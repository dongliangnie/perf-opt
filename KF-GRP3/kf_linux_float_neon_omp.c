// kf_linux_float_neon_omp.c —— 全 float 实现 + NEON + OpenMP
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <arm_neon.h>
#include <omp.h>
#include "kf_linux.h"
#include "timestamp.h"

/* ---------- 类型定义与接口说明 ---------- */
typedef float real_t;

/* 健壮性提醒：该文件为 float 专用版本，输出结构体 KfOutput_float 与标准
 * KfOutput (double* x) 内存布局不同。调用者必须确保传递正确的结构体，
 * 否则会导致内存破坏。建议在项目中统一使用一种接口，或通过编译开关切换。
 */
typedef struct {
    real_t *x;
} KfOutput_float;

/* ---------- 矩阵加法 ---------- */
static void mat_add(real_t* A, real_t* B, real_t* R, size_t n, size_t m) {
    #pragma omp parallel for
    for (size_t i = 0; i < n * m; i++)
        R[i] = A[i] + B[i];
}

/* ---------- 矩阵减法 ---------- */
static void mat_sub(real_t* A, real_t* B, real_t* R, size_t n, size_t m) {
    #pragma omp parallel for
    for (size_t i = 0; i < n * m; i++)
        R[i] = A[i] - B[i];
}

/* ---------- NEON向量化 + 分块 + OpenMP 并行的 float 矩阵乘法 ---------- */
static void mat_mul(real_t* A, real_t* B, real_t* R, size_t n, size_t k, size_t m) {
    /* 健壮性建议：检查 A,B,R 是否 16 字节对齐（同上） */
    memset(R, 0, sizeof(real_t) * n * m);
    #define BLK 64
    #pragma omp parallel for collapse(2) schedule(static)
    for (size_t i0 = 0; i0 < n; i0 += BLK) {
        for (size_t j0 = 0; j0 < m; j0 += BLK) {
            for (size_t l0 = 0; l0 < k; l0 += BLK) {
                size_t i_end = i0 + BLK < n ? i0 + BLK : n;
                size_t j_end = j0 + BLK < m ? j0 + BLK : m;
                size_t l_end = l0 + BLK < k ? l0 + BLK : k;

                for (size_t i = i0; i < i_end; i++) {
                    for (size_t j = j0; j < j_end; j++) {
                        float32x4_t sum_v = vdupq_n_f32(0.0f);
                        size_t l;
                        for (l = l0; l + 3 < l_end; l += 4) {
                            float32x4_t A_v = vld1q_f32(&A[i * k + l]);
                            /* 已修复：显式初始化 B_v 避免未定义行为 */
                            float32x4_t B_v = vdupq_n_f32(0.0f);
                            B_v = vsetq_lane_f32(B[l * m + j], B_v, 0);
                            B_v = vsetq_lane_f32(B[(l+1) * m + j], B_v, 1);
                            B_v = vsetq_lane_f32(B[(l+2) * m + j], B_v, 2);
                            B_v = vsetq_lane_f32(B[(l+3) * m + j], B_v, 3);
                            sum_v = vfmaq_f32(sum_v, A_v, B_v);
                        }
                        /* 健壮性提醒：vaddvq_f32 仅在 AArch64 下可用。
                         * 若需兼容 ARMv7，可使用手动水平相加：
                         * float sum = vgetq_lane_f32(sum_v,0) +
                         *              vgetq_lane_f32(sum_v,1) +
                         *              vgetq_lane_f32(sum_v,2) +
                         *              vgetq_lane_f32(sum_v,3);
                         */
                        float sum = vaddvq_f32(sum_v);
                        for (; l < l_end; l++)
                            sum += A[i * k + l] * B[l * m + j];
                        R[i * m + j] += sum;
                    }
                }
            }
        }
    }
    #undef BLK
}

/* ---------- 生成单位矩阵 ---------- */
static void mat_identity(real_t* I, size_t dim) {
    memset(I, 0, sizeof(real_t) * dim * dim);
    for (size_t i = 0; i < dim; i++)
        I[i * dim + i] = 1.0f;
}

/* ---------- 对角矩阵求逆 ---------- */
static void mat_inv_diag(real_t* A, real_t* R, size_t dim) {
    memset(R, 0, sizeof(real_t) * dim * dim);
    #pragma omp parallel for
    for (size_t i = 0; i < dim; i++) {
        real_t diag = A[i * dim + i];
        /* 数值保护可在此加入钳位，同 double 版 */
        if (diag != 0.0f)
            R[i * dim + i] = 1.0f / diag;
    }
}

/* ---------- 卡尔曼核心（全部使用 float） ---------- */
static int kf_core_impl(size_t dim, real_t* x_out) {
    if (dim == 0) return -1;

    size_t mat_size = dim * dim;
    size_t vec_size = dim;

    real_t* x       = (real_t*)calloc(vec_size, sizeof(real_t));
    real_t* P       = (real_t*)calloc(mat_size, sizeof(real_t));
    real_t* Q       = (real_t*)calloc(mat_size, sizeof(real_t));
    real_t* R       = (real_t*)calloc(mat_size, sizeof(real_t));
    real_t* K       = (real_t*)calloc(mat_size, sizeof(real_t));
    real_t* I_mat   = (real_t*)calloc(mat_size, sizeof(real_t));
    real_t* z       = (real_t*)calloc(vec_size, sizeof(real_t));
    real_t* true_state = (real_t*)calloc(vec_size, sizeof(real_t));
    real_t* temp1   = (real_t*)calloc(mat_size, sizeof(real_t));
    real_t* inv_temp= (real_t*)calloc(mat_size, sizeof(real_t));
    real_t* temp_I_minus_K = (real_t*)calloc(mat_size, sizeof(real_t));
    real_t* P_old   = (real_t*)calloc(mat_size, sizeof(real_t));
    real_t* z_minus_x = (real_t*)calloc(vec_size, sizeof(real_t));
    real_t* K_mult  = (real_t*)calloc(vec_size, sizeof(real_t));

    if (!x || !P || !Q || !R || !K || !I_mat || !z || !true_state ||
        !temp1 || !inv_temp || !temp_I_minus_K || !P_old || !z_minus_x || !K_mult) {
        free(x); free(P); free(Q); free(R); free(K); free(I_mat);
        free(z); free(true_state);
        free(temp1); free(inv_temp); free(temp_I_minus_K); free(P_old);
        free(z_minus_x); free(K_mult);
        return -2;
    }

    for (size_t i = 0; i < dim; i++) {
        P[i * dim + i] = 1.0f;
        Q[i * dim + i] = 0.5f;
        R[i * dim + i] = 0.5f;
        x[i] = 0.0f;
    }
    mat_identity(I_mat, dim);
    srand(12345);

    int step = 0;
    for (size_t i = 0; i < dim; i++) {
        true_state[i] = (real_t)(step + 1 + i);
        real_t noise = ((rand() / (RAND_MAX / 2.0f)) - 1.0f) * sqrtf(R[i * dim + i]);
        z[i] = true_state[i] + noise;
    }

    mat_add(P, Q, P, dim, dim);
    mat_add(P, R, temp1, dim, dim);
    mat_inv_diag(temp1, inv_temp, dim);
    mat_mul(P, inv_temp, K, dim, dim, dim);

    for (size_t i = 0; i < dim; i++)
        z_minus_x[i] = z[i] - x[i];
    memset(K_mult, 0, vec_size * sizeof(real_t));
    for (size_t i = 0; i < dim; i++)
        for (size_t j = 0; j < dim; j++)
            K_mult[i] += K[i * dim + j] * z_minus_x[j];
    for (size_t i = 0; i < dim; i++)
        x[i] += K_mult[i];

    memcpy(P_old, P, mat_size * sizeof(real_t));
    mat_sub(I_mat, K, temp_I_minus_K, dim, dim);
    mat_mul(temp_I_minus_K, P_old, P, dim, dim, dim);

    memcpy(x_out, x, vec_size * sizeof(real_t));

    free(x); free(P); free(Q); free(R); free(K); free(I_mat);
    free(z); free(true_state);
    free(temp1); free(inv_temp); free(temp_I_minus_K); free(P_old);
    free(z_minus_x); free(K_mult);
    return 0;
}

/* ---------- 接口函数 ---------- */
RetCode kf_linux_iopointer(int dim, void* input, void* output) {
    if (dim <= 0 || !input || !output) return K_RET_INVALID_PARAM;
    KfInput*  in  = (KfInput*)input;
    KfOutput_float* out = (KfOutput_float*)output;  /* 注意：必须传入 KfOutput_float* */
    if (in->dim != dim || !out->x) return K_RET_INVALID_PARAM;
    int ret = kf_core_impl((size_t)dim, out->x);
    return (ret == 0) ? K_RET_OK : K_RET_UNK_ERROR;
}

RetCode kf_linux_ioself_profiling(int dim) {
    if (dim <= 0) return K_RET_INVALID_PARAM;
    KfInput  input;
    KfOutput_float output;
    input.dim = dim;
    output.x  = (real_t*)malloc(dim * sizeof(real_t));
    if (!output.x) return K_RET_UNK_ERROR;
    uint64_t t0 = get_timestamp();
    RetCode rc = kf_linux_iopointer(dim, &input, &output);
    uint64_t t1 = get_timestamp();
    if (rc == K_RET_OK) {
        int64_t ns = timestamp_diff(t0, t1);
        printf("dim=%d, time=%ld ns\n", dim, ns);
    }
    free(output.x);
    return rc;
}

RetCode kf_linux_ioself(int dim) {
    if (dim <= 0) return K_RET_INVALID_PARAM;
    KfInput  input;
    KfOutput_float output;
    input.dim = dim;
    output.x  = (real_t*)malloc(dim * sizeof(real_t));
    if (!output.x) return K_RET_UNK_ERROR;
    RetCode rc = kf_linux_iopointer(dim, &input, &output);
    free(output.x);
    return rc;
}