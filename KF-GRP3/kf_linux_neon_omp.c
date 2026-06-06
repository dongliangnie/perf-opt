// kf_linux_neon_omp.c —— NEON向量化 + OpenMP多核并行，保留 double 精度
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>      // 健壮性: 断言支持
#include <arm_neon.h>
#include <omp.h>
#include "kf_linux.h"
#include "timestamp.h"

/* ---------- 矩阵加法 ---------- */
static void mat_add(double* A, double* B, double* R, size_t n, size_t m) {
    /* 健壮性建议：A,B,R 允许重叠，但不能部分重叠导致数据竞争 */
    #pragma omp parallel for
    for (size_t i = 0; i < n * m; i++)
        R[i] = A[i] + B[i];
}

/* ---------- 矩阵减法 ---------- */
static void mat_sub(double* A, double* B, double* R, size_t n, size_t m) {
    #pragma omp parallel for
    for (size_t i = 0; i < n * m; i++)
        R[i] = A[i] - B[i];
}

/* ---------- NEON向量化 + 循环分块 + OpenMP并行的矩阵乘法 ---------- */
static void mat_mul(double* A, double* B, double* R, size_t n, size_t k, size_t m) {
    /* 健壮性建议：NEON 加载要求 16 字节对齐，若性能敏感请保证对齐
     * assert(((uintptr_t)A & 0xF) == 0);
     * assert(((uintptr_t)B & 0xF) == 0);
     * assert(((uintptr_t)R & 0xF) == 0);
     */
    memset(R, 0, sizeof(double) * n * m);
    #define BLK 64
    /* 并行策略：每个 (i0,j0) 块更新不同的 R 区域，无数据竞争 */
    #pragma omp parallel for collapse(2) schedule(static)
    for (size_t i0 = 0; i0 < n; i0 += BLK) {
        for (size_t j0 = 0; j0 < m; j0 += BLK) {
            for (size_t l0 = 0; l0 < k; l0 += BLK) {
                size_t i_end = i0 + BLK < n ? i0 + BLK : n;
                size_t j_end = j0 + BLK < m ? j0 + BLK : m;
                size_t l_end = l0 + BLK < k ? l0 + BLK : k;

                for (size_t i = i0; i < i_end; i++) {
                    for (size_t j = j0; j < j_end; j++) {
                        float64x2_t sum_v = vdupq_n_f64(0.0);
                        size_t l;
                        for (l = l0; l + 1 < l_end; l += 2) {
                            float64x2_t A_v = vld1q_f64(&A[i * k + l]);
                            /* 已修复：显式初始化 B_v 避免未定义行为 */
                            float64x2_t B_v = vdupq_n_f64(0.0);
                            B_v = vsetq_lane_f64(B[l * m + j], B_v, 0);
                            B_v = vsetq_lane_f64(B[(l+1) * m + j], B_v, 1);
                            sum_v = vfmaq_f64(sum_v, A_v, B_v);
                        }
                        double sum = vaddvq_f64(sum_v);
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
static void mat_identity(double* I, size_t dim) {
    memset(I, 0, sizeof(double) * dim * dim);
    for (size_t i = 0; i < dim; i++)
        I[i * dim + i] = 1.0;
}

/* ---------- 对角矩阵求逆（仅处理对角线） ---------- */
static void mat_inv_diag(double* A, double* R, size_t dim) {
    /* 健壮性建议：检查 A 是否确实为对角阵（仅对角线非零） */
    memset(R, 0, sizeof(double) * dim * dim);
    #pragma omp parallel for
    for (size_t i = 0; i < dim; i++) {
        double diag = A[i * dim + i];
        /* 数值保护：可启用极小分母钳位
         * #ifdef NUMERIC_ROBUST
         * if (fabs(diag) < 1e-15) diag = 1e-15;
         * #endif
         */
        if (diag != 0.0)
            R[i * dim + i] = 1.0 / diag;
    }
}

/* ---------- 卡尔曼核心计算（一次迭代） ---------- */
static int kf_core_impl(size_t dim, double* x_out) {
    if (dim == 0) return -1;

    /* 健壮性建议：防止超大维度导致分配失败或溢出
     * #define MAX_DIM 4096
     * if (dim > MAX_DIM) return -1;
     */

    size_t mat_size = dim * dim;
    size_t vec_size = dim;

    double* x       = (double*)calloc(vec_size, sizeof(double));
    double* P       = (double*)calloc(mat_size, sizeof(double));
    double* Q       = (double*)calloc(mat_size, sizeof(double));
    double* R       = (double*)calloc(mat_size, sizeof(double));
    double* K       = (double*)calloc(mat_size, sizeof(double));
    double* I_mat   = (double*)calloc(mat_size, sizeof(double));
    double* z       = (double*)calloc(vec_size, sizeof(double));
    double* true_state = (double*)calloc(vec_size, sizeof(double));
    double* temp1   = (double*)calloc(mat_size, sizeof(double));
    double* inv_temp= (double*)calloc(mat_size, sizeof(double));
    double* temp_I_minus_K = (double*)calloc(mat_size, sizeof(double));
    double* P_old   = (double*)calloc(mat_size, sizeof(double));
    double* z_minus_x = (double*)calloc(vec_size, sizeof(double));
    double* K_mult  = (double*)calloc(vec_size, sizeof(double));

    if (!x || !P || !Q || !R || !K || !I_mat || !z || !true_state ||
        !temp1 || !inv_temp || !temp_I_minus_K || !P_old || !z_minus_x || !K_mult) {
        free(x); free(P); free(Q); free(R); free(K); free(I_mat);
        free(z); free(true_state);
        free(temp1); free(inv_temp); free(temp_I_minus_K); free(P_old);
        free(z_minus_x); free(K_mult);
        return -2;
    }

    // 初始化（对角矩阵）
    for (size_t i = 0; i < dim; i++) {
        P[i * dim + i] = 1.0;
        Q[i * dim + i] = 0.5;
        R[i * dim + i] = 0.5;
        x[i] = 0.0;
    }
    mat_identity(I_mat, dim);
    srand(12345);   // 固定种子，保证可重复

    int step = 0;
    // 生成真实状态与测量（串行执行，无并发问题）
    for (size_t i = 0; i < dim; i++) {
        true_state[i] = (double)(step + 1 + i);
        double noise = ((rand() / (RAND_MAX / 2.0)) - 1.0) * sqrt(R[i * dim + i]);
        z[i] = true_state[i] + noise;
    }

    // 卡尔曼滤波步骤
    mat_add(P, Q, P, dim, dim);                 // P = P + Q
    mat_add(P, R, temp1, dim, dim);             // temp1 = P + R
    mat_inv_diag(temp1, inv_temp, dim);         // inv_temp = (P+R)^{-1}
    mat_mul(P, inv_temp, K, dim, dim, dim);     // K = P * inv_temp

    for (size_t i = 0; i < dim; i++)
        z_minus_x[i] = z[i] - x[i];
    memset(K_mult, 0, vec_size * sizeof(double));
    for (size_t i = 0; i < dim; i++)
        for (size_t j = 0; j < dim; j++)
            K_mult[i] += K[i * dim + j] * z_minus_x[j];
    for (size_t i = 0; i < dim; i++)
        x[i] += K_mult[i];

    memcpy(P_old, P, mat_size * sizeof(double));
    mat_sub(I_mat, K, temp_I_minus_K, dim, dim);
    mat_mul(temp_I_minus_K, P_old, P, dim, dim, dim);

    memcpy(x_out, x, vec_size * sizeof(double));

    /* 健壮性建议：检查输出是否包含 NaN/Inf
     * #ifdef ENABLE_NAN_CHECK
     * for (size_t i = 0; i < dim; i++) {
     *     assert(!isnan(x_out[i]) && "NaN in state");
     *     assert(!isinf(x_out[i]) && "Inf in state");
     * }
     * #endif
     */

    free(x); free(P); free(Q); free(R); free(K); free(I_mat);
    free(z); free(true_state);
    free(temp1); free(inv_temp); free(temp_I_minus_K); free(P_old);
    free(z_minus_x); free(K_mult);
    return 0;
}

/* ---------- 接口函数 ---------- */
RetCode kf_linux_iopointer(int dim, void* input, void* output) {
    /* 健壮性建议：增强参数检查，例如检查 dim 上限 */
    if (dim <= 0 || !input || !output) return K_RET_INVALID_PARAM;
    KfInput*  in  = (KfInput*)input;
    KfOutput* out = (KfOutput*)output;
    if (in->dim != dim || !out->x) return K_RET_INVALID_PARAM;
    int ret = kf_core_impl((size_t)dim, out->x);
    return (ret == 0) ? K_RET_OK : K_RET_UNK_ERROR;
}

RetCode kf_linux_ioself_profiling(int dim) {
    if (dim <= 0) return K_RET_INVALID_PARAM;
    KfInput  input;
    KfOutput output;
    input.dim = dim;
    output.x  = (double*)malloc(dim * sizeof(double));
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
    KfOutput output;
    input.dim = dim;
    output.x  = (double*)malloc(dim * sizeof(double));
    if (!output.x) return K_RET_UNK_ERROR;
    RetCode rc = kf_linux_iopointer(dim, &input, &output);
    free(output.x);
    return rc;
}