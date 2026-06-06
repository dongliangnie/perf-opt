// kf_linux_neon.c —— NEON 向量化矩阵乘法 + 循环分块
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>      // 健壮性: 断言
#include <arm_neon.h>
#include "kf_linux.h"
#include "timestamp.h"

/* ---------- 基本矩阵运算 ---------- */
static void mat_add(double* A, double* B, double* R, size_t n, size_t m) {
    for (size_t i = 0; i < n * m; i++)
        R[i] = A[i] + B[i];
}

static void mat_sub(double* A, double* B, double* R, size_t n, size_t m) {
    for (size_t i = 0; i < n * m; i++)
        R[i] = A[i] - B[i];
}

static void mat_mul(double* A, double* B, double* R, size_t n, size_t k, size_t m) {
    /* 健壮性建议 0：确保矩阵地址满足 NEON 对齐要求
     * vld1q_f64 要求 16 字节对齐。若指针未对齐，在部分平台会触发总线错误，
     * 或在支持非对齐访问的平台造成性能下降。
     * 建议：使用 aligned_alloc(16, ...) 或 posix_memalign 分配矩阵内存，
     * 并在函数入口处添加断言：
     *   assert(((uintptr_t)A & 0xF) == 0 && "A not 16-byte aligned");
     *   assert(((uintptr_t)B & 0xF) == 0);
     *   assert(((uintptr_t)R & 0xF) == 0);
     * 此处假定已满足对齐（如通过编译器选项或分配器保证）。
     */

    memset(R, 0, sizeof(double) * n * m);
    #define BLK 64
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
                            // B 的列访问不连续，手动构造向量
                            /* 健壮性建议 1：初始化 B_v 以避免未定义行为
                             * 原代码：float64x2_t B_v; 
                             *          B_v = vsetq_lane_f64(..., B_v, 0);  // B_v 未初始化！
                             * 这属于未定义行为。应显式初始化为零或使用其他安全方式。
                             * 修复如下：
                             */
                            float64x2_t B_v = vdupq_n_f64(0.0); // 显式初始化
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

static void mat_identity(double* I, size_t dim) {
    memset(I, 0, sizeof(double) * dim * dim);
    for (size_t i = 0; i < dim; i++)
        I[i * dim + i] = 1.0;
}

static void mat_inv_diag(double* A, double* R, size_t dim) {
    /* 健壮性建议 2：明确对角假设并保护除零
     * 该函数仅对对角线元素求倒数，隐含假设 A 为对角阵（或仅对角线有意义）。
     * 若调用方误传非对角阵，会得到错误结果且无警告。
     * 建议在调试版添加断言，并引入最小阈值避免除以过小数。
     */
    memset(R, 0, sizeof(double) * dim * dim);
    for (size_t i = 0; i < dim; i++) {
        double diag = A[i * dim + i];
        #ifdef NUMERIC_ROBUST
        if (fabs(diag) < 1e-15) {
            diag = 1e-15;   // 或直接报错
        }
        #endif
        if (diag != 0.0)
            R[i * dim + i] = 1.0 / diag;
    }
}

/* ---------- 卡尔曼核心 ---------- */
static int kf_core_impl(size_t dim, double* x_out) {
    if (dim == 0) return -1;

    size_t mat_size = dim * dim;
    size_t vec_size = dim;

    /* 健壮性建议 3：分配前检查 dim 上限，避免 size_t 溢出 */
    #ifdef SIZE_CHECK
    if (dim > 4096) return -1; // 示例上限
    #endif

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

    // 初始化（仅对角线非零，符合对角假设）
    for (size_t i = 0; i < dim; i++) {
        P[i * dim + i] = 1.0;
        Q[i * dim + i] = 0.5;
        R[i * dim + i] = 0.5;
        x[i] = 0.0;
    }
    mat_identity(I_mat, dim);
    srand(12345);

    int step = 0;
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

    /* 健壮性建议 4：检查输出是否包含 NaN/Inf */
    #ifdef ENABLE_NAN_CHECK
    for (size_t i = 0; i < dim; i++) {
        assert(!isnan(x_out[i]) && "NaN in state");
        assert(!isinf(x_out[i]) && "Inf in state");
    }
    #endif

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