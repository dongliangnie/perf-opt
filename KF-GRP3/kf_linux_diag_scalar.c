// kf_linux_diag_scalar.c —— 极致优化版卡尔曼滤波器（对角阵 + 常量折叠 + 标量化）
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>    // 健壮性: 用于断言
#include <float.h>     // 健壮性: DBL_MIN
#include "kf_linux.h"
#include "timestamp.h"

/* ---------- 内部核心计算（一次迭代） ---------- */
static int kf_core_impl(size_t dim, double* x_out) {
    if (dim == 0) return -1;

    const double sigma = 0.7071067811865475;    // sqrt(0.5)
    const double r_max = RAND_MAX / 2.0;

    double* x = (double*)calloc(dim, sizeof(double));
    double* P = (double*)malloc(dim * sizeof(double));
    if (!x || !P) {
        free(x); free(P);
        return -2;
    }

    // 初始化协方差 P = I (对角线全1)
    for (size_t i = 0; i < dim; i++)
        P[i] = 1.0;

    srand(12345);   // 固定种子，保证可重复

    // 单步迭代（step=0），所有步骤合并在一个循环中
    for (size_t i = 0; i < dim; i++) {
        // 1. 模拟测量值 z（即时使用，不存储数组）
        double true_val = (double)(i + 1);
        double noise = ((rand() / r_max) - 1.0) * sigma;
        double z = true_val + noise;

        // 2. 预测：P = P + 0.5 （过程噪声 Q=0.5）
        P[i] += 0.5;

        // 3. 卡尔曼增益 K = P / (P + 0.5) （测量噪声 R=0.5）
        /* 健壮性建议 1: 分母保护
         * 当前 P[i] 初始为 1 且 Q=0.5>0，分母恒 ≥ 0.5，不存在除零风险。
         * 但若将来修改 Q/R 为可变参数且可能为 0 或负数，则必须保护：
         *
         *   double denom = P[i] + 0.5;
         *   #ifdef NUMERIC_ROBUST
         *   if (denom < 1e-15) denom = 1e-15;   // 或 DBL_MIN
         *   #endif
         *   double K = P[i] / denom;
         */
        double K = P[i] / (P[i] + 0.5);

        // 4. 状态更新 x = x + K*(z - x)
        x[i] += K * (z - x[i]);

        // 5. 协方差更新 P = (1 - K) * P
        P[i] = (1.0 - K) * P[i];

        /* 健壮性建议 2: 协方差非负保护
         * 由于浮点舍入，(1-K)*P 可能产生极微小的负数（如 -1e-17）。
         * 若后续需要对 P 开方（如计算 sigma），负值将导致域错误。
         * 建议在调试或关键场合钳位至 0.0：
         *
         * #ifdef NUMERIC_ROBUST
         * if (P[i] < 0.0) P[i] = 0.0;
         * #endif
         */
    }

    /* 健壮性建议 3: 输出值的合理性检查（调试）
     * 可在此循环后加入断言检查 x_out 是否出现 NaN 或 Inf：
     *
     * #ifdef ENABLE_NAN_CHECK
     * for (size_t i = 0; i < dim; i++) {
     *     assert(!isnan(x[i]) && "NaN detected in state");
     *     assert(!isinf(x[i]) && "Inf detected in state");
     * }
     * #endif
     */

    // 输出最终状态
    memcpy(x_out, x, dim * sizeof(double));

    free(x);
    free(P);
    return 0;
}

/* ---------- 接口函数（类型名与 kf_linux.h 严格一致） ---------- */
RetCode kf_linux_iopointer(int dim, void* input, void* output) {
    /* 健壮性建议 4: 扩展参数校验
     * 当前仅检查 dim<=0 和空指针。可增加上限检查防止意外超大分配：
     *
     * #define MAX_DIM 4096
     * if (dim > MAX_DIM) return K_RET_INVALID_PARAM;
     */
    if (dim <= 0 || !input || !output)
        return K_RET_INVALID_PARAM;

    KfInput*  in  = (KfInput*)input;
    KfOutput* out = (KfOutput*)output;

    if (in->dim != dim || !out->x)
        return K_RET_INVALID_PARAM;

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