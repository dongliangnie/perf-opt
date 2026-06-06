// kf_linux_float_diag_scalar.c —— 对角阵 + 常量折叠 + 标量化 + float
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>    // 健壮性
#include <float.h>     // 健壮性: FLT_MIN, FLT_EPSILON
#include "kf_linux.h"
#include "timestamp.h"

/* ---------- 内部核心计算（一次迭代），全部使用 float ---------- */
static int kf_core_impl(size_t dim, double* x_out) {
    if (dim == 0) return -1;

    const float sigma = 0.7071067811865475f;    // sqrt(0.5)
    const float r_max = RAND_MAX / 2.0f;

    float* x = (float*)calloc(dim, sizeof(float));
    float* P = (float*)malloc(dim * sizeof(float));
    if (!x || !P) {
        free(x); free(P);
        return -2;
    }

    // 初始化协方差 P = I (对角线全1)
    for (size_t i = 0; i < dim; i++)
        P[i] = 1.0f;

    srand(12345);   // 固定种子

    // 单步迭代，所有步骤合并在一个循环中
    for (size_t i = 0; i < dim; i++) {
        // 1. 模拟测量值 z
        float true_val = (float)(i + 1);
        /* 健壮性建议 1: 噪声生成注意 float 精度
         * rand() 返回 int，r_max 为 float，计算中可能产生微小的截断误差。
         * 对于高维或长时间运行，这些误差累积可能影响滤波器一致性。
         * 建议对关键系统使用 double 进行噪声生成再截断，或记录噪声统计。
         */
        float noise = ((rand() / r_max) - 1.0f) * sigma;
        float z = true_val + noise;

        // 2. 预测：P = P + 0.5
        P[i] += 0.5f;

        // 3. 卡尔曼增益 K = P / (P + 0.5)
        /* 健壮性建议 2: 分母保护 (float 版)
         * 虽然当前参数下分母安全，但若改为可变参数，注意除零。
         * float 的极小数可能导致 K 计算为 nan。
         */
        float K = P[i] / (P[i] + 0.5f);

        // 4. 状态更新 x = x + K*(z - x)
        x[i] += K * (z - x[i]);

        // 5. 协方差更新 P = (1 - K) * P
        P[i] = (1.0f - K) * P[i];

        /* 健壮性建议 3: 协方差非负保护 (float 版)
         * float 尾数仅 24 位，减法相消更容易产生负零或微小负数。
         * 强烈建议启用钳位，尤其在迭代次数多时。
         *
         * #ifdef NUMERIC_ROBUST
         * if (P[i] < 0.0f) P[i] = 0.0f;
         * #endif
         */
    }

    /* 健壮性建议 4: float→double 转换与精度检查
     * 从 float 转为 double 不会提升精度，但可以避免后续运算中的额外舍入。
     * 若下游期望 double 精度，应确保 float 版本的计算误差在可接受范围内。
     * 建议在调试阶段对比 double 版本的状态输出差异。
     */
    for (size_t i = 0; i < dim; i++)
        x_out[i] = (double)x[i];

    /* 健壮性建议 5: NaN/Inf 检测
     * 同 double 版，可在输出前循环检查 x 数组。
     */
    #ifdef ENABLE_NAN_CHECK
    for (size_t i = 0; i < dim; i++) {
        assert(!isnan(x[i]) && "NaN detected in state (float)");
        assert(!isinf(x[i]) && "Inf detected in state (float)");
    }
    #endif

    free(x);
    free(P);
    return 0;
}

/* ---------- 接口函数（与之前完全相同） ---------- */
RetCode kf_linux_iopointer(int dim, void* input, void* output) {
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