// kf_linux_diag.c —— 对角阵优化版（类型名与 kf_linux.h 对齐）
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>   // 健壮性: 用于静态/动态断言
#include <float.h>    // 健壮性: 用于浮点常量
#include "kf_linux.h"
#include "timestamp.h"

/* ---------- 卡尔曼核心计算（一次迭代），仅使用对角阵性质 ---------- */
static int kf_core_impl(size_t dim, double* x_out) {
    if (dim == 0) return -1;

    // 原先 dim×dim 的矩阵现在只用长度为 dim 的向量存储对角线
    double* x       = (double*)calloc(dim, sizeof(double));  // 状态估计
    double* P       = (double*)calloc(dim, sizeof(double));  // 协方差对角线
    double* Q       = (double*)calloc(dim, sizeof(double));  // 过程噪声对角线
    double* R       = (double*)calloc(dim, sizeof(double));  // 测量噪声对角线
    double* K       = (double*)calloc(dim, sizeof(double));  // 卡尔曼增益对角线
    double* z       = (double*)calloc(dim, sizeof(double));  // 测量值
    double* true_state = (double*)calloc(dim, sizeof(double)); // 真实状态（模拟用）

    if (!x || !P || !Q || !R || !K || !z || !true_state) {
        free(x); free(P); free(Q); free(R); free(K); free(z); free(true_state);
        return -2;
    }

    // 初始化（与原来完全一致，只是不写非对角线元素）
    for (size_t i = 0; i < dim; i++) {
        P[i] = 1.0;   // P = I
        Q[i] = 0.5;
        R[i] = 0.5;
        x[i] = 0.0;
    }
    srand(12345);  // 固定种子

    // 模拟一步
    int step = 0;
    // 1. 生成真实状态与测量（与原来完全相同）
    for (size_t i = 0; i < dim; i++) {
        true_state[i] = (double)(step + 1 + i);
        double noise = ((rand() / (RAND_MAX / 2.0)) - 1.0) * sqrt(R[i]);
        z[i] = true_state[i] + noise;
    }

    // --- 卡尔曼滤波核心步骤（对角阵版本）---

    /* 健壮性建议 1: 运行时检查对角性前提（调试阶段启用）
     * 如果将来有人错误使用了非对角矩阵初始化，这里可以尽早暴露问题。
     * 注意：本实现已默认只存对角线，因此该检查在当前代码中并非必需，
     * 但若你的外部接口仍可能传入完整矩阵，则强烈建议加入。
     */
    #ifdef ENABLE_DIAGONAL_CHECK
    {
        // 示例：如果是从外部完整矩阵转换而来，可以这样检查非对角线元素
        // for (size_t i = 0; i < dim; i++) {
        //     for (size_t j = 0; j < dim; j++) {
        //         if (i != j) {
        //             assert(fabs(matrix[i*dim + j]) < 1e-15 && "Non-diagonal element detected!");
        //         }
        //     }
        // }
    }
    #endif

    // 预测：P = P + Q
    for (size_t i = 0; i < dim; i++)
        P[i] += Q[i];

    // 卡尔曼增益：K = P * inv(P + R)，对角阵运算退化为逐元素
    for (size_t i = 0; i < dim; i++) {
        double temp = P[i] + R[i];
        /* 健壮性建议 2: 分母保护与数值稳定性
         * 避免除零和极小值导致的数值问题。
         * 建议用 DBL_MIN 或一个合理的下限（如 1e-15）钳位。
         */
        // 原始行:
        // double inv_temp = (temp != 0.0) ? 1.0 / temp : 0.0;
        // 健壮版本:
        #ifdef NUMERIC_ROBUST
        const double MIN_DENOM = 1e-15;   // 根据实际噪声水平调整
        if (temp < MIN_DENOM) {
            temp = MIN_DENOM;
        }
        #endif
        double inv_temp = (temp != 0.0) ? 1.0 / temp : 0.0;
        K[i] = P[i] * inv_temp;
    }

    // 状态更新：x = x + K * (z - x) ，对角阵 K 乘向量退化为逐元素乘
    for (size_t i = 0; i < dim; i++)
        x[i] += K[i] * (z[i] - x[i]);

    // 协方差更新：P = (I - K) * P，对角阵版本
    for (size_t i = 0; i < dim; i++) {
        P[i] = (1.0 - K[i]) * P[i];
        /* 健壮性建议 3: 协方差非负保护
         * 由于浮点舍入，P[i] 可能变为极小的负数，导致后续 sqrt 出错或滤波器锁死。
         * 建议钳位至 0.0。
         */
        #ifdef NUMERIC_ROBUST
        if (P[i] < 0.0) {
            P[i] = 0.0;
        }
        #endif
    }

    // 输出最终状态
    memcpy(x_out, x, dim * sizeof(double));

    // 释放内存
    free(x); free(P); free(Q); free(R); free(K); free(z); free(true_state);
    return 0;
}

/* ---------- 接口函数实现（类型名与 kf_linux.h 保持一致）---------- */
RetCode kf_linux_iopointer(int dim, void *input, void *output) {
    /* 健壮性建议 4: 接口参数强校验
     * 现有的 dim <= 0 和空指针检查已较完善，可进一步检查 dim 是否超过合理上限，
     * 避免栈/堆溢出（如果后续改为固定长度数组）。
     */
    if (dim <= 0 || !input || !output)
        return K_RET_INVALID_PARAM;

    // 可扩展: 检查 dim 是否在安全范围内
    // #define MAX_DIM 2048
    // if (dim > MAX_DIM) return K_RET_INVALID_PARAM;

    KfInput  *in  = (KfInput*)input;
    KfOutput *out = (KfOutput*)output;

    if (in->dim != dim || !out->x)
        return K_RET_INVALID_PARAM;

    /* 健壮性建议 5: 输入一致性检查
     * 如果 KfInput 包含协方差矩阵等，需确保它们为对角阵。此处仅在开发阶段启用断言。
     */
    #ifdef ENABLE_DIAGONAL_CHECK
    // 假设输入中可能携带原始矩阵，可在此检查
    // assert(is_diagonal(in->P, dim) && "Input P must be diagonal");
    #endif

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