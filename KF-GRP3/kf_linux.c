#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "kf_linux.h"
#include "timestamp.h"

/* ---------- 内部矩阵运算（蛇形命名，参数全小写）---------- */
static void mat_add(double *a, double *b, double *result, size_t n, size_t m) {
    for (size_t i = 0; i < n * m; i++)
        result[i] = a[i] + b[i];
}

static void mat_sub(double *a, double *b, double *result, size_t n, size_t m) {
    for (size_t i = 0; i < n * m; i++)
        result[i] = a[i] - b[i];
}

static void mat_mul(double *a, double *b, double *result, size_t n, size_t k, size_t m) {
    memset(result, 0, sizeof(double) * n * m);
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < m; j++) {
            double sum = 0.0;
            for (size_t l = 0; l < k; l++)
                sum += a[i * k + l] * b[l * m + j];
            result[i * m + j] = sum;
        }
}

static void mat_identity(double *imat, size_t dim) {
    memset(imat, 0, sizeof(double) * dim * dim);
    for (size_t i = 0; i < dim; i++)
        imat[i * dim + i] = 1.0;
}

static void mat_inv_diag(double *a, double *result, size_t dim) {
    memset(result, 0, sizeof(double) * dim * dim);
    for (size_t i = 0; i < dim; i++)
        if (a[i * dim + i] != 0.0)
            result[i * dim + i] = 1.0 / a[i * dim + i];
}

/* ---------- 卡尔曼核心实现（变量全改为蛇形小写）---------- */
static int kf_core_impl(size_t dim, double *x_out) {
    if (dim == 0) return -1;

    double *x       = calloc(dim, sizeof(double));
    double *p       = calloc(dim * dim, sizeof(double));
    double *q       = calloc(dim * dim, sizeof(double));
    double *r       = calloc(dim * dim, sizeof(double));
    double *k       = calloc(dim * dim, sizeof(double));
    double *i_mat   = calloc(dim * dim, sizeof(double));
    double *z       = calloc(dim, sizeof(double));
    double *true_state = calloc(dim, sizeof(double));
    double *temp_p_plus_r = calloc(dim * dim, sizeof(double));
    double *inv_temp = calloc(dim * dim, sizeof(double));
    double *temp_i_minus_k = calloc(dim * dim, sizeof(double));
    double *p_old    = calloc(dim * dim, sizeof(double));
    double *z_minus_x = calloc(dim, sizeof(double));
    double *k_mult   = calloc(dim, sizeof(double));

    if (!x || !p || !q || !r || !k || !i_mat || !z || !true_state ||
        !temp_p_plus_r || !inv_temp || !temp_i_minus_k || !p_old ||
        !z_minus_x || !k_mult) {
        free(x); free(p); free(q); free(r); free(k); free(i_mat); free(z);
        free(true_state);
        free(temp_p_plus_r); free(inv_temp); free(temp_i_minus_k); free(p_old);
        free(z_minus_x); free(k_mult);
        return -2;   /* 内存分配失败 */
    }

    /* 初始化矩阵 */
    for (size_t i = 0; i < dim; i++) {
        p[i * dim + i] = 1.0;
        q[i * dim + i] = 0.5;
        r[i * dim + i] = 0.5;
        x[i] = 0.0;
    }
    mat_identity(i_mat, dim);
    srand(12345);   /* 固定随机种子 */

    int step = 0;
    for (size_t i = 0; i < dim; i++) {
        true_state[i] = (double)(step + 1 + i);
        z[i] = true_state[i] + ((rand() / (RAND_MAX / 2.0)) - 1.0) * sqrt(r[i * dim + i]);
    }

    /* 预测：p = p + q */
    mat_add(p, q, p, dim, dim);
    /* 更新：k = p * inv(p + r) */
    mat_add(p, r, temp_p_plus_r, dim, dim);
    mat_inv_diag(temp_p_plus_r, inv_temp, dim);
    mat_mul(p, inv_temp, k, dim, dim, dim);
    /* 状态更新：x = x + k * (z - x) */
    for (size_t i = 0; i < dim; i++)
        z_minus_x[i] = z[i] - x[i];
    memset(k_mult, 0, dim * sizeof(double));
    for (size_t i = 0; i < dim; i++)
        for (size_t j = 0; j < dim; j++)
            k_mult[i] += k[i * dim + j] * z_minus_x[j];
    for (size_t i = 0; i < dim; i++)
        x[i] += k_mult[i];
    /* 协方差更新：p = (i_mat - k) * p */
    memcpy(p_old, p, dim * dim * sizeof(double));
    mat_sub(i_mat, k, temp_i_minus_k, dim, dim);
    mat_mul(temp_i_minus_k, p_old, p, dim, dim, dim);

    memcpy(x_out, x, dim * sizeof(double));

    free(x); free(p); free(q); free(r); free(k); free(i_mat); free(z);
    free(true_state);
    free(temp_p_plus_r); free(inv_temp); free(temp_i_minus_k); free(p_old);
    free(z_minus_x); free(k_mult);
    return 0;
}

/* ---------- 接口函数（类型名使用大驼峰，变量蛇形）---------- */
RetCode kf_linux_iopointer(int dim, void *input, void *output) {
    if (dim <= 0 || !input || !output)
        return K_RET_INVALID_PARAM;

    KfInput  *in  = (KfInput *)input;
    KfOutput *out = (KfOutput *)output;

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
    output.x  = malloc(dim * sizeof(double));
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
    output.x  = malloc(dim * sizeof(double));
    if (!output.x) return K_RET_UNK_ERROR;

    RetCode rc = kf_linux_iopointer(dim, &input, &output);

    free(output.x);
    return rc;
}