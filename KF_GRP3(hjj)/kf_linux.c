#include "kf_linux.h"
#include "timestamp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void locality_optimized_mat_mul(const double *A, const double *B, double *result, size_t n, size_t k_common, size_t m)
{
    memset(result, 0, sizeof(double) * n * m);

    for (size_t i = 0; i < n; i++) {
        double *result_row = result + (i * m);
        const double *a_row = A + (i * k_common);

        for (size_t l = 0; l < k_common; l++) {
            double a_val = a_row[l];
            const double *b_row = B + (l * m);

            for (size_t j = 0; j < m; j++) {
                result_row[j] += a_val * b_row[j];
            }
        }
    }
}

static void optimized_mat_inv_diag(const double *A, double *result, size_t dim)
{
    memset(result, 0, sizeof(double) * dim * dim);
    for (size_t i = 0; i < dim; i++) {
        if (A[i * dim + i] != 0.0) {
            result[i * dim + i] = 1.0 / A[i * dim + i];
        }
    }
}

static void mat_add_inplace(double *A, const double *B, size_t n, size_t m)
{
    for (size_t i = 0; i < n * m; i++) {
        A[i] += B[i];
    }
}

static void mat_add(const double *A, const double *B, double *result, size_t n, size_t m)
{
    for (size_t i = 0; i < n * m; i++) {
        result[i] = A[i] + B[i];
    }
}

static void mat_sub(const double *A, const double *B, double *result, size_t n, size_t m)
{
    for (size_t i = 0; i < n * m; i++) {
        result[i] = A[i] - B[i];
    }
}

static void mat_copy(const double *src, double *dst, size_t n, size_t m)
{
    memcpy(dst, src, sizeof(double) * n * m);
}

static void mat_identity(double *I, size_t dim)
{
    memset(I, 0, sizeof(double) * dim * dim);
    for (size_t i = 0; i < dim; i++) {
        I[i * dim + i] = 1.0;
    }
}

static void mat_vec_mul(const double *mat, const double *vec, double *result, size_t dim)
{
    memset(result, 0, sizeof(double) * dim);
    for (size_t i = 0; i < dim; i++) {
        const double *mat_row = mat + (i * dim);
        double sum = 0.0;
        for (size_t j = 0; j < dim; j++) {
            sum += mat_row[j] * vec[j];
        }
        result[i] = sum;
    }
}

static int kf_core_compute(
    size_t dim,
    int steps,
    double *x,
    double *P,
    const double *Q,
    const double *R,
    double *K,
    const double *I_mat,
    double *z,
    double *temp_p_plus_r,
    double *inv_temp_p_plus_r,
    double *temp_i_minus_k,
    double *p_old,
    double *z_minus_x,
    double *k_mult_z_minus_x)
{
    if (dim == 0 || steps == 0) {
        return K_RET_INVALID_PARAM;
    }

    double *true_state_sim = (double *)calloc(dim, sizeof(double));
    if (!true_state_sim) {
        return K_RET_UNK_ERROR;
    }

    for (int step = 0; step < steps; step++) {
        for (size_t i = 0; i < dim; i++) {
            true_state_sim[i] = (double)(step + 1 + i);
            z[i] = true_state_sim[i] + ((rand() / (RAND_MAX / 2.0)) - 1.0) * sqrt(0.5);
        }

        mat_add_inplace(P, Q, dim, dim);

        mat_add(P, R, temp_p_plus_r, dim, dim);

        optimized_mat_inv_diag(temp_p_plus_r, inv_temp_p_plus_r, dim);

        locality_optimized_mat_mul(P, inv_temp_p_plus_r, K, dim, dim, dim);

        for (size_t i = 0; i < dim; i++) {
            z_minus_x[i] = z[i] - x[i];
        }

        mat_vec_mul(K, z_minus_x, k_mult_z_minus_x, dim);

        for (size_t i = 0; i < dim; i++) {
            x[i] += k_mult_z_minus_x[i];
        }

        mat_copy(P, p_old, dim, dim);
        mat_sub(I_mat, K, temp_i_minus_k, dim, dim);
        locality_optimized_mat_mul(temp_i_minus_k, p_old, P, dim, dim, dim);
    }

    free(true_state_sim);
    return K_RET_OK;
}

static RetCode run_one_profile(int horizon)
{
    size_t dim = (size_t)horizon;
    int steps = 1;

    double *x = (double *)calloc(dim, sizeof(double));
    double *P = (double *)calloc(dim * dim, sizeof(double));
    double *Q = (double *)calloc(dim * dim, sizeof(double));
    double *R = (double *)calloc(dim * dim, sizeof(double));
    double *K = (double *)calloc(dim * dim, sizeof(double));
    double *I_mat = (double *)calloc(dim * dim, sizeof(double));
    double *z = (double *)calloc(dim, sizeof(double));

    double *temp_p_plus_r = (double *)calloc(dim * dim, sizeof(double));
    double *inv_temp_p_plus_r = (double *)calloc(dim * dim, sizeof(double));
    double *temp_i_minus_k = (double *)calloc(dim * dim, sizeof(double));
    double *p_old = (double *)calloc(dim * dim, sizeof(double));
    double *z_minus_x = (double *)calloc(dim, sizeof(double));
    double *k_mult_z_minus_x = (double *)calloc(dim, sizeof(double));

    if (!x || !P || !Q || !R || !K || !I_mat || !z ||
        !temp_p_plus_r || !inv_temp_p_plus_r || !temp_i_minus_k || !p_old ||
        !z_minus_x || !k_mult_z_minus_x) {
        free(x); free(P); free(Q); free(R); free(K); free(I_mat); free(z);
        free(temp_p_plus_r); free(inv_temp_p_plus_r); free(temp_i_minus_k);
        free(p_old); free(z_minus_x); free(k_mult_z_minus_x);
        return K_RET_UNK_ERROR;
    }

    for (size_t i = 0; i < dim; i++) {
        P[i * dim + i] = 1.0;
        Q[i * dim + i] = 0.5;
        R[i * dim + i] = 0.5;
        I_mat[i * dim + i] = 1.0;
    }
    srand(12345);

    TimeStamp start = timestamp();

    int ret = kf_core_compute(dim, steps, x, P, Q, R, K, I_mat, z,
                               temp_p_plus_r, inv_temp_p_plus_r, temp_i_minus_k,
                               p_old, z_minus_x, k_mult_z_minus_x);

    TimeStamp end = timestamp();
    int64_t elapsed_ns = timestamp_diff(start, end);

    printf("%d,%.6f,A72,KF\n", horizon, (double)elapsed_ns / 1000000.0);

    free(x); free(P); free(Q); free(R); free(K); free(I_mat); free(z);
    free(temp_p_plus_r); free(inv_temp_p_plus_r); free(temp_i_minus_k);
    free(p_old); free(z_minus_x); free(k_mult_z_minus_x);

    return (RetCode)ret;
}

static RetCode kf_self_test(int horizon, int enable_timing)
{
    if (horizon <= 0 || horizon > KF_MAX_DIM) {
        return K_RET_INVALID_PARAM;
    }

    size_t dim = (size_t)horizon;
    int steps = 1;

    double *x = (double *)calloc(dim, sizeof(double));
    double *P = (double *)calloc(dim * dim, sizeof(double));
    double *Q = (double *)calloc(dim * dim, sizeof(double));
    double *R = (double *)calloc(dim * dim, sizeof(double));
    double *K = (double *)calloc(dim * dim, sizeof(double));
    double *I_mat = (double *)calloc(dim * dim, sizeof(double));
    double *z = (double *)calloc(dim, sizeof(double));

    double *temp_p_plus_r = (double *)calloc(dim * dim, sizeof(double));
    double *inv_temp_p_plus_r = (double *)calloc(dim * dim, sizeof(double));
    double *temp_i_minus_k = (double *)calloc(dim * dim, sizeof(double));
    double *p_old = (double *)calloc(dim * dim, sizeof(double));
    double *z_minus_x = (double *)calloc(dim, sizeof(double));
    double *k_mult_z_minus_x = (double *)calloc(dim, sizeof(double));

    if (!x || !P || !Q || !R || !K || !I_mat || !z ||
        !temp_p_plus_r || !inv_temp_p_plus_r || !temp_i_minus_k || !p_old ||
        !z_minus_x || !k_mult_z_minus_x) {
        free(x); free(P); free(Q); free(R); free(K); free(I_mat); free(z);
        free(temp_p_plus_r); free(inv_temp_p_plus_r); free(temp_i_minus_k);
        free(p_old); free(z_minus_x); free(k_mult_z_minus_x);
        return K_RET_UNK_ERROR;
    }

    for (size_t i = 0; i < dim; i++) {
        P[i * dim + i] = 1.0;
        Q[i * dim + i] = 0.5;
        R[i * dim + i] = 0.5;
        I_mat[i * dim + i] = 1.0;
    }
    srand(12345);

    RetCode ret = K_RET_OK;

    if (enable_timing) {
        TimeStamp start = timestamp();
        ret = (RetCode)kf_core_compute(dim, steps, x, P, Q, R, K, I_mat, z,
                                        temp_p_plus_r, inv_temp_p_plus_r,
                                        temp_i_minus_k, p_old, z_minus_x,
                                        k_mult_z_minus_x);
        TimeStamp end = timestamp();
        int64_t elapsed_ns = timestamp_diff(start, end);
        printf("%d,%.6f,A72,KF\n", horizon, (double)elapsed_ns / 1000000.0);
    } else {
        ret = (RetCode)kf_core_compute(dim, steps, x, P, Q, R, K, I_mat, z,
                                        temp_p_plus_r, inv_temp_p_plus_r,
                                        temp_i_minus_k, p_old, z_minus_x,
                                        k_mult_z_minus_x);
    }

    free(x); free(P); free(Q); free(R); free(K); free(I_mat); free(z);
    free(temp_p_plus_r); free(inv_temp_p_plus_r); free(temp_i_minus_k);
    free(p_old); free(z_minus_x); free(k_mult_z_minus_x);

    return ret;
}

RetCode kf_linux_iopointer(int horizon, void *input, void *output)
{
    if (horizon <= 0 || horizon > KF_MAX_DIM) {
        return K_RET_INVALID_PARAM;
    }
    if (input == NULL || output == NULL) {
        return K_RET_INVALID_PARAM;
    }

    KfInput *kf_in = (KfInput *)input;
    KfOutput *kf_out = (KfOutput *)output;

    if (kf_in->dim != horizon) {
        return K_RET_INVALID_PARAM;
    }

    if (kf_in->mode == KF_MODE_SOA_NEON) {
        if (!kf_in->x || !kf_in->p || !kf_in->q || !kf_in->r || !kf_in->k ||
            !kf_in->measurements || !kf_in->i_mat) {
            return K_RET_INVALID_PARAM;
        }

        int ret = kf_core_compute(
            (size_t)horizon,
            kf_in->steps_as_iterations,
            kf_in->x,
            kf_in->p,
            kf_in->q,
            kf_in->r,
            kf_in->k,
            kf_in->i_mat,
            kf_in->measurements,
            kf_in->temp_p_plus_r,
            kf_in->inv_temp_p_plus_r,
            kf_in->temp_i_minus_k,
            kf_in->p_old,
            kf_in->z_minus_x,
            kf_in->k_mult_z_minus_x);

        double final_sum = 0.0;
        for (int i = 0; i < horizon; i++) {
            final_sum += kf_in->x[i];
        }
        kf_out->final_sum = final_sum;

        return (RetCode)ret;
    } else {
        return K_RET_UNK_ERROR;
    }
}

RetCode kf_linux_ioself_profiling(int horizon)
{
    return run_one_profile(horizon);
}

RetCode kf_linux_ioself(int horizon)
{
    return kf_self_test(horizon, 0);
}