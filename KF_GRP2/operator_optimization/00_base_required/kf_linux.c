#include "kf_linux.h"
#include "timestamp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static RetCode_t prepare_kf_input(int horizon, KfInput *input);
static void release_kf_input(KfInput *input);
static RetCode_t check_kf_input(int horizon, const KfInput *input, const KfOutput *output);

static void mat_add(const double *a, const double *b, double *result, size_t dim_n, size_t dim_m);
static void mat_sub(const double *a, const double *b, double *result, size_t dim_n, size_t dim_m);
static void mat_mul(const double *a, const double *b, double *result, size_t n, size_t k_common, size_t m);
static void mat_identity(double *identity, size_t dim);
static void mat_inv_diag(const double *a, double *result, size_t dim);

RetCode_t kf_linux_iopointer(int horizon, void *input, void *output) {
    KfInput *kf_input = (KfInput *)input;
    KfOutput *kf_output = (KfOutput *)output;
    RetCode_t check_ret;
    size_t dim;
    size_t matrix_count;
    double final_sum = 0.0;

    double *x = NULL;
    double *p = NULL;
    double *q = NULL;
    double *r = NULL;
    double *k = NULL;
    double *identity = NULL;
    double *temp_p_plus_r = NULL;
    double *inv_temp_p_plus_r = NULL;
    double *temp_i_minus_k = NULL;
    double *p_old = NULL;
    double *z_minus_x = NULL;
    double *k_mult_z_minus_x = NULL;

    check_ret = check_kf_input(horizon, kf_input, kf_output);
    if (check_ret != K_RET_OK) {
        return check_ret;
    }

    dim = (size_t)kf_input->dim;
    matrix_count = dim * dim;
    kf_output->final_sum = 0.0;

    /*
     * Step 0 baseline:
     * Keep teacher's dense matrix KF logic.
     * Keep workspace allocation inside iopointer().
     */
    x = (double *)calloc(dim, sizeof(double));
    p = (double *)calloc(matrix_count, sizeof(double));
    q = (double *)calloc(matrix_count, sizeof(double));
    r = (double *)calloc(matrix_count, sizeof(double));
    k = (double *)calloc(matrix_count, sizeof(double));
    identity = (double *)calloc(matrix_count, sizeof(double));

    temp_p_plus_r = (double *)calloc(matrix_count, sizeof(double));
    inv_temp_p_plus_r = (double *)calloc(matrix_count, sizeof(double));
    temp_i_minus_k = (double *)calloc(matrix_count, sizeof(double));
    p_old = (double *)calloc(matrix_count, sizeof(double));
    z_minus_x = (double *)calloc(dim, sizeof(double));
    k_mult_z_minus_x = (double *)calloc(dim, sizeof(double));

    if (x == NULL ||
        p == NULL ||
        q == NULL ||
        r == NULL ||
        k == NULL ||
        identity == NULL ||
        temp_p_plus_r == NULL ||
        inv_temp_p_plus_r == NULL ||
        temp_i_minus_k == NULL ||
        p_old == NULL ||
        z_minus_x == NULL ||
        k_mult_z_minus_x == NULL) {
        fprintf(stderr, "KF workspace allocation failed, horizon = %d\n", horizon);

        free(x);
        free(p);
        free(q);
        free(r);
        free(k);
        free(identity);
        free(temp_p_plus_r);
        free(inv_temp_p_plus_r);
        free(temp_i_minus_k);
        free(p_old);
        free(z_minus_x);
        free(k_mult_z_minus_x);

        return K_RET_UNK_ERROR;
    }

    for (size_t i = 0; i < dim; i++) {
        x[i] = 0.0;
        p[i * dim + i] = 1.0;
        q[i * dim + i] = 0.5;
        r[i * dim + i] = 0.5;
    }

    mat_identity(identity, dim);

    for (int step = 0; step < kf_input->steps_as_iterations; step++) {
        const double *z = kf_input->measurements + (size_t)step * dim;

        /*
         * Teacher's simplified KF:
         * P = P + Q
         * K = P * inv(P + R)
         * x = x + K * (z - x)
         * P = (I - K) * P
         */
        mat_add(p, q, p, dim, dim);

        mat_add(p, r, temp_p_plus_r, dim, dim);
        mat_inv_diag(temp_p_plus_r, inv_temp_p_plus_r, dim);
        mat_mul(p, inv_temp_p_plus_r, k, dim, dim, dim);

        for (size_t i = 0; i < dim; i++) {
            z_minus_x[i] = z[i] - x[i];
        }

        memset(k_mult_z_minus_x, 0, dim * sizeof(double));

        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j < dim; j++) {
                k_mult_z_minus_x[i] += k[i * dim + j] * z_minus_x[j];
            }
        }

        for (size_t i = 0; i < dim; i++) {
            x[i] += k_mult_z_minus_x[i];
        }

        memcpy(p_old, p, matrix_count * sizeof(double));
        mat_sub(identity, k, temp_i_minus_k, dim, dim);
        mat_mul(temp_i_minus_k, p_old, p, dim, dim, dim);
    }

    for (size_t i = 0; i < dim; i++) {
        final_sum += x[i];
    }

    kf_output->final_sum = final_sum;

    free(x);
    free(p);
    free(q);
    free(r);
    free(k);
    free(identity);
    free(temp_p_plus_r);
    free(inv_temp_p_plus_r);
    free(temp_i_minus_k);
    free(p_old);
    free(z_minus_x);
    free(k_mult_z_minus_x);

    return K_RET_OK;
}

RetCode_t kf_linux_ioself_profiling(int horizon) {
    KfInput input;
    KfOutput output;
    TimeStamp start;
    TimeStamp end;
    int64_t elapsed_ns;
    RetCode_t ret;

    ret = prepare_kf_input(horizon, &input);
    if (ret != K_RET_OK) {
        return ret;
    }

    output.final_sum = 0.0;

    start = timestamp();
    ret = kf_linux_iopointer(horizon, &input, &output);
    end = timestamp();

    elapsed_ns = timestamp_diff(start, end);

    if (ret == K_RET_OK) {
        printf("%d,%lld,%.6f,ARM Linux,KF,OK,%.6f\n",
               horizon,
               (long long)elapsed_ns,
               (double)elapsed_ns / 1000000.0,
               output.final_sum);
    } else {
        printf("%d,%lld,%.6f,ARM Linux,KF,ERROR,%.6f\n",
               horizon,
               (long long)elapsed_ns,
               (double)elapsed_ns / 1000000.0,
               output.final_sum);
    }

    release_kf_input(&input);
    return ret;
}

RetCode_t kf_linux_ioself(int horizon) {
    KfInput input;
    KfOutput output;
    RetCode_t ret;

    ret = prepare_kf_input(horizon, &input);
    if (ret != K_RET_OK) {
        return ret;
    }

    output.final_sum = 0.0;
    ret = kf_linux_iopointer(horizon, &input, &output);

    release_kf_input(&input);
    return ret;
}

static RetCode_t prepare_kf_input(int horizon, KfInput *input) {
    size_t dim;
    size_t measurement_count;

    if (input == NULL || horizon <= 0 || horizon > KF_MAX_DIM) {
        return K_RET_INVALID_PARAM;
    }

    memset(input, 0, sizeof(*input));

    input->dim = horizon;
    input->steps_as_iterations = KF_DEFAULT_STEPS;

    dim = (size_t)input->dim;
    measurement_count = dim * (size_t)input->steps_as_iterations;

    input->measurements = (double *)calloc(measurement_count, sizeof(double));
    if (input->measurements == NULL) {
        fprintf(stderr, "KF input allocation failed, horizon = %d\n", horizon);
        return K_RET_UNK_ERROR;
    }

    srand(12345);

    for (int step = 0; step < input->steps_as_iterations; step++) {
        double *z = input->measurements + (size_t)step * dim;

        for (size_t i = 0; i < dim; i++) {
            double true_state = (double)(step + 1 + i);
            double measurement_noise = ((rand() / (RAND_MAX / 2.0)) - 1.0) * sqrt(0.5);
            z[i] = true_state + measurement_noise;
        }
    }

    return K_RET_OK;
}

static void release_kf_input(KfInput *input) {
    if (input == NULL) {
        return;
    }

    free(input->measurements);
    memset(input, 0, sizeof(*input));
}

static RetCode_t check_kf_input(int horizon, const KfInput *input, const KfOutput *output) {
    if (horizon <= 0 || horizon > KF_MAX_DIM || input == NULL || output == NULL) {
        return K_RET_INVALID_PARAM;
    }

    if (input->dim != horizon ||
        input->dim <= 0 ||
        input->dim > KF_MAX_DIM ||
        input->steps_as_iterations <= 0 ||
        input->measurements == NULL) {
        return K_RET_INVALID_PARAM;
    }

    return K_RET_OK;
}

static void mat_identity(double *identity, size_t dim) {
    memset(identity, 0, sizeof(double) * dim * dim);

    for (size_t i = 0; i < dim; i++) {
        identity[i * dim + i] = 1.0;
    }
}

static void mat_add(const double *a, const double *b, double *result, size_t dim_n, size_t dim_m) {
    size_t total_count = dim_n * dim_m;

    for (size_t i = 0; i < total_count; i++) {
        result[i] = a[i] + b[i];
    }
}

static void mat_sub(const double *a, const double *b, double *result, size_t dim_n, size_t dim_m) {
    size_t total_count = dim_n * dim_m;

    for (size_t i = 0; i < total_count; i++) {
        result[i] = a[i] - b[i];
    }
}

static void mat_mul(const double *a,
                    const double *b,
                    double *result,
                    size_t n,
                    size_t k_common,
                    size_t m) {
    memset(result, 0, sizeof(double) * n * m);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            double sum_value = 0.0;

            for (size_t l = 0; l < k_common; l++) {
                sum_value += a[i * k_common + l] * b[l * m + j];
            }

            result[i * m + j] = sum_value;
        }
    }
}

static void mat_inv_diag(const double *a, double *result, size_t dim) {
    memset(result, 0, sizeof(double) * dim * dim);

    for (size_t i = 0; i < dim; i++) {
        double diag_value = a[i * dim + i];

        if (diag_value != 0.0) {
            result[i * dim + i] = 1.0 / diag_value;
        }
    }
}
