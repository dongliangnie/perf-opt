#include "kf_linux.h"
#include "timestamp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static RetCode_t prepare_kf_input(int horizon, KfInput *input);
static void release_kf_input(KfInput *input);
static RetCode_t check_kf_input(int horizon, const KfInput *input, const KfOutput *output);

static RetCode_t create_workspace(size_t dim, KfWorkspace **workspace_out);
static void release_workspace(KfWorkspace *workspace);
static RetCode_t init_workspace(KfWorkspace *workspace, size_t dim);

static void mat_add(const double *a, const double *b, double *result, size_t dim_n, size_t dim_m);
static void mat_sub(const double *a, const double *b, double *result, size_t dim_n, size_t dim_m);
static void mat_mul(const double *a, const double *b, double *result, size_t n, size_t k_common, size_t m);
static void mat_identity(double *identity, size_t dim);
static void mat_inv_diag(const double *a, double *result, size_t dim);

RetCode_t kf_linux_iopointer(int horizon, void *input, void *output) {
    KfInput *kf_input = (KfInput *)input;
    KfOutput *kf_output = (KfOutput *)output;
    KfWorkspace *workspace;
    RetCode_t check_ret;
    size_t dim;
    size_t matrix_count;
    double final_sum = 0.0;

    check_ret = check_kf_input(horizon, kf_input, kf_output);
    if (check_ret != K_RET_OK) {
        return check_ret;
    }

    dim = (size_t)kf_input->dim;
    matrix_count = dim * dim;
    workspace = kf_input->workspace;
    kf_output->final_sum = 0.0;

    for (int step = 0; step < kf_input->steps_as_iterations; step++) {
        const double *z = kf_input->measurements + (size_t)step * dim;

        mat_add(workspace->p, workspace->q, workspace->p, dim, dim);

        mat_add(workspace->p, workspace->r, workspace->temp_p_plus_r, dim, dim);
        mat_inv_diag(workspace->temp_p_plus_r, workspace->inv_temp_p_plus_r, dim);
        mat_mul(workspace->p, workspace->inv_temp_p_plus_r, workspace->k, dim, dim, dim);

        {
            size_t i = 0;
            for (; i + 3 < dim; i += 4) {
                workspace->z_minus_x[i] = z[i] - workspace->x[i];
                workspace->z_minus_x[i + 1] = z[i + 1] - workspace->x[i + 1];
                workspace->z_minus_x[i + 2] = z[i + 2] - workspace->x[i + 2];
                workspace->z_minus_x[i + 3] = z[i + 3] - workspace->x[i + 3];
            }
            for (; i < dim; i++) {
                workspace->z_minus_x[i] = z[i] - workspace->x[i];
            }
        }

        memset(workspace->k_mult_z_minus_x, 0, dim * sizeof(double));

        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j < dim; j++) {
                workspace->k_mult_z_minus_x[i] +=
                    workspace->k[i * dim + j] * workspace->z_minus_x[j];
            }
        }

        {
            size_t i = 0;
            for (; i + 3 < dim; i += 4) {
                workspace->x[i] += workspace->k_mult_z_minus_x[i];
                workspace->x[i + 1] += workspace->k_mult_z_minus_x[i + 1];
                workspace->x[i + 2] += workspace->k_mult_z_minus_x[i + 2];
                workspace->x[i + 3] += workspace->k_mult_z_minus_x[i + 3];
            }
            for (; i < dim; i++) {
                workspace->x[i] += workspace->k_mult_z_minus_x[i];
            }
        }

        memcpy(workspace->p_old, workspace->p, matrix_count * sizeof(double));
        mat_sub(workspace->identity, workspace->k, workspace->temp_i_minus_k, dim, dim);
        mat_mul(workspace->temp_i_minus_k, workspace->p_old, workspace->p, dim, dim, dim);
    }

    {
        size_t i = 0;
        for (; i + 3 < dim; i += 4) {
            final_sum += workspace->x[i];
            final_sum += workspace->x[i + 1];
            final_sum += workspace->x[i + 2];
            final_sum += workspace->x[i + 3];
        }
        for (; i < dim; i++) {
            final_sum += workspace->x[i];
        }
    }

    kf_output->final_sum = final_sum;
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
    RetCode_t ret;

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

    ret = create_workspace(dim, &input->workspace);
    if (ret != K_RET_OK) {
        release_kf_input(input);
        return ret;
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
    release_workspace(input->workspace);

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
        input->measurements == NULL ||
        input->workspace == NULL) {
        return K_RET_INVALID_PARAM;
    }

    if (input->workspace->x == NULL ||
        input->workspace->p == NULL ||
        input->workspace->q == NULL ||
        input->workspace->r == NULL ||
        input->workspace->k == NULL ||
        input->workspace->identity == NULL ||
        input->workspace->temp_p_plus_r == NULL ||
        input->workspace->inv_temp_p_plus_r == NULL ||
        input->workspace->temp_i_minus_k == NULL ||
        input->workspace->p_old == NULL ||
        input->workspace->z_minus_x == NULL ||
        input->workspace->k_mult_z_minus_x == NULL) {
        return K_RET_INVALID_PARAM;
    }

    return K_RET_OK;
}

static RetCode_t create_workspace(size_t dim, KfWorkspace **workspace_out) {
    KfWorkspace *workspace;
    size_t matrix_count;

    if (workspace_out == NULL || dim == 0 || dim > KF_MAX_DIM) {
        return K_RET_INVALID_PARAM;
    }

    *workspace_out = NULL;
    matrix_count = dim * dim;

    workspace = (KfWorkspace *)calloc(1, sizeof(KfWorkspace));
    if (workspace == NULL) {
        return K_RET_UNK_ERROR;
    }

    workspace->x = (double *)calloc(dim, sizeof(double));
    workspace->p = (double *)calloc(matrix_count, sizeof(double));
    workspace->q = (double *)calloc(matrix_count, sizeof(double));
    workspace->r = (double *)calloc(matrix_count, sizeof(double));
    workspace->k = (double *)calloc(matrix_count, sizeof(double));
    workspace->identity = (double *)calloc(matrix_count, sizeof(double));

    workspace->temp_p_plus_r = (double *)calloc(matrix_count, sizeof(double));
    workspace->inv_temp_p_plus_r = (double *)calloc(matrix_count, sizeof(double));
    workspace->temp_i_minus_k = (double *)calloc(matrix_count, sizeof(double));
    workspace->p_old = (double *)calloc(matrix_count, sizeof(double));
    workspace->z_minus_x = (double *)calloc(dim, sizeof(double));
    workspace->k_mult_z_minus_x = (double *)calloc(dim, sizeof(double));

    if (workspace->x == NULL ||
        workspace->p == NULL ||
        workspace->q == NULL ||
        workspace->r == NULL ||
        workspace->k == NULL ||
        workspace->identity == NULL ||
        workspace->temp_p_plus_r == NULL ||
        workspace->inv_temp_p_plus_r == NULL ||
        workspace->temp_i_minus_k == NULL ||
        workspace->p_old == NULL ||
        workspace->z_minus_x == NULL ||
        workspace->k_mult_z_minus_x == NULL) {
        release_workspace(workspace);
        return K_RET_UNK_ERROR;
    }

    if (init_workspace(workspace, dim) != K_RET_OK) {
        release_workspace(workspace);
        return K_RET_UNK_ERROR;
    }

    *workspace_out = workspace;
    return K_RET_OK;
}

static void release_workspace(KfWorkspace *workspace) {
    if (workspace == NULL) {
        return;
    }

    free(workspace->x);
    free(workspace->p);
    free(workspace->q);
    free(workspace->r);
    free(workspace->k);
    free(workspace->identity);

    free(workspace->temp_p_plus_r);
    free(workspace->inv_temp_p_plus_r);
    free(workspace->temp_i_minus_k);
    free(workspace->p_old);
    free(workspace->z_minus_x);
    free(workspace->k_mult_z_minus_x);

    free(workspace);
}

static RetCode_t init_workspace(KfWorkspace *workspace, size_t dim) {
    if (workspace == NULL || dim == 0 || dim > KF_MAX_DIM) {
        return K_RET_INVALID_PARAM;
    }

    for (size_t i = 0; i < dim; i++) {
        workspace->x[i] = 0.0;
        workspace->p[i * dim + i] = 1.0;
        workspace->q[i * dim + i] = 0.5;
        workspace->r[i * dim + i] = 0.5;
    }

    mat_identity(workspace->identity, dim);
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
    size_t i = 0;

    for (; i + 3 < total_count; i += 4) {
        result[i] = a[i] + b[i];
        result[i + 1] = a[i + 1] + b[i + 1];
        result[i + 2] = a[i + 2] + b[i + 2];
        result[i + 3] = a[i + 3] + b[i + 3];
    }

    for (; i < total_count; i++) {
        result[i] = a[i] + b[i];
    }
}

static void mat_sub(const double *a, const double *b, double *result, size_t dim_n, size_t dim_m) {
    size_t total_count = dim_n * dim_m;
    size_t i = 0;

    for (; i + 3 < total_count; i += 4) {
        result[i] = a[i] - b[i];
        result[i + 1] = a[i + 1] - b[i + 1];
        result[i + 2] = a[i + 2] - b[i + 2];
        result[i + 3] = a[i + 3] - b[i + 3];
    }

    for (; i < total_count; i++) {
        result[i] = a[i] - b[i];
    }
}

static void mat_mul(const double *a,
                    const double *b,
                    double *result,
                    size_t n,
                    size_t k_common,
                    size_t m) {
    /*
     * Matrix multiplication loop-order optimization.
     * Change i-j-l to i-l-j so b_row[j] and result_row[j]
     * are accessed continuously in the innermost loop.
     */
    memset(result, 0, sizeof(double) * n * m);

    for (size_t i = 0; i < n; i++) {
        double *result_row = result + i * m;

        for (size_t l = 0; l < k_common; l++) {
            double a_value = a[i * k_common + l];
            const double *b_row = b + l * m;

            for (size_t j = 0; j < m; j++) {
                result_row[j] += a_value * b_row[j];
            }
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
