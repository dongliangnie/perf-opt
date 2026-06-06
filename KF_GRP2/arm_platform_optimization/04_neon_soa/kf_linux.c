#include "kf_linux.h"
#include "timestamp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__aarch64__)
#include <arm_neon.h>
#define KF_USE_NEON 1
#else
#define KF_USE_NEON 0
#endif

static RetCode_t run_one_profile(int horizon, KfMode mode, const char *method_name);
static RetCode_t prepare_kf_input(int horizon, KfMode mode, KfInput *input);
static void release_kf_input(KfInput *input);
static RetCode_t check_kf_input(int horizon, const KfInput *input, const KfOutput *output);
static RetCode_t kf_core_aos_scalar(KfInput *input, KfOutput *output);
static RetCode_t kf_core_soa_neon(KfInput *input, KfOutput *output);

RetCode_t kf_linux_iopointer(int horizon, void *input, void *output) {
    KfInput *kf_input = (KfInput *)input;
    KfOutput *kf_output = (KfOutput *)output;
    RetCode_t check_ret = check_kf_input(horizon, kf_input, kf_output);

    if (check_ret != K_RET_OK) {
        return check_ret;
    }

    if (kf_input->mode == KF_MODE_AOS_SCALAR) {
        return kf_core_aos_scalar(kf_input, kf_output);
    }

    return kf_core_soa_neon(kf_input, kf_output);
}

RetCode_t kf_linux_ioself_profiling(int horizon) {
    RetCode_t ret_aos = run_one_profile(horizon, KF_MODE_AOS_SCALAR, "AoS scalar");
    RetCode_t ret_neon = run_one_profile(horizon, KF_MODE_SOA_NEON, "SoA + NEON");

    if (ret_aos != K_RET_OK) {
        return ret_aos;
    }
    return ret_neon;
}

RetCode_t kf_linux_ioself(int horizon) {
    KfInput input;
    KfOutput output;
    RetCode_t ret = prepare_kf_input(horizon, KF_MODE_SOA_NEON, &input);

    if (ret != K_RET_OK) {
        return ret;
    }

    output.final_sum = 0.0;
    ret = kf_linux_iopointer(horizon, &input, &output);
    release_kf_input(&input);
    return ret;
}

static RetCode_t run_one_profile(int horizon, KfMode mode, const char *method_name) {
    KfInput input;
    KfOutput output;
    TimeStamp start;
    TimeStamp end;
    int64_t elapsed_ns;
    RetCode_t ret = prepare_kf_input(horizon, mode, &input);

    if (ret != K_RET_OK) {
        return ret;
    }

    output.final_sum = 0.0;
    start = timestamp();
    ret = kf_linux_iopointer(horizon, &input, &output);
    end = timestamp();
    elapsed_ns = timestamp_diff(start, end);

    if (ret == K_RET_OK) {
        printf("%d,%d,%lld,%.6f,ARM Linux,KF,%s,OK,%.6f\n",
               horizon,
               input.steps_as_iterations,
               (long long)elapsed_ns,
               (double)elapsed_ns / 1000000.0,
               method_name,
               output.final_sum);
    } else {
        printf("%d,%d,%lld,%.6f,ARM Linux,KF,%s,ERROR,%.6f\n",
               horizon,
               input.steps_as_iterations,
               (long long)elapsed_ns,
               (double)elapsed_ns / 1000000.0,
               method_name,
               output.final_sum);
    }

    release_kf_input(&input);
    return ret;
}

static RetCode_t prepare_kf_input(int horizon, KfMode mode, KfInput *input) {
    size_t dim;
    size_t measurement_count;

    if (input == NULL || horizon <= 0 || horizon > KF_MAX_DIM) {
        return K_RET_INVALID_PARAM;
    }

    memset(input, 0, sizeof(*input));
    input->dim = horizon;
    input->steps_as_iterations = KF_DEFAULT_STEPS;
    input->mode = mode;
    dim = (size_t)input->dim;
    measurement_count = dim * (size_t)input->steps_as_iterations;

    input->measurements = (double *)calloc(measurement_count, sizeof(double));
    if (input->measurements == NULL) {
        release_kf_input(input);
        return K_RET_UNK_ERROR;
    }

    if (mode == KF_MODE_AOS_SCALAR) {
        input->aos_states = (KfAosState *)calloc(dim, sizeof(KfAosState));
        if (input->aos_states == NULL) {
            release_kf_input(input);
            return K_RET_UNK_ERROR;
        }

        for (size_t i = 0; i < dim; i++) {
            input->aos_states[i].x = 0.0;
            input->aos_states[i].p = 1.0;
            input->aos_states[i].q = 0.5;
            input->aos_states[i].r = 0.5;
            input->aos_states[i].k = 0.0;
        }
    } else {
        input->x = (double *)calloc(dim, sizeof(double));
        input->p = (double *)calloc(dim, sizeof(double));
        input->q = (double *)calloc(dim, sizeof(double));
        input->r = (double *)calloc(dim, sizeof(double));
        input->k = (double *)calloc(dim, sizeof(double));

        if (input->x == NULL ||
            input->p == NULL ||
            input->q == NULL ||
            input->r == NULL ||
            input->k == NULL) {
            release_kf_input(input);
            return K_RET_UNK_ERROR;
        }

        for (size_t i = 0; i < dim; i++) {
            input->x[i] = 0.0;
            input->p[i] = 1.0;
            input->q[i] = 0.5;
            input->r[i] = 0.5;
            input->k[i] = 0.0;
        }
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
    free(input->aos_states);
    free(input->x);
    free(input->p);
    free(input->q);
    free(input->r);
    free(input->k);
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

    if (input->mode == KF_MODE_AOS_SCALAR) {
        if (input->aos_states == NULL) {
            return K_RET_INVALID_PARAM;
        }
    } else {
        if (input->x == NULL ||
            input->p == NULL ||
            input->q == NULL ||
            input->r == NULL ||
            input->k == NULL) {
            return K_RET_INVALID_PARAM;
        }
    }

    return K_RET_OK;
}

static RetCode_t kf_core_aos_scalar(KfInput *input, KfOutput *output) {
    size_t dim = (size_t)input->dim;
    double final_sum = 0.0;

    for (int step = 0; step < input->steps_as_iterations; step++) {
        const double *z = input->measurements + (size_t)step * dim;

        for (size_t i = 0; i < dim; i++) {
            KfAosState *state = &input->aos_states[i];
            double denominator;
            double innovation;

            state->p = state->p + state->q;
            denominator = state->p + state->r;
            if (denominator == 0.0) {
                return K_RET_INVALID_PARAM;
            }
            state->k = state->p / denominator;
            innovation = z[i] - state->x;
            state->x = state->x + state->k * innovation;
            state->p = (1.0 - state->k) * state->p;
        }
    }

    for (size_t i = 0; i < dim; i++) {
        final_sum += input->aos_states[i].x;
    }

    output->final_sum = final_sum;
    return K_RET_OK;
}

static RetCode_t kf_core_soa_neon(KfInput *input, KfOutput *output) {
    size_t dim = (size_t)input->dim;
    double final_sum = 0.0;

    for (int step = 0; step < input->steps_as_iterations; step++) {
        const double *z = input->measurements + (size_t)step * dim;
        size_t i = 0;

#if KF_USE_NEON
        float64x2_t one_vec = vdupq_n_f64(1.0);

        for (; i + 1 < dim; i += 2) {
            float64x2_t p_vec = vld1q_f64(input->p + i);
            float64x2_t q_vec = vld1q_f64(input->q + i);
            float64x2_t r_vec = vld1q_f64(input->r + i);
            float64x2_t x_vec = vld1q_f64(input->x + i);
            float64x2_t z_vec = vld1q_f64(z + i);
            float64x2_t denominator_vec;
            float64x2_t k_vec;
            float64x2_t innovation_vec;

            p_vec = vaddq_f64(p_vec, q_vec);
            denominator_vec = vaddq_f64(p_vec, r_vec);
            k_vec = vdivq_f64(p_vec, denominator_vec);
            innovation_vec = vsubq_f64(z_vec, x_vec);
            x_vec = vaddq_f64(x_vec, vmulq_f64(k_vec, innovation_vec));
            p_vec = vmulq_f64(vsubq_f64(one_vec, k_vec), p_vec);

            vst1q_f64(input->p + i, p_vec);
            vst1q_f64(input->k + i, k_vec);
            vst1q_f64(input->x + i, x_vec);
        }
#endif

        for (; i < dim; i++) {
            double denominator;
            double innovation;

            input->p[i] = input->p[i] + input->q[i];
            denominator = input->p[i] + input->r[i];
            if (denominator == 0.0) {
                return K_RET_INVALID_PARAM;
            }
            input->k[i] = input->p[i] / denominator;
            innovation = z[i] - input->x[i];
            input->x[i] = input->x[i] + input->k[i] * innovation;
            input->p[i] = (1.0 - input->k[i]) * input->p[i];
        }
    }

    for (size_t i = 0; i < dim; i++) {
        final_sum += input->x[i];
    }

    output->final_sum = final_sum;
    return K_RET_OK;
}
