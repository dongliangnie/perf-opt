#ifndef KF_LINUX_H
#define KF_LINUX_H

#include <stddef.h>

#define KF_MAX_DIM 300

typedef enum {
    K_RET_OK = 0,
    K_RET_INVALID_PARAM = 1,
    K_RET_UNK_ERROR = 2
} RetCode;

typedef enum {
    KF_MODE_AOS_SCALAR = 0,
    KF_MODE_SOA_NEON = 1
} KfMode;

typedef struct {
    double x;
    double p;
    double q;
    double r;
    double k;
} KfAosState;

typedef struct {
    int dim;
    int steps_as_iterations;
    KfMode mode;
    double *measurements;
    KfAosState *aos_states;
    double *x;
    double *p;
    double *q;
    double *r;
    double *k;
    double *temp_p_plus_r;
    double *inv_temp_p_plus_r;
    double *temp_i_minus_k;
    double *p_old;
    double *z_minus_x;
    double *k_mult_z_minus_x;
    double *i_mat;
} KfInput;

typedef struct {
    double final_sum;
} KfOutput;

RetCode kf_linux_iopointer(int horizon, void *input, void *output);
RetCode kf_linux_ioself_profiling(int horizon);
RetCode kf_linux_ioself(int horizon);

#endif