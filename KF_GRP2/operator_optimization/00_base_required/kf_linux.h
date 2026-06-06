#ifndef KF_LINUX_H
#define KF_LINUX_H

#include <stdint.h>

#define KF_MAX_DIM 300
#define KF_DEFAULT_STEPS 1

typedef enum RetCode_e {
    K_RET_UNK_ERROR = -1,
    K_RET_OK = 0,
    K_RET_TIMEOUT,
    K_RET_INVALID_PARAM,
    K_RET_BUSY,
} RetCode_t;

typedef struct KfInput_s {
    int dim;
    int steps_as_iterations;
    double *measurements;
} KfInput;

typedef struct KfOutput_s {
    double final_sum;
} KfOutput;

RetCode_t kf_linux_iopointer(int horizon, void *input, void *output);
RetCode_t kf_linux_ioself_profiling(int horizon);
RetCode_t kf_linux_ioself(int horizon);

#endif
