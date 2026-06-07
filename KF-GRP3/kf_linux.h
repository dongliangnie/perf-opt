#ifndef KF_LINUX_H
#define KF_LINUX_H

typedef enum {
    K_RET_OK = 0,
    K_RET_UNK_ERROR = -1,
    K_RET_TIMEOUT,
    K_RET_INVALID_PARAM,
    K_RET_BUSY
} RetCode;

typedef struct {
    int dim;
} KfInput;

typedef struct {
    double *x;
} KfOutput;

RetCode kf_linux_iopointer(int dim, void *input, void *output);
RetCode kf_linux_ioself_profiling(int dim);
RetCode kf_linux_ioself(int dim);

#endif