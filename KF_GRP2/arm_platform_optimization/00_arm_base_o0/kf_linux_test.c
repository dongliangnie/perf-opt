#include "kf_linux.h"

#include <stdio.h>

int main(void) {
    printf("horizon,time_ns,time_ms,platform,operator,status,final_sum\n");

    for (int horizon = 10; horizon <= KF_MAX_DIM; horizon += 9) {
        RetCode_t ret = kf_linux_ioself_profiling(horizon);

        if (ret != K_RET_OK) {
            fprintf(stderr,
                    "kf_linux_ioself_profiling failed, horizon = %d, ret = %d\n",
                    horizon,
                    ret);
        }
    }

    return 0;
}
