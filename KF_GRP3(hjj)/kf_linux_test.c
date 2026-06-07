#include "kf_linux.h"
#include <stdio.h>

int main(void)
{
    for (int horizon = 10; horizon <= 300; horizon += 9) {
        kf_linux_ioself_profiling(horizon);
    }
    return 0;
}