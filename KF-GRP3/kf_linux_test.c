#include <stdio.h>
#include "kf_linux.h"

int main(void) {
    // 维度从 10 到 300，步长 9，覆盖不同计算规模
    for (int dim = 10; dim <= 300; dim += 9) {
        kf_linux_ioself_profiling(dim);
    }
    return 0;
}