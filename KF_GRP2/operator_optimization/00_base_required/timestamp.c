#define _POSIX_C_SOURCE 200809L

#include "timestamp.h"

TimeStamp timestamp(void) {
    TimeStamp now;
    clock_gettime(CLOCK_MONOTONIC, &now.value);
    return now;
}

int64_t timestamp_diff(TimeStamp start, TimeStamp end) {
    int64_t sec_diff = (int64_t)(end.value.tv_sec - start.value.tv_sec);
    int64_t nsec_diff = (int64_t)(end.value.tv_nsec - start.value.tv_nsec);
    return sec_diff * 1000000000LL + nsec_diff;
}
