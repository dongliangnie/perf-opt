#include "timestamp.h"
#include <time.h>

uint64_t get_timestamp(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
}

int64_t timestamp_diff(uint64_t start, uint64_t end) {
    return (int64_t)(end - start);
}