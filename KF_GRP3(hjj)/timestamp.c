#include "timestamp.h"

TimeStamp timestamp(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts;
}

int64_t timestamp_diff(TimeStamp start, TimeStamp end)
{
    int64_t seconds = end.tv_sec - start.tv_sec;
    int64_t nanoseconds = end.tv_nsec - start.tv_nsec;
    return seconds * 1000000000LL + nanoseconds;
}