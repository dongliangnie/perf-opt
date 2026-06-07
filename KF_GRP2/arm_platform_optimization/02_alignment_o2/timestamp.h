#ifndef TIMESTAMP_H
#define TIMESTAMP_H

#include <stdint.h>
#include <time.h>

typedef struct TimeStamp_s {
    struct timespec value;
} TimeStamp;

TimeStamp timestamp(void);
int64_t timestamp_diff(TimeStamp start, TimeStamp end);

#endif
