#ifndef TIMESTAMP_H
#define TIMESTAMP_H

#include <time.h>
#include <stdint.h>

typedef struct timespec TimeStamp;

TimeStamp timestamp(void);
int64_t timestamp_diff(TimeStamp start, TimeStamp end);

#endif