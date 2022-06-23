#include <stdint.h>
#include <stdio.h>

extern "C" int64_t dataenter(double * state[]) {
    int64_t res;

    printf("Inside of dataenter\n");
    res = 0;
    return res;
}
