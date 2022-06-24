#include <stdint.h>
#include <stdio.h>

extern "C" int64_t dataenter(int64_t dtype, int64_t * offsets, int64_t * size, double *** state) {
    int64_t res;

    printf("Inside of dataenter: %p\n", state);
    //printf("Inside of dataenter: %f\n", state[0][0]);
    res = 0;
    return res;
}
