#include <stdint.h>
#include <stdio.h>

extern "C" int64_t kernel() {
    int64_t res;

    printf("Inside of kernel\n");

    res = 0;
    return res;
}
