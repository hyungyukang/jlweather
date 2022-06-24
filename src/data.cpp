#include <stdint.h>
#include <stdio.h>

extern "C" int64_t dataenter(int64_t nargs, int64_t dtype, int64_t offsets[3], int64_t size[3], double state[4][54][104]) {
    int64_t res;

    printf("Inside of dataenter: nargs = %lld\n", nargs);

    for (int i=0; i<nargs; i++) {
        printf("Inside of dataenter, arg %d: dtype = %lld\n", i, dtype);
        printf("Inside of dataenter, arg %d: offsets = %lld, %lld, %lld\n", i, offsets[0], offsets[1], offsets[2]);
        printf("Inside of dataenter, arg %d: size = %lld, %lld, %lld\n", i, size[0], size[1], size[2]);
        printf("Inside of dataenter, arg %d: state[0][1][2] = %lf\n", i, state[0][1][2]);

        state[0][1][2] = 200.0;
    }
    res = 0;
    return res;
}
