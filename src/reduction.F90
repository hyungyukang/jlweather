module reducemod
USE, INTRINSIC :: ISO_C_BINDING

public launch

contains

INTEGER (C_INT64_T) FUNCTION launch(nargs, dtype, offsets, size, state) BIND(C, name="launch")
    USE, INTRINSIC :: ISO_C_BINDING

    INTEGER (C_INT64_T), INTENT(IN) :: nargs
    INTEGER (C_INT64_T), INTENT(IN) :: dtype
    INTEGER (C_INT64_T), INTENT(IN), DIMENSION(3) :: offsets
    INTEGER (C_INT64_T), INTENT(IN), DIMENSION(3) :: size
    REAL    (C_DOUBLE), INTENT(INOUT), DIMENSION(-1:102, -1:52, 1:4) :: state

    INTEGER (C_INT64_T) :: res, i

    res = 0

    print *, "Inside of kernel"

    do i=1, nargs

        print *, "Inside of kernel, arg ", i, ": dtype = ", dtype
        print *, "Inside of kernel, arg ", i, ": offsets = ", offsets
        print *, "Inside of kernel, arg ", i, ": size = ", size
        print *, "Inside of kernel, arg ", i, ": state(1, 0, 1) = ", state(1, 0, 1)

        state(1, 0, 1) = 200.0
    end do

    launch = res

END FUNCTION

end module
