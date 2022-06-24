module datamod
USE, INTRINSIC :: ISO_C_BINDING

public dataenter

contains

INTEGER (C_INT64_T) FUNCTION dataenter(nargs, dtype, offsets, size, state) BIND(C, name="dataenter")
    USE, INTRINSIC :: ISO_C_BINDING

    INTEGER (C_INT64_T), INTENT(IN) :: nargs
    INTEGER (C_INT64_T), INTENT(IN) :: dtype
    INTEGER (C_INT64_T), INTENT(IN), DIMENSION(3) :: offsets
    INTEGER (C_INT64_T), INTENT(IN), DIMENSION(3) :: size
    REAL    (C_DOUBLE), INTENT(INOUT), DIMENSION(-1:102, -1:52, 1:4) :: state

    INTEGER (C_INT64_T) :: res, i

    res = 0

    print *, "Inside of dataenter: nargs = ", nargs

    do i=1, nargs

        print *, "Inside of dataenter, arg ", i, ": dtype = ", dtype
        print *, "Inside of dataenter, arg ", i, ": offsets = ", offsets
        print *, "Inside of dataenter, arg ", i, ": size = ", size
        print *, "Inside of dataenter, arg ", i, ": state(1, 0, 1) = ", state(1, 0, 1)

        state(1, 0, 1) = 200.0
    end do

    dataenter = res

END FUNCTION

end module
