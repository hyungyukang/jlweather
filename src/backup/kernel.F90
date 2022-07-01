module kernelmod
USE, INTRINSIC :: ISO_C_BINDING

public kernel

contains

INTEGER (C_INT64_T) FUNCTION kernel()
    USE, INTRINSIC :: ISO_C_BINDING

    INTEGER (C_INT64_T) :: res

    res = 0

    print *, "Inside of kernel"

    kernel = res

END FUNCTION

end module
