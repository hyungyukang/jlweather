module exitdatamod
USE, INTRINSIC :: ISO_C_BINDING

public dataexit

contains

INTEGER (C_INT64_T) FUNCTION dataexit() BIND(C, name="dataexit")
    USE, INTRINSIC :: ISO_C_BINDING

    INTEGER (C_INT64_T) :: res

    res = 0

    dataexit = res

END FUNCTION

end module
