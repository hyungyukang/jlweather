module FortranInterface

enterdata = """
module enterdatamod
USE, INTRINSIC :: ISO_C_BINDING

public dataenter

contains

INTEGER (C_INT64_T) FUNCTION dataenter() BIND(C, name="dataenter")
    USE, INTRINSIC :: ISO_C_BINDING

    INTEGER (C_INT64_T) :: res

    res = 0

    dataenter = res

END FUNCTION

end module
"""

end
