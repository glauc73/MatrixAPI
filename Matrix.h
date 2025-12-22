#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdio.h>




#ifndef _MATRIX_INT_
    #define _MATRIX_INT_
    #undef GENERICS_NAME
    #undef GENERICS_TYPE
    #define GENERICS_NAME int
    #define GENERICS_TYPE int
    #include "genericMatrix.h"
#endif

#ifndef _MATRIX_CHAR_
    #define _MATRIX_CHAR_
    #undef GENERICS_NAME
    #undef GENERICS_TYPE
    #define GENERICS_NAME char
    #define GENERICS_TYPE char
    #include "genericMatrix.h"
#endif

#ifndef _MATRIX_UINT_
    #define _MATRIX_UINT_
    #undef GENERICS_NAME
    #undef GENERICS_TYPE
    #define GENERICS_NAME uint
    #define GENERICS_TYPE unsigned
    #include "genericMatrix.h"
#endif

#ifndef _MATRIX_FLOAT_
    #define _MATRIX_FLOAT_
    #undef GENERICS_NAME
    #undef GENERICS_TYPE
    #define GENERICS_NAME float
    #define GENERICS_TYPE float    
    #include "genericMatrix.h"
#endif

#ifndef _MATRIX_DOUBLE_
    #define _MATRIX_DOUBLE_   
    #undef GENERICS_NAME
    #undef GENERICS_TYPE
    #define GENERICS_NAME double
    #define GENERICS_TYPE double
    #include "genericMatrix.h"
#endif

#undef MATRIX_TYPE

#endif