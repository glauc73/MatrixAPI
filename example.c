#include <stdio.h>
#include "include/Matrix.h"

typedef struct{
    float a, b;
}params;

float fx(float x, void* p){
    params* ps = (params*)p;
    return x*ps->a + ps->b;
}

int main(){
    matrix_float* A = new_matrix_float(2, 2);
    matrix_float* B = new_matrix_float(2, 2);
    matrix_float* C = new_matrix_float(2, 2);
    
    float A_2d[2][2] = {
        {3, -7}, 
        {0, 4}
    };
    matrix_float_set_array_2d(A, A_2d);
    A->pow(&A, A, 2);
    B->inv(&B, A);

    params p = {2, -1};
    
    B->print(B); 
    B->apply(&B, fx, &p);
    B->print(B);

    A->free(&A);
    B->free(&B);
    C->free(&C);
    return 0;
}
