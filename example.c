#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "include/Matrix.h"

typedef struct{
    float a, b;
}params;

float fx(float x, void* p){
    params* ps = (params*)p;
    return (x)*ps->a + ps->b;
}

void free_all(matrix_float* mtrzs[], int n){
    if(!mtrzs) return;
    for(int i = 0; i < n; i++){
        if(mtrzs[i] != NULL){
            M(mtrzs[i], free); 
            mtrzs[i] = NULL;
        }
    }
}

#define true 1

int main(){
    srand(time(NULL));
    matrix_float* A = new_matrix_float(3, 3);
    matrix_float* A_inv = new_matrix_float(3, 3);
    matrix_float* C = new_matrix_float(3, 3);
    matrix_float* b = new_matrix_float(3, 1);
    matrix_float* x = new_matrix_float(3, 1);
      
    M(A, fill_rand, -10, 10, true);   
    M(A, print);
    M(C, fill_rand, -10, 10, true);
    M(C, print);
    printf("dot(A, C) = %.0f\n", M(A, dot, C));

    M(C, fill, 2);
    M(C, print);

    matrix_float* mtzs[] = {A, C, A_inv, b, x};
    free_all(mtzs, 5);
    return 0;
}
