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
    
    matrix_float* A = new_matrix_float(3, 6);
    matrix_float* b = new_matrix_float(3, 1);
    matrix_float* x = new_matrix_float(3, 1);

    M(x, fill_rand, -10, 10, true);
    M(b, fill_rand, -10, 10, true);
    { 
        matrix_float* xt = new_matrix_float(3, 1);
        M(xt, transp, x);
        M(b, print);
        M(xt, print);
        M(A, mult, b, xt);
        M(xt, free);
    }
    M(A, print);
    M(A, diag, A);
    M(A, print);
    matrix_float* mtzs[] = {A, b, x};
    free_all(mtzs, 3);
    return 0;
}
