#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "include/cMatrix.h"

typedef struct{
    float a, b;
}params;

float fx(float x, void* p){
    params* ps = (params*)p;
    return (x)*ps->a + ps->b;
}

void free_all(cmatrix_float* mtrzs[], int n){
    if(!mtrzs) return;
    for(int i = 0; i < n; i++){
        if(mtrzs[i] != NULL){
            cM(mtrzs[i], free); 
            mtrzs[i] = NULL;
        }
    }
}

#define true 1

int main(){
    srand(time(NULL));
    
    cmatrix_float* A = new_cmatrix_float(3, 6);
    cmatrix_float* b = new_cmatrix_float(3, 1);

    float x_2d[3][1] = {1, 2, 3};
    cmatrix_float* x = new_cmatrix_float(3, 1);
    cmatrix_float_set2D(&x, 3, 1, x_2d);
    
    cM(b, fill_rand, -10, 10, true);
    { 
        cmatrix_float* xt = new_cmatrix_float(3, 1);
        cM(xt, transp, x);
        cM(b, print);
        cM(xt, print);
        cM(A, mult, b, xt);
        cM(xt, free);
    }
    cM(A, print);
    cM(A, diag, A);
    cM(A, print);
    cmatrix_float* mtzs[] = {A, b, x};
    free_all(mtzs, 3);
    return 0;
}
