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
    cmatrix_float* Ab = get_cmatrix_float(5, 8);
    int n = Ab->cols - 1;
    cmatrix_float* b = get_cmatrix_float(Ab->rows, 1);
    cM(b, getCol, Ab, n);
    cmatrix_float* v[n];

    for(int i = 0; i < n; i++) {
        v[i] = new_cmatrix_float(Ab->rows, 1);
        cM(v[i], getCol, Ab, i);
    }

    cM(Ab, free);
    return 0;
}












