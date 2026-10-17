#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "include/cMatrix.h"

void free_all(cmatrix_float* mtrzs[], int n){
    if(!mtrzs) return;
    for(int i = 0; i < n; i++){
        if(mtrzs[i] != NULL){
            cM(mtrzs[i], free); 
            mtrzs[i] = NULL;
        }
    }
}

int64_t next_combination(int64_t x){
    int64_t c = x & -x;
    int64_t r = x + c;
    return (((x^r) >> 2) / c) | r;
}

int comb(int n, int m){
    if(n < m) return 0;
    float c = 1;
    for(int i = 0; i < m; i++){
        c *= ((float)(n - i)) / (m - i);
    }
    return  roundf(c);
}

void print_bits(long x, int n){
    for(int i = n - 1; i >= 0; i--){
        printf("%lld", (x >> i) & 1ll);
    }
    putchar('\n');
}

#define true 1

float diag_mult(cmatrix_float* A){
    float diag_det = 1;
    for(size_t i = 0; i < A->rows; i++){
        size_t idx = i*A->cols + i;
        diag_det *= A->data[idx];
    }
    return diag_det;
}

int main(){
    size_t m, n;
    int vbn;
    scanf("%d", &vbn);
    scanf(" %zu %zu", &m, &n);
    cmatrix_float* Ab = get_cmatrix_float(m, n); 
    --n;
    cmatrix_float* b = new_cmatrix_float(Ab->rows, 1);
    cM(b, getCol, Ab, n);
    cmatrix_float* v[n];

    for(int i = 0; i < n; i++) {
        v[i] = new_cmatrix_float(Ab->rows, 1);
        cM(v[i], getCol, Ab, i);
    }

    int n_combs = comb(n, m);
    int64_t mask = (1 << m) - 1;

    cmatrix_float* Bb = new_cmatrix_float(m, m + 1);
    cmatrix_float* x = new_cmatrix_float(n, 1);
    int cntbasics = 0;
    for(int i = 0; i < n_combs; i++){
        cM(Bb, setCol, b, m);
        int cnt = 0;
        for(int j = 0; j < n; j++){
            if((mask >> j) & 1)    
                cM(Bb, setCol, v[j], cnt++);
        }

        cM(Bb, rref, Bb);
        if(fabs(diag_mult(Bb)) > 1e-6){
            int top = 0;
            for(int j = 0; j < n; j++){
                x->data[j] = ((mask >> j) & 1) ? Bb->data[top*Bb->cols + m] : 0;
                top += (mask >> j) & 1;
            }
            if((cM(x, atLeast, 0) && vbn) || !vbn){
                cM(x, print);
                cntbasics++;
            }
        }
        mask = next_combination(mask);
    }
    printf("%d\n", cntbasics);
    cM(Ab, free);
    return 0;
}

