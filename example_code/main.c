#include "Matrix.h"
#include <math.h>
#include <stdio.h>

typedef struct{
    int x, y;
}parameters;

double comb_linear(double a, void* args){
    if(args == NULL)  //se nao passar nenhum parametro
        return log10(a);
        
    parameters* p = (parameters*) args;
    return a*p->x + a*a*p->y;
}

int main(){
    size_t m, n;

    parameters p = {
        .x = 2,
        .y = 3
    };

    scanf("%zu %zu", &m, &n);
    matrix_double* B = get_matrix_double(m, n);

    B->apply(&B, comb_linear, &p);
    B->print(B);
    
    putchar('\n'); 
    
    B->rref(B);
    B->print(B);
    
    B->free(&B);
    return 0;
}