#ifndef GENERICS_BASE_H
#define GENERICS_BASE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>


//nome
#define _CONCAT(X,Y) X##Y
#define CONCAT(X,Y) _CONCAT(X,Y)

#define abs_gen(x) x > 0 ? x : -x

//comparação
#ifndef eps
#define eps 1e-8
#endif

static inline int cmp_double(double a, double b){
    if(fabs(a - b) < eps) return 0;
    return (a < b) ? -1 : 1;
}

static inline int cmp_int64(long long a, long long b){
    if(a == b) return 0;
    return (a < b) ? -1 : 1;
}

static inline int cmp_uint64(unsigned long long a, unsigned long long b){
    if(a == b) return 0;
    return (a < b) ? -1 : 1;
}

#define numcmp(x, y) \
_Generic((x),\
    long long: cmp_int64,\
    int: cmp_int64,\
    unsigned: cmp_uint64,\
    char: cmp_uint64,\
    unsigned long long: cmp_uint64, \
    float: cmp_double, \
    double: cmp_double\
)((x), (y))\

#define equal(x, y) (numcmp((x), (y)) == 0)


//leitura
static inline void input_ll(long long* ptr){
    scanf(" %lld", ptr);
}

static inline void input_int(int* ptr){
    scanf(" %d", ptr);
}

static inline void input_uint(unsigned int* ptr){
    scanf(" %u", ptr);
}

static inline void input_ull(uint64_t *ptr){
    scanf(" %zu", ptr);
}

static inline void input_f(float* ptr){
    scanf(" %f", ptr);
}

static inline void input_dbl(double *ptr){
    scanf(" %lf", ptr);
}

static inline void input_ch(char* ptr){
    scanf(" %c", ptr);
}

#define input(x) \
_Generic((x),\
    long long: input_ll,\
    int: input_int,\
    char: input_ch,\
    unsigned int: input_uint,\
    unsigned long long: input_ull, \
    float: input_f, \
    double: input_dbl\
)(&(x))


//escrita
static inline void print_ll(long long x){
    printf("%lld", x);
}

static inline void print_ull(unsigned long long x){
    printf("%llu", x);
}

static inline void print_f(float x){
    printf("%.2f", x);
}

static inline void print_dbl(double x){
    printf("%.2lf", x);
}

static inline void print_ch(char x){
    printf("%c", x);
}

#define printg(x) \
_Generic((x),\
    long long: print_ll,\
    unsigned int: print_ull,\
    int: print_ll,\
    char: print_ch,\
    unsigned long long: print_ull, \
    float: print_f, \
    double: print_dbl\
)(x)

#endif 