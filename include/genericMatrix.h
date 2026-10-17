#include "genericsBase.h"

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#ifndef GENERICS_TYPE
#define GENERICS_TYPE int
#endif

#ifndef GENERICS_NAME
#define GENERICS_NAME int
#endif

#ifndef tostr
#define tostr(x) #x
#endif

#ifndef eps
#define eps 1e-8
#endif

#define MATRIX_TYPE CONCAT(cmatrix_, GENERICS_NAME)

#ifndef CMATRIX_ERROR_DEFINED
#define CMATRIX_ERROR_DEFINED
typedef enum {
    CMATRIX_OK = 0,           // operação realizada com sucesso
    CMATRIX_NULL_POINTER,     // ponteiro nulo passado como argumento
    CMATRIX_ALLOC_ERROR,      // falha de alocação de memória
    CMATRIX_DIMENSION_MISMATCH, // dimensões incompatíveis entre as matrizes
    CMATRIX_NOT_SQUARE,       // operação exige matriz quadrada
    CMATRIX_SINGULAR_MATRIX,  // matriz singular (não possui inversa)
    CMATRIX_INCOMPATIBLE_RESHAPE, // novo shape tem número de elementos diferente
    CMATRIX_INVALID_DIAG_INPUT,   // diag exige matriz quadrada ou vetor linha/coluna
    CMATRIX_INDEX_OUT_OF_RANGE,  // índice de linha ou coluna fora do intervalo válido
} cmatrix_error;
#endif

typedef struct CONCAT(VTable_, GENERICS_NAME) CONCAT(VTable_, GENERICS_NAME);

#define cM(obj, method, ...) (obj)->vptr->method(&(obj), ##__VA_ARGS__)

typedef struct MATRIX_TYPE{
    GENERICS_TYPE* data;
    GENERICS_TYPE** at;
    size_t rows, cols;
    CONCAT(VTable_, GENERICS_NAME)* vptr;

    #if !defined(__INTELLISENSE__) && !defined(__clangd__) //Atributos privados
    GENERICS_TYPE det_cache; //atributo que guarda o determinante calculado para evitar recalculo
    int8_t was_modified; //flag alterada cada vez que um metodo modifica a matriz para recalcular o det
    #endif
}MATRIX_TYPE;

typedef struct CONCAT(VTable_, GENERICS_NAME){ 
    cmatrix_error (*mult)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*, const struct MATRIX_TYPE*);
    cmatrix_error (*hadmult)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*, const struct MATRIX_TYPE*);
    cmatrix_error (*add)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*, const struct MATRIX_TYPE*);
    cmatrix_error (*sub)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*, const struct MATRIX_TYPE*);
    cmatrix_error (*pow)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*, uint16_t);
    cmatrix_error (*scale)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*, GENERICS_TYPE);
    cmatrix_error (*cpy)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*);
    cmatrix_error (*free)(struct MATRIX_TYPE**);
    cmatrix_error (*ref)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*);
    cmatrix_error (*identity)(struct MATRIX_TYPE**);
    cmatrix_error (*print)(struct MATRIX_TYPE**);
    cmatrix_error (*hstack)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*, const struct MATRIX_TYPE*);
    cmatrix_error (*apply)(struct MATRIX_TYPE**, GENERICS_TYPE (*)(GENERICS_TYPE, void*), void*);
    cmatrix_error (*transp)(struct MATRIX_TYPE**, struct MATRIX_TYPE*);
    cmatrix_error (*diag)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*); 
    cmatrix_error (*fill)(struct MATRIX_TYPE**, GENERICS_TYPE);
    cmatrix_error (*reshape)(struct MATRIX_TYPE**, size_t, size_t);
    cmatrix_error (*swapRows)(struct MATRIX_TYPE**, size_t, size_t);
    cmatrix_error (*swapCols)(struct MATRIX_TYPE**, size_t, size_t);
    cmatrix_error (*slice)(struct MATRIX_TYPE**, struct MATRIX_TYPE*, size_t, size_t, size_t, size_t);
    cmatrix_error (*getCol)(struct MATRIX_TYPE**, struct MATRIX_TYPE*, size_t);
    cmatrix_error (*getRow)(struct MATRIX_TYPE**, struct MATRIX_TYPE*, size_t);
    cmatrix_error (*setCol)(struct MATRIX_TYPE**, struct MATRIX_TYPE*, size_t);
    cmatrix_error (*setRow)(struct MATRIX_TYPE**, struct MATRIX_TYPE*, size_t);
    GENERICS_TYPE (*det)(struct MATRIX_TYPE**);
    int (*atLeast)(struct MATRIX_TYPE**, GENERICS_TYPE);
   
    
    #if !defined(_MATRIX_DOUBLE_) && !defined(_MATRIX_FLOAT_)
    cmatrix_error (*fill_rand)(struct MATRIX_TYPE**, GENERICS_TYPE, GENERICS_TYPE);
    #endif

    #if defined(_MATRIX_DOUBLE_) || defined(_MATRIX_FLOAT_) //metodos que só fazem sentido para double
    GENERICS_TYPE (*dot)(struct MATRIX_TYPE**, const MATRIX_TYPE*);
    cmatrix_error (*rref)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*);
    cmatrix_error (*inv)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*);
    cmatrix_error (*fill_rand)(struct MATRIX_TYPE**, GENERICS_TYPE, GENERICS_TYPE, int);
    #endif
}CONCAT(VTable, GENERICS_NAME);

#define NEW_MATRIX_TYPE CONCAT(new_, MATRIX_TYPE)
static inline MATRIX_TYPE* NEW_MATRIX_TYPE(size_t rows, size_t cols);

static inline cmatrix_error CONCAT(swap_rows_, MATRIX_TYPE)(MATRIX_TYPE** mat, size_t r1, size_t r2){
    if(!mat || !*mat) return CMATRIX_NULL_POINTER;
    if(r1 >= (*mat)->rows || r2 >= (*mat)->rows) return CMATRIX_INDEX_OUT_OF_RANGE;
    if(r1 == r2) return CMATRIX_OK;
    GENERICS_TYPE* a = (*mat)->at[r1];
    GENERICS_TYPE* b = (*mat)->at[r2];
    for(size_t j = 0; j < (*mat)->cols; j++){
        GENERICS_TYPE t = a[j]; a[j] = b[j]; b[j] = t;
    }
    (*mat)->was_modified = 1;
    return CMATRIX_OK;
}

#if !defined(__INTELLISENSE__) && !defined(__clangd__)

static inline cmatrix_error CONCAT(delete_, MATRIX_TYPE)(MATRIX_TYPE** self){ //destrói corretamente a matriz
    if(self == NULL || *self == NULL) return CMATRIX_NULL_POINTER;
    free((*self)->data);
    free((*self)->at);
    free((*self));
    *self = NULL;
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(copy_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* B){
    if(!self || !*self || !B) return CMATRIX_NULL_POINTER;

    if((*self)->rows != B->rows || (*self)->cols != B->cols){   
        cM(*self, free);
        *self = NEW_MATRIX_TYPE(B->rows, B->cols);
        if(!*self) return CMATRIX_ALLOC_ERROR;
    } 
    memcpy((*self)->data, B->data, B->rows * B->cols * sizeof(GENERICS_TYPE));

    (*self)->was_modified = B->was_modified;

    if(!(*self)->was_modified){
        (*self)->det_cache = B->det_cache;
    }
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(reshape_, MATRIX_TYPE)(MATRIX_TYPE** self, size_t m, size_t n){
    if(!self || !*self) return CMATRIX_NULL_POINTER;

    if((*self)->rows*(*self)->cols != m*n){
        return CMATRIX_INCOMPATIBLE_RESHAPE;
    }
    if((*self)->rows != m){
        GENERICS_TYPE** new_at = (GENERICS_TYPE**) realloc((*self)->at, m*sizeof(GENERICS_TYPE*));
        if(!new_at) return CMATRIX_ALLOC_ERROR;
        (*self)->at = new_at;
    }
    (*self)->rows = m;
    (*self)->cols = n;

    for(int i = 0; i < m; i++){
        (*self)->at[i] = (*self)->data + (i*n);
    }

    (*self)->was_modified = 1;
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(fill_, MATRIX_TYPE)(MATRIX_TYPE** self, GENERICS_TYPE value){
    if(!self || !*self) return CMATRIX_NULL_POINTER;

    size_t n = (*self)->rows * (*self)->cols;
    for(size_t i = 0; i < n; i++){
        (*self)->data[i] = value;
    }
    (*self)->was_modified = 1;
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(add_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A, const MATRIX_TYPE* B){
    if(!self || !*self || !A || !B) return CMATRIX_NULL_POINTER;
    if(A->rows != B->rows || A->cols != B->cols) return CMATRIX_DIMENSION_MISMATCH;
    if((*self)->rows != A->rows || (*self)->cols != A->cols){
        cM(*self, free);
        *self = NEW_MATRIX_TYPE(A->rows, A->cols);
        if(!*self) return CMATRIX_ALLOC_ERROR;
    }
    //soma elemento a elemento e atribui em self
    int n = A->rows * A->cols;
    for(int k = 0; k < n; k++)
        (*self)->data[k] = A->data[k] + B->data[k];
    
    (*self)->was_modified = 1;
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(sub_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A, const MATRIX_TYPE* B){
    if(!self || !*self || !A || !B) return CMATRIX_NULL_POINTER;
    if(A->rows != B->rows || A->cols != B->cols) return CMATRIX_DIMENSION_MISMATCH;

    if((*self)->rows != A->rows || (*self)->cols != A->cols){
        cM(*self, free);
        *self = NEW_MATRIX_TYPE(A->rows, A->cols);
        if(!*self) return CMATRIX_ALLOC_ERROR;
    }

    int n = A->rows * A->cols;
    for(int k = 0; k < n; k++)
        (*self)->data[k] = A->data[k] - B->data[k];
   
    (*self)->was_modified = 1;
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(scale_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A, GENERICS_TYPE alpha){
    if(!self || !*self || !A) return CMATRIX_NULL_POINTER;

    if((*self)->rows != A->rows || (*self)->cols != A->cols){
        cM(*self, free);
        *self = NEW_MATRIX_TYPE(A->rows, A->cols);
        if(!*self) return CMATRIX_ALLOC_ERROR;
    }
     
    int n = A->rows * A->cols;
    for(int k = 0; k < n; k++)
        (*self)->data[k] = A->data[k]*alpha;
    
    (*self)->was_modified = 1;
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(mul_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A, const MATRIX_TYPE* B){
    if(!self || !*self || !A || !B) return CMATRIX_NULL_POINTER;
    if(A->cols != B->rows) return CMATRIX_DIMENSION_MISMATCH;
    
    MATRIX_TYPE* TMP = NEW_MATRIX_TYPE(A->rows, B->cols);
    if(!TMP) return CMATRIX_ALLOC_ERROR;

    //aplica o algoritmo de multiplicação de matrizes;
    for(size_t i = 0; i < TMP->rows; i++){
        for(size_t k = 0; k < A->cols; k++){
            GENERICS_TYPE aik = A->at[i][k];
            for(size_t j = 0; j < TMP->cols; j++){
                TMP->at[i][j] += aik*B->at[k][j];
            }
        }
    }

    //usa temp para permitir in-place
    cM(*self, cpy, TMP);
    cM(TMP, free); 

    (*self)->was_modified = 1;
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(identity_, MATRIX_TYPE)(MATRIX_TYPE** self){
    if(!self || !*self) return CMATRIX_NULL_POINTER;
    if((*self)->rows != (*self)->cols) return CMATRIX_NOT_SQUARE;
    for(size_t i = 0; i < (*self)->rows; i++){
        for(size_t j = 0; j < (*self)->cols; j++){
            (*self)->at[i][j] = (i == j);
        }
    }
    (*self)->det_cache = 1;
    (*self)->was_modified = 0;
    return CMATRIX_OK;
}

static inline int CONCAT(at_least_, MATRIX_TYPE)(MATRIX_TYPE** self, GENERICS_TYPE value){
    size_t n = (*self)->cols * (*self)->rows;
    for(int i = 0; i < n; i++){
        if((*self)->data[i] < value) return 0;
    }
    return 1;
}

static inline cmatrix_error CONCAT(hadamard_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A, const MATRIX_TYPE* B){
    if(!self || !*self || !A || !B) return CMATRIX_NULL_POINTER;
    if(A->cols != B->cols || A->rows != B->rows) return CMATRIX_DIMENSION_MISMATCH;
    if((*self)->cols != A->cols || (*self)->rows != A->rows){
        cM(*self, free);
        *self = NEW_MATRIX_TYPE(A->rows, A->cols);
        if(!*self) return CMATRIX_ALLOC_ERROR;
    }
    
    int r = A->rows, c = A->cols;
    
    int n = A->rows * A->cols;
    for(int k = 0; k < n; k++)
        (*self)->data[k] = A->data[k] * B->data[k];
   
    (*self)->was_modified = 1;
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(transp_, MATRIX_TYPE)(MATRIX_TYPE** self, MATRIX_TYPE* A){
    if(!self || !*self || !A) return CMATRIX_NULL_POINTER;
    MATRIX_TYPE* tmp = NEW_MATRIX_TYPE(A->cols, A->rows);
    if(!tmp) return CMATRIX_ALLOC_ERROR;

    for(int i = 0; i < A->rows; i++){
        for(int j = 0; j < A->cols; j++){
            tmp->at[j][i] = A->at[i][j];
        }
    }
    cM(*self, cpy, tmp);
    (*self)->was_modified = A->was_modified;
    cM(tmp, free);
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(swap_cols_, MATRIX_TYPE)(MATRIX_TYPE** self, size_t c1, size_t c2){
    if(!self || !*self) return CMATRIX_NULL_POINTER;
    if(c1 >= (*self)->cols || c2 >= (*self)->cols) return CMATRIX_INDEX_OUT_OF_RANGE;
    if(c1 == c2) return CMATRIX_OK;

    GENERICS_TYPE* pc1 = (*self)->data + c1*(*self)->rows;
    GENERICS_TYPE* pc2 = (*self)->data + c2*(*self)->rows;
    for(size_t i = 0; i < (*self)->rows; i++){
        GENERICS_TYPE tmp = pc1[i]; pc1[i] = pc2[i]; pc2[i] = tmp;
    }
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(slice_, MATRIX_TYPE)(MATRIX_TYPE** self, MATRIX_TYPE* A, size_t i0, size_t i1, size_t j0, size_t j1){
    if(!A) return CMATRIX_NULL_POINTER;
    size_t m = A->rows;   size_t n = A->cols;
    
    if(i1 < i0) {size_t tmp = i0, i0 = i1, i1 = tmp;}
    if(j1 < j0) {size_t tmp = j0, j0 = j1, j1 = tmp;}
    
    if(i1 >= m || j1 >= n) return CMATRIX_INDEX_OUT_OF_RANGE;
    if(!self || !*self) return CMATRIX_NULL_POINTER;
    
    size_t slice_m = i1 - i0 + 1;
    size_t slice_n = j1 - j0 + 1;

    if((*self)->rows != slice_m || (*self)->cols != slice_n){
        cM(*self, free);
        *self = NEW_MATRIX_TYPE(slice_m, slice_n);
    }
    
    for(size_t i = 0; i < slice_m; i++){
        for(size_t j = 0; j < slice_n; j++){
            (*self)->data[i*slice_n + j] = A->data[(i + i0)*n + j + j0];
        }
    }
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(get_col_, MATRIX_TYPE)(MATRIX_TYPE** self, MATRIX_TYPE* A, size_t col){
    if(!self || !*self || !A) return CMATRIX_NULL_POINTER;
    if(col >= A->cols) return CMATRIX_INDEX_OUT_OF_RANGE;

    cM(*self, slice, A, 0, A->rows - 1, col, col);
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(get_row_, MATRIX_TYPE)(MATRIX_TYPE** self, MATRIX_TYPE* A, size_t row){
    if(!self || !*self || !A) return CMATRIX_NULL_POINTER;
    if(row >= A->rows) return CMATRIX_INDEX_OUT_OF_RANGE;


    cM(*self, slice, A, row, row, 0, A->cols - 1);
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(set_col_, MATRIX_TYPE)(MATRIX_TYPE** self, MATRIX_TYPE* A, size_t col){
    if(!self || !*self || !A) return CMATRIX_NULL_POINTER;
    if(col >= (*self)->cols) return CMATRIX_INDEX_OUT_OF_RANGE;
    if(A->cols > 1) return CMATRIX_DIMENSION_MISMATCH;    
    if((*self)->rows < A->rows) return CMATRIX_DIMENSION_MISMATCH;
    
    for(size_t i = 0; i < A->rows; i++){
        int idx_self = i* (*self)->cols + col;
        (*self)->data[idx_self] = A->data[i];
    }
    for(size_t i = A->rows; i < (*self)->rows; i++){
        int idx = i * (*self)->cols + col;
        (*self)->data[idx] = 0;
    }

    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(set_row_, MATRIX_TYPE)(MATRIX_TYPE** self, MATRIX_TYPE* A, size_t row){
    if(!self || !*self || !A) return CMATRIX_NULL_POINTER;
    if(row >= (*self)->rows) return CMATRIX_INDEX_OUT_OF_RANGE;
    if(A->rows > 1) return CMATRIX_DIMENSION_MISMATCH;    
    if((*self)->cols < A->cols) return CMATRIX_DIMENSION_MISMATCH;
    
    for(size_t j = 0; j < A->cols; j++){
        int idx_self = row* (*self)->cols + j;
        (*self)->data[idx_self] = A->data[j];
    }
    for(size_t j = A->cols; j < (*self)->cols; j++){
        int idx = row * (*self)->cols + j;
        (*self)->data[idx] = 0;
    }

    return CMATRIX_OK;
}


static inline cmatrix_error CONCAT(diag_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A){
    if(!self || !*self || !A) return CMATRIX_NULL_POINTER;

    int is_col = 0, is_row = 0, dim;
    if(A->rows == A->cols){
        is_row = 1, is_col = 1;
        dim = A->rows;
    }else if(A->rows == 1){
        is_row = 1;
        dim = A->cols;
    }else if(A->cols == 1){
        is_col = 1;
        dim = A->rows;
    }else{
        return CMATRIX_INVALID_DIAG_INPUT;
    }
     
    if((*self)->rows != dim || (*self)->cols != dim){
        cM(*self, free);
        *self = NEW_MATRIX_TYPE(dim, dim);
        if(!*self) return CMATRIX_ALLOC_ERROR;
    }
    MATRIX_TYPE* TMP = NEW_MATRIX_TYPE(dim , dim);
    if(!TMP) return CMATRIX_ALLOC_ERROR;
    for(int i = 0; i < dim; i++){
        TMP->at[i][i] = A->at[is_col*i][is_row*i];
    }

    cM(*self, cpy, TMP);
    cM(TMP, free);
    return CMATRIX_OK;
}

static inline cmatrix_error CONCAT(pow_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A, uint16_t n){
    if(!self || !*self || !A) return CMATRIX_NULL_POINTER;
    if(A->cols != A->rows) return CMATRIX_NOT_SQUARE;
    if(n == 0){
        return cM(*self, identity);
    }
    
    if((*self)->cols != A->cols || (*self)->rows != A->rows){
        cM(*self, free);
        *self = NEW_MATRIX_TYPE(A->rows, A->cols);
        if(!*self) return CMATRIX_ALLOC_ERROR;
    }
    n--;
    MATRIX_TYPE* base = NEW_MATRIX_TYPE(A->rows, A->cols);  //base para permitir in-place
    if(!base) return CMATRIX_ALLOC_ERROR;
    cM(base, cpy, A);

    //exponenciação rapida iterativa em O(n³logn)
    while(n > 0){
        if(n & 1){
            cM(*self, mult, *self, base);
        }
        cM(base, mult, base, base);
        n >>= 1;
    }
    cM(base, free);
    (*self)->was_modified = 1;
    return CMATRIX_OK;
}

static inline void CONCAT(add_rows_, MATRIX_TYPE)(MATRIX_TYPE* mat, size_t r1, size_t r2, double k){
    if(!mat){
        printf("ponteiro Nulo\n");
    }
    for(size_t j = 0; j < mat->cols ; j++)
        mat->at[r1][j] += (GENERICS_TYPE)(double)(mat->at[r2][j] * k);
}

static inline void CONCAT(mul_row_, MATRIX_TYPE)(MATRIX_TYPE* mat, size_t row, GENERICS_TYPE k){
    //multiplica uma coluna por um escalar
    if(!mat){
        printf("ponteiro Nulo\n");
        return;
    }
    for(size_t j = 0; j < mat->cols ; j++)
        mat->at[row][j] *= k;
    mat->was_modified = 1;
}

//verifica se uma submatriz de mat delimitada por i0, n_rows, j0 e jmax é nula
static inline int CONCAT(sub_iszero_, MATRIX_TYPE)(MATRIX_TYPE* mat, size_t i0, size_t j0, size_t jmax){
    if(j0 > jmax || j0 >= mat->cols) return 1;
    if(jmax >= mat->cols) jmax = mat->cols - 1;
    for(size_t i = i0; i < mat->rows; i++){
        for(size_t j = j0; j <= jmax; j++){
            if(!equal(mat->at[i][j], 0)) return 0;
        }
    }
    return 1;
}

//reduz a matriz a sua forma Escalonada Reduzida por linhas
static inline cmatrix_error CONCAT(REF_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* mat){
    if(!self || !*self || !mat) return CMATRIX_NULL_POINTER;
    size_t i, j, i0 = 0;
    cM(*self, cpy, mat);

    for(j = 0; j < (*self)->cols; j++){
        i = i0; //inicializa a primeira linha para buscar
        while(i < (*self)->rows && equal((*self)->at[i][j], 0)) i++; //procura o primeiro elemento nao nulo na linha
        if(i == (*self)->rows) //se nao encontou, continua
            continue;
        if(i != i0){
            cM(*self, swapRows, i0, i);
        }
    
        for(i = i0 + 1; i < (*self)->rows; i++){
            double pivot = (*self)->at[i0][j]; //atribui o elemento nao nulo ao pivo
            if(!equal((*self)->at[i][j], 0)){ 
                double aij = (*self)->at[i][j];

                CONCAT(mul_row_, MATRIX_TYPE)(*self, i, pivot); 
                CONCAT(add_rows_, MATRIX_TYPE)(*self, i, i0, -aij);
                //ao encontrar um elemento nao nuloa partir do pivo na linha, 
                //faz as combinações lineares necessarias para zerar os elementos abaixo do pivo
            }
        }
        i0++;
    }
    
    (*self)->was_modified = 1;
    return CMATRIX_OK;
}



static inline GENERICS_TYPE CONCAT(det_, MATRIX_TYPE)(MATRIX_TYPE** self){
    if((*self)->cols != (*self)->rows){
        printf("Determinante so e definido para matrizes quadradas\n");
        return 0;
    }
    if(!(*self)->was_modified) return (*self)->det_cache; //se nenhum metodo modificou a matriz, o detrminante ainda é o mesmo e é retornado em O(1)

    MATRIX_TYPE* TMP = NEW_MATRIX_TYPE((*self)->rows, (*self)->cols);
    TMP->vptr->cpy(&TMP, *self);

    GENERICS_TYPE det_fact = 1;
    (*self)->det_cache = 1;

    size_t i, j, i0 = 0;
    for(j = 0; j < TMP->cols; j++){
        i = i0;
        while(i < (*self)->rows && equal(TMP->at[i][j], 0)) i++;
        if(i == (*self)->rows) //procura o primeiro elemento nao nulo em uma coluna
            continue;
        if(i != i0){
            cM(TMP, swapRows, i0, i);
            (*self)->det_cache *= -1; //trocar linhas equivale a multiplicar o determinante por -1
        }
    
        for(i = i0 + 1; i < TMP->rows; i++){
            double pivot = TMP->at[i0][j]; //pivot recebe o primeiro elemento nao nulo
            if(!equal(TMP->at[i][j], 0)){
                double aij = TMP->at[i][j];

                CONCAT(mul_row_, MATRIX_TYPE)(TMP, i, pivot); // multoplica uma linha pelo pivot
                det_fact *= pivot; //custo é multiplicar o pivot pelo denominador

                CONCAT(add_rows_, MATRIX_TYPE)(TMP, i, i0, -aij); //adiciona a linha i pela linha do pivo por -aij para zerar as linhas abaixo
            }
        }
        i0++;
    }
    
    for(size_t i = 0; i < TMP->rows; i++){
        (*self)->det_cache *= TMP->at[i][i]; //apos o escalonamento o determinante é o produto da diagonal principal
    }
    (*self)->det_cache /= det_fact; //divide pelos fatores acumulados nas transformações
    
    cM(TMP, free);
    (*self)->was_modified = 0; //agora a matriz esta com o determinante correto, entao fica comom nao modificada

    return (*self)->det_cache;
}

//concatena duas matrizes de mesmo numero de linhas horizontalmente
static inline cmatrix_error CONCAT(concat_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A, const MATRIX_TYPE* B){
    if(!self || !*self || !A || !B) return CMATRIX_NULL_POINTER;
    if(A->rows != B->rows) return CMATRIX_DIMENSION_MISMATCH;

    size_t new_cols = A->cols + B->cols;

    MATRIX_TYPE* new_matrix = NEW_MATRIX_TYPE(A->rows, new_cols);
    if(!new_matrix) return CMATRIX_ALLOC_ERROR;

    for(size_t i = 0; i < A->rows; i++){
        for(size_t j = 0; j < new_cols; j++){
            if(j < A->cols){
                new_matrix->at[i][j] = A->at[i][j];
            } else {
                new_matrix->at[i][j] = B->at[i][j - A->cols];
            }
        }
    }

    cM(*self, free);
    *self = new_matrix;
    return CMATRIX_OK;
}


#if !defined(_MATRIX_DOUBLE_) && !defined(_MATRIX_FLOAT_)
static inline cmatrix_error CONCAT(rand_fill_, MATRIX_TYPE)(MATRIX_TYPE** self, GENERICS_TYPE min, GENERICS_TYPE max){
    if(!self || !*self) return CMATRIX_NULL_POINTER;
    if(min > max){
        GENERICS_TYPE t = min; min = max; max = t;
    }

    GENERICS_TYPE r = max - min + 1;
    size_t n = (*self)->rows * (*self)->cols;
    for(size_t k = 0; k < n; k++){
        (*self)->data[k] = min + rand() % r;
    }
    (*self)->was_modified = 1;
    return CMATRIX_OK;
}
#endif

#if defined(_MATRIX_DOUBLE_) || defined(_MATRIX_FLOAT_)

static inline cmatrix_error CONCAT(rand_fill_, MATRIX_TYPE)(MATRIX_TYPE** self, GENERICS_TYPE min, GENERICS_TYPE max, int only_int){
    if(!self || !*self) return CMATRIX_NULL_POINTER;
    if(min > max){
        GENERICS_TYPE t = min; min = max; max = t;
    }

    size_t r = max - min + 1;
    size_t n = (*self)->rows * (*self)->cols;
    for(size_t k = 0; k < n; k++){
        double scale = (double)rand() / (double) RAND_MAX;
        (*self)->data[k] = only_int ? min + rand() % r : (min + scale * (max - min));
    }
    (*self)->was_modified = 1;
    return CMATRIX_OK;
}

//coloca a matriz na sua Forma Escalonada Reduzida por Linhas
static inline cmatrix_error CONCAT(RREF_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* mat){
    if(!self || !*self || !mat) return CMATRIX_NULL_POINTER;
    cM(*self, cpy, mat);
    cM(*self, ref, *self);

    size_t j0 = 0, i = 0;
    //aqui aplica o mesmo raciocinio do escalonamento, mas procurando o pivo na linha e zerando todos acima e abaixo dele
    for(i = 0; i < (*self)->rows; i++){
        size_t j = j0;
        while(j < (*self)->cols && equal((*self)->at[i][j], 0)) j++;
        if(j == (*self)->cols) continue;

        double pivot = (*self)->at[i][j];
        CONCAT(mul_row_, MATRIX_TYPE)(*self, i, 1/pivot); //transforma cada pivo em 1 dividindo a linha por pivot

        for(int k = i - 1; k >= 0; k--){
            if(!equal((*self)->at[k][j], 0)){
                CONCAT(add_rows_, MATRIX_TYPE)(*self, k, i, -(*self)->at[k][j]);
            }
        }
        j0 = j + 1; // como um pivo na coluna seguinte só pode acontecer mais a direita, procuramos a partir da direita de j;
    }
    (*self)->was_modified = 1;
    return CMATRIX_OK;
}

//transforma a matriz na sua inversa
static inline cmatrix_error CONCAT(inverse_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A) {
    if (!self || !*self || !A) return 0;
    if (A->cols != A->rows) return CMATRIX_NOT_SQUARE; 

    if((*self)->cols != A->cols || (*self)->rows != A->rows){
        cM(*self, free);
        *self = NEW_MATRIX_TYPE(A->rows, A->cols);
    }
    
    MATRIX_TYPE* I = NEW_MATRIX_TYPE(A->rows, A->cols);
    cM(I, identity);

    MATRIX_TYPE* TMP = NEW_MATRIX_TYPE(A->rows, A->cols);

    cM(TMP, cpy, A);
    cM(TMP, hstack, TMP, I);
    cM(TMP, rref, TMP);

    for (size_t i = 0; i < A->rows; i++) {
        for (size_t j = 0; j < A->cols; j++) {
            GENERICS_TYPE esperado = (i == j) ? (GENERICS_TYPE)1 : (GENERICS_TYPE)0;
            if (!equal(TMP->at[i][j], esperado)) {
                cM(I, free);
                cM(TMP, free);
                return CMATRIX_SINGULAR_MATRIX;
            }
        }
    }

    for (size_t i = 0; i < I->rows; i++) {
        for (size_t j = 0; j < I->cols; j++) {
            I->at[i][j] = TMP->at[i][j + I->cols];
        }
    }

    cM(TMP, free);
    cM(*self, cpy, I);
    cM(I, free);
    return CMATRIX_OK;
}

static inline GENERICS_TYPE CONCAT(dot_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A){
    if(!self || !*self || !A)   return NAN;
    size_t n = (*self)->rows * (*self)->cols;
    if(n != A->cols * A->rows)  return NAN;

    GENERICS_TYPE dot = 0;

    for(size_t k = 0; k < n; k++){
        dot += (*self)->data[k] * A->data[k];
    }
    return dot;
}

void CONCAT(solve_, MATRIX_TYPE)(MATRIX_TYPE** x, const MATRIX_TYPE* A, const MATRIX_TYPE* b){
    MATRIX_TYPE* Ab = NEW_MATRIX_TYPE(A->rows, A->cols + 1);
    MATRIX_TYPE* Ab_s = NEW_MATRIX_TYPE(A->rows, A->cols + 1);
    MATRIX_TYPE* A_s = NEW_MATRIX_TYPE(A->rows, A->cols);
    cM(Ab, hstack, A, b);
    
    cM(Ab_s, rref, Ab);
    cM(A_s, rref, A);

    int rank_Ab = 0, rank_A = 0;
    for(int i = 0; i < A->rows; i++){
        int A_is_z = 1, Ab_is_z = 1;
        for(int j = 0; j < Ab->cols; j++){
            if(j >= A->cols){
                if(Ab_s->at[i][j] != 0) {Ab_is_z = 0; break;}
            }else{
                if(Ab_s->at[i][j] != 0) {Ab_is_z = 0; break;}
                if(A_s->at[i][j] != 0) {A_is_z = 0; break;}
            }
        }
        rank_A += !A_is_z;
        rank_Ab += !Ab_is_z;
    }
    
    if(rank_A != rank_Ab) return;
    if(rank_A == A->cols){
        
    }

}

#endif

static inline cmatrix_error CONCAT(print_, MATRIX_TYPE)(MATRIX_TYPE** self){
    if(!self || !*self) return CMATRIX_NULL_POINTER;
    
    for(size_t i = 0; i < (*self)->rows; i++){
        for(size_t j = 0; j < (*self)->cols; j++){
            if(fabs((*self)->at[i][j]) < eps) 
                (*self)->at[i][j] = 0; // resolve o problema do -0.00
            printf("%f", (*self)->at[i][j]); //printgeneric
            putchar('\t');
        }
        putchar('\n');
    } 
    putchar('\n');
    return CMATRIX_OK;
}

#endif

//aplica a função func com parametros params a cada elemento de self
static inline cmatrix_error CONCAT(apply_, MATRIX_TYPE)(MATRIX_TYPE** self, GENERICS_TYPE (*func)(GENERICS_TYPE, void*), void* params){
    if(!self || !*self) return CMATRIX_NULL_POINTER;
    if(!func) return CMATRIX_NULL_POINTER;

    size_t quant = (*self)->rows * (*self)->cols;
    GENERICS_TYPE* ptr = (*self)->data;
    for(size_t i = 0; i < quant; i++){
        ptr[i] = func(ptr[i], params);
    }
    return CMATRIX_OK;
}


static CONCAT(VTable_, GENERICS_NAME) CONCAT(vtable_, GENERICS_NAME) = {
    .cpy = CONCAT(copy_, MATRIX_TYPE),
    .scale = CONCAT(scale_, MATRIX_TYPE),
    .diag = CONCAT(diag_, MATRIX_TYPE),
    .add = CONCAT(add_, MATRIX_TYPE),
    .sub =CONCAT(sub_, MATRIX_TYPE),
    .mult = CONCAT(mul_, MATRIX_TYPE),
    .pow = CONCAT(pow_, MATRIX_TYPE),
    .free = CONCAT(delete_, MATRIX_TYPE),
    .hadmult = CONCAT(hadamard_, MATRIX_TYPE),
    .ref = CONCAT(REF_, MATRIX_TYPE),
    .det = CONCAT(det_, MATRIX_TYPE),
    .print = CONCAT(print_, MATRIX_TYPE),
    .fill = CONCAT(fill_, MATRIX_TYPE),
    .fill_rand = CONCAT(rand_fill_, MATRIX_TYPE),
    .identity = CONCAT(identity_, MATRIX_TYPE),
    .hstack = CONCAT(concat_, MATRIX_TYPE),
    .apply = CONCAT(apply_, MATRIX_TYPE),
    .swapRows = CONCAT(swap_rows_, MATRIX_TYPE),
    .swapCols = CONCAT(swap_cols_, MATRIX_TYPE),
    .reshape = CONCAT(reshape_, MATRIX_TYPE),
    .slice = CONCAT(slice_, MATRIX_TYPE),
    .getCol = CONCAT(get_col_, MATRIX_TYPE),
    .getRow = CONCAT(get_row_, MATRIX_TYPE),
    .setCol = CONCAT(set_col_, MATRIX_TYPE),
    .setRow = CONCAT(set_row_, MATRIX_TYPE),
    .atLeast = CONCAT(at_least_, MATRIX_TYPE),

    #if defined(_MATRIX_DOUBLE_) || defined(_MATRIX_FLOAT_)
    .rref = CONCAT(RREF_, MATRIX_TYPE),
    .inv = CONCAT(inverse_, MATRIX_TYPE),
    .dot = CONCAT(dot_, MATRIX_TYPE),
    #endif
    
    .transp = CONCAT(transp_, MATRIX_TYPE)
};


static inline MATRIX_TYPE* NEW_MATRIX_TYPE(size_t rows, size_t cols){
    MATRIX_TYPE* mat = (MATRIX_TYPE*) malloc(sizeof(MATRIX_TYPE));
    if(!mat){
        printf("memoria insuficiente\n");
        return NULL;
    }

    mat->data = (GENERICS_TYPE*) calloc(rows*cols, sizeof(GENERICS_TYPE));
    if(!mat->data) {
        printf("memoria insuficiente\n");
        free(mat);
        return NULL;
    }

    mat->at = (GENERICS_TYPE**) malloc(rows*sizeof(GENERICS_TYPE*));
    if(!mat->at){
        printf("Memoria insuficiente\n");
        free(mat->data);
        free(mat);
        return NULL;
    }

    mat->rows = rows;
    mat->cols = cols;

    for(size_t i = 0; i < mat->rows; i++)
        mat->at[i] = mat->data + i*mat->cols;
        
    mat->vptr = &CONCAT(vtable_, GENERICS_NAME);
    mat->was_modified = 1;

    return mat;
}

//constroi a matriz e inicializa com dados do teclado
static inline MATRIX_TYPE* CONCAT(get_, MATRIX_TYPE)(size_t rows, size_t cols){
    MATRIX_TYPE* mat = NEW_MATRIX_TYPE(rows, cols);
    if(!mat){
        printf("Memoria insuficiente\n");
        return NULL;
    }
    for(size_t i = 0; i < mat->rows; i++){
        for(size_t j = 0; j < mat->cols; j++){
            input(mat->at[i][j]); // entrada generica com _Generic
        }
    }
    return mat;
}

//transforma um array 2d em uma matrix
static inline MATRIX_TYPE* CONCAT(MATRIX_TYPE,_init2D)(size_t m, size_t n, const GENERICS_TYPE mat[m][n]){
    MATRIX_TYPE* A = NEW_MATRIX_TYPE(m, n);
    if(!A){
        printf("Memoria insuficiente!\n");
        return NULL;
    }
    memcpy(A->data, mat, n*m*sizeof(GENERICS_TYPE));
    return A;
}

static inline cmatrix_error CONCAT(MATRIX_TYPE, _set2D)(MATRIX_TYPE** A, size_t m, size_t n, const GENERICS_TYPE mat[m][n]){

    if(!A || !*A || !mat) return CMATRIX_NULL_POINTER;

    if((*A)->rows != m || (*A)->cols != n){
        *A =  CONCAT(MATRIX_TYPE,_init2D)(m, n, mat);
        return CMATRIX_OK;
    }
    memcpy((*A)->data, mat, m*n*sizeof(GENERICS_TYPE));
    return CMATRIX_OK;
}


