#include "genericsBase.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

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

#define MATRIX_TYPE CONCAT(matrix_, GENERICS_NAME)


typedef struct MATRIX_TYPE{
    GENERICS_TYPE* data;
    GENERICS_TYPE** at;
    size_t rows, cols;
    //metodos in-place. eles nao retornam as matrizes para evitar memory leak
    void (*mult)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*, const struct MATRIX_TYPE*);
    void (*add)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*, const struct MATRIX_TYPE*);
    void (*sub)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*, const struct MATRIX_TYPE*);
    void (*pow)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*, uint16_t);
    void (*scale)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*, GENERICS_TYPE);
    void (*cpy)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*);
    void (*free)(struct MATRIX_TYPE**);
    void (*ref)(struct MATRIX_TYPE*);
    void (*identity)(struct MATRIX_TYPE*);
    void (*print)(struct MATRIX_TYPE*);
    void (*hstack)(struct MATRIX_TYPE**, const struct MATRIX_TYPE*);
    void (*apply)(struct MATRIX_TYPE**, GENERICS_TYPE (*)(GENERICS_TYPE, void*), void*);
    GENERICS_TYPE (*det)(struct MATRIX_TYPE*);

    #if defined(_MATRIX_DOUBLE_) || defined(_MATRIX_FLOAT_) //metodos que só fazem sentido para double
    void (*rref)(struct MATRIX_TYPE*);
    int (*inv)(struct MATRIX_TYPE**);
    #endif

    #if !defined(__INTELLISENSE__) && !defined(__clangd__) //Atributos privados
    GENERICS_TYPE det_cache; //atributo que guarda o determinante calculado para evitar recalculo
    int8_t was_modified; //flag alterada cada vez que um metodo modifica a matriz para recalcular o det
    #endif
}MATRIX_TYPE;


#define NEW_MATRIX_TYPE CONCAT(new_, MATRIX_TYPE)
static inline MATRIX_TYPE* NEW_MATRIX_TYPE(size_t rows, size_t cols);

static inline void CONCAT(swap_rows_, MATRIX_TYPE)(MATRIX_TYPE* mat, size_t r1, size_t r2){//troca linhas em O(1)
    if(r1 >= mat->rows || r2 >= mat->rows) return;
    GENERICS_TYPE* temp = mat->at[r1];
    mat->at[r1] = mat->at[r2];
    mat->at[r2] = temp;
}

#if !defined(__INTELLISENSE__) && !defined(__clangd__)


static inline void CONCAT(delete_, MATRIX_TYPE)(MATRIX_TYPE** self){ //destrói corretamente a matriz
    if(self == NULL || *self == NULL) return;
    free((*self)->data);
    free((*self)->at);
    free((*self));
    *self = NULL;
}

//copia a matriz B para self
static inline void CONCAT(copy_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* B){
    if((*self)->rows != B->rows || (*self)->cols != B->cols){
        printf("Matrizes precisam ter dimensões iguais para serem copiadas!\n");
        return;
    } 
    //como data aponta para um bloco contiguo, memcpy copia tudo de uma vez
    memcpy((*self)->data, B->data, B->rows * B->cols * sizeof(GENERICS_TYPE));

    (*self)->was_modified = B->was_modified;

    if(!(*self)->was_modified){
        (*self)->det_cache = B->det_cache;
    }
}

static inline void CONCAT(add_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A, const MATRIX_TYPE* B){
    if(A->rows != B->rows || A->cols != B->cols){
        printf("O numero de linhas e colunas das matrizes precisam ser iguais\n");
        return;
    }
    //soma elemento a elemento e atribui em self
    for(size_t i = 0; i < (*self)->rows; i++){
        for(size_t j = 0; j < (*self)->cols; j++){
            (*self)->at[i][j] = A->at[i][j] + B->at[i][j];
        }
    }
    (*self)->was_modified = 1;
}

static inline void CONCAT(sub_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A, const MATRIX_TYPE* B){
    if(A->rows != B->rows || A->cols != B->cols){
        printf("O numero de linhas e colunas das matrizes precisam ser iguais\n");
        return;
    }

    for(size_t i = 0; i < (*self)->rows; i++){
        for(size_t j = 0; j < (*self)->cols; j++){
            (*self)->at[i][j] = A->at[i][j] - B->at[i][j];
        }
    }
    (*self)->was_modified = 1;
}

static inline void CONCAT(scale_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A, GENERICS_TYPE k){
    if(!self || !A){
        printf("ponteiro nulo!\n");
        return;
    }

    for(size_t i = 0; i < (*self)->rows; i++){
        for(size_t j = 0; j < (*self)->cols; j++){
            (*self)->at[i][j] = k * A->at[i][j];
        }
    }
    (*self)->was_modified = 1;

}

static inline void CONCAT(mul_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A, const MATRIX_TYPE* B){
    if(!self || !A || !B){
        printf("ponteiro nulo!\n");
        return;
    }
    if(A->cols != B->rows){
        printf("O numero de colunas de A(%zu x %zu) deve ser igual ao numero de linhas de B(%zu x %zu)\n",
             A->rows, A->cols, B->rows, B->cols
        );
        return;
    }
    
    MATRIX_TYPE* TMP = NEW_MATRIX_TYPE((*self)->rows, (*self)->cols);

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
    (*self)->cpy(self, TMP);
    
    TMP->free(&TMP);

    (*self)->was_modified = 1;
}

static inline void CONCAT(identity_, MATRIX_TYPE)(MATRIX_TYPE* self){
    for(size_t i = 0; i < self->rows; i++){
        for(size_t j = 0; j < self->cols; j++){
            self->at[i][j] = (i == j);
        }
    }
    self->det_cache = 1;
}

static inline void CONCAT(pow_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A, uint16_t n){
    if(!self || !A){
        printf("ponteiro nulo!\n");
        return;
    }
    if(A->cols != A->rows) {
        printf("Potenciacao so funciona em matrizes quadradas\n");
        return ;
    }

    MATRIX_TYPE* base = NEW_MATRIX_TYPE(A->rows, A->cols);  //base para permitir in-place
    base->cpy(&base, A);

    //exponenciação rapida iterativa em O(n³logn)
    while(n > 0){
        if(n % 2 == 1){
           (*self)->mult(self, *self, base);
        }
        base->mult(&base, base, base);
        n /= 2;
    }
    CONCAT(delete_, MATRIX_TYPE)(&base);

    (*self)->was_modified = 1;
}

static inline void CONCAT(add_rows_, MATRIX_TYPE)(MATRIX_TYPE* mat, size_t r1, size_t r2, double k){
    if(!mat){
        printf("ponteiro Nulo\n");
    }
    //adiciona uma linha ao produto de outra por um escalar
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
static inline void CONCAT(REF_, MATRIX_TYPE)(MATRIX_TYPE* mat){
    size_t i, j, i0 = 0;
    for(j = 0; j < mat->cols; j++){
        i = i0; //inicializa a primeira linha para buscar
        while(i < mat->rows && equal(mat->at[i][j], 0)) i++; //procura o primeiro elemento nao nulo na linha
        if(i == mat->rows) //se nao encontou, continua
            continue;
        if(i != i0){
            CONCAT(swap_rows_, MATRIX_TYPE)(mat, i0, i); //se encontrar, troca com a primeira linha dessa iteração
        }
    
        for(i = i0 + 1; i < mat->rows; i++){
            double pivot = mat->at[i0][j]; //atribui o elemento nao nulo ao pivo
            if(!equal(mat->at[i][j], 0)){ 
                double aij = mat->at[i][j];

                CONCAT(mul_row_, MATRIX_TYPE)(mat, i, pivot); 
                CONCAT(add_rows_, MATRIX_TYPE)(mat, i, i0, -aij);
                //ao encontrar um elemento nao nuloa partir do pivo na linha, 
                //faz as combinações lineares necessarias para zerar os elementos abaixo do pivo
            }
        }
        i0++;
    }
    
    mat->was_modified = 1;
}

//ultiliza o escalonamento para triangularizar a matriz e tornar o determinante um simples produto da diagonal principal
//calcula em O(n³) e Θ(1)
static inline GENERICS_TYPE CONCAT(det_, MATRIX_TYPE)(MATRIX_TYPE* self){
    if(self->cols != self->rows){
        printf("Determinante so e definido para matrizes quadradas\n");
        return 0;
    }
    if(!self->was_modified) return self->det_cache; //se nenhum metodo modificou a matriz, o detrminante ainda é o mesmo e é retornado em O(1)

    MATRIX_TYPE* TMP = NEW_MATRIX_TYPE(self->rows, self->cols);
    TMP->cpy(&TMP, self);

    GENERICS_TYPE det_fact = 1;
    self->det_cache = 1;

    size_t i, j, i0 = 0;
    for(j = 0; j < TMP->cols; j++){
        i = i0;
        while(i < self->rows && equal(TMP->at[i][j], 0)) i++;
        if(i == self->rows) //procura o primeiro elemento nao nulo em uma coluna
            continue;
        if(i != i0){
            CONCAT(swap_rows_, MATRIX_TYPE)(TMP, i0, i); //caso encontre, troca se necessario
            self->det_cache *= -1; //trocar linhas equivale a multiplicar o determinante por -1
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
        self->det_cache *= TMP->at[i][i]; //apos o escalonamento o determinante é o produto da diagonal principal
    }
    self->det_cache /= det_fact; //divide pelos fatores acumulados nas transformações
    
    TMP->free(&TMP);
    self->was_modified = 0; //agora a matriz esta com o determinante correto, entao fica comom nao modificada

    return self->det_cache;
}

//concatena duas matrizes de mesmo numero de linhas horizontalmente
static inline void CONCAT(concat_, MATRIX_TYPE)(MATRIX_TYPE** self, const MATRIX_TYPE* A){
    if(A->rows != (*self)->rows){
        printf("A concatenacao horizontal de matrizes exige o mesmo numero de linhas");
        return;
    }

    size_t new_cols = (*self)->cols + A->cols;

    MATRIX_TYPE* new_matrix = NEW_MATRIX_TYPE((*self)->rows, new_cols);

    for(size_t i = 0; i < (*self)->rows; i++){
        for(size_t j = 0; j < new_cols; j++){
            if(j < (*self)->cols){
                new_matrix->at[i][j] = (*self)->at[i][j];
            } else {
                new_matrix->at[i][j] = A->at[i][j - (*self)->cols];
            }
        }
    }

    (*self)->free(self);
    *self = new_matrix;
}


#if defined(_MATRIX_DOUBLE_) || defined(_MATRIX_FLOAT_)

//coloca a matriz na sua Forma Escalonada Reduzida por Linhas
static inline void CONCAT(RREF_, MATRIX_TYPE)(MATRIX_TYPE* mat){
    mat->ref(mat);
    size_t j0 = 0, i = 0;
    //aqui aplica o mesmo raciocinio do escalonamento, mas procurando o pivo na linha e zerando todos acima e abaixo dele
    for(i = 0; i < mat->rows; i++){
        size_t j = j0;
        while(j < mat->cols && equal(mat->at[i][j], 0)) j++;
        if(j == mat->cols) continue;

        double pivot = mat->at[i][j];
        CONCAT(mul_row_, MATRIX_TYPE)(mat, i, 1/pivot); //transforma cada pivo em 1 dividindo a linha por pivot

        for(int k = i - 1; k >= 0; k--){
            if(!equal(mat->at[k][j], 0)){
                CONCAT(add_rows_, MATRIX_TYPE)(mat, k, i, -mat->at[k][j]);
            }
        }
        j0 = j + 1; // como um pivo na coluna seguinte só pode acontecer mais a direita, procuramos a partir da direita de j;
    }
    mat->was_modified = 1;
}

//transforma a matriz na sua inversa
static inline int CONCAT(inverse_, MATRIX_TYPE)(MATRIX_TYPE** self){
    if((*self)->cols != (*self)->rows){
        printf("Erro: A matriz nao é quadrada, portanto nao possui inversa\n");
        return 0;
    }
    //concatena a matriz identidade em self e aplica rref em ambas com os mesmos parametros
    //a matriz identidade vai se transformar na inversa, caso a matriz seja invertivel
    MATRIX_TYPE* I = NEW_MATRIX_TYPE((*self)->rows, (*self)->cols);
    CONCAT(identity_, MATRIX_TYPE)(I);
    
    MATRIX_TYPE* TMP = NEW_MATRIX_TYPE((*self)->rows, (*self)->cols);
    TMP->cpy(&TMP, *self);

    CONCAT(concat_, MATRIX_TYPE)(&TMP, I);
    TMP->rref(TMP);

    for(size_t i = 0; i < TMP->rows; i++){
        if(!equal(TMP->at[i][i], 1)){
            printf("Erro: Nao possui matriz inversa!\n"); 
            //se tiver algum elemento diferente de 1 na diagonal, entao ela nao admite inversa,
            //pois indice que tem uma linha nula, logo det = 0
            I->free(&I);
            TMP->free(&TMP);
            return 0;
        }
    }

    //atribui a I apenas a parte direita de TMP, que é a inversa
    for(size_t i = 0; i < I->rows; i++){
        for(size_t j = 0; j < I->cols; j++){
            I->at[i][j] = TMP->at[i][j + I->cols];
        }
    }

    (*self)->free(self);
    TMP->free(&TMP);

    *self = I; //aponta self para I
    return 1;
}

#endif

static inline void CONCAT(print_, MATRIX_TYPE)(MATRIX_TYPE* self){
    if(!self) {
        printf("Ponteiro Nulo");
        return;
    }

    for(size_t i = 0; i < self->rows; i++){
        for(size_t j = 0; j < self->cols; j++){
            if(abs_gen(self->at[i][j]) < eps) 
                self->at[i][j] = 0; // resolve o problema do -0.00
            printg(self->at[i][j]); //printgeneric
            putchar('\t');
        }
        putchar('\n');
    }
    
}

#endif

//aplica a função func com parametros params a cada elemento de self
static inline void CONCAT(apply_, MATRIX_TYPE)(MATRIX_TYPE** self, GENERICS_TYPE (*func)(GENERICS_TYPE, void*), void* params){
    if(self == NULL || *self == NULL) {
        printf("Ponteiro Nulo\n");
        return;
    }

    size_t quant = (*self)->rows * (*self)->cols;
    GENERICS_TYPE* ptr = (*self)->data;
    for(size_t i = 0; i < quant; i++){
        ptr[i] = func(ptr[i], params);
    }
}

static inline MATRIX_TYPE* NEW_MATRIX_TYPE(size_t rows, size_t cols){ constructor da matriz
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
        
    mat->cpy = CONCAT(copy_, MATRIX_TYPE);
    mat->scale = CONCAT(scale_, MATRIX_TYPE);
    mat->add = CONCAT(add_, MATRIX_TYPE);
    mat->sub = CONCAT(sub_, MATRIX_TYPE);
    mat->mult = CONCAT(mul_, MATRIX_TYPE);
    mat->pow = CONCAT(pow_, MATRIX_TYPE);
    mat->free = CONCAT(delete_, MATRIX_TYPE);
    mat->ref = CONCAT(REF_, MATRIX_TYPE);
    mat->det = CONCAT(det_, MATRIX_TYPE);
    mat->print = CONCAT(print_, MATRIX_TYPE);
    mat->identity = CONCAT(identity_, MATRIX_TYPE);
    mat->hstack = CONCAT(concat_, MATRIX_TYPE);
    mat->apply = CONCAT(apply_, MATRIX_TYPE);

    #if defined(_MATRIX_DOUBLE_) || defined(_MATRIX_FLOAT_)
    mat->rref = CONCAT(RREF_, MATRIX_TYPE);
    mat->inv = CONCAT(inverse_, MATRIX_TYPE);
    #endif
    
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

//constroi a matriz e inicializa ela com uma matriz A
static inline MATRIX_TYPE* CONCAT(init_, MATRIX_TYPE)(MATRIX_TYPE* A){
    if(!A){
        printf("Ponteiro Nulo\n");
        return NULL;
    }
    MATRIX_TYPE* TMP = NEW_MATRIX_TYPE(A->rows, A->cols);
    if(!TMP){
        printf("Memoria insuficiente!\n");
        return NULL;
    }
    TMP->cpy(&TMP, A);
    return TMP;
}

//transforma um array 2d em uma matrix
static inline MATRIX_TYPE* CONCAT(matrix_from_array_2d)(size_t m, size_t n, const GENERICS_TYPE mat[m][n]){
    MATRIX_TYPE* A = NEW_MATRIX_TYPE(m, n);
    if(!A){
        printf("Memoria insuficiente!\n");
        return NULL;
    }
    memcpy(A->data, mat, n*m*sizeof(GENERICS_TYPE));
    return A;
}