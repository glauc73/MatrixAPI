# cmatrix

Uma biblioteca *header-only* de matrizes genéricas em **C11 puro**, com interface orientada a objetos via macro `cM` e simulação de virtual tables.

```bash
gcc -std=c11 meu_programa.c -Iinclude -lm -o meu_programa
```

---

## Estrutura de arquivos

```
cmatrix/
├── include/
│   ├── genericsBase.h    # infraestrutura de tipos genéricos e macros base
│   ├── genericMatrix.h   # implementação do TAD e de todas as funções
│   └── cMatrix.h         # ponto de entrada — inclua só este
├── example.c             # demonstração completa (soluções básicas de PPL)
├── build.sh              # compilação Linux/macOS
└── build.bat             # compilação Windows
```

---

## Uso básico

```c
#include "include/cMatrix.h"
```

### Tipos disponíveis

```c
cmatrix_int*    A = new_cmatrix_int(3, 3);
cmatrix_float*  B = new_cmatrix_float(3, 3);
cmatrix_double* C = new_cmatrix_double(3, 3);
```

Construtores auxiliares:

```c
// lê os valores do stdin
cmatrix_float* A = get_cmatrix_float(3, 3);

// constrói a partir de um array 2D existente
float arr[2][3] = {{1,2,3},{4,5,6}};
cmatrix_float* B = cmatrix_float_init2D(2, 3, arr);
```

### A macro `cM`

Toda operação é feita via `cM(matriz, método, ...)`:

```c
cM(A, print);
cM(C, mult, A, B);   // C = A · B
```

---

## Funcionalidades

| Operação | Macro `cM` | O que faz | Expressão aritmética |
|----------|------------|-----------|----------------------|
| Imprimir | `cM(A, print)` | Exibe a matriz formatada no terminal | — |
| Copiar | `cM(B, cpy, A)` | Copia o conteúdo de A em B | B ← A |
| Preencher constante | `cM(A, fill, v)` | Preenche todos os elementos com `v` | aᵢⱼ = v |
| Preencher aleatório (`int`) | `cM(A, fill_rand, min, max)` | Preenche com inteiros aleatórios em `[min, max]` | aᵢⱼ ~ U(min, max) ∩ ℤ |
| Preencher aleatório (`float`/`double`) | `cM(A, fill_rand, min, max, inteiros)` | Preenche aleatoriamente; `inteiros=1` restringe a ℤ | aᵢⱼ ~ U(min, max) |
| Identidade | `cM(A, identity)` | Transforma A na identidade (exige matriz quadrada) | A = I |
| Transposta | `cM(B, transp, A)` | Armazena a transposta de A em B | B = Aᵀ |
| Adição | `cM(C, add, A, B)` | Soma elemento a elemento | C = A + B |
| Subtração | `cM(C, sub, A, B)` | Subtrai elemento a elemento | C = A − B |
| Multiplicação | `cM(C, mult, A, B)` | Produto matricial | C = A · B |
| Produto de Hadamard | `cM(C, hadmult, A, B)` | Multiplica elemento a elemento | cᵢⱼ = aᵢⱼ · bᵢⱼ |
| Escalonamento | `cM(B, scale, A, alpha)` | Multiplica todos os elementos por `alpha` | B = α · A |
| Potência | `cM(B, pow, A, n)` | Eleva A à potência inteira n | B = Aⁿ |
| Produto interno | `float d = cM(A, dot, B)` | Soma dos produtos elemento a elemento | d = Σ aᵢⱼ · bᵢⱼ |
| Diagonal | `cM(B, diag, A)` | Extrai diagonal de A ou monta diagonal de vetor | B = diag(a₁₁, a₂₂, …) |
| REF | `cM(B, ref, A)` | Forma escalonada por linhas | B = ref(A) |
| RREF | `cM(B, rref, A)` | Forma escalonada reduzida (pivôs = 1) | B = rref(A) |
| Inversa | `cM(B, inv, A)` | Inversa via Gauss-Jordan | B = A⁻¹ |
| Determinante | `float d = cM(A, det)` | Determinante em O(n³), com cache | d = det(A) |
| Concatenação horizontal | `cM(C, hstack, A, B)` | Junta A e B lado a lado | C = [A \| B] |
| Aplicar função | `cM(A, apply, f, params)` | Aplica `f(aᵢⱼ, params)` a cada elemento | aᵢⱼ ← f(aᵢⱼ) |
| Redimensionar | `cM(A, reshape, m, n)` | Muda dimensões mantendo os dados | — |
| Trocar linhas | `cM(A, swapRows, i, j)` | Troca as linhas i e j em O(1) | Lᵢ ↔ Lⱼ |
| Trocar colunas | `cM(A, swapCols, i, j)` | Troca as colunas i e j | Cᵢ ↔ Cⱼ |
| Fatiar | `cM(B, slice, A, i0, i1, j0, j1)` | Extrai submatriz de A | B = A[i0:i1, j0:j1] |
| Obter coluna | `cM(B, getCol, A, j)` | Extrai a coluna j de A como vetor coluna | B = A[:, j] |
| Obter linha | `cM(B, getRow, A, i)` | Extrai a linha i de A como vetor linha | B = A[i, :] |
| Definir coluna | `cM(A, setCol, v, j)` | Substitui a coluna j de A pelo vetor v | A[:, j] ← v |
| Definir linha | `cM(A, setRow, v, i)` | Substitui a linha i de A pelo vetor v | A[i, :] ← v |
| Todos ≥ escalar | `int ok = cM(A, atLeast, v)` | Retorna 1 se todos os elementos de A são ≥ v | ∀ aᵢⱼ : aᵢⱼ ≥ v |
| Liberar | `cM(A, free)` | Libera a memória alocada | — |

> `rref`, `inv`, `dot` e `atLeast` só existem para `cmatrix_float` e `cmatrix_double`.  
> `fill_rand` para `int` não recebe o parâmetro `inteiros`.

---

## Tratamento de erros

Todas as funções que retornam `void` foram convertidas para retornar `cmatrix_error`:

| Código | Valor | Quando ocorre |
|--------|-------|---------------|
| `CMATRIX_OK` | 0 | Operação realizada com sucesso |
| `CMATRIX_NULL_POINTER` | 1 | Um ponteiro nulo foi passado como argumento |
| `CMATRIX_ALLOC_ERROR` | 2 | Falha na alocação de memória |
| `CMATRIX_DIMENSION_MISMATCH` | 3 | Dimensões incompatíveis entre as matrizes |
| `CMATRIX_NOT_SQUARE` | 4 | A operação exige matriz quadrada |
| `CMATRIX_SINGULAR_MATRIX` | 5 | A matriz não possui inversa |
| `CMATRIX_INCOMPATIBLE_RESHAPE` | 6 | O novo shape tem número de elementos diferente |
| `CMATRIX_INVALID_DIAG_INPUT` | 7 | `diag` recebeu matriz que não é quadrada nem vetor |
| `CMATRIX_INDEX_OUT_OF_RANGE` | 8 | Índice de linha ou coluna fora do intervalo válido |

As funções com retorno próprio — `det`, `dot`, `inv` e `atLeast` — não foram alteradas.

### Exemplo de uso

```c
cmatrix_error err;

cmatrix_float* A = new_cmatrix_float(3, 3);
cmatrix_float* B = new_cmatrix_float(2, 4);
cmatrix_float* C = new_cmatrix_float(3, 3);

cM(A, fill_rand, -5, 5, 1);
cM(B, fill_rand, -5, 5, 1);

err = cM(C, add, A, B);
if (err != CMATRIX_OK) {
    printf("Erro %d\n", err);  // CMATRIX_DIMENSION_MISMATCH
}

err = cM(C, mult, A, A);
if (err != CMATRIX_OK) {
    cM(A, free); cM(B, free); cM(C, free);
    return err;
}

cM(A, free); cM(B, free); cM(C, free);
```

---

## Exemplo completo

O `example.c` demonstra a geração de todas as soluções básicas (e básicas viáveis) de um problema de programação linear a partir de uma matriz aumentada `[A | b]`:

```c
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "include/cMatrix.h"

typedef enum {
    BASIC_SOLUTION = 0,
    BASIC_FEASIBLE_SOLUTION
} PPLSolution_t;

// Gera todas as soluções básicas (ou básicas viáveis) de Ab = [A | b]
cmatrix_error generate_all_BS(cmatrix_float*** solutions, size_t* n_solutions,
                               cmatrix_float* Ab, PPLSolution_t bfs) {
    if (!Ab) return CMATRIX_NULL_POINTER;

    size_t n = Ab->cols - 1, m = Ab->rows;

    cmatrix_float* b = new_cmatrix_float(m, 1);
    cM(b, getCol, Ab, n);

    cmatrix_float* v[n];
    for (int i = 0; i < n; i++) {
        v[i] = new_cmatrix_float(m, 1);
        cM(v[i], getCol, Ab, i);
    }

    size_t n_combs = comb(n, m);
    *solutions = malloc(n_combs * sizeof(cmatrix_float*));
    for (size_t i = 0; i < n_combs; i++)
        (*solutions)[i] = new_cmatrix_float(m, 1);

    int64_t mask = (1LL << m) - 1;
    cmatrix_float* Bb = new_cmatrix_float(m, m + 1);
    cmatrix_float* x  = new_cmatrix_float(n, 1);
    int cnt = 0;

    for (int i = 0; i < n_combs; i++) {
        // monta submatriz aumentada [B | b] para a combinação de colunas dada pela máscara
        cM(Bb, setCol, b, m);
        int col = 0;
        for (int j = 0; j < n; j++)
            if ((mask >> j) & 1) cM(Bb, setCol, v[j], col++);

        cM(Bb, rref, Bb);

        // verifica se a submatriz tem posto completo (nenhum pivô nulo na diagonal)
        int full_rank = 1;
        for (size_t r = 0; r < Bb->rows; r++)
            if (fabsf(Bb->at[r][r]) < 1e-6) { full_rank = 0; break; }

        if (full_rank) {
            // reconstrói x completo (variáveis não-básicas = 0)
            int top = 0;
            for (int j = 0; j < n; j++) {
                x->data[j] = ((mask >> j) & 1)
                             ? Bb->data[top++ * Bb->cols + m]
                             : 0;
            }
            // se bfs=BASIC_FEASIBLE_SOLUTION, filtra apenas x >= 0
            if (!bfs || cM(x, atLeast, 0)) {
                cM((*solutions)[cnt], cpy, x);
                cnt++;
            }
        }
        mask = next_combination(mask);
    }

    if (cnt == 0) return CMATRIX_SINGULAR_MATRIX;

    *solutions = realloc(*solutions, cnt * sizeof(cmatrix_float*));
    *n_solutions = cnt;

    cM(Bb, free); cM(b, free); cM(x, free);
    for (size_t i = 0; i < n; i++) cM(v[i], free);
    return CMATRIX_OK;
}

int main() {
    size_t m, n;
    PPLSolution_t tipo;
    scanf("%d", &tipo);          // 0 = todas as SB, 1 = apenas SBV (x >= 0)
    scanf(" %zu %zu", &m, &n);

    cmatrix_float* Ab = get_cmatrix_float(m, n);

    cmatrix_float** x = NULL;
    size_t n_x;
    cmatrix_error err = generate_all_BS(&x, &n_x, Ab, tipo);

    if (err) { printf("Erro: %d\n", err); goto fim; }

    printf("%zu solucoes encontradas:\n", n_x);
    for (int i = 0; i < n_x; i++) cM(x[i], print);

    for (int i = 0; i < n_x; i++) cM(x[i], free);
    free(x);

fim:
    cM(Ab, free);
    return err;
}
```

**Compilar e rodar:**

```bash
gcc -std=c11 example.c -Iinclude -lm -o main
echo "0\n5 8\n1 0 1 0 0 0 0 40\n0 1 0 1 0 0 0 60\n3 2 0 0 1 0 0 180\n0 1 0 0 0 -1 0 10\n1 1 0 0 0 0 -1 20" | ./main
```

---

## Compilação

```bash
# Linux / macOS
gcc -std=c11 meu_programa.c -Iinclude -lm -o meu_programa

# Windows (MinGW)
gcc -std=c11 meu_programa.c -Iinclude -lm -o meu_programa.exe
```

Ou use os scripts incluídos:

```bash
./build.sh meu_programa.c   # Linux/macOS
build.bat  meu_programa.c   # Windows
```

---

*por [glauc73](https://github.com/glauc73)*
