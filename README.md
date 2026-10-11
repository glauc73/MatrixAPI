# cmatrix

Biblioteca *header-only* de matrizes genéricas em **C11 puro**, com interface orientada a objetos via macro `cM`.

```bash
gcc -std=c11 meu_programa.c -Iinclude -lm -o meu_programa
```

---

## Estrutura de arquivos

```
cmatrix/
├── include/
│   ├── genericsBase.h   # infraestrutura de tipos genéricos
│   ├── genericMatrix.h  # implementação do TAD e das funções
│   └── Matrix.h         # ponto de entrada — inclua só este
├── example.c
└── in.txt
```

---

## Uso básico

```c
#include "include/Matrix.h"
```

### Tipos disponíveis

```c
cmatrix_int*    A = new_cmatrix_int(3, 3);
cmatrix_float*  B = new_cmatrix_float(3, 3);
cmatrix_double* C = new_cmatrix_double(3, 3);
```

Também há construtores auxiliares:

```c
// lê os valores do stdin
cmatrix_float* A = get_cmatrix_float(3, 3);

// cria uma cópia de outra matriz
cmatrix_float* B = init_cmatrix_float(A);

// constrói a partir de um array 2D existente
float arr[2][3] = {{1,2,3},{4,5,6}};
cmatrix_float* C = cmatrix_float_from_array_2d(2, 3, arr);
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
| Copiar | `cM(B, cpy, A)` | Copia o conteúdo de A em B, redimensionando se necessário | B ← A |
| Preencher constante | `cM(A, fill, v)` | Preenche todos os elementos com o valor `v` | aᵢⱼ = v |
| Preencher aleatório (`int`) | `cM(A, fill_rand, min, max)` | Preenche com inteiros aleatórios em `[min, max]` | aᵢⱼ ~ U(min, max) ∩ ℤ |
| Preencher aleatório (`float`/`double`) | `cM(A, fill_rand, min, max, inteiros)` | Preenche com valores aleatórios; `inteiros=1` restringe a ℤ | aᵢⱼ ~ U(min, max) |
| Identidade | `cM(A, identity)` | Transforma A na matriz identidade (deve ser quadrada) | A = I |
| Transposta | `cM(B, transp, A)` | Armazena a transposta de A em B | B = Aᵀ |
| Adição | `cM(C, add, A, B)` | Soma elemento a elemento | C = A + B |
| Subtração | `cM(C, sub, A, B)` | Subtrai elemento a elemento | C = A − B |
| Multiplicação | `cM(C, mult, A, B)` | Produto matricial (exige colunas de A = linhas de B) | C = A · B |
| Produto de Hadamard | `cM(C, hadmult, A, B)` | Multiplica elemento a elemento (exige mesmas dimensões) | cᵢⱼ = aᵢⱼ · bᵢⱼ |
| Escalonamento | `cM(B, scale, A, alpha)` | Multiplica todos os elementos por um escalar `alpha` | B = α · A |
| Potência | `cM(B, pow, A, n)` | Eleva A à potência inteira n (exige matriz quadrada; n=0 retorna I) | B = Aⁿ |
| Produto interno | `float d = cM(A, dot, B)` | Soma dos produtos elemento a elemento (mesmas dimensões) — só `float`/`double` | d = Σ aᵢⱼ · bᵢⱼ |
| Diagonal | `cM(B, diag, A)` | Se A é quadrada: extrai diagonal; se A é vetor linha/coluna: monta matriz diagonal | B = diag(a₁₁, a₂₂, …, aₙₙ) |
| REF | `cM(B, ref, A)` | Forma escalonada por linhas via eliminação gaussiana (pivôs não normalizados) | B = ref(A) |
| RREF | `cM(B, rref, A)` | Forma escalonada reduzida (pivôs = 1, zeros acima e abaixo) — só `float`/`double` | B = rref(A) |
| Inversa | `int ok = cM(B, inv, A)` | Inversa via Gauss-Jordan; retorna 1 se existe, 0 se singular — só `float`/`double` | B = A⁻¹, A · B = I |
| Determinante | `float d = cM(A, det)` | Determinante via escalonamento gaussiano O(n³); cacheado enquanto A não mudar | d = det(A) |
| Concatenação horizontal | `cM(C, hstack, A, B)` | Junta A e B lado a lado (exige mesmo número de linhas) | C = [A \| B] |
| Aplicar função | `cM(A, apply, f, params)` | Aplica `f(aᵢⱼ, params)` a cada elemento de A in-place | aᵢⱼ ← f(aᵢⱼ) |
| Redimensionar | `cM(A, reshape, m, n)` | Muda as dimensões mantendo os dados (exige m×n = linhas×colunas originais) | — |
| Trocar linhas | `cM(A, swapRows, i, j)` | Troca as linhas i e j em O(1) | Lᵢ ↔ Lⱼ |
| Trocar colunas | `cM(A, swapCols, i, j)` | Troca as colunas i e j | Cᵢ ↔ Cⱼ |
| Liberar | `cM(A, free)` | Libera a memória e anula o ponteiro | — |

> `rref`, `inv` e `dot` só existem para `cmatrix_float` e `cmatrix_double`.  
> `fill_rand` para `int` não recebe o parâmetro `inteiros`.

---

## Exemplos por operação

### Adição e subtração

```c
cmatrix_float* A = new_cmatrix_float(3, 3);
cmatrix_float* B = new_cmatrix_float(3, 3);
cmatrix_float* C = new_cmatrix_float(3, 3);

cM(A, fill_rand, -5, 5, 1);
cM(B, fill_rand, -5, 5, 1);

cM(C, add, A, B);   // C = A + B
cM(C, sub, A, B);   // C = A - B
```

### Escalonamento (multiplicação por escalar)

```c
cmatrix_float* B = new_cmatrix_float(3, 3);
cM(B, scale, A, 2.5f);   // B = 2.5 * A
```

### Potência de matriz

```c
cmatrix_float* B = new_cmatrix_float(3, 3);
cM(B, pow, A, 4);   // B = A⁴
// pow com n=0 retorna a identidade
```

### Produto de Hadamard

```c
cmatrix_float* C = new_cmatrix_float(3, 3);
cM(C, hadmult, A, B);   // cᵢⱼ = aᵢⱼ * bᵢⱼ
```

### Produto interno (dot)

```c
// A e B podem ter qualquer forma, desde que tenham o mesmo total de elementos
float d = cM(A, dot, B);   // d = Σ aᵢⱼ * bᵢⱼ
printf("dot = %.4f\n", d);
```

### Identidade

```c
cmatrix_float* I = new_cmatrix_float(4, 4);
cM(I, identity);   // I vira a identidade 4×4
```

### Diagonal

```c
// Extraindo a diagonal de uma matriz quadrada
cmatrix_float* D = new_cmatrix_float(3, 3);
cM(D, diag, A);   // D = diag(a₁₁, a₂₂, a₃₃)

// Construindo matriz diagonal a partir de vetor
cmatrix_float* v = new_cmatrix_float(1, 3);
cM(v, fill_rand, 1, 9, 1);
cmatrix_float* D2 = new_cmatrix_float(3, 3);
cM(D2, diag, v);   // D2 = diag(v[0], v[1], v[2])
```

### Concatenação horizontal

```c
// A: 3×2, B: 3×4 → C: 3×6
cmatrix_float* C = new_cmatrix_float(3, 2);
cM(C, hstack, A, B);
```

### Apply (aplicar função elemento a elemento)

```c
typedef struct { float a, b; } Params;

float linear(float x, void* p) {
    Params* ps = (Params*)p;
    return ps->a * x + ps->b;
}

Params ps = {2.0f, 1.0f};
cM(A, apply, linear, &ps);   // aᵢⱼ ← 2*aᵢⱼ + 1
```

### Reshape

```c
// Transforma uma matriz 2×6 em 3×4 (mesmo total de elementos)
cM(A, reshape, 3, 4);
```

### Trocar linhas e colunas

```c
cM(A, swapRows, 0, 2);   // troca linha 0 com linha 2
cM(A, swapCols, 1, 3);   // troca coluna 1 com coluna 3
```

### REF e RREF

```c
cmatrix_float* R  = new_cmatrix_float(3, 3);
cmatrix_float* RR = new_cmatrix_float(3, 3);

cM(R,  ref,  A);   // forma escalonada
cM(RR, rref, A);   // forma escalonada reduzida
```

### Inversa

```c
cmatrix_float* Ainv = new_cmatrix_float(3, 3);
int ok = cM(Ainv, inv, A);
if (!ok) printf("A nao e invertivel\n");
```

### Determinante

```c
float d = cM(A, det);
printf("det(A) = %.4f\n", d);
// O resultado fica em cache: chamar det() novamente sem modificar A é O(1)
```

---

## Liberação de memória

```c
cM(A, free);
```

Para liberar vários de uma vez:

```c
void free_all(cmatrix_float* ms[], int n) {
    for (int i = 0; i < n; i++)
        if (ms[i]) { cM(ms[i], free); ms[i] = NULL; }
}

cmatrix_float* mats[] = {A, B, C};
free_all(mats, 3);
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
