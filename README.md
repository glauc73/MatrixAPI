# MatrixC — Biblioteca Genérica de Matrizes em C Puro

> Álgebra linear com estilo orientado a objetos, sem depender de C++.

[![Linguagem](https://img.shields.io/badge/Linguagem-C11-blue.svg)](https://en.cppreference.com/w/c/11)
[![Header-only](https://img.shields.io/badge/Header--only-sim-green.svg)]()
[![Licença](https://img.shields.io/badge/Licença-MIT-yellow.svg)]()

---

## Sumário

- [Sobre o projeto](#sobre-o-projeto)
- [Nome sugerido](#nome-sugerido)
- [Funcionalidades](#funcionalidades)
- [Estrutura de arquivos](#estrutura-de-arquivos)
- [Pré-requisitos](#pré-requisitos)
- [Instalação](#instalação)
- [Compilação](#compilação)
  - [Linux / macOS](#linux--macos)
  - [Windows (MinGW / MSYS2)](#windows-mingw--msys2)
  - [Usando os scripts automáticos](#usando-os-scripts-automáticos)
- [Como usar a biblioteca](#como-usar-a-biblioteca)
  - [Incluindo no projeto](#incluindo-no-projeto)
  - [Criando e liberando matrizes](#criando-e-liberando-matrizes)
  - [Preenchimento](#preenchimento)
  - [Impressão](#impressão)
  - [Operações aritméticas](#operações-aritméticas)
  - [Operações de álgebra linear](#operações-de-álgebra-linear)
  - [Liberação de memória](#liberação-de-memória)
- [Exemplo completo](#exemplo-completo)
- [Referência da API](#referência-da-api)
- [FAQ](#faq)
- [Contribuindo](#contribuindo)

---

## Sobre o projeto

**MatrixC** é uma biblioteca *header-only* escrita em **C11 puro** que oferece:

- Um tipo abstrato de dado (TAD) `matrix_<tipo>` genérico para `int`, `float` e `double`
- Interface orientada a objetos usando a macro `M(matriz, método, ...)`, sem C++
- Algoritmos de álgebra linear clássicos implementados com eliminação gaussiana
- Zero dependências externas além da libc padrão

É ideal para projetos acadêmicos, sistemas embarcados e qualquer situação em que C++ não esteja disponível mas você ainda precise de operações matriciais robustas.

---

## Nome sugerido

O nome original **MatrixAPI** descreve o que a biblioteca é, mas não o que ela faz nem em qual linguagem. Sugestões melhores:

| Nome | Justificativa |
|------|--------------|
| **MatrixC** ✅ | Direto: matriz em C. Fácil de pesquisar, sem ambiguidade. |
| **cmatrix** | Segue convenção Unix (`cjson`, `curl`). Minúsculo, limpo. |
| **gmatrix** | De *generic matrix*. Destaca o diferencial genérico. |
| **linalgc** | Reforça que é álgebra linear em C. |

**Recomendação: `MatrixC`** — simples, memorável e imediatamente claro para quem pesquisa "matrix library C".

---

## Funcionalidades

| Operação | Descrição |
|----------|-----------|
| `new_matrix_<T>` | Aloca uma matriz M×N do tipo T |
| `fill_rand` | Preenche com valores aleatórios em um intervalo |
| `print` | Exibe a matriz formatada no terminal |
| `transp` | Transposta |
| `mult` | Multiplicação de matrizes |
| `diag` | Extrai / monta a diagonal principal |
| `REF` | Forma escalonada por linhas (*Row Echelon Form*) |
| `RREF` | Forma escalonada reduzida (*Reduced REF*) via eliminação gaussiana |
| `inv` | Matriz inversa |
| `det` | Determinante (eliminação gaussiana, O(n³)) |
| `free` | Libera a memória alocada |

---

## Estrutura de arquivos

```
MatrixC/
├── include/
│   ├── genericsBase.h    # Macros genéricas (infraestrutura de tipos)
│   ├── genericMatrix.h   # Implementação do TAD e das funções genéricas
│   └── Matrix.h          # Ponto de entrada único — inclua este no seu código
├── example.c             # Demonstração de uso
├── in.txt                # Entrada de exemplo para testar via stdin
├── build.sh              # Script de compilação (Linux/macOS)
├── build.bat             # Script de compilação (Windows)
└── README.md
```

> **Regra de ouro:** inclua apenas `Matrix.h` no seu código. Os outros headers são dependências internas.

---

## Pré-requisitos

### Linux / macOS
- GCC 7+ ou Clang 6+ com suporte a C11
- `make` (opcional, para projetos maiores)

```bash
# Ubuntu / Debian
sudo apt install gcc

# Fedora / RHEL
sudo dnf install gcc

# macOS (via Homebrew)
brew install gcc
```

### Windows
- [MinGW-w64](https://www.mingw-w64.org/) — recomendado via [MSYS2](https://www.msys2.org/)
- Ou [TDM-GCC](https://jmeubank.github.io/tdm-gcc/)

```powershell
# Via MSYS2 (dentro do terminal MSYS2)
pacman -S mingw-w64-x86_64-gcc
```

---

## Instalação

A biblioteca é *header-only*: basta copiar a pasta `include/` para dentro do seu projeto.

```bash
# Clone o repositório
git clone https://github.com/glauc73/MatrixAPI.git
cd MatrixAPI
```

---

## Compilação

### Linux / macOS

```bash
gcc -std=c11 example.c -Iinclude -lm -o main
./main < in.txt
```

> A flag `-lm` é necessária se seu código usar funções de `<math.h>`.

### Windows (MinGW / MSYS2)

```cmd
gcc -std=c11 example.c -Iinclude -lm -o main.exe
main.exe < in.txt
```

### Usando os scripts automáticos

Os scripts detectam o compilador disponível, compilam e executam o exemplo automaticamente.

**Linux / macOS:**
```bash
chmod +x build.sh
./build.sh
```

**Windows:**
```cmd
build.bat
```

---

## Como usar a biblioteca

### Incluindo no projeto

```c
#include "include/Matrix.h"
```

Isso é tudo. Não há arquivos `.c` para compilar separadamente.

### Criando e liberando matrizes

A convenção de nomenclatura é `new_matrix_<tipo>(linhas, colunas)`.

```c
// Tipos disponíveis: int, float, double
matrix_float* A  = new_matrix_float(3, 3);   // 3×3 float
matrix_int*   B  = new_matrix_int(2, 4);     // 2×4 int
matrix_double* C = new_matrix_double(4, 4);  // 4×4 double
```

### Preenchimento

```c
// Preenche com valores aleatórios entre -10 e 10
// Último parâmetro: 1 = inteiros, 0 = floats
M(A, fill_rand, -10, 10, 1);

// Preenche com zeros (comportamento padrão de new_matrix)
// As matrizes já são alocadas zeradas
```

### Impressão

```c
M(A, print);
```

Saída típica:
```
[ -3.00   7.00   1.00 ]
[  5.00  -2.00   9.00 ]
[  0.00   4.00  -6.00 ]
```

### Operações aritméticas

```c
matrix_float* resultado = new_matrix_float(3, 3);

// Multiplicação: resultado = B * C
M(resultado, mult, B, C);

// Transposta de A (resultado é armazenado no próprio A)
M(A, transp, A);
```

### Operações de álgebra linear

```c
// Forma escalonada (REF)
M(A, REF, A);

// Forma escalonada reduzida (RREF)
M(A, RREF, A);

// Inversa
matrix_float* inv_A = new_matrix_float(3, 3);
M(inv_A, inv, A);

// Determinante (retorna float/double/int)
float d = M(A, det);
printf("det(A) = %.4f\n", d);

// Diagonal
matrix_float* diag_A = new_matrix_float(3, 3);
M(diag_A, diag, A);
```

### Liberação de memória

Sempre libere a memória após o uso:

```c
M(A, free);
M(B, free);
M(C, free);
```

Ou use o padrão helper do exemplo:

```c
void free_all(matrix_float* mtrzs[], int n) {
    for (int i = 0; i < n; i++) {
        if (mtrzs[i] != NULL) {
            M(mtrzs[i], free);
            mtrzs[i] = NULL;
        }
    }
}

matrix_float* mtzs[] = {A, B, C};
free_all(mtzs, 3);
```

---

## Exemplo completo

O arquivo `example.c` demonstra o fluxo completo:

```c
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "include/Matrix.h"

#define true 1

void free_all(matrix_float* mtrzs[], int n) {
    if (!mtrzs) return;
    for (int i = 0; i < n; i++) {
        if (mtrzs[i] != NULL) {
            M(mtrzs[i], free);
            mtrzs[i] = NULL;
        }
    }
}

int main() {
    srand(time(NULL));

    // Cria matrizes
    matrix_float* A = new_matrix_float(3, 6);
    matrix_float* b = new_matrix_float(3, 1);
    matrix_float* x = new_matrix_float(3, 1);

    // Preenche x e b com valores aleatórios inteiros entre -10 e 10
    M(x, fill_rand, -10, 10, true);
    M(b, fill_rand, -10, 10, true);

    // Imprime b e x transposta
    M(b, print);

    {
        matrix_float* xt = new_matrix_float(3, 1);
        M(xt, transp, x);
        M(xt, print);

        // A = b * xᵀ  (produto externo)
        M(A, mult, b, xt);
        M(xt, free);
    }

    // Imprime A e sua diagonal
    M(A, print);
    M(A, diag, A);
    M(A, print);

    // Libera tudo
    matrix_float* mtzs[] = {A, b, x};
    free_all(mtzs, 3);

    return 0;
}
```

**Compilar e rodar:**
```bash
gcc -std=c11 example.c -Iinclude -lm -o main
./main < in.txt
```

---

## Referência da API

### Tipos disponíveis

```c
matrix_int    // Matriz de inteiros
matrix_float  // Matriz de floats
matrix_double // Matriz de doubles
```

### Construtores / Destrutores

| Chamada | Descrição |
|---------|-----------|
| `new_matrix_float(m, n)` | Aloca matriz float M×N |
| `new_matrix_int(m, n)` | Aloca matriz int M×N |
| `new_matrix_double(m, n)` | Aloca matriz double M×N |
| `M(mat, free)` | Libera a memória |

### Métodos (via macro `M`)

| Método | Parâmetros extras | Descrição |
|--------|-------------------|-----------|
| `print` | — | Imprime a matriz |
| `fill_rand` | `min, max, inteiros` | Preenche aleatoriamente |
| `transp` | `src` | Transposta de `src` → `dest` |
| `mult` | `A, B` | Multiplicação: `dest = A * B` |
| `diag` | `src` | Diagonal principal de `src` → `dest` |
| `REF` | `src` | Forma escalonada de `src` |
| `RREF` | `src` | Forma escalonada reduzida de `src` |
| `inv` | `src` | Inversa de `src` → `dest` |
| `det` | — | Retorna o determinante |

### A macro `M`

A macro `M(matriz, metodo, ...)` é o coração da interface OO:

```c
// Equivalente a: matriz->metodo(matriz, ...)
M(A, print);
M(A, mult, B, C);   // A = B * C
```

---

## FAQ

**P: Posso usar tipos além de `int`, `float` e `double`?**
R: Sim. O sistema de macros em `genericsBase.h` permite adicionar novos tipos instanciando as funções genéricas de `genericMatrix.h` para o novo tipo. Veja os exemplos existentes como modelo.

**P: A biblioteca é thread-safe?**
R: As operações em matrizes diferentes são independentes. Operações na *mesma* matriz a partir de múltiplas threads precisam de sincronização externa.

**P: Preciso linkar alguma biblioteca extra?**
R: Apenas `-lm` (libmath) se seu código usar `<math.h>` diretamente. A biblioteca em si depende apenas de `stdlib.h` e `string.h`.

**P: Qual padrão C é necessário?**
R: C11 (`-std=c11`). Versões mais antigas podem funcionar, mas não são garantidas.

---

## Contribuindo

1. Faça um fork do repositório
2. Crie um branch: `git checkout -b minha-feature`
3. Commit suas mudanças: `git commit -m 'Adiciona feature X'`
4. Push: `git push origin minha-feature`
5. Abra um Pull Request

---

*Desenvolvido por [glauc73](https://github.com/glauc73)*
