**Include** guarda as funções e suas dependencias organizadas da seguinte maneira:

**genericsBase.h**
contem as dependencias necessarias para criar funções para varios tipos

**genericMatrix.h**
contem a implementação das funções genericas e do TAD matrix

**Matrix.h**
é o arquivo que deve ser incluido no seu codigo para usar as implementações.
nele possui o Guard do projeto e as definições dos tipos

**example_code** tem um exemplo de como usar a biblioteca, e deve ser compilado da seguinte maneira:

**gcc -std=c11 main.c -Iinclude -o main**

