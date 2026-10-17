#!/bin/bash

gcc -std=c11 -g $1 -Iinclude -lm -o main

if [ $? -ne 0 ]; then
    echo "Erro: compilacao falhou."
    exit 1
fi

echo "Compilado com sucesso!"
./main
