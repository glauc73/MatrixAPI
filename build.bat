@echo off

gcc -std=c11 %1 -Iinclude -lm -o main.exe

if errorlevel 1 (
    echo Erro: compilacao falhou.
    exit /b 1
)

echo Compilado com sucesso! Executando...
main.exe
