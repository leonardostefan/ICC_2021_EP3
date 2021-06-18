# Exercício EP-03
## Objetivo

O objetivo deste exercício é implementar programas para a solução de Equações Diferenciais Ordinárias e Parciais.

## Especificação
Implementar programas para resolver as equações diferenciais abaixo, considerando as condições de contorno indicadas:



Em todos os casos, gere resultados para n = 5 e 10 (malha x). Nos casos (b) e (d) considere também m =  3 (malha y).

Seu programa deve, para cada caso, exibir a matriz aumentada do sistema linear resultante (mariz completa OU diagonais e termos independentes apenas), sua solução (calculada com o método de Gauss-Seidel), a norma L2 do resíduo e o tempo gasto para a solução.

Use como critério de parada do método como sendo o máximo de 50 iterações.

O código que implementa Gauss-Seidel deve estar adaptado para matrizes k−diagonais, considerando que são implementadas com vetores representando apenas as diagonais e termos independentes.

## O que deve ser entregue
Devem ser entregues os códigos-fonte do programa (em linguagem C), makefile e outros arquivos que sejam necessários à execução do programa (se houver)  usados pelo aluno. Todos estes itens devem ser compactados em um arquivo do tipo .tar , .tgz ou .zip.

Em particular, o arquivo Makefile deve ter os targets de limpeza, como no Exercício EP-02.
