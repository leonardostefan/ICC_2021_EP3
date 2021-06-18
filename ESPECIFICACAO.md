# Exercício EP-02
## Objetivo
O objetivo deste exercício é implementar uma biblioteca para a resolução de sistemas lineares.

## Especificação
As funções a serem implementadas estão declaradas no arquivo SistemasLineares.h que não pode ser alterado.

- Método da Eliminação de Gauss (com pivoteamento parcial )
- Método de Jacobi
- Método de Gauss-Seidel
    - Valor inicial: o vetor x deve ter valor inicial nulo
    - erro: erro máximo da norma em x
    - Critério de parada é a norma máxima do erro absoluto aproximado em x:


- Método de Refinamento

Além das funções, o arquivo SistemasLineares.h também define:

- MAXIT: Número máximo de iterações antes de abortar os métodos iterativos;
- real_t: tipo de dados que pode ser alterado em tempo de compilação para float ou double;
- struct SistLinear_t: estrutura para definição/armazenamento dos coeficientes e termos independentes do Sistema Linear, bem como sua dimensão e valor para critério de parada.
A alocação e liberação de memória são feitas respectivamente pelas funções alocaSistLinear(...) e  liberaSistLinear(...), definidas no arquivo SistemasLineares.c

O arquivo sistemas.dat contém o conjunto Sistemas Lineares (SL) que deverão ser usados neste exercício. Este arquivo contém diversos blocos   separados por uma linha em branco. Cada bloco tem a ordem do sistema linear na 1a linha, o valor a ser usado no critério de parada na 2ª linha, e nas linhas restantes do bloco estão a matriz de coeficientes, seguida do vetor de termos independentes.

Uso da biblioteca
Todas funções da biblioteca devem "falhar graciosamente". Em caso de erro, o retorno das funções deve ser negativo, e uma mensagem de erro deve ser impressa na saída de erros (stderr), sem no entanto encerrar a execução do programa. CUIDADO: Nenhuma função deve imprimir na saída padrão, apenas seu programa principal deve imprimir resultados.

Tarefa
O objetivo é comparar o comportamento dos métodos em diferentes sistemas lineares. Deve ser feito um programa de nome labSisLin que, usando as funções da biblioteca indicada na especificação e tendo como entrada cada SL no arquivo sistemas.dat, faça:

Encontrar a solução do SL com os métodos de Eliminação de Gauss (com pivotamento), Jacobi e Gauss-Seidel, indicando  (para os dois últimos) o número de iterações até solução;
Medir o tempo (em milisegundos) de cada método até a solução utilizando a função timestamp() (vide arquivos utils.*);
Calcular a norma L2 (ou norma euclidiana) do resíduo de cada uma das soluções;
Aplicar o método de Refinamento à solução obtida em cada um dos outros métodos caso a norma L2 do resíduo seja maior que 5.0;
Exemplo de execução do programa:

Entrada:
4
0.05
4 -1 0 -1
-1 4 -1 0
0 -1 4 -1
-1 0 -1 4
3 6 1 12

4
0.001
31.0153    23.8669     17.7260    4.4339
8.7305      27.9141     26.8869     1.6694
22.6307     6.4646      26.4602     21.2894
2.4314     17.8153      19.4762     26.2457
9.7445      29.9868      18.3525    15.0411


Saída:

***** Sistema 1 --> n = 4, erro: 0.05
===> Eliminação Gauss: 0.0009765625 ms
     --> X: 2.416667 2.583333 1.916667 4.083333
     --> Norma L2 do residuo: 4.768372e-07

===> Jacobi: 0.0009765625 ms --> 7 iteracoes
     --> X: 2.390625 2.566406 1.890625 4.066406
     --> Norma L2 do residuo: 0.1018625

===> Gauss-Seidel: 0 ms --> 5 iteracoes
     --> X: 2.405941 2.57554 1.910909 4.079212
     --> Norma L2 do residuo: 0.03605055

***** Sistema 2 --> n = 4, erro: 0.001
===> Eliminação Gauss: 0 ms
     --> X: -0.3530174 -0.006185221 1.256357 -0.3223174
     --> Norma L2 do residuo: 9.536743e-07

===> Jacobi: 0.004150391 ms --> 50 iteracoes
     --> X: -2.821267e+09 -2.740078e+09 -3.423687e+09 -2.932857e+09
     --> Norma L2 do residuo: 4.305108e+11

===> Refinamento: 0.0009765625 ms --> 2 iteracoes
     --> X: -0.3530048 -0.006214148 1.256385 -0.3223572
     --> Norma L2 do residuo: 0.0009798451

===> Gauss-Seidel: 0.004882812 ms --> 50 iteracoes
     --> X: -2.129541 -10.23473 8.774721 1.206116
     --> Norma L2 do residuo: 224.2306

===> Refinamento: 0.0009765625 ms --> 1 iteracoes
     --> X: -0.3530164 -0.006184578 1.256355 -0.3223165
      --> Norma L2 do residuo: 2.778778e-05

Acesso aos arquivos

Os arquivos para este exercício encontram-se no git em

https://gitlab.c3sl.ufpr.br/nicolui/ci1164_2020-sislin.git
ou via linha de comando:

git clone git@gitlab.c3sl.ufpr.br:nicolui/ci1164_2020-sislin.git
Status de envio