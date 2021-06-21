#include "odeSolver.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Constantes utilizadas
#define N1 5
#define N2 10
#define M 7
#define L 6
#define W 8
#define M_PI 3.14159265359
#define M_PI_2 1.57079632679

// Declaração das struckts equações dos exercicios (terminados em 1 = usa N1;  terminados em 2 = usa N2)
Edo *a1;
Edo *a2;
Edo *c1;
Edo *c2;

Edp *b1;
Edp *b2;
Edp *d1;
Edp *d2;

// Definição dos termos independentes(da função) da equação
double zero(double x) { return 0; }
double zeroxy(double x, double y) { return 0; }

double const1(double x) { return 1; }
double const1xy(double x, double y) { return 1; }
double constsub1xy(double x, double y) { return -1; }
double rA(double x) { return (6 * x) - (0.5 * x * x); }
double rB(double x, double y)
{
    double s = sin(x);
    return s * s;
}
// double rC(double x) { return 0; }// mesmo que zero()
double rD(double x, double y) { return (-cos(x + y) + cos(x - y)); }

double *Y;
double **Z;

double res;
int time = 0;
int main()
{

    Y = calloc(N1, sizeof(double));
    a1 = alocaEdo(N1, 0, 12, 0, 0, zero, zero, rA);
    printf("\e[92m\n***** item (a)\e[39m: n = %d, H = %lf\n", N1, (a1->b - a1->a) / (a1->n + 1));
    time = gaussSeidelEDO(a1, Y);
    printY(&Y, N1);
    printf("Norma L2: %lf, Tempo: %d ms\n", normaL2Edo(a1, Y, &res), time);
    free(Y);
    free(a1);

    Y = calloc(N2, sizeof(double));
    a2 = alocaEdo(N2, 0, 12, 0, 0, zero, zero, rA);
    printf("\e[92m\n***** item (a)\e[39m: n = %d, H = %lf\n", N2, (a2->b - a2->a) / (a2->n + 1));
    time = gaussSeidelEDO(a2, Y);
    printY(&Y, N2);
    printf("Norma L2: %lf, Tempo: %d ms\n", normaL2Edo(a2, Y, &res), time);
    free(Y);
    free(a2);

    Y = calloc(N1, sizeof(double));
    c1 = alocaEdo(N1, 0, 1, 0, 1, zero, const1, zero);
    printf("\e[92m\n***** item (c)\e[39m: n = %d, H = %lf\n", N1, (c1->b - c1->a) / (c1->n + 1));
    time = gaussSeidelEDO(c1, Y);
    printY(&Y, N1);
    printf("Norma L2: %lf, Tempo: %d ms\n", normaL2Edo(c1, Y, &res), time);
    free(Y);
    free(c1);

    Y = calloc(N2, sizeof(double));
    c2 = alocaEdo(N2, 0, 1, 0, 1, zero, const1, zero);
    printf("\e[92m\n***** item (c)\e[39m: n = %d, H = %lf\n", N2, (c2->b - c2->a) / (c2->n + 1));
    time = gaussSeidelEDO(c2, Y);
    printY(&Y, N2);
    printf("Norma L2: %lf, Tempo: %d ms\n", normaL2Edo(c2, Y, &res), time);
    free(Y);
    free(c2);

    b1 = alocaEdp(N1, M, 0, L, 0, W, 20, 45, 0, 100, constsub1xy, const1xy, const1xy, rB);

    Z = alocaZ(N1, M, 20, 45, 0, 100);
    dBug("aquio3");
    printf("\e[92m\n***** item (b)\e[39m: L = %d, W = %d, n = %d, m = %d,  Hx = %lf,  Hy = %lf\n", L, W, N1, M, (b1->bx - b1->ax) / (b1->n + 1), (b1->by - b1->ay) / (b1->m + 1));
    dBug("aquio4");

    time = gaussSeidelEDP(b1, Z);
    dBug("aquio5");

    printZ(Z, N1, M);
    printf("Norma L2: %lf, Tempo: %d ms\n", normaL2Edp(b1, Z, &res), time);
    freeZ(Z, N1, M);
    free(b1);

    b2 = alocaEdp(N2, M, 0, L, 0, W, 20, 45, 0, 100, constsub1xy, const1xy, const1xy, rB);
    Z = alocaZ(N2, M, 20, 45, 0, 100);
    printf("\e[92m\n***** item (b)\e[39m: L = %d, W = %d, n = %d, m = %d,  Hx = %lf,  Hy = %lf\n", L, W, N2, M, (b2->bx - b2->ax) / (b2->n + 1), (b2->by - b2->ay) / (b2->m + 1));
    time = gaussSeidelEDP(b2, Z);
    printZ(Z, N2, M);
    printf("Norma L2: %lf, Tempo: %d ms\n", normaL2Edp(b2, Z, &res), time);
    freeZ(Z, N2, M);
    free(b2);

    d1 = alocaEdp(N1, M, 0, M_PI, 0, M_PI_2, 0, 0, 0, 0, constsub1xy, const1xy, zeroxy, rD);
    Z = alocaZ(N1, M, 0, 45, 0, 100);
    printf("\e[92m\n***** item (d)\e[39m: L = %d, W = %d, n = %d, m = %d,  Hx = %lf,  Hy = %lf\n", L, W, N1, M, (d1->bx - d1->ax) / (d1->n + 1), (d1->by - d1->ay) / (d1->m + 1));
    time = gaussSeidelEDP(d1, Z);
    printZ(Z, N1, M);
    printf("Norma L2: %lf, Tempo: %d ms\n", normaL2Edp(d1, Z, &res), time);
    freeZ(Z, N1, M);
    free(d1);

    d2 = alocaEdp(N2, M, 0, M_PI, 0, M_PI_2, 0, 0, 0, 0, constsub1xy, const1xy, zeroxy, rD);
    Z = alocaZ(N2, M, 20, 45, 0, 100);

    double hx = (d2->bx - d2->ax) / (d2->n + 1), hy = (d2->by - d2->ay) / (d2->m + 1);

    for (int i = 0; i < N2; i++)
    {
        Z[i][0] = cos(i * hx);
        Z[i][M - 1] = 0;
    }
    for (int j = 01; j < M - 1; j++)
    {
        Z[0][j] = cos(j * hy);
        Z[N2 - 1][j] = -cos(j * hy);
    }
    printf("\e[92m\n***** item (d)\e[39m: L = %d, W = %d, n = %d, m = %d,  Hx = %lf,  Hy = %lf\n", L, W, N2, M, hx, hy);
    time = gaussSeidelEDP(d2, Z);
    printZ(Z, N2, M);
    printf("Norma L2: %lf, Tempo: %d ms\n", normaL2Edp(d2, Z, &res), time);
    freeZ(Z, N2, M);
    free(d2);

    // printf("Teste C1\n");
    // b1 = alocaEdo(N2, 0, 1, 0, 1, zero, const1, zero);
    // gaussSeidelEDO(c2, Y2);
    // printY(&Y2, N2);
    // b1 = alocaEdp(N1, M, 0, L, 0, W, 20, 45, 0, 100, constsub1xy, const1xy, const1xy, rB);
    // Z1 = alocaZ(N1, M, 20, 45, 0, 100);
    // gaussSeidelEDP(b1, Z1);
    // printZ(Z1, N1, M);
    // b1= alocaEdp()
}