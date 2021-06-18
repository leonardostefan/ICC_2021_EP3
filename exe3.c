#include "odeSolver.h"
#include <stdio.h>
#include <stdlib.h>

#define N1 5
#define N2 10
#define M 7

Edo *a1;
Edo *a2;
Edp *b1;
Edp *b2;
Edo *c1;
Edo *c2;
Edp *d1;
Edp *d2;
Edo *e1;
Edo *e2;
void dBug(const char *_msg_)
{
  printf("\e[33mDEBUG:\e[39m");
  printf(_msg_);
  printf("\n");
  fflush(NULL);
}
double zero(double x)
{
    return 0;
}
double const1(double x)
{
    return 1;
}
double rA(double x)
{
    return (6 * x) - (0.5 * x * x);
}
double *Y1;
double *Y2;

void printY(double **Y, int n)
{

    for (int i = 0; i < n; i++)
    {
        printf("Y[%d]: %lf", i, (*Y)[i]);
    }
    printf("\n");
}

int main()
{
    Y1= calloc (N1, sizeof(double));
    Y2= calloc (N2, sizeof(double));
    printf("Teste A1\n");
    dBug("1");
    a1 = alocaEdo(N1, 0, 12, 0, 0, zero, zero, rA);
    dBug("2");
    gaussSeidelEDO(a1, Y1);
    dBug("3");
    printY(&Y1, N1);
    dBug("4");
    printf("Teste A2\n");
    a2 = alocaEdo(N2, 0, 12, 0, 0, zero, zero, rA);
    gaussSeidelEDO(a2, Y2);
    printY(&Y2, N2);

    printf("Teste C1\n");

    c1 = alocaEdo(N1, 0, 1, 0, 1, zero, const1, zero);
    gaussSeidelEDO(c1, Y1);
    printY(&Y1, N1);

    printf("Teste C1\n");
    c2 = alocaEdo(N2, 0, 1, 0, 1, zero, const1, zero);
    gaussSeidelEDO(c2, Y2);
    printY(&Y2, N2);

    // b1= alocaEdp()
}