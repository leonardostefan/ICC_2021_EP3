#include <stdlib.h>
#include "odeSolver.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

void dBug(const char *_msg_)
{
    printf("\e[33mDEBUG:\e[39m");
    printf(_msg_);
    printf("\n");
    fflush(NULL);
}
int gaussSeidelEDO(Edo *edo, double *Y)
{
    double initTime = timestamp();
    int n = edo->n, k, i;
    double h, xi, bi, yi, d, di, ds;

    h = (edo->b - edo->a) / (n + 1); // Largura do passo da malha
    // Constantes baseadas na largura do passo da malha
    double hpot2 = h * h, hdiv2 = h / 2.0;

    for (k = 0; k < 50; ++k)
    { //9 +10*n +10 flops por iteração, antes eram 24*n

        // (i == 0)
        xi = edo->a + h;             // valor xi da malha: 1 FLOP
        bi = hpot2 * edo->r(xi);     // termo independente: 1 FLOP
        di = 1 - hdiv2 * edo->p(xi); // diagonal inferior: 2 FLOP
        d = -2 + hpot2 * edo->q(xi); // diagonal principal: 2 FLOP
        ds = 1 + hdiv2 * edo->p(xi); // diagonal superior: 2 FLOP
        bi -= ds * Y[1] + edo->ya * (1 - hdiv2 * edo->p(edo->a + h));
        Y[0] = bi / d; // Calcula incógnita: 1 FLOP

        // 10 FLOP por iteração do método
        for (i = 1; i < n - 1; ++i)
        {                                        // Para cada equação do SL
            xi = edo->a + (i + 1) * h;           // valor xi da malha: 1 FLOP
            bi = hpot2 * edo->r(xi);             // termo independente: 2 FLOP
            di = 1 - hdiv2 * edo->p(xi);         // diagonal inferior: 2 FLOP
            d = -2 + hpot2 * edo->q(xi);         // diagonal principal: 2 FLOP
            ds = 1 + hdiv2 * edo->p(xi);         // diagonal superior: 2 FLOP
            bi -= ds * Y[i + 1] + di * Y[i - 1]; // 4 FLOP
            Y[i] = bi / d;                       // Calcula incógnita: 1 FLOP
        }
        // (i == n - 1)
        xi = edo->a + (n - 1 + 1) * h;               // valor xi da malha: 1 FLOP
        bi = hpot2 * edo->r(xi);                     // termo independente: 2 FLOP
        di = 1 - hdiv2 * edo->p(xi);                 // diagonal inferior: 2 FLOP
        d = -2 + hpot2 * edo->q(xi);                 // diagonal principal: 2 FLOP
        ds = 1 + hdiv2 * edo->p(xi);                 // diagonal superior: 2 FLOP
        bi -= ds * Y[n - 1 + 1] + di * Y[n - 1 - 1]; // 4 FLOP
        Y[n - 1] = bi / d;                           // Calcula incógnita: 1 FLOP
        bi -= di * Y[n - 1 - 1] + edo->yb * (1 + hdiv2 * edo->p(edo->b - h));
    }
    return timestamp() - initTime;
}

int gaussSeidelEDP(Edp *edp, double **Z)
{
    double initTime = timestamp();

    int n = edp->n, k, i;
    int m = edp->m;
    double hx, hy, xi, yj, bij;

    double d0, d1, d2, d3, d4;

    double d0z, d1z, d2z, d3z, d4z;

    hx = (edp->bx - edp->ax) / (n + 1); // Largura do passo da malha
    hy = (edp->by - edp->ay) / (m + 1); // Largura do passo da malha

    double hxpot2 = hx * hx, hxdiv2 = hx / 2.0;
    double hypot2 = hy * hy, hydiv2 = hy / 2.0;

    double px, py;
    for (k = 0; k < MAXIT; ++k)
    { // 23 FLOP por iteração do método
        for (i = 1; i < n - 1; ++i)
        {
            xi = edp->ax + (i + 1) * hx;
            for (int j = 1; j < m - 1; j++)
            {
                yj = edp->ay + (j + 1) * hy;
                px = edp->px(xi, yj);
                py = edp->py(xi, yj);

                d0 = hxpot2 * py;
                d1 = hypot2 * px;
                d2 = -2 * (hxpot2 * py + hypot2 * px) + edp->p(xi, yj);
                d3 = hypot2 * px;
                d4 = hxpot2 * py;
                bij = hxpot2 * hypot2 * edp->r(xi, yj);

                if ((i == 1))
                {
                    // d0z = d0 * edp->zxa;
                    d0z = 0;
                }
                else if ((j == 1))
                {
                    // d0z = d0 * edp->zya;
                    d0z = 0;
                }
                else
                    d0z = d0 * Z[i - 2][j - 2];
                if (i == n - 2)
                {
                    // d4z = d4 * edp->zxb;
                    d4z = 0;
                }
                else if (j == n - 2)
                {
                    // d4z = d4 * edp->zyb;
                    d4z = 0;
                }
                else
                    d4z = d4 * Z[i + 2][j + 2];

                bij -= d0z + d1 * Z[i - 1][j - 1] + d3 * Z[i + 1][j + 1] + d4z;
                Z[i][j] = bij / d2; // Calcula incógnita: 1 FLOP}
            }
        }
    }
    return timestamp() - initTime;
}

Edo *alocaEdo(unsigned int n, double a, double b, double ya, double yb, double (*p)(double), double (*q)(double), double (*r)(double))
{
    Edo *alocado = calloc(1, sizeof(Edo));
    alocado->n = n;
    alocado->a = a;
    alocado->b = b;
    alocado->ya = ya;
    alocado->yb = yb;
    alocado->p = p;
    alocado->q = q;
    alocado->r = r;
    return alocado;
}

Edp *alocaEdp(unsigned int n, unsigned int m, double ax, double bx, double ay, double by, double zxa, double zxb, double zya, double zyb, double (*p)(double, double), double (*px)(double, double), double (*py)(double, double), double (*r)(double, double))
{
    Edp *alocado = calloc(1, sizeof(Edp));
    alocado->n = n;
    alocado->m = m;
    alocado->ax = ax;
    alocado->bx = bx;
    alocado->ay = ay;
    alocado->by = by;
    alocado->zxa = zxa;
    alocado->zxb = zxb;
    alocado->zya = zya;
    alocado->zyb = zyb;
    alocado->p = p;
    alocado->px = px;
    alocado->py = py;
    alocado->r = r;
    return alocado;
}

double **alocaZ(int n, int m, double zxa, double zxb, double zya, double zyb)
{
    double **z;
    z = calloc(n, sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        z[i] = calloc(m, sizeof(double));
    }

    for (int i = 0; i < n; i++)
    {
        z[i][0] = zya;
        z[i][m - 1] = zyb;
    }
    for (int j = 01; j < m - 1; j++)
    {
        z[0][j] = zxa;
        z[n - 1][j] = zxb;
    }

    return z;
}

void printZ(double **z, int n, int m)
{
    printf("Malha de resultados \e[1mcom contorno\e[0m utilizado:\n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            printf("%.3lf\t", z[i][j]);
        }
        printf("\n");
    }
}

void printY(double **Y, int n)
{
    printf("Y: ");
    for (int i = 0; i < n; i++)
    {
        printf("%.3lf; ", (*Y)[i]);
    }
    printf("\n");
}

void freeZ(double **Z, int n, int m)
{
    for (int i = n - 1; i > 0; i--)
    {
        free(Z[i]);
    }
    free(Z);
}

double normaL2Edo(Edo *SL, double *Y, double *res)
{
    // *res = 0;
    // double bxp = 0;
    // double prod = 0;
    // double hx, xi = SL->a, bi;
    // int n = SL->n;
    // hx = (SL->b - SL->a) / (n + 1); // Largura do passo da malha

    // prod += SL->p(xi) * Y[0];
    // prod += SL->p(xi) * Y[0];
    // bxp = SL->r(xi) - prod;
    // (*res) += bxp * bxp;

    // for (int i = 1; i < n - 1; i++, prod = 0, xi += hx)
    // {
    //     prod += (1) * Y[i - 1];
    //     prod += SL->p(xi) * Y[i];
    //     prod += SL->p(xi) * Y[i + 1];
    //     bxp = SL->r(xi) - prod;
    //     (*res) += bxp * bxp;
    // }
    // prod += (1) * Y[n - 2];
    // prod += SL->p(xi) * Y[n - 1];
    // bxp = SL->r(xi) - prod;
    // (*res) += bxp * bxp;

    // return (sqrtf(*res));
    return INFINITY;

}

double normaL2Edp(Edp *SL, double **Z, double *res)
{
    // *res = 0;
    // double bxp = 0;
    // double prod = 0;
    // double hx, hy, xi = SL->ax, yj = SL->ay, bij;
    // int n = SL->n;
    // int m = SL->m;

    // hx = (SL->bx - SL->ax) / (n + 1); // Largura do passo da malha

    // hy = (SL->by - SL->ay) / (m + 1); // Largura do passo da malha

    // for (int i = 0; i < n; i++, prod = 0, xi += hx)
    // {
    //     for (int j = 0; j < n; j++, yj += hy)
    //     {
    //         prod += (1) * Z[i - 1][j - 1];
    //         prod += SL->p(xi,yj) * Z[i][j - 1];
    //         prod += SL->p(xi,yj) * Z[i + 1][j - 1];
    //         bxp = SL->r(xi,yj) - prod;
    //         (*res) += bxp * bxp;
    //     }
    // }

    // return (sqrtf(*res));
    return INFINITY;
}
