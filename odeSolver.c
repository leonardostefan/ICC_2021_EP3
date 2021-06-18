#include "odeSolver.h"

void gaussSeidelEDO(Edo *edo, double *Y)
{
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

Edp *alocaEdp(unsigned int n, double a, double b, double ya, double yb, double (*px)(double), double (*pxx)(double), double (*pxy)(double), double (*pyy)(double), double (*py)(double))
{
    // Edp *alocado = calloc(1, sizeof(Edp));
    // alocado->n = n;
    // alocado->a = a;
    // alocado->b = b;
    // alocado->ya = ya;
    // alocado->yb = yb;
    // alocado->pxx = pxx;
    // alocado->px = px;
    // alocado->pxy = pxy;
    // alocado->pyy = pyy;
    // alocado->py = py;
    // return alocado;
}