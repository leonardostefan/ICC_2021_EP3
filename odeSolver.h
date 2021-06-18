#ifndef __ODESOLVER__
#define __ODESOLVER__

// Parâmetros para teste de convergência
#define MAXIT 50 // Número máximo de iterações em métodos iterativos

typedef float real_t;

// Matriz tridiagonal
typedef struct
{
    double *D, *Di, *Ds, *B;
    int n;
} SL_Tridiag;

// Equação Diferencial Ordinária para MDF
typedef struct
{
    int n;         // número de pontos internos na malha
    double a, b;   // intervalo
    double ya, yb; // condições contorno
    double (*p)(double), (*q)(double), (*r)(double);
} Edo;

// Equação Diferencial Ordinária para MDF
typedef struct
{
    int n; // número de pontos internos na malha

    double ax, bx; // intervalo
    double ay, by;
    double zya, zyb; // condições contorno
    double zxa, zxb; // condições contorno
    double (*px)(double), (*pxx)(double), (*pxy)(double), (*pyy)(double), (*py)(double);
} Edp;

// Alocaçao e desalocação de memória
Edo *alocaEdo(unsigned int n, double a, double b, double ya, double yb, double (*p)(double), double (*q)(double), double (*r)(double));
Edp *alocaEdp(unsigned int n, double a, double b, double ya, double yb, double (*p0)(double), double (*p1)(double), double (*p2)(double), double (*p3)(double), double (*p4)(double));

// Leitura e impressão de sistemas lineares
void prnSistLinear3(SL_Tridiag *SL);
void prnVetor(real_t *vet, unsigned int n);

// Retorna a normaL2 do resíduo. Parâmetro 'res' deve ter o resíduo.
real_t normaL2Residuo(SL_Tridiag *SL, real_t *x, real_t *res);

void gaussSeidelEDO(Edo *edo, double *Y);

void gaussSeidelEDP(Edp *edp, double *Z);

// Método de Refinamento. Valor inicial e resultado no parâmetro 'x'
int refinamento(SL_Tridiag *SL, real_t *x, double *tTotal);

#endif // __SISLINEAR_H__
