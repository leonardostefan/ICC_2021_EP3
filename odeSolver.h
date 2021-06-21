#ifndef __ODESOLVER__
#define __ODESOLVER__

// Parâmetros para teste de convergência
#define MAXIT 5000 // Número máximo de iterações em métodos iterativos

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
    int n; // número de pontos x internos na malha
    int m; // número de pontos y internos na malha

    double ax, bx; // intervalo
    double ay, by;
    double zxa, zxb, zya, zyb;
    double (*p)(double, double), (*px)(double, double), (*py)(double, double);
    double (*r)(double, double);
} Edp;

// Alocaçao e desalocação de memória
Edo *alocaEdo(unsigned int n, double a, double b, double ya, double yb, double (*p)(double), double (*q)(double), double (*r)(double));
Edp *alocaEdp(unsigned int n, unsigned int m, double ax, double bx, double ay, double by, double zxa, double zxb, double zya, double zyb, double (*p)(double, double), double (*px)(double, double), double (*py)(double, double), double (*r)(double, double));
double **alocaZ(int n, int m, double zxa, double zxb, double zya, double zyb);
void freeZ(double **Z, int n, int m);

// Leitura e impressão de sistemas lineares
// void prnSistLinear3(SL_Tridiag *SL);
void printZ(double **z, int n, int m);
void printY(double **Y, int n);

// Retorna a normaL2 do resíduo. Parâmetro 'res' deve ter o resíduo.
real_t normaL2Residuo(SL_Tridiag *SL, real_t *x, real_t *res);

int gaussSeidelEDO(Edo *edo, double *Y);
int gaussSeidelEDP(Edp *edp, double **Z);

// Método de Refinamento. Valor inicial e resultado no parâmetro 'x'
int refinamento(SL_Tridiag *SL, real_t *x, double *tTotal);

double normaL2Edo(Edo *SL, double *Y, double *res);
double normaL2Edp(Edp *SL, double **Z, double *res);


void dBug(const char *_msg_);
#endif //__ODESOLVER__
