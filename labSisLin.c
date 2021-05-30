#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"



int main ()
{
    real_t *x;
    real_t *res;
    double *tTotal;
    SistLinear_t *SL = lerSistLinear();
    x = (real_t *)calloc(SL->n, sizeof(real_t));
    res = (real_t *)calloc(SL->n, sizeof(real_t));
    prnSistLinear(SL);
    gaussSeidel(SL, x, tTotal);
    prnVetor(x, SL->n);
    refinamento(SL, x, tTotal);
    prnVetor(x, SL->n);
}

