#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"



int main ()
{
    real_t *x;
    double *tTotal;
    SistLinear_t *SL = lerSistLinear();
    x = (real_t *)calloc(SL->n, sizeof(real_t));
    prnSistLinear(SL);
    printf("%d", gaussSeidel(SL, x, 0));
    

}

