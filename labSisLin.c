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
    prnSistLinear(SL);
    x = multiMatrix(SL, SL->b);
    prnVetor(x, 4);

}

