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
    eliminacaoGauss(SL, x, tTotal);
    prnVetor(x, 4);
    refinamento(SL, x, tTotal);
    prnVetor(x, 4);

}

