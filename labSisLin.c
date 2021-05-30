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
    SistLinear_t *newSL = alocaSistLinear(SL->n);

    x = (real_t *)calloc(SL->n, sizeof(real_t));
    res = (real_t *)calloc(SL->n, sizeof(real_t));
    printf("Entrada\n");
    prnSistLinear(SL);
    
    //***salvando sis
    for(int i=0; i<SL->n; i++)
        for(int j=0; j<SL->n; j++)
            newSL->A[i][j] = SL->A[i][j];
    for(int i=0; i<SL->n; i++)
        newSL->b[i] = SL->b[i];//???
    newSL->n = SL->n;
    newSL->erro = SL->erro;
    //*****
    
    eliminacaoGauss(SL, x, tTotal);
    printf("\nDepois do gauss\n\n");
    prnSistLinear(SL);
    printf("\nVetor solução");
    prnVetor(x, SL->n);
    normaL2Residuo(newSL, x, res);

}

