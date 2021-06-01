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
    long int pos;
    char c;
    

    while(1){
        SistLinear_t *SL = lerSistLinear();
        SistLinear_t *newSL = alocaSistLinear(SL->n);
        x = (real_t *)calloc(SL->n, sizeof(real_t));
        res = (real_t *)calloc(SL->n, sizeof(real_t));
        
        //***salvando sistema original
        for(int i=0; i<SL->n; i++)
            for(int j=0; j<SL->n; j++)
                newSL->A[i][j] = SL->A[i][j];
        for(int i=0; i<SL->n; i++)
            newSL->b[i] = SL->b[i];
        newSL->n = SL->n;
        newSL->erro = SL->erro;
        //*****

        
        eliminacaoGauss(SL, x, tTotal);
        printf("***** Sistema x --> n = %d, erro: %f\n", SL->n, SL->erro);
        printf("===> Eliminação Gauss: tempo\n --> X: ");
        prnVetor(x, SL->n);
        printf(" --> Norma L2 do residuo: %1.8e\n", normaL2Residuo(newSL, x, res));

        gaussJacobi(newSL, x, tTotal);//newSL pois SL foi modificado
        printf("\n===> Jacobi: tempo\n --> X: ");
        prnVetor(x, SL->n);
        printf(" --> Norma L2 do residuo: %1.8e\n", normaL2Residuo(newSL, x, res));
        
        gaussSeidel(newSL, x, tTotal);//newSL pois SL foi modificado
        printf("\n===> Gauss-Seidel: tempo\n --> X: ");
        prnVetor(x, SL->n);
        printf(" --> Norma L2 do residuo: %1.8e\n\n", normaL2Residuo(newSL, x, res));
        
        pos = ftell(stdin);
        for(c=fgetc(stdin); c == '\n' || c == ' '; c=fgetc(stdin))//searching for new linearSystems
            pos = ftell(stdin);//saving position of pointer
        if(c == EOF) break;
        else fseek (stdin, pos, SEEK_SET);//Number found!!
    
    }
    
    

}

