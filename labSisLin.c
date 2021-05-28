#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int encontraMax(SistLinear_t *SL, int i){
    double maior = SL->A[i][i];
    int maiorLinha = i;
    for(int linha=i+1; linha<SL->n; linha++)
        if(fabs(SL->A[linha][i]) > fabs(maior)){
            maior = SL->A[linha][i];
            maiorLinha = linha;
        }
    return maiorLinha;
}

void trocaLinha(SistLinear_t *SL, int i, int iPivo){
    double aux;
    for(int coluna=0; coluna<SL->n; coluna++){
        aux = SL->A[iPivo][coluna];
        SL->A[iPivo][coluna] = SL->A[i][coluna];
        SL->A[i][coluna] = aux;
    }

    aux = SL->b[iPivo];
    SL->b[iPivo] = SL->b[i];
    SL->b[i] = aux;

}

void retrossubs(SistLinear_t *SL, real_t *x, int n){
    for(int i=n-1; i>=0; --i){
        x[i] = SL->b[i];
        for(int j = i+1; j<n; ++j)
            x[i] -= SL->A[i][j]*x[j];
        x[i] /= SL->A[i][i];
    }
}



// int eliminacaoGauss(SistLinear_t *SL, real_t *x, double *tTotal) {
// /* para cada linha a partir da primeira */   
//     for (int i=0; i < SL->n; ++i) {      
// 	    int iPivo = encontraMax(SL, i);      
//         if ( i != iPivo )         
//             trocaLinha(SL, i, iPivo);
//         for(int k=i+1; k < SL->n; ++k) {         
//             double m = SL->A[k][i] / SL->A[i][i];         
//             SL->A[k][i] = 0.0;         
//             for(int j=i+1; j < SL->n; ++j)            
//                 SL->A[k][j] -= SL->A[i][j] * m;         
//             SL->b[k] -= SL->b[i] * m;      
//         }   
// 	}
//     retrossubs(SL, x, SL->n);
// }


int main ()
{


    SistLinear_t *SL = lerSistLinear();
    prnSistLinear(SL);

}

