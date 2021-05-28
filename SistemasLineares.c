#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "SistemasLineares.h"

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param res Valor do resíduo

  \return Norma L2 do resíduo.
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res)
{

}


/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTriangulariza tempo gasto na triangularização
  \param tRetroSubst tempo gasto na retrosubstituição

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal)
{


}

/*!
  \brief Método de Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo gasto pelo método
  \param tIteração tempo gasto em cada iteração

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
*/
int gaussJacobi (SistLinear_t *SL, real_t *x, double *tTotal)
{


}

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo gasto pelo método
  \param tIteração tempo gasto em cada iteração

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, double *tTotal)
{


}


/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial para início do refinamento
  \param erro menor erro aproximado para encerrar as iterações

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int refinamento (SistLinear_t *SL, real_t *x, double *tTotal)
{


}

/*!
  \brief Alocaçao de memória 

  \param n tamanho do SL

  \return ponteiro para SL. NULL se houve erro de alocação
  */
SistLinear_t* alocaSistLinear (unsigned int n)
{
  SistLinear_t *SL = (SistLinear_t*)malloc(sizeof(SistLinear_t));
  SL->A = (real_t **)malloc(n*sizeof(real_t*));
  for(int i=0; i<n; i++) SL->A[i] = (real_t *) malloc(n*sizeof(real_t));
  SL->b = (real_t *)malloc(n*sizeof(real_t));

  if(SL == NULL || SL->A == NULL || SL->b == NULL)
    return NULL;
  else
    return SL;

}

/*!
  \brief Liberaçao de memória 

  \param sistema linear SL
  */
void liberaSistLinear (SistLinear_t *SL)
{
  for(int i=0; i<SL->n; i++) free(SL->A[i]);
  free(SL->A);
  free(SL->b);
  free(SL);

}

/*!
  \brief Leitura de SL a partir de Entrada padrão (stdin).

  \return sistema linear SL. NULL se houve erro (leitura ou alocação)
  */
SistLinear_t *lerSistLinear ()
{
  //*******falta checar erro de leitura
  int num; real_t erro;
  scanf("%d %f", &num, &erro);
  while((getchar()) != '\n');//limpando buffer

  SistLinear_t *SL = alocaSistLinear(num);
  if(SL == NULL) return NULL;
  SL->n = num;
  SL->erro = erro;

  for(int i=0; i<SL->n; i++){
    for(int j=0; j<SL->n; j++)
      scanf("%f", &SL->A[i][j]);
    while((getchar()) != '\n');//limpando buffer
  }

  for(int i=0; i<SL->n; i++)
    scanf("%f", &SL->b[i]);

  return SL;

}


// Exibe SL na saída padrão
void prnSistLinear (SistLinear_t *SL)
{
  for(int i=0; i<SL->n; i++){
    for(int j=0; j<SL->n; j++){
      printf("%f ", SL->A[i][j]);
    }
    printf("\n");

  }
  for(int i=0; i<SL->n; i++)
    printf("%f ", SL->b[i]);
  printf("\n");

}

// Exibe um vetor na saída padrão
void prnVetor (real_t *v, unsigned int n)
{
  printf("\n");
  for(int i=0; i<n; i++)
    printf("%f ", v[i]);
  printf("\n");
}

