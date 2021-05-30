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
/* para cada linha a partir da primeira */   
    for (int i=0; i < SL->n; ++i) {      
	    int iPivo = encontraMax(SL, i);      
        if ( i != iPivo )         
            trocaLinha(SL, i, iPivo);
        for(int k=i+1; k < SL->n; ++k) {         
            double m = SL->A[k][i] / SL->A[i][i];         
            SL->A[k][i] = 0.0;         
            for(int j=i+1; j < SL->n; ++j)            
                SL->A[k][j] -= SL->A[i][j] * m;         
            SL->b[k] -= SL->b[i] * m;      
        }   
	}
    retrossubs(SL, x, SL->n);
}

real_t maxError(real_t * nextX, real_t *x, int n){
  real_t max = fabs(nextX[0] - x[0]);
  for(int i=1; i<n; i++){
    if(max < fabs(nextX[i] - x[i]))
      max = fabs(nextX[i] - x[i]);
  }
  return max;
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
  real_t *nextX;
  int num=0;
  nextX = (real_t *)calloc(SL->n,sizeof(real_t));  
  for(int i=0; i<SL->n; i++)
    x[i] = 0;

  while(1){  
    
    for(int i=0; i<SL->n; i++){
      for(int j=0; j<SL->n; j++){
        if(i == j) continue;
        nextX[i] += -SL->A[i][j]*x[j];
      }
      nextX[i] += SL->b[i];
      nextX[i] /= SL->A[i][i]; 
    }
    num++;
    prnVetor(nextX, 4);
    if(maxError(nextX, x, SL->n) <= (SL->erro)){//critério de parada
      for(int i=0; i<SL->n; i++)
        x[i] = nextX[i];
      return num;//modify to number of interations
    }
    
    for(int i=0; i<SL->n; i++){
      x[i] = nextX[i];
      nextX[i] = 0;
    }

  }

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
  real_t *nextX;
  int num=0;
  nextX = (real_t *)calloc(SL->n,sizeof(real_t));  
  for(int i=0; i<SL->n; i++)
    x[i] = 0;

  while(1){  
    
    for(int i=0; i<SL->n; i++){
      for(int j=0; j<SL->n; j++){
        if(i == j) continue;
        else if(j < i)
          nextX[i] += -SL->A[i][j]*nextX[j];//Using variables already written
        else 
          nextX[i] += -SL->A[i][j]*x[j];
      }
      nextX[i] += SL->b[i];
      nextX[i] /= SL->A[i][i]; 
    }
    num++;
    prnVetor(nextX, 4);
    if(maxError(nextX, x, SL->n) <= (SL->erro)){//critério de parada
      for(int i=0; i<SL->n; i++)
        x[i] = nextX[i];
      return num;//modify to number of interations
    }
    
    for(int i=0; i<SL->n; i++){
      x[i] = nextX[i];
      nextX[i] = 0;
    }

  }
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


void multiMatrix(SistLinear_t *SL, real_t *x, real_t *sol){
  //sol = [A].[X]
  for(int i=0; i<SL->n; i++)
    sol[i] = 0;

  for (int i = 0; i < SL->n; ++i) {
    for (int k = 0; k < SL->n; ++k)
      sol[i] += SL->A[i][k] * x[k];
  }

}   
  
int refinamento (SistLinear_t *SL, real_t *x, double *tTotal)
{
  //calculando resíduo
  real_t *aux;//A*X
  real_t *w;
  real_t *r;
  aux = (real_t *)malloc(SL->n*sizeof(real_t));
  w = (real_t *)malloc(SL->n*sizeof(real_t));//residuo
  r = (real_t *)malloc(SL->n*sizeof(real_t));//residuo
  
  multiMatrix(SL, x, aux);
  for(int i=0; i<SL->n; i++)
    r[i] = SL->b[i] - aux[i];
  
  SistLinear_t *newSL = alocaSistLinear(SL->n);
  for(int i=0; i<SL->n; i++)
    for(int j=0; j<SL->n; j++)
      newSL->A[i][j] = SL->A[i][j];
  newSL->b = r;//???
  eliminacaoGauss(newSL, w, tTotal); //Solving A*w = r 
  for(int i=0; i<SL->n; i++)
    x[i] += w[i];


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

