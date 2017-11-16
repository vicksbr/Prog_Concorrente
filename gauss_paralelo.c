#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include <string.h>
#include <unistd.h>


//vars globais

int rank,size;
// map[i] == quantas linhas o processador i possui 

// depois uso isso para encontrar o indice das linhas referentes a cada processador 

void print_matriz( double *Matriz, int y, int x ) {  
  
    int i,j;
    for(i=0; i<y; i++) {
        for( j=0; j<x; j++) {
            printf("\t%.1lf",Matriz[i*x+j]);            
        }
    printf("\n");
    }
}

void meu_print_matriz(double *M,int nlinhas,int ncolunas)  { 
    int i,j;
    
    for (i=0; i < nlinhas; i++) {                
        for(j=0; j < ncolunas; j++) {
            printf("\t%.1f",M[i*nlinhas+j]);
        }
        printf("\n");
    }    
}

void substituicao_reversa_com_print( double *A, double *y, int n) {
  
  int i, j;  
  double x[n];

  for(i =n-1; i >= 0; i--) {
    x[i] = y[i];
    for(j = n-1; j > i; j--) {
      x[i] = x[i] - x[j] * A[i*n+j];        
    }
  }
  
  //printa o resultado
  for(i = 0; i < n; i=i+4) {
    for(j = i; j<i+4 && j<n; j++ ) { 
        printf("x[%d]=%.2f ",j,x[j]);        
    }
    printf("\n");        
  }
}


// map[i] => quantos linhas o processo i possui
void mapearProcessos_linhas(int *map, int nproc,int nlinhas) { 

    int i;
    for (i = 0; i < nlinhas; i++) 
        map[i] = i%nproc;
}


void mapearProcesssos_Linhas(int nlinhas,int *map) { 
    
    int i,j; 
    
    for (j = 0; j < size; j++)
        map[j] = 0;    
    
    for (i = size-1; i>=0; i--)
        for (j = 0; j < nlinhas; j++)
            if (j%size == i) { 
                map[i] += 1;
            }
}

void gerar_matriz_inteiros(double *A,double *b,int nlinhas,int ncolunas) { 
  
   int i,j,indice;
   
    for(i=0; i<nlinhas; i++ ) {	
        b[i] = 0;			
        for(j=0; j<ncolunas; j++ ) {	            
            indice = (i*nlinhas+j);
            A[i*nlinhas+j] = indice+1;  //A[i][j] = numero
            b[i] = b[i] + j*indice; 
        }   
    }    
}

void gerar_matriz(double *A,double *b,int n) { 
  
   int i,j,num_ram;
   
    for(i=0; i<n; i++ ) {	
        b[i] = 0.0;			
        for(j=0; j<n; j++ ) {	
            num_ram = rand();		
            A[i*n+j] = num_ram;  //A[i][j] = numero
            b[i] = b[i] + j*num_ram; 
        }   
    }    
}


void distribuirMatrizes(double *A,double *matrizLocal,int nlinhas,int ncolunas,int *map) { 
 
    MPI_Status status;

    int n,x,i,j;
    
    int procAtual;

    if (rank == 0) {    
        for(procAtual = size-1; procAtual >= 0; procAtual--) { 
            for (i = 0; i < map[procAtual]; i++) {                                     
                //printf("linha %d para processo %d\n",(procAtual + i*size),procAtual);
                for (j = 0; j < ncolunas; j++) { 
                    //printf("rank %d para %d: matrizLocal[%d] = A[%d] = %f\n",rank,procAtual,i*ncolunas +j,(procAtual + i*size)*ncolunas+j,A[i*ncolunas +j,(procAtual + i*size)*ncolunas+j]);
                    matrizLocal[i*ncolunas+j] = A[(procAtual+i*size)*ncolunas+j];                    
                }
            }
            if(procAtual != 0) {                                                
                MPI_Send(matrizLocal,map[procAtual]*ncolunas, MPI_DOUBLE, procAtual, 10, MPI_COMM_WORLD );                            
            }
        }
    }
    else {        
        MPI_Recv(matrizLocal,map[rank]*ncolunas, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, &status );
    }
}


void juntar_matriz(double *A, double *matrizLocal, int nlinhas, int ncolunas,int *map) { 
    
    MPI_Status status;

    int i, j, proc;    
    //Processo 0 recebe de todos os outros as matrizes locais    

    double *novaMatrizLocal;
    
    if (rank == 0) { 
        for(proc = 0; proc < size; proc++) {
                        
            //precisamos disso pois se a matrizlocal dos outros processos tiver tamanho diferente 
            //da matrizlocal do processo 0, segmentation fault, buffer overflow  
            novaMatrizLocal = (double*)malloc(map[proc]*ncolunas*sizeof(double));            
            
            if (proc != 0) {                                                               
                MPI_Recv(novaMatrizLocal,map[proc]*ncolunas, MPI_DOUBLE, proc, 10, MPI_COMM_WORLD, &status );
                // printf("rank %d recebendo\n",rank);
                // print_matriz(n,map[proc],ncolunas);                 
            }            
            for(i = 0; i < map[proc]; i++) {                 
                for(j = 0; j < ncolunas; j++) {                                                             
                    //printf("A[%d] = matrizLocal[%d] do processo: %d\n",(proc+i*size)*ncolunas+j,i*ncolunas+j,proc);
                    if (proc == 0)
                        A[(proc+i*size)*ncolunas+j] = matrizLocal[i*ncolunas+j];
                    else 
                        A[(proc+i*size)*ncolunas+j] = novaMatrizLocal[i*ncolunas+j];
                }             
            
            }                          
            
        }    
    }
    else { 
        MPI_Send(matrizLocal,map[rank]*ncolunas, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
        // printf("enviando do %d com tamanho %d\n",rank,map[rank]*ncolunas);
        // print_matriz(matrizLocal,map[rank],ncolunas);            
        // printf("\n");
    }
}

void gaus_paralelo(double *A,double *b,double *y,int nlinhas,int ncolunas,int *map) { 

    MPI_Status status;
    MPI_Request request;
    MPI_Comm com = MPI_COMM_WORLD;
        
    double *matrizLocal;    
    double *bLocal;
    double *yLocal;
    double *linhapivo;    

    int nlocais = map[rank]; // numero de linhas para cada processo    
    
    matrizLocal = (double*)malloc(nlocais*ncolunas*sizeof(double));
    bLocal = (double*)malloc(nlocais*sizeof(double));
    yLocal = (double*)malloc(nlocais*sizeof(double));    
    linhapivo = (double*)malloc(ncolunas*sizeof(double)+1); //usado para enviar as linhas pivo entre os processos

    int rank_predecessor = (size+(rank-1)) % size; 
    int rank_sucessor = (rank+1) % size;

    distribuirMatrizes(A,matrizLocal,nlinhas,ncolunas,map);    
    distribuirMatrizes(b,bLocal,nlinhas,1,map);

    
    // print_matriz(matrizLocal,nlocais,ncolunas);
    // printf("\n");
    // printf("\n");
    // printf("\n");
    int i,j;
    int k; //k controla a linha pivo
    int klocal;

    for (k = 0;k < nlinhas; k++) { 
        if (rank == k%size) { 
            //printf("calcular a linha pivo rank %d k = %d\n",rank,k);
            klocal = (k/size)*nlinhas;
            
            //printf("rank %d elemento pivo %.2f\n",rank,A[klocal]+k);      
            //printf("rank %d klocal %d\n",rank,klocal);
            for (j = k+1; j < ncolunas; j++) {
                matrizLocal[klocal+j] = matrizLocal[klocal+j] / matrizLocal[klocal+k];                
                //printf("matrizLocal[%d] = matrizlocal[%d])/matrizlocal[%d]\n",klocal+j,klocal+j,klocal+k);
            }
            yLocal[k/size] = bLocal[k/size] / matrizLocal[klocal+k];            
            //printf("rank %d yLocal[%d] = bLocal[%d] / matrizLocal[%d]\n",rank,k/size,k/size,klocal+k);
            //printf("rank %d %2.f = %2.f / %2.f\n",rank,yLocal[k/size],bLocal[k/size],matrizLocal[klocal+k]);
            matrizLocal[klocal+k] = 1.0;

            for (j = 0; j < ncolunas; j++) { 
                linhapivo[j] = matrizLocal[klocal + j];            
                //printf("rank %d linhapivo[%d] = matrizLocal[%d]\n",rank,j,klocal+j);
            }
            linhapivo[ncolunas] = yLocal[k/size];                             

            MPI_Send(linhapivo, ncolunas+1, MPI_DOUBLE, rank_sucessor, 20, com );	//manda pivo pro sucessor 

        }
        else {
            //recebe a linha pivo do processo anterior 
            MPI_Recv(linhapivo, ncolunas+1, MPI_DOUBLE, rank_predecessor, 20, com, &status );
            
            
            //se manda a linha pivo recebida do processo anterior pro processo sucessor
            //menos no caso dele ter recebido a linha pivo do próprio processo sucessor            
            //k%size representa o numero do processo de quem veio a linha
            
            if (k%size != rank_sucessor) {                
                MPI_Send(linhapivo, ncolunas+1, MPI_DOUBLE, rank_sucessor, 20, com );    
            }
            
        }
    
        //nesse momento cada processo tem a linha pivo necessária pra realizar seus calculos 
        //abaixo da linha diagonal local calculada
        
        
        //se o meu rank for menor do rank de quem me mandou a linha        
        
         int comeco;
         int in,ink; 
        // comeco = (int)k/size;                            
                
        //printf("k = %d rank %d começar a eliminação a partir da linha %d\n",k,rank,comeco);
        
        if (rank <= k%size) {
            comeco = (int)k/size+1;            
            printf("rank %d k = %d comeco(0) = %d/%d = %d\n",rank,k,k,size,k/size);
        }
        else {
            comeco = (int)k/size;           
            printf("rank %d k = %d comeco(1) = %d/%d = %d\n",rank,k,k,size,k/size);

        }

        
        for (i = comeco; i < map[rank]; i++) {             
            
            in = i*ncolunas; 
            ink = in+k;
            
            for (j = k+1; j < ncolunas; j++) {                 
                matrizLocal[in+j] = matrizLocal[in+j] - matrizLocal[ink]*linhapivo[j];                                    
            }
            bLocal[i] = bLocal[i] - matrizLocal[ink]*linhapivo[ncolunas];
            matrizLocal[ink] = 0.0;
        }           
    }
    
    
    
    juntar_matriz(A,matrizLocal,nlocais,ncolunas,map);        
    juntar_matriz(y,yLocal,nlinhas,1,map);     

}


int main () { 

    int nlinhas = 5;
    int ncolunas = 5;                
    
    int *map;
    
    MPI_Init (0,0);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);    
    
    //if (rank == 0) {
    
        double *A;
        double *b;
        double *y; //variaveis onde serão alocadas as matrizes
        
        A = (double*)malloc(nlinhas*ncolunas*sizeof(double));
        b = (double*)malloc(nlinhas*sizeof(double));
        y = (double*)malloc(nlinhas*sizeof(double));
    
        gerar_matriz(A,b,nlinhas);
    
        map = (int*)malloc(size*sizeof(int));                
        mapearProcesssos_Linhas(nlinhas,map);        
        
        
    //Matriz teste
    // double A[] = {33,47,50,30,18,38,13,59,17,1,60,35,40,87,4,32,90,5,37,91,88,67,31,85,68};
    // double b[] = {11,18,5,12,40};


    if (rank == 0) { 
        print_matriz(A,nlinhas,ncolunas);
    }
    
    gaus_paralelo(A,b,y,nlinhas,ncolunas,map);    
    
    MPI_Barrier( MPI_COMM_WORLD );

    if (rank == 0) { 
        printf("\n");
        print_matriz(A,nlinhas,ncolunas);
        printf("\n");
        print_matriz(b,1,ncolunas);
        printf("\n");
        print_matriz(y,1,ncolunas);
        printf("\n");
        substituicao_reversa_com_print(A,y,nlinhas);
        printf("\nprocesso 0 dizendo adeus\n");
    }
    
    MPI_Finalize();
    
    return 0;
}