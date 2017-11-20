#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int num_threads;


int getDimensaoMatriz() { 

    FILE *fp = fopen("matriz.txt","r");    
    
    int linesize = 1024;        
    char line[linesize];            

    char *token = NULL;
    const char *espaco = " "; //delimitador

    int nele = 0; //auxiliar pra contar o numero de elementos
    
    fgets(line, sizeof(line), fp);
    token = strtok(line, espaco);            
    
    while(token != NULL) {
        //printf("%s\n",token);
        token=strtok(NULL,espaco);                
        nele++;
    }                                    

    fclose(fp);
    return nele;
}   

//retorna a matriz(vetor..) double lida do arquivo texto
double *ler_matriz_arquivo(int dimensao) { 

    float num;
    
    int linesize = 1024;        
    char line[linesize];            
    
    int i=0;    

    char *token = NULL;
    
    const char *espaco = " ";
    
    double *A = malloc(dimensao*dimensao*sizeof(double));

    
    FILE *fp = fopen("matriz.txt","r");    
    
    if( fp == NULL ) {                       
        printf("problemas ao abrir o arquivo");        
        return 0;
    }

    while (!feof (fp)) {
        if (fgets(line, sizeof (line), fp)) {            
            //printf("%s",line);            
            token = strtok(line, espaco);            
            while(token != NULL) {                
                num = atof(token);                                                
                //printf("%.3f\n",num);
                A[i] = num;                                                                
                token=strtok(NULL,espaco);                
                i++;
            }                                    
        }
    }    
    fclose(fp);        
    return A;
}

double *ler_vetor_arquivo(int dimensao) { 

    float num;
    
    int linesize = 1024;        
    char line[linesize];            
    
    int i=0;    

    char *token = NULL;
    
    const char *espaco = " ";
    
    double *A = malloc(dimensao*sizeof(double));

    
    FILE *fp = fopen("vetor.txt","r");    
    
    if( fp == NULL ) {                       
        printf("problemas ao abrir o arquivo");        
        return 0;
    }

    while (!feof (fp)) {
        if (fgets(line, sizeof (line), fp)) {            
            //printf("%s",line);            
            token = strtok(line, espaco);            
            while(token != NULL) {                
                num = atof(token);                                                
                //printf("%.3f\n",num);
                A[i] = num;                                                                
                token=strtok(NULL,espaco);                
                i++;
            }                                    
        }
    }
    
    fclose(fp);
        
    return A;
}

int escrever_y_arquivo(double *b,int dim) { 
    
    FILE *fp = fopen("resposta.txt","w+");    
           
    for(int j=0;j<dim;j++) {
        //printf("%2.f\n",array[j]);
         fprintf(fp,"%f\n",b[j]);         
    }        
    return 0;
}


void print_matriz(double *M,int nlinhas,int ncolunas)  { 
    int i,j;
    
    for (i=0; i < nlinhas; i++) {                
        for(j=0; j < ncolunas; j++) {
            printf("\t%.1f",M[i*ncolunas+j]);
        }
        printf("\n");
    }    
}


void substituicao_reversa_com_print(double *A, double *y, int n) {
  
    int i, j;  
    double x[n];

    for(i =n-1; i >= 0; i--) {
        x[i] = y[i];
        for(j = n-1; j > i; j--) {
            x[i] = x[i] - x[j] * A[i*n+j];        
        }
        y[i] = x[i];
    }
  
  //printa o resultado
    for(i = 0; i < n; i=i+4) {
        for(j = i; j<i+4 && j<n; j++ ) { 
            printf("x[%d]=%.3f ",j,x[j]);        
        }
        printf("\n");        
    }

}

void substituicao_reversa(double *A, double *y, int n) {
  
    int i, j;  
    double x[n];

    for(i =n-1; i >= 0; i--) {
        x[i] = y[i];
        for(j = n-1; j > i; j--) {
            x[i] = x[i] - x[j] * A[i*n+j];        
        }
        y[i] = x[i];
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

void gauss_serial(double *A,double *b,double *y,int n) { 

    int i,j,lp;

    for(lp = 0; lp < n; lp++) {
        for (j=lp+1; j < n; j++) { 
            A[lp*n+j] = A[lp*n+j] / A[lp*n+lp]; 
            //printf("oi\n");
        }
        y[lp] = b[lp] / A[lp*n+lp];
        A[lp*n+lp] = 1.0;    
    
        //calculada a linha pivo; fazer para as proximas lp-n linhas        
        #pragma omp parallel default(none) shared(A,b,y,lp,n) private(i,j) num_threads(num_threads)
        {        
            int rank = omp_get_thread_num();
            #pragma omp for schedule(auto)
            for (i=lp+1; i < n; i++) { 
                for (j=lp+1; j < n; j++) { 
                    A[i*n+j] = A[i*n+j] - A[i*n+lp] * A[lp*n+j];
                    //printf("rank %d - A[%d] = A[%d] - A[%d] * A[%d]\n",rank,i*n+j,i*n+j,i*n+lp,lp*n+j);
                }
                b[i] = b[i] - A[i*n+lp] * y[lp];
                //printf("rank %d - b[%d] = b[%d] - A[%d] * y[%d]\n",rank,i,i,i*n+lp,lp);
                A[i*n+lp] = 0.0;
            }        
        }
    
    //printf("\n");
    }
} 




int main (int argc, char *argv[]) {     

    double *A,*b,*y;
    int dimensao;
    double tfinal,tinicial;
        
    if (argc > 1) {
        num_threads = atoi(argv[1]);
    }
    else {         
        printf("uso: ./executavel numero_threads\n");
        return 1;
    }

    dimensao = getDimensaoMatriz();
        
    A = ler_matriz_arquivo(dimensao);
    b = ler_vetor_arquivo(dimensao);
    y = (double*)malloc(dimensao*sizeof(double));
                
    // testes 
    // dimensao=2000;
    
    // A = (double*)malloc(dimensao*dimensao*sizeof(double));
    // b = (double*)malloc(dimensao*sizeof(double));
    // y = (double*)malloc(dimensao*sizeof(double));

    // gerar_matriz(A,b,dimensao);
    
    //print_matriz(A,dimensao,dimensao);

    tinicial = omp_get_wtime();    
    gauss_serial(A,b,y,dimensao);
    tfinal = omp_get_wtime();
        
    //printf("\n");
    //print_matriz(A,dimensao,dimensao);
    //printf("\n");
    substituicao_reversa(A,y,dimensao);
    escrever_y_arquivo(y,dimensao);
    printf("num_threads = %d dimensao = %d tempo = %f\n",num_threads,dimensao,tfinal-tinicial);
    
    return 0;

}