//EP2 - Metodos Numericos - MAP3121

// Gabriel da Cunha Rodrigues - No USP: 8992930 - Turma 2

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "fftpack4.h"

/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Declaracao de funcoes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

/* >>>>>>>>>>> Funcoes basicas <<<<<<<<<<< */
double* criar_vetor(int N);
int* criar_vetor_int(int N);
double complex* criar_vetor_complexo(int N);
void imprimir_vetor(double* vetor, int N);
void imprimir_complexo(double complex *c, int N);

/* >>>>>>>>>>> Funcoes de Transformada de Fourier <<<<<<<<<<< */
int tamanho_arquivo(char *nome_arquivo, int *canais);
int ler_arquivo_dat(char *nome_arquivo, int n, double complex *x, double complex *f, double complex *f2, double *f_linha, double *f_linha2);
void escrever_arquivo_dat(char *nome_arquivo, int sample_rate, int canais, int n, double complex *x, double complex *f, double complex *f2);
void fourier(double complex *c, double complex *f, double complex *x, int n);
void anti_fourier(double complex *c, double complex *f, double complex *x, int n);
void fftrec(double complex *c, double complex *f, int n, bool dir);
void comprimir_sinal(double complex *c, double taxa_minima, int N);
void passa_altas(double complex *c, int N, int freq_corte );
void passa_baixas(double complex *c, int N, int freq_corte);
void passa_bandas(double complex *c, int N, int freq1, int freq2);

/* ------------------------------------------------------------------------------------- */

int main() {
	// Variaveis globais
	int n; // numero total de amostras. Usado para a alocacao de vetores
	double complex *c; // Vetor de coeficientes
	double complex *x; // Vetor de amostras 
	double complex *f; // Vetor de valores da função em cada amostra
	int tipo_problema, tipo_transformada;
	bool direta; // Usado na fftrec

	// Variaveis especificas dos arquivos de audio
	int escolha, tipo_filtro, tipo_compressao;
	char nome_arquivo[128];
	double complex *c2; // Vetor de coeficientes para o segundo canal
    double complex *f2; // Vetor de valores da função para o segundo canal
    int canais, sample_rate; // Parametros fornecidos pelo arquivo
    int K, K1, K2; // Parametros de corte utilizados nos filtros
    int S; // Parametro de compressao

	// Variaveis especificas do fftpack4
	double *a, *b, *a2, *b2;
	double a0, a02; 
	int *ifac, *ifac2;
	double *wsave, *wsave2;
	double *f_linha, *f_linha2; 

    printf("EP2 - Numerico\n\n");
    printf("1 - Teste a\n");
    printf("2 - Teste b\n");
    printf("3 - Teste c\n");
    printf("4 - Analise de arquivo\n\n");
    printf("Digite o numero do problema a ser resolvido: ");
    scanf("%d", &tipo_problema);
    printf("\n");

    switch(tipo_problema) {

        case 1: // Teste a
        	n = 2;

			x = criar_vetor_complexo(2*n);
			x[0] = 0;
   			x[1] = M_PI / 2;
    		x[2] = M_PI;
    		x[3] = 3 * M_PI / 2;

    		f = criar_vetor_complexo(2*n);
    		f[0] = 5;
		    f[1] = -1;
		    f[2] = 3;
		    f[3] = 1;

		    f_linha = criar_vetor(2*n);  // No FFTPACK4, trabalha-se com valores em double e nao complexos:
		    f_linha[0] = 5;
		    f_linha[1] = -1;
		    f_linha[2] = 3;
		    f_linha[3] = 1;
            break;

        case 2: // Teste b
        	n = 4;

			x = criar_vetor_complexo(2*n);
			x[0] = 0;
   			x[1] = M_PI / 4;
    		x[2] = M_PI / 2;
    		x[3] = 3 * M_PI / 4;
    		x[4] = M_PI;
    		x[5] = 5 * M_PI / 4;
    		x[6] = 3 * M_PI / 2;
    		x[7] = 7 * M_PI / 4;

			f = criar_vetor_complexo(2*n);
			f[0] = 6;
		    f[1] = 2;
		    f[2] = 5;
		    f[3] = 2;
		    f[4] = 11;
		    f[5] = 2;
		    f[6] = 8;
		    f[7] = 8;

		    f_linha = criar_vetor(2*n); // No FFTPACK4, trabalha-se com valores em double e nao complexos
	    	f_linha[0] = 6;
			f_linha[1] = 2;
			f_linha[2] = 5;
			f_linha[3] = 2;
			f_linha[4] = 11;
			f_linha[5] = 2;
			f_linha[6] = 8;
			f_linha[7] = 8;
            break;

        case 3: // Teste c
        	n = 512;
			x = criar_vetor_complexo(2*n);
			f = criar_vetor_complexo(2*n);
			f_linha = criar_vetor(2*n); // No FFTPACK4, trabalha-se com valores em double e nao complexos
			
			for(int i = 0; i < 2*n; i++) {
		        x[i] = (i * M_PI) / 512; 

		        f[i] = 10*sin(x[i]) + 7*cos(30*x[i]) + 11*sin(352*x[i]) - 8*cos(711*x[i]);
		        f_linha[i] = 10*sin(x[i]) + 7*cos(30*x[i]) + 11*sin(352*x[i]) - 8*cos(711*x[i]);
		    }
            break;

        case 4: // Arquivos de audio
		    printf("Digite o nome do arquivo a ser analizado (com a terminacao .dat): ");
		    scanf("%s", nome_arquivo);
			n = tamanho_arquivo(nome_arquivo, &canais); // Define o numero de amostras com base no arquivo
			printf("\n---------------Analise de Fourier----------------\n");

			// Alocacao de vetores
			x = criar_vetor_complexo(2*n);
			f = criar_vetor_complexo(2*n);
			f2 = criar_vetor_complexo(2*n);
			// No FFTPACK4, trabalha-se com valores em double e nao complexos
			f_linha = criar_vetor(2*n); 
			f_linha2 = criar_vetor(2*n);
			// Preenchimento dos vetores de acordo com o arquivo
			sample_rate = ler_arquivo_dat(nome_arquivo, n, x, f, f2, f_linha, f_linha2);
            break;

        default:
            printf("Opcao invalida!\n");
            break;
    }

	if(tipo_problema != 4) { // Testes a, b, c:
		printf("\nO sinal analisado tem %d amostras. \n", n);	
		c = criar_vetor_complexo(2*n);

		printf("\n------------- Forma Direta ---------------\n\n");
		printf("Amostra original do sinal:\n");
		imprimir_complexo(f, 2*n);

		printf("Transformada de Fourier - Vetor de coeficientes:\n");
		fourier(c, f, x, n);
		imprimir_complexo(c, 2*n);

	    printf("Antitransformada de Fourier - Sinal recuperado:\n");
	    anti_fourier(c, f, x, n);
		imprimir_complexo(f, 2*n);

		printf("\n------------- FFT Recursiva ---------------\n\n");
		printf("Amostra original do sinal:\n");
	    imprimir_complexo(f, 2*n);

		printf("Transformada de Fourier - Vetor de coeficientes:\n");
		fftrec(c, f, n, true);  // Transformada direta pela fftrec

		// Como indicado no algoritmo, faz-se necessaria a divisao dos coeficientes por 2N
		for(int i = 0; i < 2*n; i++){
	        c[i] = c[i]/(2*n);
	    }
	    imprimir_complexo(c, 2*n);

	    printf("Antitransformada de Fourier - Sinal recuperado:\n");
		fftrec(c, f, n, false);  // Anti-transformada pela fftrec	
	    imprimir_complexo(f, 2*n);

		printf("\n------------- FFTPACK4 ---------------\n\n");
	    printf("Amostra original do sinal:\n");
	    imprimir_complexo(f, 2*n);
		
		wsave = criar_vetor(3 * n + 15);
		ifac = criar_vetor_int(8);
		a = criar_vetor(n/2);
		b = criar_vetor(n/2);

		n = 2*n; // Dobra-se n para o uso nas funções do fftpack4
		ezffti(&n, wsave, ifac);  // inicializacao da fftpack4
		ezfftf(&n, f_linha, &a0, a, b, wsave, ifac);  // transformada direta de fourier
		n = n/2; // Retorna-se para o valor de n original

		// Conversao de valores do tipo a*cos() + b*sen() para coeficientes complexos do tipo ck
		c[0] = a0;
		for(int i = 1; i < n; i++) {
			c[i] = (a[i-1] - (I * b[i-1]))/2;	
		}
		c[n] = a[n-1] + I*b[n-1];
		for(int i = 1; i < n; i++) {
			c[2*n-i] = (a[i-1] + (I * b[i-1]))/2;
		}

		printf("Transformada de Fourier - Vetor de coeficientes:\n");
	    imprimir_complexo(c, 2*n);

	    printf("Antitransformada de Fourier - Sinal recuperado:\n");
		n = 2*n; // Dobra-se n para o uso nas funções do fftpack4
		ezfftb(&n, f_linha, &a0, a, b, wsave, ifac);  // Antitransformada de fourier
	    imprimir_vetor(f_linha, n);
	    n = n/2; // Retorna-se para o valor de n original
	    return 0;

        
    } else { // Para os arquivos de audio
    	printf("\nO sinal analisado tem %d amostras. \n", n);	
    	c = criar_vetor_complexo(2*n);
    	if(canais == 2) {
			c2 = criar_vetor_complexo(2*n);
		}
    	
		printf("1 - Forma direta\n");
		printf("2 - FFT Recursiva\n");
		printf("3 - FFTPACK4\n");
		printf("Selecione o tipo de transformada desejada: ");
		scanf("%d", &tipo_transformada);


		// Transformada de Fourier

		switch(tipo_transformada) {
			case 1:

				fourier(c, f, x, n);
				if(canais == 2) {
					fourier(c2, f2, x, n);
				}
				printf("\nTransformada de Fourier realizada.\n\n");
				break;

			case 2:
				fftrec(c, f, n, true);
				if(canais == 2) {
					fftrec(c2, f2, n, true);
				}
				// Como indicado no algoritmo, faz-se necessaria a divisao dos coeficientes por 2N
				for(int i = 0; i < 2*n; i++){
			        c[i] = c[i]/(2*n);
			    }
			    if(canais == 2) {
			    	for(int i = 0; i < 2*n; i++){
			        	c2[i] = c[i]/(2*n);
			   		}
			    }
				printf("\nTransformada de Fourier (Forma Recursiva) realizada.\n\n");
			    break;

			case 3:
				
				wsave = criar_vetor(3 * n + 15);
				ifac = criar_vetor_int(8);
				a = criar_vetor(n/2);
				b = criar_vetor(n/2);
				wsave2 = criar_vetor(3 * n + 15);
				ifac2 = criar_vetor_int(8);
				a2 = criar_vetor(n/2);
				b2 = criar_vetor(n/2);

				n = 2*n; // Dobra-se n para o uso nas funções do fftpack4
				ezffti(&n, wsave, ifac);  // inicializacao da fftpack4
				ezfftf(&n, f_linha, &a0, a, b, wsave, ifac);  // transformada
				n = n/2; // Retorna-se para o valor de n original
			
				// Conversao de valores do tipo a*cos() + b*sen() para coeficientes complexos do tipo ck
				c[0] = a0;
				for(int i = 1; i < n; i++) {
					c[i] = (a[i-1] - (I * b[i-1]))/2;	
				}
				c[n] = a[n-1] + I*b[n-1];
				for(int i = 1; i < n; i++) {
					c[2*n-i] = (a[i-1] + (I * b[i-1]))/2;
				}

				if(canais == 2) {
					n = 2*n; // Dobra-se n para o uso nas funções do fftpack4
					ezffti(&n, wsave2, ifac2);  // inicializacao da fftpack4
					ezfftf(&n, f_linha2, &a02, a2, b2, wsave2, ifac2);  // transformada direta de fourier
					n = n/2; // Retorna-se para o valor de n original

					// Conversao de valores do tipo a*cos() + b*sen() para coeficientes complexos do tipo ck
					c2[0] = a02;
					for(int i = 1; i < n; i++) {
						c2[i] = (a2[i-1] - (I * b2[i-1]))/2;	
					}
					c2[n] = a2[n-1] + I*b2[n-1];
					for(int i = 1; i < n; i++) {
						c2[2*n-i] = (a2[i-1] + (I * b2[i-1]))/2;
					}
				}
				break;
			default:
				printf("Opcao invalida!\n");
            	break;
        }


        // Aplicacao de Filtros

        printf("Deseja aplicar um filtro no sinal?\n");
        printf("1 - Sim\n");
	    printf("2 - Nao\n");
	    printf("Resposta: ");
	    scanf("%d", &escolha);
	    if (escolha == 1) {
	    	printf("\n1 - Passa-altas\n");
	    	printf("2 - Passa-baixas\n");
	    	printf("3 - Passa-bandas\n");
	    	printf("Digite o filtro desejado: ");
	    	scanf("%d", &tipo_filtro);
	    	printf("O sinal analisado tem %d amostras. \n", n);	    
	    	switch(tipo_filtro) {
	    		case 1:
	    			printf("Digite o parametro de corte K: ");
	    			scanf("%d", &K);
	    			passa_altas(c, n, K);
	    			if (canais == 2) {
	    				passa_altas(c2, n, K);
	    			}
	    			printf("\nFiltro passa-altas aplicado.\n\n");
	    			break;
	    		case 2:
	    			printf("Digite o parametro de corte K: ");
	    			scanf("%d", &K);
	    			passa_baixas(c, n, K);
	    			if (canais == 2) {
	    				passa_baixas(c2, n, K);
	    			}
	    			printf("\nFiltro passa-baixas aplicado.\n\n");
	    			break;
	    		case 3:
	    			printf("Digite o parametro de corte K1: ");
	    			scanf("%d", &K1);
	    			printf("Digite o parametro de corte K2: ");
	    			scanf("%d", &K2);
	    			passa_bandas(c, n, K1, K2);
	    			if (canais == 2) {
	    				passa_bandas(c2, n, K1, K2);
	    			}
	    			printf("\nFiltro passa-bandas aplicado.\n\n");
	    			break;

	    	}

	    } else if(escolha != 2) {
	    	printf("Escolha invalida!\n");
	    }


	    // Compressao

	    printf("\nDeseja comprimir o sinal?\n");
        printf("1 - Sim\n");
	    printf("2 - Nao\n");
	    printf("Resposta: ");
	    scanf("%d", &escolha);
	    if (escolha == 1) {
	    	printf("Digite o parametro de compressao S: ");
	    	scanf("%d", &S);
	    	comprimir_sinal(c, S, n);
	    	if(canais == 2){
	    		comprimir_sinal(c2, S, n);
	    	}
	    	printf("\nSinal comprimido.\n\n");
	    } else if(escolha != 2) {
	    	printf("Escolha invalida!\n");
	    }


        // Antitransformada de Fourier

        switch(tipo_transformada) {
			case 1:
			    anti_fourier(c, f, x, n);
			    if(canais == 2){
			    	anti_fourier(c2, f2, x, n);
			    }
			    printf("\nAntitransformada de Fourier aplicada.\n\n");
				break;

			case 2:
				fftrec(c, f, n, false);
				if(canais == 2){
			    	fftrec(c2, f2, n, false);
			    }	
				printf("\nAntitransformada de Fourier (FFT Recursiva) aplicada.\n\n");
			    break;

			case 3:
				n = 2*n; // Dobra-se n para o uso nas funções do fftpack4
				ezfftb(&n, f_linha, &a0, a, b, wsave, ifac);  
			    if(canais == 2){
			    	ezfftb(&n, f_linha2, &a02, a2, b2, wsave2, ifac2);  
			    }	
			    n = n/2; // Retorna-se para o valor de n original
			    printf("\nAntitransformada de Fourier (FFTPACK4) aplicada.\n\n");
			    break;
        }


        //Gravando arquivo

		printf("Deseja gravar o resultado em um arquivo?\n");
	    printf("1 - Sim\n");
	    printf("2 - Nao\n");
	    printf("Resposta: ");
	    scanf("%d", &escolha);

	    if(escolha == 1) {
	    	printf("\nDigite o nome do arquivo a ser escrito (com a terminacao .dat): ");
	    	scanf("%s", nome_arquivo);
	    	// Escreve o resultado da analise no arquivo
	    	escrever_arquivo_dat(nome_arquivo, sample_rate, canais, n, x, f, f2);
	    	printf("\nArquivo gravado com sucesso!\n");
	    } 
    }

    // Desalocacao de memoria
    free(ifac);
	free(wsave);
	free(a);
	free(b);
	free(ifac2);
	free(wsave2);
	free(a2);
	free(b2);
	free(c);
	free(c2);
	free(x);
	free(f);
	free(f2);
	free(f_linha);
	free(f_linha2);
	return 0;
}


/* >>>>>>>>>>>>>>>>>>>>>>>> Funcoes basicas <<<<<<<<<<<<<<<<<<<<<<<< */

double* criar_vetor(int N) {
    /* Cria um vetor de N doubles. */

    double *vetor;
    vetor = (double*) calloc(N, sizeof(double));
    return vetor;
}

int* criar_vetor_int(int N) {
    /* Cria um vetor de N inteiros. */

    int *vetor;
    vetor = (int*) calloc(N, sizeof(int));
    return vetor;
}

double complex* criar_vetor_complexo(int N) {
    /* Cria um vetor de N doubles complexos. */

    double complex *vetor;
    vetor = (double complex*) calloc(N, sizeof(double complex));
    return vetor;
}

void imprimir_vetor(double* vetor, int N) {
    /* Impressao de vetor de doubles */

    for(int i = 0; i < N; i++) {
        if(vetor[i] >= 0) {
            printf("|  %.3e |\n", vetor[i]);
        }
        else {
            printf("| %.3e |\n", vetor[i]);
        }
    }
    printf("\n");
}

void imprimir_complexo(double complex *c, int N) { 
	/* Impressao de vetor de complexos */

    for(int i = 0; i < N; i++) {
    	if (creal(c[i]) > 0) {
    		printf("| %.2f %+.2fi|\n", creal(c[i]), cimag(c[i]));
    	} 
    	else {
    		printf("|%.2f %+.2fi|\n", creal(c[i]), cimag(c[i]));
    	}  
    }
    printf("\n");
}

int tamanho_arquivo(char *nome_arquivo, int *canais) {
	int n = -1;
	int var_temp1;
	double var_temp2, var_temp3, var_temp4;
	char linha[512];

	FILE *arquivo = fopen(nome_arquivo, "r");

	if(arquivo == NULL) {
        printf("\nArquivo nao encontrado\n");
        exit(EXIT_FAILURE);
    }

    // Passando pelas duas primeiras linhas
    fscanf(arquivo, "%*[^0-9]%d", &var_temp1);
    fscanf(arquivo, "%*[^0-9]%d", &var_temp1);

    *canais = var_temp1; // retornando numero de canais por meio de um ponteiro

    if(*canais == 1) {
    	while(fgets(linha, sizeof(linha), arquivo) != NULL) { // pega uma linha de até 512 caracteres. Null quando acabar as linhas 

	        sscanf(linha, "%lf %lf", &var_temp2, &var_temp3);
	        n++;  // contando o numero de dados
	        
        }
    }

    else if(*canais == 2) {
    	while(fgets(linha, sizeof(linha), arquivo) != NULL) { // pega uma linha de até 512 caracteres. Null quando acabar as linhas

	        sscanf(linha, "%lf %lf %lf", &var_temp2, &var_temp3, &var_temp4);
	        n++;
	       
	    }
    }

    else {
        printf("Numero de canais nao suportado para analise\n");
    }

    fclose(arquivo);
    return n;

}


int ler_arquivo_dat(char *nome_arquivo, int n, double complex *x, double complex *f, double complex *f2, double *f_linha, double *f_linha2) {
	int sample_rate, channels;
	double var_x, var_f, var_f2;
	char linha[512];

	FILE *arquivo = fopen(nome_arquivo, "r");

    if(arquivo == NULL) {
        printf("\nArquivo nao encontrado\n");
        exit(EXIT_FAILURE);
    }

    /* Solucao para fscanf encontrada em: https://stackoverflow.com/questions/19413569/can-i-use-fscanf-to-get-only-digits-from-text-that-contain-chars-and-ints */
    fscanf(arquivo, "%*[^0-9]%d", &sample_rate);  // Sample Rate
    fscanf(arquivo, "%*[^0-9]%d", &channels);  // Channels

    /* leitura de dados do arquivo */
    if(channels == 1) {

    	for (int i = 0; i < n; i++) {
	        fscanf(arquivo, "%lf", &x[i]);
	        fscanf(arquivo, "%lf", &f[i]);
	        f_linha[i] = f[i]; // amplitude do sinal - para fftpack4
    	}
    }

    else if(channels == 2) {
    	for (int i = 0; i < n; i++) {
	        fscanf(arquivo, "%lf", &x[i]);
	        fscanf(arquivo, "%lf", &f[i]);
	        fscanf(arquivo, "%lf", &f2[i]);
	        f_linha[i] = f[i]; // amplitude do sinal - para fftpack4
	        f_linha2[i] = f2[i]; // amplitude do sinal - para fftpack4
    	}
    }

    else {
        printf("Numero de canais nao suportado para analise\n");
    }

    fclose(arquivo);
    return sample_rate;
}


void escrever_arquivo_dat(char *nome_arquivo, int sample_rate, int canais, int n, double complex *x, double complex *f, double complex *f2) {
	int i;
	FILE *arquivo = fopen(nome_arquivo, "w");

	// Escreve a sample rate e o numero de canais no mesmo formato dos arquivos .dat fornecidos
	fprintf(arquivo, "; Sample Rate %d\n", sample_rate);
	fprintf(arquivo, "; Channels %d\n", canais);

	// Se houver apenas 1 canal, havera 2 colunas de dados
	if(canais == 1) {
		for(i = 0; i < n; i++) {
			fprintf(arquivo, "%.13lf %.13lf\n", creal(x[i]), creal(f[i]));
		}
	}

	// Se houverem 2 canais, havera 3 colunas de dados
	else if(canais == 2) {
		for(i = 0; i < n; i++) {
			fprintf(arquivo, "%.13lf %.13lf %.13lf\n", creal(x[i]), creal(f[i]), creal(f2[i]));
		}
	}
	else {
		printf("Numero de canais nao suportado.\n");
	}
	fclose(arquivo);
}


/* >>>>>>>>>>>>>>>>>>>>>>>> Funcoes de Transformada de Fourier <<<<<<<<<<<<<<<<<<<<<<<< */

void fourier(double complex *c, double complex *f, double complex *x, int n){ 
	/*Obtem um vetor de coeficientes c como resultado da Transformada de Fourier 
	direta de f(x)*/

    double complex somatorio;
    
    for (int k = 0; k < 2*n; k++){
        somatorio = 0;
        for (int j = 0; j < 2*n; j++){
            somatorio += f[j]* cexp(-I * k * x[j]);
        }
        c[k] = somatorio * (double)1 / (2*n);
    }
}


void anti_fourier(double complex *c, double complex *f, double complex *x, int n){
	/*Obtem um vetor f com a funcao antitransformada a partir dos coeficientes do vetor c
	e do vetor x*/

    double complex somatorio;

    for (int j = 0; j < 2*n; j++){
        somatorio = 0;
        for (int k = 0; k < 2*n; k++){
            somatorio += c[k] * cexp(I * k * x[j]);
        }
        f[j] = somatorio;
    }
}


void fftrec(double complex *c, double complex *f, int n, bool dir) {
	/* Funcao feita com base no pseudo-codigo do enunciado do EP2 */

	double complex eij;
	double complex *even, *odd, *fe, *fo;

	even = criar_vetor_complexo(n);
	odd = criar_vetor_complexo(n);
	fe = criar_vetor_complexo(n);
	fo = criar_vetor_complexo(n);

	if(n == 1) {
		c[0] = f[0] + f[1];
		c[1] = f[0] - f[1];
	}
	else {
		for(int j = 0; j < n; j++) {
			fe[j] = f[2 * j];
			fo[j] = f[2 * j + 1];
		}

		fftrec(even, fe, n/2, dir);
		fftrec(odd, fo, n/2, dir);

		for(int j = 0; j < n; j++) {
			if(dir) {
				eij = cexp(- I * j * M_PI / n);  // https://stackoverflow.com/questions/2834865/computing-e-j-in-c
			}
			else {
				eij = cexp(I * j * M_PI / n);
			}

			c[j] = even[j] + eij * odd[j];
			c[j+n] = even[j] - eij * odd[j];
		}
	}
}


void comprimir_sinal(double complex *c, double taxa_minima, int N) {
	/* Zera as frequencias menos significativas do sinal, definidas pelo
	dobro da amplitude de cada frequencia comparada a uma taxa minima */

    double amplitude = 0;

    for (int k = 0; k < N; k++){
        amplitude = 2 * cabs(c[k]);
        if (amplitude < taxa_minima){
            c[k] = 0;
        }
    }
}


void passa_altas(double complex *c, int N, int freq_corte ) {
    /* Filtro que zera todas as frequencias abaixo da frequencia de corte */

    for (int k = freq_corte-1; k > 0; k--) {
        c[k] = 0;
    } 
}


void passa_baixas(double complex *c, int N, int freq_corte) {
	/* Filtro que zera todas as frequencias acima da frequencia de corte */

    for (int k = freq_corte+1; k < N; k++) {
        c[k] = 0;
    } 
}


void passa_bandas(double complex *c, int N, int freq1, int freq2) {
	/* Filtro que zera todas as frequencias fora de uma faixa de frequencias */

    for (int k = freq1-1; k > 0; k--) {
        c[k] = 0;
    }

    for (int k = freq2+1; k < N; k++) {
        c[k] = 0;
    }
}

