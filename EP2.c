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
void imprimir_complexo_R(double complex *c, int N);

/* >>>>>>>>>>> Funcoes de Transformada de Fourier <<<<<<<<<<< */
int contar_amostras(char *nome_arquivo, int *sample_rate, int *canais, int *tratar_dados, int *delta);
void ler_arquivo(char *nome_arquivo, int n, double complex *x, double complex *f, double complex *f2, double *f_linha, double *f2_linha, int tratar_dados, int delta);
void escrever_arquivo(char *nome_arquivo, int sample_rate, int canais, int n, double complex *x, double complex *f, double complex *f2);
void escrever_arquivo_linha(char *nome_arquivo, int sample_rate, int canais, int n, double complex *x, double *f_linha, double *f2_linha);
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
	double complex *f_rec; // Vetor de valores da função reconstituida pela antitransformada
	int tipo_problema, tipo_transformada;
	bool direta; // Usado na fftrec
	clock_t tempo[2]; // Usado para medir o tempo de execucao
	double tempo_execucao;

	// Variaveis especificas dos arquivos de audio
	int escolha, tipo_filtro, tipo_compressao;
	char nome_arquivo[128];
	double complex *c2; // Vetor de coeficientes para o segundo canal
    double complex *f2; // Vetor de valores da função para o segundo canal
    double complex *f2_rec; // Vetor de valores da função reconstituida pela antitransformada
    int canais, sample_rate; // Parametros fornecidos pelo arquivo
    int K, K1, K2; // Parametros de corte utilizados nos filtros
    int S; // Parametro de compressao
    int tratar_dados; // Para o caso do numero de amostras nao ser uma potencia de 2
    int delta; // Parametro de tratamento de dados

	// Variaveis especificas do fftpack4
	double *a, *b, *a2, *b2;
	double a0, a02; 
	int *ifac, *ifac2;
	double *wsave, *wsave2;
	double *f_linha, *f2_linha; 

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
        	n = 4;

			x = criar_vetor_complexo(n);
			x[0] = 0;
   			x[1] = M_PI / 2;
    		x[2] = M_PI;
    		x[3] = 3 * M_PI / 2;

    		f = criar_vetor_complexo(n);
    		f[0] = 5;
		    f[1] = -1;
		    f[2] = 3;
		    f[3] = 1;

		    f_linha = criar_vetor(n);  // No FFTPACK4, trabalha-se com valores em double e nao complexos:
		    f_linha[0] = 5;
		    f_linha[1] = -1;
		    f_linha[2] = 3;
		    f_linha[3] = 1;
            break;

        case 2: // Teste b
        	n = 8;

			x = criar_vetor_complexo(n);
			x[0] = 0;
   			x[1] = M_PI / 4;
    		x[2] = M_PI / 2;
    		x[3] = 3 * M_PI / 4;
    		x[4] = M_PI;
    		x[5] = 5 * M_PI / 4;
    		x[6] = 3 * M_PI / 2;
    		x[7] = 7 * M_PI / 4;

			f = criar_vetor_complexo(n);
			f[0] = 6;
		    f[1] = 2;
		    f[2] = 5;
		    f[3] = 2;
		    f[4] = 11;
		    f[5] = 2;
		    f[6] = 8;
		    f[7] = 8;

		    f_linha = criar_vetor(n); // No FFTPACK4, trabalha-se com valores em double e nao complexos
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
        	n = 1024;
			x = criar_vetor_complexo(n);
			f = criar_vetor_complexo(n);
			f_linha = criar_vetor(n); // No FFTPACK4, trabalha-se com valores em double e nao complexos
			
			for(int i = 0; i < n; i++) {
		        x[i] = (i * M_PI) / 512; 

		        f[i] = 10*sin(x[i]) + 7*cos(30*x[i]) + 11*sin(352*x[i]) - 8*cos(711*x[i]);
		        f_linha[i] = 10*sin(x[i]) + 7*cos(30*x[i]) + 11*sin(352*x[i]) - 8*cos(711*x[i]);
		    }
            break;

        case 4: // Arquivos de audio
        	//Destroi os ponteiros que nao serao utilizados
        	f_rec = NULL;
        	f2_rec = NULL;
		    printf("Digite o nome do arquivo a ser analizado (com a terminacao .dat): ");
		    scanf("%s", nome_arquivo);
		    n = contar_amostras(nome_arquivo, &sample_rate, &canais, &tratar_dados, &delta); // Define o numero de amostras, canais e o tratamento dos dados

		    switch(tratar_dados){
				case 0:
					printf("\nO sinal analisado tem %d amostras, que corresponde a uma potencia de 2.\n", n);
					break;
				case 1:
					printf("\nO sinal analisado tem %d amostras, que nao corresponde a uma potencia de 2, tendo %d amostras a mais da potencia mais proxima.\n", n, delta);
					break;
				case 2:
				printf("\nO sinal analisado tem %d amostras, que nao corresponde a uma potencia de 2, tendo %d amostras a menos da proxima potencia.\n", n, delta);
					break;
				default:
					break;
			}

		    printf("\n---------------Analise de Fourier----------------\n");
			printf("1 - Forma direta\n");
			printf("2 - FFT Recursiva\n");
			printf("3 - FFTPACK4\n");
			printf("Selecione o tipo de transformada desejada: ");
			scanf("%d", &tipo_transformada);

			// Tratamento de dados
			if(tipo_transformada != 3 && tratar_dados == 1){
				n -= delta; // As amostras em excesso serao cortadas
			} else if(tipo_transformada != 3 && tratar_dados == 2) {
				n += delta; // As amostras faltantes serao preenchidas com a media das demais
			}

			// Caso n nao seja uma potencia de 2, mas o metodo seja fftpack4, nao ha tratamento de dados
			if(tratar_dados != 0 && tipo_transformada == 3){
				tratar_dados = 0;
			}

			// Alocacao de vetores
			x = criar_vetor_complexo(n);
			f = criar_vetor_complexo(n);
			f2 = criar_vetor_complexo(n);
			// No FFTPACK4, trabalha-se com valores em double e nao complexos
			f_linha = criar_vetor(n); 
			f2_linha = criar_vetor(n);
	    	c = criar_vetor_complexo(n);
	    	if(canais == 2) {
				c2 = criar_vetor_complexo(n);
			} else {
				//Desaloca-se os ponteiros que nao serao utilizados
				c2 = NULL;
				f2 = NULL;
			}
			// Preenchimento dos vetores de acordo com o arquivo
			ler_arquivo(nome_arquivo, n, x, f, f2, f_linha, f2_linha, tratar_dados, delta);
			
			printf("\nO sinal analisado tem %d amostras. \n", n);	
            break;

        default:
            printf("Opcao invalida!\n");
            exit(EXIT_FAILURE);
            break;
    }

	if(tipo_problema != 4) { // Testes a, b, c:
		//Destroi os ponteiros que nao serao utilizados
		c2 = NULL;
		f2 = NULL;
		f2_linha = NULL;
		ifac2 = NULL;
		wsave2 = NULL;
		a2 = NULL;
		b2 = NULL;
		printf("\nO sinal analisado tem %d amostras. \n", n);	
		c = criar_vetor_complexo(n);


		printf("\n------------- Forma Direta ---------------\n\n");
		printf("Amostra original do sinal:\n");
		imprimir_complexo_R(f, n);

		printf("Transformada de Fourier - Vetor de coeficientes:\n");
		tempo[0] = clock();
		fourier(c, f, x, n);
		tempo[1] = clock();
		tempo_execucao = (tempo[1] - tempo[0]) * 1000.0 / CLOCKS_PER_SEC;
		imprimir_complexo(c, n);
		printf("Tempo gasto: %g ms.\n\n", tempo_execucao);

	    f_rec = criar_vetor_complexo(n); // Inicializa vetor para o sinal reconstituido
	    printf("Antitransformada de Fourier - Sinal recuperado:\n");
	    tempo[0] = clock();
	    anti_fourier(c, f_rec, x, n);
	    tempo[1] = clock();
	    tempo_execucao = (tempo[1] - tempo[0]) * 1000.0 / CLOCKS_PER_SEC;
		imprimir_complexo_R(f_rec, n);

		// Libera a memoria para o uso na FFT Recursiva
		free(f_rec);
		free(c);
		printf("Tempo gasto: %g ms.\n", tempo_execucao);

		printf("\n------------- FFT Recursiva ---------------\n\n");
		printf("Amostra original do sinal:\n");
	    imprimir_complexo_R(f, n);

		printf("Transformada de Fourier - Vetor de coeficientes:\n");
		c = criar_vetor_complexo(n); //Reinicializa o vetor de coeficientes
		tempo[0] = clock();
		fftrec(c, f, n/2, true);  // Transformada direta pela fftrec
		
		// Como indicado no algoritmo, faz-se necessaria a divisao dos coeficientes por N
		for(int i = 0; i < n; i++){
	        c[i] = c[i]/n;
	    }
	    tempo[1] = clock();
	    tempo_execucao = (tempo[1] - tempo[0]) * 1000.0 / CLOCKS_PER_SEC;
	    imprimir_complexo(c, n);
	    printf("Tempo gasto: %g ms.\n\n", tempo_execucao);
	    
	    printf("Antitransformada de Fourier - Sinal recuperado:\n");
	    tempo[0] = clock();
	    f_rec = criar_vetor_complexo(n); // Reinicializa vetor para o sinal reconstituido
		fftrec(c, f_rec, n/2, false);  // Anti-transformada pela fftrec
		tempo[1] = clock();	
		tempo_execucao = (tempo[1] - tempo[0]) * 1000.0 / CLOCKS_PER_SEC;
	    imprimir_complexo_R(f_rec, n);
	    free(c); // Desaloca para o uso no FFTPACK4
	    // Nao serao mais utilizados
	    free(f_rec);
	    free(f);
	    printf("Tempo gasto: %g ms.\n", tempo_execucao);
	  
		printf("\n------------- FFTPACK4 ---------------\n\n");
	    printf("Amostra original do sinal:\n");
	    imprimir_vetor(f_linha, n);
		tempo[0] = clock();
		wsave = criar_vetor(3 * n + 15);
		ifac = criar_vetor_int(8);
		a = criar_vetor(n/2);
		b = criar_vetor(n/2);

		ezffti(&n, wsave, ifac);  // inicializacao da fftpack4
		ezfftf(&n, f_linha, &a0, a, b, wsave, ifac);  // transformada direta de fourier

		// Conversao de valores do tipo a*cos() + b*sen() para coeficientes complexos do tipo ck
		c = criar_vetor_complexo(n); //Reinicializa o vetor de coeficientes
		c[0] = a0;
		for(int i = 1; i < (n/2); i++) {
			c[i] = (a[i-1] - (I * b[i-1]))/2;	
		}
		c[n/2] = a[(n/2)-1] + I*b[(n/2)-1];
		for(int i = 1; i < (n/2); i++) {
			c[n-i] = (a[i-1] + I * b[i-1])/2;
		}
		tempo[1] = clock();
		tempo_execucao = (tempo[1] - tempo[0]) * 1000.0 / CLOCKS_PER_SEC;
		printf("Transformada de Fourier - Vetor de coeficientes:\n");
	    imprimir_complexo(c, n);
	    printf("Tempo gasto: %g ms.\n\n", tempo_execucao);

	    free(f_linha); // Limpando o vetor para usar novamente na antitransformada
	    f_linha = criar_vetor(n); // Reinicializando o vetor
	    printf("Antitransformada de Fourier - Sinal recuperado:\n");
	    tempo[0] = clock();
	    // Conversao de valores ck para ak, bk
	    a0 = c[0];
		for(int k = 1; k < (n/2); k++) {
			a[k-1] = c[k] + c[n-k];
			b[k-1] = I*(c[k] - c[n-k]);
		}
		ezfftb(&n, f_linha, &a0, a, b, wsave, ifac);  // Antitransformada de fourier
	    tempo[1] = clock();
	    tempo_execucao = (tempo[1] - tempo[0]) * 1000.0 / CLOCKS_PER_SEC;
	    imprimir_vetor(f_linha, n);
	    printf("Tempo gasto: %g ms.\n", tempo_execucao);

	    // Desalocacao de memoria
	    free(ifac);
		free(wsave);
		free(a);
		free(b);
		free(c);
		free(x);
		free(f_linha);
		ifac = NULL;
		wsave = NULL;
		a = NULL;
		b = NULL;
		c = NULL;
		x = NULL;
		f = NULL;
		f_linha = NULL;
		f_rec = NULL;
	    return 0;

        
    } else { // Para os arquivos de audio
    	
		// Transformada de Fourier

		switch(tipo_transformada) {
			case 1:
				//Desaloca os vetores que nao serao utilizados
				free(f_linha);
				free(f2_linha);

				//Destroi os ponteiros que nao serao utilizados
				ifac = NULL;
				wsave = NULL;
				a = NULL;
				b = NULL;
				ifac2 = NULL;
				wsave2 = NULL;
				a2 = NULL;
				b2 = NULL;

				tempo[0] = clock();
				fourier(c, f, x, n);
				free(f); // Limpa o vetor para usar novamente na antitransformada
				if(canais == 2) {
					fourier(c2, f2, x, n);
					free(f2); // Limpa o vetor para usar novamente na antitransformada
				}
				tempo[1] = clock();
				tempo_execucao = (tempo[1] - tempo[0]) * 1000.0 / CLOCKS_PER_SEC;
				printf("\nTransformada de Fourier realizada em %g ms.\n\n", tempo_execucao);
				break;

			case 2:
				//Desaloca os vetores que nao serao utilizados
				free(f_linha);
				free(f2_linha);

				//Destroi os ponteiros que nao serao utilizados
				ifac = NULL;
				wsave = NULL;
				a = NULL;
				b = NULL;
				ifac2 = NULL;
				wsave2 = NULL;
				a2 = NULL;
				b2 = NULL;

				tempo[0] = clock();
				fftrec(c, f, n/2, true);
				free(f); // Limpa o vetor para usar novamente na antitransformada
				if(canais == 2) {
					fftrec(c2, f2, n/2, true);
					free(f2); // Limpa o vetor para usar novamente na antitransformada
				}
				// Como indicado no algoritmo, faz-se necessaria a divisao dos coeficientes por 2N
				for(int i = 0; i < n; i++){
			        c[i] = c[i]/n;
			    }
			    if(canais == 2) {
			    	for(int i = 0; i < n; i++){
			        	c2[i] = c[i]/n;
			   		}
			    }
			    tempo[1] = clock();
			    tempo_execucao = (tempo[1] - tempo[0]) * 1000.0 / CLOCKS_PER_SEC;
				printf("\nTransformada de Fourier (Forma Recursiva) realizada em %g ms.\n\n", tempo_execucao);
			    break;

			case 3:
				//Desaloca-se os vetores que nao serao utilizados
				free(f);
				free(f2);

				tempo[0] = clock();
				wsave = criar_vetor(3 * n + 15);
				ifac = criar_vetor_int(8);
				a = criar_vetor(n/2);
				b = criar_vetor(n/2);
				if(canais == 2) {
					wsave2 = criar_vetor(3 * n + 15);
					ifac2 = criar_vetor_int(8);
					a2 = criar_vetor(n/2);
					b2 = criar_vetor(n/2);
				} else {
					//Desaloca-se os ponteiros que nao serao utilizados
					ifac2 = NULL;
					wsave2 = NULL;
					a2 = NULL;
					b2 = NULL;
				}
				
				ezffti(&n, wsave, ifac);  // inicializacao da fftpack4
				ezfftf(&n, f_linha, &a0, a, b, wsave, ifac);  // transformada
				free(f_linha); // Limpa o vetor para usar novamente na antitransformada
				
				// Conversao de valores do tipo a*cos() + b*sen() para coeficientes complexos do tipo ck
				c[0] = a0;
				for(int i = 1; i < (n/2); i++) {
					c[i] = (a[i-1] - (I * b[i-1]))/2;	
				}
				c[n/2] = a[(n/2)-1] + I*b[(n/2)-1];
				for(int i = 1; i < (n/2); i++) {
					c[n-i] = (a[i-1] + I * b[i-1])/2;
				}

				if(canais == 2) {
			
					ezffti(&n, wsave2, ifac2);  // inicializacao da fftpack4
					ezfftf(&n, f2_linha, &a02, a2, b2, wsave2, ifac2);  // transformada direta de fourier
					free(f2_linha); // Limpa o vetor para usar novamente na antitransformada
					// Conversao de valores do tipo a*cos() + b*sen() para coeficientes complexos do tipo ck
					c2[0] = a02;
					for(int i = 1; i < (n/2); i++) {
						c2[i] = (a2[i-1] - (I * b2[i-1]))/2;	
					}
					c2[n/2] = a2[(n/2)-1] + I*b2[(n/2)-1];
					for(int i = 1; i < (n/2); i++) {
						c2[n-i] = (a2[i-1] + I * b2[i-1])/2;
					}
				}
				tempo[1] = clock();
				tempo_execucao = (tempo[1] - tempo[0]) * 1000.0 / CLOCKS_PER_SEC;
				printf("\nTransformada de Fourier (FFTPACK4) realizada em %g ms.\n\n", tempo_execucao);
				break;
			default:
				printf("Opcao invalida!\n");
				exit(EXIT_FAILURE);
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
				tempo[0] = clock();
				f = criar_vetor_complexo(n); // Reinicializa o vetor para obter o sinal reconstituido
			    anti_fourier(c, f, x, n);
			    //Desalocando memoria
			    free(c);
			    c = NULL;
			    if(canais == 2){
			    	f2 = criar_vetor_complexo(n); // Reinicializa o vetor para obter o sinal reconstituido
			    	anti_fourier(c2, f2, x, n);
			    	//Desalocando memoria
			    	free(c2);
			    	c2 = NULL;
			    }
			    tempo[1] = clock();
			    tempo_execucao = (tempo[1] - tempo[0]) * 1000.0 / CLOCKS_PER_SEC;
			    printf("\nAntitransformada de Fourier aplicada em %g ms.\n\n", tempo_execucao);
				break;

			case 2:
				tempo[0] = clock();
				f = criar_vetor_complexo(n); // Reinicializa o vetor para obter o sinal reconstituido
				fftrec(c, f, n/2, false);
				//Desalocando memoria
			    free(c);
			    c = NULL;
				if(canais == 2){
					f2 = criar_vetor_complexo(n); // Reinicializa o vetor para obter o sinal reconstituido
			    	fftrec(c2, f2, n/2, false);
			    	//Desalocando memoria
			    	free(c2);
			    	c2 = NULL;
			    }
			    tempo[1] = clock();	
			    tempo_execucao = (tempo[1] - tempo[0]) * 1000.0 / CLOCKS_PER_SEC;
				printf("\nAntitransformada de Fourier (FFT Recursiva) aplicada em %g ms.\n\n", tempo_execucao);
			    break;

			case 3:
				tempo[0] = clock();
				f_linha = criar_vetor(n); // Reinicializa o vetor para obter o sinal reconstituido
				// Conversao do vetor de coeficientes c para a0, a b:
				a0 = c[0];
				for(int k = 1; k < (n/2); k++) {
					a[k-1] = c[k] + c[n-k];
					b[k-1] = I*(c[k] - c[n-k]);
				}

				//Desalocando memoria
			    free(c);
			    c = NULL;

				ezfftb(&n, f_linha, &a0, a, b, wsave, ifac); 
				//Desalocando memoria
				free(a);
				free(b);
				free(wsave);
				free(ifac);
				a = NULL;
				b = NULL;
				wsave = NULL;
				ifac = NULL;

			    if(canais == 2){
			    	f2_linha = criar_vetor(n); // Reinicializa o vetor para obter o sinal reconstituido
			    	// Conversao do vetor de coeficientes c para a0, a b:
					a02 = c2[0];
					for(int k = 1; k < (n/2); k++) {
						a2[k-1] = c2[k] + c2[n-k];
						b2[k-1] = I*(c2[k] - c2[n-k]);
					}
			    	ezfftb(&n, f2_linha, &a02, a2, b2, wsave2, ifac2);
			    	//Desalocando memoria
					free(a2);
					free(b2);
					free(wsave2);
					free(ifac2);
					free(c2);
			    	c2 = NULL;
					a2 = NULL;
					b2 = NULL;
					wsave2 = NULL;
					ifac2 = NULL;  
			    }	
			    tempo[1] = clock();
			    tempo_execucao = (tempo[1] - tempo[0]) * 1000.0 / CLOCKS_PER_SEC;
			    printf("\nAntitransformada de Fourier (FFTPACK4) aplicada em %g ms.\n\n", tempo_execucao);
			    break;
        }


        //Gravando arquivo
    	printf("\nDigite o nome do arquivo para gravar o resultado (com a terminacao .dat): ");
    	scanf("%s", nome_arquivo);
    	// Escreve o resultado da analise no arquivo
    	if(tipo_transformada != 3) {
    		escrever_arquivo(nome_arquivo, sample_rate, canais, n, x, f, f2);
    		// Desalocacao de memoria
    		free(f);
    		free(x);
    		if(canais == 2) {
    			free(f2);
    		}
    	}
    	else {
    		escrever_arquivo_linha(nome_arquivo, sample_rate, canais, n, x, f_linha, f2_linha);
    		free(f_linha);
    		free(x);
    		if(canais == 2) {
    			free(f2_linha);
    		}
    	}
    	printf("\nArquivo gravado com sucesso!\n");
	    
    }

    // Desalocacao de ponteiros
	x = NULL;
	f = NULL;
	f2 = NULL;
	f_linha = NULL;
	f2_linha = NULL;
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
            printf("| %.2e |\n", vetor[i]);
        }
        else {
            printf("|%.2e |\n", vetor[i]);
        }
    }
    printf("\n");
}

void imprimir_complexo_R(double complex *c, int N) { 
	/* Impressao de vetor de complexos */

    for(int i = 0; i < N; i++) {
    	if (creal(c[i]) > 0) {
    		printf("| %.2e |\n", creal(c[i]));
    	} 
    	else {
    		printf("|%.2e |\n", creal(c[i]));
    	}  
    }
    printf("\n");
}

void imprimir_complexo(double complex *c, int N) { 
	 // Impressao de vetor de complexos 

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

int contar_amostras(char *nome_arquivo, int *sample_rate, int *canais, int *tratar_dados, int *delta) {
	int n = 0;
	int var_temp1;
	double var_temp2, var_temp3, var_temp4;
	char linha[512];
	int i, potencia, delta_prev, delta_next;

	FILE *arquivo = fopen(nome_arquivo, "r");

	if(arquivo == NULL) {
        printf("\nArquivo nao encontrado\n");
        exit(EXIT_FAILURE);
    }

    // Passando pelas duas primeiras linhas
    fscanf(arquivo, "%*[^0-9]%d", &var_temp1);
    *sample_rate = var_temp1;
    fscanf(arquivo, "%*[^0-9]%d", &var_temp1);
    *canais = var_temp1; 

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
    n-=1;

    //Verifica se o numero de amostras e uma potencia de 2
    i = 0;
    potencia = pow(2, i);
    while (potencia < n){ // Itera a potencia de 2 ate que seja maior ou igual a n
    	i++;
    	potencia = pow(2, i);
    }

    if(potencia != n) {
    	delta_next = potencia - n; // Distancia de n ate a proxima potencia de 2
    	delta_prev = n - pow(2, i-1); // Distancia de n ate a potencia de 2 anterior
    	if (delta_prev < 50) { // Valor arbitrado para o numero maximo de amostras que podem ser descartadas
    		*tratar_dados = 1; // Descartar as amostras que estao sobrando
    		*delta = delta_prev;
    	} else {
    		*tratar_dados = 2; // Completar as amostras com a media dos valores ate n = proxima potencia de 2
    		*delta = delta_next;
    	}
    } else {
    	*tratar_dados = 0; //Nenhum tratamento de dado necessario, pois n = potencia de 2
    }
	
    return n;
}


void ler_arquivo(char *nome_arquivo, int n, double complex *x, double complex *f, double complex *f2, double *f_linha, double *f2_linha, int tratar_dados, int delta) {
	int sample_rate, channels;
	double soma, soma2, media, media2, delta_x;
	char linha[512];
	int i;

	FILE *arquivo = fopen(nome_arquivo, "r");

    if(arquivo == NULL) {
        printf("\nArquivo nao encontrado\n");
        exit(EXIT_FAILURE);
    }

    /* Solucao para fscanf encontrada em: https://stackoverflow.com/questions/19413569/can-i-use-fscanf-to-get-only-digits-from-text-that-contain-chars-and-ints */
    fscanf(arquivo, "%*[^0-9]%d", &sample_rate);  // Sample Rate
    fscanf(arquivo, "%*[^0-9]%d", &channels);  // Channels

	if(tratar_dados == 2) { // Completa-se o vetor com a media das amostras
		 
		if(channels == 1) {
			soma = 0;
	    	for (i = 0; i < (n-delta); i++) {
		        fscanf(arquivo, "%lf", &x[i]);
		        fscanf(arquivo, "%lf", &f[i]);
		        f_linha[i] = f[i]; // amplitude do sinal - para fftpack4
		        soma += f[i];
	    	}
	    	media = soma/(n-delta); // Media dos valores das amostras;
	    	delta_x = x[i] - x[i-1]; // Armazena a distancia entre duas amostras

	    	for (i = (n-delta); i < n; i++) {
	    		x[i] = x[i-1] + delta_x; // Preenche o vetor de amostras com o passo definido por delta_x
	    		f[i] = media;
	    		f_linha[i] = media;
	    	}
	    }

	    else if(channels == 2) {
	    	soma = 0;
	    	soma2 = 0;
	    	for (i = 0; i < (n-delta); i++) {
		        fscanf(arquivo, "%lf", &x[i]);
		        fscanf(arquivo, "%lf", &f[i]);
		        fscanf(arquivo, "%lf", &f2[i]);
		        f_linha[i] = f[i]; // amplitude do sinal - para fftpack4
		        f2_linha[i] = f2[i]; // amplitude do sinal - para fftpack4
		        soma += f[i];
		        soma2 += f2[i];
	    	}

	    	media = soma/(n-delta); // Media dos valores das amostras;
	    	media2 = soma2/(n-delta); // Media dos valores das amostras para o canal 2;
	    	delta_x = x[i] - x[i-1]; // Armazena a distancia entre duas amostras

	    	for (i = (n-delta); i < n; i++) {
	    		x[i] = x[i-1] + delta_x; // Preenche o vetor de amostras com o passo definido por delta_x
	    		f[i] = media;
	    		f_linha[i] = media;
	    		f2[i] = media2;
	    		f2_linha[i] = media2;
	    	}
	    }

	    else {
	        printf("Numero de canais nao suportado para analise\n");
	        exit(EXIT_FAILURE);
	    }

	} else { // Sem tratamento de dados
		
	    if(channels == 1) {

	    	for (i = 0; i < n; i++) {
		        fscanf(arquivo, "%lf", &x[i]);
		        fscanf(arquivo, "%lf", &f[i]);
		        f_linha[i] = f[i]; // amplitude do sinal - para fftpack4
	    	}
	    }

	    else if(channels == 2) {
	    	for (i = 0; i < n; i++) {
		        fscanf(arquivo, "%lf", &x[i]);
		        fscanf(arquivo, "%lf", &f[i]);
		        fscanf(arquivo, "%lf", &f2[i]);
		        f_linha[i] = f[i]; // amplitude do sinal - para fftpack4
		        f2_linha[i] = f2[i]; // amplitude do sinal - para fftpack4
	    	}
	    }

	    else {
	        printf("Numero de canais nao suportado para analise\n");
	        exit(EXIT_FAILURE);
	    }
	}

    fclose(arquivo);
}


void escrever_arquivo(char *nome_arquivo, int sample_rate, int canais, int n, double complex *x, double complex *f, double complex *f2) {
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

void escrever_arquivo_linha(char *nome_arquivo, int sample_rate, int canais, int n, double complex *x, double *f_linha, double *f2_linha) {
	int i;
	FILE *arquivo = fopen(nome_arquivo, "w");

	// Escreve a sample rate e o numero de canais no mesmo formato dos arquivos .dat fornecidos
	fprintf(arquivo, "; Sample Rate %d\n", sample_rate);
	fprintf(arquivo, "; Channels %d\n", canais);

	// Se houver apenas 1 canal, havera 2 colunas de dados
	if(canais == 1) {
		for(i = 0; i < n; i++) {
			fprintf(arquivo, "%.13lf %.13lf\n", creal(x[i]), f_linha[i]);
		}
	}

	// Se houver 2 canais, havera 3 colunas de dados
	else if(canais == 2) {
		for(i = 0; i < n; i++) {
			fprintf(arquivo, "%.13lf %.13lf %.13lf\n", creal(x[i]), f_linha[i], f2_linha[i]);
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
    
    for (int k = 0; k < n; k++){
        somatorio = 0;
        for (int j = 0; j < n; j++){
            somatorio += f[j]* cexp(-I * k * x[j]);
        }
        c[k] = somatorio / n;
    }
}


void anti_fourier(double complex *c, double complex *f, double complex *x, int n){
	/*Obtem um vetor f com a funcao antitransformada a partir dos coeficientes do vetor c
	e do vetor x*/

    double complex somatorio;

    for (int j = 0; j < n; j++){
        somatorio = 0;
        for (int k = 0; k < n; k++){
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

	if (dir) {
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
				eij = cexp(- I * j * M_PI / n);
				c[j] = even[j] + eij * odd[j];
				c[j+n] = even[j] - eij * odd[j];
			}
		}
	} else {
		if(n == 1) {
			f[0] = c[0] + c[1];
			f[1] = c[0] - c[1];
		}
		else {
			for(int j = 0; j < n; j++) {
				fe[j] = c[2 * j];
				fo[j] = c[2 * j + 1];
			}

			fftrec(fe, even, n/2, dir);
			fftrec(fo, odd, n/2, dir);

			for(int j = 0; j < n; j++) {
				eij = cexp(I * j * M_PI / n);
				f[j] = even[j] + eij * odd[j];
				f[j+n] = even[j] - eij * odd[j];
			}
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

    for (int k = freq_corte-1; k >= 0; k--) {
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

