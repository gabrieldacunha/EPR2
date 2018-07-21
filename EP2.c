/* EP2 - Metodos Numericos - MAP3121
 *
 * Gabriel da Cunha Rodrigues - No USP: 8992930 - Turma 2
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#include "fftpack4.h"

/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Declaracao de funcoes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

/* >>>>>>>>>>> Funcoes basicas <<<<<<<<<<< */
double* criarVetor(int N);
int* criarVetorInt(int N);
double complex* criarVetorComplexo(int N);
void teste_complexos();
int tamanho_arquivo(char *nome_arquivo, int *channels);
int ler_arquivo_dat(char *nome_arquivo, int N, double complex *tempo, double complex *f, double complex *f2);
void escrever_arquivo_dat(char *nome_arquivo, int sample_rate, int channels, int N, double complex *tempo, double complex *f, double complex *f2);

/* >>>>>>>>>>> Funcoes de Transformada de Fourier <<<<<<<<<<< */
void fftrec(double complex *c, double complex *f, int n, bool dir);
void filtro(double complex *c, int N, double freq_corte, char *tipo);

/* >>>>>>>>>>> Testes iniciais <<<<<<<<<<< */
void teste_inicial_a();
void teste_inicial_b();
void teste_inicial_c();

/* ------------------------------------------------------------------------------------- */

int main() {
	// Declaracao de variaveis
	int N, tipo_problema, tipo_transformada, escrever_arquivo;
	char nome_arquivo[128];
	bool direta;
	int sample_rate, channels;
	double complex *c;
	double complex *tempo;
	double complex *f;
	double complex *f2;


	//teste_complexos();

	/* Execucao do codigo */
    printf("EP2 - Numerico\n\n");
    printf("1 - Teste a\n");
    printf("2 - Teste b\n");
    printf("3 - Teste c\n");
    printf("4 - Analise de arquivo\n\n");
    printf("Digite o numero do problema a ser resolvido: ");
    scanf("%d", &tipo_problema);
    printf("\n");

    switch(tipo_problema) {

        case 1:
			teste_inicial_a();
            break;

        case 2:
        	teste_inicial_b();
            break;

        case 3:
        	teste_inicial_c();
            break;

        case 4:
		    printf("Digite o nome do arquivo a ser analizado (com a terminacao .dat): ");
		    scanf("%s", nome_arquivo);
		    printf("\n");
		    printf("Transformada direta ou inversa?\n");
		    printf("1 - Direta\n");
		    printf("2 - Inversa\n");
		    printf("Tipo de transformada: ");
		    scanf("%d", &tipo_transformada);
		    printf("\n");

		    if(tipo_transformada == 2) {
		    	direta = false;
		    }
		    else {
		    	direta = true;
		    }

		    // N = o numero total de amostras. Usado para a alocacao de vetores
			N = tamanho_arquivo(nome_arquivo, &channels);

			printf("N = %d\n", N);

			// Alocacao de vetores
			c = criarVetorComplexo(2*N);
			tempo = criarVetorComplexo(2*N);
			f = criarVetorComplexo(2*N);
			f2 = criarVetorComplexo(2*N);

			// Leitura de fato do arquivo
			sample_rate = ler_arquivo_dat(nome_arquivo, N, tempo, f, f2);

			// Se houver 2 canais, a analise eh feita em ambos
			if(channels == 1) {
				fftrec(c, f, N, direta);
			}
			else if(channels == 2) {
				fftrec(c, f, N, direta);
				fftrec(c, f2, N, direta);
			}

			printf("Escrever resultado em arquivo?\n");
		    printf("1 - Sim\n");
		    printf("2 - Nao\n");
		    printf("Resposta: ");
		    scanf("%d", &escrever_arquivo);

		    if(escrever_arquivo == 1) {
		    	printf("Digite o nome do arquivo a ser escrito (com a terminacao .dat): ");
		    	scanf("%s", nome_arquivo);

		    	// Escrever resultado da analise no arquivo
		    	escrever_arquivo_dat(nome_arquivo, sample_rate, channels, N, tempo, f, f2);
		    }

		    // Desalocacao de memoria
			free(c);
			free(tempo);
			free(f);
			free(f2);
            break;

        default:
            printf("Opcao invalida!\n");
            break;
    }


	return 0;
}


/* >>>>>>>>>>>>>>>>>>>>>>>> Funcoes basicas <<<<<<<<<<<<<<<<<<<<<<<< */

double* criarVetor(int N) {
    /* Cria um vetor de N doubles. */
    double *vetor;

    vetor = (double*) calloc(N, sizeof(double));

    return vetor;
}

int* criarVetorInt(int N) {
    /* Cria um vetor de N inteiros. */
    int *vetor;

    vetor = (int*) calloc(N, sizeof(int));

    return vetor;
}

double complex* criarVetorComplexo(int N) {
    /* Cria um vetor de N doubles complexos. */
    double complex *vetor;

    vetor = (double complex*) malloc(N * sizeof(double complex));

    return vetor;
}

void teste_complexos() {
	/* Como usar numeros complexos em C: https://stackoverflow.com/questions/6418807/how-to-work-with-complex-numbers-in-c# */
	double complex z1 = 1.0 + 3.0 * I;
	double complex z2 = 1.0 - 4.0 * I;

	printf("Working with complex numbers:\n\v");

    printf("Starting values: Z1 = %.2f + %.2fi\tZ2 = %.2f %+.2fi\n", creal(z1), cimag(z1), creal(z2), cimag(z2));

    double complex sum = z1 + z2;
    printf("The sum: Z1 + Z2 = %.2f %+.2fi\n", creal(sum), cimag(sum));

    double complex difference = z1 - z2;
    printf("The difference: Z1 - Z2 = %.2f %+.2fi\n", creal(difference), cimag(difference));

    double complex product = z1 * z2;
    printf("The product: Z1 x Z2 = %.2f %+.2fi\n", creal(product), cimag(product));

    double complex quotient = z1 / z2;
    printf("The quotient: Z1 / Z2 = %.2f %+.2fi\n", creal(quotient), cimag(quotient));

    double complex conjugate = conj(z1);
    printf("The conjugate of Z1 = %.2f %+.2fi\n", creal(conjugate), cimag(conjugate));
}


int tamanho_arquivo(char *nome_arquivo, int *channels) {
	int n = 0;
	int var_temp1;
	double var_temp2, var_temp3, var_temp4;  // estas variaveis nao tem funcao de fato
	char linha[512];

	FILE *arquivo = fopen(nome_arquivo, "r");

	if(arquivo == NULL) {
        printf("\nArquivo nao encontrado\n");
        exit(EXIT_FAILURE);
    }

    // Passando pelas duas primeiras linhas
    fscanf(arquivo, "%*[^0-9]%d", &var_temp1);
    fscanf(arquivo, "%*[^0-9]%d", &var_temp1);

    *channels = var_temp1; // retornando numero de canais por meio de um ponteiro

    if(*channels == 1) {
    	while(fgets(linha, sizeof(linha), arquivo) != NULL) { /* pega uma linha de até 512 caracteres. Null quando acabar as linhas */

	        sscanf(linha, "%lf %lf", &var_temp2, &var_temp3);

	        if(var_temp2 != 0) {
	        	n++;  // contando o numero de dados
	        }
        }
    }

    else if(*channels == 2) {
    	while(fgets(linha, sizeof(linha), arquivo) != NULL) { /* pega uma linha de até 512 caracteres. Null quando acabar as linhas */

	        sscanf(linha, "%lf %lf %lf", &var_temp2, &var_temp3, &var_temp4);

	        if(var_temp2 != 0) {
	        	n++;
	        }
	    }
    }

    else {
        printf("Numero de canais nao suportado para analise\n");
    }

    return n;
}


int ler_arquivo_dat(char *nome_arquivo, int N, double complex *tempo, double complex *f, double complex *f2) {
	int i = 0;
	int sample_rate, channels;
	double var_tempo, var_f, var_f2;
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
    	while(fgets(linha, sizeof(linha), arquivo) != NULL) { /* pega uma linha de até 512 caracteres. Null quando acabar as linhas */
	        sscanf(linha, "%lf %lf", &var_tempo, &var_f);

	        tempo[i] = var_tempo;  // vetor de tempo na amostragem
	        f[i] = var_f;  // amplitude do sinal

	        i++;  // posicao de alocacao no vetor acompanha a passagem de linhas
        }
    }

    else if(channels == 2) {
    	while(fgets(linha, sizeof(linha), arquivo) != NULL) { /* pega uma linha de até 512 caracteres. Null quando acabar as linhas */
	        sscanf(linha, "%lf %lf %lf", &var_tempo, &var_f, &var_f2);

	        tempo[i] = var_tempo;  // vetor de tempo na amostragem
	        f[i] = var_f;  // amplitude do sinal do canal 1
	        f2[i] = var_f2;  // amplitude do sinal do canal 2

	        i++;  // posicao de alocacao no vetor acompanha a passagem de linhas
	    }
    }

    else {
        printf("Numero de canais nao suportado para analise\n");
    }

    fclose(arquivo);

    return sample_rate;
}


void escrever_arquivo_dat(char *nome_arquivo, int sample_rate, int channels, int N, double complex *tempo, double complex *f, double complex *f2) {
	int i;

	FILE *arquivo = fopen(nome_arquivo, "r");

	// Escreve a sample rate e o numero de canais no mesmo formato dos arquivos .dat fornecidos
	fprintf(arquivo, "; Sample Rate %d\n", sample_rate);
	fprintf(arquivo, "; Channels %d\n", channels);

	// Se houver apenas 1 canal, havera 2 colunas de dados
	if(channels == 1) {
		for(i = 0; i < N; i++) {
			fprintf(arquivo, " %lf %lf\n", tempo[i], f[i]);
		}
	}

	// Se houverem 2 canais, havera 3 colunas de dados
	else if(channels == 2) {
		for(i = 0; i < N; i++) {
			fprintf(arquivo, " %lf %lf %lf\n", tempo[i], f[i], f2[i]);
		}
	}
	else {
		printf("Ops! Numero de canais nao suportado\n");
	}
	fclose(arquivo);
}


/* >>>>>>>>>>>>>>>>>>>>>>>> Funcoes de Transformada de Fourier <<<<<<<<<<<<<<<<<<<<<<<< */
void fftrec(double complex *c, double complex *f, int n, bool dir) {
	/* Funcao feita com base no pseudo-codigo do enunciado do EP2 */

	double complex eij;
	double complex *even, *odd, *fe, *fo;
	int j; /* Variavel auxiliar*/

	even = criarVetorComplexo(n);
	odd = criarVetorComplexo(n);
	fe = criarVetorComplexo(n);
	fo = criarVetorComplexo(n);

	if(n == 1) {
		c[0] = f[0] + f[1];
		c[1] = f[0] - f[1];
	}
	else {
		for(j = 0; j < n; j++) {
			fe[j] = f[2 * j];
			fo[j] = f[2 * j + 1];
		}

		fftrec(even, fe, n/2, dir);
		fftrec(odd, fo, n/2, dir);

		for(j = 0; j < n; j++) {

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


void filtro(double complex *c, int N, double freq_corte, char *tipo) {
	/* Filtro passa-altas ou passa-baixas, dependendo do input do usuario
	* Deixa passar apenas as frequencias "permitidas", zerando o valor das outras
	*/

	if(tipo == "passa baixas") {
		for(int j = 0; j < N; j++) {
			if(j > freq_corte) {
				c[j] = 0;
			}
		}
	}
	else {
		for(int j = 0; j < N; j++) {
			if(j < freq_corte) {
				c[j] = 0;
			}
		}
	}
}


/* >>>>>>>>>>>>>>>>>>>>>>>> Testes iniciais <<<<<<<<<<<<<<<<<<<<<<<< */
void teste_inicial_a() {
	/* Documentacao: http://people.sc.fsu.edu/~jburkardt/c_src/fftpack4/fftpack4.html
	* Baseado no teste 4 em: http://people.sc.fsu.edu/~jburkardt/c_src/fftpack4/fftpack4_prb.c
	*/
	int n = 4;
	int nh;
	int i;

	double *a;
	double a0;
	double *b;
	double complex *c;
	double complex *c_linha;

	double complex *amostras;
	double complex *amostras_valores;

	// double amostras[n-1];
	// double amostras_valores[n-1];

	int *ifac;
	double *wsave;

	amostras = criarVetorComplexo(n);
	amostras_valores = criarVetorComplexo(n);

	amostras[0] = 0;
	amostras[1] = M_PI / 2;
	amostras[2] = M_PI;
	amostras[3] = 3 * M_PI / 2;

	amostras_valores[0] = 5;
	amostras_valores[1] = -1;
	amostras_valores[2] = 3;
	amostras_valores[3] = 1;

	c = criarVetorComplexo(n);

	printf("fftrec:\n");

	// printf("amostra original = (");
	// for(i = 0; i < n; i++) {
	// 	printf("%.2f", amostras_valores[i]);
	// 	if(i != n - 1) {
	// 		printf(", ");
	// 	}
	// }
	// printf(")\n");

	fftrec(c, amostras_valores, n, true);  // Transformada direta pela fftrec

	// prints
	printf("c = (");
	for(i = 0; i < n; i++) {
		printf("%.2f %+.2fi", creal(c[i]) / n, cimag(c[i]) / n);
		if(i != n - 1) {
			printf(", ");
		}
	}
	printf(")\n");

	fftrec(c, amostras_valores, n, false);  // anti-transformada pela fftrec

	// prints
	printf("amostra anti-transformada = (");
	for(i = 0; i < n; i++) {
		printf("%.2f", amostras_valores[i]);
		if(i != n - 1) {
			printf(", ");
		}
	}
	printf(")\n\n");

	// FFTPACK4:

	wsave = criarVetor(3 * n + 15);
	ifac = criarVetorInt(8);

	ezffti(&n, wsave, ifac);  // inicializacao da fftpack4

	nh = n / 2;
	a = criarVetor(nh);
	b = criarVetor(nh);
	c_linha = criarVetorComplexo(n);

	ezfftf(&n, amostras_valores, &a0, a, b, wsave, ifac);  //transformada direta de fourier

	// Conversao de valores do tipo a*cos() + b*sen() para coeficientes complexos do tipo ck
	c_linha[0] = a0 + 0 * I;
	c_linha[nh] = a[nh-1];

	for(i = 1; i < nh; i++) {
		c_linha[i] = (a[i-1] - (I * b[i-1]))/2;
		c_linha[n-i] = (a[i-1] + (I * b[i-1]))/2;
	}

	printf("FFTPACK4:\n");

	// prints
	printf("amostra original = (");
	for(i = 0; i < n; i++) {
		printf("%.2f", amostras_valores[i]);
		if(i != n - 1) {
			printf(", ");
		}
	}
	printf(")\n");

	// prints
	printf("c = (");
	for(i = 0; i < n; i++) {
		printf("%.2f %+.2fi", creal(c_linha[i]), cimag(c_linha[i]));
		if(i != n - 1) {
			printf(", ");
		}
	}
	printf(")\n");

	ezfftb(&n, amostras_valores, &a0, a, b, wsave, ifac);  //anti-transformada de fourier

	// prints
	printf("amostra anti-transformada = (");
	for(i = 0; i < n; i++) {
		printf("%.2f", amostras_valores[i]);
		if(i != n - 1) {
			printf(", ");
		}
	}
	printf(")\n");

	// desalocacao de vetores
	free(a);
	free(b);
	free(c);
	free(c_linha);
	free(amostras);
	free(amostras_valores);
	free(ifac);
	free(wsave);
}


void teste_inicial_b() {
	/* Documentacao: http://people.sc.fsu.edu/~jburkardt/c_src/fftpack4/fftpack4.html
	* Baseado no teste 4 em: http://people.sc.fsu.edu/~jburkardt/c_src/fftpack4/fftpack4_prb.c
	*/
	int n = 8;
	int nh;
	int i;

	double *a;
	double a0;       
	double *b;
	double complex *c;
	double complex *c_linha;

	double complex *amostras;
	double complex *amostras_valores;

	// double amostras[n-1];
	// double amostras_valores[n-1];

	int *ifac;
	double *wsave;

	amostras = criarVetorComplexo(n);
	amostras_valores = criarVetorComplexo(n);

	// amostras[0] = 0;
	// amostras[1] = M_PI / 2;
	// amostras[2] = M_PI;
	// amostras[3] = 3 * M_PI / 2;

	amostras_valores[0] = 6;
	amostras_valores[1] = 2;
	amostras_valores[2] = 5;
	amostras_valores[3] = 2;
	amostras_valores[4] = 11;
	amostras_valores[5] = 2;
	amostras_valores[6] = 8;
	amostras_valores[7] = 8;

	c = criarVetorComplexo(n);

	printf("fftrec:\n");

	// printf("amostra original = (");
	// for(i = 0; i < n; i++) {
	// 	printf("%.2f", amostras_valores[i]);
	// 	if(i != n - 1) {
	// 		printf(", ");
	// 	}
	// }
	// printf(")\n");

	fftrec(c, amostras_valores, n, true);  // Transformada direta pela fftrec

	printf("c = (");
	for(i = 0; i < n; i++) {
		printf("%.2f %+.2fi", creal(c[i]) / n, cimag(c[i]) / n);
		if(i != n - 1) {
			printf(", ");
		}
	}
	printf(")\n");

	fftrec(c, amostras_valores, n, false);  // Anti-transformada pela fftrec

	printf("amostra anti-transformada = (");
	for(i = 0; i < n; i++) {
		printf("%.2f", amostras_valores[i]);
		if(i != n - 1) {
			printf(", ");
		}
	}
	printf(")\n\n");

	// FFTPACK4:

	wsave = criarVetor(3 * n + 15);
	ifac = criarVetorInt(8);

	ezffti(&n, wsave, ifac);  // inicializacao da fftpack4

	nh = n / 2;
	a = criarVetor(nh);
	b = criarVetor(nh);
	c_linha = criarVetorComplexo(n);

	ezfftf(&n, amostras_valores, &a0, a, b, wsave, ifac);  // transformada direta de fourier

	// Conversao de valores do tipo a*cos() + b*sen() para coeficientes complexos do tipo ck
	c_linha[0] = a0 + 0 * I;
	c_linha[nh] = a[nh-1];

	for(i = 1; i < nh; i++) {
		c_linha[i] = (a[i-1] - (I * b[i-1]))/2;
		c_linha[n-i] = (a[i-1] + (I * b[i-1]))/2;
	}

	printf("FFTPACK4:\n");

	printf("amostra original = (");
	for(i = 0; i < n; i++) {
		printf("%.2f", amostras_valores[i]);
		if(i != n - 1) {
			printf(", ");
		}
	}
	printf(")\n");

	printf("c = (");
	for(i = 0; i < n; i++) {
		printf("%.2f %+.2fi", creal(c_linha[i]), cimag(c_linha[i]));
		if(i != n - 1) {
			printf(", ");
		}
	}
	printf(")\n");

	ezfftb(&n, amostras_valores, &a0, a, b, wsave, ifac);  // anti-transformada de fourier

	printf("amostra anti-transformada = (");
	for(i = 0; i < n; i++) {
		printf("%.2f", amostras_valores[i]);
		if(i != n - 1) {
			printf(", ");
		}
	}
	printf(")\n");

	free(a);
	free(b);
	free(c);
	free(c_linha);
	free(amostras);
	free(amostras_valores);
	free(ifac);
	free(wsave);
}


void teste_inicial_c() {
	/* Documentacao: http://people.sc.fsu.edu/~jburkardt/c_src/fftpack4/fftpack4.html
	* Baseado no teste 4 em: http://people.sc.fsu.edu/~jburkardt/c_src/fftpack4/fftpack4_prb.c
	*/
	int n = 1024;
	int nh;
	int i;
	int x;
	double funcao;

	double *a;
	double a0;
	double *b;
	double complex *c;
	double complex *c_linha;

	double complex *amostras;
	double complex *amostras_valores;

	// double amostras[n-1];
	// double amostras_valores[n-1];

	int *ifac;
	double *wsave;

	amostras_valores = criarVetorComplexo(n);
	
	for(i = 0; i < n; i++) {
		x = (i * M_PI) / 512;

		funcao = 0;
		funcao = 10 * sin(x);
		funcao += 7 * cos(30 * x);
		funcao += 11 * sin(352 * x);
		funcao -= 8 * cos(711 * x);
		amostras_valores[i] = funcao;
	}

	c = criarVetorComplexo(n);

	printf("******************* fftrec: *******************\n");

	// printf("amostra original = (");
	// for(i = 0; i < n; i++) {
	// 	printf("%.2f", amostras_valores[i]);
	// 	if(i != n - 1) {
	// 		printf(", ");
	// 	}
	// }
	// printf(")\n");

	fftrec(c, amostras_valores, n, true);  // Transformada direta pela fftrec

	printf(">>>>>>>>>>>>>>>>>>> c:\n");
	for(i = 0; i < n; i++) {
		printf("%d: %.2f %+.2fi\n", i + 1, creal(c[i]) / n, cimag(c[i]) / n);
	}
	printf("\n");

	fftrec(c, amostras_valores, n, false);  // Anti-transformada direta pela fftrec

	printf(">>>>>>>>>>>>>>>>>>> amostra anti-transformada:\n");
	for(i = 0; i < n; i++) {
		printf("%d: %.2f\n", i + 1, amostras_valores[i]);
	}
	printf("\n\n");

	// FFTPACK4:

	wsave = criarVetor(3 * n + 15);
	ifac = criarVetorInt(8);

	ezffti(&n, wsave, ifac);  // inicializacao da fftpack4

	nh = n / 2;
	a = criarVetor(nh);
	b = criarVetor(nh);
	c_linha = criarVetorComplexo(n);

	ezfftf(&n, amostras_valores, &a0, a, b, wsave, ifac);  // transformada direta de fourier

	// Conversao de valores do tipo a*cos() + b*sen() para coeficientes complexos do tipo ck
	c_linha[0] = a0 + 0 * I;
	c_linha[nh] = a[nh-1];

	for(i = 1; i < nh; i++) {
		c_linha[i] = (a[i-1] - (I * b[i-1]))/2;
		c_linha[n-i] = (a[i-1] + (I * b[i-1]))/2;
	}

	printf("******************* FFTPACK4: *******************\n");

	// prints
	printf(">>>>>>>>>>>>>>>>>>> amostra original:\n");
	for(i = 0; i < n; i++) {
		printf("%d: %.2f\n", i + 1, amostras_valores[i]);
	}
	printf("\n");

	// prints
	printf(">>>>>>>>>>>>>>>>>>>c:\n");
	for(i = 0; i < n; i++) {
		printf("%d: %.2f %+.2fi\n", i + 1, creal(c_linha[i]), cimag(c_linha[i]));
	}
	printf("\n");

	ezfftb(&n, amostras_valores, &a0, a, b, wsave, ifac);  // anti-transformada de fourier

	// prints
	printf(">>>>>>>>>>>>>>>>>>>amostra anti-transformada\n");
	for(i = 0; i < n; i++) {
		printf("%d: %.2f\n", i + 1, amostras_valores[i]);
	}
	printf("\n");

	// desalocacao de vetores
	free(a);
	free(b);
	free(c);
	free(c_linha);
	free(amostras_valores);
	free(ifac);
	free(wsave);
}