/* kcr:
 *   k es el parametro de acoplamiento
 *   c es la conectividad de red (fraccion de conexiones existentes
 *     de las n*(n-1) posibles [considerando conexiones dirigidas])
 *   r es el parametro de orden de sincronizacion
 *     ($ |\sum_{i=1}^N \exp(i \theta_i)|$, donde $\theta_i$=fase)
 */

#define DEBUG 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "utils.h"
#include "redutils.h"
#include "redes.h"
#include "circadiano.h"
#include "orden.h"

int main(int argc, char *argv[]) {
  int num; // numero de neuronas
  int i,j,k,l;

  neurona *n; // puntero para array de neuronas
  float *m; // puntero para matriz de conectividad
	float coef_k=0.8;
	float con;
	
  float t, h=PASO_T; // tiempo
	float q;

	int sprom;

	// para el calculo de fases y el parametro de orden
	int *nregs, *nregact; //cantidad de registros de picos para cada neurona
	float *x0,*x1;
	float **reg;
	float r;

	int indexparametro=-1;
	int ncorridas;
	int nparametros;
	int indexred;
	char *etiqred=NULL;
	// k0, kf y dk
	float k_valores[3];
	float *parametros=NULL; //muy importante!
	float parametro_valores[3];
	void (*funcionred)(int num, float *m, float param[]);

	float *auxper,*auxperfase;

	char filename[200]; //supongo que es suficiente...
	char charaux[100];
	FILE *fout;

	leer_argumentos_iniciales(argc,argv,&num,k_valores,&ncorridas,&nparametros);
	#if DEBUG
		printf("argumentos iniciales ok\n");
	#endif

	leer_argumentos_redes(argc,argv,nparametros,&parametros,&indexred,etiqred,&indexparametro,parametro_valores);
	#if DEBUG
		printf("argumentos redes ok\n");
	#endif

	redes(FUNCION,&indexred,NULL,&nparametros,&funcionred);

	#if DEBUG
		printf("ok. empezamos...\n");
	#endif

	// seteamos nombre de archivo
	sprintf(filename,"data_n");
	for(i=1;i<argc;i++) {
		sprintf(charaux,"%s",argv[i]);
		if(i<argc-1) strcat(charaux,"_");
		strcat(filename,charaux);
	}
	// abrimos el archivo donde guardamos los datos
	fout=fopen(filename,"w");

	// barremos sobre valores de k (intensidad de acoplamiento)
	for(coef_k=k_valores[0];coef_k<=k_valores[1];coef_k+=k_valores[2]) {
		fprintf(fout,"#k=%g\n",coef_k);
		for(parametros[indexparametro]=parametro_valores[0];
				parametros[indexparametro]<=parametro_valores[1];
				parametros[indexparametro]+=parametro_valores[2]) {
			fflush(fout);
			for(j=0;j<ncorridas;j++) {
				#if DEBUG
					printf("corrida %d, coefk=%f\n",j,coef_k);
					printf("valores k = %f, %f, %f\n", k_valores[0],k_valores[1],k_valores[2]);
					printf("valores p = %f, %f, %f\n", parametro_valores[0],parametro_valores[1],parametro_valores[2]);
				#endif
				iniciar_m_n(num,&n,&m);
				iniciar_neuronas(n,num);
				(*funcionred)(num,m,parametros);
			 	iniciar_fases(num,&x0,&x1,&reg,&nregs,&nregact);

				auxper=malloc(sizeof(float)*num);
				auxperfase=malloc(sizeof(float)*num);

				t=0;
				do {
					//sw=(t>=TMIN);
					sprom=(t>=TMAV);

				  // mover x0,x1 (para el calculo de fases)
					if(sprom) {
						for(i=0;i<num;i++) {
							x0[i]=x1[i];
							x1[i]=n[i].z[0];
						}
					}

					for(i=0;i<num;i++) {
						// calcular q_i
						q=0;
						for(k=0;k<num;k++) {
							q+=M(i,k)*n[k].V; 
						}
						q*=coef_k;
						
						// paso rk4...
						dzdt(t,q,n[i].e,n[i].z,n[i].dz);
						rk4(n[i].z,n[i].dz,VARS,t,h,q,n[i].e,n[i].zout,dzdt);
					}
					// actualiza los pasos
					for(i=0;i<num;i++) for(k=0;k<VARS;k++) n[i].z[k]=n[i].zout[k];

					//fases....
					if(sprom)	{
							actualizar_registro(num, n, t, x0, x1, nregs, &reg,auxper,auxperfase);
					}

					t+=h;
				} while(t<=TMAX);

				r=parametro_de_orden(num,nregs,nregact,reg);
				
				con=conectividad(num,m);

				fprintf(fout,"%e\t%e\t%e\t",coef_k,con,r);
				for(l=0;l<nparametros;l++) {
					fprintf(fout,"%e\t",parametros[l]);
				}
				fprintf(fout,"\n");

				//liberar memoria!
				free(n);
				free(m);

				free(auxper);
				free(auxperfase);

				liberar_fases(num,x0,x1,reg,nregs,nregact);
				}
				fprintf(fout,"\n");
			}
		fprintf(fout,"\n");
	}
	fclose(fout);
				

  return 0;
}
