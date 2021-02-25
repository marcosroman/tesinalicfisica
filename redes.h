/* Este header contiene, ademas de un par de macros, la estructura 'topologia'
 * (de red) y la funcion (void) 'redes' que se utiliza para definir la forma
 * de interacccion entre las neuronas
 */

#include "funcredes.h"

// constantes para redes
#define MAXLABEL 20
#define MAXPARAM 5
#define MAXLEND 30
#define MAXLENP 35

char **etiqpar;

// estructura topologia [de red]
typedef struct {
	int numparam;
	char etiqueta[MAXLABEL];
	char desc_etiqueta[MAXLEND];
	char desc_param[MAXPARAM][MAXLENP];
	void (*ftop)(int num, float *m, float param[]);
} topologia;

// tipos de redes
#define TIPOS 10

/* esta funcion tiene multiples usos dependiendo de la tarea que se le asigne:
  [definida con el parametro 'tarea']
	 0=listar redes (LISTAR)
	 1=fijar codigos y verificar correcto numero de parametros (VERIFICAR)
	 2=fijar funcion de red correspondiente (FUNCION)
*/
void redes(
		int tarea, int *indred, char *nred, int *npars,
		void (**funcion)(int num, float *m, float param[])
		) {
	int i,j,k=0;
	topologia topo[TIPOS];

	//0
	strcpy(topo[k].etiqueta,"al");
  strcpy(topo[k].desc_etiqueta,"Aleatoria");
	topo[k].numparam=1;
	strcpy(topo[k].desc_param[0],"conectividad");
	topo[k].ftop=&f_a;
	k++;

	//1
	strcpy(topo[k].etiqueta,"alauto");
	strcpy(topo[k].desc_etiqueta,"Aleatoria autocrina");
	topo[k].numparam=1;
	strcpy(topo[k].desc_param[0],"conectividad aprox.");
	topo[k].ftop=&f_aauto;
	k++;

	//2
	strcpy(topo[k].etiqueta,"alautosim");
	strcpy(topo[k].desc_etiqueta,"Aleatoria autocrina simetrica");
	topo[k].numparam=1;
	strcpy(topo[k].desc_param[0],"p");
	topo[k].ftop=&f_aautosim;
	k++;

	//3
	strcpy(topo[k].etiqueta,"pv");
	strcpy(topo[k].desc_etiqueta,"Primeros vecinos");
	topo[k].numparam=1;
	strcpy(topo[k].desc_param[0],"num vec");
	topo[k].ftop=&f_pvnvec;
	k++;

	//4
	strcpy(topo[k].etiqueta,"rozin");
	strcpy(topo[k].desc_etiqueta,"Rozenfeld (in)");
	topo[k].numparam=3;
	strcpy(topo[k].desc_param[0],"gamma");
	strcpy(topo[k].desc_param[1],"k_min");
	strcpy(topo[k].desc_param[2],"k_max");
	topo[k].ftop=&f_rozin;
	k++;

  //5
	strcpy(topo[k].etiqueta,"rozout");
	strcpy(topo[k].desc_etiqueta,"Rozenfeld (out)");
	topo[k].numparam=3;
	strcpy(topo[k].desc_param[0],"gamma");
	strcpy(topo[k].desc_param[1],"k_min");
	strcpy(topo[k].desc_param[2],"k_max");
	topo[k].ftop=&f_rozout;
	k++;
	
	//6
	strcpy(topo[k].etiqueta,"roz");
	strcpy(topo[k].desc_etiqueta,"Rozenfeld");
	topo[k].numparam=4;
	strcpy(topo[k].desc_param[0],"gamma");
	strcpy(topo[k].desc_param[1],"k_min");
	strcpy(topo[k].desc_param[2],"k_max");
	strcpy(topo[k].desc_param[2],"A");
	topo[k].ftop=&f_roz;
	k++;

	//7
	strcpy(topo[k].etiqueta,"rozrand");
	strcpy(topo[k].desc_etiqueta,"Rozenfeld rand");
	topo[k].numparam=4;
	strcpy(topo[k].desc_param[0],"gamma");
	strcpy(topo[k].desc_param[1],"k_min");
	strcpy(topo[k].desc_param[2],"k_max");
	strcpy(topo[k].desc_param[2],"A");
	topo[k].ftop=&f_rozrand;
	k++;

	//8
	strcpy(topo[k].etiqueta,"ba");
	strcpy(topo[k].desc_etiqueta,"Barabasi-Albert");
	topo[k].numparam=2;
	strcpy(topo[k].desc_param[0],"n0");
	strcpy(topo[k].desc_param[1],"p0");
	topo[k].ftop=&f_ba;
	k++;

	//9	
	strcpy(topo[k].etiqueta,"pvcp");
	strcpy(topo[k].desc_etiqueta,"PV cuad. c/ CPC");
	topo[k].numparam=1;
	strcpy(topo[k].desc_param[0],"d_max");
	topo[k].ftop=&f_pvcp;
	k++;
	
	switch (tarea) {
		case LISTAR:
			for(i=0;i<TIPOS;i++) {
				//printf("%d: %s  (%s) {",i, topo[i].etiqueta, topo[i].desc_etiqueta);
				printf("  %s  (%s) {",topo[i].etiqueta, topo[i].desc_etiqueta);
				for(j=0;j<topo[i].numparam;j++) {
					printf("%s",topo[i].desc_param[j]);
					if(j!=(topo[i].numparam-1)) {
						printf(",");
					}
				}
				printf("}\n");
			}
			//listar.... indice, etiqueta, nombre completo, parametros
			break;
		case VERIFICAR:
			if(*indred>=0) {
				if(*indred>=TIPOS) {
					printf("no existe la red indred= %d\n",*indred);
					exit(1);
				}
				sprintf(nred,"%s",topo[*indred].etiqueta);
			} else {
				j=0;
				while(j<TIPOS) {
					if(strcmp(nred,topo[j].etiqueta)==0) break;
					j++;
				}
				if(j>=TIPOS) {
					printf("no existe la red %s\n",nred);
					exit(1);
				} else {
					*indred=j;
				}
			}
			if(topo[*indred].numparam!=*npars) {
				printf("error en numero de parametros!\n");
				exit(1);
			}
			break;
		case FUNCION:
			*funcion=topo[*indred].ftop;
			break;
		case GETNUMPARS:
			j=0;
			while(j<TIPOS) {
				if(strcmp(nred,topo[j].etiqueta)==0) break;
				j++;
			}
			if(j>=TIPOS) {
				printf("no existe la red %s\n",nred);
				exit(1);
			} else {
				*indred=j;
			}
			*npars=topo[*indred].numparam;
			break;
		default:
			printf("algun error...\n");
			exit(1);
	}
}
