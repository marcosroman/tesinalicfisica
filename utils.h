/* Macros, definiciones y funciones utiles varias para correr el modelo
 */

// tiempo
#define PASO_T 0.1
#define TMAX 24*30 //sacar el 2 dsp!!!
#define TMIN 0
#define TMAV (TMAX-24*5)
#define TPOR (TMAX-24*4)

// el sistema tiene 10 variables
#define VARS 10
// estamos en un plano
#define DIM  2

// geometria euclidiana XD
#define PI   3.14159265

// para simplificar la lectura, especialmente con la matriz
#define Y(i) z[i-1]
#define V    z[7]
#define X(i) z[7+i]
#define M(x,y)  m[x*num+y]

#define A(x,y) a[x*n+y]

#define K_INICIAL 0.1
#define K_FINAL 1.05
#define K_PASO 0.1

//tareas
#define LISTAR 0
#define VERIFICAR 1
#define FUNCION 2
#define GETNUMPARS 3
#define ETIQUETAS 4

// estructura neurona
typedef struct {
  float z[VARS];
  float zout[VARS];
  float dz[VARS];
  float e;
} neurona;

// parametros para generar nros aleatorios
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

void redes(int tarea, int *indred, char *nred, int *npars,void (**funcion)(int num, float *m, float param[]));
void dzdt(float t, float q, float e, float z[], float dz[]);
void rk4(float z[], float dzdt[], int n, float t, float h, float q, float e, float zout[], void (*derivs)(float, float, float, float [], float []));
void iniciar_neuronas(neurona *n, int num);

// ran1, de Numerical Recipes in C
float ran1(long *idum) {
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

void iniciar_m_n(int num, neurona **n, float **m) {
  *n=malloc(sizeof(neurona)*(num));
  *m=calloc((num)*(num),sizeof(float));

  if(*m==NULL||*n==NULL) {
    fprintf(stderr, "No hay memoria para tantas neuronas!\n");
    free(*m);free(*n);
    exit(EXIT_FAILURE);
  }
}

// gasdev, de Numerical Recipes in C
float gasdev(long *idum) {
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if(*idum<0) iset=0;
	if(iset==0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

void imprimir_PerCry(int num, neurona *n) {
  int i,j;
	int sqrn=(int)sqrt(num);

	for(i=0;i<sqrn;i++) {
		for(j=0;j<sqrn;j++) {
			printf("%f ",n[j*sqrn+i].z[0]);
		}
		printf("\n");
	}
}

void plotear_matriz_cuadrada(int n, float *a, FILE *f) {
	int x,y;

	for(x=0;x<n;x++) {
		for(y=0;y<n;y++) {
			fprintf(f,"%d %d %f\n",x,y,A(x,y));
			fprintf(f,"%d %d %f\n",x,y+1,A(x,y));
		}
		fprintf(f,"\n");

		for(y=0;y<n;y++) {
			fprintf(f,"%d %d %f\n",x+1,y,A(x,y));
			fprintf(f,"%d %d %f\n",x+1,y+1,A(x,y));
		}
		fprintf(f,"\n");
	}
}

void imprimir_flechas(int num, float *m, FILE *f) {
	int sqrn=(int)sqrt(num);
	int i,j,k;
	int xi,xj,yi,yj;

	for(k=0;k<num*num;k++) {
		if(m[k]>0) {
			//j=k%num;
			//i=(k-j)/num;
			i=k%num;
			j=(k-i)/num;

			//yi=i%sqrn;
			//xi=(i-yi)/sqrn;
			xi=i%sqrn;
			yi=(i-xi)/sqrn;


			//yj=j%sqrn;
			//xj=(j-yj)/sqrn;
			xj=j%sqrn;
			yj=(j-xj)/sqrn;


			fprintf(f==NULL?stdout:f,"set arrow from %f,%f to %f,%f lw 0.1\n",
					(int)xj+0.5,(int)yj+0.5,(int)xi+0.5,(int)yi+0.5);
		}
	}
}

void imprimir_valores(int n, float **pv, FILE *f) {
	int i;
	float *a;

	a=malloc(sizeof(float)*n*n);
	if(a==NULL) {
		exit(1);
	}

	for(i=0;i<n*n;i++) {
		a[i] = *(pv[i]);
	}

	plotear_matriz_cuadrada(n,a,f);

	free(a);
}

void iniciar_R(int num, float **s, float **s2, float *sp, float *s2p, float *pt) {
	int i;

	*s=calloc(num,sizeof(float));
	*s2=calloc(num,sizeof(float));
	if(*s==NULL || *s2==NULL) exit(1);
	
	for(i=0;i<num;i++) {
		(*s)[i]=0;
		(*s2)[i]=0;
	}

	*pt=0;
	*s2p=0;
	*sp=0;
}

void actualizar_R(int num, neurona n[], float *s, float *s2, float *sp, float *s2p, float *pt) {
	int i;
	float sx=0;

	for(i=0;i<num;i++) {
		s[i]+=n[i].Y(1);
		s2[i]+=n[i].Y(1)*n[i].Y(1);
		sx+=n[i].Y(1);
	}
	sx/=(float)num;
	*sp+=sx;
	*s2p+=sx*sx;

	(*pt)++;
}

void terminar_R(int num, float *suma, float s[], float s2[], float *sp, float *s2p, float pt) {
	int i;

	*suma = 0;
	for(i=0;i<num;i++) {
		s[i]/=pt; s2[i]/=pt;
		*suma += (s2[i]-s[i]*s[i]);
	}

	*sp /= pt;
	*s2p /= pt;
}

float devolver_R(int num, float sp, float s2p, float suma) {
	return (((float)num)*(s2p-sp*sp))/suma;
}

void liberar_R(float *s, float *s2) {
	free(s);
	free(s2);
}

float parametro_de_orden(int num, int nregs[], int nregact[], float *reg[]) {
	int i;
	float t,pt,sumpar,sumsen,sumcos;
	float fase;
	
	// calculo de parametro de orden
	t=TPOR;
	pt=0;
	sumpar=0;
	do {
		sumsen=0;sumcos=0;
		for(i=0;i<num;i++) {
			if(nregs[i]<2) continue; // (no tener en cuenta las que no oscilan...)

			if(t<reg[i][0]) {
				//pendiente de los dos sgtes puntos
				fase=fmod(reg[i][1]+(t-reg[i][0])/(reg[i][2]-reg[i][0]),2); //pend. de 0 y 1
			}	else if(t>=reg[i][2*(nregs[i]-1)]) {
				//seguir con la ultima pendiente...
				fase=fmod(reg[i][1+2*(nregs[i]-2)]+
						(t-reg[i][2*(nregs[i]-2)])/(reg[i][2*(nregs[i]-1)]-reg[i][2*(nregs[i]-2)]),2);
			}
			else { 
				if(t>=reg[i][2*(nregact[i]+1)]) nregact[i]++;
			
				fase=fmod(reg[i][1+2*nregact[i]]+(t-reg[i][2*nregact[i]])/(reg[i][2*(nregact[i]+1)]-reg[i][2*nregact[i]]),2);
			}

			sumsen+=sin(PI*fase);
			sumcos+=cos(PI*fase);
		}

		pt++;
		sumpar+=sqrt(sumsen*sumsen+sumcos*sumcos)/((float)num);

		t+=PASO_T;
	} while(t<=TMAX);

	return sumpar/pt;
}

void iniciar_fases(int num, float **x0, float **x1, float ***reg, int **nregs, int **nregact) {
	int i;

	*x0=calloc(num,sizeof(float));
	*x1=calloc(num,sizeof(float));
	*reg=calloc(num,sizeof(float*));
	*nregs=calloc(num,sizeof(int));
	*nregact=calloc(num,sizeof(int));
	for(i=0;i<num;i++) (*reg)[i]=calloc(2,sizeof(float));
}

void liberar_fases(int num, float *x0, float *x1, float **reg, int *nregs, int *nregact) {
	int i;
	
	for(i=0;i<num;i++) {
		free(reg[i]);
	}
	free(x0);
	free(x1);
	free(reg);
	free(nregs);
	free(nregact);
}

/* lee un argumento (arg) y lo guarda (en el 3-vector float p)
   el vector p se interpreta como p_inicial (p[0]), p_final (p[1]), deltap (p[2])
   hay 2 casos|tipos:
	 1. (float0)a(float1)p(float2) -> guarda (0),(1),(2)
	 	-> lo cual se lee: de float0 a float1, con paso float2
	 2. (float0) -> guarda (float0),(float0),1
*/
void leer_parametro(char arg[], float p[]) {
	int i,l;
	l=strlen(arg);

	// vamos a ver si arg es tipo 1 o tipo 2
	// (escaneando la cadena y viendo si llego al final sin encontrar una 'a')
	i=0;
	while(arg[i]!='a' && i<l) i++;
	if(l==i) {
		// si llegamos al final...
		p[0]=p[1]=atof(arg);
		p[2]=1;
	} else {
		sscanf(arg,"%fa%fp%f",&p[0],&p[1],&p[2]);	
	}

	if(p[1]<p[0]) {
		printf("error, pmin>pmax ! \n");
		exit(1);
	}
}

// esta funcion recibe argumentos y se encarga de que todo este en orden,
// sin tocar los parametros para la red;
// fija num neuronas, k y nro de corridas para cada set de parametros
// como argumentos debemos tener:
//   0.(nombre de programa)
//   1.(nro de osciladores)
//   2.(tipo de red), definido como cadena
//   (3-(ac-3)).(parametros de red, separados por espacios)
//   ...
//   (ac-2). k|'k0'a'kf'p'dp' (sin las comillas)
//      -> acoplamiento: unico numero | rango, donde dp es el step
//   (ac-1).(nro de corridas)
// donde
//   ac,av=argc,argv
//   k=[k0,kf,deltak]
void leer_argumentos_iniciales(int ac, char *av[], int *nneuronas, float k[],
	int *ncorridas, int *nparametros) {
	int i;

	/*
	if(ac>=2) {
	 	if(0==strcmp(av[1],"?")) {
			printf("Uso: %s",av[0]);
			redes(LISTAR,NULL,NULL,0,NULL);
			exit(0);
		}
	} else {
		printf("argumentos insuficientes. (? para ver lista)\n");
		exit(1);
	}

	if(ac<5) {
		printf("argumentos incompletos\n");
		exit(1);
	} 
	*/

	if(ac<3) {
		printf("Uso:\n%s <numero de neuronas> <red> {parametros de red} <acoplamiento:K(unico valor)|K0aKfdDK(rango)> <nro_corridas>\n",av[0]);
		printf("\ndonde las opciones de redes son:\n");
		redes(LISTAR,NULL,NULL,0,NULL);
		exit(1);
	}

	*nneuronas=atoi(av[1]);
	*nparametros=ac-5;
	leer_parametro(av[ac-2],k);

	if(!(*ncorridas=atoi(av[ac-1]))) {
		printf("nro corridas?\n");
		exit(1);
	}
}

// si params=NULL, alocar memoria... si no, liberar y alocar
void leer_argumentos_redes(
		int ac, char *av[], 
		int npars, float **params, int *indred, 
		char *nred, int *indpar, float *valspar) {
	int i,j=3;
	//char *red=av[2];
	float paux[3];

	if(*params!=NULL) {
		free(*params);
	}
	*params=calloc(npars,sizeof(float));
	if(*params==NULL) {
		printf("No se pudo alocar memoria\n");
		exit(1);
	}

	*indpar=-1; //por si no hay argumentos (?)

	// av[2] nos da la red, como cadena o digito
	// fija uno de los dos, codigo string de red o codigo numerico... lo q no se fija es NULL o -1
	if((*indred)=atoi(av[2])) {
		//ESTA SECCION NO FUNCIONA TODAVIA!!!!
		printf("funciona, indred=%d\n",*indred);
		getchar();
		nred=NULL;
	} else {
		*indred=-1;
		nred=av[2];
	}

	redes(VERIFICAR,indred,nred,&npars,NULL);

	for(j=3;j<ac-2;j++) {
		leer_parametro(av[j],paux);
		(*params)[j-3]=paux[0];
		if(paux[0]!=paux[1]) {
			for(i=0;i<3;i++) valspar[i]=paux[i];
			*indpar=j-3;
		}
	}
	if(*indpar<0 && npars>0) {
		valspar[0]=valspar[1]=paux[0];
		valspar[2]=1;
	}
}

void imprimir_corrida(int num, float *m, float coef_k, FILE *f) {
	int i,k;
	float t=0,s, q, h=PASO_T;
	neurona n[num];

	iniciar_neuronas(n,num);

	do {
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
		
		fprintf(f,"%e ",t);

		s=0;
		for(i=0;i<num;i++) s+=n[i].z[0];
		s/=(float)num;

		fprintf(f,"%e ",s);


		for(i=0;i<num;i++) {
			fprintf(f,"%e ",n[i].z[0]);
		}
		fprintf(f,"\n");

		// actualiza los pasos
		for(i=0;i<num;i++) for(k=0;k<VARS;k++) n[i].z[k]=n[i].zout[k];
		t+=PASO_T;
	} while(t<=TMAX);
}

