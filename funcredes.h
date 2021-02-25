/* Funciones para pasar a la funcion (void) redes (ver redutils.h)
 * y fijar las conexiones entre neuronas.
 */

#include "rozenf.h"

// Rozenfeld sin A
// parametros: g (gamma), kmin, kmax
void f_rozin(int num, float *m, float param[]) {
	int sqrn=(int)sqrt(num);
	int x,y,d[2],xx,yy;
	int i,j,k,n;
  int pos[num][DIM];
	float g=param[0];
	int kmin=(int)param[1],
			kmax=(int)param[2];
	
	for(i=0;i<num;i++) {
    pos[i][0]=i%sqrn;
    pos[i][1]=(i-pos[i][0])/sqrn;
  }

	// (i,j)>0 si va de j a i
	// la neurona n corresponde a x=n%sqrn,y=(n-x)/sqrn
	// (x,y)->n=y*sqrn+x

	/* como aca el orden no importa ya que estamos haciendo
		 las conexiones en una sola direccion (downstream)
		 seleccionamos las posiciones secuencialmente; para cada
		 una sorteamos k y conectamos con los primeros k vecinos
		 mas cercanos, es decir...
		 llegamos a la neurona i, con coordenadas xi,yi
		 conectamos con los vecinos cercanos... cada uno de estos tiene
		 coordenadas x'i,y'i que determinan una neurona j, asi que ponemos
		 (j,i)=1 y cuando acabamos normalizamos
		 */

	if(fmod(sqrt((float)num),1.)!=0) {
		fprintf(stderr,"El numero no es cuadrado!\n");
		exit(EXIT_FAILURE);
	}
	
	srand(time(NULL));

  iniciar_pares();
	for(i=0;i<num;i++) {
		ip=0;
		pp=p0;
		iniciar_numeros();

	  // la neurona n corresponde a x=n%sqrn,y=(n-x)/sqrn
		x=i%sqrn;
		y=(i-x)/sqrn;

		k=sortear_k(g,kmin,kmax);
		n=0;
		while(n<k) {
			mostrar_dupla(d);
			n++;
			pp->rest--;

			//coordenadas de la otra neurona, CPC
			xx=x+d[0]>=0?(x+d[0])%sqrn:(x+d[0]+sqrn)%sqrn;
			yy=y+d[1]>=0?(y+d[1])%sqrn:(y+d[1]+sqrn)%sqrn;

			j=yy*sqrn+xx;
			
			M(j,i)=1;

			siguiente();
		}
	}

	normalizar(num,m);
	free(n0);
	free(p0);
	p0=NULL;
  
 	npar=1; //numero de pares
	pa=0; //par actual (indice de ultimo elemento en array)
	ua=0; //ultimo a
	ip=0; //indice de pp

/*
	for(i=0;i<num;i++) {
		for(j=0;j<num;j++) {
			M(i,j)*=f(mindist2(sqrn,pos[i],pos[j]),0.25*(float)num);
		}
	}
	*/
}

// Rozenfeld sin A
// parametros: g (gamma), kmin, kmax
void f_rozout(int num, float *m, float param[]) {
	int sqrn=(int)sqrt(num);
	int x,y,d[2],xx,yy;
	int i,j,k,n;
  int pos[num][DIM];
	float g=param[0];
	int kmin=(int)param[1],
			kmax=(int)param[2];
	
	if(fmod(sqrt((float)num),1.)!=0) {
		fprintf(stderr,"El numero no es cuadrado!\n");
		exit(EXIT_FAILURE);
	}
	for(i=0;i<num;i++) {
    pos[i][0]=i%sqrn;
    pos[i][1]=(i-pos[i][0])/sqrn;
  }

	// (i,j)>0 si va de j a i
	// la neurona n corresponde a x=n%sqrn,y=(n-x)/sqrn
	// (x,y)->n=y*sqrn+x

	/* como aca el orden no importa ya que estamos haciendo
		 las conexiones en una sola direccion (downstream)
		 seleccionamos las posiciones secuencialmente; para cada
		 una sorteamos k y conectamos con los primeros k vecinos
		 mas cercanos, es decir...
		 llegamos a la neurona i, con coordenadas xi,yi
		 conectamos con los vecinos cercanos... cada uno de estos tiene
		 coordenadas x'i,y'i que determinan una neurona j, asi que ponemos
		 (j,i)=1 y cuando acabamos normalizamos
		 */

	srand(time(NULL));

  iniciar_pares();
	for(i=0;i<num;i++) {
		ip=0;
		pp=p0;
		iniciar_numeros();

	  // la neurona n corresponde a x=n%sqrn,y=(n-x)/sqrn
		x=i%sqrn;
		y=(i-x)/sqrn;

		k=sortear_k(g,kmin,kmax);
		n=0;
		while(n<k) {
			mostrar_dupla(d);
			n++;
			pp->rest--;

			//coordenadas de la otra neurona, CPC
			xx=x+d[0]>=0?(x+d[0])%sqrn:(x+d[0]+sqrn)%sqrn;
			yy=y+d[1]>=0?(y+d[1])%sqrn:(y+d[1]+sqrn)%sqrn;

			j=yy*sqrn+xx;
			
			M(i,j)=1;

			siguiente();
		}
	}

	normalizar(num,m);
	free(n0);
	free(p0);
	p0=NULL;
  
 	npar=1; //numero de pares
	pa=0; //par actual (indice de ultimo elemento en array)
	ua=0; //ultimo a
	ip=0; //indice de pp

}

void f_roz(int num, float *m, float param[]) {
	int sqrn=(int)sqrt(num);
	int x,y,d[2],xx,yy;
	int i,j,k,n;
  int pos[num][DIM];
	float g=param[0];
	int kmin=(int)param[1],
			kmax=(int)param[2];
	float A=param[3];
	int r,temp;

	int ks[num], // k sorteado
		 	ka[num]; // k actual (todos empiezan con 0)
	int is[num]; // indices desordenados

	if(fmod(sqrt((float)num),1.)!=0) {
		fprintf(stderr,"El numero no es cuadrado!\n");
		exit(EXIT_FAILURE);
	}

	for(i=0;i<num;i++) {
		ka[i]=0;
		// indices ordenados
		is[i]=i;
	}

	//desordenar indices...
	srand(time(NULL));
	for(i=0;i<(num-1);i++) {
		r = i + (rand()%(num-i));
		temp = is[i];
		is[i]=is[r]; 
		is[r]=temp;
	}

	//posiciones...
	for(i=0;i<num;i++) {
    pos[i][0]=i%sqrn;
    pos[i][1]=(i-pos[i][0])/sqrn;
  }

	//sortear k para cada sitio
	for(i=0;i<num;i++) {
		ks[i]=sortear_k(g,kmin,kmax);
	}

	iniciar_pares();
	for(i=0;i<num;i++) {
		// (el indice es en realidad is[i])
		ip=0;	pp=p0;
		iniciar_numeros();

	  // la neurona n corresponde a x=n%sqrn,y=(n-x)/sqrn
		x=is[i]%sqrn;
		y=(is[i]-x)/sqrn;

		//n=0;
		while(/*n<k*/ ka[is[i]]<ks[is[i]]) {
			mostrar_dupla(d);

			// si se pasa el radio permitido NO
			if( (d[0]*d[0]+d[1]*d[1]) > A*A*((float)ks[i]) ) break;

			//coordenadas de la otra neurona, CPC
			xx=x+d[0]>=0?(x+d[0])%sqrn:(x+d[0]+sqrn*sqrn)%sqrn;
			yy=y+d[1]>=0?(y+d[1])%sqrn:(y+d[1]+sqrn*sqrn)%sqrn;

			j=yy*sqrn+xx;
			//printf("x=%d, y=%d, j=%d\n",xx,yy,j);

			// si j ya tiene si cuota llena
			if(ka[j]==ks[j]) {
				pp->rest--;
				siguiente();
				continue; // pasar al siguiente nodo... 
			}

			ka[ is[i] ]++;
			if(j!=is[i]) ka[j]++; //para evitar que para d=(0,0) => ka[i]+=2;
			
			M(j,(is[i]))=1;
			M((is[i]),j)=1;

			pp->rest--;
			siguiente();
		}
	}

	free(n0);
	free(p0);
	p0=NULL;
  
 	npar=1; //numero de pares
	pa=0; //par actual (indice de ultimo elemento en array)
	ua=0; //ultimo a
	ip=0; //indice de pp
	
	normalizar(num,m);
}

/* randomizacion de red de rozenfeld
	(se asume que k_min>1)
 */
void f_rozrand(int num, float *m, float param[]) {
	int conteo, veces, j;
	int i1, j1, *li1=NULL, *nli1=NULL, tamli1;
	int i2, j2, *li2=NULL, *nli2=NULL, tamli2;

	f_roz(num,m,param);
	denormalizar(num,m);

	srand(time(NULL));

	veces=0;
	while(veces<4*num) {
		// tomamos cualquier nodo como el primero
		i1=rand()%num;
		// creamos listas de conectados a i1 (li1), no conectados (nli1) y cantidad de conectados (tamli1).
		// no contamos las autointeracciones; entonces si li1 tiene tamli1 cantidad de elementos, nli1 tendra (num-1-tamli1)
		crear_listas(i1,&tamli1,&li1,&nli1,num,m);

		j1=li1[rand()%tamli1];
		
		conteo=0;
		do {
			// queremos i2 no conectado ni con i1 ni con j1
			i2=nli1[rand()%(num-1-tamli1)];
			conteo++;
		} while(M(j1,i2)!=0||conteo<num);
		if(conteo==num) continue;

		// hasta el momento se tiene i1 conectado a j1 y a la vez estos no estan conectados a i2... necesitamos ahora j2 conectado a j1 pero no conectado a i1 ni j1
		crear_listas(i2,&tamli2,&li2,&nli2,num,m);

		conteo=0;
		do {
			// queremos j2 no conectado ni con i1 ni con j1 pero si con i2
			j2=li2[rand()%(tamli2)];
			conteo++;
		} while((M(i1,j2)!=0&&M(j1,j2)!=0)||conteo<num);
		if(conteo==num) continue;

		/*interseccion(li2,nli1,tamli2,num-1-tamli1,&inter,&taminter);
		if(taminter==0) continue;

		j2=inter[rand()%taminter];
		*/

		// cruzamos cables...
		M(i1,j1)=0; M(j1,i1)=0;
		M(i2,j2)=0; M(j2,i2)=0;
		
		M(i1,j2)=1; M(j2,i1)=1;
		M(i2,j1)=1; M(j1,i2)=1;

		veces++;
	}

	normalizar(num,m);
}

// primeros vecinos c/ nro de vecinos como parametro
void f_pvnvec(int num, float *m, float param[]) {
	int sqrn=(int)sqrt(num);
	int x,y,d[2],xx,yy;
	int i,j,k,n;
  int pos[num][DIM];
	float g=2;
	int kmin=(int)param[0],
			kmax=(int)param[0];
	
	for(i=0;i<num;i++) {
    pos[i][0]=i%sqrn;
    pos[i][1]=(i-pos[i][0])/sqrn;
  }

	// (i,j)>0 si va de j a i
	// la neurona n corresponde a x=n%sqrn,y=(n-x)/sqrn
	// (x,y)->n=y*sqrn+x

	/* como aca el orden no importa ya que estamos haciendo
		 las conexiones en una sola direccion (downstream)
		 seleccionamos las posiciones secuencialmente; para cada
		 una sorteamos k y conectamos con los primeros k vecinos
		 mas cercanos, es decir...
		 llegamos a la neurona i, con coordenadas xi,yi
		 conectamos con los vecinos cercanos... cada uno de estos tiene
		 coordenadas x'i,y'i que determinan una neurona j, asi que ponemos
		 (j,i)=1 y cuando acabamos normalizamos
		 */

	if(fmod(sqrt((float)num),1.)!=0) {
		fprintf(stderr,"El numero no es cuadrado!\n");
		exit(EXIT_FAILURE);
	}
	
	srand(time(NULL));

  iniciar_pares();
	for(i=0;i<num;i++) {
		ip=0;
		pp=p0;
		iniciar_numeros();

	  // la neurona n corresponde a x=n%sqrn,y=(n-x)/sqrn
		x=i%sqrn;
		y=(i-x)/sqrn;

		k=sortear_k(g,kmin,kmax);
		n=0;
		while(n<k) {
			mostrar_dupla(d);
			n++;
			pp->rest--;

			//coordenadas de la otra neurona, CPC
			xx=x+d[0]>=0?(x+d[0])%sqrn:(x+d[0]+sqrn)%sqrn;
			yy=y+d[1]>=0?(y+d[1])%sqrn:(y+d[1]+sqrn)%sqrn;

			j=yy*sqrn+xx;
			
			M(j,i)=1;

			siguiente();
		}
	}

	normalizar(num,m);
	free(n0);
	free(p0);
	p0=NULL;
  
 	npar=1; //numero de pares
	pa=0; //par actual (indice de ultimo elemento en array)
	ua=0; //ultimo a
	ip=0; //indice de pp

/*
	for(i=0;i<num;i++) {
		for(j=0;j<num;j++) {
			M(i,j)*=f(mindist2(sqrn,pos[i],pos[j]),0.25*(float)num);
		}
	}
	*/
}


// Rozenfeld sin A no autocrina
// parametros: g (gamma), kmin, kmax
void f_rozna(int num, float *m, float param[]) {
	int sqrn=(int)sqrt(num);
	int x,y,d[2],xx,yy;
	int i,j,k,n;
  int pos[num][DIM];
	float g=param[0];
	int kmin=(int)param[1],
			kmax=(int)param[2];
	
	for(i=0;i<num;i++) {
    pos[i][0]=i%sqrn;
    pos[i][1]=(i-pos[i][0])/sqrn;
  }

	// (i,j)>0 si va de j a i
	// la neurona n corresponde a x=n%sqrn,y=(n-x)/sqrn
	// (x,y)->n=y*sqrn+x

	/* como aca el orden no importa ya que estamos haciendo
		 las conexiones en una sola direccion (downstream)
		 seleccionamos las posiciones secuencialmente; para cada
		 una sorteamos k y conectamos con los primeros k vecinos
		 mas cercanos, es decir...
		 llegamos a la neurona i, con coordenadas xi,yi
		 conectamos con los vecinos cercanos... cada uno de estos tiene
		 coordenadas x'i,y'i que determinan una neurona j, asi que ponemos
		 (j,i)=1 y cuando acabamos normalizamos
		 */

	if(fmod(sqrt((float)num),1.)!=0) {
		fprintf(stderr,"El numero no es cuadrado!\n");
		exit(EXIT_FAILURE);
	}
	
	srand(time(NULL));

  iniciar_pares();
	for(i=0;i<num;i++) {
		ip=0;
		pp=p0;
		iniciar_numeros();

	  // la neurona n corresponde a x=n%sqrn,y=(n-x)/sqrn
		x=i%sqrn;
		y=(i-x)/sqrn;

		k=sortear_k(g,kmin,kmax);
		n=0;
		while(n<k) {
			mostrar_dupla(d);
			n++;
			pp->rest--;

			//coordenadas de la otra neurona, CPC
			xx=x+d[0]>=0?(x+d[0])%sqrn:(x+d[0]+sqrn)%sqrn;
			yy=y+d[1]>=0?(y+d[1])%sqrn:(y+d[1]+sqrn)%sqrn;

			j=yy*sqrn+xx;
			
			M(j,i)=1;

			siguiente();
		}
	}

	//eliminar 1s de la diagonal!!! (antes de normalizar)

	for(i=0;i<num;i++) M(i,i)=0;

	normalizar(num,m);
	free(n0);
	free(p0);
	p0=NULL;
  
 	npar=1; //numero de pares
	pa=0; //par actual (indice de ultimo elemento en array)
	ua=0; //ultimo a
	ip=0; //indice de pp

/*
	for(i=0;i<num;i++) {
		for(j=0;j<num;j++) {
			M(i,j)*=f(mindist2(sqrn,pos[i],pos[j]),0.25*(float)num);
		}
	}
	*/
}

// primeros vecinos cuadrado con cpc con dmax como parametro
void f_pvcp(int num, float *m, float param[]){
  int i,j;
  float fac=1/((float)num);
  float d2, v;
  int pos[num][DIM];
	float dmax2;
	int sqrn=(int)sqrt(num);

	if(fmod(sqrt((float)num),1.)!=0) {
		fprintf(stderr,"El numero no es cuadrado!\n");
		exit(EXIT_FAILURE);
	}

	dmax2=param[0]*param[0];
  
  for(i=0;i<num;i++) {
    pos[i][0]=i%sqrn;
    pos[i][1]=(i-pos[i][0])/sqrn;
  }

  for(i=0;i<num;i++) {
		for(j=0;j<num;j++) {
			M(i,j)=1;//fac;
		}
  }
	normalizar(num,m);

  for(i=0;i<num-1;i++) {
    for(j=i+1;j<num;j++) {
      d2=mindist2(sqrn,pos[i],pos[j]);
      if(d2>dmax2) {
        M(i,j)=0; 
        M(j,i)=0; 
      }/* else {
        //v=fac;//f(d2,0.25*(float)num)*fac; 
        M(i,j)=1; 
        M(j,i)=1;
      }*/
    }
  }
}

// primeros vecinos cuadrado sin cpc
void f_pvc(int num, float *m, float param[]){
  int i,j;
  float fac=1/((float)num);
  float d2, v;
  int pos[num][DIM];
	float dmax2;
	int sqrn=(int)sqrt(num);

	if(fmod(sqrt((float)num),1.)!=0) {
		fprintf(stderr,"El numero no es cuadrado!\n");
		exit(EXIT_FAILURE);
	}

	dmax2=param[0]*param[0];
  
  for(i=0;i<num;i++) {
    pos[i][0]=i%sqrn;
    pos[i][1]=(i-pos[i][0])/sqrn;
  }

  for(i=0;i<num;i++) {
    M(i,i)=1;//fac;
  }

  for(i=0;i<num-1;i++) {
    for(j=i+1;j<num;j++) {
      d2=dist2(pos[i],pos[j]);
      if(d2>dmax2) {
        M(i,j)=0; 
        M(j,i)=0; 
      } else {
        //v=fac;//f(d2,0.25*(float)num)*fac; 
        M(i,j)=1; 
        M(j,i)=1;
      }
    }
  }
	normalizar(num,m);
}


// primeros vecinos cuadrado con cpc + int. endocrina relativa
void f_pvcpr(int num, float *m, float param[]){
  int i,j;
  float fac=1/((float)num);
  float d2, v;
  int pos[num][DIM];
	float dmax2;
	float intrel;
	int sqrn=(int)sqrt(num);

	if(fmod(sqrt((float)num),1.)!=0) {
		fprintf(stderr,"El numero no es cuadrado!\n");
		exit(EXIT_FAILURE);
	}

	dmax2=param[0]*param[0];
	intrel=param[1];
  
  for(i=0;i<num;i++) {
    pos[i][0]=i%sqrn;
    pos[i][1]=(i-pos[i][0])/sqrn;
  }

  for(i=0;i<num;i++) {
    M(i,i)=intrel;//fac;
  }

  for(i=0;i<num-1;i++) {
    for(j=i+1;j<num;j++) {
      d2=mindist2(sqrn,pos[i],pos[j]);
      if(d2>dmax2) {
        M(i,j)=0; 
        M(j,i)=0; 
      } else {
        //v=fac;//f(d2,0.25*(float)num)*fac; 
        M(i,j)=1; 
        M(j,i)=1;
      }
    }
  }
	normalizar(num,m);
}


// Aleatoria
void f_a(int num, float *m, float param[]) {
  int i,j;
  long s;
	float p;

	p=param[0];

  srand(time(NULL));
  s=-(long)rand();
  for(i=0;i<num;i++) {
    for(j=0;j<num;j++) {
      M(i,j)=((ran1(&s)<=p)?1:0);
    }
  }

	normalizar(num,m);
}

// Aleatoria autocrina
void f_aauto(int num, float *m, float param[]) {
  int i,j;
  long s;
	float p;

	p=(param[0]-1/(float)num)/(1-1/(float)num);;

  srand(time(NULL));
  s=-(long)rand();
  for(i=0;i<num;i++) {
    for(j=0;j<num;j++) {
			if(i==j) M(i,i)=1;
			else {
				M(i,j)=((ran1(&s)<=p)?1:0);
			}
    }
  }

	normalizar(num,m);
}

// Aleatoria autocrina simetrica
void f_aautosim(int num, float *m, float param[]) {
  int i,j;
  long s;
	float p;
	int r;

	p=(param[0]-1/(float)num)/(1-1/(float)num);;

  srand(time(NULL));
  s=-(long)rand();

  for(i=0;i<num;i++) M(i,i)=1;

  for(i=0;i<num;i++) {
    for(j=i+1;j<num;j++) {
	r=((ran1(&s)<=p)?1:0);
	M(i,j)=r;
	M(j,i)=r;
    }
  }

	normalizar(num,m);
}


// primeros vecinos
void f_pv(int num, float *m, float param[]){
  int i,j;
  float fac=1/((float)num);
  float d2, v;
  int pos[num][DIM];
	int lx, ly;
	float dmax2;

	lx=(int)param[0];
	ly=(int)param[1];
	dmax2=param[2]*param[2];
  
  if(num>lx*ly) {
    fprintf(stderr,"n>lx*ly\n");
    exit(EXIT_FAILURE);
  }
  for(i=0;i<num;i++) {
    pos[i][0]=i%lx;
    pos[i][1]=(i-pos[i][0])/lx;
  }

  for(i=0;i<num;i++) {
    M(i,i)=fac;
  }

  for(i=0;i<num-1;i++) {
    for(j=i+1;j<num;j++) {
      d2=dist2(pos[i],pos[j]);
      if(d2>dmax2) {
        M(i,j)=0; 
        M(j,i)=0; 
      } else {
        v=f(d2,dmax2)*fac; 
        M(i,j)=v; 
        M(j,i)=v;
      }
    }
  }
}

// Aleatoria cuadrada con cpc
void f_acp(int num, float *m, float param[]) {
  int i,j;
	float d2;
  long s;
	float p;
  int pos[num][DIM];
	int sqrn=(int)sqrt(num);

	p=param[0];

	if(fmod(sqrt((float)num),1.)!=0) {
		fprintf(stderr,"El numero no es cuadrado!\n");
		exit(EXIT_FAILURE);
	}
	
	for(i=0;i<num;i++) {
    pos[i][0]=i%sqrn;
    pos[i][1]=(i-pos[i][0])/sqrn;
  }
	
  srand(time(NULL));
  s=-(long)rand();
  for(i=0;i<num;i++) {
    for(j=0;j<num;j++) {
      M(i,j)=((ran1(&s)<=p)?1:0);
    }
  }

	normalizar(num,m);
	
	for(i=0;i<num;i++) {
		for(j=0;j<num;j++) {
			d2=mindist2(sqrn,pos[i],pos[j]);
			M(i,j)*=f(d2,0.25*(float)num);
		}
	}
}

// primeros vecinos + aleatoria
void f_pva(int num, float *m, float param[]) {
	int i,j;
  float fac=1/((float)num);
  float d2, v;
	float z;
	long s;
  int pos[num][DIM];
  int lx, ly;
	float dmax2;
	float p;

	lx=(int)param[0];
	ly=(int)param[1];
	dmax2=param[2]*param[2];
  p=param[3];

  if(num>lx*ly) {
    fprintf(stderr,"n>lx*ly\n");
    exit(EXIT_FAILURE);
  }
  for(i=0;i<num;i++) {
    pos[i][0]=i%lx;
    pos[i][1]=(i-pos[i][0])/lx;
  }

  for(i=0;i<num;i++) {
    M(i,i)=fac;
  }

  for(i=0;i<num-1;i++) {
    for(j=i+1;j<num;j++) {
      d2=dist2(pos[i],pos[j]);
      if(d2>dmax2) {
        M(i,j)=0; 
        M(j,i)=0; 
      } else {
        v=1; 
        M(i,j)=v; 
        M(j,i)=v;
      }
    }
  }

	srand(time(NULL));
  s=-(long)rand();

  for(i=0;i<num;i++) {
    for(j=0;j<num;j++) {
			// no tocar lo que ya esta conectado por primeros vecinos
      if(M(i,j)!=1) M(i,j)=((ran1(&s)<=p)?1:0);
    }
  }

	// normaliza
  for(i=0;i<num;i++) {
    z=0;
    for(j=0;j<num;j++) {
      z+=M(i,j);
    }
    if(z!=0) {
      for(j=0;j<num;j++) {
        M(i,j)/=z;
      }
    }
  }
}

void f_s(int num, float *m, float param[]) {
	int i,j;

  for(i=0;i<num;i++) {
    for(j=0;j<num;j++) {
      M(i,j)=((i==j)?1:0);
    }
  }
}

// Metodo de Barabasi-Albert:
// parametros: n_0 (numero inicial de nodos), p (probabilidad de conexion aleatoria inicial)
void f_ba(int num, float *m, float param[]) {
	int i,j,l,k,suma,sumak,n0;
	float p0,r,p;
	long s;

	// numero inicial de neuronas conectadas (n0<num, p>=2/(n0+1))
	n0=(int)param[0];
	p0=param[1];

	if(n0>num) {
		fprintf(stderr, "error, n0 > numero de neuronas\n");
	/*}	else if(p0<2./(float)(n0+1)) {
		fprintf(stderr, "error, p0 = %f < 2/(n0+1) = %f\n", p0, 2./((float)(n0+1)));
	*/}

	// conectar las primeras n0 neuronas de forma simetrica (y autocrina, por ahora)
	srand(time(NULL));
  s=-(long)rand();
  for(i=0;i<n0;i++) {
    for(j=i+1;j<n0;j++) {
			r=((ran1(&s)<=p0)?1:0);

			M(i,j)=r;
			M(j,i)=r;
    }
  }
	for(i=0;i<n0;i++) M(i,i)=1; //autocrina

	/*
	//verificar que la red este conectada.... (de otra forma volver a correr(?))
	for(i=0;i<n0;i++) {
		suma=0;
		for(j=i+1;j<n0;j++) {
			suma+=M(i,j);
		}
		if(suma==0) {
			fprintf(stderr,"La red no esta conectada del todo (existen nodos aislados), intentando de nuevo...\n");
			//exit(1);
			f_ba(num, m, param);
		}
	}
	*/

	// y ahora empiezo a generar la red de barabasi albert simetrica...

	// para cada nuevo nodo
	for(i=n0;i<num;i++) {
		// nos fijamos en cada nodo previo y vemos cuantas conexiones tiene,
		// relativo al numero total de conexiones existentes...
		
		// asi que primero sumamos todos los k, antes de agregar el nuevo nodo
		// (obviamente esto se podria hacer mas eficiente, pero ahora no vale la pena
		sumak=0;
		//sumamos todos los 1's de la seccion triangular inferior-izquierda de M,
		//en la submatriz
		for(j=0;j<i;j++) {
			for(l=0;l<j;l++) {
				sumak+=(int)M(j,l);
			}
		}
		// ahora tenemos la cantidad de conexiones totales en sumak

		// ahora por cada nodo calculamos k_i, es decir, la cantidad de conexiones
		// que salen/entran (es simetrica la red)

		//(para cada nodo preexistente...)
		for(j=0;j<i;j++) {
			// queremos saber cuantas conexiones tiene el nodo j,  es decir, suma de M(j,x), x!=j
			k=0;
			for(l=j+1;l<i;l++) {
				k+=(int)M(j,l);
				//k+=(int)M(l,j);
			}

			p=((float)k)/((float)sumak);
			//printf("p(conexion(%d,%d)=%f\n",i,j,p);
	//		getchar();
			r=((ran1(&s)<=p) ? 1:0);

			M(i,j)=r;
			M(j,i)=r;
		}
		M(i,i)=1;
	}

	normalizar(num,m);
}
