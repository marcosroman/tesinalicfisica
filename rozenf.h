/* Funciones utiles para configurar una red libre de escala embebida en 2D
 * 
 * para mas info ver:
 *       Rozenfeld, Alejandro F., et al. "Scale-free networks on lattices."
 *       Physical Review Letters 89.21 (2002): 218701.
 */

typedef struct spar {
	int r[2]; // 2-vector con [0]>=[1]>=0
	int cant; // posibles combinaciones
	int rest;
//	int orden; // indice que ordena los pares
//	struct spar *ste; // siguiente par, ordenado segun modulo r2=a2+b2
} par;

typedef struct snumero {
	int n; // guarda un numero, q indica la operacion
	struct snumero *ant,*ste; // numeros anterior y siguiente
} numero;

numero *n0,*pn;
par *p0=NULL,*pp,*pf;

int npar=1, //numero de pares
		pa=0, //par actual (indice de ultimo elemento en array)
		ua=0, //ultimo a
		ip=0; //indice de pp

void mostrar_dupla(int d[]) {
	int sgn[2];
	int regla=pn->n, swap;
	
	if(regla%2==0) {swap=1;	regla/=2;}
	else {swap=0;	regla=(regla+1)/2;}

	switch(regla) {
		case 1: sgn[0]= 1; sgn[1]= 1;	break;
		case 2: sgn[0]=-1; sgn[1]=-1;	break;
		case 3: sgn[0]= 1; sgn[1]=-1;	break;
		case 4: sgn[0]=-1; sgn[1]= 1;	break;
		default: printf("hay algun error!!!\n");
	}

	d[0]=sgn[0]*(pp->r[swap%2]);
	d[1]=sgn[1]*(pp->r[(swap+1)%2]);
}

void generar_numeros() {
	int caso,cantidad;
	int i,di,k;
	int a=pp->r[0], b=pp->r[1];
	
	if(b!=0) {
		if(b!=a) {
			caso=0;//12345678
		} else if(b==a) {
			caso=1;//1357
		}
	}	else {
		caso=2;//1234
	}	

	if(caso>0) pp->cant=4;
	else pp->cant=8;

	// cantidad de operaciones que generan duplas distintas
	cantidad=pp->cant;
	n0=realloc(n0,sizeof(numero)*cantidad);
	
	// conectamos los 'numero' circularmente
	for(i=0;i<cantidad;i++) {
		n0[i].ant=&n0[(i-1)%cantidad];
		n0[i].ste=&n0[(i+1)%cantidad];
	}
	n0[0].ant=&n0[cantidad-1];
	n0[cantidad-1].ste=&n0[0];
	pp->rest=cantidad;
	pn=n0;

	// asignamos numeros
	switch(caso%2) {
		case 0:	di=1;	break;
 		case 1:	di=2; break;
	}
	k=rand()%cantidad;
	i=1;
	while(cantidad>0) {
		pn->n=i;
		i+=di;
		cantidad--;
		pn=pn->ste;
	}

	//'aleatorizar'!
	while(k>0) {
		pn=pn->ste;
		k--;
	}
}

int comparar(const void *vp1, const void *vp2) {
	par *p1=(par*)vp1,
			*p2=(par*)vp2;
	int a1=p1->r[0],a2=p2->r[0],
			b1=p1->r[1],b2=p2->r[1];
	int r12=a1*a1+b1*b1, r22=a2*a2+b2*b2;

	return r12-r22;
}

void sgte_par() {
	int a=pp->r[0];
	int aa,bb;
	int am=(int)ceil(sqrt(2)*a);
	int i,n;
	
	if(a==0) am=1;
	if(ua<am) {
		//agregar pares que faltan
		//cambiar ua, cambiar npar, cambiar pa
		//apuntar pp a el nuevo siguiente...
		n=0;
		for(aa=ua+1;aa<=am;aa++) {
			for(bb=0;bb<=aa;bb++) {
				n++;
			}
		}
		p0=realloc(p0,(npar+n)*sizeof(par));
		i=npar;
		for(aa=ua+1;aa<=am;aa++) {
			for(bb=0;bb<=aa;bb++) {
				p0[i].r[0]=aa;
				p0[i].r[1]=bb;
				i++;
			}
		}
		npar+=n;
		ua=am;
		//ordenar!
		qsort(p0,npar,sizeof(par),comparar);
	}
	ip++;
	pp=&p0[ip];
	generar_numeros();
}

void siguiente() {
	int k;

	if(pp->rest>0) {
  	//hace que el anterior apunte al siguiente
		(pn->ant)->ste=pn->ste;
		(pn->ste)->ant=pn->ant;
		//apunta el puntero al siguiente, removiendo el numero de la lista circular
		pn=pn->ste;
	  //se mueve el puntero a traves de la lista hasta una posicion aleatoria	
		k=rand()%(pp->cant);
		while(k>0) {
			pn=pn->ste;
			k--;
		}
	}
	else {
		sgte_par();
	}
}

void iniciar_pares() {
	p0=malloc(sizeof(par));
	pp=p0;
	pp->r[0]=0;	pp->r[1]=0;
	pp->cant=1;
	pp->rest=1;
	npar=1;
}

void iniciar_numeros() {
	n0=malloc(sizeof(numero));
	pn=n0;
	pn->n=1;
	pn->ant=pn;	pn->ste=pn;
}

int sortear_k(float g, int kmin, int kmax) {
	long int s=-rand();
	float u=ran1(&s);

	return 
//		  (int) nearbyint(
		(int) floor(
					kmin*pow(
						( 1 - u*(1 - pow((float)kmax/(float)kmin,1-g))),
						1/(1-g)));
}


//para rozrand...
void crear_listas(int i, int *tam, int **l, int **nl, int num, float *m) {
	int c=0; //contador
	int j;
	int ucl=0, ucnl=0;

	for(j=0;j<num;j++) {
		c+=(i!=j?M(i,j):0);
	}
	(*tam)=c;

	// TIENEN QUE EMPEZAR COMO NULL!!!!
	if(*l!=NULL) {
		free(*l);
		free(*nl);
	}
	(*l)=malloc(c*sizeof(int));
	(*nl)=malloc((num-1-c)*sizeof(int));

	for(j=0;j<num;j++) {
		if(i!=j) {
			if(M(i,j)==1) {
				(*l)[ucl]=j;
				ucl++;
			} else {
				(*nl)[ucnl]=j;
				ucnl++;
			}
		}
	}
}

