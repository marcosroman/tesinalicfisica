/*  
 */

void actualizar_registros(int num, neurona *n, float t,
	float *x0, float *x1, int *nregs, float ***reg,
	float *auxper, float *auxperfase,
	float *ppico, float *pvalle, float *upico, float *uvalle
	) {
	int i;
	float dif, fasei, dper;

	for(i=0;i<num;i++) {
		if(!(x1[i]==x0[i]&&x1[i]==n[i].z[0]&&x0[i]==n[i].z[0])&&(x1[i]!=0&&x0[i]!=0)) {
			// si los registros a tiempo t+h (z0),  (x1) t y t-h (x0) son distintos...
			dper=(x1[i]-x0[i])*(n[i].z[0]-x1[i]);
			// si cambia la derivada entre los 3 tiempos y ningun valor es igual, se cumplira que dper<0
			if(dper<0) {
				// nuevo registro, si ya hay algo registrado
				dif=x1[i]-x0[i];
				
				if(nregs[i]!=0) {
					(*reg)[i]=realloc((*reg)[i],2*sizeof(float)*(nregs[i]+1));
				}
				nregs[i]++;
				//para el registro enesimo de la neurona i, reg[i][2*(n-1)] almacena el tiempo y reg[i][2*(n-1)+1] la fase/PI
				(*reg)[i][(nregs[i]-1)*2] = t;
				(*reg)[i][(nregs[i]-1)*2+1] = (dif>0) ? 0.5 : 1.5;

				if(dif>0) {
					if(ppico[i]==0) {
						ppico[i]=x0[i];
					}
					upico[i]=x0[i];
				}
				
				if(dif<0) { // es un minimo
					if(pvalle[i]==0) {
						pvalle[i]=x0[i];
					}
					uvalle[i]=x0[i];
				}
			} else if(dper==0 && nregs[i]>0) { // la segunda condicion (nmins[i]>0) evita considerar el caso inicial (x0[i]=x1[i])
				/* si se cumple que dper==0, entonces
						x1==x0 y n.z0!=x0,x1
						o bien x1==n.z0 y x0!=x1,n.z0
					 
					 esto ocurrira siempre que tengamos valores de la forma:
					 t_(j)			v0
					 t_(j+1)		v1
					 t_(j+2)		v1
					 ...
					 t_(j+k)		v1
					 t_(j+k+1)	v2

					 definiremos el tiempo del pico como 0.5*(t_(j+k)-t_(j+1)).
					 considerando que la primera vez que se den las condiciones de arriba (dper==0) sera cuando x0=v0, x1=v1 y z0=v1, guardamos el tiempo en la variable auxper: auxper[i]=t
					 la siguiente vez sera cuando z0=v2 (ya con auxper!=0), entonces auxper=0.5*(auxper+t), y guardamos este tiempo para luego setear auxper=0 nuevamente
					 */
				if(auxper[i]==0) {
					if(x1[i]==n[i].z[0]) { // para q no hayan errores...
						auxper[i]=t;
						auxperfase[i]=(x0[i]<x1[i])?0.5:1.5;
					} /*else {
						printf("algun error en calculo de periodos? (neurona %d)\n",i);
					}*/
				} else {
					if(x1[i]==x0[i]) {
						dif=n[i].z[0]-x1[i];
						fasei = (dif<0) ? 0.5 : 1.5;
						if(auxperfase[i]==fasei) {
							auxper[i]=0.5*(auxper[i]+t);
							(*reg)[i]=realloc((*reg)[i],2*sizeof(float)*(nregs[i]+1));
							nregs[i]++;
							(*reg)[i][(nregs[i]-1)*2] = auxper[i];
							(*reg)[i][(nregs[i]-1)*2+1] = (dif<0) ? 0.5 : 1.5;
							if(dif>0) { // es un minimo
								uvalle[i]=x1[i];
							} else {
								upico[i]=x1[i];
							}
						}
						auxper[i]=0;
						auxperfase[i]=0;
					}
				}
			}
		}	
	}
}

void actualizar_registro(int num, neurona *n, float t,
	float *x0, float *x1, int *nregs, float ***reg,
	float *auxper, float *auxperfase
	) {
	int i;
	float dif, fasei, dper;

	for(i=0;i<num;i++) {
		if(!(x1[i]==x0[i]&&x1[i]==n[i].z[0]&&x0[i]==n[i].z[0])&&(x1[i]!=0&&x0[i]!=0)) {
			// si los registros a tiempo t+h (z0),  (x1) t y t-h (x0) son distintos...
			dper=(x1[i]-x0[i])*(n[i].z[0]-x1[i]);
			// si cambia la derivada entre los 3 tiempos y ningun valor es igual, se cumplira que dper<0
			if(dper<0) {
				// nuevo registro, si ya hay algo registrado
				dif=x1[i]-x0[i];
				
				if(nregs[i]!=0) {
					(*reg)[i]=realloc((*reg)[i],2*sizeof(float)*(nregs[i]+1));
				}
				nregs[i]++;
				//para el registro enesimo de la neurona i, reg[i][2*(n-1)] almacena el tiempo y reg[i][2*(n-1)+1] la fase/PI
				(*reg)[i][(nregs[i]-1)*2] = t;
				(*reg)[i][(nregs[i]-1)*2+1] = (dif>0) ? 0.5 : 1.5;
			} else if(dper==0 && nregs[i]>0) { // la segunda condicion (nmins[i]>0) evita considerar el caso inicial (x0[i]=x1[i])
			/* si se cumple que dper==0, entonces
					x1==x0 y n.z0!=x0,x1
					o bien x1==n.z0 y x0!=x1,n.z0
				 
				 esto ocurrira siempre que tengamos valores de la forma:
				 t_(j)			v0
				 t_(j+1)		v1
				 t_(j+2)		v1
				 ...
				 t_(j+k)		v1
				 t_(j+k+1)	v2

				 definiremos el tiempo del pico como 0.5*(t_(j+k)-t_(j+1)).
				 considerando que la primera vez que se den las condiciones de arriba (dper==0) sera cuando x0=v0, x1=v1 y z0=v1, guardamos el tiempo en la variable auxper: auxper[i]=t
				 la siguiente vez sera cuando z0=v2 (ya con auxper!=0), entonces auxper=0.5*(auxper+t), y guardamos este tiempo para luego setear auxper=0 nuevamente
				 */
				if(auxper[i]==0) {
					if(x1[i]==n[i].z[0]) { // para q no hayan errores...
						auxper[i]=t;
						auxperfase[i]=(x0[i]<x1[i])?0.5:1.5;
					}
				} else {
					if(x1[i]==x0[i]) {
						dif=n[i].z[0]-x1[i];
						fasei = (dif<0) ? 0.5 : 1.5;
						if(auxperfase[i]==fasei) {
							auxper[i]=0.5*(auxper[i]+t);
							(*reg)[i]=realloc((*reg)[i],2*sizeof(float)*(nregs[i]+1));
							nregs[i]++;
							(*reg)[i][(nregs[i]-1)*2] = auxper[i];
							(*reg)[i][(nregs[i]-1)*2+1] = (dif<0) ? 0.5 : 1.5;

											}

						auxper[i]=0;
						auxperfase[i]=0;
					}
				}
			}
		}	
	}
}
