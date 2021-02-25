/* Funciones utiles para modificar|manipular&visualizar la matriz de conectividad
 * (que define la red de interacciones entre 'neuronas')
 */

float f(float x2, float dm2) {
  return exp(-x2/dm2);
}

float dist2(int p1[], int p2[]) {
  float s=0, r;
  int i;

  for(i=0;i<DIM;i++) {
    r = p1[i]-p2[i];
    s += r*r;
  }

  return s;
}

float mindist2(int sqrn, int p1[], int p2[]) {
	float s=0;
	int hsn=sqrn/2;
  int i,r;

  for(i=0;i<DIM;i++) {
		r = abs(p1[i]-p2[i]);
		if(r>hsn) {
		 r-=sqrn;
		}
		s+=r*r;
	}

	return s;
}

void imprimir_red(int num, float *m) {
	int i,j;

	for(i=0;i<num;i++) {
		for(j=0;j<num;j++) {
			printf("%e ",M(i,j));
		}
		printf("\n");
	}
}

void normalizar(int num, float *m) {
	int i,j;
	float z;

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

void denormalizar(int num, float *m) {
	int i;

	for(i=0;i<num*num;i++) {
		if(m[i]>0) m[i]=1;
	}
}

void imprimir_C(int num, float *m) {
	int i,j;

	for(i=0;i<num;i++) {
		for(j=0;j<num;j++) {
			printf("%d ",M(i,j)>0?1:0);
		}
		//printf("\n");
	}
}
	
void imprimir_M(int num, float *m) {
	int i,j;

	for(i=0;i<num;i++) {
		for(j=0;j<num;j++) {
			printf("%e ",M(i,j));
		}
		printf("\n");
	}
}

float conectividad(int num, float *m) {
	int i;
	int numnum=num*num;
	int con=0;

	for(i=0;i<numnum;i++) {
		if(m[i]>0) con++;
	}

	return ((float)con)/((float)num*num);
}

// (asumimos que la red es cuadrada)
float longitud(int num, float *m) {
	int i,j;
	int sqrn=(int)sqrt(num);
  int pos[num][DIM];
	float l=0;

	for(i=0;i<num;i++) {
    pos[i][0]=i%sqrn;
    pos[i][1]=(i-pos[i][0])/sqrn;
  }

	for(i=0;i<num;i++) {
		for(j=0;j<num;j++) {
			if(M(i,j)>0) l+=sqrt(mindist2(sqrn,pos[i],pos[j]));
		}
	}

	return l;
}

void imprimir_distribucion(int num, float *m, FILE *f) {
	//M(i,j) > 0 si el nodo j -> nodo i
	float *d_in, *d_out, s;
	int i,j,k;

	// para bins...
	d_in=calloc(num+1,sizeof(float));
	d_out=calloc(num+1,sizeof(float));

	// para p_in (k)
	for(i=0;i<num;i++) {
		k=0;
		for(j=0;j<num;j++) {
			k+=M(i,j)>0?1:0;
		}
		d_in[k]++;
	}

	// para p_out (k)
	for(i=0;i<num;i++) {
	
		k=0;
		for(j=0;j<num;j++) {
			k+=M(j,i)>0?1:0;
		}
		d_out[k]++;
	}

	// normalizar
	s=0;
	for(k=0;k<=num;k++) s+=d_in[k];
	for(k=0;k<=num;k++) d_in[k]/=s;
	s=0;
	for(k=0;k<=num;k++) s+=d_out[k];
	for(k=0;k<=num;k++) d_out[k]/=s;

	// imprimir
	fprintf(f,"#in\n");
	for(k=0;k<=num;k++) {
		fprintf(f,"%d\t%e\n",k,d_in[k]);
	}
	fprintf(f,"\n\n#out\n");
	for(k=0;k<=num;k++) {
		fprintf(f,"%d\t%e\n",k,d_out[k]);
	}

	free(d_in);
	free(d_out);
}

void imprimir_distribucion_media(int num, float *m, int nsamples, FILE *f,
	 	void (*fred)(int num, float *m, float param[]), float *pars) {
		//M(i,j) > 0 si el nodo j -> nodo i
	long double *d_inind, *d_insum, s;
	int i,j,k,l;

	// para bins... desde 1 hasta num
	d_inind=calloc(num,sizeof(long double));
	d_insum=calloc(num,sizeof(long double));

	for(l=0;l<nsamples;l++) {
		if(l%100==0) printf("%g%%...",((float)l/(float)nsamples)*100);
		// nueva red con los mismos parametros...
		for(k=0;k<num*num;k++) m[k]=0;
		for(k=0;k<num;k++) d_inind[k]=0;

		(*fred)(num,m,pars);

		// para p_in (k)
		for(i=0;i<num;i++) {
			k=0;
			for(j=0;j<num;j++) {
				k+=M(i,j)>0?1:0;
			}
			//k=1 va en el primer registro, k=num en el ultimo 
			d_inind[k-1]++;
		}
		
		s=0;
		for(k=0;k<num;k++) s+=d_inind[k];
		for(k=0;k<num;k++) d_inind[k]/=s;

		for(k=0;k<num;k++) d_insum[k]+=d_inind[k];
	}
	
	s=0;
	for(k=0;k<num;k++) s+=d_insum[k];
	for(k=0;k<num;k++) d_insum[k]/=s;

		// para p_out (k)


	// normalizar

	// imprimir
	fprintf(f,"#in\n");
	for(k=0;k<num;k++) {
		fprintf(f,"%d\t%Le\n",k+1,d_insum[k]);
	}

	free(d_inind);
	free(d_insum);
	//free(d_out);
}

