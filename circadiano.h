/* Definiciones y funciones para correr el modelo de ecuaciones diferenciales
 * ordinarias de Bernard et al, utilizando el metodo Runge Kutta 4.
 *
 * El modelo fue propuesto para modelar la din\'amica de especies qu\'imicas
 * (ARN, proteinas y neuropeptidos) en neuronas del NSQ (nucleo supraquiasm\'a-
 * tico, el principal 'reloj' circadiano en mam\'iferos).
 * Se trata de 10*N variables acopladas, donde N es el n\'umero de neuronas.
 * Las 10 variables capturan, de forma simplificada, los compuestos que
 * intervienen en la generaci'on de ritmos circadianos a nivel gen\'etico.
 * La interconexi\'on entre neuronas est\'a dada mediante una matriz de
 * conectividad.
 * (El modelo fue propuesto originalmente (por Bernard et al) para sugerir y/o
 * explorar num\'ericamente la posibilidad de que oscilaciones amortiguadas
 * puedan ser sostenidas tras acoplar a los osciladores.
 * En nuestro caso, lo utilizamos para explorar el efecto de diferentes
 * topolog\'ias de red para conectar a las neuronas)
 *
 * Para mas informaci\'on, ver:
 *   Bernard, Samuel, et al. "Synchronization-induced rhythmicity of circadian
 *   oscillators in the suprachiasmatic nucleus." PLoS Comput Biol 3.4 (2007):
 *   e68.
 */

// parametros para tener oscilaciones amortiguadas con una sola celula
#define V1B  9.0
#define K1B  1.0
#define K1I  0.56
#define P    3
#define H    2
#define K1D  0.18
#define K2B  0.3
#define Q    2
#define K2D  0.1
#define K2T  0.36
#define K3T  0.02
#define K3D  0.18
#define V4B  1.0
#define K4B  2.16
#define R    3
#define K4D  1.1
#define K5B  0.24
#define K5D  0.09
#define K5T  0.45
#define K6T  0.06
#define K6D  0.18
#define K6A  0.09
#define K7A  0.003
#define K7D  0.13
#define K8   1.0
#define K8D  4.0
#define K    1.0
#define KX1  3.0
#define X1T  15.0
#define KDX1 4.0
#define KX2  0.25
#define X2T  15.0
#define KDX2 10.0
#define L0   0.2

// runge-kutta 4 (de numerical recipes in c)
void rk4(float z[], float dzdt[], int n, float t, float h, float q, float e,
		float zout[], void (*derivs)(float, float, float, float [], float [])) {
  int i;
  float th,hh,h6,*dzm,*dzt,*zt;
    
  dzm=malloc(sizeof(float)*(n));
  dzt=malloc(sizeof(float)*(n));
  zt=malloc(sizeof(float)*(n));

  hh=h*0.5;
  h6=h/6.0;
  th=t+hh;
  
  for (i=0;i<n;i++) zt[i]=z[i]+hh*dzdt[i];
  (*derivs)(th,q,e,zt,dzt);
  for (i=0;i<n;i++) zt[i]=z[i]+hh*dzt[i];
  (*derivs)(th,q,e,zt,dzm);
  for (i=0;i<n;i++) {
    zt[i]=z[i]+h*dzm[i];
    dzm[i] += dzt[i];
  }
  (*derivs)(t+h,q,e,zt,dzt);
  for (i=0;i<n;i++)
    zout[i]=z[i]+h6*(dzdt[i]+dzt[i]+2.0*dzm[i]);

  free(zt);
  free(dzt);
  free(dzm);
}

// la luz... (funcion en ecuacion 1 del modelo de Bernard et al)
float L(float t) {
	float tluz=12, tdark=12,m;

	m=fmod(t,tluz);
	if(t<24) return 0;
	else
	if (m+tluz<=tdark) {
		return L0*sin(PI*m/t);
	} else return 0;
}

// derivadas del modelo de Bernard et al
void dzdt(float t, float q, float e, float z[], float dz[]) {
  int i;
  float a1=z[6]+pow(z[9],H);
  float a2=pow(z[2]/K1I,P);
  float a3=pow(z[2],R);
  float k4br=pow(K4B,R);

  dz[0]=V1B*a1/(K1B*(1+a2)+a1)-K1D*z[0]/*+L(t)*/; // no olvidar L, p/ despues 
  dz[1]=K2B*pow(z[0],Q)-(K2D+K2T)*z[1]+K3T*z[2];
  dz[2]=K2T*z[1]-(K3T+K3D)*z[2];
  dz[3]=V4B*a3/(k4br+a3)-K4D*z[3];
  dz[4]=K5B*z[3]-(K5D+K5T)*z[4]+K6T*z[5];
  dz[5]=K5T*z[4]-(K6T+K6D+K6A)*z[5]+K7A*z[6];
  dz[6]=K6A*z[5]-(K7A+K7D)*z[6];
  // neurotransmisor VIP
  dz[7]=K8*z[1]-K8D*z[7];
  // cascada
  dz[8]=KX1*q*(X1T-z[8])-KDX1*z[8];
  dz[9]=KX2*z[8]*(X2T-z[9])-KDX2*z[9];;

  for(i=0;i<VARS;i++) dz[i]*=e;
}

// se inician los valores de las especies
// entre 0 y 2 veces su valor medio cuando oscilan
void iniciar_neuronas(neurona *n, int num) {
  int i,k;
  long s;

  float z0[] = {1.5, 2, 2, 0.5, 0.13, 0.3, 0.25, 0.3, 2.5,1};

  // iniciar generador de num aleatorios, para generar la semilla para generar
  // gaussian deviates
  srand(time(NULL));
  s=-(long)rand();

  // iniciamos concentraciones y e_i
  for(i=0;i<num;i++) {
    n[i].e=1/(1+gasdev(&s)*0.05);
    for(k=0;k<VARS;k++) {
      n[i].z[k]=z0[k]*2*ran1(&s);
    }
	}
}

