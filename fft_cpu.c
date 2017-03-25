/*****************************************************************************/
// gcc -O0 -o fft_cpu fft_cpu.c -lrt -lm

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#define GIG 1000000000
#define CPG 2.9           // Cycles per GHz -- Adjust to your computer

#define ITERS 1
#define DELTA 0
#define BASE 8

#define OPTIONS 7       // NEED TO MODIFY
#define IDENT 1.0
#define OP *

typedef float data_t;

typedef struct{
	float Re; 
	float Im;
} complex;

/*****************************************************************************/
main(int argc, char *argv[])
{
  int OPTION;

  struct timespec diff(struct timespec start, struct timespec end);
  struct timespec time1, time2;
  struct timespec time_stamp[OPTIONS][ITERS+1];
  int clock_gettime(clockid_t clk_id, struct timespec *tp);
 

  void init_complex(complex *v, int N);
  void init_array(data_t *v, int N);

  void radix2dit(complex *v, int n, complex *temp);
  void radix2dif(complex *v, int n);
  void radix4dit(complex *v, int n, complex *temp);
  void radix4dif(complex *v, int n);
  void radix8dit(complex *v, int n, complex *temp);
  void radix8dif(complex *v, int n);
  void radix2_non_recursive(float *data, int *nn);

  int i, j,q,N;
  long long int time_sec, time_ns;
  long int MAXSIZE = BASE+(ITERS+1)*DELTA;

  int *nn;
  nn = &N;

  printf("\n Hello World -- fft examples\n");

  // execute and time ... options from B&O 
  OPTION = 0;
  for (i = 0; i < ITERS; i++) {
    q = BASE+(i+1)*DELTA;
    N=(1<<q);
    complex v[N],scratch[N];
    init_complex(v,N);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    radix2dit(v, N, scratch);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
  }
  OPTION++;
  for (i = 0; i < ITERS; i++) {
    q = BASE+(i+1)*DELTA;
    N=(1<<q);
    complex v[N];
    init_complex(v,N);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    radix2dif(v, N);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
  }

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    q = BASE+(i+1)*DELTA;
    N=(1<<q);
    complex v[N],scratch[N];
    init_complex(v,N);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    radix4dit(v, N, scratch);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
  }

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    q = BASE+(i+1)*DELTA;
    N=(1<<q);
    complex v[N];
    init_complex(v,N);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    radix4dif(v, N);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
  }

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    q = BASE+(i+1)*DELTA;
    N=(1<<q);
    complex v[N],scratch[N];
    init_complex(v,N);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    radix8dit(v, N, scratch);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
  }

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    q = BASE+(i+1)*DELTA;
    N=(1<<q);
    complex v[N];
    init_complex(v,N);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    radix8dif(v, N);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
  }

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    q = BASE+(i+1)*DELTA;
    N=(1<<(q+1));
    data_t v[N];
    init_array(v,N);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    radix2_non_recursive(v, nn);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
  }

  /* output times */
  printf("\nlog2(size), r2dit, r2dif,  r4dit,  r4dif, r2nr");
  for (i = 0; i < ITERS; i++) {

    printf("\n%d, ", BASE+(i+1)*DELTA);
    for (j = 0; j < OPTIONS; j++) {
      if (j != 0) printf(", ");
      printf("%ld", (long int)((double)(CPG)*(double)
		 (GIG * time_stamp[j][i].tv_sec + time_stamp[j][i].tv_nsec)));
    }
  }

  printf("\n");
  
}/* end main */

/**********************************************/
/* Create array of specified length */

void init_complex(complex *v, int N)
{
  int k;
  for(k=0; k<N; k++) {		//loop inserts values into array of points
    v[k].Re = 0;	
    v[k].Im = 0;
  }
  v[0].Re=1;

}

void init_array(data_t *v, int N){
   int k;
   for(k=0; k<N; k++){
    v[k] = 0;
   }
}
/*************************************************/
struct timespec diff(struct timespec start, struct timespec end)
{
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}

/*************************************************/

//Cooley-Tukey FFT algorithm  
void radix2dit(complex *v, int n, complex *temp){		
  if(n>1) {			//return
    int k;	//incrementer    
    complex *vo,*ve;
    float real,imag,cosw,sinw;
    ve = temp; vo = temp+n/2;
    for(k=0; k<n/2; k++) {	//divide v into even vector (ve) and odd vector (vo)
      ve[k] = v[2*k];		
      vo[k] = v[2*k+1];
    }

    radix2dit(ve,n/2,v);		//recursively call fft for even indices of n
    radix2dit(vo,n/2,v);		//recursively call fft for odd indices of n

    for(k=0; k<n/2; k++) {

	/*when we divide into even and odd indices,
	 the odd terms are multiplied by a constant exp(-j*2*pi*k/n)
	 realconstant*vo = a*cosw + jb*sinw        imagconst*vo = b*cosw - ja*sinw	*/

	//v = ve + constant*vo
	cosw = cos(2*M_PI*k/n);
	sinw = sin(2*M_PI*k/n);
	real = vo[k].Re*cosw + vo[k].Im*sinw;
	imag = vo[k].Im*cosw - vo[k].Re*sinw;
	//v = ve + constant*vo
	v[k].Re = ve[k].Re + real;	
	v[k].Im = ve[k].Im + imag;	
	
	//v[k+n/2] = ve - w*vo
	v[k+n/2].Re = ve[k].Re - real;			
	v[k+n/2].Im = ve[k].Im - imag;			
    }
  }
  return;
}

/*************************************************/

//Cooley-Tukey FFT algorithm  
void radix2dif(complex *v, int n)
{
    int    k;
    complex w, v1[2];

    /** Do 2 Point DFT */
    for (k=0; k<n/2; k++)
    {
	/** Don't hurt the butterfly */

	w.Re=cos((double)k*2.0*M_PI/(double)n);
        w.Im=sin((double)k*2.0*M_PI/(double)n);
	v1[0].Re = (v[k].Re + v[k + n/2].Re);
	v1[0].Im = (v[k].Im + v[k + n/2].Im);
	v1[1].Re = (v[k].Re - v[k + n/2].Re) * w.Re + ((v[k].Im - v[k + n/2].Im) * w.Im); 
	v1[1].Im = (v[k].Im - v[k + n/2].Im) * w.Re - ((v[k].Re - v[k + n/2].Re) * w.Im);

	/** In-place results */
	v[k].Re = v1[0].Re;
	v[k].Im = v1[0].Im;
	v[k + n/2].Re = v1[1].Re;
	v[k + n/2].Im = v1[1].Im;

    }
    
    /** Don't recurse if we're down to one butterfly */
    if ((n/2)!=1){
	 radix2dif(&v[0], n/2);
	 radix2dif(&v[n/2], n/2);
    }
}
/*************************************************/

//Cooley-Tukey FFT algorithm  
void radix4dit(complex *v, int n, complex *temp){		
  if(n>1) {			//return
    int k;	//incrementer    
    complex *v1,*v2,*v3,*v4;
    float real,imag,cosw,sinw,t1,t1m,t2,t2m,t3,t3m,t4,t4m,cos2w,sin2w,cos3w,sin3w;
    v1 = temp; v2 = temp+n/4;	v3 = temp+n/2;	v4 = temp+3*n/4;
    for(k=0; k<n/4; k++) {	//divide v into even vector (ve) and odd vector (vo)
      v1[k] = v[4*k];		
      v2[k] = v[4*k+1];
      v3[k] = v[4*k+2];
      v4[k] = v[4*k+3];
    }

    radix4dit(v1,n/4,v);		//recursively call fft for even indices of n
    radix4dit(v2,n/4,v);		//recursively call fft for odd indices of n
    radix4dit(v3,n/4,v);		//recursively call fft for even indices of n
    radix4dit(v4,n/4,v);		//recursively call fft for odd indices of n

    for(k=0; k<n/4; k++) {

	/*when we divide into even and odd indices,
	 the odd terms are multiplied by a constant exp(-j*2*pi*k/n)
	 realconstant*vo = a*cosw + jb*sinw        imagconst*vo = b*cosw - ja*sinw	*/
	cosw = cos(2*M_PI*k/n);
	sinw = sin(2*M_PI*k/n);

	t1 = v1[k].Re + v3[k].Re*cosw + v3[k].Im*sinw;
	t1m = v1[k].Im + v3[k].Im*cosw - v3[k].Re*sinw;

	t2 =  v2[k].Re*cosw + v2[k].Im*sinw + v4[k].Re*cosw + v4[k].Re*sinw;
	t2m =  v2[k].Im*cosw - v2[k].Re*sinw + v4[k].Im*cosw - v4[k].Re*sinw;

	t3 = v1[k].Re - v3[k].Re*cosw - v3[k].Im*sinw;
	t3m = v1[k].Im - v3[k].Im*cosw + v3[k].Re*sinw;

	t4 =  v2[k].Im*cosw - v2[k].Re*sinw - v4[k].Im*cosw + v4[k].Re*sinw;
	t4m =  -v2[k].Re*cosw - v2[k].Im*sinw + v4[k].Re*cosw + v4[k].Im*sinw;

	v[k].Re = t1 + t2;
	v[k].Im = t1m + t2m;

	v[k+n/4].Re = t3 + t4;
	v[k+n/4].Im = t3m + t4m;

	v[k+n/2].Re = t1 - t2;
	v[k+n/2].Im = t1m - t2m;

	v[k+3*n/4].Re = t3 - t4;
	v[k+3*n/4].Im = t3m - t4m;

    }
  }
  return;
}

/*************************************************/
void radix4dif(complex *v, int n)
{ 
    int    k, m;
    complex w, v1[4];

    // printf("\n Here\n");
    /** Do 4 Point DFT */ 
    for (k=0; k<n/4; k++)
    {
	/** Don't hurt the butterfly */
	v1[0].Re = (v[k].Re + v[k + n/4].Re + v[k + n/2].Re + v[k + 3*n/4].Re);
	v1[0].Im = (v[k].Im + v[k + n/4].Im + v[k + n/2].Im + v[k + 3*n/4].Im);

	v1[1].Re = (v[k].Re + v[k + n/4].Re - v[k + n/2].Re - v[k + 3*n/4].Im);
	v1[1].Im = (v[k].Im - v[k + n/4].Im - v[k + n/2].Im + v[k + 3*n/4].Re);

	v1[2].Re = (v[k].Re - v[k + n/4].Re + v[k + n/2].Re - v[k + 3*n/4].Re);
	v1[2].Im = (v[k].Im - v[k + n/4].Im + v[k + n/2].Im - v[k + 3*n/4].Im);

	v1[3].Re = (v[k].Re - v[k + n/4].Im - v[k + n/2].Re + v[k + 3*n/4].Im);
	v1[3].Im = (v[k].Im + v[k + n/4].Re - v[k + n/2].Im - v[k + 3*n/4].Re);


	/** In-place results */
	for (m=0; m<4; m++)
	{
	    w.Re = cos((double)m*(double)k*2.0*M_PI/(double)n);
	    w.Im = sin((double)m*(double)k*2.0*M_PI/(double)n);
	    v[k + n/4*m].Re = v1[m].Re*w.Re + v1[m].Im*w.Im;
	    v[k + n/4*m].Im = v1[m].Im*w.Re - v1[m].Re*w.Im;

	}
    }
    
    /** Don't recurse if we're down to one butterfly */
    if ((n/4)!=1){
	    radix4dif(&v[0], n/4);
	    radix4dif(&v[n/4], n/4);
	    radix4dif(&v[n/2], n/4);
            radix4dif(&v[3*n/4], n/4);
    }
   return;
}


/*************************************************/

//Cooley-Tukey FFT algorithm  
void radix8dit(complex *v, int n, complex *temp){		
  if(n>1) {			//return
    int k;	//incrementer    
    complex *v1,*v2,*v3,*v4,*v5,*v6,*v7,*v8;
    float real,imag,cosw,sinw,t1,t1m,t2,t2m,t3,t3m,t4,t4m,t5,t5m,t6,t6m,t7,t7m,t8,t8m;
    v1 = temp; v2 = temp+n/8;	v3 = temp+n/4;	v4 = temp+3*n/8;    
    v5 = temp+n/2; v6 = temp+5*n/8;   v7 = temp+3*n/4;	v8 = temp+7*n/8;
    for(k=0; k<n/8; k++) {	//divide v into even vector (ve) and odd vector (vo)
      v1[k] = v[8*k];		
      v2[k] = v[8*k+1];
      v3[k] = v[8*k+2];
      v4[k] = v[8*k+3];
      v5[k] = v[8*k+4];
      v6[k] = v[8*k+5];
      v7[k] = v[8*k+6];
      v8[k] = v[8*k+7];
    }

    radix8dit(v1,n/8,v);		//recursively call fft for even indices of n
    radix8dit(v2,n/8,v);		//recursively call fft for odd indices of n
    radix8dit(v3,n/8,v);		//recursively call fft for even indices of n
    radix8dit(v4,n/8,v);		//recursively call fft for odd indices of n
    radix8dit(v5,n/8,v);		//recursively call fft for even indices of n
    radix8dit(v6,n/8,v);		//recursively call fft for odd indices of n
    radix8dit(v7,n/8,v);		//recursively call fft for even indices of n
    radix8dit(v8,n/8,v);		//recursively call fft for odd indices of n

    for(k=0; k<n/8; k++) {

	/*when we divide into even and odd indices,
	 the odd terms are multiplied by a constant exp(-j*2*pi*k/n)
	 realconstant*vo = a*cosw + jb*sinw        imagconst*vo = b*cosw - ja*sinw	*/
	cosw = cos(2*M_PI*k/n);
	sinw = sin(2*M_PI*k/n);


	t1 = v1[k].Re + v3[k].Re*cosw + v3[k].Im * sinw + v5[k].Re*cosw + v5[k].Im * sinw + v7[k].Re*cosw + v7[k].Im * sinw;
	t1m = v1[k].Im + v3[k].Im*cosw - v3[k].Re * sinw + v5[k].Im*cosw - v5[k].Re * sinw + v7[k].Im*cosw - v7[k].Re * sinw;

	t2 = v2[k].Re*cosw + v2[k].Im * sinw + v4[k].Re*cosw + v4[k].Im * sinw + v6[k].Re*cosw + v6[k].Im * sinw + v8[k].Re*cosw + v8[k].Im * sinw;
	t2m = v2[k].Im*cosw - v2[k].Re * sinw + v4[k].Im*cosw - v4[k].Re * sinw + v6[k].Im*cosw - v6[k].Re * sinw + v8[k].Im*cosw - v8[k].Re * sinw;
	
	t3 = v1[k].Re + v3[k].Im*cosw - v3[k].Re * sinw - v5[k].Re*cosw - v5[k].Im * sinw + v7[k].Re*sinw - v7[k].Im * cosw;
	t3m = v1[k].Im - v3[k].Re*cosw - v3[k].Im * sinw - v5[k].Im*cosw + v5[k].Re * sinw + v7[k].Re*cosw + v7[k].Im * sinw;

	t4 = v2[k].Re*cosw + v2[k].Im * sinw + v4[k].Im*cosw - v4[k].Re * sinw - v6[k].Re*cosw - v6[k].Im * sinw + v8[k].Re*sinw - v8[k].Im * cosw;
	t4m = v2[k].Im*cosw - v2[k].Re * sinw - v4[k].Re*cosw - v4[k].Im * sinw - v6[k].Im*cosw + v6[k].Re * sinw + v8[k].Re*cosw + v8[k].Im * sinw;

	t5 = v1[k].Re - v3[k].Re*cosw - v3[k].Im * sinw + v5[k].Re*cosw + v5[k].Im * sinw - v7[k].Re*cosw - v7[k].Im * sinw;
	t5m = v1[k].Im - v3[k].Im*cosw + v3[k].Re * sinw + v5[k].Im*cosw - v5[k].Re * sinw - v7[k].Im*cosw + v7[k].Re * sinw;

	t6 = v2[k].Im*cosw - v2[k].Re * sinw + v4[k].Re*sinw - v4[k].Im * cosw + v6[k].Im*cosw - v6[k].Re * sinw + v8[k].Re*sinw - v8[k].Im * cosw;
	t6m = -v2[k].Re*cosw - v2[k].Im * sinw + v4[k].Re*cosw - v4[k].Im * sinw - v6[k].Re*cosw - v6[k].Im * sinw + v8[k].Re*cosw + v8[k].Im * sinw;

	t7 = v1[k].Re + v3[k].Re*sinw - v3[k].Im * cosw - v5[k].Re*cosw - v5[k].Im * sinw + v7[k].Im*cosw - v7[k].Re * sinw;
	t7m = v1[k].Im + v3[k].Re*cosw + v3[k].Im * sinw - v5[k].Im*cosw + v5[k].Re * sinw - v7[k].Re*cosw - v7[k].Im * sinw;

	t8 = v2[k].Im*cosw - v2[k].Re * sinw + v4[k].Re*cosw + v4[k].Im * sinw + v6[k].Re*sinw - v6[k].Im * cosw - v8[k].Re*cosw - v8[k].Im * sinw;
	t8m = -v2[k].Re*cosw - v2[k].Im * sinw + v4[k].Im*cosw - v4[k].Re * sinw + v6[k].Re*cosw + v6[k].Im * sinw - v8[k].Im*cosw + v8[k].Re * sinw;

	v[k].Re = t1 + t2;
	v[k].Im = t1m + t2m;

	v[k+n/8].Re = t3 + t4;
	v[k+n/8].Im = t3m + t4m;

	v[k+n/4].Re = t5 + t6;
	v[k+n/4].Im = t5m + t6m;

	v[k+3*n/8].Re = t7 + t8;
	v[k+3*n/8].Im = t7m + t8m;

	v[k+n/2].Re = t1 - t2;
	v[k+n/2].Im = t1m - t2m;

	v[k+5*n/8].Re = t3 - t4;
	v[k+5*n/8].Im = t3 - t4m;

	v[k+3*n/4].Re = t5 - t6;
	v[k+3*n/4].Im = t5m - t6m;

	v[k+7*n/8].Re = t7 - t8;
	v[k+7*n/8].Im = t7m - t8m;


    }
  }
  return;
}

/*************************************************/
//Cooley-Tukey FFT algorithm  
void radix8dif(complex *v, int n)
{ 
    int    k, k1;
    complex w, v1[8];

 
    /** Do 4 Point DFT */ 
    for (k=0; k<n/8; k++)
    {
	/** Don't hurt the butterfly */
	v1[0].Re = (v[k].Re + v[k + n/8].Re + v[k + n/4].Re + v[k + 3*n/8].Re + v[k + n/2].Re + v[k + 5*n/8].Re + v[k + 3*n/4].Re + v[k + 7*n/8].Re);
	v1[0].Im = (v[k].Im + v[k + n/8].Im + v[k + n/4].Im + v[k + 3*n/8].Im + v[k + n/2].Im + v[k + 5*n/8].Im + v[k + 3*n/4].Im + v[k + 7*n/8].Im);

	v1[1].Re = (v[k].Re + v[k + n/8].Re + v[k + n/4].Re + v[k + 3*n/8].Re - v[k + n/2].Re - v[k + 5*n/8].Re - v[k + 3*n/4].Im - v[k + 7*n/8].Im);
	v1[1].Im = (v[k].Im + v[k + n/8].Im - v[k + n/4].Im - v[k + 3*n/8].Im - v[k + n/2].Im - v[k + 5*n/8].Im + v[k + 3*n/4].Re + v[k + 7*n/8].Re);

	v1[2].Re = (v[k].Re + v[k + n/8].Re - v[k + n/4].Re + v[k + 3*n/8].Re + v[k + n/2].Re + v[k + 5*n/8].Re - v[k + 3*n/4].Re - v[k + 7*n/8].Re);
	v1[2].Im = (v[k].Im + v[k + n/8].Im - v[k + n/4].Im + v[k + 3*n/8].Im + v[k + n/2].Im + v[k + 5*n/8].Im - v[k + 3*n/4].Im - v[k + 7*n/8].Im);

	v1[3].Re = (v[k].Re + v[k + n/8].Re - v[k + n/4].Im - v[k + 3*n/8].Im - v[k + n/2].Re - v[k + 5*n/8].Re + v[k + 3*n/4].Im + v[k + 7*n/8].Im);
	v1[3].Im = (v[k].Im + v[k + n/8].Im + v[k + n/4].Re + v[k + 3*n/8].Re - v[k + n/2].Im - v[k + 5*n/8].Im - v[k + 3*n/4].Re - v[k + 7*n/8].Re);
/////
	v1[4].Re = (v[k].Re + v[k + n/8].Re + v[k + n/4].Re + v[k + 3*n/8].Re + v[k + n/2].Re + v[k + 5*n/8].Re + v[k + 3*n/4].Re + v[k + 7*n/8].Re);
	v1[4].Im = (v[k].Im + v[k + n/8].Im + v[k + n/4].Im + v[k + 3*n/8].Im + v[k + n/2].Im + v[k + 5*n/8].Im + v[k + 3*n/4].Im + v[k + 7*n/8].Im);

	v1[5].Re = (v[k].Re + v[k + n/8].Re + v[k + n/4].Re + v[k + 3*n/8].Re - v[k + n/2].Re - v[k + 5*n/8].Re - v[k + 3*n/4].Im - v[k + 7*n/8].Im);
	v1[5].Im = (v[k].Im + v[k + n/8].Im - v[k + n/4].Im - v[k + 3*n/8].Im - v[k + n/2].Im - v[k + 5*n/8].Im + v[k + 3*n/4].Re + v[k + 7*n/8].Re);

	v1[6].Re = (v[k].Re + v[k + n/8].Re - v[k + n/4].Re + v[k + 3*n/8].Re + v[k + n/2].Re + v[k + 5*n/8].Re - v[k + 3*n/4].Re - v[k + 7*n/8].Re);
	v1[6].Im = (v[k].Im + v[k + n/8].Im - v[k + n/4].Im + v[k + 3*n/8].Im + v[k + n/2].Im + v[k + 5*n/8].Im - v[k + 3*n/4].Im - v[k + 7*n/8].Im);

	v1[7].Re = (v[k].Re + v[k + n/8].Re - v[k + n/4].Im - v[k + 3*n/8].Im - v[k + n/2].Re - v[k + 5*n/8].Re + v[k + 3*n/4].Im + v[k + 7*n/8].Im);
	v1[7].Im = (v[k].Im + v[k + n/8].Im + v[k + n/4].Re + v[k + 3*n/8].Re - v[k + n/2].Im - v[k + 5*n/8].Im - v[k + 3*n/4].Re - v[k + 7*n/8].Re);

	/** In-place results */
	for (k1=0; k1<8; k1++)
	{
	    w.Re = cos((double)k1*(double)k*2.0*M_PI/(double)n);
	    w.Im = sin((double)k1*(double)k*2.0*M_PI/(double)n);
	    v[k + n/8*k1].Re = v1[k1].Re*w.Re - v1[k1].Im*w.Im;
	    v[k + n/8*k1].Im = v1[k1].Im*w.Re + v1[k1].Re*w.Im;

	}
    }
    
    /** Don't recurse if we're down to one butterfly */
    if (n/8!=1){
	    radix8dif(&v[0], n/8);
	    radix8dif(&v[n/8], n/8);
	    radix8dif(&v[n/4], n/8);
            radix8dif(&v[3*n/8], n/8);
	    radix8dif(&v[n/2], n/8);
	    radix8dif(&v[5*n/8], n/8);
            radix8dif(&v[3*n/4], n/8);
            radix8dif(&v[7*n/8], n/8);
    }
   return;
}


/*************************************************/
 void radix2_non_recursive(float *data, int *nn)
{
  int n, mmax, m, j, i,p;
  float wtemp, wr, wpr, wpi, wi, theta, wpin,t,t1;   
  float tempr, tempi, datar, datai, data1r, data1i;
  //printf("%d\n",*nn);
  n = *nn;
  *nn/=2;
  j=0;
  for(i=0; i<n; i+=2){
      if(j>i){		//swap data[j] with data[i] and data[j+1] with data[i+1]
	t = data[j];
	t1 = data[j+1];
        data[j] = data[i];
	data[j+1] = data[i+1];
	data[i] = t;
	data[i+1] = t1;
      }
      m = *nn;
      while(m>=2 && j>=m){
	j-=m;
	m>>=1;
      }
      j+=m;
  }

  theta = -M_PI/2;
  wpin=0;	//sin(+-PI)
  for(mmax=2; n > mmax; mmax*=2){
     wpi = wpin;
     wpin = sin(theta);
     wpr = 1 - wpin * wpin - wpin * wpin;
     //cos(theta/2)
     theta /=2;
     wr = 1;
     wi = 0;
     for(m=0;m<mmax;m+=2){
	j=m+mmax;
	tempr = (float) wr * (data1r = data[j]);
	tempi = (float) wi * (data1i = data[j+1]);
	for(i=m;i<n-mmax*2;i+=mmax*2){
	    tempr -= tempi;
	    tempi = (float) wr * data1i + (float) wi * data1r;
	    data1r = data[j+mmax*2];
	    data1i = data[j+mmax*2+1];
	    data[i] = (datar = data[i]) + tempr;
	    data[i+1] = (datai = data[i+1]) + tempi;
	    data[j] = datar - tempr;
	    data[j+1] = datai - tempi;
	    tempr = (float) wr * data1r;
	    tempi = (float) wi * data1i;
	    j += mmax*2;
	}
	tempr -= tempi;
	tempi = (float) wr * data1i + (float) wi * data1r;
	data[i] = (datar = data[i]) + tempr;
	data[i+1] = (datai = data[i+1]) + tempi;
	data[j] = datar - tempr;
	data[j+1] = datai - tempi;
	wr = (wtemp = wr) * wpr - wi * wpi;
	wi = wtemp * wpi + wi * wpr;
      }
  }

  return;
}

