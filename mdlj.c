/*
   Microcanonical Molecular Dynamics simulation of a Lennard-Jones fluid
   in a periodic boundary

   Initially written for the course CHE 800-002, Molecular 
   Simulation Spring 0304, Drexel University, Department
   of Chemical Engineering, Philadelphia

   compile using "gcc -o mdlj mdlj.c -lm"

Copyright (c) 2004, Cameron F. Abrams & Board of Trustees of Drexel University
Copyright (c) 2012, Christoph Junghans
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the copyright holders nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Writes the coordinates in XYZ format to the output stream fp.
   The integer "z" is the atomic number of the particles, required
   for the XYZ format. The array ix contains the number of x-dir
   periodic boundary crossings a particle has performed; thus,
   the "unfolded" coordinate is rx[i]+ix[i]*L. */
void xyz_out (FILE * fp,
	      double * rx, double * ry, double * rz,
	      double * vx, double * vy, double * vz,
	      int * ix, int * iy, int * iz, double L,
	      int N, int z, int put_vel, int unfold) {
  int i;

  fprintf(fp,"%i %i\n\n",N,put_vel);
  for (i=0;i<N;i++) {
    fprintf(fp,"%i %.8lf %.8lf %.8lf ",z,
	    rx[i]+(unfold?(ix[i]*L):0.0),
	    ry[i]+(unfold?(iy[i]*L):0.0),
	    rz[i]+(unfold?(iz[i]*L):0.0));
    if (put_vel)
      fprintf(fp,"%.8lf %.8lf %.8lf",vx[i],vy[i],vz[i]);
    fprintf(fp,"\n");
  }
}

int xyz_in (FILE * fp, double * rx, double * ry, double * rz,
	     double * vx, double * vy, double * vz,
	     int * N) {
  int i;
  int has_vel, dum;
  fscanf(fp,"%i %i\n\n",N,&has_vel);
  for (i=0;i<(*N);i++) {
    fscanf(fp,"%i %lf %lf %lf ",&dum,&rx[i],&ry[i],&rz[i]);
    if (has_vel) { // read velocities
      fscanf(fp,"%lf %lf %lf",&vx[i],&vy[i],&vz[i]);
    }
  }
  return has_vel;
}

/** \brief save current positions
 * \param[in] rx array of x-coordinate of the position of all particles
 * \param[in] ry array of y-coordinate of the position of all particles
 * \param[in] rz array of z-coordinate of the position of all particles
 * \param[in] N number of particles
 * \param[out] rx0 array of x-coordinate of the position of all particles
 * \param[out] ry0 array of y-coordinate of the position of all particles
 * \param[out] rz0 array of z-coordinate of the position of all particles
 */
void save_positions(double *rx,double *ry,double *rz,int N,
                    double *rx0,double *ry0,double *rz0) {
  int i;
  for (i=0;i<N;i++) {
    rx0[i]=rx[i];
    ry0[i]=ry[i];
    rz0[i]=rz[i];
  }
}

/** \brief calculates the max distance a particle has moved
 * \param[in] rx array of x-coordinate of the position of all particles
 * \param[in] ry array of y-coordinate of the position of all particles
 * \param[in] rz array of z-coordinate of the position of all particles
 * \param[in] N number of particles
 * \param[in] rx0 array of x-coordinate of the position of all particles
 * \param[in] ry0 array of y-coordinate of the position of all particles
 * \param[in] rz0 array of z-coordinate of the position of all particles
 * \return square of the max distance moved
 */
double max_moved_distance2(double *rx,double *ry,double *rz, int N,
                        double *rx0,double *ry0,double *rz0) {
  double max2=0;
  int i;
  for(i=0;i<N;i++) {
    double dx2= (rx[i]-rx0[i])*(rx[i]-rx0[i]) \
               +(ry[i]-ry0[i])*(ry[i]-ry0[i]) \
               +(rz[i]-rz0[i])*(rz[i]-rz0[i]);
    if (dx2>max2) {
      max2=dx2;
    }
  }
  return max2;
}

/** \brief calculates periodic distance between two particles
 * \param[in] i index of particle i
 * \param[in] j index of particle j
 * \param[in] rx array of x-coordinate of the position of all particles
 * \param[in] ry array of y-coordinate of the position of all particles
 * \param[in] rz array of z-coordinate of the position of all particles
 * \param[out] dx difference in x-direction
 * \param[out] dy difference in y-direction
 * \param[out] dz difference in z-direction
 * \param[in] L box length
 * \return square shorted distance of all periodic images
 */
double per_dist2(int i, int j,
    double *rx, double *ry, double *rz,
    double *dx,double *dy,double *dz,double L) {
  double  r2,hL=L/2.0;

  /* Periodic boundary conditions: Apply the minimum image
     convention; note that this is *not* used to truncate the
     potential as long as there an explicit cutoff. */

  *dx  = (rx[i]-rx[j]);
  *dy  = (ry[i]-ry[j]);
  *dz  = (rz[i]-rz[j]);

  if (*dx>hL)       *dx-=L;
  else if (*dx<-hL) *dx+=L;
  if (*dy>hL)       *dy-=L;
  else if (*dy<-hL) *dy+=L;
  if (*dz>hL)       *dz-=L;
  else if (*dz<-hL) *dz+=L;
  r2 = *dx* *dx + *dy* *dy + *dz* *dz;
  return r2;
}

/** \brief calculates Lennard-Jones interaction between two particles
 * \param[in] i index of particle i
 * \param[in] j index of particle j
 * \param[in] dx array of x-coordinate of the distance vetor
 * \param[in] dy array of y-coordinate of the distance vetor
 * \param[in] dz array of z-coordinate of the distance vetor
 * \param[in] r2 square of the distance vector
 * \param[in,out] fx array of x-coordinate of the force of all particles
 * \param[in,out] fy array of y-coordinate of the force of all particles
 * \param[in,out] fz array of z-coordinate of the force of all particles
 * \param[in] L box length
 * \param[in] rc2 square of the cutoff distance
 * \param[in] ecut energy varlue at the cutoff
 * \param[in,out] vir variable to add the virial contribute to
 * \param[in,out] e variable to add the energy contribute to
 */
void calc_lj ( int i, int j, double dx, double dy, double dz, double r2,
               double * fx, double * fy, double * fz, double L,
	       double rc2, double ecut, double * vir, double *e) {
  double r6i, f;

  if (r2<rc2) {
    r6i   = 1.0/(r2*r2*r2);
    *e   += 4*(r6i*r6i - r6i) - ecut;
    f     = 48*(r6i*r6i-0.5*r6i);
    fx[i] += dx*f/r2;
    fx[j] -= dx*f/r2;
    fy[i] += dy*f/r2;
    fy[j] -= dy*f/r2;
    fz[i] += dz*f/r2;
    fz[j] -= dz*f/r2;
    *vir += f;
  }
}

/** \brief contains all information about the neighbor list */
typedef struct{
  int *head; ///< point to first neighboring particles
  int *neighbors; ///<list of all neighboring paricles
  int size; ///< size of neighbors array
  int frequence; ///<update frequenz
  int update; ///< if list should be updated
  double cutoff2; ///< square of the list cutoff
  double skin2; ///<square of the list skin (list cutoff - lj cutoff)
} nblist_t;

/** \brief algorithm for computing forces and potential energy.
 * algorithm for computing forces and potential energy. The virial
 *  is also computed and returned in *vir.
 */
double total_e ( double * rx, double * ry, double * rz,
		 double * fx, double * fy, double * fz,
		 int N, double L, nblist_t *nblist,
		 double rc2, double ecor, double ecut, double * vir ) {
   int i,j,k;
   double e = 0.0;

   double dx,dy,dz,r2;

   *vir=0.0;

   /* Zero the forces */
   for (i=0;i<N;i++) {
     fx[i]=fy[i]=fz[i]=0.0;
   }

   if (!nblist) {
     /* do N square calculation */
     for (i=0;i<(N-1);i++) {
       for (j=i+1;j<N;j++) {
	 r2=per_dist2(i,j,rx,ry,rz,&dx,&dy,&dz,L);
         calc_lj(i,j,dx,dy,dz,r2,fx,fy,fz,L,rc2,ecut,vir,&e);
       }
     }
   }
   else if (nblist->update) {
     /* update neighbor list and calc lj */
     k=0;
     for (i=0;i<(N-1);i++) {
       nblist->head[i]=k;
       for (j=i+1;j<N;j++) {
	 r2=per_dist2(i,j,rx,ry,rz,&dx,&dy,&dz,L);
	 if (r2<nblist->cutoff2) {
	   if(k == nblist->size ) {
	     /* double the size */
	     nblist->size*=2;
	     nblist->neighbors=(int *)realloc(nblist->neighbors,nblist->size*sizeof(int));
	   }
	   nblist->neighbors[k]=j;
	   k++;
           calc_lj(i,j,dx,dy,dz,r2,fx,fy,fz,L,rc2,ecut,vir,&e);
	 }
       }
     }
     nblist->head[N-1]=k;
     nblist->update=0;
   }
   else {
     /* use existing neighbor list */
     for (i=0;i<(N-1);i++) {
       for (k=nblist->head[i];k<nblist->head[i+1];k++) {
	 j=nblist->neighbors[k];
	 r2=per_dist2(i,j,rx,ry,rz,&dx,&dy,&dz,L);
         calc_lj(i,j,dx,dy,dz,r2,fx,fy,fz,L,rc2,ecut,vir,&e);
       }
     }
   }
   return e+N*ecor;
}

/** \brief generate gaussian random number
 * generate a gaussian random number using the
 * Box-Muller method
 * \param[in] mean mean value of the gaussian
 * \param[in sigma standard deviation of the gaussian
 * \return gaussian random number
 */
double drand_gaussian(double mean, double sigma) {
  static int i=0;
  static double x2;
  /* http://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform */
  double r,theta,u1,u2,x1;
  if (i==0) {
    i=1-i;
    u1=1-drand48();
    u2=1-drand48();
    r=sqrt(-2*log(u1))*sigma;
    theta=2*M_PI*u2;
    x1 = r*cos(theta) + mean;
    x2 = r*sin(theta) + mean ;
    return(x1);
  } else {
    i=1-i;
    return(x2);
  }
}

/* Initialize particle positions by assigning them
   on a cubic grid, then scaling positions
   to achieve a given box size and thereby, volume,
   and density */
void init ( double * rx, double * ry, double * rz,
	    double * vx, double * vy, double * vz,
	    int * ix, int * iy, int * iz,
	    int * N, double L, double T0,
	    double * KE, char * icf) {
  int i,iix,iiy,iiz;
  double cmvx=0.0,cmvy=0.0,cmvz=0.0;
  double T, fac;
  int n3=2;
  int vel_ok=0;

  /* If icf has a value, assume it is the name of a file containing
     the input configuration in XYZ format */
  if (icf) {
    FILE * fp = fopen(icf,"r");
    if (fp) {
      vel_ok = xyz_in(fp,rx,ry,rz,vx,vy,vz,N);
    }
    else {
      fprintf(stderr,"# error: could not read %s\n",icf);
      exit(-1);
    }
  }
  /* Assign particles on a cubic lattice */
  else {

    /* Find the lowest perfect cube, n3, greater than or equal to the
       number of particles */
    while ((n3*n3*n3)<*N) n3++;

    iix=iiy=iiz=0;
    /* Assign particle positions */
    for (i=0;i<*N;i++) {
      rx[i] = ((double)iix+0.5)*L/n3;
      ry[i] = ((double)iiy+0.5)*L/n3;
      rz[i] = ((double)iiz+0.5)*L/n3;
      iix++;
      if (iix==n3) {
	iix=0;
	iiy++;
	if (iiy==n3) {
	  iiy=0;
	  iiz++;
	}
      }
    }
  }
  /* If no velocities yet assigned, randomly pick some */
  if (!vel_ok) {
    for (i=0;i<*N;i++) {
      vx[i]=drand_gaussian(0.0,1.0);
      vy[i]=drand_gaussian(0.0,1.0);
      vx[i]=drand_gaussian(0.0,1.0);
    }
  }
  /* Take away any center-of-mass drift; compute initial KE */
  for (i=0;i<*N;i++) {
    cmvx+=vx[i];
    cmvy+=vy[i];
    cmvz+=vz[i];
  }
  (*KE)=0;
  for (i=0;i<*N;i++) {
    vx[i]-=cmvx/(*N);
    vy[i]-=cmvy/(*N);
    vz[i]-=cmvz/(*N);
    (*KE)+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  (*KE)*=0.5;
  /* if T0 is specified, scale velocities */
  if (T0>0.0) {
    T=(*KE)/(*N)*2./3.;
    fac=sqrt(T0/T);
    (*KE)=0;
    for (i=0;i<*N;i++) {
      vx[i]*=fac;
      vy[i]*=fac;
      vz[i]*=fac;
      (*KE)+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
    (*KE)*=0.5;
  }
  /* Initialize periodic boundary crossing counter arrays */
  for (i=0;i<*N;i++) {
    ix[i]=0;
    iy[i]=0;
    iz[i]=0;
  }
}

int main ( int argc, char * argv[] ) {

  double * rx, * ry, * rz;
  double * vx, * vy, * vz;
  double * fx, * fy, * fz;
  int * ix, * iy, * iz;
  int N=216;
  double L=0.0;
  double rho=pow(1./1.05,3), rc2 = 1.12246, vir, vir_old, V;
  double PE, KE, TE, ecor, ecut, T0=1.0, TE0;
  double rr3,dt=0.001, dt2;
  int i,s;
  int nSteps = 10, fSamp=100;
  int use_e_corr=0;
  int unfold = 0;

  char fn[20];
  FILE * out;
  char * wrt_code_str = "w";
  char * init_cfg_file = NULL;
  double rlist=0;
  int nblist_frequenz=0;

  unsigned long int Seed = 23410981;

  /* Here we parse the command line arguments;  If
   you add an option, document it in the usage() function! */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dt")) dt=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rc")) rc2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-ns")) nSteps = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-T0")) T0=atof(argv[++i]);
    else if (!strcmp(argv[i],"-fs")) fSamp=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-sf")) wrt_code_str = argv[++i];
    else if (!strcmp(argv[i],"-icf")) init_cfg_file = argv[++i];
    else if (!strcmp(argv[i],"-ecorr")) use_e_corr = 1;
    else if (!strcmp(argv[i],"-seed")) Seed = (unsigned long)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-rlist")) rlist = atof(argv[++i]);
    else if (!strcmp(argv[i],"-uplist")) nblist_frequenz = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-uf")) unfold = 1;
    else if ((!strcmp(argv[i],"-h"))||(!strcmp(argv[i],"--help"))) {
      fprintf(stdout,"mdlj usage:\n");
      fprintf(stdout,"mdlj [options]\n\n");
      fprintf(stdout,"Options:\n");
      fprintf(stdout,"\t -N [integer]\t\tNumber of particles (default %i)\n",N);
      fprintf(stdout,"\t -rho [real]\t\tNumber density (default %f)\n",rho);
      fprintf(stdout,"\t -dt [real]\t\tTime step (default %f)\n",dt);
      fprintf(stdout,"\t -rc [real]\t\tCutoff radius (default %f)\n",rc2);
      fprintf(stdout,"\t -ns [real]\t\tNumber of integration steps (default %i)\n",nSteps);
      fprintf(stdout,"\t -T0 [real]\t\tInitial temperature (default %f)\n",T0);
      fprintf(stdout,"\t -fs [integer]\t\tSample frequency (default %i)\n",fSamp);
      fprintf(stdout,"\t -sf [a|w]\t\tAppend or write config output file (default %s)\n",wrt_code_str);
      fprintf(stdout,"\t -icf [string]\t\tInitial configuration file (default %s)\n",init_cfg_file);
      fprintf(stdout,"\t -seed [integer]\tRandom number generator seed (default %li)\n",Seed);
      fprintf(stdout,"\t -rlist [real]\t\tUse Verlet lists with this cutoff (default %f)\n",rlist);
      fprintf(stdout,"\t -uplist [int]\t\tVerlet lists update frequence with 0=auto (default %i)\n",nblist_frequenz);
      fprintf(stdout,"\t -uf          \t\tPrint unfolded coordinates in output files (default %i)\n",unfold);
      fprintf(stdout,"\t -h           \t\tPrint this info\n");
      exit(0);
    }
    else {
      fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n",
	      argv[i]);
      exit(-1);
    }
  }

  /* Compute the side-length */
  L = pow((V=N/rho),0.3333333);

  /* Compute the tail-corrections; assumes sigma and epsilon are both 1 */
  rr3 = 1.0/(rc2*rc2*rc2);
  ecor = use_e_corr?8*M_PI*rho*(rr3*rr3*rr3/9.0-rr3/3.0):0.0;
  ecut = 4*(rr3*rr3*rr3*rr3-rr3*rr3);

  /* Compute the *squared* cutoff, reusing the variable rc2 */
  rc2*=rc2;

  /* compute the squared time step */
  dt2=dt*dt;

  /* Output some initial information */
  fprintf(stdout,"# NVE MD Simulation of a Lennard-Jones fluid\n");
  fprintf(stdout,"# L = %.5lf; rho = %.5lf; N = %i; rc = %.5lf\n",
	  L,rho,N,sqrt(rc2));
  fprintf(stdout,"# nSteps %i, seed %li, dt %.5lf\n",
nSteps,Seed,dt);

  nblist_t *nblist=NULL;
  double *rx0,*ry0,*rz0;
  /* verlet list stuff */
  if (rlist > 0) {
    nblist=(nblist_t *)malloc(sizeof(nblist_t));
    nblist->cutoff2 = rlist*rlist;
    if (nblist->cutoff2 < rc2) {
      fprintf(stderr,"-rlist must be bigger than -rc\n");
      exit(1);
    }
    nblist->skin2 = pow(rlist-sqrt(rc2),2);
    fprintf(stdout,"# using Verlet lists with cutoff %.5f generated by %s search with update frequence %i\n",rlist,"grid",nblist_frequenz);
    /* assume 25 neighbors and realloc later if needed */
    nblist->size= 25*N;
    nblist->neighbors = (int *)malloc(nblist->size*sizeof(int));
    nblist->head = (int *)malloc(N*sizeof(int));
    /* always update nblist one time */
    nblist->update=1;
    /* old position saver*/
    rx0 = (double*)malloc(N*sizeof(double));
    ry0 = (double*)malloc(N*sizeof(double));
    rz0 = (double*)malloc(N*sizeof(double));
  }

  /* Seed the random number generator */
  srand48(Seed);

  /* Allocate the position arrays */
  rx = (double*)malloc(N*sizeof(double));
  ry = (double*)malloc(N*sizeof(double));
  rz = (double*)malloc(N*sizeof(double));

  /* Allocate the boundary crossing counter arrays */
  ix = (int*)malloc(N*sizeof(int));
  iy = (int*)malloc(N*sizeof(int));
  iz = (int*)malloc(N*sizeof(int));

  /* Allocate the velocity arrays */
  vx = (double*)malloc(N*sizeof(double));
  vy = (double*)malloc(N*sizeof(double));
  vz = (double*)malloc(N*sizeof(double));

  /* Allocate the force arrays */
  fx = (double*)malloc(N*sizeof(double));
  fy = (double*)malloc(N*sizeof(double));
  fz = (double*)malloc(N*sizeof(double));

  /* Generate initial positions on a cubic grid,
     and measure initial energy */
  init(rx,ry,rz,vx,vy,vz,ix,iy,iz,&N,L,T0,&KE,init_cfg_file);
  sprintf(fn,"%i.xyz",0);
  out=fopen(fn,"w");
  xyz_out(out,rx,ry,rz,vx,vy,vz,ix,iy,iz,L,N,16,1,unfold);
  fclose(out);

  if (nblist) {
    save_positions(rx,ry,rz,N,rx0,ry0,rz0);
  }

  PE = total_e(rx,ry,rz,fx,fy,fz,N,L,nblist,rc2,ecor,ecut,&vir_old);
  TE0=PE+KE;

  fprintf(stdout,"# step PE KE TE drift T P\n");

  for (s=0;s<nSteps;s++) {

    /* First integration half-step */
    for (i=0;i<N;i++) {
      rx[i]+=vx[i]*dt+0.5*dt2*fx[i];
      ry[i]+=vy[i]*dt+0.5*dt2*fy[i];
      rz[i]+=vz[i]*dt+0.5*dt2*fz[i];
      vx[i]+=0.5*dt*fx[i];
      vy[i]+=0.5*dt*fy[i];
      vz[i]+=0.5*dt*fz[i];
      /* Apply periodic boundary conditions */
      if (rx[i]<0.0) { rx[i]+=L; ix[i]--; }
      if (rx[i]>L)   { rx[i]-=L; ix[i]++; }
      if (ry[i]<0.0) { ry[i]+=L; iy[i]--; }
      if (ry[i]>L)   { ry[i]-=L; iy[i]++; }
      if (rz[i]<0.0) { rz[i]+=L; iz[i]--; }
      if (rz[i]>L)   { rz[i]-=L; iz[i]++; }
    }

    if (nblist) {
      /* auto nblist update */
      if (nblist->frequence==0) {
	/*2* dx > skin */
	if (4*max_moved_distance2(rx,ry,rz,N,rx0,ry0,rz0) > nblist->skin2) {
	  nblist->update=1;
          save_positions(rx,ry,rz,N,rx0,ry0,rz0);
	}
      }
      else if (s%nblist->frequence==0) {
        nblist->update=1;
      }
    }

    /* Calculate forces */
    PE = total_e(rx,ry,rz,fx,fy,fz,N,L,nblist,rc2,ecor,ecut,&vir);

    /* Second integration half-step */
    KE = 0.0;
    for (i=0;i<N;i++) {
      vx[i]+=0.5*dt*fx[i];
      vy[i]+=0.5*dt*fy[i];
      vz[i]+=0.5*dt*fz[i];
      KE+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
    KE*=0.5;
    TE=PE+KE;
    fprintf(stdout,"%i %.5lf %.5lf %.5lf %.5lf %.5le %.5lf %.5lf\n",
	    s,s*dt,PE,KE,TE,(TE-TE0)/TE0,KE*2/3./N,rho*KE*2./3./N+vir/3.0/V);
    if (!(s%fSamp)) {
      sprintf(fn,"%i.xyz",!strcmp(wrt_code_str,"a")?0:s);
      out=fopen(fn,wrt_code_str);
      xyz_out(out,rx,ry,rz,vx,vy,vz,ix,iy,iz,L,N,16,1,unfold);
      fclose(out);
    }
  }
  return(0);
}
