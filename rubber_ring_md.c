// Two Rubber rings hitting each other Problem

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NY      100
#define NX      100
#define DY      1.00
#define N       5116
#define NB      2558
#define RHO0    1.1547
#define MASSP   1.0
#define MASSB   1.0
#define VRING   0.10//-04.0
#define SIGMA   0.80
#define EPS     100.0
#define RADIUS  6.0
#define MAX     400.00
#define DT      0.01
#define r       1.4142135
#define rad     30.00
#define M       2.0
#define NPRINT  10
#define PI      3.143
#define ri      30.00
#define ro      40.00


double X[N], Y[N], VX[N], VY[N], MASS[N];
double XN[N], YN[N], VXN[N], VYN[N];
double E[N];
double KE[N], PE[N];
double XDOT[N], YDOT[N], VXDOT[N], VYDOT[N];
double VXDOTOLD[N], VYDOTOLD[N];
double FX[N], FY[N];
double SXX[N], SXY[N], SYX[N], SYY[N];

FILE *fp;
FILE *fp5;

void initialize();
void derivatives();
void Update(double dt);
void callprint(int loopcounter);
void Energy();

int main()  {
    double t=0.0;
    int i, loopcounter=0;
    fp = fopen("Trajectory.lammpstrj","w");
    fp5 = fopen("Energy.output","w");

    initialize();
    while(t<MAX)    {

        //Incorporating Velocity Verlet
        for(i=0;i<N;i++)     {
            X[i] += VX[i]*DT + 0.5*VXDOT[i]*DT*DT;
            Y[i] += VY[i]*DT + 0.5*VYDOT[i]*DT*DT;
            VXDOTOLD[i] = VXDOT[i];
            VYDOTOLD[i] = VYDOT[i];
        }
        derivatives();

        for(i=0;i<N;i++)     {
            VX[i] += 0.5*(VXDOTOLD[i] + VXDOT[i])*DT;
            VY[i] += 0.5*(VYDOTOLD[i] + VYDOT[i])*DT;
        }

        Energy();

        if(loopcounter%NPRINT==0) {
            double KE1 = 0.0,PE1 = 0.0;
            printf("Time(ms): %lf ", (double) DT*loopcounter);
           for(i=0; i<N; i++) {
                KE1 += 0.5 * MASS[i] * ( VX[i]*VX[i] + VY[i]*VY[i] );
                PE1 +=  PE[i];
            }
            fprintf(fp5,"%lf %lf %lf %lf\n",(double) DT*loopcounter,KE1,PE1,KE1+PE1);

            callprint(loopcounter);
        }
        t+=DT;
        loopcounter++;
    }
    fclose(fp5);
    return(0);
}


void initialize()   {

    int n,i,j,count=0;
    double dy = DY, dx = sqrt(3.0)/2.0*dy, displacement = 0.5*DY;
    double xlim = (double) NY*dx, ylim = (double) NX*dy*0.30;
    double tot_momentum = (double) MASSB*VRING*N, tot_mass = 0.0;
    double temp,x, y, val, radi, rado;
    int tempint;

    FILE *fp2;
    fp2 = fopen("Initial.dat","w");


    //FILE *fp4;
    //fp4 = fopen("Consistency_Check.dat","r");

    //This is for left ring Particles
    for(j=0;j<NY;j++)   {
        for(i=0;i<NX;i++)   {
            if(j%2==0)  {
                x = (double) j*dx - 0.5*xlim;
                y = -(double) i*dy + 40.00;
            }
            if(j%2==1)  {
                x = (double) j*dx - 0.5*xlim;
                y = -(double) i*dy - displacement + 40.00;
            }
            val = x*x + y*y;
            radi = ri*ri;
            rado = ro*ro;
            if(val >= radi && val <= rado){
               VX[count] = VRING;
               VY[count] = 0.0;
               E[count] = 0.0;
               MASS[count] = MASSP;
               X[count] = x;
               Y[count] = y;
               //printf("count:%d \n",count);
               count++;
            }
        }
    }
    //This is for right ring Particle
 for(j=0;j<NY;j++)   {
        for(i=0;i<NX;i++)   {
            if(j%2==0)  {
                x = (double) j*dx - 0.5*xlim;
                y = -(double) i*dy + 40.00;
            }
            if(j%2==1)  {
                x = (double) j*dx - 0.5*xlim;
                y = -(double) i*dy - displacement + 40.00;
            }
            val = x*x + y*y;
            radi = ri*ri;
            rado = ro*ro;
            if(val >= radi && val <= rado){
               VX[count] = -VRING;
               VY[count] = 0.0;
               E[count] = 0.0;
               MASS[count] = MASSP;
               X[count] = 87.00 + x;
               Y[count] = y;
               //fprintf(fp2,"%lf %lf \n",X[count],Y[count]);
               count++;
            }
        }
    }
    //Read Mass from SPH data file constency_check.dat
    //for(i=0;i<N;i++)  fscanf(fp4, "%d %lf %lf %lf %lf", & tempint, &temp, &temp, &MASS[i], &temp);
    tot_mass = tot_momentum = 0.0;
    for(i=0;i<N;i++) tot_mass += MASS[i];
    for(i=0;i<N;i++) tot_momentum += MASS[i]*VX[i];
    printf("Total Momentum is: %lf \n",tot_momentum);
    printf("Total Mass is: %lf \n",tot_mass);

    //Initialize Velocity of the Particles so that the net momentum of the system is zero
     //for(i=0;i<N-NB;i++) VY[i] = (double) -1.0*tot_momentum/(tot_mass);
    fclose(fp2);
}

void derivatives()  {
    double xij, yij, vxij, vyij, dist, rsq = r*r;
    double term1, term12, term2,term3,term4,term5,fpair;
    int i,j;

    for(i=0;i<N;i++)  {
    	VXDOT[i]    = VYDOT[i] = 0.0;
    	SXX[i] = 0.0;
        SXY[i] = 0.0;
        SYX[i] = 0.0;
        SYY[i] = 0.0;
    }


	for(i=0;i<N;i++)	{  	
		for(j=0;j<N;j++)	{
			if(i!=j)	{
				//Ring1 with Ring1
				if((i<N-NB)&&(j<N-NB))	{
					xij  = X[i] - X[j];
                    yij  = Y[i] - Y[j];
                    vxij = VX[i] - VX[j];
                    vyij = VY[i] - VY[j];

                    dist = xij*xij + yij*yij;
                    if(dist<rsq)   {
                        term1       = 4.0*M;
                        term2       = 2.0-dist;
                        term3       = pow(term2,2.0*M-1.0);
                        term4       = pow(term2,M-1.0);
                        term5       = term3 - term4;
                        fpair       = term1*term5;

                        VXDOT[i]    += fpair*xij/MASS[i]; //Computing the Force on ith particle due to jth particle along x axis
                        VYDOT[i]    += fpair*yij/MASS[i];

                  		SXX[i] -= 0.5*(X[i]*fpair*xij - X[j]*fpair*xij);
                  		SXY[i] -= 0.5*(X[i]*fpair*yij - X[j]*fpair*yij);
                  		SYX[i] -= 0.5*(Y[i]*fpair*xij - Y[j]*fpair*xij);
                  		SYY[i] -= 0.5*(Y[i]*fpair*yij - Y[j]*fpair*yij);  
                  	}
				}
				//Ring2 with Ring2
				if((i>=N-NB)&&(j>=N-NB))	{
					xij  = X[i] - X[j];
                    yij  = Y[i] - Y[j];
                    vxij = VX[i] - VX[j];
                    vyij = VY[i] - VY[j];

                    dist = xij*xij + yij*yij;
                    if(dist<rsq)   {
                        term1       = 4.0*M;
                        term2       = 2.0-dist;
                        term3       = pow(term2,2.0*M-1.0);
                        term4       = pow(term2,M-1.0);
                        term5       = term3 - term4;
                        fpair       = term1*term5;

                        VXDOT[i]    += fpair*xij/MASS[i]; //Computing the Force on ith particle due to jth particle along x axis
                        VYDOT[i]    += fpair*yij/MASS[i];

                  		SXX[i] -= 0.5*(X[i]*fpair*xij - X[j]*fpair*xij);
                  		SXY[i] -= 0.5*(X[i]*fpair*yij - X[j]*fpair*yij);
                  		SYX[i] -= 0.5*(Y[i]*fpair*xij - Y[j]*fpair*xij);
                  		SYY[i] -= 0.5*(Y[i]*fpair*yij - Y[j]*fpair*yij);  
                  	}
				}				
				//Ring1 with Ring2
				if((i<N-NB)&&(j>=N-NB))	{
					xij  = X[i] - X[j];
                    yij  = Y[i] - Y[j];
                    vxij = VX[i] - VX[j];
                    vyij = VY[i] - VY[j];
                    dist = xij*xij + yij*yij;
                    dist = sqrt(dist);
                    term1 = (dist - RADIUS)/SIGMA;
                    if(term1 < 1.0) {
                        term12 = term1*term1;
                        term2 = 8.0*EPS/SIGMA*term1*(1-term12)*(1-term12)*(1-term12);
                        VXDOT[i] += term2*xij/dist/MASS[i]; //VXDOT[N-1]
                        VYDOT[i] += term2*yij/dist/MASS[i]; //VXDOT[N-1]

            			SXX[i] -= 0.5*(X[i]*term2*xij/dist - X[j]*term2*xij/dist);
	            		SXY[i] -= 0.5*(X[i]*term2*yij/dist - X[j]*term2*yij/dist);
    	        		SYX[i] -= 0.5*(Y[i]*term2*xij/dist - Y[j]*term2*xij/dist);
	    				SYY[i] -= 0.5*(Y[i]*term2*yij/dist - Y[j]*term2*yij/dist);	    			
                    }
				}
				//Ring2 with Ring1
				if((i>=N-NB)&&(j<N-NB))	{
					xij  = X[i] - X[j];
                    yij  = Y[i] - Y[j];
                    vxij = VX[i] - VX[j];
                    vyij = VY[i] - VY[j];
                    dist = xij*xij + yij*yij;
                    dist = sqrt(dist);
                    term1 = (dist - RADIUS)/SIGMA;
                    if(term1 < 1.0) {
                        term12 = term1*term1;
                        term2 = 8.0*EPS/SIGMA*term1*(1-term12)*(1-term12)*(1-term12);
                        VXDOT[i] += term2*xij/dist/MASS[i]; //VXDOT[N-1]
                        VYDOT[i] += term2*yij/dist/MASS[i]; //VXDOT[N-1]

            			SXX[i] -= 0.5*(X[i]*term2*xij/dist - X[j]*term2*xij/dist);
	            		SXY[i] -= 0.5*(X[i]*term2*yij/dist - X[j]*term2*yij/dist);
    	        		SYX[i] -= 0.5*(Y[i]*term2*xij/dist - Y[j]*term2*xij/dist);
	    				SYY[i] -= 0.5*(Y[i]*term2*yij/dist - Y[j]*term2*yij/dist);	    			
                    }
				}				
								
			}
		}
		XDOT[i] = VX[i];
        YDOT[i] = VY[i];
	}
}

void Update(double dt)  {

    int i;

    for(i=0;i<N;i++)   {
        XN[i]   = X[i] + 0.5*dt*XDOT[i];
        YN[i]   = Y[i] + 0.5*dt*YDOT[i];
        VXN[i]  = VX[i] + 0.5*dt*VXDOT[i];
        VYN[i]  = VY[i] + 0.5*dt*VYDOT[i];
    }
}

void callprint(int loop) {

    int i,j;
    double TotE = 0.0, TotKE = 0.0, TotPE = 0.0;
    //For trajectory visualization in LAMMPS
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",loop);
    fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,"%d\n",N);
    fprintf(fp,"ITEM: BOX BOUNDS ss ss ss\n");
    fprintf(fp,"0.00000 1.0\n");
    fprintf(fp,"0.00000 1.0\n");
    fprintf(fp,"0.00000 0.00000\n");
    fprintf(fp,"ITEM: ATOMS id type xs ys zs Pxx Pxy Pyx Pyy\n");
    for(i=0; i<N; i++) {
        E[i] = KE[i] + PE[i];
        TotE += E[i];
        TotKE += KE[i];
        TotPE += PE[i];
        fprintf(fp,"%d %d %lf %lf 0.000 %lf %lf %lf %lf \n", i, 0, X[i], Y[i], SXX[i], SXY[i], SYX[i], SYY[i] );
    }

    printf("Total KE: %lf Total PE: %lf Total Energy is :%lf\n",TotKE,TotPE,TotE);
}

void Energy()   {
    int i,j;
    double xij,yij,dist;
    double term1, term12, term2,term3,term4,term5,fpair, rsq = r*r;

    for(i=0;i<N;i++)  KE[i]    = 0.5*MASS[i]*(VX[i]*VX[i] + VY[i]*VY[i]);
    for(i=0;i<N;i++)  PE[i]    = 0.0;

    for(i=0;i<N;i++)    {
        for(j=0;j<N-NB;j++)    {
            if(i!=j)    {
                //For Ball-Plate Interaction
                if(i>=N-NB)    {
                    xij  = X[i] - X[j];
                    yij  = Y[i] - Y[j];
                    dist = xij*xij + yij*yij;
                    dist = sqrt(dist);
                    term1 = (dist - RADIUS)/SIGMA;
                    if(term1 < 1.0) {
                        term12 = 1.0 - term1*term1;
                        term2 = 0.5*EPS*term12*term12*term12*term12;
                        PE[i] += term2;
                        PE[j] += term2;
                    }
                }
                //For Plate-Plate Interaction
                if(i<N-NB)  {
                    xij  = X[i] - X[j];
                    yij  = Y[i] - Y[j];

                    dist = xij*xij + yij*yij;
                    if(dist<rsq)   {
                        term1       = 2.0 - dist;
                        term2       = pow(term1,2.0*M)*M/(2*M-M);
                        term3       = pow(term1,M)*2.0*M/(2*M-M);
                        fpair       = term2-term3;

                        PE[i]    += 0.5*fpair; //Computing the Force on ith particle due to jth particle along x axis
                    }
                }
            }
        }
    }
    //For Ball-Ball Interaction
    for(i=N-NB;i<N;i++)    {
        for(j=N-NB;j<N;j++)    {
            if(i!=j)    {
                xij  = X[i] - X[j];
                yij  = Y[i] - Y[j];
                dist = xij*xij + yij*yij;
                if(dist<rsq)   {
                    term1       = 2.0 - dist;
                    term2       = pow(term1,2.0*M)*M/(2*M-M);
                    term3       = pow(term1,M)*2.0*M/(2*M-M);
                    fpair       = term2-term3;

                    PE[i]    += 0.5*fpair; //Computing the Force on ith particle due to jth particle along x axis
                }
            }
        }
    }
}

