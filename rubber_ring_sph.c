// Mass has been modified to account for boundaries -- changes in the initialize and the derivative
// Added a function for artificial viscosity
// Value of density changed
// Added a function to test consistency
// Changed mass of particles to ensure consistency
// Changed Rhobar in the derivative function
// Changed Fracture Criteria
// Altered Monaghan's correction in VXBAR as well
// Using B-Spline Kernel

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NX      100
#define NY      100
#define DY      1.0
#define N       5116
#define NB      2558
#define RHO0    1.1547
#define MASSP   1.00 //Comes from consistency
#define MASSB   1.00
#define VRING   0.10//-04.0
#define SIGMA   0.80
#define EPS     100.0
#define RADIUS  6.0
#define MAX     400.0
#define DT      0.01
#define H       3.0
#define M       2.0
#define EMOD    36.950
#define G       13.856
#define B       27.713
#define YIELD   1.360
#define ULTI    4.270
#define NPRINT  10
#define alpha   0.50
#define beta    0.50
#define ART_VISCOSITY 1
#define MONA_CORR     1
#define MONA_CORR2    2
#define JAUMANN       1
#define TENSILE       1
#define GRADCORR      1
#define ri      30.00
#define ro      40.00


double X[N], Y[N], VX[N], VY[N], RHO[N];
double E[N], SXX[N], SXY[N], SYY[N];
double XO[N], YO[N], VXO[N], VYO[N], RHOO[N]; //For Storing Old variables
double EO[N], SXXO[N], SXYO[N], SYYO[N];
double XN[N], YN[N], VXN[N], VYN[N], RHON[N]; //For Storing New variables
double EN[N], SXXN[N], SXYN[N], SYYN[N];
double P[N], KE[N], PE[N], MASS[N];
double XDOT[N], YDOT[N], VXDOT[N], VYDOT[N], RHODOT[N], EDOT[N], SXXDOT[N], SXYDOT[N], SYYDOT[N];
double DVXDX[N], DVXDY[N], DVYDX[N], DVYDY[N], FX[N], FY[N], VXBAR[N], VYBAR[N], DeltaWP[N];
double EINT, WP;
double DAMAGE[N];
double wfgrad[N][5];

FILE *fp;
FILE *fp5;

void initialize();
void derivatives();
void kernel(double wf[], double dist );
void Update(double dt);
void callprint(int loopcounter);
void plasticity();
void viscosity(double visc[], int i, int j, double xij, double yij, double vxij, double vyij, double dist, double wf[], double wf_grad, double wfdx, double wfdy);
void consistency();
void artificial_pressure(double Pi, double Pj, double rho2i, double rho2j, double wf[], double xij, double yij, double dist, int i, int j, double vxij, double vyij, double wf_grad, double wfdx, double wfdy);
void monacorr(int i, int j, double wf[], double vxij, double vyij, double wf_grad, double wfdx, double wfdy);
void basicsph(int i,int j,double wf[], double xij, double yij, double vxij, double vyij, double rho2i, double rho2j, double dist, double sigxxi, double sigyyi, double sigxyi, double sigxxj, double sigyyj, double sigxyj, double drhobar, double wf_grad, double wfdx, double wfdy);

int main()  {
    double t=0.00;
    int i, loopcounter=0;
    fp = fopen("Trajectory.lammpstrj","w");
    fp5 = fopen("Energy.output","w");

    initialize();

    while(t<MAX)    {

        for(i=0;i<N;i++)    {
            XO[i]   = X[i];     YO[i]   = Y[i];
            VXO[i]  = VX[i];   VYO[i]  = VY[i];
            RHOO[i] = RHO[i];   EO[i]   = E[i];
            SXXO[i] = SXX[i];   SXYO[i] = SXY[i];  SYYO[i] = SYY[i];
        }
        derivatives();
         Update(DT);

        for(i=0;i<N;i++)    {
            X[i]   = XN[i];     Y[i]   = YN[i];
            VX[i]  = VXN[i];    VY[i]  = VYN[i];
            RHO[i] = RHON[i];   E[i]   = EN[i];
            SXX[i] = SXXN[i];   SXY[i] = SXYN[i];  SYY[i] = SYYN[i];
        }
        derivatives(); //Rate Computation is with respect to XN
        for(i=0;i<N;i++)    {
            X[i]   = XO[i];     Y[i]   = YO[i];
            VX[i]  = VXO[i];    VY[i]  = VYO[i];
            RHO[i] = RHOO[i];   E[i]   = EO[i];
            SXX[i] = SXXO[i];   SXY[i] = SXYO[i];  SYY[i] = SYYO[i];
        }
        Update(DT); //Update is with respect to XO;

        for(i=0;i<N;i++)   {
            X[i]        = 2.0 * XN[i] - XO[i];
            Y[i]        = 2.0 * YN[i] - YO[i];
            VX[i]       = 2.0 * VXN[i] - VXO[i];
            VY[i]       = 2.0 * VYN[i] - VYO[i];
            SXX[i]      = 2.0 * SXXN[i] - SXXO[i];
            SXY[i]      = 2.0 * SXYN[i] - SXYO[i];
            SYY[i]      = 2.0 * SYYN[i] - SYYO[i];
            RHO[i]      = 2.0 * RHON[i] - RHOO[i];
            E[i]        = 2.0 * EN[i] - EO[i];
        }

        //plasticity();

        if(loopcounter%NPRINT==0) {
            double KE1 = 0.0, IE1 = 0.0, KE2 = 0.0, IE2 = 0.0, WP1 = 0.0;
            printf("Time(ms): %lf ", (double) DT*loopcounter);
            for(i=0; i<N; i++) {
                if(i< N-NB) {
                    KE1 += 0.5 * MASS[i] * ( VX[i]*VX[i] + VY[i]*VY[i] );
                    IE1 += MASS[i]* E[i];
                }
                if(i>= N-NB)    {
                    KE2 += 0.5 * MASS[i] * ( VX[i]*VX[i] + VY[i]*VY[i] );
                    IE2 += MASS[i]* E[i];
                }
                WP1 += DeltaWP[i];
            }
            fprintf(fp5,"%lf %lf %lf %lf %lf %lf %lf \n",(double) DT*loopcounter,KE1,IE1,KE2,IE2,WP1,KE1+IE1+KE2+IE2+EINT);
            callprint(loopcounter);
        }
        t+=DT;
        loopcounter++;
    }
    fclose(fp5);
    return(0);
}


void derivatives()  {
    double xij, yij, vxij, vyij, dist, hsq = H*H;
    double rhobar,drhobar;
    double smxx, smxy, smyx, smyy, rho0rho;
    double wf[2], lambda=B - G, eta=G, trace;
    double sigxxi, sigxyi, sigyyi, sigxxj, sigxyj, sigyyj;
    double rho2i, rho2j;
    double term1, term12, term2;
    double visc[1], rxydot[N];
    int i,j;
    double tensterm = 0.0, tensgamma = 0.3, distinit, xinitij, yinitij, wfinit[2];
    double mata, matb, matc, matd, wfdx, wfdy, wf_grad;

    for(i=0;i<N;i++)    {
        RHODOT[i]   = 0.0;
        VXDOT[i]    = VYDOT[i] = 0.0;
        DVXDX[i]    = DVXDY[i] = DVYDX[i] = DVYDY[i] = 0.0;
        VXBAR[i]    = 0.0;
        VYBAR[i]    = 0.0;
        FX[i]       = FY[i] = 0.0;
        EDOT[i]     = 0.0;
        rho0rho     = RHO0/RHO[i];
        P[i]        = -6.0*M*RHO0*( pow(2.0-rho0rho,2*M-1) - pow(2.0-rho0rho,M-1) );
        rxydot[i]   = 0.0;
        EINT        = 0.0;
    }
    
    //Kernel Computation for Gradient Correction
    for(i=0;i<N;i++)	{
    	for(j=0;j<5;j++)	{
    		wfgrad[i][j] = 0.0;
    	}
    }
    if(GRADCORR == 0)	{
    	for(i=0;i<N;i++)	{
    		wfgrad[i][0] = 1.0;
    		wfgrad[i][1] = 1.0;
    		wfgrad[i][2] = 0.0;
    		wfgrad[i][3] = 0.0;
    		wfgrad[i][4] = 1.0;    		    		    		    		
    	}
    }
    if(GRADCORR == 1)	{
    	for(i=0;i<N-NB;i++)	{
    		for(j=0;j<N-NB;j++)	{
				if(i!=j)	{
					xij  = X[i] - X[j];
                	yij  = Y[i] - Y[j];
					dist = xij*xij + yij*yij;
                	if(dist<=hsq)   {
                		dist = sqrt(dist);
                    	kernel(wf, dist);
                    	wfgrad[i][0] += wf[0]*MASS[j]/RHO[j];
                    	wfgrad[i][1] += -1.0*xij*wf[1]*xij/dist*MASS[j]/RHO[j];
                    	wfgrad[i][2] += -1.0*yij*wf[1]*xij/dist*MASS[j]/RHO[j];
                    	wfgrad[i][3] += -1.0*xij*wf[1]*yij/dist*MASS[j]/RHO[j];
                    	wfgrad[i][4] += -1.0*yij*wf[1]*yij/dist*MASS[j]/RHO[j];	
                	}
				}
    		}
    	}
    }
    if(GRADCORR == 1)	{
    	for(i=N-NB;i<N;i++)	{
    		for(j=N-NB;j<N;j++)	{
				if(i!=j)	{
					xij  = X[i] - X[j];
                	yij  = Y[i] - Y[j];
					dist = xij*xij + yij*yij;
                	if(dist<=hsq)   {
                		dist = sqrt(dist);
                    	kernel(wf, dist);
                    	wfgrad[i][0] += wf[0]*MASS[j]/RHO[j];
                    	wfgrad[i][1] += -1.0*xij*wf[1]*xij/dist*MASS[j]/RHO[j];
                    	wfgrad[i][2] += -1.0*yij*wf[1]*xij/dist*MASS[j]/RHO[j];
                    	wfgrad[i][3] += -1.0*xij*wf[1]*yij/dist*MASS[j]/RHO[j];
                    	wfgrad[i][4] += -1.0*yij*wf[1]*yij/dist*MASS[j]/RHO[j];	
                	}
				}
    		}
    	}
    }
    
//Actual Computation
    for(i=0;i<N;i++)    {
        sigxxi = SXX[i] + P[i];
        sigyyi = SYY[i] + P[i];
        sigxyi = SXY[i];
        rho2i  = 1.0/(RHO[i]*RHO[i]);
        for(j=0;j<N-NB;j++)    { //Excluding the ball particle
            if(i!=j)    {
                //For Ball-Plate Interaction
                if(i>=N-NB)    {
                    xij  = X[i] - X[j];
                    yij  = Y[i] - Y[j];
                    dist = xij*xij + yij*yij;
                    dist = sqrt(dist);
                    term1 = (dist - RADIUS)/SIGMA;
                    if(term1 < 1.0) {
                        term12 = term1*term1;
                        term2 = 8.0*EPS/SIGMA*term1*(1-term12)*(1-term12)*(1-term12);
                        VXDOT[i] += term2*xij/dist/MASS[i];
                        VYDOT[i] += term2*yij/dist/MASS[i];
                        VXDOT[j] -= term2*xij/dist/MASS[j];
                        VYDOT[j] -= term2*yij/dist/MASS[j];
                        EINT     += EPS*(1.0-term12)*(1.0-term12)*(1.0-term12)*(1.0-term12);
                    }
                }
                //For Plate-Plate Interaction
                if(i<N-NB)  {
                    xij  = X[i] - X[j];
                    yij  = Y[i] - Y[j];
                    vxij = VX[i] - VX[j];
                    vyij = VY[i] - VY[j];

                    dist = xij*xij + yij*yij;
                    if(dist<=hsq)   {
                        dist        = sqrt(dist);
                        mata = wfgrad[i][1];
                        matb = wfgrad[i][2];
                        matc = wfgrad[i][3];
                        matd = wfgrad[i][4];
                        
                        kernel(wf, dist);
                        wf_grad 	= wf[0]*MASS[j]/RHO[j]/wfgrad[i][0];
                        //printf("%d %d %lf \n", i, j, 1.0/(mata*matd-matb*matc));
                    	wfdx		= 1.0/(mata*matd-matb*matc)*(wfgrad[i][4]*wf[1]*xij/dist - wfgrad[i][2]*wf[1]*yij/dist);
                    	wfdy		= 1.0/(mata*matd-matb*matc)*(wfgrad[i][1]*wf[1]*yij/dist - wfgrad[i][3]*wf[1]*xij/dist);
                    	//wfdx		= 0.50/(mata*matd-matb*matc)*(wfgrad[i][4]*wf[1]*xij/dist - wfgrad[i][2]*wf[1]*yij/dist);
                    	//wfdy		= 0.50/(mata*matd-matb*matc)*(wfgrad[i][1]*wf[1]*yij/dist - wfgrad[i][3]*wf[1]*xij/dist);

                        //mata = wfgrad[j][1];
                        //matb = wfgrad[j][2];
                        //matc = wfgrad[j][3];
                        //matd = wfgrad[j][4];
                    	//wfdx		-= 0.50/(mata*matd-matb*matc)*(-1.0*wfgrad[j][4]*wf[1]*xij/dist + wfgrad[j][2]*wf[1]*yij/dist);
                    	//wfdy		-= 0.50/(mata*matd-matb*matc)*(-1.0*wfgrad[j][1]*wf[1]*yij/dist + wfgrad[j][3]*wf[1]*xij/dist);
                        
                        rhobar      = RHO[j];
                        drhobar     = dist*rhobar;
                    	//RHODOT[i]   += wf[1]*MASS[j]*(vxij*xij + vyij*yij)/dist;
                    	RHODOT[i]   += MASS[j]*(vxij*wfdx + vyij*wfdy);
                    
                        sigxxj      = SXX[j] + P[j];
                        sigyyj      = SYY[j] + P[j];
                        sigxyj      = SXY[j];
                        rho2j       = 1.0/(RHO[j]*RHO[j]);

                        if(MONA_CORR == 1)  monacorr(i,j,wf,vxij,vyij,wf_grad,wfdx,wfdy);
                        if(JAUMANN == 1) rxydot[i] += -0.5 * MASS[j]/RHO[j]*(vxij*wfdy - vyij*wfdx);
                        
                        basicsph(i,j,wf,xij,yij,vxij,vyij,rho2i,rho2j,dist,sigxxi,sigyyi,sigxyi,sigxxj,sigyyj,sigxyj,drhobar,wf_grad,wfdx,wfdy);
                        
                        if(ART_VISCOSITY == 1)  viscosity(visc,i,j,xij,yij,vxij,vyij,dist,wf,wf_grad,wfdx,wfdy);
                        if(TENSILE == 1) artificial_pressure(P[i],P[j],rho2i,rho2j,wf, xij, yij, dist, i, j, vxij, vyij,wf_grad,wfdx,wfdy);
                    }
                }
            }
        }
        trace       = -1.0/3.0*(DVXDX[i] + DVYDY[i]);
        SXXDOT[i]   = 2.0*G*(DVXDX[i] + trace) + 2.0*SXY[i]*rxydot[i];
        SYYDOT[i]   = 2.0*G*(DVYDY[i] + trace) - 2.0*SXY[i]*rxydot[i];
        SXYDOT[i]   = 1.0*G*(DVXDY[i] + DVYDX[i]) - rxydot[i]*(SXX[i] - SYY[i]);
        XDOT[i] = VX[i] + VXBAR[i];
        YDOT[i] = VY[i] + VYBAR[i];
        EDOT[i] = -0.5*EDOT[i];
    }

    //For Ball-Ball Interaction
    for(i=N-NB;i<N;i++)    {
        sigxxi = SXX[i] + P[i];
        sigyyi = SYY[i] + P[i];
        sigxyi = SXY[i];
        rho2i  = 1.0/(RHO[i]*RHO[i]);
        for(j=N-NB;j<N;j++)    {
            if(i!=j)    {
                xij  = X[i] - X[j];
                yij  = Y[i] - Y[j];
                vxij = VX[i] - VX[j];
                vyij = VY[i] - VY[j];
                dist = xij*xij + yij*yij;
                if(dist<=hsq)   {
                    dist        = sqrt(dist);
                        mata = wfgrad[i][1];
                        matb = wfgrad[i][2];
                        matc = wfgrad[i][3];
                        matd = wfgrad[i][4];                    
                    kernel(wf, dist);
                    wf_grad 	= wf[0]*MASS[j]/RHO[j]/wfgrad[i][0];
                    wfdx		= 1.0/(mata*matd-matb*matc)*(wfgrad[i][4]*wf[1]*xij/dist - wfgrad[i][2]*wf[1]*yij/dist);
                    wfdy		= 1.0/(mata*matd-matb*matc)*(wfgrad[i][1]*wf[1]*yij/dist - wfgrad[i][3]*wf[1]*xij/dist);
                    //wfdx		= 0.50/(mata*matd-matb*matc)*(wfgrad[i][4]*wf[1]*xij/dist - wfgrad[i][2]*wf[1]*yij/dist);
                    //wfdy		= 0.50/(mata*matd-matb*matc)*(wfgrad[i][1]*wf[1]*yij/dist - wfgrad[i][3]*wf[1]*xij/dist);

                    //mata = wfgrad[j][1];
                    //matb = wfgrad[j][2];
                    //matc = wfgrad[j][3];
                    //matd = wfgrad[j][4];
                    //wfdx		-= 0.50/(mata*matd-matb*matc)*(-1.0*wfgrad[j][4]*wf[1]*xij/dist + wfgrad[j][2]*wf[1]*yij/dist);
                    //wfdy		-= 0.50/(mata*matd-matb*matc)*(-1.0*wfgrad[j][1]*wf[1]*yij/dist + wfgrad[j][3]*wf[1]*xij/dist);
                    	                    
                    rhobar      = RHO[j];
                    drhobar     = dist*rhobar;
                    //RHODOT[i]   += wf[1]*MASS[j]*(vxij*xij + vyij*yij)/dist;
                    RHODOT[i]   += MASS[j]*(vxij*wfdx + vyij*wfdy);
                    
                    sigxxj      = SXX[j] + P[j];
                    sigyyj      = SYY[j] + P[j];
                    sigxyj      = SXY[j];
                    rho2j       = 1.0/(RHO[j]*RHO[j]);

                    if(MONA_CORR == 1)  monacorr(i,j,wf,vxij,vyij,wf_grad,wfdx,wfdy);
                    if(JAUMANN == 1) rxydot[i] += -0.5 * MASS[j]/RHO[j]*(vxij*wfdy - vyij*wfdx);
                    
                    basicsph(i,j,wf,xij,yij,vxij,vyij,rho2i,rho2j,dist,sigxxi,sigyyi,sigxyi,sigxxj,sigyyj,sigxyj,drhobar,wf_grad,wfdx,wfdy);
                    
                    if(ART_VISCOSITY == 1)  viscosity(visc,i,j,xij,yij,vxij,vyij,dist,wf,wf_grad,wfdx,wfdy);
                    if(TENSILE == 1)  artificial_pressure(P[i],P[j],rho2i,rho2j,wf, xij, yij, dist, i, j, vxij, vyij,wf_grad,wfdx,wfdy);
                }
            }
        }
        trace       = -1.0/3.0*(DVXDX[i] + DVYDY[i]);
        SXXDOT[i]   = 2.0*G*(DVXDX[i] + trace) + 2.0*SXY[i]*rxydot[i];
        SYYDOT[i]   = 2.0*G*(DVYDY[i] + trace) - 2.0*SXY[i]*rxydot[i];
        SXYDOT[i]   = 1.0*G*(DVXDY[i] + DVYDX[i]) - rxydot[i]*(SXX[i] - SYY[i]);
        XDOT[i] = VX[i] + VXBAR[i];
        YDOT[i] = VY[i] + VYBAR[i];
        EDOT[i] = -0.5*EDOT[i];
    }
    
    for(i=0;i<N;i++)	{
    	if(MONA_CORR2 == 1)     {
    		VX[i]   = VX[i] + VXBAR[i];
            VY[i]   = VY[i] + VYBAR[i];
    	}
    }
}

void basicsph(int i,int j,double wf[], double xij, double yij, double vxij, double vyij, double rho2i, double rho2j, double dist, double sigxxi, double sigyyi, double sigxyi, double sigxxj, double sigyyj, double sigxyj, double drhobar, double wf_grad, double wfdx, double wfdy)   {
    double smxx, smyx, smxy, smyy;
	double rhobar = drhobar/dist;
	
    smxx        = -wfdx*vxij/rhobar;
    smxy		= -wfdy*vxij/rhobar;
    smyx        = -wfdx*vyij/rhobar;
    smyy        = -wfdy*vyij/rhobar;
    
//    smxx        = -wf[1]*(vxij * xij)/drhobar;
//    smxy        = -wf[1]*(vxij * yij)/drhobar;
//    smyx        = -wf[1]*(vyij * xij)/drhobar;
//    smyy        = -wf[1]*(vyij * yij)/drhobar;
    DVXDX[i]    += MASS[j]*smxx;
    DVXDY[i]    += MASS[j]*smxy;
    DVYDX[i]    += MASS[j]*smyx;
    DVYDY[i]    += MASS[j]*smyy;

    VXDOT[i]    += MASS[j]*(sigxxi*rho2i + sigxxj*rho2j)*wfdx + MASS[j]*(sigxyi*rho2i + sigxyj*rho2j)*wfdy;
    VYDOT[i]    += MASS[j]*(sigyyi*rho2i + sigyyj*rho2j)*wfdy + MASS[j]*(sigxyi*rho2i + sigxyj*rho2j)*wfdx;
    
    EDOT[i]     += MASS[j]*(vxij*(sigxxi*rho2i + sigxxj*rho2j)*wfdx + vxij*(sigxyi*rho2i + sigxyj*rho2j)*wfdy);
    EDOT[i]     += MASS[j]*(vyij*(sigyyi*rho2i + sigyyj*rho2j)*wfdy + vyij*(sigxyi*rho2i + sigxyj*rho2j)*wfdx);
}

void monacorr(int i, int j, double wf[], double vxij, double vyij, double wf_grad, double wfdx, double wfdy)  {
    VXBAR[i]    -= 0.50*MASS[j]*wf[0]*vxij/(RHO[i] + RHO[j]);
    VYBAR[i]    -= 0.50*MASS[j]*wf[0]*vyij/(RHO[i] + RHO[j]);
}

void artificial_pressure(double Pi, double Pj, double rho2i, double rho2j, double wf[], double xij, double yij, double dist, int i, int j, double vxij, double vyij, double wf_grad, double wfdx, double wfdy) {
    double tensterm1=0.0, tensterm2=0.0, tensterm, wfinit[2];
    double tensgamma = 0.50;
    double temp;
    kernel(wfinit, DY);

    temp = wf[0]/wfinit[0];
    temp = temp*temp;
    temp = temp*temp;
    //temp = temp*temp;
    
    if(Pi > 0.0) tensterm1 = Pi*rho2i;
    if(Pj > 0.0) tensterm2 = Pj*rho2j;
    tensterm = tensgamma*(tensterm1 + tensterm2)*temp;
    VXDOT[i]    += -1.0*MASS[j]*tensterm*wfdx;
    VYDOT[i]    += -1.0*MASS[j]*tensterm*wfdy;
    EDOT[i]     += -1.0*MASS[j]*vxij*tensterm*wfdx;
    EDOT[i]     += -1.0*MASS[j]*vyij*tensterm*wfdy;
}

void kernel(double wf[], double dist )   {	
    double h = H;	
    double q = dist/h;	
    double par = 5.0/(3.14159265359*h*h);	
    wf[0] = par*(1.0+3.0*q)*(1.0-q)*(1.0-q)*(1.0-q);	
    wf[1] = -par*12/h*q*(1-q)*(1-q);	
}

void Update(double dt)  {

    int i;

    for(i=0;i<N;i++)   {
        RHON[i] = RHO[i] + 0.5*dt*RHODOT[i];
        EN[i]   = E[i] + 0.5*dt*EDOT[i];
        XN[i]   = X[i] + 0.5*dt*XDOT[i];
        YN[i]   = Y[i] + 0.5*dt*YDOT[i];
        VXN[i]  = VX[i] + 0.5*dt*VXDOT[i];
        VYN[i]  = VY[i] + 0.5*dt*VYDOT[i];
        SXXN[i] = SXX[i] + 0.5*dt*SXXDOT[i];
        SXYN[i] = SXY[i] + 0.5*dt*SXYDOT[i];
        SYYN[i] = SYY[i] + 0.5*dt*SYYDOT[i];
    }
}

void plasticity()   {
    int i;
    double shear, sigxx, sigxy, sigyy, factor, sum, diff;

    for(i=0;i<N;i++)  {
        sigxx = SXX[i] + P[i];
        sigyy = SYY[i] + P[i];
        sigxy = SXY[i];
        shear = sqrt(sigxy*sigxy + 0.25*(sigxx - sigyy)*(sigxx - sigyy) );

        if(shear > YIELD)   {
            factor = YIELD/shear;
            sum = sigxx + sigyy;
            diff = sigxx - sigyy;
            sigxy  = factor*sigxy;
            sigxx = 0.5*sum + 0.5*factor*diff;
            sigyy = 0.5*sum - 0.5*factor*diff;
            DeltaWP[i] = (1.0 - factor)*factor*(2.0*SXY[i]*SXY[i] + SXX[i]*SXX[i] + SYY[i]*SYY[i])*MASS[i]/3.00/G/RHO[i];
        }

        if(0.5*(sigxx + sigyy) > ULTI)  {
            sigxy = 0.0;
            sigxx = P[i];
            sigyy = P[i];
            RHO[i] = RHO0;
            DAMAGE[i] = 0.0;
        }
        SXX[i] = sigxx - P[i];
        SYY[i] = sigyy - P[i];
        SXY[i] = sigxy;
    }
}

void callprint(int loop) {

    int i;
    double KE1 = 0.0,IE1 = 0.0, KE2 = 0.0, IE2 = 0.0;
    //For trajectory visualization in LAMMPS
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",loop);
    fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,"%d\n",N);
    fprintf(fp,"ITEM: BOX BOUNDS ss ss ss\n");
    fprintf(fp,"0.00000 1.0\n");
    fprintf(fp,"0.00000 1.0\n");
    fprintf(fp,"0.00000 0.00000\n");
    fprintf(fp,"ITEM: ATOMS id type xs ys zs sxx sxy syy\n");
    for(i=0; i<N; i++) {
        fprintf(fp,"%d %d %lf %lf 0.000 %lf %lf %lf \n", i, 0, X[i], Y[i],(SXX[i] +P[i])*MASS[i]/RHO[i],(SXY[i] +P[i])*MASS[i]/RHO[i],(SYY[i] +P[i])*MASS[i]/RHO[i] );
        //Total Energy calculation
        if(i< N-NB) {
            KE1 += 0.5 * MASS[i] * ( VX[i]*VX[i] + VY[i]*VY[i] );
            IE1 += MASS[i]* E[i];
        }
        if(i>= N-NB)    {
            KE2 += 0.5 * MASS[i] * ( VX[i]*VX[i] + VY[i]*VY[i] );
            IE2 += MASS[i]* E[i];
        }
        WP += DeltaWP[i];
    }
    printf("KE1:%lf IE1:%lf KE2:%lf IE2:%lf WP:%lf Total Energy:%lf Interaction Energy: %lf Total Energy: %lf \n",KE1,IE1,KE2,IE2,WP,KE1+IE1+KE2+IE2, EINT, KE1+IE1+KE2+IE2+EINT);

    //For Velocity verification in LAMMPS
    //for(i=0; i<Nt; i++) fprintf(output,"%d %Lf %Lf %Lf %Lf %Lf %Lf, %Lf, %Lf, %Lf\n", i, PosX[i], PosY[i], VelX[i], VelY[i], Pressure[i], SXX[i], SXY[i], SYX[i], SYY[i])
}

void viscosity(double visc[], int i, int j, double xij, double yij, double vxij, double vyij, double dist, double wf[], double wf_grad, double wfdx, double wfdy)    {

    double delVdelR, muij, ci = sqrt(EMOD/RHO[i]), cj = sqrt(EMOD/RHO[j]);

    delVdelR    = xij*vxij + yij*vyij;
    muij        = (H*delVdelR)/(dist*dist + 0.01*H*H);
    visc[0]     = 0.0; //Force part
    if(delVdelR < 0.0) visc[0] = ( -0.5*(ci+cj)*muij*alpha + beta*muij*muij )/ ( 0.5*(RHO[i] + RHO[j]) );
    VXDOT[i]    += -1.0*MASS[j]*visc[0]*wfdx;
    VYDOT[i]    += -1.0*MASS[j]*visc[0]*wfdy;
    EDOT[i]     += -1.0*MASS[j]*vxij*visc[0]*wfdx;
    EDOT[i]     += -1.0*MASS[j]*vyij*visc[0]*wfdy;
}

void consistency()  {
    int i,j;
    double xij,yij,dist,wf[2],sumKernel[N];
    double hsq = H*H;
    FILE *fp3;

    fp3 = fopen("Consistency_Check.dat","w");

    for(i=0;i<N;i++)    {
        sumKernel[i] = 0.0;
        for(j=0;j<N;j++)    {
            //if(i!=j)    {
            xij  = X[i] - X[j];
            yij  = Y[i] - Y[j];
            dist = xij*xij + yij*yij;
            if(dist<=hsq)   {
                dist        = sqrt(dist);
                kernel(wf, dist);
                sumKernel[i] += wf[0]*MASS[j]/RHO[j];
            }
        }
        //fprintf(fp3,"%d KernelSum: %lf \n",i,sumKernel[i]);
    }
    for(i=0;i<N;i++)    {
        MASS[i] = MASS[i]/sumKernel[i];
        fprintf(fp3,"%d %lf %lf %lf %lf \n",i,X[i],Y[i],MASS[i],sumKernel[i]);
    }
    fclose(fp3);
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
RHO[count]  = RHO0;
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
RHO[count]  = RHO0;

               //fprintf(fp2,"%lf %lf \n",X[count],Y[count]);
               count++;
            }
        }
    }
    //Read Mass from SPH data file constency_check.dat
    //for(i=0;i<N;i++)  fscanf(fp4, "%d %lf %lf %lf %lf", &tempint, &temp, &temp, &MASS[i], &temp);
    tot_mass = tot_momentum = 0.0;
    for(i=0;i<N;i++) tot_mass += MASS[i];
    for(i=0;i<N;i++) tot_momentum += 0.50*MASS[i]*VX[i]*VX[i];
    printf("Total Count is: %d \n",count);
    printf("Total Momentum is: %lf \n",tot_momentum);
    printf("Total Mass is: %lf \n",tot_mass);

    //Initialize Velocity of the Particles so that the net momentum of the system is zero
     //for(i=0;i<N-NB;i++) VY[i] = (double) -1.0*tot_momentum/(tot_mass);
    fclose(fp2);
}
