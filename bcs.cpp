#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

//////////////////////Boundary conditions on the top, bottom and side walls 

void bcs(vertex * node, fval * fvar)
{	
	int ind; 

	tdma2x(fvar,0);  //d2u/dx2  //Computing the second derivatives needed for pressure BCs
	tdma2y(fvar,1);  //d2v/dy2 

	//Top wall 

	for(int i=1;i<nx;i++)
	{
		ind = i + ny*str_x; 					
		
		fvar[ind].u[0] = 1.0;     //The velocity of the moving top plate 
		fvar[ind].u[1] = 0.0; 	
		fvar[ind].u[2] =  -(1./B1)*dy*(1./Re)*fvar[ind].uyy[1] - B2*fvar[ind-str_x].u[2] - B3*fvar[ind-2*str_x].u[2] - B4*fvar[ind-3*str_x].u[2];
	}

	//Right wall 

	for(int j=0;j<=ny;j++)
	{
		ind = nx + j*str_x; 			

		fvar[ind].u[0] = 0.0; 
		fvar[ind].u[1] = 0.0; 
		fvar[ind].u[2] = -(1./B1)*dx*(1./Re)*fvar[ind].uxx[0] - B2*fvar[ind-1].u[2] - B3*fvar[ind-2].u[2] - B4*fvar[ind-3].u[2]; 
	}

	//Bottom wall 

	for(int i=1;i<=nx;i++)
	{			
		fvar[i].u[0] = 0.0;
		fvar[i].u[1] = 0.0;
		fvar[i].u[2]=  (1./B1)*dy*(1./Re)*fvar[i].uyy[1] - (B2*fvar[i+str_x].u[2] + B3*fvar[i+2*str_x].u[2] + B4*fvar[i+3*str_x].u[2]);			
	}

	//Left wall 

	for(int j=0;j<=ny;j++)
	{
	 	ind = j*str_x;  	 		
	
		fvar[ind].u[0] = 0.0;	
		fvar[ind].u[1] = 0.0;	
		fvar[ind].u[2] = (1./B1)*dx*(1./Re)*fvar[ind].uxx[0] - ( B2*fvar[ind+1].u[2] + B3*fvar[ind+2].u[2] + B4*fvar[ind+3].u[2]) ;
	}

}
