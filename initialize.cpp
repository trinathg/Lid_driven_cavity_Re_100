#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

void initialize(vertex * node, fval * fvar)
{  
	int ind; 

	for(int j=0;j<=ny;j++)
	{
		for(int i=0;i<=nx;i++) 
		{   
			ind = i + j*str_x;   			
    						    
			fvar[ind].F = 0.0;  	    //The RHS of the pressure Poisson equation 			 			
			fvar[ind].res = 0.0; 
			fvar[ind].err = 0.0; 
			fvar[ind].rhs = 0.0; 
			fvar[ind].div = 0.0; 
			
			for(int var=0;var<NUM_VAR;var++)
			{
			  fvar[ind].u[var] = 0.0; 
			  fvar[ind].ux[var] = 0.0; 
			  fvar[ind].uy[var] = 0.0; 
			  fvar[ind].uxx[var] = 0.0; 
			  fvar[ind].uyy[var] = 0.0; 
			  fvar[ind].u0[var] = 0.0; 
			}
			
		} 
	}
	
}

/**********************************************************************/
//Initializing the values of the RHS, i.e, F and others for various levels of multigrid 
void mg_initialize(mg_grid * level)
{
	int ind,lev;
	int sp; 

	for(lev=0;lev<mg_levels;lev++)
	{	
		sp = pow(2,lev); 

		for(int j=0;j<=ny;j=j+sp)
		{
	  		for(int i=0;i<=nx;i=i+sp)
	  		{
			  	ind = i+j*str_x;				  			
	       			    			
	       			level[lev].phi_s[ind] = 0.0; 			 					
				level[lev].F[ind] = 0.0; 									
	  			level[lev].rhs[ind] =0.0;
	  			level[lev].res[ind] = 0.0;
	  			level[lev].cor_rhs[ind] = 0.0;
	  			level[lev].coeff[ind]=0.0; 
	  		}	
		}
	} 
}
/**********************************************************************/
