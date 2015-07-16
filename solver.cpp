#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

///////////////////Solver function////////////////////

void solver(vertex * node, fval * fvar, mg_grid * level, pbcs& pbc, int restart_count)
{		
	double max_div; 
	ofstream fp_res; //Filestream to write the residuals as per time for the 3 flow variables
	
	fp_res.open("global_residuals.dat"); 
	
	cout<<"The time step is "<<dt<<"\n"; 
	cout<<"The no.of timesteps is "<<ite<<"\n"; 
	cout<<"The spacing is dx = "<<dx<< " dy= "<<dy<<"\n"; 

	for(int t=restart_count+1;t<=ite;t++)
	{
		explicit_solver(node,fvar,level,pbc);
		max_div = div_calc(node,fvar);
		  
		cout<<"----------- "<<"ite "<<t<<" "<<max_div<<"------------"<<"\n"; 	
		global_residuals(fvar,fp_res,t);
		
		if(t%file_freq==0) write_to_file(node,fvar,t); 		
		if(t%r_file_freq==0) write_restart(fvar,t);
	}
	
	fp_res.close(); 
}

/**********************************************************************************************/

void explicit_solver(vertex * node, fval * fvar, mg_grid * level, pbcs& pbc)
{ 
	int ind; 
			
	for(int i=0;i<tot_p;i++)
	{
		fvar[i].u0[0] = fvar[i].u[0];                                    //U1  stage-1 mod-RK scheme 
		fvar[i].u0[1] = fvar[i].u[1]; 					 //V1 
		fvar[i].u0[2] = fvar[i].u[2]; 					 //Pressure for the time residual
	}

	bcs(node,fvar);
	mg_clear_levels(level); 
	poisson_source(fvar,level,4.0); 
	/*compute_pbc(fvar,pbc,1); 		
	mg_bicgstab_corn(fvar, level, 0, 15000, 1, 0, pbc);
	mg_final(level,fvar,0);*/
        mg_poisson_solver(level,fvar,pbc);
	tdma1x(fvar,2); 
	tdma1y(fvar,2); 

	for(int j=1;j<ny;j++)
	{
		for(int i=1;i<nx;i++) 
		{
			ind = i + j*str_x; 
	
			fvar[ind].u[0] = fvar[ind].u0[0] + (dt/4.)*( fvar[ind].u[3] - fvar[ind].ux[2]);     //U2 stage-2 mod-RK scheme 
			fvar[ind].u[1] = fvar[ind].u0[1] + (dt/4.)*( fvar[ind].u[4] - fvar[ind].uy[2]);     //V2 
		}
	}												

	bcs(node,fvar);
	mg_clear_levels(level);
	poisson_source(fvar,level,3.0);
	/*compute_pbc(fvar,pbc,1);
	mg_bicgstab_corn(fvar, level, 0, 15000, 1, 0, pbc);
	mg_final(level,fvar,0);*/
        mg_poisson_solver(level,fvar,pbc);
	tdma1x(fvar,2); 
	tdma1y(fvar,2); 

	for(int j=1;j<ny;j++)
	{
		for(int i=1;i<nx;i++) 
		{
			ind = i + j*str_x; 
	
			fvar[ind].u[0] = fvar[ind].u0[0] + (dt/3.)*( fvar[ind].u[3] - fvar[ind].ux[2]);    //U3 
			fvar[ind].u[1] = fvar[ind].u0[1] + (dt/3.)*( fvar[ind].u[4] - fvar[ind].uy[2]);	   //V3 
		}
	}

	bcs(node,fvar);
	mg_clear_levels(level);	
	poisson_source(fvar,level,2.0);
	/*compute_pbc(fvar,pbc,1);
	mg_bicgstab_corn(fvar, level, 0, 15000, 1, 0, pbc);
	mg_final(level,fvar,0);*/
        mg_poisson_solver(level,fvar,pbc);
	tdma1x(fvar,2); 
	tdma1y(fvar,2); 

	for(int j=1;j<ny;j++)
	{
		for(int i=1;i<nx;i++) 
		{
			ind = i + j*str_x; 
	
			fvar[ind].u[0] = fvar[ind].u0[0] + (dt/2.)*( fvar[ind].u[3] - fvar[ind].ux[2]);    //U4 
			fvar[ind].u[1] = fvar[ind].u0[1] + (dt/2.)*( fvar[ind].u[4] - fvar[ind].uy[2]);	   //V4 
		}
	}

	bcs(node,fvar);	
	mg_clear_levels(level	);
	poisson_source(fvar,level,1.0);
	/*compute_pbc(fvar,pbc,1);
	mg_bicgstab_corn(fvar, level, 0, 15000, 1, 0, pbc);
	mg_final(level,fvar,0);*/
        mg_poisson_solver(level,fvar,pbc);
	tdma1x(fvar,2); 
	tdma1y(fvar,2); 

	for(int j=1;j<ny;j++)
	{
		for(int i=1;i<nx;i++)
		{
			ind = i + j*str_x; 
	
			fvar[ind].u[0] = fvar[ind].u0[0] + dt*( fvar[ind].u[3] - fvar[ind].ux[2]);    	//U(n+1)
			fvar[ind].u[1] = fvar[ind].u0[1] + dt*( fvar[ind].u[4] - fvar[ind].uy[2]);    	//V(n+1) 
		}
	}	
}
/*************************************************************************************/
