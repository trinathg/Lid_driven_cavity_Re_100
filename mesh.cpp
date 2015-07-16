#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

void make_mesh()
{

	int ind;  // The global index
	double start_time, end_time;
	int precond_count=0, restart_count; 

	cout<<"The no. of MG levels "<<mg_levels<<"\n";

	/*vector<vertex> node (tot_p);
	vector<fval> fvar (tot_p); 
	vector<mg_grid> level(mg_levels);*/
	
	vertex * node = new vertex[tot_p]; 
	fval * fvar = new fval[tot_p]; 
	mg_grid * level = new mg_grid[mg_levels];
	
	pbcs pbc; 

	for(int j=0;j<=ny;j++) 		//Grid generation 
	{
		for(int i=0;i<=nx;i++)
		{
			ind = i + j*str_x; 
			node[ind].x[0] = i*dx;
			node[ind].x[1] = j*dy;   
		}
	} 
	
	omp_set_num_threads(num_threads); 
  		
	if(restart==0)
        {
	  initialize(node,fvar);
	  mg_initialize(level);			
          restart_count=0;
        }
        else
        {
          read_restart(fvar,restart_count);
          mg_restart(fvar,level);
        }
        
        /***Standalone Direct solver***/
         
	/*int lev=0; 
	int nx_level = nx/pow(2,lev)-1; 
	int ny_level = ny/pow(2,lev)-1; 
	
	mg_coeff();  	      			//Coefficients for smoothing operator at all levels
	mg_read_coeff(level); 			//Reading the null space vectors for all grid sizes
	compute_pbc(pbc,1);
	mg_evaluate_rhs(level,lev,pbc); 
	mg_compute_rhs_tot(level,lev); 	
	double rhs_sum = mg_eval_rhs_sum(level,lev); 	
	level[lev].point_correc = rhs_sum/(nx_level*ny_level); 	  
	
	cout<<"The RHS sum = "<<rhs_sum<<"\n";
	cout<<"The point correction = "<<level[lev].point_correc<<"\n";
	
	//arma_direct_solve(level,lev);
	arma_direct_trunc(level,lev);
	mg_bcs_neu(level,0,pbc); 
	mg_final(level,fvar,0);	
       	write_to_file(node,fvar);       	
       	cout<<"The error norm: "<<error_calc(fvar,0)<<"\n";
       	
  	write_restart(fvar,0);*/
	
	/***Standalone CG Solver****/ 
	
	/*mg_coeff();  	      // Coefficients for smoothing operator at all levels
	mg_read_coeff(level); // Reading the null space vectors for all grid sizes
	compute_pbc(pbc,1);	
	mg_interpolate(level,3,0); 
	mg_conjugate_gradient(level,0, 100000, 0, 0, pbc); 
	mg_bcs_neu(level,0,pbc);
	mg_final(level,fvar,0); 
	write_to_file(node,fvar);
	write_restart(fvar,0);

	cout<<"The error norm: "<<error_calc(fvar,0)<<"\n";*/
	
	
	/***Stand alone interpolation + level:2 testing + arma_direct***/ 
	
	/*mg_coeff();  	      // Coefficients for smoothing operator at all levels
	mg_read_coeff(level); // Reading the null space vectors for all grid sizes			
	compute_pbc(pbc,1);	
	mg_interpolate(level,4);
	mg_interpolate(level,3);  
	mg_conjugate_gradient(level,2, 4000, 0, 2, pbc); 
	mg_bcs_neu(level,0,pbc);
	mg_final(level,fvar,0); 
	write_to_file(node,fvar);
	write_restart(fvar,0);

	cout<<"The error norm: "<<error_calc(fvar)<<"\n";*/
	
	/***V-cycle and FMG solvers*****/

	/*mg_poisson_solver(level,fvar,pbc);  
	write_to_file(node,fvar);
	write_restart(fvar,0);

        cout<<"The error norm: "<<error_calc(fvar,0)<<"\n";*/

	
	/***Stand alone Bi-CGSTAB Solver***/ 
	
	/*mg_coeff();
	mg_read_coeff(level); // Reading the null space vectors for all grid sizes
	compute_pbc(pbc,1); 
	mg_bicgstab_corn(level,0, 1000, 0, 0, pbc); 
	mg_bcs_neu(level,0,pbc);
	mg_final(level,fvar,0); 
	write_to_file(node,fvar);
	write_restart(fvar,0);
	
	cout<<"The error norm: "<<error_calc(fvar,0)<<"\n";*/
	
	/***Computing the truncation error field***/ 
	//The initial condition is an exact solution here
	
	/*int lev=0; 	
	mg_coeff();
	mg_read_coeff(level); // Reading the null space vectors for all grid sizes
	compute_pbc(pbc,1); 
	trunc_error_map(level,lev,pbc);
	mg_bcs_neu(level,lev,pbc);
	mg_final(level,fvar,lev); 
	write_to_file(node,fvar);
	write_restart(fvar,lev);*/	
	
	/***Computation of truncation error from Direct solver***/
         
	/*int lev=0; 
	int nx_level = nx/pow(2,lev)-1; 
	int ny_level = ny/pow(2,lev)-1; 
	
	mg_coeff();  	      			// Coefficients for smoothing operator at all levels
	mg_read_coeff(level); 			// Reading the null space vectors for all grid sizes
	compute_pbc(pbc,1);
	mg_evaluate_rhs(level,lev,pbc); 
	mg_compute_rhs_tot(level,lev); 	
	arma_direct_trunc(level, lev);		
	mg_bcs_neu(level,0,pbc); 
	mg_final(level,fvar,0);	
       	write_to_file(node,fvar);
       	
       	cout<<"The error norm: "<<error_calc(fvar,0)<<"\n";*/
       	
       	/******Computing for the flow********/ 
       	
        solver(node, fvar, level, pbc, restart_count); 	       	
       	
	/*bcs(node,fvar);
	mg_clear_levels(level); 
	poisson_source(fvar,level,4.0); 
	compute_pbc(fvar,pbc,1); 
	
	//mg_coeff();  	      // Coefficients for smoothing operator at all levels
	//mg_read_coeff(level); // Reading the null space vectors for all grid sizes     		
	//mg_bicgstab_corn(fvar, level, 0, 1, 1, 0, pbc);
	
	mg_poisson_solver(level,fvar,pbc);
	mg_final(level,fvar,0);	        
	tdma1x(fvar,2); 
	tdma1y(fvar,2); 
	
	write_to_file(node,fvar,0);*/
}

/***************************************************/

double get_wall_time()
{
    struct timeval time;

    if(gettimeofday(&time,NULL))
    {
        return 0; //Handle error
    }

    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

/***************************************************/
