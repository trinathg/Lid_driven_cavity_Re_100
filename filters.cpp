#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

///////////////////Filtering////////////////////
//8th order filter 

/*void filter_gaitonde_8(vector<vertex>& node, vector<fval>& fvar)
{
	vector<double> MPf(nx+1),MWf(nx+1),MEf(nx+1); 	
	
////////Filter coefficients 
	for(int i=0;i<=nx;i++)
	{
		if(i==0)
		{
			MPf[i] = 1.0; 
			MEf[i] = alpha_f;
			MWf[i] = 0.0; 				
		}
		else if(i==1)
		{
			MPf[i] = 1.0; 
			MWf[i] = alpha_f; 
			MEf[i] = alpha_f; 
		}
		else if(i==2)
		{
			MPf[i] = 1.0;
			MWf[i] = alpha_f; 
			MEf[i] = alpha_f; 
		}
		else if(i==3)
		{
			MPf[i] = 1.0; 
			MEf[i] = alpha_f; 
			MWf[i] = alpha_f;
		}
		else if(i==4)
		{
			MPf[i] = 1.0; 
			MEf[i] = alpha_f; 
			MWf[i] = alpha_f; 
		}
		else if(i==nx-4)
		{
			MPf[i] = 1.0; 
			MEf[i] = alpha_f; 
			MWf[i] = alpha_f; 
		}
		else if(i==nx-3 )   
		{
			MPf[i] = 1.0; 
			MEf[i] = alpha_f; 
			MWf[i] = alpha_f; 
		}
		else if(i==nx-2)
		{  
			MPf[i] = 1.0; 
			MEf[i] = alpha_f; 
			MWf[i] = alpha_f; 
		} 
		else if(i==nx-1)
		{
			MPf[i] = 1.0; 
			MWf[i] = alpha_f; 
			MEf[i] = alpha_f; 
		}
		else if(i==nx)
		{
			MPf[i] = 1.0; 
			MWf[i] = alpha_f; 
			MEf[i] = 0.0;
		}
		else
		{
			MWf[i]=alpha_f; 
			MEf[i]=alpha_f; 
			MPf[i]=1.;	
		}
	}*/ 
	////////////////////////////////////////

	/*double A0_i = (93. + 70.*alpha_f)/128. ;   //Verified for the correctness of the coefficients
	double A1_i = (7. + 18.*alpha_f)/16. ;
	double A2_i = (-7. + 14.*alpha_f)/32. ; 
	double A3_i = (1./16.) - (alpha_f/8.); 
	double A4_i = (-1./128.) + (alpha_f/64.);

	double a_4 = 1./256. - alpha_f/128.; 	 //Verified for the correctness of the coefficients
	double b_4 = -1./32. + alpha_f/16.; 
	double c_4 = 7./64. + (25.*alpha_f)/32. ; 
	double d_4 = 25./32. + (7.*alpha_f)/16. ; 
	double e_4 = 35./128. + (29.*alpha_f)/64. ; 
	double f_4 = -7./32. + (7.*alpha_f)/16. ; 
	double g_4 = 7./64. - (7.*alpha_f)/32. ; 
	double h_4 = -1./32. + alpha_f/16. ; 
	double i_4 = 1./256. - alpha_f/128. ; 

	double a_3 = -1./256. + alpha_f/128. ;      //Verified for the correctness of the coefficients
	double b_3 = 1./32. + (15.*alpha_f)/16. ; 
	double c_3 = 57./64. + (7.*alpha_f)/32. ; 
	double d_3 = 7./32. + (9.*alpha_f)/16. ; 
	double e_3 = (7./128.)*(-5. + 10.*alpha_f); 
	double f_3 = 7./32. - (7.*alpha_f)/16. ; 
	double g_3 = -7./64. + (7.*alpha_f)/32. ; 
	double h_3 = 1./32. - alpha_f/16. ; 
	double i_3 = -1./256. + alpha_f/128. ; 

	double a_2 = 1./256. + (127.*alpha_f)/128.; //Verified for the correctness of the coefficients      
	double b_2 = 31./32. + alpha_f/16.; 
	double c_2 = 7./64. + (25.*alpha_f)/32. ; 
	double d_2 = -7./32. + (7.*alpha_f)/16. ; 
	double e_2 = (7./128.)*(5.-10.*alpha_f); 
	double f_2 = -7./32. + (7.*alpha_f)/16.;
	double g_2 = 7./64. - (7.*alpha_f)/32. ; 
	double h_2 = -1./32. + alpha_f/16. ; 
	double i_2 = 1./256. - alpha_f/128. ; 

	double a_1 = 255./256. + alpha_f/256. ;   //Verified for the correctness of the coefficients
	double b_1 = 1./32. + 31.*alpha_f/ 32. ; 
	double c_1 = -7./64. + 7.*alpha_f/64. ; 
	double d_1 = 7./32. - 7.*alpha_f/32. ; 
	double e_1 = (7./128.)*(-5. + 5.*alpha_f); 
	double f_1 = 7./32. - (7.*alpha_f)/32.; 
	double g_1 = 7.*(-1.+alpha_f)/64. ; 
	double h_1 = (1.-alpha_f)/32. ; 
	double i_1 = (alpha_f-1.)/256. ; 

//cout<<a_1<<" "<<b_1<<" "<<c_1<<" "<<d_1<<" "<<e_1<<" "<<f_1<<" "<<g_1<<" "<<h_1<<" "<<i_1<<" "<<"\n";


	vector<double> Qx(nx+1), Qstarx(nx+1);
 
	vector<double> phi(nx+1), phi_bar(nx+1);*/
 

	/*******Initializing the array Phi*/

	/*for(int i=0;i<=nx;i++)
	{
		phi[i]=fvar[i].u[0];
		phi_bar[i]=0.0; 
	}*/

	/*******Populating the matrix*/

	/*for(int i=0;i<=nx;i++)
	{

		if(i==0)
		{				
			Qx[i] = a_1*phi[i] + b_1*phi[i+1] + c_1*phi[i+2] + d_1*phi[i+3] + e_1*phi[i+4] + f_1*phi[i+5] + g_1*phi[i+6] + h_1*phi[i+7] + i_1*phi[i+8] ;
		}  
		else if(i==1)
		{
			Qx[i] = a_2*phi[i-1] + b_2*phi[i] + c_2*phi[i+1] + d_2*phi[i+2] + e_2*phi[i+3] + f_2*phi[i+4] + g_2*phi[i+5] + h_2*phi[i+6] + i_2*phi[i+7]; 
		}
		else if(i==2)
		{
			Qx[i] = a_3*phi[i-2] + b_3*phi[i-1] + c_3*phi[i] + d_3*phi[i+1] + e_3*phi[i+2] + f_3*phi[i+3] + g_3*phi[i+4] + h_3*phi[i+5] + i_3*phi[i+6]; 
		}
		else if(i==3)
		{
			Qx[i] =  a_4*phi[i-3] + b_4*phi[i-2] + c_4*phi[i-1] + d_4*phi[i] + e_4*phi[i+1] + f_4*phi[i+2] + g_4*phi[i+3] + h_4*phi[i+4] + i_4*phi[i+5];  
		}
		else if(i==nx-3 )   
		{
			Qx[i] = a_4*phi[i+3] + b_4*phi[i+2] + c_4*phi[i+1] + d_4*phi[i] + e_4*phi[i-1] + f_4*phi[i-2] + g_4*phi[i-3] + h_4*phi[i-4] + i_4*phi[i-5]; 
		}
		else if(i==nx-2)
		{  
			Qx[i] = a_3*phi[i+2] + b_3*phi[i+1] + c_3*phi[i] + d_3*phi[i-1] + e_3*phi[i-2] + f_3*phi[i-3] + g_3*phi[i-4] + h_3*phi[i-5] + i_3*phi[i-6]; 
		} 
		else if(i==nx-1)
		{
			Qx[i] = a_2*phi[i+1] + b_2*phi[i] + c_2*phi[i-1] + d_2*phi[i-2] + e_2*phi[i-3] + f_2*phi[i-4] + g_2*phi[i-5] + h_2*phi[i-6] + i_2*phi[i-7] ;   
		}
		else if(i==nx)
		{
			Qx[i] = a_1*phi[i] + b_1*phi[i-1] + c_1*phi[i-2] + d_1*phi[i-3] + e_1*phi[i-4] + f_1*phi[i-5] + g_1*phi[i-6] + h_1*phi[i-7] + i_1*phi[i-8]; 
		}
		else
		{
			Qx[i] =   A0_i*phi[i] + A1_i*(phi[i+1] + phi[i-1] )*0.5 + A2_i*(phi[i+2] + phi[i-2])*0.5 + A3_i*(phi[i+3]+phi[i-3])*0.5 + A4_i*(phi[i+4] + phi[i-4] )*0.5 ; 
     
		}

	}

// Thomas Algorithm for X-derivative 

//Forward Elimination

	for(int i=0;i<=nx;i++)
	{
		if(i==0)
		{
			Qstarx[i]=Qx[i];
		}
		else 
		{
			MPf[i] =  MPf[i]- (MWf[i]*MEf[i-1])/(MPf[i-1]) ; 
			Qstarx[i] = Qx[i] - (MWf[i]*Qstarx[i-1])/(MPf[i-1]) ; 
		}
	}

//Backward Substitution
	
	for(int i=nx;i>=0;i--)
	{					
		if(i==nx)
		{
			phi_bar[i] = Qstarx[i]/MPf[i]; 
		}	
		else 
		{
			phi_bar[i] = ( Qstarx[i] - MEf[i]*phi_bar[i+1] )/MPf[i]; 
		}
		
		fvar[i].u[0] = phi_bar[i];		
	}

}*/
///////////////////////////////////////////////////////////////////////////////////
