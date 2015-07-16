const int DIM=2;
const int NUM_VAR=5;     // U, V, P,H0,H1 corresponding to the indices 0, 1, 2, 3, 4; Here 3 & 4 correspond to advective + viscous fluxes in both the directions 

const int nx=128, ny=128, nz=20; 
const int str_x = nx+1;  
 
const int file_freq=500; // Frequency of writing down files for data visualization
const int r_file_freq = 3*file_freq; //Frequency for wriitng the restart file

const double l_y = 1.0;    //The height of the 2D channel  
const double l_x = 1.0; //Non dimensional units 
const double delta_p = 2.0; //This is the difference between the outlet and inlet pressure. The flow is from left to right

const double Re = 100.0; //Reynolds number based on the size of the cavity 
const double RRe = 1./Re;  

const double dx = l_x/nx; //The grid spacing in X and Y-directions 
const double dy = l_y/ny; 

const double dt = 5.0e-4; //CFL_d*(dx)*(dx)*Re; // 1/Re acts as a diffusion coefficient here 

const double tot_time = 20.0; //Some-non-dimenisonal time 

const int ite = tot_time/dt;
//const int ite = 100; //Giving a small no. of iterations for testing

const int poi_ite=100;  

const double alpha_f = 0.3;  

const int restart=0; //0 --> No restart; 1->Restart

typedef struct coord          // Defining the structure for the coordinate system
{
	double x[DIM];               // X[0]: for X-coordinate // X[1]: for Y-coordinate //X[2]: for Z-coordinate 
}vertex; 

typedef struct field          // Defining the structure for the field variables
{
	double u[NUM_VAR], uxx[NUM_VAR], uyy[NUM_VAR];  // P for pressure. Doesn't matter here thermodynamic or whatever. Only grad matters
	double ux[NUM_VAR], uy[NUM_VAR]; 
	double u0[NUM_VAR]; 
	double F,div,res,dp,rhs,err; 		
}fval; 

typedef struct bc
{
	double p1x[ny+1], p2x[ny+1], p1y[nx+1], p2y[nx+1]; 
}pbcs; 

const int tot_p = (nx+1)*(ny+1);    //Total no.of points in the domain 
const int tot_inp = (nx-1)*(ny-1);  //Total no.of interior points

extern std::vector<double> MWx, MEx, MPx, MWy, MEy, MPy, MW2x, ME2x, MP2x, MW2y, ME2y, MP2y;
extern std::vector<double> MWf, MEf, MPf;

const double td_ax = 14./9. , td_bx = 1./9.;
const double td_ay = 14./9. , td_by = 1./9.;

const double td_ax2 = 12./11. , td_bx2 = 3./11.;
const double td_ay2 = 12./11. , td_by2 = 3./11.;
 
const double a_nb = 12./10., alpha_nb = 1./10.; //The coefficients at the near boundary points 
const double a = 12./11., b = (1./4.)*(3./11.), alpha = 2./11.; // Coeficients in the interior 

const double alpha_b = 11.0;
//const double hx2 = dx*dx, hy2 = dy*dy;  

const double B1 = -11./6.;           //One sided coefficients normalized with B1
const double B2 = 3./B1;
const double B3 = (-3./2.)/B1; 
const double B4 = (1./3.)/B1; 

const double e1 = 5.0, e2 = -10.0, e3 = 10.0, e4 = -5.0, e5 = 1.0;

/*********************************************************************************************/ 
//OMP variables 
const int numprocs = omp_get_num_procs();	  
const int num_threads = 1; 
//const int chunksize = (ny+1)/num_threads;

/*********************************************************************************************/
/*********************************************************************************************/
//Data structures and indices for multigrid 

const int min_stencil_req = 8; 

const int mg_levels_x = log2(nx) - log2(min_stencil_req) + 1; // Calculating the no. of multi grid levels needed in x-direction 
const int mg_levels_y = log2(ny) - log2(min_stencil_req) + 1; // Calculating the no. of multi grid levels needed in y-direction 

const int mg_levels = min(mg_levels_x,mg_levels_y); 

const int mg_full_x = log2(nx); 
const int mg_full_y = log2(ny); 

const int mg_full = min(mg_full_x, mg_full_y);  

//The "-1" in the levels above is after consideration of the stencil which requries 4 grid spacings to be there for the coarsest grid to preserve the (3,4,6,4,3) accuracy. 

//If the no.of spacings (nx and ny) are different from each other then common sense instructs us to take the min(mg_levels_x,mg_levels_y). This is because beyond the coarsest of the levels in one direction there aren't any points left for the stencil. Another thing to look for is to have enough no.of points available for the scheme that we are using. We can't work with absolutely just 1 point in the domain and because the scheme demands a minimum 20 point stencil. The coarsest grid should be able to accommodate the scheme in full without any loss of the stencil. All the Coefficients required for the operator are going to be changed os have to take that into account too. 

typedef struct mg
{
	double phi_s[tot_p], F[tot_p], res[tot_p], coeff[tot_p], cor_rhs[tot_p], rhs[tot_p], p1x[ny+1], p2x[ny+1], p1y[nx+1], p2y[nx+1],err[tot_p]; 
	double point_correc; 
}mg_grid; 

extern std::vector<double> hx, hy, hx2, hy2, hxx, hyy; //The grid sizes at various levels of the multigrid cycle
extern std::vector<double> X1,X2,X3,X4,X5,X6; 
extern std::vector<double> Y1,Y2,Y3,Y4,Y5,Y6,Y7;
extern std::vector<double> XA,YA,XH,YH,YB;
extern std::vector<double> Z1,Z2; 

extern std::vector<double> C11,C12,C13,C14,C15,C16; 
extern std::vector<double> C21,C22,C23,C24,C25,C26,C27,C28,C29,C210,C211,C212; 
extern std::vector<double> C31,C32,C33,C34,C35,C36,C37,C38,C39,C310;
extern std::vector<double> C41,C42,C43,C44,C45,C46,C47,C48; 

/*Case-1*/
extern std::vector<double> dc_b1_11, dc_b1_12, dc_b1_13, dc_nb1_21, dc_nb1_22, dc_nb1_23, dc_nb1_24, dc_i1_32, dc_i1_33, dc_i1_34, dc_i1_35, dc_i1_36; 
extern std::vector<double> dc_b2_11, dc_b2_12, dc_b2_13, dc_nb2_21, dc_nb2_22, dc_nb2_23, dc_nb2_24, dc_i2_32, dc_i2_33, dc_i2_34, dc_i2_35, dc_i2_36; 
extern std::vector<double> dc_b3_11, dc_b3_12, dc_b3_13, dc_nb3_21, dc_nb3_22, dc_nb3_23, dc_nb3_24, dc_i3_32, dc_i3_33, dc_i3_34, dc_i3_35, dc_i3_36; 
 
/*Case-2*/
extern std::vector<double> dc_b4_11, dc_b4_12, dc_b4_13, dc_nb4_21, dc_nb4_22, dc_nb4_23, dc_nb4_24, dc_i4_32, dc_i4_33, dc_i4_34, dc_i4_35, dc_i4_36; 
extern std::vector<double> dc_b5_11, dc_b5_12, dc_b5_13, dc_nb5_21, dc_nb5_22, dc_nb5_23, dc_nb5_24, dc_i5_32, dc_i5_33, dc_i5_34, dc_i5_35, dc_i5_36; 
extern std::vector<double> dc_b6_11, dc_b6_12, dc_b6_13, dc_nb6_21, dc_nb6_22, dc_nb6_23, dc_nb6_24, dc_i6_32, dc_i6_33, dc_i6_34, dc_i6_35, dc_i6_36; 
extern std::vector<double> dc_b7_11, dc_b7_12, dc_b7_13, dc_nb7_21, dc_nb7_22, dc_nb7_23, dc_nb7_24, dc_i7_32, dc_i7_33, dc_i7_34, dc_i7_35, dc_i7_36; 

/*Case-3*/
extern std::vector<double> dc_b8_11, dc_b8_12, dc_b8_13, dc_nb8_21, dc_nb8_22, dc_nb8_23, dc_nb8_24, dc_i8_32, dc_i8_33, dc_i8_34, dc_i8_35, dc_i8_36; 
extern std::vector<double> dc_b9_11, dc_b9_12, dc_b9_13, dc_nb9_21, dc_nb9_22, dc_nb9_23, dc_nb9_24, dc_i9_32, dc_i9_33, dc_i9_34, dc_i9_35, dc_i9_36; 
extern std::vector<double> dc_b10_11, dc_b10_12, dc_b10_13, dc_nb10_21, dc_nb10_22, dc_nb10_23, dc_nb10_24, dc_i10_32, dc_i10_33, dc_i10_34, dc_i10_35, dc_i10_36; 
extern std::vector<double> dc_b11_11, dc_b11_12, dc_b11_13, dc_nb11_21, dc_nb11_22, dc_nb11_23, dc_nb11_24, dc_i11_32, dc_i11_33, dc_i11_34, dc_i11_35, dc_i11_36; 
extern std::vector<double> dc_b12_11, dc_b12_12, dc_b12_13, dc_nb12_21, dc_nb12_22, dc_nb12_23, dc_nb12_24, dc_i12_32, dc_i12_33, dc_i12_34, dc_i12_35, dc_i12_36; 

/*Coefficients for the RHS side*/

extern std::vector<double> dc_b13_11, dc_b13_12, dc_b13_13, dc_nb13_21, dc_nb13_22, dc_nb13_23, dc_nb13_24, dc_nb13_25; 
extern std::vector<double> dc_b14_11, dc_b14_12, dc_b14_13, dc_nb14_21, dc_nb14_22, dc_nb14_23, dc_nb14_24, dc_nb14_25; 
