extern void make_mesh(void); 
extern void initialize(vertex *,fval *);
extern void write_to_file(vertex *, fval *);  
extern void write_to_file(vertex *, fval *,int); //Function overloading to write files at certain frequency 

extern void tdma1x(fval *,int);
extern void tdma1y(fval *,int); 
 
extern void tdma2x(fval *,int); 
extern void tdma2y(fval *,int); 

extern void filter_gaitonde_8(vertex *, fval *); 
extern void solver(vertex *, fval *, mg_grid *, pbcs&, int);

extern void explicit_solver(vertex *, fval *, mg_grid *, pbcs&);
extern void bcs(vertex *,fval *); 
extern void compute_pbc(fval *, pbcs&,int);
extern void compute_pbc(pbcs&, int);
extern double div_calc(vertex *,fval *);

extern double residual_calc(fval *); 
extern void global_residuals(fval *, ofstream & , int);
extern double error_calc(fval *);

extern double get_wall_time();

extern void bcs_neu(fval *); 
extern void poisson_Gauss_Seidel6(fval *);
extern void poisson_source(fval * , mg_grid *, double);
extern void read_coeff(double *, double *, double &);

extern void evaluate_rhs6(fval *);
extern void evaluate_rhs_total(fval *,double *); 
extern void evlauate_F(fval*); 

extern void write_restart(fval *,int); 
extern void read_restart(fval *, int); 
extern void mg_restart(fval *, mg_grid *); 

/****************************************************************************/
/***************************Multigrid functions******************************/ 

extern void mg_coeff(); 
extern void mg_interpolate(mg_grid *, int, int);

extern void mg_restrict(mg_grid *, int);
extern void mg_final(mg_grid *, fval *,int); 

extern void mg_read_coeff(mg_grid *); 
extern void mg_evaluate_rhs(mg_grid *, int);
extern void mg_evaluate_rhs(mg_grid *, int, pbcs&, int);

extern void mg_residual_neu(mg_grid *, int);

extern void mg_compute_rhs_tot(mg_grid *, int); 
extern double mg_eval_rhs_sum(mg_grid *, int);

extern void mg_gauss_seidel6(mg_grid *, int, int, int, int, pbcs&); //Overloaded function for FMG
extern double mg_res_norm(mg_grid *, int);

extern void mg_poisson_solver(mg_grid *, fval *, pbcs&);
extern void v_cycle(fval * , mg_grid *, int, pbcs&); //Overloaded function for FMG
extern void fmg(fval*, mg_grid *, pbcs&); //Function for FMG 
extern void fmg(mg_grid *, pbcs&);
extern void mg_bcs_neu(mg_grid * ,int , pbcs&);

extern void mg_bcs_neu(mg_grid *, int);
extern void mg_initialize(mg_grid *);

extern void impose_zero_mean(mg_grid *, int);
extern void arma_direct_solve(mg_grid *, int);
extern void arma_direct_trunc(mg_grid *, int);
extern void mg_clear_levels(mg_grid *);

extern void mg_conjugate_gradient(mg_grid *, int, int, int, int, pbcs& pbc);

extern void mg_bicgstab(mg_grid *, int, int, int, int, pbcs& pbc);
extern void mg_bicgstab_corn(fval *, mg_grid *, int, int, int, int, pbcs& pbc);

extern vec mulvec(vec,mg_grid *);
extern vec mulvec(vec,mg_grid *,int);
extern vec error_from_res(vec &, mg_grid *); 
extern double error_calc(fval *, int);
extern void trunc_error_map(mg_grid *, int, pbcs&);

extern double extra_pol(int, int, mg_grid *, int);

/****************************************************************************/	
