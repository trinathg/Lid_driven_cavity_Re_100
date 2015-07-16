#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

/**************************************************************************************************************************************/
/*Specifying and initializing the coefficients globally*/

vector<double> hx(mg_levels), hy(mg_levels), hX2(mg_levels), hY2(mg_levels), hxx(mg_levels), hyy(mg_levels); //The grid sizes at various levels of the multigrid cycle
vector<double> X1(mg_levels),X2(mg_levels),X3(mg_levels),X4(mg_levels),X5(mg_levels),X6(mg_levels); 
vector<double> Y1(mg_levels),Y2(mg_levels),Y3(mg_levels),Y4(mg_levels),Y5(mg_levels),Y6(mg_levels),Y7(mg_levels);
vector<double> XA(mg_levels),YA(mg_levels),XH(mg_levels),YH(mg_levels),YB(mg_levels);
vector<double> Z1(mg_levels),Z2(mg_levels); 

vector<double> C11(mg_levels),C12(mg_levels),C13(mg_levels),C14(mg_levels),C15(mg_levels),C16(mg_levels); 
vector<double> C21(mg_levels),C22(mg_levels),C23(mg_levels),C24(mg_levels),C25(mg_levels),C26(mg_levels),C27(mg_levels),C28(mg_levels),C29(mg_levels),C210(mg_levels),C211(mg_levels),C212(mg_levels); 
vector<double> C31(mg_levels),C32(mg_levels),C33(mg_levels),C34(mg_levels),C35(mg_levels),C36(mg_levels),C37(mg_levels),C38(mg_levels),C39(mg_levels),C310(mg_levels);
vector<double> C41(mg_levels),C42(mg_levels),C43(mg_levels),C44(mg_levels),C45(mg_levels),C46(mg_levels),C47(mg_levels),C48(mg_levels); 

/*Case-1*/
vector<double> dc_b1_11(mg_levels), dc_b1_12(mg_levels), dc_b1_13(mg_levels), dc_nb1_21(mg_levels), dc_nb1_22(mg_levels), dc_nb1_23(mg_levels), dc_nb1_24(mg_levels), dc_i1_32(mg_levels), dc_i1_33(mg_levels), dc_i1_34(mg_levels), dc_i1_35(mg_levels), dc_i1_36(mg_levels); 
vector<double> dc_b2_11(mg_levels), dc_b2_12(mg_levels), dc_b2_13(mg_levels), dc_nb2_21(mg_levels), dc_nb2_22(mg_levels), dc_nb2_23(mg_levels), dc_nb2_24(mg_levels), dc_i2_32(mg_levels), dc_i2_33(mg_levels), dc_i2_34(mg_levels), dc_i2_35(mg_levels), dc_i2_36(mg_levels); 
vector<double> dc_b3_11(mg_levels), dc_b3_12(mg_levels), dc_b3_13(mg_levels), dc_nb3_21(mg_levels), dc_nb3_22(mg_levels), dc_nb3_23(mg_levels), dc_nb3_24(mg_levels), dc_i3_32(mg_levels), dc_i3_33(mg_levels), dc_i3_34(mg_levels), dc_i3_35(mg_levels), dc_i3_36(mg_levels); 
 
/*Case-2*/
vector<double> dc_b4_11(mg_levels), dc_b4_12(mg_levels), dc_b4_13(mg_levels), dc_nb4_21(mg_levels), dc_nb4_22(mg_levels), dc_nb4_23(mg_levels), dc_nb4_24(mg_levels), dc_i4_32(mg_levels), dc_i4_33(mg_levels), dc_i4_34(mg_levels), dc_i4_35(mg_levels), dc_i4_36(mg_levels); 
vector<double> dc_b5_11(mg_levels), dc_b5_12(mg_levels), dc_b5_13(mg_levels), dc_nb5_21(mg_levels), dc_nb5_22(mg_levels), dc_nb5_23(mg_levels), dc_nb5_24(mg_levels), dc_i5_32(mg_levels), dc_i5_33(mg_levels), dc_i5_34(mg_levels), dc_i5_35(mg_levels), dc_i5_36(mg_levels); 
vector<double> dc_b6_11(mg_levels), dc_b6_12(mg_levels), dc_b6_13(mg_levels), dc_nb6_21(mg_levels), dc_nb6_22(mg_levels), dc_nb6_23(mg_levels), dc_nb6_24(mg_levels), dc_i6_32(mg_levels), dc_i6_33(mg_levels), dc_i6_34(mg_levels), dc_i6_35(mg_levels), dc_i6_36(mg_levels); 
vector<double> dc_b7_11(mg_levels), dc_b7_12(mg_levels), dc_b7_13(mg_levels), dc_nb7_21(mg_levels), dc_nb7_22(mg_levels), dc_nb7_23(mg_levels), dc_nb7_24(mg_levels), dc_i7_32(mg_levels), dc_i7_33(mg_levels), dc_i7_34(mg_levels), dc_i7_35(mg_levels), dc_i7_36(mg_levels); 

/*Case-3*/
vector<double> dc_b8_11(mg_levels), dc_b8_12(mg_levels), dc_b8_13(mg_levels), dc_nb8_21(mg_levels), dc_nb8_22(mg_levels), dc_nb8_23(mg_levels), dc_nb8_24(mg_levels), dc_i8_32(mg_levels), dc_i8_33(mg_levels), dc_i8_34(mg_levels), dc_i8_35(mg_levels), dc_i8_36(mg_levels); 
vector<double> dc_b9_11(mg_levels), dc_b9_12(mg_levels), dc_b9_13(mg_levels), dc_nb9_21(mg_levels), dc_nb9_22(mg_levels), dc_nb9_23(mg_levels), dc_nb9_24(mg_levels), dc_i9_32(mg_levels), dc_i9_33(mg_levels), dc_i9_34(mg_levels), dc_i9_35(mg_levels), dc_i9_36(mg_levels); 
vector<double> dc_b10_11(mg_levels), dc_b10_12(mg_levels), dc_b10_13(mg_levels), dc_nb10_21(mg_levels), dc_nb10_22(mg_levels), dc_nb10_23(mg_levels), dc_nb10_24(mg_levels), dc_i10_32(mg_levels), dc_i10_33(mg_levels), dc_i10_34(mg_levels), dc_i10_35(mg_levels), dc_i10_36(mg_levels); 
vector<double> dc_b11_11(mg_levels), dc_b11_12(mg_levels), dc_b11_13(mg_levels), dc_nb11_21(mg_levels), dc_nb11_22(mg_levels), dc_nb11_23(mg_levels), dc_nb11_24(mg_levels), dc_i11_32(mg_levels), dc_i11_33(mg_levels), dc_i11_34(mg_levels), dc_i11_35(mg_levels), dc_i11_36(mg_levels); 
vector<double> dc_b12_11(mg_levels), dc_b12_12(mg_levels), dc_b12_13(mg_levels), dc_nb12_21(mg_levels), dc_nb12_22(mg_levels), dc_nb12_23(mg_levels), dc_nb12_24(mg_levels), dc_i12_32(mg_levels), dc_i12_33(mg_levels), dc_i12_34(mg_levels), dc_i12_35(mg_levels), dc_i12_36(mg_levels); 

/*Coefficients for the RHS side*/

vector<double> dc_b13_11(mg_levels), dc_b13_12(mg_levels), dc_b13_13(mg_levels), dc_nb13_21(mg_levels), dc_nb13_22(mg_levels), dc_nb13_23(mg_levels), dc_nb13_24(mg_levels), dc_nb13_25(mg_levels); 
vector<double> dc_b14_11(mg_levels), dc_b14_12(mg_levels), dc_b14_13(mg_levels), dc_nb14_21(mg_levels), dc_nb14_22(mg_levels), dc_nb14_23(mg_levels), dc_nb14_24(mg_levels), dc_nb14_25(mg_levels); 

/****************************************************************************************************************************************/
/*******************************************Multi grid coefficients**********************************************************************/

void mg_coeff()
{
	ofstream coeff_c, coeff_d; 
	
	/*coeff_c.open("coeff_c.dat");
	coeff_d.open("coeff_d.dat");*/ 
	
	for(int lev=0;lev<mg_levels;lev++)
	{	
		hx[lev] = dx*pow(2,lev); 
		hy[lev] = dy*pow(2,lev);
		
		hX2[lev] = hx[lev]*hx[lev];
		hY2[lev] = hy[lev]*hy[lev]; 
		
		hxx[lev] = 1./(hx[lev]*hx[lev]); 
		hyy[lev] = 1./(hy[lev]*hy[lev]); 
				
		/******Dirichlet coeffs for all levels*******/
		
		C11[lev] = -2.0*a_nb*( 1./hX2[lev] + 1./hY2[lev] ); 		
		C12[lev] = (a_nb/hX2[lev]) - (2.*alpha_nb*a_nb/hY2[lev]) ; 
		C13[lev] = a_nb*alpha_nb*( 1./hX2[lev] + 1./hY2[lev] ); 
		C14[lev] = a_nb/hY2[lev] - (2.*a_nb*alpha_nb)/hX2[lev];
		C15[lev] = alpha_nb*alpha_nb; 
		C16[lev] = alpha_nb; 
				
		C21[lev] = -2.0*(a+b)/hX2[lev] -(2.*a_nb)/hY2[lev] ; 
		C22[lev] = (a/hX2[lev]) - (2.*alpha*a_nb)/hY2[lev]; 
		C23[lev] = (a/hX2[lev]) - (2.*alpha*a_nb)/hY2[lev]; 
		C24[lev] = b/hX2[lev]; 
		C25[lev] = alpha_nb*b/hX2[lev]; 
		C26[lev] = (a_nb*alpha)/hY2[lev]  + (alpha_nb*a)/hX2[lev]; 
		C27[lev] = a_nb/hY2[lev] - (2.*alpha_nb*(a+b))/hX2[lev] ; 
		C28[lev] = (a_nb*alpha)/hY2[lev] + (alpha_nb*a)/hX2[lev] ; 
		C29[lev] = (a_nb/hY2[lev]) - (2*(a+b))*alpha_nb/hX2[lev] ; 
		C210[lev] = b*alpha_nb/hX2[lev]; 
		C211[lev] = (a*alpha_nb)/hX2[lev] + (alpha*a_nb)/hY2[lev] ; 
		C212[lev] = (a*alpha_nb)/hX2[lev] + (alpha*a_nb)/hY2[lev] ;

		C31[lev] = (-2.*a_nb)/hX2[lev] -(2.*(a+b))/hY2[lev]; 
		C32[lev] = a_nb/hX2[lev] - (2.*alpha_nb*(a+b))/hY2[lev] ; 
		C33[lev] = (alpha*a_nb)/hX2[lev] + (a*alpha_nb)/hY2[lev] ; 
		C34[lev] = a/hY2[lev] - (2.*a_nb*alpha)/hX2[lev]; 
		C35[lev] = (alpha*a_nb)/hX2[lev] + (alpha_nb*a)/hY2[lev]; 
		C36[lev] = b/hY2[lev]; 
		C37[lev] = (a*alpha_nb)/hY2[lev] + (alpha*a_nb)/hX2[lev]; 
		C38[lev] = a/hY2[lev] - (2.*a_nb*alpha)/hX2[lev]; 
		C39[lev] = (a*alpha_nb)/hY2[lev] + (alpha*a_nb)/hX2[lev]; 
		C310[lev] = a_nb/hX2[lev] - (2.*(a+b)*alpha_nb/hY2[lev]);
		
		C41[lev] = -2.0*(a + b)*( 1./hX2[lev] + 1./hY2[lev]); 
		C42[lev] = (-2.*alpha/hX2[lev])*(a+b) + a/hY2[lev];
		C43[lev] = (-2.*alpha/hY2[lev])*(a+b) + a/hX2[lev] ; 
		C44[lev] = b/hY2[lev]; 
		C45[lev] = b/hX2[lev]; 
		C46[lev] = b*alpha/hX2[lev]; 
		C47[lev] = b*alpha/hY2[lev]; 
		C48[lev] = a*alpha*(1./hX2[lev] + 1./hY2[lev]); 
		
		/*coeff_c<<C11[lev]<<"\n"; 		
		coeff_c<<C12[lev]<<"\n"; 
		coeff_c<<C13[lev]<<"\n"; 
		coeff_c<<C14[lev]<<"\n";
		coeff_c<<C15[lev]<<"\n"; 
		coeff_c<<C16[lev]<<"\n"; 
		
		coeff_c<<C21[lev]<<"\n"; 
		coeff_c<<C22[lev]<<"\n"; 
		coeff_c<<C23[lev]<<"\n"; 
		coeff_c<<C24[lev]<<"\n"; 
		coeff_c<<C25[lev]<<"\n"; 
		coeff_c<<C26[lev]<<"\n"; 
		coeff_c<<C27[lev]<<"\n"; 
		coeff_c<<C28[lev]<<"\n"; 
		coeff_c<<C29[lev]<<"\n"; 
		coeff_c<<C210[lev]<<"\n"; 
		coeff_c<<C211[lev]<<"\n"; 
		coeff_c<<C212[lev]<<"\n";

		coeff_c<<C31[lev]<<"\n"; 
		coeff_c<<C32[lev]<<"\n"; 
		coeff_c<<C33[lev]<<"\n"; 
		coeff_c<<C34[lev]<<"\n"; 
		coeff_c<<C35[lev]<<"\n"; 
		coeff_c<<C36[lev]<<"\n"; 
		coeff_c<<C37[lev]<<"\n"; 
		coeff_c<<C38[lev]<<"\n"; 
		coeff_c<<C39[lev]<<"\n"; 
		coeff_c<<C310[lev]<<"\n";		
		
		coeff_c<<C41[lev]<<"\n"; 
		coeff_c<<C42[lev]<<"\n"; 
		coeff_c<<C43[lev]<<"\n"; 
		coeff_c<<C44[lev]<<"\n"; 
		coeff_c<<C45[lev]<<"\n"; 
		coeff_c<<C46[lev]<<"\n"; 
		coeff_c<<C47[lev]<<"\n"; 
		coeff_c<<C48[lev]<<"\n";*/ 
			
		/****Neumann Coefficients after BCs have been accommodated***/
		
		X1[lev] = hxx[lev]*(1.-B2*alpha_nb); 
		X2[lev] = -a_nb*hyy[lev]*(2.+ B2); 
		
		X3[lev] = (1.-B3)*alpha_nb*hxx[lev];
		X4[lev] = a_nb*hyy[lev]*(1.-B3); 
		
		X5[lev] = -B4*alpha_nb*hxx[lev];
		X6[lev] = -B4*a_nb*hyy[lev]; 

		Y1[lev] = alpha*hxx[lev];
		Y2[lev] = hyy[lev]*(a - B2*b); 
		
		Y3[lev] = hxx[lev]; 
		Y4[lev] = -(2.*(a+b) + B3*b)*hyy[lev]; 
		
		Y5[lev] = alpha*hxx[lev]; 
		Y6[lev] = (a - B4*b)*hyy[lev]; 
		
		Y7[lev] = b*hyy[lev]; 

		XA[lev] = alpha*hxx[lev];
		YA[lev] = a*hyy[lev];
		
		XH[lev] = hxx[lev];
		YH[lev] = -2.*(a+b)*hyy[lev];
		
		YB[lev] = b*hyy[lev]; 

		Z1[lev] = alpha_nb*hxx[lev]; 
		Z2[lev] = a_nb*hyy[lev]; 
		
		/*coeff_c<<X1[lev]<<"\n"; 
		coeff_c<<X2[lev]<<"\n"; 
		coeff_c<<X3[lev]<<"\n"; 
		coeff_c<<X4[lev]<<"\n"; 
		coeff_c<<X5[lev]<<"\n"; 
		coeff_c<<X6[lev]<<"\n"; 
		
		coeff_c<<Y1[lev]<<"\n"; 
		coeff_c<<Y2[lev]<<"\n"; 
		coeff_c<<Y3[lev]<<"\n"; 
		coeff_c<<Y4[lev]<<"\n"; 
		coeff_c<<Y5[lev]<<"\n"; 
		coeff_c<<Y6[lev]<<"\n";
		coeff_c<<Y7[lev]<<"\n";
		
		coeff_c<<XA[lev]<<"\n"; 
		coeff_c<<YA[lev]<<"\n"; 
		coeff_c<<XH[lev]<<"\n"; 
		coeff_c<<YH[lev]<<"\n"; 
		coeff_c<<YB[lev]<<"\n"; 
		
		coeff_c<<Z1[lev]<<"\n"; 
		coeff_c<<Z2[lev]<<"\n";*/
		

		
		//Case-1 (j=1) 
		//j 

		dc_b1_11[lev] = -2.*X1[lev]*a_nb + X2[lev] - B2*(X1[lev]*a_nb + X2[lev]*alpha_nb); 
		dc_b1_12[lev] = X1[lev]*a_nb + X2[lev]*alpha_nb - B3*(X1[lev]*a_nb + X2[lev]*alpha_nb); 
		dc_b1_13[lev] = -B4*( X1[lev]*a_nb + X2[lev]*alpha_nb ); 

		dc_nb1_21[lev] = X1[lev]*a + X2[lev]*alpha - B2*X1[lev]*b; 
		dc_nb1_22[lev] = -2.*(a+b)*X1[lev] + X2[lev] - B3*X1[lev]*b; 
		dc_nb1_23[lev] = X1[lev]*a + X2[lev]*alpha - B4*X1[lev]*b; 
		dc_nb1_24[lev] = X1[lev]*b; 

		dc_i1_32[lev] = X1[lev]*b;    
		dc_i1_33[lev] = X1[lev]*a + X2[lev]*alpha; 
		dc_i1_34[lev] = -2.*(a+b)*X1[lev] + X2[lev]; 
		dc_i1_35[lev] = X1[lev]*a + X2[lev]*alpha;
		dc_i1_36[lev] = X1[lev]*b; 
		
		/*coeff_d<<dc_b1_11[lev]<<"\n"; 
		coeff_d<<dc_b1_12[lev]<<"\n"; 
		coeff_d<<dc_b1_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb1_21[lev]<<"\n"; 
		coeff_d<<dc_nb1_22[lev]<<"\n"; 
		coeff_d<<dc_nb1_23[lev]<<"\n"; 
		coeff_d<<dc_nb1_24[lev]<<"\n"; 
		
		coeff_d<<dc_i1_32[lev]<<"\n"; 
		coeff_d<<dc_i1_33[lev]<<"\n"; 
		coeff_d<<dc_i1_34[lev]<<"\n"; 
		coeff_d<<dc_i1_35[lev]<<"\n";
		coeff_d<<dc_i1_36[lev]<<"\n";*/
		 
		//j+1

		dc_b2_11[lev] = -2.*X3[lev]*a_nb + X4[lev] - B2*(X3[lev]*a_nb + X4[lev]*alpha_nb); 
		dc_b2_12[lev] = X3[lev]*a_nb + X4[lev]*alpha_nb - B3*(X3[lev]*a_nb + X4[lev]*alpha_nb);
		dc_b2_13[lev] = -B4*( X3[lev]*a_nb + X4[lev]*alpha_nb );

		dc_nb2_21[lev] = X3[lev]*a + X4[lev]*alpha - B2*X3[lev]*b;
		dc_nb2_22[lev] = -2.*(a+b)*X3[lev] + X4[lev] - B3*X3[lev]*b;
		dc_nb2_23[lev] = X3[lev]*a + X4[lev]*alpha - B4*X3[lev]*b;
		dc_nb2_24[lev] = X3[lev]*b;

		dc_i2_32[lev] = X3[lev]*b;
		dc_i2_33[lev] = X3[lev]*a + X4[lev]*alpha;
		dc_i2_34[lev] = -2.*(a+b)*X3[lev] + X4[lev];
		dc_i2_35[lev] = X3[lev]*a + X4[lev]*alpha;
		dc_i2_36[lev] = X3[lev]*b;
		
		/*coeff_d<<dc_b2_11[lev]<<"\n"; 
		coeff_d<<dc_b2_12[lev]<<"\n"; 
		coeff_d<<dc_b2_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb2_21[lev]<<"\n"; 
		coeff_d<<dc_nb2_22[lev]<<"\n"; 
		coeff_d<<dc_nb2_23[lev]<<"\n"; 
		coeff_d<<dc_nb2_24[lev]<<"\n"; 
		
		coeff_d<<dc_i2_32[lev]<<"\n"; 
		coeff_d<<dc_i2_33[lev]<<"\n"; 
		coeff_d<<dc_i2_34[lev]<<"\n"; 
		coeff_d<<dc_i2_35[lev]<<"\n";
		coeff_d<<dc_i2_36[lev]<<"\n";*/

		//j+2 

		dc_b3_11[lev] = -2.*X5[lev]*a_nb + X6[lev] - B2*(X5[lev]*a_nb + X6[lev]*alpha_nb); 
		dc_b3_12[lev] = X5[lev]*a_nb + X6[lev]*alpha_nb - B3*(X5[lev]*a_nb + X6[lev]*alpha_nb);
		dc_b3_13[lev] = -B4*( X5[lev]*a_nb + X6[lev]*alpha_nb );

		dc_nb3_21[lev] = X5[lev]*a + X6[lev]*alpha - B2*X5[lev]*b;
		dc_nb3_22[lev] = -2.*(a+b)*X5[lev] + X6[lev] - B3*X5[lev]*b;
		dc_nb3_23[lev] = X5[lev]*a + X6[lev]*alpha - B4*X5[lev]*b;
		dc_nb3_24[lev] = X5[lev]*b;

		dc_i3_32[lev] = X5[lev]*b;
		dc_i3_33[lev] = X5[lev]*a + X6[lev]*alpha;
		dc_i3_34[lev] = -2.*(a+b)*X5[lev] + X6[lev];
		dc_i3_35[lev] = X5[lev]*a + X6[lev]*alpha;
		dc_i3_36[lev] = X5[lev]*b;		
		
		/*coeff_d<<dc_b3_11[lev]<<"\n"; 
		coeff_d<<dc_b3_12[lev]<<"\n"; 
		coeff_d<<dc_b3_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb3_21[lev]<<"\n"; 
		coeff_d<<dc_nb3_22[lev]<<"\n"; 
		coeff_d<<dc_nb3_23[lev]<<"\n"; 
		coeff_d<<dc_nb3_24[lev]<<"\n"; 
		
		coeff_d<<dc_i3_32[lev]<<"\n"; 
		coeff_d<<dc_i3_33[lev]<<"\n"; 
		coeff_d<<dc_i3_34[lev]<<"\n"; 
		coeff_d<<dc_i3_35[lev]<<"\n";
		coeff_d<<dc_i3_36[lev]<<"\n";*/
				
		//Case-2 (j=2 in code here) 

		//j-1 

		dc_b4_11[lev] = -2.*Y1[lev]*a_nb + Y2[lev] - B2*(Y1[lev]*a_nb + Y2[lev]*alpha_nb); 
		dc_b4_12[lev] = Y1[lev]*a_nb + Y2[lev]*alpha_nb - B3*(Y1[lev]*a_nb + Y2[lev]*alpha_nb); 
		dc_b4_13[lev] = -B4*( Y1[lev]*a_nb + Y2[lev]*alpha_nb ); 

		dc_nb4_21[lev] = Y1[lev]*a + Y2[lev]*alpha - B2*Y1[lev]*b; 
		dc_nb4_22[lev] = -2.*(a+b)*Y1[lev] + Y2[lev] - B3*Y1[lev]*b; 
		dc_nb4_23[lev] = Y1[lev]*a + Y2[lev]*alpha - B4*Y1[lev]*b; 
		dc_nb4_24[lev] = Y1[lev]*b; 

		dc_i4_32[lev] = Y1[lev]*b;    
		dc_i4_33[lev] = Y1[lev]*a + Y2[lev]*alpha; 
		dc_i4_34[lev] = -2.*(a+b)*Y1[lev] + Y2[lev]; 
		dc_i4_35[lev] = Y1[lev]*a + Y2[lev]*alpha;
		dc_i4_36[lev] = Y1[lev]*b; 
		
		/*coeff_d<<dc_b4_11[lev]<<"\n"; 
		coeff_d<<dc_b4_12[lev]<<"\n"; 
		coeff_d<<dc_b4_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb4_21[lev]<<"\n"; 
		coeff_d<<dc_nb4_22[lev]<<"\n"; 
		coeff_d<<dc_nb4_23[lev]<<"\n"; 
		coeff_d<<dc_nb4_24[lev]<<"\n"; 
		
		coeff_d<<dc_i4_32[lev]<<"\n"; 
		coeff_d<<dc_i4_33[lev]<<"\n"; 
		coeff_d<<dc_i4_34[lev]<<"\n"; 
		coeff_d<<dc_i4_35[lev]<<"\n";
		coeff_d<<dc_i4_36[lev]<<"\n";*/

		//j 

		dc_b5_11[lev] = -2.*Y3[lev]*a_nb + Y4[lev] - B2*(Y3[lev]*a_nb + Y4[lev]*alpha_nb); 
		dc_b5_12[lev] = Y3[lev]*a_nb + Y4[lev]*alpha_nb - B3*(Y3[lev]*a_nb + Y4[lev]*alpha_nb); 
		dc_b5_13[lev] = -B4*( Y3[lev]*a_nb + Y4[lev]*alpha_nb ); 

		dc_nb5_21[lev] = Y3[lev]*a + Y4[lev]*alpha - B2*Y3[lev]*b; 
		dc_nb5_22[lev] = -2.*(a+b)*Y3[lev] + Y4[lev] - B3*Y3[lev]*b; 
		dc_nb5_23[lev] = Y3[lev]*a + Y4[lev]*alpha - B4*Y3[lev]*b; 
		dc_nb5_24[lev] = Y3[lev]*b; 

		dc_i5_32[lev] = Y3[lev]*b;    
		dc_i5_33[lev] = Y3[lev]*a + Y4[lev]*alpha; 
		dc_i5_34[lev] = -2.*(a+b)*Y3[lev] + Y4[lev]; 
		dc_i5_35[lev] = Y3[lev]*a + Y4[lev]*alpha;
		dc_i5_36[lev] = Y3[lev]*b;
		
		/*coeff_d<<dc_b5_11[lev]<<"\n"; 
		coeff_d<<dc_b5_12[lev]<<"\n"; 
		coeff_d<<dc_b5_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb5_21[lev]<<"\n"; 
		coeff_d<<dc_nb5_22[lev]<<"\n"; 
		coeff_d<<dc_nb5_23[lev]<<"\n"; 
		coeff_d<<dc_nb5_24[lev]<<"\n"; 
		
		coeff_d<<dc_i5_32[lev]<<"\n"; 
		coeff_d<<dc_i5_33[lev]<<"\n"; 
		coeff_d<<dc_i5_34[lev]<<"\n"; 
		coeff_d<<dc_i5_35[lev]<<"\n";
		coeff_d<<dc_i5_36[lev]<<"\n";*/
		

		//j+1

		dc_b6_11[lev] = -2.*Y5[lev]*a_nb + Y6[lev] - B2*(Y5[lev]*a_nb + Y6[lev]*alpha_nb); 
		dc_b6_12[lev] = Y5[lev]*a_nb + Y6[lev]*alpha_nb - B3*(Y5[lev]*a_nb + Y6[lev]*alpha_nb); 
		dc_b6_13[lev] = -B4*( Y5[lev]*a_nb + Y6[lev]*alpha_nb ); 

		dc_nb6_21[lev] = Y5[lev]*a + Y6[lev]*alpha - B2*Y5[lev]*b; 
		dc_nb6_22[lev] = -2.*(a+b)*Y5[lev] + Y6[lev] - B3*Y5[lev]*b; 
		dc_nb6_23[lev] = Y5[lev]*a + Y6[lev]*alpha - B4*Y5[lev]*b; 
		dc_nb6_24[lev] = Y5[lev]*b; 

		dc_i6_32[lev] = Y5[lev]*b;    
		dc_i6_33[lev] = Y5[lev]*a + Y6[lev]*alpha; 
		dc_i6_34[lev] = -2.*(a+b)*Y5[lev] + Y6[lev]; 
		dc_i6_35[lev] = Y5[lev]*a + Y6[lev]*alpha;
		dc_i6_36[lev] = Y5[lev]*b;
		
		/*coeff_d<<dc_b6_11[lev]<<"\n"; 
		coeff_d<<dc_b6_12[lev]<<"\n"; 
		coeff_d<<dc_b6_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb6_21[lev]<<"\n"; 
		coeff_d<<dc_nb6_22[lev]<<"\n"; 
		coeff_d<<dc_nb6_23[lev]<<"\n"; 
		coeff_d<<dc_nb6_24[lev]<<"\n"; 
		
		coeff_d<<dc_i6_32[lev]<<"\n"; 
		coeff_d<<dc_i6_33[lev]<<"\n"; 
		coeff_d<<dc_i6_34[lev]<<"\n"; 
		coeff_d<<dc_i6_35[lev]<<"\n";
		coeff_d<<dc_i6_36[lev]<<"\n";*/

		//j+2 

		dc_b7_11[lev] =  Y7[lev] - B2*(Y7[lev]*alpha_nb); 
		dc_b7_12[lev] =  Y7[lev]*alpha_nb - B3*(Y7[lev]*alpha_nb); 
		dc_b7_13[lev] = -B4*(Y7[lev]*alpha_nb ); 

		dc_nb7_21[lev] = Y7[lev]*alpha ; 
		dc_nb7_22[lev] = Y7[lev]; 
		dc_nb7_23[lev] = Y7[lev]*alpha; 
		dc_nb7_24[lev] = 0.0; 

		dc_i7_32[lev] = 0.0;    
		dc_i7_33[lev] = Y7[lev]*alpha; 
		dc_i7_34[lev] = Y7[lev]; 
		dc_i7_35[lev] = Y7[lev]*alpha;
		dc_i7_36[lev] = 0.0;
		
		/*coeff_d<<dc_b7_11[lev]<<"\n"; 
		coeff_d<<dc_b7_12[lev]<<"\n"; 
		coeff_d<<dc_b7_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb7_21[lev]<<"\n"; 
		coeff_d<<dc_nb7_22[lev]<<"\n"; 
		coeff_d<<dc_nb7_23[lev]<<"\n"; 
		coeff_d<<dc_nb7_24[lev]<<"\n"; 
		
		coeff_d<<dc_i7_32[lev]<<"\n"; 
		coeff_d<<dc_i7_33[lev]<<"\n"; 
		coeff_d<<dc_i7_34[lev]<<"\n"; 
		coeff_d<<dc_i7_35[lev]<<"\n";
		coeff_d<<dc_i7_36[lev]<<"\n";*/
		
		//Case-3 (j=3 in code here)

		//j-2 
		dc_b8_11[lev] =  YB[lev] - B2*(YB[lev]*alpha_nb); 
		dc_b8_12[lev] =  YB[lev]*alpha_nb - B3*(YB[lev]*alpha_nb); 
		dc_b8_13[lev] = -B4*(YB[lev]*alpha_nb); 

		dc_nb8_21[lev] = YB[lev]*alpha ; 
		dc_nb8_22[lev] = YB[lev]; 
		dc_nb8_23[lev] = YB[lev]*alpha; 
		dc_nb8_24[lev] = 0.0; 

		dc_i8_32[lev] = 0.0;    
		dc_i8_33[lev] = YB[lev]*alpha; 
		dc_i8_34[lev] = YB[lev]; 
		dc_i8_35[lev] = YB[lev]*alpha;
		dc_i8_36[lev] = 0.0;
		
		/*coeff_d<<dc_b8_11[lev]<<"\n"; 
		coeff_d<<dc_b8_12[lev]<<"\n"; 
		coeff_d<<dc_b8_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb8_21[lev]<<"\n"; 
		coeff_d<<dc_nb8_22[lev]<<"\n"; 
		coeff_d<<dc_nb8_23[lev]<<"\n"; 
		coeff_d<<dc_nb8_24[lev]<<"\n"; 
		
		coeff_d<<dc_i8_32[lev]<<"\n"; 
		coeff_d<<dc_i8_33[lev]<<"\n"; 
		coeff_d<<dc_i8_34[lev]<<"\n"; 
		coeff_d<<dc_i8_35[lev]<<"\n";
		coeff_d<<dc_i8_36[lev]<<"\n";*/		

		//j-1

		dc_b9_11[lev] = -2.*XA[lev]*a_nb + YA[lev] - B2*(XA[lev]*a_nb + YA[lev]*alpha_nb); 
		dc_b9_12[lev] = XA[lev]*a_nb + YA[lev]*alpha_nb - B3*(XA[lev]*a_nb + YA[lev]*alpha_nb); 
		dc_b9_13[lev] = -B4*( XA[lev]*a_nb + YA[lev]*alpha_nb ); 

		dc_nb9_21[lev] = XA[lev]*a + YA[lev]*alpha - B2*XA[lev]*b; 
		dc_nb9_22[lev] = -2.*(a+b)*XA[lev] + YA[lev] - B3*XA[lev]*b; 
		dc_nb9_23[lev] = XA[lev]*a + YA[lev]*alpha - B4*XA[lev]*b; 
		dc_nb9_24[lev] = XA[lev]*b; 

		dc_i9_32[lev] = XA[lev]*b;    
		dc_i9_33[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_i9_34[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_i9_35[lev] = XA[lev]*a + YA[lev]*alpha;
		dc_i9_36[lev] = XA[lev]*b;
		
		/*coeff_d<<dc_b9_11[lev]<<"\n"; 
		coeff_d<<dc_b9_12[lev]<<"\n"; 
		coeff_d<<dc_b9_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb9_21[lev]<<"\n"; 
		coeff_d<<dc_nb9_22[lev]<<"\n"; 
		coeff_d<<dc_nb9_23[lev]<<"\n"; 
		coeff_d<<dc_nb9_24[lev]<<"\n"; 
		
		coeff_d<<dc_i9_32[lev]<<"\n"; 
		coeff_d<<dc_i9_33[lev]<<"\n"; 
		coeff_d<<dc_i9_34[lev]<<"\n"; 
		coeff_d<<dc_i9_35[lev]<<"\n";
		coeff_d<<dc_i9_36[lev]<<"\n";*/

		//j 

		dc_b10_11[lev] = -2.*XH[lev]*a_nb + YH[lev] - B2*(XH[lev]*a_nb + YH[lev]*alpha_nb); 
		dc_b10_12[lev] = XH[lev]*a_nb + YH[lev]*alpha_nb - B3*(XH[lev]*a_nb + YH[lev]*alpha_nb); 
		dc_b10_13[lev] = -B4*( XH[lev]*a_nb + YH[lev]*alpha_nb ); 

		dc_nb10_21[lev] = XH[lev]*a + YH[lev]*alpha - B2*XH[lev]*b; 
		dc_nb10_22[lev] = -2.*(a+b)*XH[lev] + YH[lev] - B3*XH[lev]*b; 
		dc_nb10_23[lev] = XH[lev]*a + YH[lev]*alpha - B4*XH[lev]*b; 
		dc_nb10_24[lev] = XH[lev]*b; 

		dc_i10_32[lev] = XH[lev]*b;    
		dc_i10_33[lev] = XH[lev]*a + YH[lev]*alpha; 
		dc_i10_34[lev] = -2.*(a+b)*XH[lev] + YH[lev]; 
		dc_i10_35[lev] = XH[lev]*a + YH[lev]*alpha;
		dc_i10_36[lev] = XH[lev]*b;
		
		/*coeff_d<<dc_b10_11[lev]<<"\n"; 
		coeff_d<<dc_b10_12[lev]<<"\n"; 
		coeff_d<<dc_b10_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb10_21[lev]<<"\n"; 
		coeff_d<<dc_nb10_22[lev]<<"\n"; 
		coeff_d<<dc_nb10_23[lev]<<"\n"; 
		coeff_d<<dc_nb10_24[lev]<<"\n"; 
		
		coeff_d<<dc_i10_32[lev]<<"\n"; 
		coeff_d<<dc_i10_33[lev]<<"\n"; 
		coeff_d<<dc_i10_34[lev]<<"\n"; 
		coeff_d<<dc_i10_35[lev]<<"\n";
		coeff_d<<dc_i10_36[lev]<<"\n";*/

		//j+1 

		dc_b11_11[lev] = -2.*XA[lev]*a_nb + YA[lev] - B2*(XA[lev]*a_nb + YA[lev]*alpha_nb); 
		dc_b11_12[lev] = XA[lev]*a_nb + YA[lev]*alpha_nb - B3*(XA[lev]*a_nb + YA[lev]*alpha_nb); 
		dc_b11_13[lev] = -B4*( XA[lev]*a_nb + YA[lev]*alpha_nb ); 

		dc_nb11_21[lev] = XA[lev]*a + YA[lev]*alpha - B2*XA[lev]*b; 
		dc_nb11_22[lev] = -2.*(a+b)*XA[lev] + YA[lev] - B3*XA[lev]*b; 
		dc_nb11_23[lev] = XA[lev]*a + YA[lev]*alpha - B4*XA[lev]*b; 
		dc_nb11_24[lev] = XA[lev]*b; 

	        dc_i11_32[lev] = XA[lev]*b;    
		dc_i11_33[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_i11_34[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_i11_35[lev] = XA[lev]*a + YA[lev]*alpha;
		dc_i11_36[lev] = XA[lev]*b;
		
		/*coeff_d<<dc_b11_11[lev]<<"\n"; 
		coeff_d<<dc_b11_12[lev]<<"\n"; 
		coeff_d<<dc_b11_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb11_21[lev]<<"\n"; 
		coeff_d<<dc_nb11_22[lev]<<"\n"; 
		coeff_d<<dc_nb11_23[lev]<<"\n"; 
		coeff_d<<dc_nb11_24[lev]<<"\n"; 
		
		coeff_d<<dc_i11_32[lev]<<"\n"; 
		coeff_d<<dc_i11_33[lev]<<"\n"; 
		coeff_d<<dc_i11_34[lev]<<"\n"; 
		coeff_d<<dc_i11_35[lev]<<"\n";
		coeff_d<<dc_i11_36[lev]<<"\n";*/

		//j+2 

		dc_b12_11[lev] =  YB[lev] - B2*(YB[lev]*alpha_nb); 
		dc_b12_12[lev] =  YB[lev]*alpha_nb - B3*(YB[lev]*alpha_nb); 
		dc_b12_13[lev] = -B4*(YB[lev]*alpha_nb ); 

		dc_nb12_21[lev] = YB[lev]*alpha ; 
		dc_nb12_22[lev] = YB[lev]; 
		dc_nb12_23[lev] = YB[lev]*alpha; 
		dc_nb12_24[lev] = 0.0; 

		dc_i12_32[lev] = 0.0;    
		dc_i12_33[lev] = YB[lev]*alpha; 
		dc_i12_34[lev] = YB[lev]; 
		dc_i12_35[lev] = YB[lev]*alpha;
		dc_i12_36[lev] = 0.0;		
		
		/*coeff_d<<dc_b12_11[lev]<<"\n"; 
		coeff_d<<dc_b12_12[lev]<<"\n"; 
		coeff_d<<dc_b12_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb12_21[lev]<<"\n"; 
		coeff_d<<dc_nb12_22[lev]<<"\n"; 
		coeff_d<<dc_nb12_23[lev]<<"\n"; 
		coeff_d<<dc_nb12_24[lev]<<"\n"; 
		
		coeff_d<<dc_i12_32[lev]<<"\n"; 
		coeff_d<<dc_i12_33[lev]<<"\n"; 
		coeff_d<<dc_i12_34[lev]<<"\n"; 
		coeff_d<<dc_i12_35[lev]<<"\n";
		coeff_d<<dc_i12_36[lev]<<"\n";*/		
		
		/**********************************************************/
		//Coefficients for the RHS side 

		//j=1  and j=ny-1

		dc_b13_11[lev] = Z1[lev]*a_nb + Z2[lev]*alpha_nb; 
		dc_b13_12[lev] = -2.*Z1[lev]*a_nb + Z2[lev];
		dc_b13_13[lev] = Z1[lev]*a_nb + Z2[lev]*alpha_nb; 

		dc_nb13_21[lev] = Z1[lev]*b;  
		dc_nb13_22[lev] = Z1[lev]*a + Z2[lev]*alpha; 
		dc_nb13_23[lev] = -2.*(a+b)*Z1[lev] + Z2[lev]; 
		dc_nb13_24[lev] = Z1[lev]*a + Z2[lev]*alpha; 
		dc_nb13_25[lev] = Z1[lev]*b; 
		
		/*coeff_d<<dc_b13_11[lev]<<"\n"; 
		coeff_d<<dc_b13_12[lev]<<"\n"; 
		coeff_d<<dc_b13_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb13_21[lev]<<"\n"; 
		coeff_d<<dc_nb13_22[lev]<<"\n"; 
		coeff_d<<dc_nb13_23[lev]<<"\n"; 
		coeff_d<<dc_nb13_24[lev]<<"\n"; 
		coeff_d<<dc_nb13_25[lev]<<"\n"; */
		
		//j=2 and j=ny-2

		dc_b14_11[lev] = Y7[lev]*alpha_nb; 
		dc_b14_12[lev] = Y7[lev]; 
		dc_b14_13[lev] = Y7[lev]*alpha_nb; 

		dc_nb14_21[lev] = 0.0; 
		dc_nb14_22[lev] = Y7[lev]*alpha; 
		dc_nb14_23[lev] = Y7[lev]; 
		dc_nb14_24[lev] = Y7[lev]*alpha; 
		dc_nb14_25[lev] = 0.0; 						
		
		/*coeff_d<<dc_b14_11[lev]<<"\n"; 
		coeff_d<<dc_b14_12[lev]<<"\n"; 
		coeff_d<<dc_b14_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb14_21[lev]<<"\n"; 
		coeff_d<<dc_nb14_22[lev]<<"\n"; 
		coeff_d<<dc_nb14_23[lev]<<"\n"; 
		coeff_d<<dc_nb14_24[lev]<<"\n"; 
		coeff_d<<dc_nb14_25[lev]<<"\n";*/ 		
	}
	
		//coeff_c.close(); 
		//coeff_d.close();
}

/**************************************************************************************************/
/******************************Multigrid interpolating function************************************/
/**************************************************************************************************/

//The levels go like this  level[0]:the finest , level[mg_levels-1]: the coarsest 

void mg_interpolate(mg_grid * level,int lev,int fmg_bool)   //This is always acted upon from coarse to fine grids. So when the integer is passed it is an indication to pass the data onto the finer grid 
{
	//cout<<"In the interpolate function\n";  
	int ind; 

	int sp = pow(2,lev);     // Spacing for the level "lev"
	int sp1 = pow(2,lev-1);  
	
	int en_inx = nx-sp1, en_iny = ny-sp1; 

        /*Actual loop*/
        
        if(fmg_bool==0)
        {        
		for(int j=sp1;j<=en_iny;j=j+sp1)
		{
			for(int i=sp1;i<=en_inx;i=i+sp1)
			{
				ind = i + j*str_x; 

				if( (i%sp!=0)&&(j%sp!=0) ) level[lev-1].phi_s[ind] = level[lev-1].phi_s[ind] + 0.25*( level[lev].phi_s[ind+sp1] + level[lev].phi_s[ind-sp1] + level[lev].phi_s[ind+sp1*str_x] + level[lev].phi_s[ind-sp1*str_x] ); 
	
				if( (i%sp!=0)&&(j%sp==0) ) level[lev-1].phi_s[ind] = level[lev-1].phi_s[ind] + 0.5*( level[lev].phi_s[ind+sp1] + level[lev].phi_s[ind-sp1] ); 

				if( (i%sp==0)&&(j%sp!=0) ) level[lev-1].phi_s[ind] = level[lev-1].phi_s[ind] + 0.5*( level[lev].phi_s[ind+sp1*str_x] + level[lev].phi_s[ind-sp1*str_x] );

				if( (i%sp==0)&&(j%sp==0) ) level[lev-1].phi_s[ind] = level[lev-1].phi_s[ind] + level[lev].phi_s[ind];		
			}
		}	
	}
	else
	{
		for(int j=sp1;j<=en_iny;j=j+sp1)
		{
			for(int i=sp1;i<=en_inx;i=i+sp1)
			{
				ind = i + j*str_x; 

				if( (i%sp!=0)&&(j%sp!=0) ) level[lev-1].phi_s[ind] = 0.25*( level[lev].phi_s[ind+sp1] + level[lev].phi_s[ind-sp1] + level[lev].phi_s[ind+sp1*str_x] + level[lev].phi_s[ind-sp1*str_x] ); 

				if( (i%sp!=0)&&(j%sp==0) ) level[lev-1].phi_s[ind] = 0.5*( level[lev].phi_s[ind+sp1] + level[lev].phi_s[ind-sp1] ); 

				if( (i%sp==0)&&(j%sp!=0) ) level[lev-1].phi_s[ind] = 0.5*( level[lev].phi_s[ind+sp1*str_x] + level[lev].phi_s[ind-sp1*str_x] );

				if( (i%sp==0)&&(j%sp==0) ) level[lev-1].phi_s[ind] = level[lev].phi_s[ind];
			}
		}		
	}		
}

/**************************************************************************************************/
/********************************Multigrid Restricting function************************************/
//Restricts field values from fine grid to coarse grid 
void mg_restrict(mg_grid * level, int lev)
{
	//cout<<"Transfering data from grid level "<<lev<<" to "<<lev+1<<"\n"; 
	
	int sp = pow(2,lev); 
	int sp1 = pow(2,lev+1); //Restricting the field data on a finer grid to a coarser grid (higher level no.) 
	
	int en_inx = nx-sp, en_iny = ny-sp; 
	
	int ind; 
	
	//Full weighted operator for the interior cells of the coarse grid 
	
	for(int j=sp;j<=en_iny;j=j+sp)
	{
		for(int i=sp;i<=en_inx;i=i+sp)
		{
			ind = i+j*str_x; 
			
			if( (i%sp1==0) && (j%sp1==0) ) 
			{				
				level[lev+1].rhs[ind] = ( level[lev].res[ind-sp-sp*str_x] + level[lev].res[ind-sp+sp*str_x] + level[lev].res[ind+sp-sp*str_x] + level[lev].res[ind+sp+sp*str_x] + 2.*( level[lev].res[ind-sp*str_x] + level[lev].res[ind+sp*str_x] + level[lev].res[ind-sp] + level[lev].res[ind+sp] ) + 4.*( level[lev].res[ind] ) )*(1./16.) ; //Full weighted operator 
				
				//level[lev+1].rhs[ind] = level[lev].res[ind]; 	//Injection operator 		
				
				//level[lev+1].rhs[ind] = 0.0;  						
			}				
		}	
	}
	
	//Full weighting operator for testing purposes 
	
	/*for(int j=sp;j<=en_iny;j=j+sp)
	{
		for(int i=sp;i<=en_inx;i=i+sp)
		{
			ind = i+j*str_x; 
			
			if( (i%sp1==0) && (j%sp1==0) ) 
			{				
				level[lev+1].phi_s[ind] = ( level[lev].phi_s[ind-sp-sp*str_x] + level[lev].phi_s[ind-sp+sp*str_x] + level[lev].phi_s[ind+sp-sp*str_x] + level[lev].phi_s[ind+sp+sp*str_x] + 2.*( level[lev].phi_s[ind-sp*str_x] + level[lev].phi_s[ind+sp*str_x] + level[lev].phi_s[ind-sp] + level[lev].phi_s[ind+sp] ) + 4.*( level[lev].phi_s[ind] ) )*(1./16.) ;				 						
			}				
		}	
	}*/
					
}

/****************************************************************************/
/*Updating the pressures after relaxations are done in a multi-grid fashion*/
/****************************************************************************/
void mg_final(mg_grid * level, fval * fvar, int lev)  
{	
	for(int i=0;i<tot_p;i++)
	{
		fvar[i].u[2] = level[lev].phi_s[i]; 		
		fvar[i].F = level[lev].F[i]; 
		fvar[i].res = level[lev].res[i]; 
		fvar[i].err = level[lev].err[i]; //Error computed from residual of pressure from the formula Ae = r to compute the exact error
		fvar[i].rhs = level[lev].rhs[i]; 
	}
}

/********************************************************************/
/*Reading correction coefficients for various grids 
  from the null space vector of MATLAB*/ 
/********************************************************************/
void mg_read_coeff(mg_grid * level)
{
  int ind,st_inx,st_iny,step_x,step_y;   
  int en_inx, en_iny, sp; 
  ifstream coeff_input[mg_levels];
  ofstream coeff_check; 
  
  coeff_check.open("coeff_check.dat"); 
  
  string base_file("coeff_from_lab_"),file; 
  stringstream tag;             
      
  for(int lev=0;lev<mg_levels;lev++)
  {
     	//cout<<"Reading coefficients at level "<<lev<<"\n";
     
  	st_inx = pow(2,lev); 
        st_iny = pow(2,lev);
        
        sp = pow(2,lev); 
        
        en_inx = nx-sp; 
        en_iny = ny-sp; 
        
        step_x = pow(2,lev); 
  	step_y = pow(2,lev);
  	
  	tag<<lev+1;  	
  	file = base_file + tag.str() + ".dat";   
  		
	coeff_input[lev].open(file.c_str()); 
	
	for(int j=st_inx;j<=en_iny;j=j+step_x)
  	{
	    for(int i=st_iny;i<=en_inx;i=i+step_y)
    	    {
		ind = i + j*str_x;       
      		coeff_input[lev]>>level[lev].coeff[ind];		      		      		   		

      		if(lev==0) coeff_check<<i<<"	"<<j<<"		"<<level[lev].coeff[ind]<<"\n"; 
    	    }  
  	}
  	
  	tag.str("");
  	coeff_input[lev].close();   
  }
  
  	coeff_check.close(); 
}

/**************************************************************************************************************/
/***********************Evaluating RHS of the smoother for various levels of multigrid*************************/ 
/************************This is computed only at a particular level*******************************************/ 

void mg_evaluate_rhs(mg_grid * level, int lev, pbcs& pbc, int up)
{   
  vector<double> p1y(nx+1), p2y(nx+1), p1x(ny+1), p2x(ny+1); 
  
  int ind,j;
  
  int step_x = pow(2,lev);
  int step_y = pow(2,lev); 
  
  int sp=pow(2,lev); 
  
  int st_inx = pow(2,lev);
  int st_iny=pow(2,lev);
  
  int en_inx = nx-sp;
  int en_iny = ny-sp;   
  
  int ind1, ind2, ind3, ind4, ind5; 
  
  double extra=0.0; 
    
  for(int i=0;i<=nx;i++) //Getting the pressure gradient in the top and bottom faces of the domain
  {
     p1y[i] = pbc.p1y[i]; 
     p2y[i] = pbc.p2y[i]; 
  }
  
  for(int j=0;j<=ny;j++) //Getting the pressure gradient at the left and right faces of the domain 
  {
     p1x[j] = pbc.p1x[j]; 
     p2x[j] = pbc.p2x[j]; 
  } 
  
  /******************Assigning dp/dy at various j's*****************/
  
  /******************Case1 j=1****************/
  j=st_iny;
    
  for(int i=st_inx;i<=en_inx;i=i+sp)
  {
    ind = i + st_iny*str_x;  
    
    if(i==st_inx)
    {    
    	ind1 = sp + sp*str_x;  //Indices for extrapolation
    	ind2 = 2*ind1; 
    	ind3 = 3*ind1;
    	ind4 = 4*ind1; 
    	ind5 = 5*ind1; 
    
      //level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_b13_11[lev]*p1y[i-sp] + dc_b13_12[lev]*p1y[i] + dc_b13_13[lev]*p1y[i+sp] );     
      
      //level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_b13_12[lev]*p1y[i] + dc_b13_13[lev]*p1y[i+sp] ) - dc_b13_11[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );       
      
      //extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x] ) ;           
      
      extra = extra_pol(0, 0,level,lev);        
      
      level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_b13_12[lev]*p1y[i] + dc_b13_13[lev]*p1y[i+sp] ) - dc_b13_11[lev]*up*(extra); 
    }
    else if(i>=st_inx+sp && i<=en_inx-sp)
    {        
      //level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_nb13_21[lev]*p1y[i-2*sp] + dc_nb13_22[lev]*p1y[i-sp] + dc_nb13_23[lev]*p1y[i] + dc_nb13_24[lev]*p1y[i+sp] + dc_nb13_25[lev]*p1y[i+2*sp] );
      
      if(i==st_inx+sp)
      {
      	ind1 = sp + sp*str_x;  //Indices for extrapolation
    	ind2 = 2*ind1; 
    	ind3 = 3*ind1;
    	ind4 = 4*ind1; 
    	ind5 = 5*ind1;
      
      	//level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_nb13_22[lev]*p1y[i-sp] + dc_nb13_23[lev]*p1y[i] + dc_nb13_24[lev]*p1y[i+sp] + dc_nb13_25[lev]*p1y[i+2*sp] ) - dc_nb13_21[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );
      	
      	//extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x]); 
      	
  	extra = extra_pol(0, 0,level,lev);        
      	
      	level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_nb13_22[lev]*p1y[i-sp] + dc_nb13_23[lev]*p1y[i] + dc_nb13_24[lev]*p1y[i+sp] + dc_nb13_25[lev]*p1y[i+2*sp] ) - dc_nb13_21[lev]*up*(extra);       	
      }
      else if(i==en_inx-sp)
      {
      
      	ind1 = en_inx + sp*str_x; 
      	ind2 = 2*ind1; 
      	ind3 = 2*ind1; 
      	ind4 = 2*ind1; 
      	ind5 = 2*ind1; 
      
      	//level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_nb13_21[lev]*p1y[i-2*sp] + dc_nb13_22[lev]*p1y[i-sp] + dc_nb13_23[lev]*p1y[i] + dc_nb13_24[lev]*p1y[i+sp] ) - dc_nb13_25[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );
      	
      	//extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x] + level[lev].phi_s[ind1-sp-sp*str_x]);  //Mistake here
      	
     	extra = extra_pol(nx-3*sp, 0,level,lev);        
      	
      	level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_nb13_21[lev]*p1y[i-2*sp] + dc_nb13_22[lev]*p1y[i-sp] + dc_nb13_23[lev]*p1y[i] + dc_nb13_24[lev]*p1y[i+sp] ) - dc_nb13_25[lev]*up*(extra);       	
      	      	
      }
      else
      {
      	level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_nb13_21[lev]*p1y[i-2*sp] + dc_nb13_22[lev]*p1y[i-sp] + dc_nb13_23[lev]*p1y[i] + dc_nb13_24[lev]*p1y[i+sp] + dc_nb13_25[lev]*p1y[i+2*sp] );      
      }                                
    }
    else 
    {
    	ind1 = en_inx + sp*str_x; 
      	ind2 = 2*ind1; 
      	ind3 = 2*ind1; 
      	ind4 = 2*ind1; 
      	ind5 = 2*ind1;   	
    
        //level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_b13_11[lev]*p1y[i-sp] + dc_b13_12[lev]*p1y[i] + dc_b13_13[lev]*p1y[i+sp] );     
          	    
      	//level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_b13_11[lev]*p1y[i-sp] + dc_b13_12[lev]*p1y[i] ) - dc_b13_13[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5]  ) ;
      	
      	//extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x] + level[lev].phi_s[ind1-sp-sp*str_x]); //The earlier mistake propagated here
      	
     	extra = extra_pol(nx-3*sp, 0,level,lev);        
      	
      	level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_b13_11[lev]*p1y[i-sp] + dc_b13_12[lev]*p1y[i] ) - dc_b13_13[lev]*up*(extra);    	
    }   
    
     //mg_rhsy<<i<<"	"<<j<<"		"<<level[lev].cor_rhs[ind]<<"\n"; 
  }
  
  /***************Case2 j=2**************************/ 
  
  j = st_iny + sp; 
  
  for(int i=st_inx;i<=en_inx;i=i+sp) 
  {
    ind = i + j*str_x; 
    
    if(i==st_inx)
    {
      	ind1 = sp + sp*str_x;  //Indices for extrapolation
    	ind2 = 2*ind1; 
    	ind3 = 3*ind1;
    	ind4 = 4*ind1; 
    	ind5 = 5*ind1;  	
    
      	//level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_b14_11[lev]*p1y[i-sp] + dc_b14_12[lev]*p1y[i] + dc_b14_13[lev]*p1y[i+sp] );                
      	
      	//level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_b14_12[lev]*p1y[i] + dc_b14_13[lev]*p1y[i+sp] ) - dc_b14_11[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );      	
      	
	//extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x]); //Propagated mistake here
	
        extra = extra_pol(0, 0,level,lev);        
	
	level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_b14_12[lev]*p1y[i] + dc_b14_13[lev]*p1y[i+sp] ) - dc_b14_11[lev]*up*(extra);       	
    }
    else if(i>=st_inx+sp && i<=en_inx-sp)
    {
      //level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_nb14_21[lev]*p1y[i-2*sp] + dc_nb14_22[lev]*p1y[i-sp] + dc_nb14_23[lev]*p1y[i] + dc_nb14_24[lev]*p1y[i+sp] + dc_nb14_25[lev]*p1y[i+2*sp] );       
      
      if(i==st_inx+sp)
      {
      	ind1 = sp + sp*str_x;  //Indices for extrapolation
    	ind2 = 2*ind1; 
    	ind3 = 3*ind1;
    	ind4 = 4*ind1; 
    	ind5 = 5*ind1;    
    	     
      	//level[lev].cor_rhs[ind] = (-hy[lev]/B1)*(dc_nb14_22[lev]*p1y[i-sp] + dc_nb14_23[lev]*p1y[i] + dc_nb14_24[lev]*p1y[i+sp] + dc_nb14_25[lev]*p1y[i+2*sp] ) - dc_nb14_21[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );   
      	
      	//extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x] ); //Propagated mistake here
      	
        extra = extra_pol(0, 0,level,lev);        
      	
      	level[lev].cor_rhs[ind] = (-hy[lev]/B1)*(dc_nb14_22[lev]*p1y[i-sp] + dc_nb14_23[lev]*p1y[i] + dc_nb14_24[lev]*p1y[i+sp] + dc_nb14_25[lev]*p1y[i+2*sp] ) - dc_nb14_21[lev]*up*(extra);       	
      }
      else if(i==en_inx-sp)
      {
      	ind1 = en_inx + sp*str_x; 
      	ind2 = 2*ind1; 
      	ind3 = 2*ind1; 
      	ind4 = 2*ind1; 
      	ind5 = 2*ind1;       
      
      	//level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_nb14_21[lev]*p1y[i-2*sp] + dc_nb14_22[lev]*p1y[i-sp] + dc_nb14_23[lev]*p1y[i] + dc_nb14_24[lev]*p1y[i+sp] ) - dc_nb14_25[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );  
      	
      	//extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x] + level[lev].phi_s[ind1-sp-sp*str_x]); 
      	
        extra = extra_pol(nx-3*sp, 0,level,lev);        
      	
      	level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_nb14_21[lev]*p1y[i-2*sp] + dc_nb14_22[lev]*p1y[i-sp] + dc_nb14_23[lev]*p1y[i] + dc_nb14_24[lev]*p1y[i+sp] ) - dc_nb14_25[lev]*up*(extra); 	           }
      else      
      {            
      	level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_nb14_21[lev]*p1y[i-2*sp] + dc_nb14_22[lev]*p1y[i-sp] + dc_nb14_23[lev]*p1y[i] + dc_nb14_24[lev]*p1y[i+sp] + dc_nb14_25[lev]*p1y[i+2*sp] );           
      }           
    }
    else
    {
      //level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_b14_11[lev]*p1y[i-sp] + dc_b14_12[lev]*p1y[i] + dc_b14_13[lev]*p1y[i+sp] );                
      
      ind1 = en_inx + sp*str_x; 
      ind2 = 2*ind1; 
      ind3 = 2*ind1; 
      ind4 = 2*ind1; 
      ind5 = 2*ind1;
      
      //level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_b14_11[lev]*p1y[i-sp] + dc_b14_12[lev]*p1y[i] ) - dc_b14_13[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );            
      
      //extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x] + level[lev].phi_s[ind1-sp-sp*str_x] ); //Propagated mistake
      
      extra = extra_pol(nx-3*sp, 0,level,lev);        
      
      level[lev].cor_rhs[ind] = (-hy[lev]/B1)*( dc_b14_11[lev]*p1y[i-sp] + dc_b14_12[lev]*p1y[i] ) - dc_b14_13[lev]*up*(extra);       
    } 
    
    //mg_rhsy<<i<<"	"<<j<<"		"<<level[lev].cor_rhs[ind]<<"\n"; 
  }
  
  /***************Case3 j>2 and j<ny-2***************/ 
  
  for(int j=st_iny+2*sp;j<=en_iny-2*sp;j=j+sp)
  {
    for(int i=st_inx;i<=en_inx;i=i+sp)
    {
      ind = i + j*str_x; 
      
      level[lev].cor_rhs[ind] = 0.0;  //No Y-derivative effect here    
      
      //mg_rhsy<<i<<"	"<<j<<"		"<<level[lev].cor_rhs[ind]<<"\n"; 
    }  
  }  
  
  /***************Case4 j=ny-2***************/ 
  j = en_iny-sp; 
  
  for(int i=st_inx;i<=en_inx;i=i+sp) 
  {
    ind = i + j*str_x; 
    
    if(i==st_inx)
    {
        ind1 = sp + en_iny*str_x;  //Indices for extrapolation
    	ind2 = 2*ind1; 
    	ind3 = 3*ind1;
    	ind4 = 4*ind1; 
    	ind5 = 5*ind1;	
    
        //level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_b14_11[lev]*p2y[i-sp] + dc_b14_12[lev]*p2y[i] + dc_b14_13[lev]*p2y[i+sp] );                
         
        //level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_b14_12[lev]*p1y[i] + dc_b14_13[lev]*p1y[i+sp] ) - dc_b14_11[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );         
        
        //extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x]); //Propagated mistake here
        
        extra = extra_pol(0, ny-3*sp,level,lev);        
        
        level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_b14_12[lev]*p2y[i] + dc_b14_13[lev]*p2y[i+sp] ) - dc_b14_11[lev]*up*(extra);                   
    }
    else if(i>=st_inx+sp && i<=en_inx-sp)
    {
      //level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_nb14_21[lev]*p2y[i-2*sp] + dc_nb14_22[lev]*p2y[i-sp] + dc_nb14_23[lev]*p2y[i] + dc_nb14_24[lev]*p2y[i+sp] + dc_nb14_25[lev]*p2y[i+2*sp] ); 
            
      if(i==st_inx+sp)
      {
      	ind1 = sp + en_iny*str_x;  //Indices for extrapolation
    	ind2 = 2*ind1; 
    	ind3 = 3*ind1;
    	ind4 = 4*ind1; 
    	ind5 = 5*ind1;    
    	     
      	//level[lev].cor_rhs[ind] = (hy[lev]/B1)*(dc_nb14_22[lev]*p1y[i-sp] + dc_nb14_23[lev]*p1y[i] + dc_nb14_24[lev]*p1y[i+sp] + dc_nb14_25[lev]*p1y[i+2*sp] ) - dc_nb14_21[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );  
      	
      	//extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x] ); //Propagated mistake here
      	
	extra = extra_pol(0, ny-3*sp,level,lev);        
      	
      	level[lev].cor_rhs[ind] = (hy[lev]/B1)*(dc_nb14_22[lev]*p2y[i-sp] + dc_nb14_23[lev]*p2y[i] + dc_nb14_24[lev]*p2y[i+sp] + dc_nb14_25[lev]*p2y[i+2*sp] ) - dc_nb14_21[lev]*up*(extra);       	        	
      }
      else if(i==en_inx-sp)
      {
      	ind1 = en_inx + en_iny*str_x; 
      	ind2 = 2*ind1; 
      	ind3 = 2*ind1; 
      	ind4 = 2*ind1; 
      	ind5 = 2*ind1;      
      	
      	//level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_nb14_21[lev]*p1y[i-2*sp] + dc_nb14_22[lev]*p1y[i-sp] + dc_nb14_23[lev]*p1y[i] + dc_nb14_24[lev]*p1y[i+sp] ) - dc_nb14_25[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );  
      	
      	//extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x] + level[lev].phi_s[ind1-sp-sp*str_x]);
      	
	extra = extra_pol(nx-3*sp, ny-3*sp,level,lev);         
      	
      	level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_nb14_21[lev]*p2y[i-2*sp] + dc_nb14_22[lev]*p2y[i-sp] + dc_nb14_23[lev]*p2y[i] + dc_nb14_24[lev]*p2y[i+sp] ) - dc_nb14_25[lev]*up*(extra);      	   }
      else
      {       	
      	level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_nb14_21[lev]*p2y[i-2*sp] + dc_nb14_22[lev]*p2y[i-sp] + dc_nb14_23[lev]*p2y[i] + dc_nb14_24[lev]*p2y[i+sp] + dc_nb14_25[lev]*p2y[i+2*sp] );      
      }    
    }
    else
    {
      //level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_b14_11[lev]*p2y[i-sp] + dc_b14_12[lev]*p2y[i] + dc_b14_13[lev]*p2y[i+sp] );   
      
      ind1 = en_inx + en_iny*str_x; 
      ind2 = 2*ind1; 
      ind3 = 2*ind1; 
      ind4 = 2*ind1; 
      ind5 = 2*ind1;
      
      //level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_b14_11[lev]*p1y[i-sp] + dc_b14_12[lev]*p1y[i] ) - dc_b14_13[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );                
      
      //extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x] + level[lev].phi_s[ind1-sp-sp*str_x] ); 
      
      extra = extra_pol(nx-3*sp, ny-3*sp,level,lev);        
      
      level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_b14_11[lev]*p2y[i-sp] + dc_b14_12[lev]*p2y[i] ) - dc_b14_13[lev]*up*(extra); 
      
    }
    
    //mg_rhsy<<i<<"	"<<j<<"		"<<level[lev].cor_rhs[ind]<<"\n";  
  }
   
  /***************Case5 j=ny-1***************/   
  j = en_iny; 
  
  for(int i=st_inx;i<=en_inx;i=i+sp)
  {
    ind = i + j*str_x; 
    
    if(i==st_inx) 
    {
      //level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_b13_11[lev]*p2y[i-sp] + dc_b13_12[lev]*p2y[i] + dc_b13_13[lev]*p2y[i+sp] );     
      
        ind1 = sp + j*str_x;  //Indices for extrapolation
    	ind2 = 2*ind1; 
    	ind3 = 3*ind1;
    	ind4 = 4*ind1; 
    	ind5 = 5*ind1;           
      
        //level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_b13_12[lev]*p1y[i] + dc_b13_13[lev]*p1y[i+sp] ) - dc_b13_11[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );      
        
        //extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x]); //Propagated mistake here
        
        extra = extra_pol(0, ny-3*sp,level,lev);        
        
        level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_b13_12[lev]*p2y[i] + dc_b13_13[lev]*p2y[i+sp] ) - dc_b13_11[lev]*up*(extra);  
    }
    else if(i>=st_inx+sp && i<=en_inx-sp)
    {
      //level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_nb13_21[lev]*p2y[i-2*sp] + dc_nb13_22[lev]*p2y[i-sp] + dc_nb13_23[lev]*p2y[i] + dc_nb13_24[lev]*p2y[i+sp] + dc_nb13_25[lev]*p2y[i+2*sp] );    
      
      if(i==st_inx+sp)
      {
      	ind1 = sp + j*str_x;  //Indices for extrapolation
    	ind2 = 2*ind1; 
    	ind3 = 3*ind1;
    	ind4 = 4*ind1; 
    	ind5 = 5*ind1;
      
      	//level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_nb13_22[lev]*p1y[i-sp] + dc_nb13_23[lev]*p1y[i] + dc_nb13_24[lev]*p1y[i+sp] + dc_nb13_25[lev]*p1y[i+2*sp] ) - dc_nb13_21[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );
      	
      	//extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x]); //Propagated mistake here
      	
	extra = extra_pol(0, ny-3*sp, level, lev);        
      	
      	level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_nb13_22[lev]*p2y[i-sp] + dc_nb13_23[lev]*p2y[i] + dc_nb13_24[lev]*p2y[i+sp] + dc_nb13_25[lev]*p2y[i+2*sp] ) - dc_nb13_21[lev]*up*(extra);         
      }
      else if(i==en_inx-sp)
      {
      	ind1 = en_inx + j*str_x;       	
	ind2 = 2*ind1; 
      	ind3 = 2*ind1; 
      	ind4 = 2*ind1; 
      	ind5 = 2*ind1; 
      
      	//level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_nb13_21[lev]*p1y[i-2*sp] + dc_nb13_22[lev]*p1y[i-sp] + dc_nb13_23[lev]*p1y[i] + dc_nb13_24[lev]*p1y[i+sp] ) - dc_nb13_25[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5] );
      	
      	//extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x] + level[lev].phi_s[ind1-sp-sp*str_x] ); 
      	
	extra = extra_pol(nx-3*sp, ny-3*sp,level,lev);        
      	
      	level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_nb13_21[lev]*p2y[i-2*sp] + dc_nb13_22[lev]*p2y[i-sp] + dc_nb13_23[lev]*p2y[i] + dc_nb13_24[lev]*p2y[i+sp] ) - dc_nb13_25[lev]*up*( extra);    
      }
      else
      {
      	level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_nb13_21[lev]*p2y[i-2*sp] + dc_nb13_22[lev]*p2y[i-sp] + dc_nb13_23[lev]*p2y[i] + dc_nb13_24[lev]*p2y[i+sp] + dc_nb13_25[lev]*p2y[i+2*sp] );      
      }       
    }
    else 
    {
      	//level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_b13_11[lev]*p2y[i-sp] + dc_b13_12[lev]*p2y[i] + dc_b13_13[lev]*p2y[i+sp] );  
      
    	ind1 = en_inx + j*str_x; 
      	ind2 = 2*ind1; 
      	ind3 = 2*ind1; 
      	ind4 = 2*ind1; 
      	ind5 = 2*ind1;   	
              	    
      	//level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_b13_11[lev]*p1y[i-sp] + dc_b13_12[lev]*p1y[i] ) - dc_b13_13[lev]*up*( e1*level[lev].phi_s[ind1] + e2*level[lev].phi_s[ind2] + e3*level[lev].phi_s[ind3] + e4*level[lev].phi_s[ind4] + e5*level[lev].phi_s[ind5]  ) ;   
      	
	//extra = 4.0*( 3.*level[lev].phi_s[ind1] - 0.5*( level[lev].phi_s[ind1+sp*str_x] +  level[lev].phi_s[ind1-sp*str_x] +  level[lev].phi_s[ind1+sp] +  level[lev].phi_s[ind1-sp] ) ) - ( level[lev].phi_s[ind1-sp+sp*str_x] + level[lev].phi_s[ind1+sp-sp*str_x] + level[lev].phi_s[ind1-sp-sp*str_x]);
	
	extra = extra_pol(nx-3*sp, ny-3*sp,level,lev);        
	
	level[lev].cor_rhs[ind] = (hy[lev]/B1)*( dc_b13_11[lev]*p2y[i-sp] + dc_b13_12[lev]*p2y[i] ) - dc_b13_13[lev]*up*( extra );  	
    }    
    
    //mg_rhsy<<i<<"	"<<j<<"		"<<level[lev].cor_rhs[ind]<<"\n"; 
  }
  
  /*******************************************************************************************/  
  /*******************************Assigning dp/dx at various j's******************************/
 
  /*************Case1 j=st_iny******************/ 
  
  j = st_iny;
  
  for(int i=st_inx;i<=en_inx;i=i+sp)
  {
    ind = i + st_iny*str_x; 
    
    if(i==st_inx) 
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] - (X1[lev]*a_nb + X2[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j] - (X3[lev]*a_nb + X4[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j+sp] - (X5[lev]*a_nb + X6[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j+2*sp];     
    }
    
    if(i==st_inx+sp)
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] - X1[lev]*b*(hx[lev]/B1)*p1x[j] - X3[lev]*b*(hx[lev]/B1)*p1x[j+sp] - X5[lev]*b*(hx[lev]/B1)*p1x[j+2*sp]; 
    }
    
    if(i==en_inx-sp)
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] + X1[lev]*b*(hx[lev]/B1)*p2x[j] + X3[lev]*b*(hx[lev]/B1)*p2x[j+sp] + X5[lev]*b*(hx[lev]/B1)*p2x[j+2*sp]; 
    }
    
    if(i==en_inx)
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] + (X1[lev]*a_nb + X2[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j] + (X3[lev]*a_nb + X4[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j+sp] + (X5[lev]*a_nb + X6[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j+2*sp]  ;     
    } 
        
    //mg_rhsx<<i<<"	"<<j<<"		"<<level[lev].cor_rhs[ind]<<"\n"; 
  }
  
  /*************Case2 j=st_iny+sp******************/ 
  j = st_iny+sp; 
  
  for(int i=st_inx;i<=en_inx;i=i+sp)
  {
    ind = i + j*str_x; 
    
    if(i==st_inx) 
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] - (Y1[lev]*a_nb + Y2[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j-sp] - (Y3[lev]*a_nb + Y4[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j] - (Y5[lev]*a_nb + Y6[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j+sp] - Y7[lev]*alpha_nb*(hx[lev]/B1)*p1x[j+2*sp];     
    }
    
    if(i==st_inx+sp)
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] - Y1[lev]*b*(hx[lev]/B1)*p1x[j-sp] - Y3[lev]*b*(hx[lev]/B1)*p1x[j] - Y5[lev]*b*(hx[lev]/B1)*p1x[j+sp]; 
    }
    
    if(i==en_inx-sp)
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] + Y1[lev]*b*(hx[lev]/B1)*p2x[j-sp] + Y3[lev]*b*(hx[lev]/B1)*p2x[j] + Y5[lev]*b*(hx[lev]/B1)*p2x[j+sp]; 
    }
    
    if(i==en_inx)
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] + (Y1[lev]*a_nb + Y2[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j-sp] + (Y3[lev]*a_nb + Y4[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j] + (Y5[lev]*a_nb + Y6[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j+sp] + Y7[lev]*alpha_nb*(hx[lev]/B1)*p2x[j+2*sp] ;     
    }
    
    //mg_rhsx<<i<<"	"<<j<<"		"<<level[lev].cor_rhs[ind]<<"\n";  
  }
  
  /***********Case-3 for j>2 and j<ny-2*************/
  
  for(int j=st_iny+2*sp;j<=en_iny-2*sp;j=j+sp)
  {
    for(int i=st_inx;i<=en_inx;i=i+sp)
    {      
      ind = i + j*str_x; 
    
      if(i==st_inx) 
      {
	level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] - (YB[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j-2*sp] - (XA[lev]*a_nb + YA[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j-sp] - (XH[lev]*a_nb + YH[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j] - (XA[lev]*a_nb + YA[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j+sp] - (YB[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j+2*sp]  ;     
      }
    
      if(i==st_inx+sp)
      {
	level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] - XA[lev]*b*(hx[lev]/B1)*p1x[j-sp] - XH[lev]*b*(hx[lev]/B1)*p1x[j] - XA[lev]*b*(hx[lev]/B1)*p1x[j+sp]; 
      }
    
      if(i==en_inx-sp)
      {
	level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] + XA[lev]*b*(hx[lev]/B1)*p2x[j-sp] + XH[lev]*b*(hx[lev]/B1)*p2x[j] + XA[lev]*b*(hx[lev]/B1)*p2x[j+sp]; 
      }
    
      if(i==en_inx)
      {
	level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] + (YB[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j-2*sp] + (XA[lev]*a_nb + YA[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j-sp] + (XH[lev]*a_nb + YH[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j] + (XA[lev]*a_nb + YA[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j+sp] + (YB[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j+2*sp]  ;    
      }
      
      //mg_rhsx<<i<<"	"<<j<<"		"<<level[lev].cor_rhs[ind]<<"\n";        
    }
  }
  
  /***************Case4 j=ny-2***************/ 
  j=en_iny-sp; 
  
  for(int i=st_inx;i<=en_inx;i=i+sp)
  {
    ind = i + j*str_x; 
    
    if(i==st_inx) 
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] - (Y1[lev]*a_nb + Y2[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j+sp] - (Y3[lev]*a_nb + Y4[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j] - (Y5[lev]*a_nb + Y6[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j-sp] - Y7[lev]*alpha_nb*(hx[lev]/B1)*p1x[j-2*sp];     
    }
    
    if(i==st_inx+sp)
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] - Y1[lev]*b*(hx[lev]/B1)*p1x[j-sp] - Y3[lev]*b*(hx[lev]/B1)*p1x[j] - Y5[lev]*b*(hx[lev]/B1)*p1x[j+sp]; 
    }
    
    if(i==en_inx-sp)
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] + Y1[lev]*b*(hx[lev]/B1)*p2x[j-sp] + Y3[lev]*b*(hx[lev]/B1)*p2x[j] + Y5[lev]*b*(hx[lev]/B1)*p2x[j+sp]; 
      
      //cout<<"The cor_rhs value after correcting for dp/dx is "<<j<<"		"<<i<<"		"<<phi_s[ind].cor_rhs<<"\n"; 
    }
    
    if(i==en_inx)
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] + (Y1[lev]*a_nb + Y2[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j+sp] + (Y3[lev]*a_nb + Y4[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j] + (Y5[lev]*a_nb + Y6[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j-sp] + Y7[lev]*alpha_nb*(hx[lev]/B1)*p2x[j-2*sp];              
    }  
    
    //mg_rhsx<<i<<"	"<<j<<"		"<<level[lev].cor_rhs[ind]<<"\n";  
  }
  
  /***************Case5 j=ny-1***************/ 
  j = en_iny; 
  
  for(int i=st_inx;i<=en_inx;i=i+sp)
  {
    ind = i + j*str_x; 
    
    if(i==st_inx) 
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] - (X1[lev]*a_nb + X2[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j] - (X3[lev]*a_nb + X4[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j-sp] - (X5[lev]*a_nb + X6[lev]*alpha_nb)*(hx[lev]/B1)*p1x[j-2*sp]  ;     
    }
    
    if(i==st_inx+sp)
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] - X1[lev]*b*(hx[lev]/B1)*p1x[j] - X3[lev]*b*(hx[lev]/B1)*p1x[j-sp] - X5[lev]*b*(hx[lev]/B1)*p1x[j-2*sp]; 
    }
    
    if(i==en_inx-sp)
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] + X1[lev]*b*(hx[lev]/B1)*p2x[j] + X3[lev]*b*(hx[lev]/B1)*p2x[j-sp] + X5[lev]*b*(hx[lev]/B1)*p2x[j-2*sp]; 
    }
    
    if(i==en_inx)
    {
      level[lev].cor_rhs[ind] = level[lev].cor_rhs[ind] + (X1[lev]*a_nb + X2[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j] + (X3[lev]*a_nb + X4[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j-sp] + (X5[lev]*a_nb + X6[lev]*alpha_nb)*(hx[lev]/B1)*p2x[j-2*sp]  ;     
    }    
    
    //mg_rhsx<<i<<"	"<<j<<"		"<<level[lev].cor_rhs[ind]<<"\n";  
  }
  
  //mg_rhsx.close();    
  //mg_rhsy.close(); 
}

/****************************************************************************/
/***************Computing residual for just the interior grid****************/
void mg_residual_neu(mg_grid * level, int lev)
{
   int sp = pow(2,lev), ind; 
   
   int st_iny = sp, en_iny = ny-sp; 
   int st_inx = sp, en_inx = nx-sp; 
   
   double LHS, RHS;    	
   
    for(int j=st_iny;j<=en_iny;j=j+sp)
    {
      for(int i=st_inx;i<=en_inx;i=i+sp)
      {	        
	ind = i + j*str_x; 
	
       /*******************Case-1 j=1***********************/ 
	if(j==st_iny)
	{
	  if(i==st_inx)
	  {	        
	    RHS = dc_b1_12[lev]*level[lev].phi_s[ind+sp] + dc_b1_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b2_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b3_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_b3_13[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x]; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - level[lev].phi_s[ind]*dc_b1_11[lev]*level[lev].coeff[ind] ; 
	  }
	  else if(i==st_inx+sp)
	  {    
	    RHS = dc_nb1_21[lev]*level[lev].phi_s[ind-sp] + dc_nb1_23[lev]*level[lev].phi_s[ind+sp] + dc_nb1_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb2_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb3_21[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb3_24[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x];
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_nb1_22[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind]; 	  
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
	  {	    
	    RHS = dc_i1_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i1_33[lev]*level[lev].phi_s[ind-sp] + dc_i1_35[lev]*level[lev].phi_s[ind+sp] + dc_i1_36[lev]*level[lev].phi_s[ind+2*sp] + dc_i2_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i2_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i2_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i2_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i2_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_i3_32[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x] + dc_i3_33[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_i3_34[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_i3_35[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_i3_36[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x]; 
	    
	    level[lev].res[ind] =   level[lev].coeff[ind]*(level[lev].rhs[ind] -RHS) - level[lev].point_correc - dc_i1_34[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];
	  }
	  else if(i==en_inx-sp)
	  {	    	    
	    RHS = dc_nb1_21[lev]*level[lev].phi_s[ind+sp] + dc_nb1_23[lev]*level[lev].phi_s[ind-sp] + dc_nb1_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb2_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb3_21[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb3_24[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x];
	    
	    level[lev].res[ind] =  level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_nb1_22[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];	    
	  }
	  else 
	  {	    
	    RHS = dc_b1_12[lev]*level[lev].phi_s[ind-sp] + dc_b1_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b2_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b3_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_b3_13[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x]; 
	    
	    level[lev].res[ind]  = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_b1_11[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind] ;  	    	  
	  }	  	
	}	
	
	/*****************Case-2 (j=2)**********************/ 
	
	if(j==st_iny+sp)
	{
	  if(i==st_inx)
	  {	    
	    RHS = dc_b4_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b5_12[lev]*level[lev].phi_s[ind+sp] + dc_b5_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b6_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b7_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_b7_13[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x]; 
	    
	    level[lev].res[ind] =  level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_b5_11[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];   	  
	  }
	  else if(i==st_inx+sp) 
	  {	    
	    RHS = dc_nb4_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb5_21[lev]*level[lev].phi_s[ind-sp] + dc_nb5_23[lev]*level[lev].phi_s[ind+sp] + dc_nb5_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb6_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb7_21[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb7_24[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x]; 
	    
	    level[lev].res[ind]  =  level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - dc_nb5_22[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];   
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
	  {	    
	    RHS =  dc_i4_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i4_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i4_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i4_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i4_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_i5_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i5_33[lev]*level[lev].phi_s[ind-sp] + dc_i5_35[lev]*level[lev].phi_s[ind+sp] + dc_i5_36[lev]*level[lev].phi_s[ind+2*sp] + dc_i6_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i6_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i6_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i6_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i6_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_i7_32[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x] + dc_i7_33[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_i7_34[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_i7_35[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_i7_36[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x];
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_i5_34[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind]; 	  
	  }
	  else if(i==en_inx-sp)
	  {	    
	    RHS = dc_nb4_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb5_21[lev]*level[lev].phi_s[ind+sp] + dc_nb5_23[lev]*level[lev].phi_s[ind-sp] + dc_nb5_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb6_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb7_21[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb7_24[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x]; 
	    
	    level[lev].res[ind]  =  level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_nb5_22[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];	  
	  }
	  else 
	  {	    
	    RHS = dc_b4_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b5_12[lev]*level[lev].phi_s[ind-sp] + dc_b5_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b6_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b7_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_b7_13[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x]; 
	    
	    level[lev].res[ind]  = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_b5_11[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];	  
	  }  
	}
	
	/*****************Case-3 (j>=3 and j<=ny-3)**********************/ 
	
	if(j>=st_iny+2*sp && j<=en_iny-2*sp)
	{
	  if(i==st_inx)
	  {	    
	    RHS = dc_b8_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b8_12[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_b8_13[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x] + dc_b9_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b9_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b9_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b10_12[lev]*level[lev].phi_s[ind+sp] + dc_b10_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b11_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b11_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b11_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b12_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b12_12[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_b12_13[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x]; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - dc_b10_11[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind]; 	    
	  }
	  else if(i==st_inx+sp) 
	  {	    
	    RHS = dc_nb8_21[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb8_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb8_23[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb8_24[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x] + dc_nb9_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb9_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb9_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb9_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb10_21[lev]*level[lev].phi_s[ind-sp] + dc_nb10_23[lev]*level[lev].phi_s[ind+sp] + dc_nb10_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb11_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb11_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb11_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb11_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb12_21[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb12_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb12_23[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb12_24[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x]; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - dc_nb10_22[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind]; 
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
	  {    
	    RHS = dc_i8_33[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_i8_34[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_i8_35[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_i9_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i9_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i9_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i9_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i9_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_i10_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i10_33[lev]*level[lev].phi_s[ind-sp] +  dc_i10_35[lev]*level[lev].phi_s[ind+sp] + dc_i10_36[lev]*level[lev].phi_s[ind+2*sp] + dc_i11_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i11_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i11_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i11_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i11_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_i12_33[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_i12_34[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_i12_35[lev]*level[lev].phi_s[ind+sp+2*sp*str_x];
	    
	    level[lev].res[ind] =  level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_i10_34[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];  	    	   
	  }
	  else if(i==en_inx-sp)
	  {	    
	    RHS = dc_nb8_21[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb8_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb8_23[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb8_24[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x] + dc_nb9_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb9_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb9_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb9_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb10_21[lev]*level[lev].phi_s[ind+sp] + dc_nb10_23[lev]*level[lev].phi_s[ind-sp] + dc_nb10_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb11_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb11_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb11_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb11_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb12_21[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb12_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb12_23[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb12_24[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x]; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind]  - RHS) - level[lev].point_correc - dc_nb10_22[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind]; 	  
	  }
	  else 
	  {    
	    RHS = dc_b8_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b8_12[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_b8_13[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x] + dc_b9_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b9_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b9_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b10_12[lev]*level[lev].phi_s[ind-sp] + dc_b10_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b11_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b11_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b11_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b12_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b12_12[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_b12_13[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x]; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_b10_11[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];  
	  }	  	  
	}
	
	/*****************Case-4(j=ny-2)**********************/ 
	
	if(j==en_iny-sp)
	{
	  if(i==st_inx)
	  {	    
	    RHS = dc_b4_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b5_12[lev]*level[lev].phi_s[ind+sp] + dc_b5_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b6_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b7_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_b7_13[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x]; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_b5_11[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];   	  
	  }
	  else if(i==st_inx+sp) 
	  {	    
	    RHS = dc_nb4_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb5_21[lev]*level[lev].phi_s[ind-sp] + dc_nb5_23[lev]*level[lev].phi_s[ind+sp] + dc_nb5_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb6_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb7_21[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb7_24[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x]; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_nb5_22[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];   
	  }
	  else if(i>st_inx+sp && i<en_inx-sp)
	  {	    
	    RHS = dc_i4_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i4_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i4_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i4_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i4_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_i5_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i5_33[lev]*level[lev].phi_s[ind-sp] + dc_i5_35[lev]*level[lev].phi_s[ind+sp] + dc_i5_36[lev]*level[lev].phi_s[ind+2*sp] + dc_i6_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i6_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i6_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i6_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i6_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_i7_32[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x] + dc_i7_33[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_i7_34[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_i7_35[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_i7_36[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x]; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_i5_34[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind]; 	  
	  }
	  else if(i==en_inx-sp)
	  {	    
	    RHS = dc_nb4_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb5_21[lev]*level[lev].phi_s[ind+sp] + dc_nb5_23[lev]*level[lev].phi_s[ind-sp] + dc_nb5_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb6_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb7_21[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb7_24[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x]; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_nb5_22[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind] ;	  	    
	  }
	  else 
	  {	    
	    RHS = dc_b4_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b5_12[lev]*level[lev].phi_s[ind-sp] + dc_b5_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b6_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b7_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_b7_13[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x]; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_b5_11[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind] ;    	    
	  }  
	}
		
	/*****************Case-5(j=ny-1)**********************/ 
	
	if(j==en_iny)
	{
	  if(i==st_inx)
	  {	    
	    RHS = dc_b1_12[lev]*level[lev].phi_s[ind+sp] + dc_b1_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b2_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b3_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_b3_13[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x]; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - dc_b1_11[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];  	    
	  }
	  else if(i==st_inx+sp)
	  {	    
	    RHS = dc_nb1_21[lev]*level[lev].phi_s[ind-sp] + dc_nb1_23[lev]*level[lev].phi_s[ind+sp] + dc_nb1_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb2_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb3_21[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb3_24[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x] ; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - dc_nb1_22[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind] ; 	  
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
	  {	    
	    RHS = dc_i1_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i1_33[lev]*level[lev].phi_s[ind-sp] + dc_i1_35[lev]*level[lev].phi_s[ind+sp] + dc_i1_36[lev]*level[lev].phi_s[ind+2*sp] + dc_i2_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i2_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i2_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i2_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i2_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_i3_32[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x] + dc_i3_33[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_i3_34[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_i3_35[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_i3_36[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x] ; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - dc_i1_34[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];
	  }
	  else if(i==en_inx-sp)
	  {    
	    RHS = dc_nb1_21[lev]*level[lev].phi_s[ind+sp] + dc_nb1_23[lev]*level[lev].phi_s[ind-sp] + dc_nb1_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb2_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb3_21[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb3_24[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x]; 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - dc_nb1_22[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];	    
	  }
	  else 
	  {   	
	  
	     RHS = dc_b1_12[lev]*level[lev].phi_s[ind-sp] + dc_b1_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b2_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b3_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_b3_13[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x]; 
	    
	     level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - dc_b1_11[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind] ;	    	    	    		    		    
	  }	  	
	}
		
      }    
    }

}
/****************************************************************************/
/************Computing the total RHS sum for Neumann correction**************/
void mg_compute_rhs_tot(mg_grid * level, int lev)
{
     int ind; 
     
     int sp = pow(2,lev);      
     
     int st_inx = sp; 
     int en_inx = nx-sp;      
     
     int st_iny = sp; 
     int en_iny = ny-sp;   
     
     /*ofstream mg_rhs_tot; 
     mg_rhs_tot.open("mg_rhs_tot.dat");
     
     cout<<"In rhs to for level "<<lev<<"\n"; */
     
     for(int j=st_iny;j<=en_iny;j=j+sp)
     {
    	for(int i=st_inx;i<=en_inx;i=i+sp)
    	{
	      ind = i + j*str_x; 
	      
	      if(j==st_iny)
	      {
		if(i==st_inx)
		{
		  level[lev].rhs[ind] = alpha_nb*( alpha_nb*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha_nb*level[lev].F[ind+sp-sp*str_x] ) + ( alpha_nb*level[lev].F[ind-sp] + level[lev].F[ind] + alpha_nb*level[lev].F[ind+sp] ) + alpha_nb*( alpha_nb*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha_nb*level[lev].F[ind+sp+sp*str_x] ) + level[lev].cor_rhs[ind]; 
		}
		else if(i==st_inx+sp)
		{
		  level[lev].rhs[ind] = alpha_nb*( alpha*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x] ) + ( alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp] ) + alpha_nb*( alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x] ) + level[lev].cor_rhs[ind];  
		}
		else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
		{	  
		  level[lev].rhs[ind] = alpha_nb*( alpha*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x] ) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha_nb*( alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x] ) + level[lev].cor_rhs[ind]; 	
		}
		else if(i==en_inx-sp)
		{
		  level[lev].rhs[ind] = alpha_nb*( alpha*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x] ) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha_nb*( alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x] ) + level[lev].cor_rhs[ind];	
		}
		else 
		{
		  level[lev].rhs[ind] = alpha_nb*( alpha_nb*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha_nb*level[lev].F[ind+sp-sp*str_x] ) + ( alpha_nb*level[lev].F[ind-sp] + level[lev].F[ind] + alpha_nb*level[lev].F[ind+sp] ) + alpha_nb*( alpha_nb*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha_nb*level[lev].F[ind+sp+sp*str_x] ) + level[lev].cor_rhs[ind];	
		}	      
	      }      
	      else if(j==st_iny+sp) 
	      {
		if(i==st_inx)
		{
		  level[lev].rhs[ind] = alpha*(alpha_nb*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha_nb*level[lev].F[ind+sp-sp*str_x]) + ( alpha_nb*level[lev].F[ind-sp] + level[lev].F[ind] + alpha_nb*level[lev].F[ind+sp] ) + alpha*(alpha_nb*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha_nb*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind]; 	
		}
		else if(i==st_inx+sp)
		{
		  level[lev].rhs[ind] = alpha*(alpha*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x]) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha*(alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind]; 	
		}
		else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
		{
		  level[lev].rhs[ind] = alpha*(alpha*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x]) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha*(alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind]; 	
		}
		else if(i==en_inx-sp)
		{
		  level[lev].rhs[ind] = alpha*(alpha*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x]) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha*(alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind];	
		}
		else 
		{
		  level[lev].rhs[ind] = alpha*(alpha_nb*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha_nb*level[lev].F[ind+sp-sp*str_x]) + (alpha_nb*level[lev].F[ind-sp] + level[lev].F[ind] + alpha_nb*level[lev].F[ind+sp]) + alpha*(alpha_nb*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha_nb*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind]; 	  	
		}      
	      }
	      else if(j>=st_iny+2*sp && j<=en_iny-2*sp)
	      {
		if(i==st_inx)
		{
		  level[lev].rhs[ind] = alpha*(alpha_nb*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha_nb*level[lev].F[ind+sp-sp*str_x]) + ( alpha_nb*level[lev].F[ind-sp] + level[lev].F[ind] + alpha_nb*level[lev].F[ind+sp] ) + alpha*( alpha_nb*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha_nb*level[lev].F[ind+sp+sp*str_x] ) + level[lev].cor_rhs[ind]; 	
		}
		else if(i==st_inx+sp)
		{
		  level[lev].rhs[ind] = alpha*(alpha*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x]) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha*(alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind]; 	
		}
		else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
		{
		  level[lev].rhs[ind] = alpha*(alpha*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x]) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha*(alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind]; 	
		}
		else if(i==en_inx-sp)
		{
		  level[lev].rhs[ind] = alpha*(alpha*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x]) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha*(alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind]; 	
		}
		else 
		{
		  level[lev].rhs[ind] = alpha*(alpha_nb*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha_nb*level[lev].F[ind+sp-sp*str_x]) + (alpha_nb*level[lev].F[ind-sp] + level[lev].F[ind] + alpha_nb*level[lev].F[ind+sp]) + alpha*(alpha_nb*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha_nb*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind]; 
		}	      
	      }
	      else if(j==en_iny-sp)
	      {
		if(i==st_inx)
		{	  
		  level[lev].rhs[ind] = alpha*( alpha_nb*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha_nb*level[lev].F[ind+sp-sp*str_x] ) + (alpha_nb*level[lev].F[ind-sp] + level[lev].F[ind] + alpha_nb*level[lev].F[ind+sp]) + alpha*(alpha_nb*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha_nb*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind];		
		}
		else if(i==st_inx+sp)
		{
		  level[lev].rhs[ind] = alpha*(alpha*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x]) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha*(alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind]; 	
		}
		else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
		{
		  level[lev].rhs[ind] = alpha*(alpha*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x]) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha*(alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind];	
		}
		else if(i==en_inx-sp)
		{
		  level[lev].rhs[ind] = alpha*(alpha*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x]) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha*(alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind];
		}
		else 
		{
		  level[lev].rhs[ind] = alpha*(alpha_nb*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha_nb*level[lev].F[ind+sp-sp*str_x]) + (alpha_nb*level[lev].F[ind-sp] + level[lev].F[ind] + alpha_nb*level[lev].F[ind+sp]) + alpha*(alpha_nb*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha_nb*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind]; 	
		}     
	      }
	      else 
	      {
		if(i==st_inx)
		{
		 level[lev].rhs[ind] = alpha_nb*(alpha_nb*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha_nb*level[lev].F[ind+sp-sp*str_x]) + ( alpha_nb*level[lev].F[ind-sp] + level[lev].F[ind] + alpha_nb*level[lev].F[ind+sp] ) + alpha_nb*( alpha_nb*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha_nb*level[lev].F[ind+sp+sp*str_x] ) + level[lev].cor_rhs[ind]; 	
		}
		else if(i==st_inx+sp)
		{
		  level[lev].rhs[ind] = alpha_nb*(alpha*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x]) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha_nb*(alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind];	
		}
		else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
		{
		  level[lev].rhs[ind] = alpha_nb*(alpha*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x]) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha_nb*(alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind];	
		}
		else if(i==en_inx-sp)
		{
		  level[lev].rhs[ind] = alpha_nb*(alpha*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x]) + (alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp]) + alpha_nb*(alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind];	
		}
		else 
		{
		  level[lev].rhs[ind] = alpha_nb*(alpha_nb*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha_nb*level[lev].F[ind+sp-sp*str_x]) + (alpha_nb*level[lev].F[ind-sp] + level[lev].F[ind] + alpha_nb*level[lev].F[ind+sp]) + alpha_nb*(alpha_nb*level[lev].F[ind-sp+ sp*str_x] + level[lev].F[ind+sp*str_x] + alpha_nb*level[lev].F[ind+sp+sp*str_x]) + level[lev].cor_rhs[ind]; 	
		}      
	      }            
	    		    	
	    	//mg_rhs_tot<<i<<"	"<<j<<"		"<<level[lev].rhs[ind]<<"\n";  
	    } 
	  }
	  
	  //mg_rhs_tot.close(); 	        
} 

/************************************************************************/
//Gauss Seidel SOR smoothing at level "lev" and with no. of iterations specified as "ite_num"
/*void mg_gauss_seidel6(vector<mg_grid>& level, int lev, int ite_num, int up, int v_level, pbcs& pbc)
{
  int ind;       
  int nx_level,ny_level; 
  
  ofstream mg_gs; 
  mg_gs.open("mg_gs_imple_check.dat");
  
  double RHS, point_correc,res_norm,err_norm=1.0; 
  double rhs_sum=0.0,tol=1.0e-2; 
  double omega=1.0; //Relaxation factor 
  
  int i,j;
  int sp = pow(2,lev);      
  int st_inx = sp; 
  int en_inx = nx - sp;      
  int st_iny = sp; 
  int en_iny = ny - sp;  
   
  //point_correc=-0.0019159/((ny-1)*(nx-1));  //50X50  
  //point_correc = -0.277475/((ny-1)*(nx-1)); //500X250   
  //point_correc = -0.581476/((ny-1)*(nx-1)); //500X500   
  //point_correc = -0.275681/((ny-1)*(nx-1)); //499X249
  //point_correc = -0.000477791/((ny-1)*(nx-1)); //50X30
  //point_correc = -0.338275/((ny-1)*(nx-1)); //500X300  
  //point_correc = -0.0286014/((ny-1)*(nx-1)); //128X128  

   //Keying in the coefficients
   
  nx_level = (nx/pow(2,lev))+1; //No. of points in X-direction at each level 
  ny_level = (ny/pow(2,lev))+1; //No. of points in Y-direction at each level
 
  if(lev==v_level) //As both part of the ascending and descending parts of the cycle 
  {     
      mg_evaluate_rhs(level,lev,pbc,1);    
      mg_compute_rhs_tot(level,lev); 
      rhs_sum = mg_eval_rhs_sum(level,lev); 
      
      //cout<<"Smoothing at level "<<lev<<" and going up\n"; 
  }
  
    
  if( (lev!=v_level) && (up==0))
  {     	
      for(int j=sp;j<ny;j=j+sp)
      {
      	for(int i=sp;i<nx;i=i+sp)
      	{
      		ind = i + j*str_x;
      		
      		level[lev].phi_s[ind] = 0.0;       		      	 	      	
      	}      
      }
   
      rhs_sum = mg_eval_rhs_sum(level,lev);               
  }
   
   
  level[lev].point_correc = rhs_sum/( (nx_level-2)*(ny_level-2) ); 
    
  //cout<<"*********************************************"<<"\n"; 
  //cout<<"The RHS sum is "<<rhs_sum<<"\n"; 
  //cout<<"The point correc is "<<level[lev].point_correc<<"\n";  
  //cout<<"*********************************************"<<"\n";
    
  for(int ite=0;ite<ite_num;ite++)
  { 
    #pragma omp parallel for default(shared) private(i,j,ind,RHS) schedule(static)
    for(j=st_iny;j<=en_iny;j=j+sp)
    {
      for(i=st_inx;i<=en_inx;i=i+sp)
      {	        
	ind = i + j*str_x; 
	
	//Case-1 j=1 
	if(j==st_iny)
	{
	  if(i==st_inx)
	  {	        
	    RHS = dc_b1_12[lev]*level[lev].phi_s[ind+sp] + dc_b1_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b2_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b3_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_b3_13[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc )/(dc_b1_11[lev]*level[lev].coeff[ind]); 
	  }
	  else if(i==st_inx+sp)
	  {    
	    RHS = dc_nb1_21[lev]*level[lev].phi_s[ind-sp] + dc_nb1_23[lev]*level[lev].phi_s[ind+sp] + dc_nb1_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb2_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb3_21[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb3_24[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x];
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/(dc_nb1_22[lev]*level[lev].coeff[ind]); 	  
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
	  {	    
	    RHS = dc_i1_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i1_33[lev]*level[lev].phi_s[ind-sp] + dc_i1_35[lev]*level[lev].phi_s[ind+sp] + dc_i1_36[lev]*level[lev].phi_s[ind+2*sp] + dc_i2_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i2_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i2_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i2_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i2_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_i3_32[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x] + dc_i3_33[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_i3_34[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_i3_35[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_i3_36[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] =  (1.-omega)*level[lev].phi_s[ind] +  omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] -RHS) - level[lev].point_correc )/(dc_i1_34[lev]*level[lev].coeff[ind]);
	  }
	  else if(i==en_inx-sp)
	  {	    	    
	    RHS = dc_nb1_21[lev]*level[lev].phi_s[ind+sp] + dc_nb1_23[lev]*level[lev].phi_s[ind-sp] + dc_nb1_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb2_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb3_21[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb3_24[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x];
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/( dc_nb1_22[lev]*level[lev].coeff[ind] );	    
	  }
	  else 
	  {	    
	    RHS = dc_b1_12[lev]*level[lev].phi_s[ind-sp] + dc_b1_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b2_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b3_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_b3_13[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x]; 
	    
	    level[lev].phi_s[ind]  = (1.-omega)*level[lev].phi_s[ind] + omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc)/(dc_b1_11[lev]*level[lev].coeff[ind]) ;  	    	  
	  }	  	
	}	
	
	//Case-2 (j=2)
	
	if(j==st_iny+sp)
	{
	  if(i==st_inx)
	  {	    
	    RHS = dc_b4_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b5_12[lev]*level[lev].phi_s[ind+sp] + dc_b5_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b6_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b7_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_b7_13[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] =  (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/( dc_b5_11[lev]*level[lev].coeff[ind] );   	  
	  }
	  else if(i==st_inx+sp) 
	  {	    
	    RHS = dc_nb4_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb5_21[lev]*level[lev].phi_s[ind-sp] + dc_nb5_23[lev]*level[lev].phi_s[ind+sp] + dc_nb5_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb6_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb7_21[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb7_24[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x]; 
	    
	    level[lev].phi_s[ind]  =  (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc )/(dc_nb5_22[lev]*level[lev].coeff[ind]);   
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
	  {	    
	    RHS =  dc_i4_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i4_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i4_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i4_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i4_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_i5_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i5_33[lev]*level[lev].phi_s[ind-sp] + dc_i5_35[lev]*level[lev].phi_s[ind+sp] + dc_i5_36[lev]*level[lev].phi_s[ind+2*sp] + dc_i6_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i6_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i6_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i6_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i6_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_i7_32[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x] + dc_i7_33[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_i7_34[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_i7_35[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_i7_36[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x];
	    
	    level[lev].phi_s[ind] =  (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] -RHS) - level[lev].point_correc )/(dc_i5_34[lev]*level[lev].coeff[ind]); 	  
	  }
	  else if(i==en_inx-sp)
	  {	    
	    RHS = dc_nb4_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb5_21[lev]*level[lev].phi_s[ind+sp] + dc_nb5_23[lev]*level[lev].phi_s[ind-sp] + dc_nb5_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb6_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb7_21[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb7_24[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x]; 
	    
	    level[lev].phi_s[ind]  =  (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/(dc_nb5_22[lev]*level[lev].coeff[ind]);	  
	  }
	  else 
	  {	    
	    RHS = dc_b4_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b5_12[lev]*level[lev].phi_s[ind-sp] + dc_b5_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b6_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b7_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_b7_13[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x]; 
	    
	    level[lev].phi_s[ind]  = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/(dc_b5_11[lev]*level[lev].coeff[ind]);	  
	  }  
	}
	
	//Case-3 (j>=3 and j<=ny-3)
	
	if(j>=st_iny+2*sp && j<=en_iny-2*sp)
	{
	  if(i==st_inx)
	  {	    
	    RHS = dc_b8_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b8_12[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_b8_13[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x] + dc_b9_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b9_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b9_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b10_12[lev]*level[lev].phi_s[ind+sp] + dc_b10_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b11_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b11_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b11_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b12_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b12_12[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_b12_13[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc)/(dc_b10_11[lev]*level[lev].coeff[ind]); 	    
	  }
	  else if(i==st_inx+sp) 
	  {	    
	    RHS = dc_nb8_21[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb8_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb8_23[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb8_24[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x] + dc_nb9_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb9_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb9_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb9_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb10_21[lev]*level[lev].phi_s[ind-sp] + dc_nb10_23[lev]*level[lev].phi_s[ind+sp] + dc_nb10_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb11_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb11_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb11_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb11_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb12_21[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb12_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb12_23[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb12_24[lev]*level[lev].phi_s[ind+2*sp+2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc)/(dc_nb10_22[lev]*level[lev].coeff[ind]); 
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
	  {    
	    RHS = dc_i8_33[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_i8_34[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_i8_35[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_i9_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i9_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i9_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i9_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i9_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_i10_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i10_33[lev]*level[lev].phi_s[ind-sp] +  dc_i10_35[lev]*level[lev].phi_s[ind+sp] + dc_i10_36[lev]*level[lev].phi_s[ind+2*sp] + dc_i11_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i11_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i11_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i11_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i11_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_i12_33[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_i12_34[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_i12_35[lev]*level[lev].phi_s[ind+sp+2*sp*str_x];
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/ (dc_i10_34[lev]*level[lev].coeff[ind]);  	    	   
	  }
	  else if(i==en_inx-sp)
	  {	    
	    RHS = dc_nb8_21[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb8_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb8_23[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb8_24[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x] + dc_nb9_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb9_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb9_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb9_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb10_21[lev]*level[lev].phi_s[ind+sp] + dc_nb10_23[lev]*level[lev].phi_s[ind-sp] + dc_nb10_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb11_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb11_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb11_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb11_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb12_21[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb12_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb12_23[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb12_24[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind]  - RHS) - level[lev].point_correc )/(dc_nb10_22[lev]*level[lev].coeff[ind]); 	  
	  }
	  else 
	  {    
	    RHS = dc_b8_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b8_12[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_b8_13[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x] + dc_b9_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b9_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b9_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b10_12[lev]*level[lev].phi_s[ind-sp] + dc_b10_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b11_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b11_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b11_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b12_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b12_12[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_b12_13[lev]*level[lev].phi_s[ind-2*sp+2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/(dc_b10_11[lev]*level[lev].coeff[ind]);  
	  }	  	  
	}
	
	//Case-4(j=ny-2)
	
	if(j==en_iny-sp)
	{
	  if(i==st_inx)
	  {	    
	    RHS = dc_b4_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b5_12[lev]*level[lev].phi_s[ind+sp] + dc_b5_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b6_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b7_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_b7_13[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc)/(dc_b5_11[lev]*level[lev].coeff[ind]);   	  
	  }
	  else if(i==st_inx+sp) 
	  {	    
	    RHS = dc_nb4_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb5_21[lev]*level[lev].phi_s[ind-sp] + dc_nb5_23[lev]*level[lev].phi_s[ind+sp] + dc_nb5_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb6_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb7_21[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb7_24[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc)/( dc_nb5_22[lev]*level[lev].coeff[ind] );   
	  }
	  else if(i>st_inx+sp && i<en_inx-sp)
	  {	    
	    RHS = dc_i4_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i4_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i4_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i4_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i4_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_i5_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i5_33[lev]*level[lev].phi_s[ind-sp] + dc_i5_35[lev]*level[lev].phi_s[ind+sp] + dc_i5_36[lev]*level[lev].phi_s[ind+2*sp] + dc_i6_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i6_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i6_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i6_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i6_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_i7_32[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x] + dc_i7_33[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_i7_34[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_i7_35[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_i7_36[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/(dc_i5_34[lev]*level[lev].coeff[ind]); 	  
	  }
	  else if(i==en_inx-sp)
	  {	    
	    RHS = dc_nb4_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb5_21[lev]*level[lev].phi_s[ind+sp] + dc_nb5_23[lev]*level[lev].phi_s[ind-sp] + dc_nb5_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb6_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb7_21[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb7_24[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/( dc_nb5_22[lev]*level[lev].coeff[ind] );	  	    
	  }
	  else 
	  {	    
	    RHS = dc_b4_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b5_12[lev]*level[lev].phi_s[ind-sp] + dc_b5_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b6_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b7_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_b7_13[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/(dc_b5_11[lev]*level[lev].coeff[ind]);    	    
	  }  
	}
		
	//Case-5(j=ny-1)
	
	if(j==en_iny)
	{
	  if(i==st_inx)
	  {	    
	    RHS = dc_b1_12[lev]*level[lev].phi_s[ind+sp] + dc_b1_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b2_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b3_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_b3_13[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc )/(dc_b1_11[lev]*level[lev].coeff[ind]);  	    
	  }
	  else if(i==st_inx+sp)
	  {	    
	    RHS = dc_nb1_21[lev]*level[lev].phi_s[ind-sp] + dc_nb1_23[lev]*level[lev].phi_s[ind+sp] + dc_nb1_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb2_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb3_21[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb3_24[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x] ; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS) - level[lev].point_correc )/( dc_nb1_22[lev]*level[lev].coeff[ind] ); 	  
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
	  {	    
	    RHS = dc_i1_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i1_33[lev]*level[lev].phi_s[ind-sp] + dc_i1_35[lev]*level[lev].phi_s[ind+sp] + dc_i1_36[lev]*level[lev].phi_s[ind+2*sp] + dc_i2_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i2_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i2_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i2_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i2_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_i3_32[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x] + dc_i3_33[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_i3_34[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_i3_35[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_i3_36[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x] ; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc)/(dc_i1_34[lev]*level[lev].coeff[ind]);
	  }
	  else if(i==en_inx-sp)
	  {    
	    RHS = dc_nb1_21[lev]*level[lev].phi_s[ind+sp] + dc_nb1_23[lev]*level[lev].phi_s[ind-sp] + dc_nb1_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb2_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb3_21[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb3_24[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x]; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*(level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc)/(dc_nb1_22[lev]*level[lev].coeff[ind]);	    
	  }
	  else 
	  {	  
	     RHS = dc_b1_12[lev]*level[lev].phi_s[ind-sp] + dc_b1_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b2_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b3_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_b3_13[lev]*level[lev].phi_s[ind-2*sp-2*sp*str_x]; 
	    
	     level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc )/( dc_b1_11[lev]*level[lev].coeff[ind] );	 	    	    		    		    	     	          
	     
	  }	  	
	}
		//mg_gs<<i<<"	"<<j<<"		"<<level[lev].phi_s[ind]<<"\n";
		
		//mg_gs<<i<<"	"<<j<<"		"<<RHS<<"	"<<level[lev].F[ind]<<"\n"; 	
      }    
    }                    
       
  }    
  
  impose_zero_mean(level,lev);
  
  if(lev==0)
  {  	
  	mg_bcs_neu(level,lev,pbc); 
  }  
  //mg_gs.close();
}*/

/******************************************************************************/
/*******************Computing residual norm at level lev***********************/
double mg_res_norm(mg_grid * level, int lev)
{
	int sp = pow(2,lev); 
	int ind; 
	double abs_res = 0.0; 
	
	for(int j=sp;j<ny;j=j+sp)
	{
		for(int i=sp;i<nx;i=i+sp)
		{
			ind = i + j*str_x; 
			
			abs_res = abs_res + level[lev].res[ind]*level[lev].res[ind]; 
		
		}	
	}
		
	return sqrt(abs_res);  
}

/*******************************************************************************************/
//Multigrid Poisson Solver 
void mg_poisson_solver(mg_grid * level, fval * fvar, pbcs& pbc) 
{		
	//Assigning the field variable pressure to the multi-level variable at level-0 	
	double res_norm=0.0; 
	int num_v_c = 5, v_level=0; 

	/*ofstream res_out;	
	res_out.open("mg_v_residual.dat");*/  		
	
	mg_coeff();  	      // Coefficients for smoothing operator at all levels
	mg_read_coeff(level); // Reading the null space vectors for all grid sizes     	
		       	
	/******************************************/		       	
	/*for(int v_c=0;v_c<num_v_c;v_c++)  //The basic V-cycling 
	{	
		v_cycle(fvar,level,v_level,pbc); 		
		//res_norm = mg_res_norm(level,v_level); 		
		//cout<<v_c<<"		"<<log10(res_norm)<<"\n"; 			
	}*/
	/******************************************/
	
	fmg(fvar,level,pbc); 		
	mg_bcs_neu(level,v_level,pbc);
	mg_final(level,fvar,v_level);	
	
	//cout<<"The error norm is "<<error_calc(fvar,v_level)<<"\n";	
	//cout<<"Residual norm after FMG = "<<mg_res_norm(level,v_level)<<"\n"; 
 	//res_out.close(); 
}
/*******************************************************************************************/
//Full Multi Grid (FMG) Technique 
void fmg(fval * fvar, mg_grid * level, pbcs& pbc)
{
	int fmg_vc=1; // No.of v-cycles within a valley of FMG
	
	/***********************/
      	//Exact computation at the coarsest level	
	
	int coar_lev = mg_levels-1; //The level numbering starts from "0"
	int nx_level = (nx/pow(2,coar_lev)) + 1; //No. of points in X-direction at each level 
        int ny_level = (ny/pow(2,coar_lev)) + 1; //No. of points in Y-direction at each level 
	
	//start_poisson_source(fvar, level, coar_lev); 
	compute_pbc(fvar,pbc,1); 		
	mg_evaluate_rhs(level, coar_lev, pbc,1);    
        mg_compute_rhs_tot(level, coar_lev); 
        double rhs_sum = mg_eval_rhs_sum(level, coar_lev);
	
	level[coar_lev].point_correc = rhs_sum/( (nx_level-2)*(ny_level-2) ); 
	arma_direct_solve(level,coar_lev); 
	mg_interpolate(level,coar_lev,1); 
	
	/***********************/
	//The real FMG Loop
	int v_start = coar_lev-1; 
	
	for(int i=v_start;i>=0;i--)
	{
		for(int num_v=0;num_v<fmg_vc;num_v++)
		{
			v_cycle(fvar,level,i,pbc);		 		
		}
			
		if(i>0) mg_interpolate(level,i,1);	//"If" condition to make sure that the values are not interpolated from the finest grid
	}	
	
	/***The smallest unit that mimics the behavior of FMG***/ 
	
	/*int v_start = coar_lev-1; 
	int i = v_start; 
	v_cycle(fvar,level,i,pbc);
	compute_pbc(fvar,pbc,1); 		 		
	mg_bcs_neu(level,i,pbc); 
	if(i>0) mg_interpolate(level, i,1);	//"If" condition to make sure that the values are not interpolated from the finest grid
	v_cycle(fvar,level, i-1, pbc);*/
}	

/*******************************************************************************************/
//V-cycle required for Full Multi Grid (FMG)

void v_cycle(fval * fvar, mg_grid * level, int v_level, pbcs& pbc)
{
	int up=0,neu=1000; 
	double rhs_sum;
	int nx_level, ny_level; 
		
	//cout<<"No. of mg levels is "<<mg_levels<<"\n"; 
	//cout<<"v_level: "<<v_level<<"\n";
		
	for(int i=v_level;i<mg_levels;i++)	//Descending the stair case 
	{
		//cout<<"At level: "<<i<<"\n"; 
		
		if(i<mg_levels-1)
		{		
			if(i==v_level) 
			{
				compute_pbc(fvar,pbc,1); 
			}
			else
			{
				compute_pbc(fvar,pbc,0); 			
			}
			
			mg_bicgstab_corn(fvar, level, i, neu, up, v_level,pbc);
			mg_restrict(level,i); 				
		}
		else
		{			      
		      	nx_level = (nx/pow(2,i))+1; //No. of points in X-direction at each level 
			ny_level = (ny/pow(2,i))+1; //No. of points in Y-direction at each level
			
		      	rhs_sum = mg_eval_rhs_sum(level,i);
		      	level[i].point_correc = rhs_sum/( (nx_level-2)*(ny_level-2) ); 
		      	
		      	arma_direct_solve(level,i);
		      	compute_pbc(fvar,pbc,0);  	
		      	mg_bcs_neu(level, i, pbc); 		      			 					
		}
	}			
			
	up = 1;  //Variable indicating if the cycle is going up or down. 
	
	compute_pbc(fvar,pbc,up);
			
	for(int i=mg_levels-1; i>v_level;i--)  //Ascending the stair case 
	{
		//cout<<"At level: "<<i<<"\n";
		mg_interpolate(level,i,0);										
		mg_bicgstab_corn(fvar, level, i-1, 2*neu, up, v_level, pbc);
	}
	
	//cout<<"-------------------------------------------------\n";	
	
	//mg_residual_neu(level,v_level);
	//cout<<mg_res_norm(level,v_level)<<"\n";  	
}

/**********************************************************************/
void impose_zero_mean(mg_grid * level, int lev)
{
	int ind; 
	double sum=0.0, mean=0.0;
	int sp = pow(2,lev); 
	
	int spacin_x = nx/sp, spacin_y = ny/sp;
	int tot_inp = (spacin_x-1)*(spacin_y-1); //Total no. of interior points required for imposing zero mean  	
	
	for(int j=sp;j<ny;j=j+sp)
	{
	   for(int i=sp;i<nx;i=i+sp)
	   {
	   	ind = i + j*str_x; 
	   	sum = sum + level[lev].phi_s[ind]; 	   
	   }
	}	
	
	mean = sum/tot_inp; 
	
	//cout<<"The mean of the field is "<<mean<<"\n"; 
				
	for(int j=sp;j<ny;j=j+sp)
	{
	   for(int i=sp;i<nx;i=i+sp)
	   {	
	   	ind = i + j*str_x; 
	   	level[lev].phi_s[ind] = level[lev].phi_s[ind] - mean; 
	   }
	}	
}

/******************************************************************************************************/
/****Evaluating the RHS sum to further get the point correction for applying consistency condition*****/

double mg_eval_rhs_sum(mg_grid * level, int lev)
{
	int ind;
	
	double rhs_sum=0.0; 	
	int sp = pow(2,lev); 
	
	int st_inx = sp; 
	int en_inx = nx-sp; 
	int st_iny = sp; 
	int en_iny = ny-sp; 
	
	for(int j=st_iny;j<=en_iny;j=j+sp)
	{
		for(int i=st_inx;i<=en_inx;i=i+sp)
		{
			ind = i + j*str_x; 
			
			rhs_sum = rhs_sum + (level[lev].coeff[ind] * level[lev].rhs[ind]); 							
		}	
	}
	
	return rhs_sum; 
}

/********************************************************************/

void arma_direct_solve(mg_grid * level, int lev)
{
	int ind, ind_m, i_m,j_m; 
	
	int sp = pow(2,lev);
	int st_inx=sp, st_iny=sp;  
	int en_inx = nx-sp, en_iny = ny-sp; 
		
	int nx_sol = (nx/sp)-1, ny_sol = (ny/sp)-1; //No.of points that solved for directly. Only interior points are mentioned. 
	int tot_p_sol = nx_sol*ny_sol; 
	int str_m = nx_sol; 
	
	//cout<<"strm is "<<str_m<<"\n";
	//cout<<"total points are "<<tot_p_sol<<"\n"; 
	
	mat A = zeros<mat>(tot_p_sol,tot_p_sol);  	
	mat B = zeros<mat>(tot_p_sol,1);	
	mat trim_A = zeros<mat>(tot_p_sol-1,tot_p_sol-1);
	mat trim_B = zeros<mat>(tot_p_sol-1,1);	
	
	//vec uni = ones<vec>(tot_p_sol-1); 
		
	for(int j=sp;j<ny;j=j+sp)
	{
		for(int i=sp;i<nx;i=i+sp)
		{
			ind = i + j*str_x; 	
			
			i_m = (i/sp) - 1; 
			j_m = (j/sp) - 1; 		

			ind_m = i_m + j_m*str_m; 
			
			//cout<<"i_m= "<<i_m<<"j_m= "<<j_m<<"\n"; 
			//cout<<"Act-ind-m:   "<<ind_m<<"\n"; 
			
			/*The various cases would start now*/
			
			/********Case-1*******/
			
			/*******************Case-1 j=st_iny***********************/ 
			if(j==st_iny)
			{
			  if(i==st_inx)
			  { 			   			    
			    A(ind_m,ind_m) = dc_b1_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b1_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b1_13[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b2_11[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_b2_12[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_b2_13[lev]; 
			    
			    A(ind_m,ind_m+2*str_m) = dc_b3_11[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_b3_12[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_b3_13[lev]; 			    
			    
			  }
			  else if(i==st_inx+sp)
			  {     			    			    
			    A(ind_m,ind_m-1) = dc_nb1_21[lev]; 
			    A(ind_m,ind_m) = dc_nb1_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb1_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb1_24[lev]; 
			    
			    A(ind_m,ind_m-1+str_m) = dc_nb2_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb2_22[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_nb2_23[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_nb2_24[lev]; 
			    
			    A(ind_m,ind_m-1+2*str_m) = dc_nb3_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb3_22[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_nb3_23[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_nb3_24[lev]; 			    
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
			  {			    
			    A(ind_m,ind_m-2) = dc_i1_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i1_33[lev]; 
			    A(ind_m,ind_m) = dc_i1_34[lev];
			    A(ind_m,ind_m+1) = dc_i1_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i1_36[lev]; 
			    
			    A(ind_m,ind_m-2+str_m) = dc_i2_32[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_i2_33[lev]; 
			    A(ind_m,ind_m+str_m) = dc_i2_34[lev];
			    A(ind_m,ind_m+1+str_m) = dc_i2_35[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_i2_36[lev];
			    
			    A(ind_m,ind_m-2+2*str_m) = dc_i3_32[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_i3_33[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_i3_34[lev];
			    A(ind_m,ind_m+1+2*str_m) = dc_i3_35[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_i3_36[lev]; 			    
			  }
			  else if(i==en_inx-sp)
			  {    			     			    
			    A(ind_m,ind_m-2) = dc_nb1_24[lev];
    			    A(ind_m,ind_m-1) = dc_nb1_23[lev]; 
			    A(ind_m,ind_m) = dc_nb1_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb1_21[lev]; 
			    
			    A(ind_m,ind_m-2+str_m) = dc_nb2_24[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_nb2_23[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb2_22[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_nb2_21[lev];	    
			    
			    A(ind_m,ind_m-2+2*str_m) = dc_nb3_24[lev];	  
			    A(ind_m,ind_m-1+2*str_m) = dc_nb3_23[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb3_22[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_nb3_21[lev]; 		    
			  }
			  else 
			  {		    				    
			    A(ind_m,ind_m) = dc_b1_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b1_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b1_13[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b2_11[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_b2_12[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_b2_13[lev]; 
			    
			    A(ind_m,ind_m+2*str_m) = dc_b3_11[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_b3_12[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_b3_13[lev];		    
			  }	  	
			}
			
			/*******************Case-2 j=st_iny+sp***********************/ 
			
			if(j==st_iny+sp)
			{
			  if(i==st_inx)
			  {			    
			    A(ind_m,ind_m-str_m) = dc_b4_11[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_b4_12[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_b4_13[lev]; 
			    
			    A(ind_m,ind_m) = dc_b5_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b5_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b5_13[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b6_11[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_b6_12[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_b6_13[lev]; 
			    
			    A(ind_m,ind_m+2*str_m) = dc_b7_11[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_b7_12[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_b7_13[lev];			    
			  }
			  else if(i==st_inx+sp) 
			  {	    			    		    
			    A(ind_m,ind_m-1-str_m) = dc_nb4_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb4_22[lev]; 
			    A(ind_m,ind_m+1-str_m)= dc_nb4_23[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_nb4_24[lev]; 
			    
			    A(ind_m,ind_m-1) = dc_nb5_21[lev]; 
			    A(ind_m,ind_m) = dc_nb5_22[lev]; 
			    A(ind_m,ind_m+1)= dc_nb5_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb5_24[lev]; 
			    	    	    
			    A(ind_m,ind_m-1+str_m) = dc_nb6_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb6_22[lev]; 
			    A(ind_m,ind_m+1+str_m)= dc_nb6_23[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_nb6_24[lev];
			    
			    A(ind_m,ind_m-1+2*str_m) = dc_nb7_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb7_22[lev]; 
			    A(ind_m,ind_m+1+2*str_m)= dc_nb7_23[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_nb7_24[lev];			    
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
			  {			    
			    A(ind_m,ind_m-2-str_m) = dc_i4_32[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_i4_33[lev]; 
			    A(ind_m,ind_m-str_m) = dc_i4_34[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_i4_35[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_i4_36[lev];			    	  
			    
			    A(ind_m,ind_m-2) = dc_i5_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i5_33[lev]; 
			    A(ind_m,ind_m)   = dc_i5_34[lev]; 
			    A(ind_m,ind_m+1) = dc_i5_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i5_36[lev]	; 
			    
			    A(ind_m,ind_m-2+str_m) = dc_i6_32[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_i6_33[lev]; 
			    A(ind_m,ind_m+str_m)   = dc_i6_34[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_i6_35[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_i6_36[lev]; 
			    
			    A(ind_m,ind_m-2+2*str_m) = dc_i7_32[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_i7_33[lev]; 
			    A(ind_m,ind_m+2*str_m)   = dc_i7_34[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_i7_35[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_i7_36[lev];   	     
			  }
			  else if(i==en_inx-sp)
			  {	    			       
			    A(ind_m,ind_m+1-str_m) = dc_nb4_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb4_22[lev]; 
			    A(ind_m,ind_m-1-str_m)= dc_nb4_23[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_nb4_24[lev]; 
			    
			    A(ind_m,ind_m+1) = dc_nb5_21[lev]; 
			    A(ind_m,ind_m) = dc_nb5_22[lev]; 
			    A(ind_m,ind_m-1)= dc_nb5_23[lev]; 
			    A(ind_m,ind_m-2) = dc_nb5_24[lev]; 
			    	    	    
			    A(ind_m,ind_m+1+str_m) = dc_nb6_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb6_22[lev]; 
			    A(ind_m,ind_m-1+str_m)= dc_nb6_23[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_nb6_24[lev];
			    
			    A(ind_m,ind_m+1+2*str_m) = dc_nb7_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb7_22[lev]; 
			    A(ind_m,ind_m-1+2*str_m)= dc_nb7_23[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_nb7_24[lev];		      
			  }
			  else 
			  {	    
			    A(ind_m,ind_m-str_m) = dc_b4_11[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_b4_12[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_b4_13[lev]; 
			    
			    A(ind_m,ind_m) = dc_b5_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b5_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b5_13[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b6_11[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_b6_12[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_b6_13[lev]; 
			    
			    A(ind_m,ind_m+2*str_m) = dc_b7_11[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_b7_12[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_b7_13[lev];	       	  
			  }  
			}
			
			/*******************Case-3 j>=st_iny+2*sp && j<=en_iny-2*sp*************/
			
			if(j>=st_iny+2*sp && j<=en_iny-2*sp)
			{
			  if(i==st_inx)
			  {	    			    	    
			    A(ind_m,ind_m-2*str_m) = dc_b8_11[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_b8_12[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_b8_13[lev];
			    
			    A(ind_m,ind_m-str_m) = dc_b9_11[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_b9_12[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_b9_13[lev];
			    
			    A(ind_m,ind_m) = dc_b10_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b10_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b10_13[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b11_11[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_b11_12[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_b11_13[lev]; 
			    
			    A(ind_m,ind_m+2*str_m) = dc_b12_11[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_b12_12[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_b12_13[lev];  
			  }
			  else if(i==st_inx+sp) 
			  {	    			    
			    
			    A(ind_m,ind_m-1-2*str_m) = dc_nb8_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb8_22[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_nb8_23[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_nb8_24[lev]; 
			    
			    A(ind_m,ind_m-1-str_m) = dc_nb9_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb9_22[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_nb9_23[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_nb9_24[lev];	
			    
			    A(ind_m,ind_m-1) = dc_nb10_21[lev]; 
			    A(ind_m,ind_m) = dc_nb10_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb10_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb10_24[lev];
			    
			    A(ind_m,ind_m-1+str_m) = dc_nb11_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb11_22[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_nb11_23[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_nb11_24[lev];
			    
			    A(ind_m,ind_m-1+2*str_m) = dc_nb12_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb12_22[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_nb12_23[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_nb12_24[lev];   		    
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
			  {    			   			    
			    A(ind_m,ind_m-1-2*str_m) = dc_i8_33[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_i8_34[lev];
			    A(ind_m,ind_m+1-2*str_m) = dc_i8_35[lev];
			    
			    A(ind_m,ind_m-2-str_m) = dc_i9_32[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_i9_33[lev]; 
			    A(ind_m,ind_m-str_m) = dc_i9_34[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_i9_35[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_i9_36[lev];
			    
			    A(ind_m,ind_m-2) = dc_i10_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i10_33[lev]; 
			    A(ind_m,ind_m) = dc_i10_34[lev]; 
			    A(ind_m,ind_m+1) = dc_i10_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i10_36[lev];
			    
			    A(ind_m,ind_m-2+str_m) = dc_i11_32[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_i11_33[lev]; 
			    A(ind_m,ind_m+str_m) = dc_i11_34[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_i11_35[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_i11_36[lev];
			    
			    A(ind_m,ind_m-1+2*str_m) = dc_i12_33[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_i12_34[lev];
			    A(ind_m,ind_m+1+2*str_m) = dc_i12_35[lev];			    
			  }
			  else if(i==en_inx-sp)
			  {	    			     	  			    
			    A(ind_m,ind_m+1-2*str_m) = dc_nb8_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb8_22[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_nb8_23[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_nb8_24[lev]; 
			    
			    A(ind_m,ind_m+1-str_m) = dc_nb9_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb9_22[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_nb9_23[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_nb9_24[lev];	
			    
			    A(ind_m,ind_m+1) = dc_nb10_21[lev]; 
			    A(ind_m,ind_m) = dc_nb10_22[lev]; 
			    A(ind_m,ind_m-1) = dc_nb10_23[lev]; 
			    A(ind_m,ind_m-2) = dc_nb10_24[lev];
			    
			    A(ind_m,ind_m+1+str_m) = dc_nb11_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb11_22[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_nb11_23[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_nb11_24[lev];
			    
			    A(ind_m,ind_m+1+2*str_m) = dc_nb12_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb12_22[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_nb12_23[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_nb12_24[lev];			    
			  }
			  else 
			  {    			    		    
			    A(ind_m,ind_m-2*str_m) = dc_b8_11[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_b8_12[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_b8_13[lev];
			    
			    A(ind_m,ind_m-str_m) = dc_b9_11[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_b9_12[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_b9_13[lev];
			    
			    A(ind_m,ind_m) = dc_b10_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b10_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b10_13[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b11_11[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_b11_12[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_b11_13[lev]; 
			    
			    A(ind_m,ind_m+2*str_m) = dc_b12_11[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_b12_12[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_b12_13[lev];		      
			  }	  	  
			}	
					
			/************************Case-4 j==en_iny-sp*************************************/
			if(j==en_iny-sp)
			{
			  if(i==st_inx)
			  {	    			    			    
			    A(ind_m,ind_m+str_m) = dc_b4_11[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_b4_12[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_b4_13[lev]; 
			    
			    A(ind_m,ind_m) = dc_b5_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b5_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b5_13[lev]; 
			    
			    A(ind_m,ind_m-str_m) = dc_b6_11[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_b6_12[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_b6_13[lev]; 
			    
			    A(ind_m,ind_m-2*str_m) = dc_b7_11[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_b7_12[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_b7_13[lev];			      
			  }
			  else if(i==st_inx+sp) 
			  {			       			    
			    A(ind_m,ind_m-1+str_m) = dc_nb4_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb4_22[lev]; 
			    A(ind_m,ind_m+1+str_m)= dc_nb4_23[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_nb4_24[lev]; 
			    
			    A(ind_m,ind_m-1) = dc_nb5_21[lev]; 
			    A(ind_m,ind_m) = dc_nb5_22[lev]; 
			    A(ind_m,ind_m+1)= dc_nb5_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb5_24[lev]; 
			    	    	    
			    A(ind_m,ind_m-1-str_m) = dc_nb6_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb6_22[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_nb6_23[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_nb6_24[lev];
			    
			    A(ind_m,ind_m-1-2*str_m) = dc_nb7_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb7_22[lev]; 
			    A(ind_m,ind_m+1-2*str_m)= dc_nb7_23[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_nb7_24[lev];   	    
			  }
			  else if(i>st_inx+sp && i<en_inx-sp)
			  {  	  			    
			    A(ind_m,ind_m-2+str_m) = dc_i4_32[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_i4_33[lev]; 
			    A(ind_m,ind_m+str_m) = dc_i4_34[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_i4_35[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_i4_36[lev];			    	  
			    
			    A(ind_m,ind_m-2) = dc_i5_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i5_33[lev]; 
			    A(ind_m,ind_m)   = dc_i5_34[lev]; 
			    A(ind_m,ind_m+1) = dc_i5_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i5_36[lev]; 
			    
			    A(ind_m,ind_m-2-str_m) = dc_i6_32[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_i6_33[lev]; 
			    A(ind_m,ind_m-str_m)   = dc_i6_34[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_i6_35[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_i6_36[lev]; 
			    
			    A(ind_m,ind_m-2-2*str_m) = dc_i7_32[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_i7_33[lev]; 
			    A(ind_m,ind_m-2*str_m)   = dc_i7_34[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_i7_35[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_i7_36[lev];			    
			  }
			  else if(i==en_inx-sp)
			  {	    			    	  	    			    			    
			    A(ind_m,ind_m+1+str_m) = dc_nb4_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb4_22[lev]; 
			    A(ind_m,ind_m-1+str_m)= dc_nb4_23[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_nb4_24[lev]; 
			    
			    A(ind_m,ind_m+1) = dc_nb5_21[lev]; 
			    A(ind_m,ind_m) = dc_nb5_22[lev]; 
			    A(ind_m,ind_m-1)= dc_nb5_23[lev]; 
			    A(ind_m,ind_m-2) = dc_nb5_24[lev]; 
			    	    	    
			    A(ind_m,ind_m+1-str_m) = dc_nb6_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb6_22[lev]; 
			    A(ind_m,ind_m-1-str_m)= dc_nb6_23[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_nb6_24[lev];
			    
			    A(ind_m,ind_m+1-2*str_m) = dc_nb7_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb7_22[lev]; 
			    A(ind_m,ind_m-1-2*str_m)= dc_nb7_23[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_nb7_24[lev];			    
			  }
			  else 
			  {    			      	    			    
			    A(ind_m,ind_m+str_m) = dc_b4_11[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_b4_12[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_b4_13[lev]; 
			    
			    A(ind_m,ind_m) = dc_b5_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b5_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b5_13[lev]; 
			    
			    A(ind_m,ind_m-str_m) = dc_b6_11[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_b6_12[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_b6_13[lev]; 
			    
			    A(ind_m,ind_m-2*str_m) = dc_b7_11[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_b7_12[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_b7_13[lev];		    
			  }  
			}			
			
			/******************Case-5(j=en_iny)**********************/ 
	
			if(j==en_iny)
			{
			  if(i==st_inx)
			  {	    			        			    
			    A(ind_m,ind_m) = dc_b1_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b1_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b1_13[lev]; 
			    
			    A(ind_m,ind_m-str_m) = dc_b2_11[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_b2_12[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_b2_13[lev]; 
			    
			    A(ind_m,ind_m-2*str_m) = dc_b3_11[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_b3_12[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_b3_13[lev];			    
			  }
			  else if(i==st_inx+sp)
			  {	    			    			    
			    A(ind_m,ind_m-1) = dc_nb1_21[lev]; 
			    A(ind_m,ind_m) = dc_nb1_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb1_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb1_24[lev]; 
			    
			    A(ind_m,ind_m-1-str_m) = dc_nb2_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb2_22[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_nb2_23[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_nb2_24[lev]; 
			    
			    A(ind_m,ind_m-1-2*str_m) = dc_nb3_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb3_22[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_nb3_23[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_nb3_24[lev];			    
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
			  {	    			    			    			    
			    A(ind_m,ind_m-2) = dc_i1_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i1_33[lev]; 
			    A(ind_m,ind_m) = dc_i1_34[lev];
			    A(ind_m,ind_m+1) = dc_i1_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i1_36[lev]; 
			    
			    A(ind_m,ind_m-2-str_m) = dc_i2_32[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_i2_33[lev]; 
			    A(ind_m,ind_m-str_m) = dc_i2_34[lev];
			    A(ind_m,ind_m+1-str_m) = dc_i2_35[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_i2_36[lev];
			    
			    A(ind_m,ind_m-2-2*str_m) = dc_i3_32[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_i3_33[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_i3_34[lev];
			    A(ind_m,ind_m+1-2*str_m) = dc_i3_35[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_i3_36[lev];
			  }
			  else if(i==en_inx-sp)
			  {    			    		    
			    A(ind_m,ind_m-2) = dc_nb1_24[lev];
    			    A(ind_m,ind_m-1) = dc_nb1_23[lev]; 
			    A(ind_m,ind_m) = dc_nb1_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb1_21[lev]; 
			    
			    A(ind_m,ind_m-2-str_m) = dc_nb2_24[lev];		    
			    A(ind_m,ind_m-1-str_m) = dc_nb2_23[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb2_22[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_nb2_21[lev];	    
			    
			    A(ind_m,ind_m-2-2*str_m) = dc_nb3_24[lev];		  
			    A(ind_m,ind_m-1-2*str_m) = dc_nb3_23[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb3_22[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_nb3_21[lev];			        
			  }
			  else 
			  {   				  	       	    	    		    		    			     
			    A(ind_m,ind_m) = dc_b1_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b1_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b1_13[lev]; 
			    
			    A(ind_m,ind_m-str_m) = dc_b2_11[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_b2_12[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_b2_13[lev]; 
			    
			    A(ind_m,ind_m-2*str_m) = dc_b3_11[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_b3_12[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_b3_13[lev]; 
			     
			  }	  	
			}		
			
			/**********************************************/
			//Populating the full rhs matrix 
			
			B[ind_m] = level[lev].rhs[ind] - (level[lev].point_correc/level[lev].coeff[ind]);  							
				
		}	
	}
	
	//cout<<"The determinant of the matrix is   "<<det(A)<<"\n";
	
	//trim_A.save("trim_A.dat", raw_ascii);
	
	/**********************************************/
	//Reducing the matrix to a non-singular matrix 
	
	trim_A = A.submat(0,0,tot_p_sol-2,tot_p_sol-2); 
	trim_B = B.submat(0,0,tot_p_sol-2,0);
	
	mat X = solve(trim_A,trim_B);	
	
	//X.save("X.dat",raw_ascii); 
	
	/*vec prodvec = trim_A*uni; 
	
	prodvec.save("prodvec.dat",raw_ascii); */
	
	for(int j=sp;j<ny;j=j+sp)
	{
		for(int i=sp;i<nx;i=i+sp)
		{
			i_m = (i/sp) - 1; 
			j_m = (j/sp) - 1; 
			
			ind = i + j*str_x; 	
			ind_m = i_m + j_m*str_m; 
			
			if(ind_m<=tot_p_sol-2) 
			{
				
				level[lev].phi_s[ind] = X(ind_m);	
			}
			else
			{
				level[lev].phi_s[ind] = 0.0; 			
			}
		}
	}	
		
}

/*************************************************************/ 
//Functions linking up Momentum and Pressure Poisson's equation 

void poisson_source(fval * fvar, mg_grid * level, double Q)
{
	int ind;

        tdma1x(fvar,0); //du/dx 
        tdma1y(fvar,0); //du/dy  

        tdma1x(fvar,1); //dv/dx 
        tdma1y(fvar,1); //dv/dy 

        tdma2x(fvar,0);//d2u/dx2
        tdma2y(fvar,0);//d2u/dy2 

        tdma2x(fvar,1);//d2v/dx2
        tdma2y(fvar,1);//d2v/dy2 

        for(int j=1;j<ny;j++)
        {
                for(int i=1;i<nx;i++)
                {
                        ind = i + j*str_x;

                        fvar[ind].u[3] = -( fvar[ind].u[0]*fvar[ind].ux[0] + fvar[ind].u[1]*fvar[ind].uy[0] ) +  (1./Re)*( fvar[ind].uxx[0] + fvar[ind].uyy[0] );

                        fvar[ind].u[4] = -( fvar[ind].u[0]*fvar[ind].ux[1] + fvar[ind].u[1]*fvar[ind].uy[1] ) +  (1./Re)*( fvar[ind].uxx[1] + fvar[ind].uyy[1] );
                }
        }

        tdma1x(fvar,3);
        tdma1y(fvar,4);

        for(int j=1;j<ny;j++)
        {
                for(int i=1;i<nx;i++)
                {
                        ind = i + j*str_x;

                        fvar[ind].F = fvar[ind].ux[3] + fvar[ind].uy[4] + (Q/dt)*( fvar[ind].ux[0] + fvar[ind].uy[1] );
                        level[0].F[ind] = fvar[ind].F;                         
                }
        }
        
       //Injecting Poisson's source at all levels 
       
       for(int lev=0;lev<mg_levels;lev++)
       {
       		for(int j=1;j<ny;j++)
        	{
        	        for(int i=1;i<nx;i++)
        	        {
        			ind = i + j*str_x; 	        
       				level[lev].F[ind] = fvar[ind].F;        				
       			}                                
        	}
        }
               
        
}

/**************************************************************/
//Function to compute Pressure BCs based on momentum equation 
void compute_pbc(fval * fvar, pbcs& pbc, int up)
{
          int ind, ind1, ind2; 
  
	  /*Modifying evaluation of RHS in line with the B.Cs for the driven cavity problem*/
	  
	  tdma2x(fvar,0);  //Evaluate d2u/dx2 at all points in the domain
	  tdma2y(fvar,1);  //Evlauate d2v/dy2 at all points in the domain
	  
	  if(up==1) 
	  { 
	  	for(int i=0;i<=nx;i++) //Getting the pressure gradient in the top and bottom walls of the domain
	  	{    
	     		ind1 = i; 				     //Bottom wall
	     		ind2 = i + ny*str_x;      		     //Top wall
	     
	     		pbc.p1y[i] = (1./Re)*fvar[ind1].uyy[1];      //Bottom wall
	     		pbc.p2y[i] = (1./Re)*fvar[ind2].uyy[1];      //Top wall	     	     			     		
	  	}
	  
	  	for(int j=0;j<=ny;j++) //Getting the pressure gradient at the left and right faces of the domain 
	  	{
	     		ind1 = j*str_x; 			    //Left wall
	     		ind2 = nx + j*str_x;			    //Right wall
	     
	     		pbc.p1x[j] = (1./Re)*fvar[ind1].uxx[0];     //Left wall
	     		pbc.p2x[j] = (1./Re)*fvar[ind2].uxx[0];     //Right wall 	     
	  	} 
	  }
	  else
	  {
	  	for(int i=0;i<=nx;i++) //Getting the pressure gradient in the top and bottom walls of the domain
	  	{    
	     		ind1 = i; 			//Bottom wall
	     		ind2 = i + ny*str_x;      	//Top wall
	     
	     		pbc.p1y[i] = 0.0; 
	     		pbc.p2y[i] = 0.0; 
	     
	  	}
	  
	  	for(int j=0;j<=ny;j++) //Getting the pressure gradient at the left and right faces of the domain 
	  	{
	     		ind1 = j*str_x; 		//Left wall
	     		ind2 = nx + j*str_x;		//Right wall 
	     
	     		pbc.p1x[j] = 0.0; 
	     		pbc.p2x[j] = 0.0; 	     
	  	}	  
	  }
}

/*****************************************************************/

void mg_bcs_neu(mg_grid * level,int lev, pbcs& pbc)
{
	int ind;
	int sp = pow(2,lev); 
	
	/*ofstream t_wall,b_wall,left,right; 
  
	t_wall.open("mg_top_wall.dat"); 
	b_wall.open("mg_bot_wall.dat"); 
	left.open("mg_left_wall.dat"); 
	right.open("mg_right_wall.dat");*/ 
  
	 //Top wall 
 	 for(int i=0;i<=nx;i=i+sp)
	 {
		ind = i + ny*str_x; 						
		
		level[lev].p2y[i]=pbc.p2y[i];  //Should be 2.0. Given "0" for testing purposes
				
		level[lev].phi_s[ind]=  -(1./B1)*hy[lev]*level[lev].p2y[i] - B2*level[lev].phi_s[ind-sp*str_x] - B3*level[lev].phi_s[ind-2*sp*str_x] - B4*level[lev].phi_s[ind-3*sp*str_x]; 			
		
		//t_wall<<i<<"	"<<level[lev].phi_s[ind]<<"\n";
	 } 

	  //Right wall 
	 for(int j=0;j<=ny;j=j+sp)
	 {
		ind = nx + j*str_x; 	
		
		level[lev].p2x[j]=pbc.p2x[j]; //Should be 2.0 
				
		level[lev].phi_s[ind] = -(1./B1)*hx[lev]*level[lev].p2x[j] - B2*level[lev].phi_s[ind-sp] - B3*level[lev].phi_s[ind-2*sp] - B4*level[lev].phi_s[ind-3*sp];  
		
		//right<<j<<"	"<<level[lev].phi_s[ind]<<"\n"; 
	 }
  
	  //Left wall 
	 for(int j=0;j<=ny;j=j+sp) //changed starting j from 1 to 0
	 {
	 	ind = j*str_x;  	
	 	
		level[lev].p1x[j]=pbc.p1x[j];		
					
		level[lev].phi_s[ind] = (1./B1)*hx[lev]*level[lev].p1x[j] - ( B2*level[lev].phi_s[ind+sp] + B3*level[lev].phi_s[ind+2*sp] + B4*level[lev].phi_s[ind+3*sp] ) ;
		
		//left<<j<<"	"<<level[lev].phi_s[ind]<<"\n"; 
	 }

	  //Bottom wall 
	 for(int i=0;i<=nx;i=i+sp)  //changed starting i from 1 to 0
	 {		
	 	level[lev].p1y[i] = pbc.p1y[i]; 
	 	
		level[lev].phi_s[i] =   (1./B1)*hy[lev]*level[lev].p1y[i] - ( B2*level[lev].phi_s[i+sp*str_x] + B3*level[lev].phi_s[i+2*sp*str_x] + B4*level[lev].phi_s[i+3*sp*str_x] );
	
		//b_wall<<i<<"	"<<level[lev].phi_s[i]<<"\n"; 
		
		//fvar[i].u[2 ]= i*dx*i*dx; 
	 }
	 
	 /*t_wall.close();
	 b_wall.close();
	 left.close();
	 right.close();*/
}

/*****************************************************************************************************************/

void mg_clear_levels(mg_grid * level)
{
	int ind; 
	
	for(int lev=1;lev<mg_levels;lev++)
	{
		for(int j=0;j<=ny;j++)
		{
			for(int i=0;i<=nx;i++)
			{
				ind = i + j*str_x; 
				level[lev].phi_s[ind] = 0.0; 
				level[lev].F[ind] = 0.0; 					
			}		
		}	
	}	
}

/******************************************************************************************************************/

/*void mg_conjugate_gradient(vector<mg_grid>& level, int lev, int sm_ite, int up, int v_level, pbcs& pbc)
{
	int sp = pow(2,lev); 
	
	int st_inx = sp, en_inx = nx-sp; 
	int st_iny = sp, en_iny = ny-sp; 
	
	int nx_sol = (nx/pow(2,lev))-1, ny_sol = (ny/pow(2,lev))-1;  //No of points to be solved on the grid at each level in each direction
	int tot_p_sol = (nx_sol)*(ny_sol); 			     //No.of points that are solved for directly. 
	
	int ind; 				   
	int str_m = nx_sol; 
	int i_m, j_m, ind_m; 					   

	vec ri(tot_p_sol-1), rip(tot_p_sol-1), pi(tot_p_sol-1), pip(tot_p_sol-1); 
	
	vec ub = zeros<vec>(tot_p_sol), ux(tot_p_sol), pcorrec_vec=ones<vec>(tot_p_sol), b(tot_p_sol-1), act_ub = zeros<vec>(tot_p_sol); 
	
	vec x=zeros<vec>(tot_p_sol-1), Api(tot_p_sol-1);
	
	vec error(tot_p_sol-1), efromr(tot_p_sol-1), diffine(tot_p_sol-1); 

	double alphai,betai; 
	double tol=1.0e-12, rhs_sum=0.0;	 		
	int count=0;
	
	double res_norm=100.0,max_norm=100.0; 
	
	ofstream res_out; 
	res_out.open("mgcg_res_out.dat");
	
	//Point correction for solvability
	
	if(lev==v_level) //As both part of the ascending and descending parts of the cycle 
  	{     
      		mg_evaluate_rhs(level,lev,pbc,up);    
      		mg_compute_rhs_tot(level,lev); 
      		rhs_sum = mg_eval_rhs_sum(level,lev);             
  	}
  	
  	  
	if( (lev!=v_level) && (up==0))
  	{     	
      		for(int j=st_iny;j<=en_iny;j=j+sp)
      		{
		      	for(int i=st_inx;i<=en_inx;i=i+sp)
      			{
		      		ind = i + j*str_x;
      		
      				level[lev].phi_s[ind] = 0.0;       		      	 	      	
      			}      
      		}
   
      		rhs_sum = mg_eval_rhs_sum(level,lev);               
  	}
 
   
  	level[lev].point_correc = rhs_sum/(tot_p_sol); 			
	
	//Initializing the RHS vector 'b'
		
	for(int j=sp;j<=en_iny;j=j+sp)
	{
		for(int i=sp;i<=en_inx;i=i+sp)
		{
			ind = i + j*str_x; 		
											
			i_m = (i/sp)-1;
			j_m = (j/sp)-1; 
			
			ind_m = i_m + j_m*str_m; 						
			
			ub(ind_m) = level[lev].rhs[ind] - (level[lev].point_correc/level[lev].coeff[ind]);	//Untrimmed RHS	(Includes the last point)	
			act_ub(ind_m) = level[lev].coeff[ind]*ub(ind_m); 
			
			ux(ind_m) = level[lev].phi_s[ind];  							
		}                      
	}		
	
	b = ub.submat(0,0,tot_p_sol-2,0);   //Trimming the column matrix so that the pinned point is eliminated
	
	//vec prodfunc = mulvec(x,level); 
	
	//prodfunc.save("prodfunc.dat", raw_ascii);
	
	x = ux.submat(0,0,tot_p_sol-2,0); 
	
	
	if(v_level==2 && lev==3 && up==1) 
	{	
		ux.save("after_interpol_lev3_ux.dat",raw_ascii);
		ub.save("after_interpol_lev3_ub.dat",raw_ascii); 
	}
	
	ri = b - mulvec(x,level,lev); 	
	pi = ri; 	
	
	//res_norm = norm(ri,2);  		
		
	while(count<sm_ite && res_norm>tol)
	//while(res_norm>tol)
	{		
		Api = mulvec(pi,level,lev); 		
		alphai = dot(ri,ri)/dot(pi,Api);
								
		x = x + alphai*pi; 

		rip = ri - alphai*Api; 			
				
		betai  = dot(rip,rip)/dot(ri,ri);		
		pip = rip + betai*pi; 
				
		pi = pip;
		ri = rip; 	
						
	        res_norm = norm(ri,2); 				

		if(v_level==2 && lev==3 && up==1)
		{
			cout<<count<<"		"<<res_norm<<"\n"; 
			res_out<<count<<"		"<<res_norm<<"\n"; 
		}
		
		count++; 
	}
	
	cout<<"V_level: "<<v_level<<"	"<<"In Level: "<<lev<<"  	Up: "<<up<<"	"<<"Count:  "<<count<<"   Residual Norm:   "<<res_norm<<"	   "<<"Compatibility Sum: "<<sum(act_ub)<<"\n"; 		
	
	//double x_mean = sum(x)/(tot_p_sol);
	
	double x_mean=0.0;
	
	for(int j=sp;j<=en_iny;j=j+sp)
	{
		for(int i=sp;i<=en_inx;i=i+sp)
		{
			i_m = (i/sp)-1; 
			j_m = (j/sp)-1; 
			
			ind = i + j*str_x; 				
			ind_m = i_m + j_m*str_m; 
						
			//cout<<"i = "<<i<<"\nj = "<<j<<"i_m = "<<i_m<<"\nj_m = "<<j_m<<"\n";
					 
			if(ind_m<=tot_p_sol-2)
			{
				level[lev].phi_s[ind] = x(ind_m)-x_mean;				
				level[lev].res[ind] = ri(ind_m);
			}
			else
			{
				level[lev].phi_s[ind] = x_mean; 			
				level[lev].res[ind] = 0.0;			
			}												
		}
	}	
	
	mg_bcs_neu(level,lev,pbc);  //Applying bcs for each level
		
	res_out.close(); 
}*/

/*****************************************************************************************************/
//A function to compute the error vector from residual vector 

vec error_from_res(vec & b, mg_grid * level)
{
	int lev= 0; 
	int sp = pow(2,lev); 
	
	int nx_sol = (nx/pow(2,lev))-sp, ny_sol = (ny/pow(2,lev))-sp; //No of spacings on the grid to be solved for
	int tot_p_sol = (nx_sol)*(ny_sol); //No.of points that are solved for directly. The top and the farthest right lines of points are 	
					       //removed due to periodicity
	int ind; 				   
	int str_m = nx_sol; 
	int i_m,j_m,ind_m; 					   

	vec ri(tot_p_sol-1), rip(tot_p_sol-1), pi(tot_p_sol-1), pip(tot_p_sol-1); 
	
	vec zi(tot_p_sol-1), zip(tot_p_sol-1); 
	
	vec x=zeros<vec>(tot_p_sol-1), Api(tot_p_sol-1);
	
	vec error(tot_p_sol-1); 

	double alphai,betai; 
	double tol=1.0e-30; 	
	int count=0; 
	double res_norm=100.0,max_norm=100.0; 
	
	/****Initializing the RHS vector 'b'*****/ 
	mg_coeff();
	mg_read_coeff(level); 
	
	/****************************************/ 
	//Correcting the RHS  	
	
	double rhs_avg = sum(b)/(b.n_rows);
	 
	for(int j=sp;j<=ny_sol;j=j+sp)
	{
		for(int i=sp;i<=nx_sol;i=i+sp)
		{
			ind = i + j*str_x; 		
											
			i_m = (i/sp)-1;
			j_m = (j/sp)-1;
	
			if(i!=nx_sol && j!=ny_sol) b(ind_m) = b(ind_m) - (rhs_avg/level[lev].coeff[ind]); 
		}
	}
		
	/****************************************/ 
	//Initialization for CG method
			
	ri = b - mulvec(x,level,lev); 	
	
	zi=ri; 
	pi = zi; 	 	
		
	//while(res_norm>tol)
	while(count<50)
	{		
		Api = mulvec(pi,level,lev); 		
		alphai = dot(ri,zi)/dot(pi,Api);
								
		x = x + alphai*pi; 
		
		rip = ri - alphai*Api; 			
		zip = rip; 
		
		betai  = dot(rip,zip)/dot(ri,zi);		
		pip = zip + betai*pi; 
				
		pi = pip;
		ri = rip; 	
		zi = zip; 
						
	        res_norm = norm(ri,2); 						
		count++; 
	}
	
	return x; 
}

/*****************************************************************************************************/
//Error Norm calculation

/*double error_calc(fval * fvar, int lev)
{
  double abs_err=0.0; 
  int ind;
  
  int sp = pow(2,lev); 
  int nx_sol = nx/sp - 1, ny_sol = ny/sp - 1; 
  
  #pragma parallel for default(shared) private(i,j,ind) reduction(+:abs_err) 
  for(int j=sp;j<ny;j=j+sp)
  {
    for(int i=sp;i<nx;i=i+sp)
    {
      ind = i + j*str_x; 
      
      //fvar[ind].err = fvar[ind].u[2] - (i*dx*i*dx + j*dy*j*dy) + (nx-1)*(nx-1)*dx*dx + (ny-1)*(ny-1)*dy*dy;     

      //fvar[ind].err = fvar[ind].u[2] - sin(M_PI*i*dx)*cos(M_PI*j*dy) + sin(M_PI*(nx-1)*dx)*cos(M_PI*(ny-1)*dy);     
      
      fvar[ind].err = fvar[ind].p - sin(M_PI*i*dx)*cos(M_PI*j*dy) + sin(M_PI*(nx-1)*dx)*cos(M_PI*(ny-1)*dy);     
      
      //fvar[ind].err = fvar[ind].u[2] - (i*i*dx*dx) - (j*j*dy*dy)  + (nx-1)*(nx-1)*dx*dx + (ny-1)*(ny-1)*dy*dy;
      
      //fvar[ind].err = fvar[ind].u[2] - pow((i*dx),3.) - pow((j*dy),3.)  + pow((nx-1)*dx,3.) + pow((ny-1)*dy,3.);
      
      //fvar[ind].err = fvar[ind].u[2] - pow((i*dx),2.) - pow((j*dy),2.)  + pow((nx-1)*dx,2.) + pow((ny-1)*dy,2.);
      
      //fvar[ind].err = fvar[ind].u[2]; //For truncation error evaluation
      
      //fvar[ind].err = fvar[ind].u[2] - ( (1./3.)*( pow(i*dx,3.) + pow(j*dy,3.) ) - (1./2.)*( pow(i*dx,2.) + pow(j*dy,2.) ) ) + ( (1./3.)*( pow((nx-1)*dx,3.) + pow((ny-1)*dy,3.) ) - (1./2.)*( pow((nx-1)*dx,2.) + pow((ny-1)*dy,2.) ) );
      
      //fvar[ind].err = fvar[ind].u[2] - pow(i*dx,3) + pow((nx-1)*dx,3);
      
      //fvar[ind].err = fvar[ind].u[2] - pow(j*dy,3) + pow((ny-1)*dy,3);       
      
      //fvar[ind].err = fvar[ind].u[2] - ( (pow(i*dx,3)/3.) - (pow(i*dx,2)/2.) ) + ( (pow((nx-1)*dx,3)/3.) - (pow((nx-1)*dx,2)/2.) ); 
      
      //fvar[ind].err = fvar[ind].u[2] - ( (pow(j*dy,3.)/3.) - (pow(j*dy,2.)/2.) ) + ( (pow((ny-1)*dy,3.)/3.) - (pow((ny-1)*dy,2.)/2.) ); 
      
      //fvar[ind].err = fvar[ind].u[2] - (pow(i*dx,7.) + pow(j*dy,7.)) + (pow((nx-1)*dx,7.) + pow((ny-1)*dy,7.)) ; 
      
      abs_err = abs_err + fvar[ind].err*fvar[ind].err; 
    } 
 }	
  
  return sqrt(abs_err)/sqrt(nx_sol*ny_sol);   
}*/

/*******************************************************************************************************/
//Bi-CGSTAB for Neumann Poisson Problem 

void mg_bicgstab(mg_grid * level, int lev, int sm_ite, int up, int v_level, pbcs& pbc)
{
	int sp = pow(2,lev); 
	
	int st_inx = sp, en_inx = nx-sp; 
	int st_iny = sp, en_iny = ny-sp; 
	
	int nx_sol = (nx/pow(2,lev))-1, ny_sol = (ny/pow(2,lev))-1;  //No of points to be solved on the grid at each level in each direction
	int tot_p_sol = (nx_sol)*(ny_sol); 			     //No.of points that are solved for directly. 
	
	int ind; 				   
	int str_m = nx_sol; 
	int i_m, j_m, ind_m; 					   

	vec r0star=ones<vec>(tot_p_sol-1), ri(tot_p_sol-1), si(tot_p_sol-1), rip(tot_p_sol-1), pi(tot_p_sol-1), pip(tot_p_sol-1); 
	
	vec ub = zeros<vec>(tot_p_sol), ux(tot_p_sol), pcorrec_vec=ones<vec>(tot_p_sol), b(tot_p_sol-1), act_ub = zeros<vec>(tot_p_sol); 
	
	vec x=zeros<vec>(tot_p_sol-1), Api(tot_p_sol-1), Asi(tot_p_sol-1);
	
	vec error(tot_p_sol-1), efromr(tot_p_sol-1), diffine(tot_p_sol-1); 

	double alphai,betai,omegai; 
	double tol=1.0e-16, rhs_sum=0.0;	 		
	int count=0;
	double res_norm=100.0,max_norm=100.0; 
	
	/********Point correction for solvability**********/ 
	if(lev==v_level) //As both part of the ascending and descending parts of the cycle 
  	{     
      		mg_evaluate_rhs(level,lev,pbc,up);    
      		mg_compute_rhs_tot(level,lev); 
      		rhs_sum = mg_eval_rhs_sum(level,lev);             
  	}
  	
  	  
	if( (lev!=v_level) && (up==0))
  	{     	
      		for(int j=st_iny;j<=en_iny;j=j+sp)
      		{
		      	for(int i=st_inx;i<=en_inx;i=i+sp)
      			{
		      		ind = i + j*str_x;
      		
      				level[lev].phi_s[ind] = 0.0;       		      	 	      	
      			}      
      		}
   
      		rhs_sum = mg_eval_rhs_sum(level,lev);               
  	}
 
   
  	level[lev].point_correc = rhs_sum/(tot_p_sol); 			
	
	/****Initializing the RHS vector 'b'*****/ 
		
	for(int j=sp;j<=en_iny;j=j+sp)
	{
		for(int i=sp;i<=en_inx;i=i+sp)
		{
			ind = i + j*str_x; 		
											
			i_m = (i/sp)-1;
			j_m = (j/sp)-1; 
			
			ind_m = i_m + j_m*str_m; 						
			
			ub(ind_m) = level[lev].rhs[ind] - (level[lev].point_correc/level[lev].coeff[ind]);	//Untrimmed RHS	(Includes the last point)	
			act_ub(ind_m) = level[lev].coeff[ind]*ub(ind_m); 
			
			ux(ind_m) = level[lev].phi_s[ind];  							
		}                      
	}
		
	/****Initialization of 'b' done*******/ 
	
	b = ub.submat(0,0,tot_p_sol-2,0);   //Trimming the column matrix so that the pinned point is eliminated			
	x = ux.submat(0,0,tot_p_sol-2,0); 
	
	ri = b - mulvec(x,level,lev); 	
	pi = ri;
	r0star = ri;  				
	
	while(count<sm_ite && res_norm>tol)
	//while(res_norm>tol)
	{
		Api = mulvec(pi,level,lev);
		alphai = dot(ri,r0star)/dot(Api,r0star); 
		si = ri - alphai*Api;
		
		Asi = mulvec(si,level,lev); 
		omegai = dot(Asi,si)/dot(Asi,Asi); 
		
		x = x + (alphai*pi) + (omegai*si);  //This "x" can be assigned to 'phi_s' to assign the corner bc
		 
		rip = si - (omegai*Asi); 
		betai = ( dot(rip,r0star)/dot(ri,r0star) )*(alphai/omegai);   
		pip = rip + betai*(pi - omegai*Api); 
		
		ri = rip; 
		pi = pip;
		
		res_norm = norm(ri,2); 				
		
		count++; 
	}		
	
	
	for(int j=sp;j<=en_iny;j=j+sp)
	{
		for(int i=sp;i<=en_inx;i=i+sp)
		{
			i_m = (i/sp)-1; 
			j_m = (j/sp)-1; 
			
			ind = i + j*str_x; 				
			ind_m = i_m + j_m*str_m; 
						
			//cout<<"i = "<<i<<"\nj = "<<j<<"i_m = "<<i_m<<"\nj_m = "<<j_m<<"\n";
					 
			if(ind_m<=tot_p_sol-2)
			{
				level[lev].phi_s[ind] = x(ind_m);				
				level[lev].res[ind] = ri(ind_m);
			}
			else
			{
				level[lev].phi_s[ind] = 0.0; 			
				level[lev].res[ind] = 0.0;			
			}												
		}
	}	
	
	//mg_bcs_neu(level,lev,pbc);  //Applying bcs for each level
		
	//res_out.close(); 
}

/********************************************************************************************/
//Function generating the truncation error map 

void trunc_error_map(mg_grid * level,int lev, pbcs& pbc)
{
	int sp = pow(2,lev); 
	
	int st_inx = sp, en_inx = nx-sp; 
	int st_iny = sp, en_iny = ny-sp; 
	
	int nx_sol = (nx/pow(2,lev))-1, ny_sol = (ny/pow(2,lev))-1;  //No of points to be solved on the grid at each level in each direction
	int tot_p_sol = (nx_sol)*(ny_sol); 			     //No.of points that are solved for directly. 
	
	int ind; 				   
	int str_m = nx_sol; 
	int i_m, j_m, ind_m; 
	
	vec ub = zeros<vec>(tot_p_sol), ux(tot_p_sol), b(tot_p_sol-1), x(tot_p_sol-1), ri(tot_p_sol-1); 
	
	mg_evaluate_rhs(level,lev,pbc,1);    
	mg_compute_rhs_tot(level,lev); 

	/****Initializing the RHS vector 'b'*****/ 
		
	for(int j=sp;j<=en_iny;j=j+sp)
	{
		for(int i=sp;i<=en_inx;i=i+sp)
		{
			ind = i + j*str_x; 		
											
			i_m = (i/sp)-1;
			j_m = (j/sp)-1; 
			
			ind_m = i_m + j_m*str_m; 						
			
			ub(ind_m) = level[lev].rhs[ind];	//Untrimmed RHS	(Includes the last point)				
			ux(ind_m) = level[lev].phi_s[ind];  							
		}                      
	}
	
	
	/****Trimming of 'b' and 'x'****/ 
	
	b = ub.submat(0,0,tot_p_sol-2,0);   //Trimming the column matrix so that the pinned point is eliminated			
	x = ux.submat(0,0,tot_p_sol-2,0); 

	/***Computation of the truncation error***/ 
	
	ri = b - mulvec(x,level,lev); 		
	
	/***Writing it out to the field***/ 
	
	for(int j=sp;j<=en_iny;j=j+sp)
	{
		for(int i=sp;i<=en_inx;i=i+sp)
		{
			i_m = (i/sp)-1; 
			j_m = (j/sp)-1; 
			
			ind = i + j*str_x; 				
			ind_m = i_m + j_m*str_m; 
						
			//cout<<"i = "<<i<<"\nj = "<<j<<"i_m = "<<i_m<<"\nj_m = "<<j_m<<"\n";
					 
			if(ind_m<=tot_p_sol-2)
			{
				level[lev].phi_s[ind] = ri(ind_m);								
			}
			else
			{
				level[lev].phi_s[ind] = 0.0; 						
			}												
		}
	}	
	
}

/********************************************************************************************/
//Arma Direct Solver for truncation error mapping
void arma_direct_trunc(mg_grid * level, int lev)
{
	int ind, ind_m, i_m,j_m; 
	
	int sp = pow(2,lev);
	
	int st_inx=sp, st_iny=sp;  
	int en_inx = nx-sp, en_iny = ny-sp; 
		
	int nx_sol = (nx/sp)-1, ny_sol = (ny/sp)-1; //No.of points that solved for directly. Only interior points are mentioned. 
	int tot_p_sol = nx_sol*ny_sol; 
	int str_m = nx_sol; 
	
	//cout<<"strm is "<<str_m<<"\n";
	//cout<<"total points are "<<tot_p_sol<<"\n"; 
	
	mat A = zeros<mat>(tot_p_sol,tot_p_sol);  	
	mat B = zeros<mat>(tot_p_sol,1), X = zeros<mat>(tot_p_sol,1), trunc=zeros<mat>(tot_p_sol,1);	
	
	//vec uni = ones<vec>(tot_p_sol-1); 
		
	for(int j=sp;j<ny;j=j+sp)
	{
		for(int i=sp;i<nx;i=i+sp)
		{
			ind = i + j*str_x; 	
			
			i_m = (i/sp) - 1; 
			j_m = (j/sp) - 1; 		

			ind_m = i_m + j_m*str_m; 
			
			//cout<<"i_m= "<<i_m<<"j_m= "<<j_m<<"\n"; 
			//cout<<"Act-ind-m:   "<<ind_m<<"\n"; 
			
			/*The various cases would start now*/
			
			/********Case-1*******/
			
			/*******************Case-1 j=st_iny***********************/ 
			if(j==st_iny)
			{
			  if(i==st_inx)
			  { 			   			    
			    A(ind_m,ind_m) = dc_b1_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b1_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b1_13[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b2_11[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_b2_12[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_b2_13[lev]; 
			    
			    A(ind_m,ind_m+2*str_m) = dc_b3_11[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_b3_12[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_b3_13[lev]; 			    
			    
			  }
			  else if(i==st_inx+sp)
			  {     			    			    
			    A(ind_m,ind_m-1) = dc_nb1_21[lev]; 
			    A(ind_m,ind_m) = dc_nb1_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb1_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb1_24[lev]; 
			    
			    A(ind_m,ind_m-1+str_m) = dc_nb2_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb2_22[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_nb2_23[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_nb2_24[lev]; 
			    
			    A(ind_m,ind_m-1+2*str_m) = dc_nb3_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb3_22[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_nb3_23[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_nb3_24[lev]; 			    
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
			  {			    
			    A(ind_m,ind_m-2) = dc_i1_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i1_33[lev]; 
			    A(ind_m,ind_m) = dc_i1_34[lev];
			    A(ind_m,ind_m+1) = dc_i1_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i1_36[lev]; 
			    
			    A(ind_m,ind_m-2+str_m) = dc_i2_32[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_i2_33[lev]; 
			    A(ind_m,ind_m+str_m) = dc_i2_34[lev];
			    A(ind_m,ind_m+1+str_m) = dc_i2_35[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_i2_36[lev];
			    
			    A(ind_m,ind_m-2+2*str_m) = dc_i3_32[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_i3_33[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_i3_34[lev];
			    A(ind_m,ind_m+1+2*str_m) = dc_i3_35[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_i3_36[lev]; 			    
			  }
			  else if(i==en_inx-sp)
			  {    			     			    
			    A(ind_m,ind_m-2) = dc_nb1_24[lev];
    			    A(ind_m,ind_m-1) = dc_nb1_23[lev]; 
			    A(ind_m,ind_m) = dc_nb1_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb1_21[lev]; 
			    
			    A(ind_m,ind_m-2+str_m) = dc_nb2_24[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_nb2_23[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb2_22[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_nb2_21[lev]; 			    
			    
			    A(ind_m,ind_m-2+2*str_m) = dc_nb3_24[lev]; 			    			  
			    A(ind_m,ind_m-1+2*str_m) = dc_nb3_23[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb3_22[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_nb3_21[lev]; 		     						    
			  }
			  else 
			  {		    				    
			    A(ind_m,ind_m) = dc_b1_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b1_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b1_13[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b2_11[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_b2_12[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_b2_13[lev]; 
			    
			    A(ind_m,ind_m+2*str_m) = dc_b3_11[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_b3_12[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_b3_13[lev];		    
			  }	  	
			}
			
			/*******************Case-2 j=st_iny+sp***********************/ 
			
			if(j==st_iny+sp)
			{
			  if(i==st_inx)
			  {			    
			    A(ind_m,ind_m-str_m) = dc_b4_11[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_b4_12[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_b4_13[lev]; 
			    
			    A(ind_m,ind_m) = dc_b5_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b5_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b5_13[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b6_11[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_b6_12[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_b6_13[lev]; 
			    
			    A(ind_m,ind_m+2*str_m) = dc_b7_11[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_b7_12[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_b7_13[lev];			    
			  }
			  else if(i==st_inx+sp) 
			  {	    			    		    
			    A(ind_m,ind_m-1-str_m) = dc_nb4_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb4_22[lev]; 
			    A(ind_m,ind_m+1-str_m)= dc_nb4_23[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_nb4_24[lev]; 
			    
			    A(ind_m,ind_m-1) = dc_nb5_21[lev]; 
			    A(ind_m,ind_m) = dc_nb5_22[lev]; 
			    A(ind_m,ind_m+1)= dc_nb5_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb5_24[lev]; 
			    	    	    
			    A(ind_m,ind_m-1+str_m) = dc_nb6_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb6_22[lev]; 
			    A(ind_m,ind_m+1+str_m)= dc_nb6_23[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_nb6_24[lev];
			    
			    A(ind_m,ind_m-1+2*str_m) = dc_nb7_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb7_22[lev]; 
			    A(ind_m,ind_m+1+2*str_m)= dc_nb7_23[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_nb7_24[lev];			    
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
			  {			    
			    A(ind_m,ind_m-2-str_m) = dc_i4_32[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_i4_33[lev]; 
			    A(ind_m,ind_m-str_m) = dc_i4_34[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_i4_35[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_i4_36[lev];			    	  
			    
			    A(ind_m,ind_m-2) = dc_i5_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i5_33[lev]; 
			    A(ind_m,ind_m)   = dc_i5_34[lev]; 
			    A(ind_m,ind_m+1) = dc_i5_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i5_36[lev]; 
			    
			    A(ind_m,ind_m-2+str_m) = dc_i6_32[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_i6_33[lev]; 
			    A(ind_m,ind_m+str_m)   = dc_i6_34[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_i6_35[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_i6_36[lev]; 
			    
			    A(ind_m,ind_m-2+2*str_m) = dc_i7_32[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_i7_33[lev]; 
			    A(ind_m,ind_m+2*str_m)   = dc_i7_34[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_i7_35[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_i7_36[lev];    		     
			  }
			  else if(i==en_inx-sp)
			  {	    			       
			    A(ind_m,ind_m+1-str_m) = dc_nb4_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb4_22[lev]; 
			    A(ind_m,ind_m-1-str_m)= dc_nb4_23[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_nb4_24[lev]; 
			    
			    A(ind_m,ind_m+1) = dc_nb5_21[lev]; 
			    A(ind_m,ind_m) = dc_nb5_22[lev]; 
			    A(ind_m,ind_m-1)= dc_nb5_23[lev]; 
			    A(ind_m,ind_m-2) = dc_nb5_24[lev]; 
			    	    	    
			    A(ind_m,ind_m+1+str_m) = dc_nb6_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb6_22[lev]; 
			    A(ind_m,ind_m-1+str_m)= dc_nb6_23[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_nb6_24[lev];
			    
			    A(ind_m,ind_m+1+2*str_m) = dc_nb7_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb7_22[lev]; 
			    A(ind_m,ind_m-1+2*str_m)= dc_nb7_23[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_nb7_24[lev];		      
			  }
			  else 
			  {	    
    
			    A(ind_m,ind_m-str_m) = dc_b4_11[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_b4_12[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_b4_13[lev]; 
			    
			    A(ind_m,ind_m) = dc_b5_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b5_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b5_13[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b6_11[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_b6_12[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_b6_13[lev]; 
			    
			    A(ind_m,ind_m+2*str_m) = dc_b7_11[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_b7_12[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_b7_13[lev];	    			    	  
			  }  
			}
			
			/*******************Case-3 j>=st_iny+2*sp && j<=en_iny-2*sp*************/
			
			if(j>=st_iny+2*sp && j<=en_iny-2*sp)
			{
			  if(i==st_inx)
			  {	    			    	    
			    A(ind_m,ind_m-2*str_m) = dc_b8_11[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_b8_12[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_b8_13[lev];
			    
			    A(ind_m,ind_m-str_m) = dc_b9_11[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_b9_12[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_b9_13[lev];
			    
			    A(ind_m,ind_m) = dc_b10_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b10_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b10_13[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b11_11[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_b11_12[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_b11_13[lev]; 
			    
			    A(ind_m,ind_m+2*str_m) = dc_b12_11[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_b12_12[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_b12_13[lev];  
			  }
			  else if(i==st_inx+sp) 
			  {	    			    
			    
			    A(ind_m,ind_m-1-2*str_m) = dc_nb8_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb8_22[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_nb8_23[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_nb8_24[lev]; 
			    
			    A(ind_m,ind_m-1-str_m) = dc_nb9_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb9_22[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_nb9_23[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_nb9_24[lev];	
			    
			    A(ind_m,ind_m-1) = dc_nb10_21[lev]; 
			    A(ind_m,ind_m) = dc_nb10_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb10_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb10_24[lev];
			    
			    A(ind_m,ind_m-1+str_m) = dc_nb11_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb11_22[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_nb11_23[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_nb11_24[lev];
			    
			    A(ind_m,ind_m-1+2*str_m) = dc_nb12_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb12_22[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_nb12_23[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_nb12_24[lev];   		    
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
			  {    			   			    
			    A(ind_m,ind_m-1-2*str_m) = dc_i8_33[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_i8_34[lev];
			    A(ind_m,ind_m+1-2*str_m) = dc_i8_35[lev];
			    
			    A(ind_m,ind_m-2-str_m) = dc_i9_32[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_i9_33[lev]; 
			    A(ind_m,ind_m-str_m) = dc_i9_34[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_i9_35[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_i9_36[lev];
			    
			    A(ind_m,ind_m-2) = dc_i10_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i10_33[lev]; 
			    A(ind_m,ind_m) = dc_i10_34[lev]; 
			    A(ind_m,ind_m+1) = dc_i10_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i10_36[lev];
			    
			    A(ind_m,ind_m-2+str_m) = dc_i11_32[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_i11_33[lev]; 
			    A(ind_m,ind_m+str_m) = dc_i11_34[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_i11_35[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_i11_36[lev];
			    
			    A(ind_m,ind_m-1+2*str_m) = dc_i12_33[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_i12_34[lev];
			    A(ind_m,ind_m+1+2*str_m) = dc_i12_35[lev];			    
			  }
			  else if(i==en_inx-sp)
			  {	    			     	  			    
			    A(ind_m,ind_m+1-2*str_m) = dc_nb8_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb8_22[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_nb8_23[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_nb8_24[lev]; 
			    
			    A(ind_m,ind_m+1-str_m) = dc_nb9_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb9_22[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_nb9_23[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_nb9_24[lev];	
			    
			    A(ind_m,ind_m+1) = dc_nb10_21[lev]; 
			    A(ind_m,ind_m) = dc_nb10_22[lev]; 
			    A(ind_m,ind_m-1) = dc_nb10_23[lev]; 
			    A(ind_m,ind_m-2) = dc_nb10_24[lev];
			    
			    A(ind_m,ind_m+1+str_m) = dc_nb11_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb11_22[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_nb11_23[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_nb11_24[lev];
			    
			    A(ind_m,ind_m+1+2*str_m) = dc_nb12_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb12_22[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_nb12_23[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_nb12_24[lev];			    
			  }
			  else 
			  {    			    		    
			    A(ind_m,ind_m-2*str_m) = dc_b8_11[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_b8_12[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_b8_13[lev];
			    
			    A(ind_m,ind_m-str_m) = dc_b9_11[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_b9_12[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_b9_13[lev];
			    
			    A(ind_m,ind_m) = dc_b10_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b10_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b10_13[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b11_11[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_b11_12[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_b11_13[lev]; 
			    
			    A(ind_m,ind_m+2*str_m) = dc_b12_11[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_b12_12[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_b12_13[lev];		      
			  }	  	  
			}	
					
			/************************Case-4 j==en_iny-sp*************************************/
			if(j==en_iny-sp)
			{
			  if(i==st_inx)
			  {	    			    			    
			    A(ind_m,ind_m+str_m) = dc_b4_11[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_b4_12[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_b4_13[lev]; 
			    
			    A(ind_m,ind_m) = dc_b5_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b5_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b5_13[lev]; 
			    
			    A(ind_m,ind_m-str_m) = dc_b6_11[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_b6_12[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_b6_13[lev]; 
			    
			    A(ind_m,ind_m-2*str_m) = dc_b7_11[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_b7_12[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_b7_13[lev];			      
			  }
			  else if(i==st_inx+sp) 
			  {			       			    
			    A(ind_m,ind_m-1+str_m) = dc_nb4_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb4_22[lev]; 
			    A(ind_m,ind_m+1+str_m)= dc_nb4_23[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_nb4_24[lev]; 
			    
			    A(ind_m,ind_m-1) = dc_nb5_21[lev]; 
			    A(ind_m,ind_m) = dc_nb5_22[lev]; 
			    A(ind_m,ind_m+1)= dc_nb5_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb5_24[lev]; 
			    	    	    
			    A(ind_m,ind_m-1-str_m) = dc_nb6_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb6_22[lev]; 
			    A(ind_m,ind_m+1-str_m)= dc_nb6_23[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_nb6_24[lev];
			    
			    A(ind_m,ind_m-1-2*str_m) = dc_nb7_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb7_22[lev]; 
			    A(ind_m,ind_m+1-2*str_m)= dc_nb7_23[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_nb7_24[lev];    		    
			  }
			  else if(i>st_inx+sp && i<en_inx-sp)
			  {  	  			    
			    A(ind_m,ind_m-2+str_m) = dc_i4_32[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_i4_33[lev]; 
			    A(ind_m,ind_m+str_m) = dc_i4_34[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_i4_35[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_i4_36[lev];			    	  
			    
			    A(ind_m,ind_m-2) = dc_i5_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i5_33[lev]; 
			    A(ind_m,ind_m)   = dc_i5_34[lev]; 
			    A(ind_m,ind_m+1) = dc_i5_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i5_36[lev]; 
			    
			    A(ind_m,ind_m-2-str_m) = dc_i6_32[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_i6_33[lev]; 
			    A(ind_m,ind_m-str_m)   = dc_i6_34[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_i6_35[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_i6_36[lev]; 
			    
			    A(ind_m,ind_m-2-2*str_m) = dc_i7_32[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_i7_33[lev]; 
			    A(ind_m,ind_m-2*str_m)   = dc_i7_34[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_i7_35[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_i7_36[lev];			    
			  }
			  else if(i==en_inx-sp)
			  {	    			    	  	    			    			    
			    A(ind_m,ind_m+1+str_m) = dc_nb4_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb4_22[lev]; 
			    A(ind_m,ind_m-1+str_m)= dc_nb4_23[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_nb4_24[lev]; 
			    
			    A(ind_m,ind_m+1) = dc_nb5_21[lev]; 
			    A(ind_m,ind_m) = dc_nb5_22[lev]; 
			    A(ind_m,ind_m-1)= dc_nb5_23[lev]; 
			    A(ind_m,ind_m-2) = dc_nb5_24[lev]; 
			    	    	    
			    A(ind_m,ind_m+1-str_m) = dc_nb6_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb6_22[lev]; 
			    A(ind_m,ind_m-1-str_m)= dc_nb6_23[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_nb6_24[lev];
			    
			    A(ind_m,ind_m+1-2*str_m) = dc_nb7_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb7_22[lev]; 
			    A(ind_m,ind_m-1-2*str_m)= dc_nb7_23[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_nb7_24[lev];			    
			  }
			  else 
			  {    			      	    			    
			    A(ind_m,ind_m+str_m) = dc_b4_11[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_b4_12[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_b4_13[lev]; 
			    
			    A(ind_m,ind_m) = dc_b5_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b5_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b5_13[lev]; 
			    
			    A(ind_m,ind_m-str_m) = dc_b6_11[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_b6_12[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_b6_13[lev]; 
			    
			    A(ind_m,ind_m-2*str_m) = dc_b7_11[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_b7_12[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_b7_13[lev];		    
			  }  
			}			
			
			/******************Case-5(j=en_iny)**********************/ 
	
			if(j==en_iny)
			{
			  if(i==st_inx)
			  {	    			        			    
			    A(ind_m,ind_m) = dc_b1_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b1_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b1_13[lev]; 
			    
			    A(ind_m,ind_m-str_m) = dc_b2_11[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_b2_12[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_b2_13[lev]; 
			    
			    A(ind_m,ind_m-2*str_m) = dc_b3_11[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_b3_12[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_b3_13[lev];			    
			  }
			  else if(i==st_inx+sp)
			  {	    			     	  			    
			    A(ind_m,ind_m-1) = dc_nb1_21[lev]; 
			    A(ind_m,ind_m) = dc_nb1_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb1_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb1_24[lev]; 
			    
			    A(ind_m,ind_m-1-str_m) = dc_nb2_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb2_22[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_nb2_23[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_nb2_24[lev]; 
			    
			    A(ind_m,ind_m-1-2*str_m) = dc_nb3_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb3_22[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_nb3_23[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_nb3_24[lev];			    
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)
			  {	    			    			    			    
			    A(ind_m,ind_m-2) = dc_i1_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i1_33[lev]; 
			    A(ind_m,ind_m) = dc_i1_34[lev];
			    A(ind_m,ind_m+1) = dc_i1_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i1_36[lev]; 
			    
			    A(ind_m,ind_m-2-str_m) = dc_i2_32[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_i2_33[lev]; 
			    A(ind_m,ind_m-str_m) = dc_i2_34[lev];
			    A(ind_m,ind_m+1-str_m) = dc_i2_35[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_i2_36[lev];
			    
			    A(ind_m,ind_m-2-2*str_m) = dc_i3_32[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_i3_33[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_i3_34[lev];
			    A(ind_m,ind_m+1-2*str_m) = dc_i3_35[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_i3_36[lev];
			  }
			  else if(i==en_inx-sp)
			  {    			    		    
			    A(ind_m,ind_m-2) = dc_nb1_24[lev];
    			    A(ind_m,ind_m-1) = dc_nb1_23[lev]; 
			    A(ind_m,ind_m) = dc_nb1_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb1_21[lev]; 
			    
			    A(ind_m,ind_m-2-str_m) = dc_nb2_24[lev]; 			    
			    A(ind_m,ind_m-1-str_m) = dc_nb2_23[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb2_22[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_nb2_21[lev]; 			    
			    
			    A(ind_m,ind_m-2-2*str_m) = dc_nb3_24[lev]; 			    			  
			    A(ind_m,ind_m-1-2*str_m) = dc_nb3_23[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb3_22[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_nb3_21[lev];			        
			  }
			  else 
			  {   				  	       	    	    		    		    			     
			    A(ind_m,ind_m) = dc_b1_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b1_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b1_13[lev]; 
			    
			    A(ind_m,ind_m-str_m) = dc_b2_11[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_b2_12[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_b2_13[lev]; 
			    
			    A(ind_m,ind_m-2*str_m) = dc_b3_11[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_b3_12[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_b3_13[lev]; 
			     
			  }	  	
			}		
			
			/**********************************************/
			//Populating the full rhs matrix 
			
			//B[ind_m] = ( level[lev].rhs[ind] - (level[lev].point_correc/level[lev].coeff[ind]) );  							
			B[ind_m] = level[lev].rhs[ind]; 
				
		}	
	}

	/**************************Code for Truncation mapping starts********************/
	//Initializing the field 
		
	for(int j=sp;j<ny;j=j+sp)
	{
		for(int i=sp;i<nx;i=i+sp)
		{
			i_m = (i/sp) - 1; 
			j_m = (j/sp) - 1; 
			
			ind = i + j*str_x; 	
			ind_m = i_m + j_m*str_m; 
			
			X(ind_m) = level[lev].phi_s[ind];	

		}
	}
	
	/*********************************************/
	//Computing the truncation error field
	
	trunc = B - A*X; 
	
	/*********************************************/
	//Reassigning the truncation error field 
	
	for(int j=sp;j<ny;j=j+sp)
	{
		for(int i=sp;i<nx;i=i+sp)
		{
			i_m = (i/sp) - 1; 
			j_m = (j/sp) - 1; 
			
			ind = i + j*str_x; 	
			ind_m = i_m + j_m*str_m; 
			
			level[lev].phi_s[ind] = trunc(ind_m); 

		}
	}		
}
 
/*************************************************************************************************/		
//BICGSTAB algorithm modified to include the corner points
 
void mg_bicgstab_corn(fval * fvar, mg_grid * level, int lev, int sm_ite, int up, int v_level, pbcs& pbc)
{
	int sp = pow(2,lev); 
	
	int st_inx = sp, en_inx = nx-sp; 
	int st_iny = sp, en_iny = ny-sp; 
	
	int nx_sol = (nx/pow(2,lev))-1, ny_sol = (ny/pow(2,lev))-1;  //No of points to be solved on the grid at each level in each direction
	int tot_p_sol = (nx_sol)*(ny_sol); 			     //No.of points that are solved for directly. 
		
	int ind; 				   
	int str_m = nx_sol; 
	int i_m, j_m, ind_m; 					   

	vec r0star=ones<vec>(tot_p_sol-1), ri(tot_p_sol-1), si(tot_p_sol-1), rip(tot_p_sol-1), pi(tot_p_sol-1), pip(tot_p_sol-1); 
	
	vec ub = zeros<vec>(tot_p_sol), ux(tot_p_sol), pcorrec_vec=ones<vec>(tot_p_sol), b(tot_p_sol-1), act_ub = zeros<vec>(tot_p_sol); 
	
	vec x=zeros<vec>(tot_p_sol-1), Api(tot_p_sol-1), Asi(tot_p_sol-1);
	
	vec error(tot_p_sol-1), efromr(tot_p_sol-1), diffine(tot_p_sol-1); 

	double alphai,betai,omegai; 
	double tol=1.0e-16, rhs_sum=0.0;	 		
	int count=0, num_count;
	double res_norm=100.0,max_norm=100.0;
	
	/****************************************/
	//ofstream res_out; 
	//res_out.open("mgcg_res_out.dat");
	if( up=0 && (lev!=v_level) )
        {
                num_count = 1;
        }
        else
        {
                num_count = 5;
        }	
	
	if( (lev!=v_level) && (up==0) )
  	{     	
      		for(int j=st_iny;j<=en_iny;j=j+sp)
      		{
		      	for(int i=st_inx;i<=en_inx;i=i+sp)
      			{
		      		ind = i + j*str_x;
      	
      				level[lev].phi_s[ind] = 0.0;       		      	 	      	
      			}      
      		}
   
   		rhs_sum = mg_eval_rhs_sum(level,lev);               
  	}	
	
	for(int glob_count=0; glob_count<num_count;glob_count++)
	{
		if(lev==v_level) //As both part of the ascending and descending parts of the cycle 
  		{     
      			mg_evaluate_rhs(level,lev,pbc,up);    
      			mg_compute_rhs_tot(level,lev); 
      			rhs_sum = mg_eval_rhs_sum(level,lev);             
  		}
  	
  	  
		/*if( (lev!=v_level) && (up==0))
  		{     	
      			for(int j=st_iny;j<=en_iny;j=j+sp)
      			{
			      	for(int i=st_inx;i<=en_inx;i=i+sp)
      				{
			      		ind = i + j*str_x;
      		
      					level[lev].phi_s[ind] = 0.0;       		      	 	      	
      				}      
      			}
   
      			rhs_sum = mg_eval_rhs_sum(level,lev);               
  		}*/
 
   
  		level[lev].point_correc = rhs_sum/(tot_p_sol); 			
	
		/****Initializing the RHS vector 'b'*****/ 
		
		for(int j=sp;j<=en_iny;j=j+sp)
		{
			for(int i=sp;i<=en_inx;i=i+sp)
			{
				ind = i + j*str_x; 		
											
				i_m = (i/sp)-1;
				j_m = (j/sp)-1; 
			
				ind_m = i_m + j_m*str_m; 						
			
				ub(ind_m) = level[lev].rhs[ind] - (level[lev].point_correc/level[lev].coeff[ind]);	//Untrimmed RHS	(Includes the last point)	
				act_ub(ind_m) = level[lev].coeff[ind]*ub(ind_m); 
			
				ux(ind_m) = level[lev].phi_s[ind]; 										
			}                      
		}
		
		/****Initialization of 'b' done*******/ 
	
		b = ub.submat(0,0,tot_p_sol-2,0);   //Trimming the column matrix so that the pinned point is eliminated			
		x = ux.submat(0,0,tot_p_sol-2,0); 
	
		ri = b - mulvec(x,level,lev); 	
		pi = ri;
		r0star = ri;  		
		
		/*if(v_level==2 && up==0) ub.save("ub_vlevel2.dat",raw_ascii); 				
		if(v_level==2 && up==0) ux.save("ux_vlevel2.dat",raw_ascii); 				
		if(v_level==2 && up==0) ri.save("ri_vlevel2.dat",raw_ascii); 				
		if(v_level==2 && up==0) cout<<"Denominator is "<<dot( mulvec(pi,level,lev),r0star)<<"\n";*/
			
		while(count<sm_ite && res_norm>tol)
		//while(res_norm>tol)
		{
			Api = mulvec(pi,level,lev);
			alphai = dot(ri,r0star)/dot(Api,r0star); 
			si = ri - alphai*Api;
		
			Asi = mulvec(si,level,lev); 
			omegai = dot(Asi,si)/dot(Asi,Asi); 
		
			x = x + (alphai*pi) + (omegai*si);  //This "x" can be assigned to 'phi_s' to assign the corner bc
			 
			rip = si - (omegai*Asi); 
			betai = ( dot(rip,r0star)/dot(ri,r0star) )*(alphai/omegai);   
			pip = rip + betai*(pi - omegai*Api); 
		
			ri = rip; 
			pi = pip;
		
			res_norm = norm(ri,2); 				

			//cout<<count<<"		"<<res_norm<<"\n"; 
			//res_out<<count<<"	"<<res_norm<<"\n"; 
		
			count++; 
		}			 		
	
		for(int j=sp;j<=en_iny;j=j+sp)
		{
			for(int i=sp;i<=en_inx;i=i+sp)
			{
				i_m = (i/sp)-1; 
				j_m = (j/sp)-1; 
			
				ind = i + j*str_x; 				
				ind_m = i_m + j_m*str_m; 
						
				//cout<<"i = "<<i<<"\nj = "<<j<<"i_m = "<<i_m<<"\nj_m = "<<j_m<<"\n";
						 
				if(ind_m<=tot_p_sol-2)
				{
					level[lev].phi_s[ind] = x(ind_m);				
					level[lev].res[ind] = ri(ind_m);
				}
				else
				{
					level[lev].phi_s[ind] = 0.0; 			
					level[lev].res[ind] = 0.0;			
				}												
			}
		}	
	}
	
	//cout<<"V_level: "<<v_level<<"	"<<"In Level: "<<lev<<"  	Up: "<<up<<"	"<<"Count:  "<<count<<"   Residual Norm:   "<<res_norm<<"	   "<<"Compatibility Sum: "<<sum(act_ub)<<"\n";
	
	//mg_bcs_neu(level,lev,pbc);  //Applying bcs for each level
		
	//res_out.close(); 
}
/************************************************************************************************/
//Function that designs an extrapolation formula using least squares method using a 4th order general polynomial in "X" and "Y" direction 
double extra_pol(int i_in, int j_in, mg_grid * level, int lev)
{
//i_in: index in x-direction of the said stencil. This is the first x-index of the extrapating stencil 
//j_in: index in y-direction of the said stencil. This is the first y-index of the extrapating stencil 

int st_size = 4; 	//Stencil size. 4X4 
int num_points = 16; 
int num_coeff = 25; 

int sp=pow(2,lev); 

int i_start = i_in;  //Start and end of the stencil indices 
int j_start = j_in; 

int i_end = i_in+3*sp; 
int j_end = j_in+3*sp; 

int i_extra, j_extra, ind, ind_extra, ind_gl; 
double x,y;

/////////

if(i_in==0 && j_in==0)
{
	i_extra = 0;
	j_extra = 0; 
}

if(i_in==0 && j_in==ny-3*sp)
{
	i_extra=0;
	j_extra=3; 
}

if(i_in==nx-3*sp && j_in==ny-3*sp)
{
	i_extra = 3; 
	j_extra = 3; 
}

if(i_in==nx-3*sp && j_in==0)
{
	i_extra = 3;
	j_extra = 0; 
}

/////////

mat coeff = zeros<mat>(num_coeff,1); 
mat inter = zeros<mat>(num_points,num_coeff); 
mat rhs = zeros<mat>(num_points,1); 

mat trunc_inter = zeros<mat>(num_points,num_coeff); 
mat trunc_rhs = zeros<mat>(num_points-1,1); 

for(int j=j_start;j<=j_end;j=j+sp)
{
	for(int i=i_start;i<=i_end;i=i+sp)
	{
		x = i*dx; 
		y = j*dy; 
		
		ind = (i-i_start)/sp + ( (j-j_start)/sp )*st_size; 
		
		ind_gl = i + j*str_x; 
		
		inter(ind,0) = pow(x,4)*pow(y,4);
        	inter(ind,1) = pow(x,3)*pow(y,4);
	        inter(ind,2) = pow(x,2)*pow(y,4);
        	inter(ind,3) = x*pow(y,4);
        	inter(ind,4) = pow(y,4);
        	inter(ind,5) = pow(x,4)*pow(y,3);
        	inter(ind,6) = pow(x,3)*pow(y,3);
        	inter(ind,7) = pow(x,2)*pow(y,3); 
        	inter(ind,8) = x*pow(y,3);                      
        	inter(ind,9) = pow(y,3);
        	inter(ind,10) = pow(x,4)*pow(y,2);
        	inter(ind,11) = pow(x,3)*pow(y,2);
        	inter(ind,12) = pow(x,2)*pow(y,2);
        	inter(ind,13) = x*pow(y,2);
        	inter(ind,14) = pow(y,2);
        	inter(ind,15) = pow(x,4)*y;
        	inter(ind,16) = pow(x,3)*y;
        	inter(ind,17) = pow(x,2)*y;
        	inter(ind,18) = x*y;
        	inter(ind,19) = y; 
        	inter(ind,20) = pow(x,4);
        	inter(ind,21) = pow(x,3);
        	inter(ind,22) = pow(x,2);
        	inter(ind,23) = x;
        	inter(ind,24) = 1.0;	
        	
        	rhs(ind) = level[lev].phi_s[ind_gl];        	        	
	}	
}

ind_extra = i_extra + j_extra*st_size; 

//cout<<"The index to be extrapolated is "<<ind_extra<<"\n";

trunc_inter = inter; 
trunc_rhs = rhs; 

trunc_inter.shed_row(ind_extra); 
trunc_rhs.shed_row(ind_extra); 

//cout<<"Size of the truncated matrix is "<<trunc_inter.size()<<"\n";

//trunc_inter.save("trunc_inter.dat",raw_ascii); 
//trunc_rhs.save("trunc_rhs.dat",raw_ascii); 

coeff = solve(trunc_inter,trunc_rhs);  

mat row_slice = inter.row(ind_extra);

return dot(row_slice,coeff);
}	
/************************************************************************************************/
