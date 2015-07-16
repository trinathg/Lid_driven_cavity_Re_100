#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

vec old_mulvec(vec X, mg_grid * level)
{ 	
	int lev=0; 
	int sp=pow(2,lev), ind; 
	
	int nx_sol = (nx/pow(2,lev))-sp, ny_sol = (ny/pow(2,lev))-sp; //No of points for the computational domain after compatibility condition 
	int tot_p_sol = (nx_sol)*(ny_sol); //No.of points that are solved for directly. The top and the farthest right lines of points are  	
					   //removed due to periodicity	
					   
	vec ans(tot_p_sol), subans(tot_p_sol-1); 				   
			
	int i_m, j_m, str_m=nx_sol, ind_m; 
	
	int st_inx = sp, en_inx = nx - sp;    
  	int st_iny = sp, en_iny = ny - sp;  
  	
	double RHS1, RHS2, RHS3, RHS4, RHS5, RHS6; 
	
	for(int j=st_iny;j<=en_iny;j=j+sp)
	{
		for(int i=st_inx;i<=en_inx;i=i+sp)
		{
			ind = i + j*str_x; 						
			
			i_m = (i/sp)-1; 
			j_m = (j/sp)-1; 
			
			ind_m = i_m + j_m*str_m; 
			
			/*******************Case-1 j=1***********************/ 
			if(j==st_iny)
			{
			  if(i==st_inx) //Checked 
			  {	        	    	    
			    RHS1 = dc_b1_11[lev]*X[ind_m] + dc_b1_12[lev]*X[ind_m+sp] + dc_b1_13[lev]*X[ind_m+2*sp] ; 
			    	    
			    RHS2 = dc_b2_11[lev]*X[ind_m+sp*str_m] + dc_b2_12[lev]*X[ind_m+sp+sp*str_m] + dc_b2_13[lev]*X[ind_m+2*sp+sp*str_m];	     	    
			   
			    RHS3 = dc_b3_11[lev]*X[ind_m+2*sp*str_m] + dc_b3_12[lev]*X[ind_m+sp+2*sp*str_m] + dc_b3_13[lev]*X[ind_m+2*sp+2*sp*str_m]; 	        			         

			   //ans(ind_m) =  level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 ); 
			    
			     ans(ind_m) =  ( RHS1 + RHS2 + RHS3 ); 
			    
			    //cout<<dc_b1_11[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_ny_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n";			      
			  }
			  else if(i==st_inx+sp) //checked 
			  {    	    	    
			    RHS1 = dc_nb1_21[lev]*X[ind_m-sp] + dc_nb1_22[lev]*X[ind_m] + dc_nb1_23[lev]*X[ind_m+sp] + dc_nb1_24[lev]*X[ind_m+2*sp] ;
			       	    
		   	    RHS2 = dc_nb2_21[lev]*X[ind_m-sp+sp*str_m] + dc_nb2_22[lev]*X[ind_m+sp*str_m] + dc_nb2_23[lev]*X[ind_m+sp+sp*str_m] + dc_nb2_24[lev]*X[ind_m+2*sp+sp*str_m] ; 
		   	    
		   	    RHS3 = dc_nb3_21[lev]*X[ind_m-sp+2*sp*str_m] + dc_nb3_22[lev]*X[ind_m+2*sp*str_m] + dc_nb3_23[lev]*X[ind_m+sp+2*sp*str_m] + dc_nb3_24[lev]*X[ind_m+2*sp+2*sp*str_m] ;     	      	    		   	   
		   	    
		   	    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3); 			  
		   	    
		   	    ans(ind_m) = (RHS1 + RHS2 + RHS3); 			  
		   	    
		   	   // cout<<dc_nb1_21[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny_21[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_24[lev]<<"\n";	    
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked 
			  {	    	    	    
			    RHS1 = dc_i1_32[lev]*X[ind_m-2*sp] + dc_i1_33[lev]*X[ind_m-sp] + dc_i1_34[lev]*X[ind_m] + dc_i1_35[lev]*X[ind_m+sp] + dc_i1_36[lev]*X[ind_m+2*sp]; 			    
			    RHS2 = dc_i2_32[lev]*X[ind_m-2*sp+sp*str_m] + dc_i2_33[lev]*X[ind_m-sp+sp*str_m] + dc_i2_34[lev]*X[ind_m+sp*str_m] + dc_i2_35[lev]*X[ind_m+sp+sp*str_m] + dc_i2_36[lev]*X[ind_m+2*sp+sp*str_m]; 
			    
			    RHS3 =  dc_i3_32[lev]*X[ind_m-2*sp+2*sp*str_m] + dc_i3_33[lev]*X[ind_m-sp+2*sp*str_m] + dc_i3_34[lev]*X[ind_m+2*sp*str_m] + dc_i3_35[lev]*X[ind_m+sp+2*sp*str_m] + dc_i3_36[lev]*X[ind_m+2*sp+2*sp*str_m];	     			    
			    
			  //  ans(ind_m) =  level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3); 			    	     	    
			    
			    ans(ind_m) =  (RHS1 + RHS2 + RHS3);
			  }
			  else if(i==en_inx-sp) //Checked 
			  {	    	    	    	    
			    RHS1 = dc_nb1_21[lev]*X[ind_m+sp] + dc_nb1_22[lev]*X[ind_m] + dc_nb1_23[lev]*X[ind_m-sp] + dc_nb1_24[lev]*X[ind_m-2*sp] ; 
			    
			    RHS2 = dc_nb2_21[lev]*X[ind_m+sp+sp*str_m] + dc_nb2_22[lev]*X[ind_m+sp*str_m] + dc_nb2_23[lev]*X[ind_m-sp+sp*str_m] + dc_nb2_24[lev]*X[ind_m-2*sp+sp*str_m] ; 
			    
			    RHS3 = dc_nb3_21[lev]*X[ind_m+sp+2*sp*str_m] + dc_nb3_22[lev]*X[ind_m+2*sp*str_m] + dc_nb3_23[lev]*X[ind_m-sp+2*sp*str_m] + dc_nb3_24[lev]*X[ind_m-2*sp+2*sp*str_m];  	    	    
		   	    
		   	    //cout<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_21[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_22[lev]<<"\n"; 
		   			   	    
		   	    //ans(ind_m) =  level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3);			   
		   	    
		   	    ans(ind_m) =  (RHS1 + RHS2 + RHS3);
			  }
			  else //Checked 
			  {	    	    	    
			    RHS1 = dc_b1_11[lev]*X[ind_m] + dc_b1_12[lev]*X[ind_m-sp] + dc_b1_13[lev]*X[ind_m-2*sp]; 
			    
			    RHS2 = dc_b2_11[lev]*X[ind_m+sp*str_m] + dc_b2_12[lev]*X[ind_m-sp+sp*str_m] + dc_b2_13[lev]*X[ind_m-2*sp+sp*str_m]; 
			    
			    RHS3 = dc_b3_11[lev]*X[ind_m+2*sp*str_m] + dc_b3_12[lev]*X[ind_m-sp+2*sp*str_m] + dc_b3_13[lev]*X[ind_m-2*sp+2*sp*str_m]; 	   	        
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3);
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3);	
			    
			    //cout<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_11[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n";			        
			  }	  	
			}	
			
			/*****************Case-2 (j=2)**********************/ 
	
			if(j==st_iny+sp)
			{
			  if(i==st_inx) //Checked 
			  {	    	    	    	    
			    RHS1 = dc_b4_11[lev]*X[ind_m-sp*str_m] + dc_b4_12[lev]*X[ind_m+sp-sp*str_m] + dc_b4_13[lev]*X[ind_m+2*sp-sp*str_m];  
			    
			    RHS2 = dc_b5_11[lev]*X[ind_m] + dc_b5_12[lev]*X[ind_m+sp] + dc_b5_13[lev]*X[ind_m+2*sp]; 
			    
			    RHS3 = dc_b6_11[lev]*X[ind_m+sp*str_m] + dc_b6_12[lev]*X[ind_m+sp+sp*str_m] + dc_b6_13[lev]*X[ind_m+2*sp+sp*str_m] ; 
			     
			    RHS4 = dc_b7_11[lev]*X[ind_m+2*sp*str_m] + dc_b7_12[lev]*X[ind_m+sp+2*sp*str_m] + dc_b7_13[lev]*X[ind_m+2*sp+2*sp*str_m]; 	    	      
			    
			    //ans(ind_m) =  level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4) ;	 
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4) ;
			    
			    //cout<<dc_b4_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_pb7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n";			       			   		  
			  }
			  else if(i==st_inx+sp) //Checked 
			  {	    	    	    
			    RHS1 = dc_nb4_21[lev]*X[ind_m-sp-sp*str_m] + dc_nb4_22[lev]*X[ind_m-sp*str_m] + dc_nb4_23[lev]*X[ind_m+sp-sp*str_m] + dc_nb4_24[lev]*X[ind_m+2*sp-sp*str_m]; 
			    
			    RHS2 = dc_nb5_21[lev]*X[ind_m-sp] + dc_nb5_22[lev]*X[ind_m] + dc_nb5_23[lev]*X[ind_m+sp] + dc_nb5_24[lev]*X[ind_m+2*sp]; 
			    
			    RHS3 = dc_nb6_21[lev]*X[ind_m-sp+sp*str_m] + dc_nb6_22[lev]*X[ind_m+sp*str_m] + dc_nb6_23[lev]*X[ind_m+sp+sp*str_m] + dc_nb6_24[lev]*X[ind_m+2*sp+sp*str_m]; 
			    
			    RHS4 = dc_nb7_21[lev]*X[ind_m-sp+2*sp*str_m] + dc_nb7_22[lev]*X[ind_m+2*sp*str_m] + dc_nb7_23[lev]*X[ind_m+sp+2*sp*str_m] + dc_nb7_24[lev]*X[ind_m+2*sp+2*sp*str_m] ;
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4);
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4); 
			    
			    //cout<<dc_nb4_21[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_pb7_21[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_pb7_24[lev]<<"\n";			    			    			    	   			  	 
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked
			  {	    	    	    
			    RHS1 = dc_i4_32[lev]*X[ind_m-2*sp-sp*str_m] + dc_i4_33[lev]*X[ind_m-sp-sp*str_m] + dc_i4_34[lev]*X[ind_m-sp*str_m] + dc_i4_35[lev]*X[ind_m+sp-sp*str_m] + dc_i4_36[lev]*X[ind_m+2*sp-sp*str_m]; 
			    
			    RHS2 = dc_i5_32[lev]*X[ind_m-2*sp] + dc_i5_33[lev]*X[ind_m-sp] + dc_i5_34[lev]*X[ind_m] + dc_i5_35[lev]*X[ind_m+sp] + dc_i5_36[lev]*X[ind_m+2*sp]; 
			    
			    RHS3 = dc_i6_32[lev]*X[ind_m-2*sp+sp*str_m] + dc_i6_33[lev]*X[ind_m-sp+sp*str_m] + dc_i6_34[lev]*X[ind_m+sp*str_m] + dc_i6_35[lev]*X[ind_m+sp+sp*str_m] + dc_i6_36[lev]*X[ind_m+2*sp+sp*str_m]; 
			    
			    RHS4 = dc_i7_32[lev]*level[lev].phi_s[ind_m-2*sp+2*sp*str_m] + dc_i7_33[lev]*X[ind_m-sp+2*sp*str_m] + dc_i7_34[lev]*X[ind_m+2*sp*str_m] + dc_i7_35[lev]*X[ind_m+sp+2*sp*str_m] + dc_i7_36[lev]*X[ind_m+2*sp+2*sp*str_m];			       	    
				    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4);
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4);	
			    
			    /*cout<<dc_i4_32[lev]<<"\n"<<dc_i4_33[lev]<<"\n"<<dc_i4_34[lev]<<"\n"<<dc_i4_35[lev]<<"\n"<<dc_i4_36[lev]<<"\n"<<dc_i5_32[lev]<<"\n"<<dc_i5_33[lev]<<"\n"<<dc_i5_34[lev]<<"\n"<<dc_i5_35[lev]<<"\n"<<dc_i5_36[lev]<<"\n"<<dc_i6_32[lev]<<"\n"<<dc_i6_33[lev]<<"\n"<<dc_i6_34[lev]<<"\n"<<dc_i6_35[lev]<<"\n"<<dc_i6_36[lev]<<"\n"<<dc_i7_33[lev]<<"\n"<<dc_i7_34[lev]<<"\n"<<dc_i7_35[lev]<<"\n"<<dc_pb7_33[lev]<<"\n"<<dc_pb7_34[lev]<<"\n"<<dc_pb7_35[lev]<<"\n"; 	
			    	
			    cout<<"--------------------------------------";*/ 			    		    	        	    
			  }
			  else if(i==en_inx-sp) //Checked 
			  {	    	    	    
			    RHS1 = dc_nb4_21[lev]*X[ind_m+sp-sp*str_m] + dc_nb4_22[lev]*X[ind_m-sp*str_m] + dc_nb4_23[lev]*X[ind_m-sp-sp*str_m] + dc_nb4_24[lev]*X[ind_m-2*sp-sp*str_m]; 
			    
			    RHS2 = dc_nb5_21[lev]*X[ind_m+sp] + dc_nb5_22[lev]*X[ind_m] + dc_nb5_23[lev]*X[ind_m-sp] + dc_nb5_24[lev]*X[ind_m-2*sp]; 
			    
			    RHS3 = dc_nb6_21[lev]*X[ind_m+sp+sp*str_m] + dc_nb6_22[lev]*X[ind_m+sp*str_m] + dc_nb6_23[lev]*X[ind_m-sp+sp*str_m] + dc_nb6_24[lev]*X[ind_m-2*sp+sp*str_m]; 
			    
			    RHS4 = dc_nb7_21[lev]*X[ind_m+sp+2*sp*str_m] + dc_nb7_22[lev]*X[ind_m+2*sp*str_m] + dc_nb7_23[lev]*X[ind_m-sp+2*sp*str_m] + dc_nb7_24[lev]*X[ind_m-2*sp+2*sp*str_m]; 			    
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4);	
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4);
			    
			    //cout<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_21[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_pb7_22[lev]<<"\n";       
			  }
			  else //checked 
			  {	    	     	    
			    RHS1 = dc_b4_11[lev]*X[ind_m-sp*str_m] + dc_b4_12[lev]*X[ind_m-sp-sp*str_m] + dc_b4_13[lev]*X[ind_m-2*sp-sp*str_m] ; 	    
			    
			    RHS2 =  dc_b5_11[lev]*X[ind_m] + dc_b5_12[lev]*X[ind_m-sp] + dc_b5_13[lev]*X[ind_m-2*sp]; 
			    
			    RHS3 = dc_b6_11[lev]*X[ind_m+sp*str_m] + dc_b6_12[lev]*X[ind_m-sp+sp*str_m] + dc_b6_13[lev]*X[ind_m-2*sp+sp*str_m];	    
			    
			    RHS4 = dc_b7_11[lev]*X[ind_m+2*sp*str_m] + dc_b7_12[lev]*X[ind_m-sp+2*sp*str_m] + dc_b7_13[lev]*X[ind_m-2*sp+2*sp*str_m]; 			    			    
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4) ;	
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4) ;
			    
			   // cout<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n";		    
			  }  
			}
	
			/*****************Case-3 (j>=3 and j<=ny-3)**********************/ 
	
			if(j>=st_iny+2*sp && j<=en_iny-2*sp)
			{
			  if(i==st_inx) //Checked 
			  {	    	     	    
			    RHS1 = dc_b8_11[lev]*X[ind_m-2*sp*str_m] + dc_b8_12[lev]*X[ind_m+sp-2*sp*str_m] + dc_b8_13[lev]*X[ind_m+2*sp-2*sp*str_m];  
			    			    
			    RHS2 = dc_b9_11[lev]*X[ind_m-sp*str_m] + dc_b9_12[lev]*X[ind_m+sp-sp*str_m] + dc_b9_13[lev]*X[ind_m+2*sp-sp*str_m];
			    
			    RHS3 = dc_b10_11[lev]*X[ind_m] + dc_b10_12[lev]*X[ind_m+sp] + dc_b10_13[lev]*X[ind_m+2*sp] ; 
			    
			    RHS4 = dc_b11_11[lev]*X[ind_m+sp*str_m] + dc_b11_12[lev]*X[ind_m+sp+sp*str_m] + dc_b11_13[lev]*X[ind_m+2*sp+sp*str_m] ; 	    	    
	    
		    	    RHS5 = dc_b12_11[lev]*X[ind_m+2*sp*str_m] + dc_b12_12[lev]*X[ind_m+sp+2*sp*str_m] + dc_b12_13[lev]*X[ind_m+2*sp+2*sp*str_m];	   	
			
			//cout<<dc_b8_11[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b9_11[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b10_11[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b11_11[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b12_11[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_12[lev]<<"\n";
			    	
			    	//cout<<"--------------------------------------------";
			    
			     
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4 + RHS5);	   			      	    
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4 + RHS5);
			    
			  }
			  else if(i==st_inx+sp)  //Checked 
			  {	    	    	    
			    RHS1 = dc_nb8_21[lev]*X[ind_m-sp-2*sp*str_m] + dc_nb8_22[lev]*X[ind_m-2*sp*str_m] + dc_nb8_23[lev]*X[ind_m+sp-2*sp*str_m] + dc_nb8_24[lev]*X[ind_m+2*sp-2*sp*str_m]; 
			    
			    RHS2 = dc_nb9_21[lev]*X[ind_m-sp-sp*str_m] + dc_nb9_22[lev]*X[ind_m-sp*str_m] + dc_nb9_23[lev]*X[ind_m+sp-sp*str_m] + dc_nb9_24[lev]*X[ind_m+2*sp-sp*str_m];
			    
			    RHS3 = dc_nb10_21[lev]*X[ind_m-sp] + dc_nb10_22[lev]*X[ind_m] + dc_nb10_23[lev]*X[ind_m+sp] + dc_nb10_24[lev]*X[ind_m+2*sp]; 
		    			    			    
			    RHS4 = dc_nb11_21[lev]*X[ind_m-sp+sp*str_m] + dc_nb11_22[lev]*X[ind_m+sp*str_m] + dc_nb11_23[lev]*X[ind_m+sp+sp*str_m] + dc_nb11_24[lev]*X[ind_m+2*sp+sp*str_m]; 
			    
			    RHS5 = dc_nb12_21[lev]*X[ind_m-sp+2*sp*str_m] + dc_nb12_22[lev]*X[ind_m+2*sp*str_m] + dc_nb12_23[lev]*X[ind_m+sp+2*sp*str_m] + dc_nb12_24[lev]*X[ind_m+2*sp+2*sp*str_m] ;
			    			    
			   // ans(ind_m) =  level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4 + RHS5);			    
			    
			    ans(ind_m) =  (RHS1 + RHS2 + RHS3 + RHS4 + RHS5);
			    
			   /* cout<<dc_nb8_21[lev]<<"\n"<<dc_nb8_22[lev]<<"\n"<<dc_nb8_23[lev]<<"\n"<<dc_nb9_21[lev]<<"\n"<<dc_nb9_22[lev]<<"\n"<<dc_nb9_23[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb10_21[lev]<<"\n"<<dc_nb10_22[lev]<<"\n"<<dc_nb10_23[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb11_21[lev]<<"\n"<<dc_nb11_22[lev]<<"\n"<<dc_nb11_23[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb12_21[lev]<<"\n"<<dc_nb12_22[lev]<<"\n"<<dc_nb12_23[lev]<<"\n";
			    	
			    cout<<"--------------------------------------------";*/
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)  //Checked 
			  {    	    
			    RHS1 = dc_i8_33[lev]*X[ind_m-sp-2*sp*str_m] + dc_i8_34[lev]*X[ind_m-2*sp*str_m] + dc_i8_35[lev]*X[ind_m+sp-2*sp*str_m]; 
			    
			    RHS2 = dc_i9_32[lev]*X[ind_m-2*sp-sp*str_m] + dc_i9_33[lev]*X[ind_m-sp-sp*str_m] + dc_i9_34[lev]*X[ind_m-sp*str_m] + dc_i9_35[lev]*X[ind_m+sp-sp*str_m] + dc_i9_36[lev]*X[ind_m+2*sp-sp*str_m]; 
			    
			    RHS3 = dc_i10_32[lev]*X[ind_m-2*sp] + dc_i10_33[lev]*X[ind_m-sp] + dc_i10_34[lev]*X[ind_m] + dc_i10_35[lev]*X[ind_m+sp] + dc_i10_36[lev]*X[ind_m+2*sp]; 			    
			    RHS4 = dc_i11_32[lev]*X[ind_m-2*sp+sp*str_m] + dc_i11_33[lev]*X[ind_m-sp+sp*str_m] + dc_i11_34[lev]*X[ind_m+sp*str_m] + dc_i11_35[lev]*X[ind_m+sp+sp*str_m] + dc_i11_36[lev]*X[ind_m+2*sp+sp*str_m]; 
			    
			    RHS5 = dc_i12_33[lev]*X[ind_m-sp+2*sp*str_m] + dc_i12_34[lev]*X[ind_m+2*sp*str_m] + dc_i12_35[lev]*X[ind_m+sp+2*sp*str_m];
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4 + RHS5);
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4 + RHS5);	
			    
			    /*cout<<dc_i8_33[lev]<<"\n"<<dc_i8_34[lev]<<"\n"<<dc_i8_35[lev]<<"\n"<<dc_i9_32[lev]<<"\n"<<dc_i9_33[lev]<<"\n"<<dc_i9_34[lev]<<"\n"<<dc_i9_35[lev]<<"\n"<<dc_i9_36[lev]<<"\n"<<dc_i10_32[lev]<<"\n"<<dc_i10_33[lev]<<"\n"<<dc_i10_34[lev]<<"\n"<<dc_i10_35[lev]<<"\n"<<dc_i10_36[lev]<<"\n"<<dc_i11_32[lev]<<"\n"<<dc_i11_33[lev]<<"\n"<<dc_i11_34[lev]<<"\n"<<dc_i11_35[lev]<<"\n"<<dc_i11_36[lev]<<"\n"<<dc_i12_33[lev]<<"\n"<<dc_i12_34[lev]<<"\n"<<dc_i12_35[lev]<<"\n"; 			    		       
			    
			    cout<<"--------------------------------------------";*/
			  }
			  else if(i==en_inx-sp) //Checked
			  {	    	     	    
			    RHS1 = dc_nb8_21[lev]*X[ind_m+sp-2*sp*str_m] + dc_nb8_22[lev]*X[ind_m-2*sp*str_m] + dc_nb8_23[lev]*X[ind_m-sp-2*sp*str_m] + dc_nb8_24[lev]*X[ind_m-2*sp-2*sp*str_m]; 
			    
			    RHS2 = dc_nb9_21[lev]*X[ind_m+sp-sp*str_m] + dc_nb9_22[lev]*X[ind_m-sp*str_m] + dc_nb9_23[lev]*X[ind_m-sp-sp*str_m] + dc_nb9_24[lev]*X[ind_m-2*sp-sp*str_m]; 
			    
			    RHS3 = dc_nb10_21[lev]*X[ind_m+sp] + dc_nb10_22[lev]*X[ind_m] + dc_nb10_23[lev]*X[ind_m-sp] + dc_nb10_24[lev]*X[ind_m-2*sp]; 
			    
			    RHS4 = dc_nb11_21[lev]*X[ind_m+sp+sp*str_m] + dc_nb11_22[lev]*X[ind_m+sp*str_m] + dc_nb11_23[lev]*X[ind_m-sp+sp*str_m] + dc_nb11_24[lev]*X[ind_m-2*sp+sp*str_m]; 	   			    
			   
			    if(j+2*sp==en_iny)
			    {			    
			    	RHS5 = dc_nb12_22[lev]*X[ind_m+2*sp*str_m] + dc_nb12_23[lev]*X[ind_m-sp+2*sp*str_m] + dc_nb12_24[lev]*X[ind_m-2*sp+2*sp*str_m];	 //lcr
			    }
			    else
			    {
			    	RHS5 = dc_nb12_21[lev]*X[ind_m+sp+2*sp*str_m] + dc_nb12_22[lev]*X[ind_m+2*sp*str_m] + dc_nb12_23[lev]*X[ind_m-sp+2*sp*str_m] + dc_nb12_24[lev]*X[ind_m-2*sp+2*sp*str_m];	 
			    }		    			     
				
                            //cout<<dc_nb8_23[lev]<<"\n"<<dc_nb8_22[lev]<<"\n"<<dc_nb8_21[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_23[lev]<<"\n"<<dc_nb9_22[lev]<<"\n"<<dc_nb9_21[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_23[lev]<<"\n"<<dc_nb10_22[lev]<<"\n"<<dc_nb10_21[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_23[lev]<<"\n"<<dc_nb11_22[lev]<<"\n"<<dc_nb11_21[lev]<<"\n"<<dc_nb12_23[lev]<<"\n"<<dc_nb12_22[lev]<<"\n"<<dc_nb12_21[lev]<<"\n"; 	
			    	
			    	//cout<<"----------------------------------------------\n";		    	
			    
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4 + RHS5); 			   
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4 + RHS5); 			    
			  }
			  else //Checked 
			  {    	    	    
			    RHS1 = dc_b8_11[lev]*X[ind_m-2*sp*str_m] + dc_b8_12[lev]*X[ind_m-sp-2*sp*str_m] + dc_b8_13[lev]*X[ind_m-2*sp-2*sp*str_m]; 
			    
			    RHS2 = dc_b9_11[lev]*X[ind_m-sp*str_m] + dc_b9_12[lev]*X[ind_m-sp-sp*str_m] + dc_b9_13[lev]*X[ind_m-2*sp-sp*str_m]; 
			    
			    RHS3 = dc_b10_11[lev]*X[ind_m] + dc_b10_12[lev]*X[ind_m-sp] + dc_b10_13[lev]*X[ind_m-2*sp]; 
			    
			    RHS4 = dc_b11_11[lev]*X[ind_m+sp*str_m] + dc_b11_12[lev]*X[ind_m-sp+sp*str_m] + dc_b11_13[lev]*X[ind_m-2*sp+sp*str_m]; 	    	    
			    
			    if(j+2*sp==en_iny)
			    {
			    	RHS5 = dc_b12_12[lev]*X[ind_m-sp+2*sp*str_m] + dc_b12_13[lev]*X[ind_m-2*sp+2*sp*str_m];  //lcr
			    }
			    else
			    {
			    	RHS5 = dc_b12_11[lev]*X[ind_m+2*sp*str_m] + dc_b12_12[lev]*X[ind_m-sp+2*sp*str_m] + dc_b12_13[lev]*X[ind_m-2*sp+2*sp*str_m];
			    }
			    
			    	//cout<<dc_b8_12[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b8_11[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b9_11[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_11[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_11[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_11[lev]<<"\n";
			    	
			    	//cout<<"---------------------------------------\n";
			    
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4 + RHS5);			   
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4 + RHS5); 			  		  
			  }	  	  
			  
			}
	
			/*****************Case-4(j=ny-2)**********************/ 
	
			if(j==en_iny-sp)
			{
			  if(i==st_inx)  //checked 
			  {	     	    
			    RHS1 = dc_b4_11[lev]*X[ind_m+sp*str_m] + dc_b4_12[lev]*X[ind_m+sp+sp*str_m] + dc_b4_13[lev]*X[ind_m+2*sp+sp*str_m] ; 

			    RHS2 = dc_b5_11[lev]*X[ind_m] + dc_b5_12[lev]*X[ind_m+sp] + dc_b5_13[lev]*X[ind_m+2*sp]; 

			    RHS3 = dc_b6_11[lev]*X[ind_m-sp*str_m] + dc_b6_12[lev]*X[ind_m+sp-sp*str_m] + dc_b6_13[lev]*X[ind_m+2*sp-sp*str_m]; 
			    
			    RHS4 = dc_b7_11[lev]*X[ind_m-2*sp*str_m] + dc_b7_12[lev]*X[ind_m+sp-2*sp*str_m] + dc_b7_13[lev]*X[ind_m+2*sp-2*sp*str_m] ;			    			    
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4);			     	            
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4);
			    			   			    
			    //cout<<dc_pb7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b4_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"; 
			    
			  }
			  else if(i==st_inx+sp) //Checked
			  {	    	    	    
			    RHS1 = dc_nb4_21[lev]*X[ind_m-sp+sp*str_m] + dc_nb4_22[lev]*X[ind_m+sp*str_m] + dc_nb4_23[lev]*X[ind_m+sp+sp*str_m] + dc_nb4_24[lev]*X[ind_m+2*sp+sp*str_m]; 
			    
			    RHS2 = dc_nb5_21[lev]*X[ind_m-sp] + dc_nb5_22[lev]*X[ind_m] + dc_nb5_23[lev]*X[ind_m+sp] + dc_nb5_24[lev]*X[ind_m+2*sp]; 
			    
			    RHS3 = dc_nb6_21[lev]*X[ind_m-sp-sp*str_m] + dc_nb6_22[lev]*X[ind_m-sp*str_m] + dc_nb6_23[lev]*X[ind_m+sp-sp*str_m] + dc_nb6_24[lev]*X[ind_m+2*sp-sp*str_m]; 
			    
			    RHS4 = dc_nb7_21[lev]*X[ind_m-sp-2*sp*str_m] + dc_nb7_22[lev]*X[ind_m-2*sp*str_m] + dc_nb7_23[lev]*X[ind_m+sp-2*sp*str_m] + dc_nb7_24[lev]*X[ind_m+2*sp-2*sp*str_m]; 	    
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4);
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4); 
			    
			    //cout<<dc_pb7_21[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb4_21[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"; 			    			     	    	    
			  }
			  else if(i>st_inx+sp && i<en_inx-sp) //Checked 
			  {	    			  		
			  		  			  	
			  	if(i+2*sp==en_inx)
			  	{			  
			  		RHS1 = dc_i4_32[lev]*X[ind_m-2*sp+sp*str_m] + dc_i4_33[lev]*X[ind_m-sp+sp*str_m] + dc_i4_34[lev]*X[ind_m+sp*str_m] + dc_i4_35[lev]*X[ind_m+sp+sp*str_m] ; //lcr 				  	
			  	}
			  	else
			  	{
			  		RHS1 = dc_i4_32[lev]*X[ind_m-2*sp+sp*str_m] + dc_i4_33[lev]*X[ind_m-sp+sp*str_m] + dc_i4_34[lev]*X[ind_m+sp*str_m] + dc_i4_35[lev]*X[ind_m+sp+sp*str_m] + dc_i4_36[lev]*X[ind_m+2*sp+sp*str_m];			  	
			  	}		  	
			  		
			  	/*cout<<dc_pb7_33[lev]<<"\n"<<dc_pb7_34[lev]<<"\n"<<dc_pb7_35[lev]<<"\n"<<dc_i7_33[lev]<<"\n"<<dc_i7_34[lev]<<"\n"<<dc_i7_35[lev]<<"\n"<<dc_i6_32[lev]<<"\n"<<dc_i6_33[lev]<<"\n"<<dc_i6_34[lev]<<"\n"<<dc_i6_35[lev]<<"\n"<<dc_i6_36[lev]<<"\n"<<dc_i5_32[lev]<<"\n"<<dc_i5_33[lev]<<"\n"<<dc_i5_34[lev]<<"\n"<<dc_i5_35[lev]<<"\n"<<dc_i5_36[lev]<<"\n"<<dc_i4_32[lev]<<"\n"<<dc_i4_33[lev]<<"\n"<<dc_i4_34[lev]<<"\n"<<dc_i4_35[lev]<<"\n"<<dc_i4_36[lev]<<"\n";		  					  	 	
			  		
			  	cout<<"ind_m= "<<ind_m<<"\n"; 
			  		
			  	cout<<"----------------------------------------------------\n"; */
			  				  	    	    			    
			    
			    	RHS2 = dc_i5_32[lev]*X[ind_m-2*sp] + dc_i5_33[lev]*X[ind_m-sp] + dc_i5_34[lev]*X[ind_m] + dc_i5_35[lev]*X[ind_m+sp] + dc_i5_36[lev]*X[ind_m+2*sp]; 
			    
			    	RHS3 = dc_i6_32[lev]*X[ind_m-2*sp-sp*str_m] + dc_i6_33[lev]*X[ind_m-sp-sp*str_m] + dc_i6_34[lev]*X[ind_m-sp*str_m] + dc_i6_35[lev]*X[ind_m+sp-sp*str_m] + dc_i6_36[lev]*X[ind_m+2*sp-sp*str_m]; 
			    
			    	RHS4 = dc_i7_32[lev]*X[ind_m-2*sp-2*sp*str_m] + dc_i7_33[lev]*X[ind_m-sp-2*sp*str_m] + dc_i7_34[lev]*X[ind_m-2*sp*str_m] + dc_i7_35[lev]*X[ind_m+sp-2*sp*str_m] + dc_i7_36[lev]*X[ind_m+2*sp-2*sp*str_m] ;			    			    
			    
			    	//ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 + RHS4 ) ;
			    	
			    	ans(ind_m) = ( RHS1 + RHS2 + RHS3 + RHS4 ) ;			     	    	 
			  }
			  else if(i==en_inx-sp) //checked 
			  {	    	    	    
			    	RHS1 =  dc_nb4_22[lev]*X[ind_m+sp*str_m] + dc_nb4_23[lev]*X[ind_m-sp+sp*str_m] + dc_nb4_24[lev]*X[ind_m-2*sp+sp*str_m]; //lcr 
			    
			    	RHS2 = dc_nb5_21[lev]*X[ind_m+sp] + dc_nb5_22[lev]*X[ind_m] + dc_nb5_23[lev]*X[ind_m-sp] + dc_nb5_24[lev]*X[ind_m-2*sp]; 
			    
			    	RHS3 = dc_nb6_21[lev]*X[ind_m+sp-sp*str_m] + dc_nb6_22[lev]*X[ind_m-sp*str_m] + dc_nb6_23[lev]*X[ind_m-sp-sp*str_m] + dc_nb6_24[lev]*X[ind_m-2*sp-sp*str_m]; 
			    
			    	RHS4 = dc_nb7_21[lev]*X[ind_m+sp-2*sp*str_m] + dc_nb7_22[lev]*X[ind_m-2*sp*str_m] + dc_nb7_23[lev]*X[ind_m-sp-2*sp*str_m] +  dc_nb7_24[lev]*X[ind_m-2*sp-2*sp*str_m];	       
			    
			    	//ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 + RHS4 ) ;
			    	
			    	ans(ind_m) = ( RHS1 + RHS2 + RHS3 + RHS4 ) ;	
			    				    
			    	//cout<<"ind_m= "<<ind_m<<"\n";
			    
			    	//cout<<dc_pb7_23[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_21[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_22[lev]<<"\n";   
			  }
			  else //Checked 
			  {	    	    	    
			    RHS1 = dc_b4_12[lev]*X[ind_m-sp+sp*str_m] + dc_b4_13[lev]*X[ind_m-2*sp+sp*str_m];   	    //lcr
			    
			    RHS2 = dc_b5_11[lev]*X[ind_m] + dc_b5_12[lev]*X[ind_m-sp] + dc_b5_13[lev]*X[ind_m-2*sp]; 
			    
			    RHS3 = dc_b6_11[lev]*X[ind_m-sp*str_m] + dc_b6_12[lev]*X[ind_m-sp-sp*str_m] + dc_b6_13[lev]*X[ind_m-2*sp-sp*str_m]; 
			    
			    RHS4 = dc_b7_11[lev]*X[ind_m-2*sp*str_m] + dc_b7_12[lev]*X[ind_m-sp-2*sp*str_m] +  dc_b7_13[lev]*X[ind_m-2*sp-2*sp*str_m] ;      				   
			    
			    //ans(ind_m) =  level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 + RHS4 ) ;
			    
			    ans(ind_m) =  ( RHS1 + RHS2 + RHS3 + RHS4 ) ;	
			    
			   // cout<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n";			    		    			    	    
			  }  
			}
		
			/****************************************Case-5(j=en_iny)******************************************/ 
		
			if(j==en_iny)
			{
			  if(i==st_inx) 
			  {			  	    	    	    
			    RHS1 = dc_b1_11[lev]*X[ind_m] + dc_b1_12[lev]*X[ind_m+sp] + dc_b1_13[lev]*X[ind_m+2*sp]; 
			    
			    RHS2 = dc_b2_11[lev]*X[ind_m-sp*str_m] + dc_b2_12[lev]*X[ind_m+sp-sp*str_m] + dc_b2_13[lev]*X[ind_m+2*sp-sp*str_m]; 
			    
			    RHS3 = dc_b3_11[lev]*X[ind_m-2*sp*str_m] + dc_b3_12[lev]*X[ind_m+sp-2*sp*str_m] + dc_b3_13[lev]*X[ind_m+2*sp-2*sp*str_m]; 	        			       	    

			    //ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 );
			    
			    ans(ind_m) = ( RHS1 + RHS2 + RHS3 );	 			    
			    
			    //cout<<dc_ny_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b1_11[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"; 			    
			  }
			  else if(i==st_inx+sp) //Checked 
			  {	    	   	    
			    RHS1 = dc_nb1_21[lev]*X[ind_m-sp] + dc_nb1_22[lev]*X[ind_m] + dc_nb1_23[lev]*X[ind_m+sp] + dc_nb1_24[lev]*X[ind_m+2*sp] ; 
			    
			    RHS2 = dc_nb2_21[lev]*X[ind_m-sp-sp*str_m] + dc_nb2_22[lev]*X[ind_m-sp*str_m] + dc_nb2_23[lev]*X[ind_m+sp-sp*str_m] + dc_nb2_24[lev]*X[ind_m+2*sp-sp*str_m]; 
			    
			    RHS3 = dc_nb3_21[lev]*X[ind_m-sp-2*sp*str_m] + dc_nb3_22[lev]*X[ind_m-2*sp*str_m] + dc_nb3_23[lev]*X[ind_m+sp-2*sp*str_m] + dc_nb3_24[lev]*X[ind_m+2*sp-2*sp*str_m];
		    	       	        	    	   	    
		   	    //ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 ) ;
		   	    
		   	    ans(ind_m) = ( RHS1 + RHS2 + RHS3 ) ;    	    
		   	    
		   	    /*cout<<"ind_m = "<<ind_m<<"\n"; 		   	       	    
		   	    
		   	    cout<<dc_ny_21[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb1_21[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_24[lev]<<"\n";*/		   	    		   	        
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked
			  {	    	    	    			   	     	
			     	if(i+2*sp==en_inx)
			     	{			     				     	
			     		RHS1 = dc_i1_32[lev]*X[ind_m-2*sp] + dc_i1_33[lev]*X[ind_m-sp] + dc_i1_34[lev]*X[ind_m] + dc_i1_35[lev]*X[ind_m+sp]; //lcr 
			     	}
			     	else
			     	{			     	
			     		RHS1 = dc_i1_32[lev]*X[ind_m-2*sp] + dc_i1_33[lev]*X[ind_m-sp] + dc_i1_34[lev]*X[ind_m] + dc_i1_35[lev]*X[ind_m+sp] + dc_i1_36[lev]*X[ind_m+2*sp];
			     	}			     	
			     	
			     	/*cout<<"ind_m= "<<ind_m<<"\n"; 
			     		
				cout<<dc_ny_32[lev]<<"\n"<<dc_ny_33[lev]<<"\n"<<dc_ny_34[lev]<<"\n"<<dc_ny_35[lev]<<"\n"<<dc_ny_36[lev]<<"\n"<<dc_ny1_33[lev]<<"\n"<<dc_ny1_34[lev]<<"\n"<<dc_ny1_35[lev]<<"\n"<<dc_i3_33[lev]<<"\n"<<dc_i3_34[lev]<<"\n"<<dc_i3_35[lev]<<"\n"<<dc_i2_32[lev]<<"\n"<<dc_i2_33[lev]<<"\n"<<dc_i2_34[lev]<<"\n"<<dc_i2_35[lev]<<"\n"<<dc_i2_36[lev]<<"\n"<<dc_i1_32[lev]<<"\n"<<dc_i1_33[lev]<<"\n"<<dc_i1_34[lev]<<"\n"<<dc_i1_35[lev]<<"\n"<<dc_i1_36[lev]<<"\n";
				
				cout<<"---------------------------------------------------\n";*/ 
			     						  			     
			     
			     RHS2 = dc_i2_32[lev]*X[ind_m-2*sp-sp*str_m] + dc_i2_33[lev]*X[ind_m-sp-sp*str_m] + dc_i2_34[lev]*X[ind_m-sp*str_m] + dc_i2_35[lev]*X[ind_m+sp-sp*str_m] + dc_i2_36[lev]*X[ind_m+2*sp-sp*str_m]; 
			     
			     RHS3 =  dc_i3_32[lev]*X[ind_m-2*sp-2*sp*str_m] + dc_i3_33[lev]*X[ind_m-sp-2*sp*str_m] + dc_i3_34[lev]*X[ind_m-2*sp*str_m] + dc_i3_35[lev]*X[ind_m+sp-2*sp*str_m] + dc_i3_36[lev]*X[ind_m+2*sp-2*sp*str_m];   			        
			     
			     //ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 );			    
			     
			     ans(ind_m) = ( RHS1 + RHS2 + RHS3 );
			  }
			  else if(i==en_inx-sp)  //Checked
			  {    	    	   	    			   	        
			    RHS1 = dc_nb1_22[lev]*X[ind_m] + dc_nb1_23[lev]*X[ind_m-sp] + dc_nb1_24[lev]*X[ind_m-2*sp]; 	//lcr 		    			    
			    
			    RHS2 = dc_nb2_21[lev]*X[ind_m+sp-sp*str_m] + dc_nb2_22[lev]*X[ind_m-sp*str_m] + dc_nb2_23[lev]*X[ind_m-sp-sp*str_m] + dc_nb2_24[lev]*X[ind_m-2*sp-sp*str_m]; 
			    
			    RHS3 = dc_nb3_21[lev]*X[ind_m+sp-2*sp*str_m] + dc_nb3_22[lev]*X[ind_m-2*sp*str_m] + dc_nb3_23[lev]*X[ind_m-sp-2*sp*str_m] + dc_nb3_24[lev]*X[ind_m-2*sp-2*sp*str_m] ;	        
		   	    
		   	    //ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 );
		   	    
		   	    ans(ind_m) = ( RHS1 + RHS2 + RHS3 );
		   	    
		   	    /*cout<<"ind_m= "<<ind_m<<"\n"; 			    	    
		   	    
		   	    cout<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_21[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"; */
			  }
			  else //Checked 
			  {	  	     	     
			     RHS1 = dc_b1_12[lev]*X[ind_m-sp] + dc_b1_13[lev]*X[ind_m-2*sp]; //lcr 
			     
			     RHS2 = dc_b2_11[lev]*X[ind_m-sp*str_m] + dc_b2_12[lev]*X[ind_m-sp-sp*str_m] + dc_b2_13[lev]*X[ind_m-2*sp-sp*str_m]; 
			     
			     RHS3 = dc_b3_11[lev]*X[ind_m-2*sp*str_m] + dc_b3_12[lev]*X[ind_m-sp-2*sp*str_m] + dc_b3_13[lev]*X[ind_m-2*sp-2*sp*str_m];    	    			     			     	    
			     //ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 ); //This mulvec part belongs to the pinned point. Hence, it is trimmed out and not involved in computation at all. 
			     
			     ans(ind_m) = ( RHS1 + RHS2 + RHS3 ); 			     			     
			  }	  	
			}											
		
		}	
	}

	subans = ans.submat(0,0,tot_p_sol-2,0); 
	return subans; 
}

/***************************************************************************************/
vec mulvec(vec X, mg_grid * level, int lev)
{ 	
	int sp=pow(2,lev), ind; 
	
	int nx_sol = (nx/sp)-1, ny_sol = (ny/sp)-1; //No of points for the computational domain after compatibility condition 
	int tot_p_sol = (nx_sol)*(ny_sol); //No.of points that are solved for directly. The top and the farthest right lines of points are  	
					   //removed due to periodicity	
					   
	vec ans(tot_p_sol), subans(tot_p_sol-1); 				   
			
	int i_m, j_m, str_m=nx_sol, ind_m; 
	
	int st_inx = sp, en_inx = nx - sp;    
  	int st_iny = sp, en_iny = ny - sp;  
  	
	double RHS1, RHS2, RHS3, RHS4, RHS5, RHS6; 
	
	for(int j=st_iny;j<=en_iny;j=j+sp)
	{
		for(int i=st_inx;i<=en_inx;i=i+sp)
		{
			ind = i + j*str_x;	
			
			i_m = i/sp-1; 
			j_m = j/sp-1; 
			
			ind_m = i_m + j_m*str_m; 
			
			/*******************Case-1 j=1***********************/ 
			if(j==st_iny)
			{
			  if(i==st_inx) //Checked 
			  {	        	    	    
			    RHS1 = dc_b1_11[lev]*X[ind_m] + dc_b1_12[lev]*X[ind_m+1] + dc_b1_13[lev]*X[ind_m+2] ; 
			    	    
			    RHS2 = dc_b2_11[lev]*X[ind_m+str_m] + dc_b2_12[lev]*X[ind_m+1+str_m] + dc_b2_13[lev]*X[ind_m+2+str_m];	     	    
			   
			    RHS3 = dc_b3_11[lev]*X[ind_m+2*str_m] + dc_b3_12[lev]*X[ind_m+1+2*str_m] + dc_b3_13[lev]*X[ind_m+2+2*str_m]; 	        			         

			   //ans(ind_m) =  level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 ); 
			    
			     ans(ind_m) =  ( RHS1 + RHS2 + RHS3 ); 
			    
			    //cout<<dc_b1_11[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_ny_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n";			      
			  }
			  else if(i==st_inx+sp) //checked 
			  {    	    	    
			    RHS1 = dc_nb1_21[lev]*X[ind_m-1] + dc_nb1_22[lev]*X[ind_m] + dc_nb1_23[lev]*X[ind_m+1] + dc_nb1_24[lev]*X[ind_m+2] ;
			       	    
		   	    RHS2 = dc_nb2_21[lev]*X[ind_m-1+str_m] + dc_nb2_22[lev]*X[ind_m+str_m] + dc_nb2_23[lev]*X[ind_m+1+str_m] + dc_nb2_24[lev]*X[ind_m+2+str_m] ; 
		   	    
		   	    RHS3 = dc_nb3_21[lev]*X[ind_m-1+2*str_m] + dc_nb3_22[lev]*X[ind_m+2*str_m] + dc_nb3_23[lev]*X[ind_m+1+2*str_m] + dc_nb3_24[lev]*X[ind_m+2+2*str_m] ;     	      	    		   	   
		   	    
		   	    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3); 			  
		   	    
		   	    ans(ind_m) = (RHS1 + RHS2 + RHS3); 			  
		   	    
		   	   // cout<<dc_nb1_21[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny_21[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_24[lev]<<"\n";	    
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked 
			  {	    	    	    
			    RHS1 = dc_i1_32[lev]*X[ind_m-2] + dc_i1_33[lev]*X[ind_m-1] + dc_i1_34[lev]*X[ind_m] + dc_i1_35[lev]*X[ind_m+1] + dc_i1_36[lev]*X[ind_m+2]; 			    
			    RHS2 = dc_i2_32[lev]*X[ind_m-2+str_m] + dc_i2_33[lev]*X[ind_m-1+str_m] + dc_i2_34[lev]*X[ind_m+str_m] + dc_i2_35[lev]*X[ind_m+1+str_m] + dc_i2_36[lev]*X[ind_m+2+str_m]; 
			    
			    RHS3 =  dc_i3_32[lev]*X[ind_m-2+2*str_m] + dc_i3_33[lev]*X[ind_m-1+2*str_m] + dc_i3_34[lev]*X[ind_m+2*str_m] + dc_i3_35[lev]*X[ind_m+1+2*str_m] + dc_i3_36[lev]*X[ind_m+2+2*str_m];	     			    
			    
			  //  ans(ind_m) =  level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3); 			    	     	    
			    
			    ans(ind_m) =  (RHS1 + RHS2 + RHS3);
			  }
			  else if(i==en_inx-sp) //Checked 
			  {	    	    	    	    
			    RHS1 = dc_nb1_21[lev]*X[ind_m+1] + dc_nb1_22[lev]*X[ind_m] + dc_nb1_23[lev]*X[ind_m-1] + dc_nb1_24[lev]*X[ind_m-2] ; 
			    
			    RHS2 = dc_nb2_21[lev]*X[ind_m+1+str_m] + dc_nb2_22[lev]*X[ind_m+str_m] + dc_nb2_23[lev]*X[ind_m-1+str_m] + dc_nb2_24[lev]*X[ind_m-2+str_m] ; 
			    
			    RHS3 = dc_nb3_21[lev]*X[ind_m+1+2*str_m] + dc_nb3_22[lev]*X[ind_m+2*str_m] + dc_nb3_23[lev]*X[ind_m-1+2*str_m] + dc_nb3_24[lev]*X[ind_m-2+2*str_m];  	    	    
		   	    
		   	    //cout<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_21[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_22[lev]<<"\n"; 
		   			   	    
		   	    //ans(ind_m) =  level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3);			   
		   	    
		   	    ans(ind_m) =  (RHS1 + RHS2 + RHS3);
			  }
			  else //Checked 
			  {	    	    	    
			    RHS1 = dc_b1_11[lev]*X[ind_m] + dc_b1_12[lev]*X[ind_m-1] + dc_b1_13[lev]*X[ind_m-2]; 
			    
			    RHS2 = dc_b2_11[lev]*X[ind_m+str_m] + dc_b2_12[lev]*X[ind_m-1+str_m] + dc_b2_13[lev]*X[ind_m-2+str_m]; 
			    
			    RHS3 = dc_b3_11[lev]*X[ind_m+2*str_m] + dc_b3_12[lev]*X[ind_m-1+2*str_m] + dc_b3_13[lev]*X[ind_m-2+2*str_m]; 	   	        
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3);
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3);	
			    
			    //cout<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_11[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n";			        
			  }	  	
			}	
			
			/*****************Case-2 (j=2)**********************/ 
	
			if(j==st_iny+sp)
			{
			  if(i==st_inx) //Checked 
			  {	    	    	    	    
			    RHS1 = dc_b4_11[lev]*X[ind_m-str_m] + dc_b4_12[lev]*X[ind_m+1-str_m] + dc_b4_13[lev]*X[ind_m+2-str_m];  
			    
			    RHS2 = dc_b5_11[lev]*X[ind_m] + dc_b5_12[lev]*X[ind_m+1] + dc_b5_13[lev]*X[ind_m+2]; 
			    
			    RHS3 = dc_b6_11[lev]*X[ind_m+str_m] + dc_b6_12[lev]*X[ind_m+1+str_m] + dc_b6_13[lev]*X[ind_m+2+str_m] ; 
			     
			    RHS4 = dc_b7_11[lev]*X[ind_m+2*str_m] + dc_b7_12[lev]*X[ind_m+1+2*str_m] + dc_b7_13[lev]*X[ind_m+2+2*str_m]; 	    	      
			    
			    //ans(ind_m) =  level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4) ;	 
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4) ;
			    
			    //cout<<dc_b4_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_pb7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n";			       			   		  
			  }
			  else if(i==st_inx+sp) //Checked 
			  {	    	    	    
			    RHS1 = dc_nb4_21[lev]*X[ind_m-1-str_m] + dc_nb4_22[lev]*X[ind_m-str_m] + dc_nb4_23[lev]*X[ind_m+1-str_m] + dc_nb4_24[lev]*X[ind_m+2-str_m]; 
			    
			    RHS2 = dc_nb5_21[lev]*X[ind_m-1] + dc_nb5_22[lev]*X[ind_m] + dc_nb5_23[lev]*X[ind_m+1] + dc_nb5_24[lev]*X[ind_m+2]; 
			    
			    RHS3 = dc_nb6_21[lev]*X[ind_m-1+str_m] + dc_nb6_22[lev]*X[ind_m+str_m] + dc_nb6_23[lev]*X[ind_m+1+str_m] + dc_nb6_24[lev]*X[ind_m+2+str_m]; 
			    
			    RHS4 = dc_nb7_21[lev]*X[ind_m-1+2*str_m] + dc_nb7_22[lev]*X[ind_m+2*str_m] + dc_nb7_23[lev]*X[ind_m+1+2*str_m] + dc_nb7_24[lev]*X[ind_m+2+2*str_m] ;
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4);
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4); 
			    
			    //cout<<dc_nb4_21[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_pb7_21[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_pb7_24[lev]<<"\n";			    			    			    	   			  	 
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked
			  {	    	    	    
			    RHS1 = dc_i4_32[lev]*X[ind_m-2-str_m] + dc_i4_33[lev]*X[ind_m-1-str_m] + dc_i4_34[lev]*X[ind_m-str_m] + dc_i4_35[lev]*X[ind_m+1-str_m] + dc_i4_36[lev]*X[ind_m+2-str_m]; 
			    
			    RHS2 = dc_i5_32[lev]*X[ind_m-2] + dc_i5_33[lev]*X[ind_m-1] + dc_i5_34[lev]*X[ind_m] + dc_i5_35[lev]*X[ind_m+1] + dc_i5_36[lev]*X[ind_m+2]; 
			    
			    RHS3 = dc_i6_32[lev]*X[ind_m-2+str_m] + dc_i6_33[lev]*X[ind_m-1+str_m] + dc_i6_34[lev]*X[ind_m+str_m] + dc_i6_35[lev]*X[ind_m+1+str_m] + dc_i6_36[lev]*X[ind_m+2+str_m]; 
			    
			    RHS4 = dc_i7_32[lev]*level[lev].phi_s[ind_m-2+2*str_m] + dc_i7_33[lev]*X[ind_m-1+2*str_m] + dc_i7_34[lev]*X[ind_m+2*str_m] + dc_i7_35[lev]*X[ind_m+1+2*str_m] + dc_i7_36[lev]*X[ind_m+2+2*str_m];			       	    
				    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4);
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4);	
			    
			    /*cout<<dc_i4_32[lev]<<"\n"<<dc_i4_33[lev]<<"\n"<<dc_i4_34[lev]<<"\n"<<dc_i4_35[lev]<<"\n"<<dc_i4_36[lev]<<"\n"<<dc_i5_32[lev]<<"\n"<<dc_i5_33[lev]<<"\n"<<dc_i5_34[lev]<<"\n"<<dc_i5_35[lev]<<"\n"<<dc_i5_36[lev]<<"\n"<<dc_i6_32[lev]<<"\n"<<dc_i6_33[lev]<<"\n"<<dc_i6_34[lev]<<"\n"<<dc_i6_35[lev]<<"\n"<<dc_i6_36[lev]<<"\n"<<dc_i7_33[lev]<<"\n"<<dc_i7_34[lev]<<"\n"<<dc_i7_35[lev]<<"\n"<<dc_pb7_33[lev]<<"\n"<<dc_pb7_34[lev]<<"\n"<<dc_pb7_35[lev]<<"\n"; 	
			    	
			    cout<<"--------------------------------------";*/ 			    		    	        	    
			  }
			  else if(i==en_inx-sp) //Checked 
			  {	    	    	    
			    RHS1 = dc_nb4_21[lev]*X[ind_m+1-str_m] + dc_nb4_22[lev]*X[ind_m-str_m] + dc_nb4_23[lev]*X[ind_m-1-str_m] + dc_nb4_24[lev]*X[ind_m-2-str_m]; 
			    
			    RHS2 = dc_nb5_21[lev]*X[ind_m+1] + dc_nb5_22[lev]*X[ind_m] + dc_nb5_23[lev]*X[ind_m-1] + dc_nb5_24[lev]*X[ind_m-2]; 
			    
			    RHS3 = dc_nb6_21[lev]*X[ind_m+1+str_m] + dc_nb6_22[lev]*X[ind_m+str_m] + dc_nb6_23[lev]*X[ind_m-1+str_m] + dc_nb6_24[lev]*X[ind_m-2+str_m]; 
			    
			    RHS4 = dc_nb7_21[lev]*X[ind_m+1+2*str_m] + dc_nb7_22[lev]*X[ind_m+2*str_m] + dc_nb7_23[lev]*X[ind_m-1+2*str_m] + dc_nb7_24[lev]*X[ind_m-2+2*str_m]; 			    
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4);	
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4);
			    
			    //cout<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_21[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_pb7_22[lev]<<"\n";       
			  }
			  else //checked 
			  {	    	     	    
			    RHS1 = dc_b4_11[lev]*X[ind_m-str_m] + dc_b4_12[lev]*X[ind_m-1-str_m] + dc_b4_13[lev]*X[ind_m-2-str_m] ; 	    
			    
			    RHS2 =  dc_b5_11[lev]*X[ind_m] + dc_b5_12[lev]*X[ind_m-1] + dc_b5_13[lev]*X[ind_m-2]; 
			    
			    RHS3 = dc_b6_11[lev]*X[ind_m+str_m] + dc_b6_12[lev]*X[ind_m-1+str_m] + dc_b6_13[lev]*X[ind_m-2+str_m];	    
			    
			    RHS4 = dc_b7_11[lev]*X[ind_m+2*str_m] + dc_b7_12[lev]*X[ind_m-1+2*str_m] + dc_b7_13[lev]*X[ind_m-2+2*str_m]; 			    			    
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4) ;	
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4) ;
			    
			   // cout<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n";		    
			  }  
			}
	
			/*****************Case-3 (j>=3 and j<=ny-3)**********************/ 
	
			if(j>=st_iny+2*sp && j<=en_iny-2*sp)
			{
			  if(i==st_inx) //Checked 
			  {	    	     	    
			    RHS1 = dc_b8_11[lev]*X[ind_m-2*str_m] + dc_b8_12[lev]*X[ind_m+1-2*str_m] + dc_b8_13[lev]*X[ind_m+2-2*str_m];  
			    			    
			    RHS2 = dc_b9_11[lev]*X[ind_m-str_m] + dc_b9_12[lev]*X[ind_m+1-str_m] + dc_b9_13[lev]*X[ind_m+2-str_m];
			    
			    RHS3 = dc_b10_11[lev]*X[ind_m] + dc_b10_12[lev]*X[ind_m+1] + dc_b10_13[lev]*X[ind_m+2] ; 
			    
			    RHS4 = dc_b11_11[lev]*X[ind_m+str_m] + dc_b11_12[lev]*X[ind_m+1+str_m] + dc_b11_13[lev]*X[ind_m+2+str_m] ; 	    	    
	    
		    	    RHS5 = dc_b12_11[lev]*X[ind_m+2*str_m] + dc_b12_12[lev]*X[ind_m+1+2*str_m] + dc_b12_13[lev]*X[ind_m+2+2*str_m];	   	
			
			//cout<<dc_b8_11[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b9_11[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b10_11[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b11_11[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b12_11[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_12[lev]<<"\n";
			    	
			    	//cout<<"--------------------------------------------";
			    
			     
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4 + RHS5);	   			      	    
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4 + RHS5);
			    
			  }
			  else if(i==st_inx+sp)  //Checked 
			  {	    	    	    
			    RHS1 = dc_nb8_21[lev]*X[ind_m-1-2*str_m] + dc_nb8_22[lev]*X[ind_m-2*str_m] + dc_nb8_23[lev]*X[ind_m+1-2*str_m] + dc_nb8_24[lev]*X[ind_m+2-2*str_m]; 
			    
			    RHS2 = dc_nb9_21[lev]*X[ind_m-1-str_m] + dc_nb9_22[lev]*X[ind_m-str_m] + dc_nb9_23[lev]*X[ind_m+1-str_m] + dc_nb9_24[lev]*X[ind_m+2-str_m];
			    
			    RHS3 = dc_nb10_21[lev]*X[ind_m-1] + dc_nb10_22[lev]*X[ind_m] + dc_nb10_23[lev]*X[ind_m+1] + dc_nb10_24[lev]*X[ind_m+2]; 
		    			    			    
			    RHS4 = dc_nb11_21[lev]*X[ind_m-1+str_m] + dc_nb11_22[lev]*X[ind_m+str_m] + dc_nb11_23[lev]*X[ind_m+1+str_m] + dc_nb11_24[lev]*X[ind_m+2+str_m]; 
			    
			    RHS5 = dc_nb12_21[lev]*X[ind_m-1+2*str_m] + dc_nb12_22[lev]*X[ind_m+2*str_m] + dc_nb12_23[lev]*X[ind_m+1+2*str_m] + dc_nb12_24[lev]*X[ind_m+2+2*str_m] ;
			    			    
			   // ans(ind_m) =  level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4 + RHS5);			    
			    
			    ans(ind_m) =  (RHS1 + RHS2 + RHS3 + RHS4 + RHS5);
			    
			   /* cout<<dc_nb8_21[lev]<<"\n"<<dc_nb8_22[lev]<<"\n"<<dc_nb8_23[lev]<<"\n"<<dc_nb9_21[lev]<<"\n"<<dc_nb9_22[lev]<<"\n"<<dc_nb9_23[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb10_21[lev]<<"\n"<<dc_nb10_22[lev]<<"\n"<<dc_nb10_23[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb11_21[lev]<<"\n"<<dc_nb11_22[lev]<<"\n"<<dc_nb11_23[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb12_21[lev]<<"\n"<<dc_nb12_22[lev]<<"\n"<<dc_nb12_23[lev]<<"\n";
			    	
			    cout<<"--------------------------------------------";*/
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)  //Checked 
			  {    	    
			    RHS1 = dc_i8_33[lev]*X[ind_m-1-2*str_m] + dc_i8_34[lev]*X[ind_m-2*str_m] + dc_i8_35[lev]*X[ind_m+1-2*str_m]; 
			    
			    RHS2 = dc_i9_32[lev]*X[ind_m-2-str_m] + dc_i9_33[lev]*X[ind_m-1-str_m] + dc_i9_34[lev]*X[ind_m-str_m] + dc_i9_35[lev]*X[ind_m+1-str_m] + dc_i9_36[lev]*X[ind_m+2-str_m]; 
			    
			    RHS3 = dc_i10_32[lev]*X[ind_m-2] + dc_i10_33[lev]*X[ind_m-1] + dc_i10_34[lev]*X[ind_m] + dc_i10_35[lev]*X[ind_m+1] + dc_i10_36[lev]*X[ind_m+2]; 			    
			    RHS4 = dc_i11_32[lev]*X[ind_m-2+str_m] + dc_i11_33[lev]*X[ind_m-1+str_m] + dc_i11_34[lev]*X[ind_m+str_m] + dc_i11_35[lev]*X[ind_m+1+str_m] + dc_i11_36[lev]*X[ind_m+2+str_m]; 
			    
			    RHS5 = dc_i12_33[lev]*X[ind_m-1+2*str_m] + dc_i12_34[lev]*X[ind_m+2*str_m] + dc_i12_35[lev]*X[ind_m+1+2*str_m];
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4 + RHS5);
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4 + RHS5);	
			    
			    /*cout<<dc_i8_33[lev]<<"\n"<<dc_i8_34[lev]<<"\n"<<dc_i8_35[lev]<<"\n"<<dc_i9_32[lev]<<"\n"<<dc_i9_33[lev]<<"\n"<<dc_i9_34[lev]<<"\n"<<dc_i9_35[lev]<<"\n"<<dc_i9_36[lev]<<"\n"<<dc_i10_32[lev]<<"\n"<<dc_i10_33[lev]<<"\n"<<dc_i10_34[lev]<<"\n"<<dc_i10_35[lev]<<"\n"<<dc_i10_36[lev]<<"\n"<<dc_i11_32[lev]<<"\n"<<dc_i11_33[lev]<<"\n"<<dc_i11_34[lev]<<"\n"<<dc_i11_35[lev]<<"\n"<<dc_i11_36[lev]<<"\n"<<dc_i12_33[lev]<<"\n"<<dc_i12_34[lev]<<"\n"<<dc_i12_35[lev]<<"\n"; 			    		       
			    
			    cout<<"--------------------------------------------";*/
			  }
			  else if(i==en_inx-sp) //Checked
			  {	    	     	    
			    RHS1 = dc_nb8_21[lev]*X[ind_m+1-2*str_m] + dc_nb8_22[lev]*X[ind_m-2*str_m] + dc_nb8_23[lev]*X[ind_m-1-2*str_m] + dc_nb8_24[lev]*X[ind_m-2-2*str_m]; 
			    
			    RHS2 = dc_nb9_21[lev]*X[ind_m+1-str_m] + dc_nb9_22[lev]*X[ind_m-str_m] + dc_nb9_23[lev]*X[ind_m-1-str_m] + dc_nb9_24[lev]*X[ind_m-2-str_m]; 
			    
			    RHS3 = dc_nb10_21[lev]*X[ind_m+1] + dc_nb10_22[lev]*X[ind_m] + dc_nb10_23[lev]*X[ind_m-1] + dc_nb10_24[lev]*X[ind_m-2]; 
			    
			    RHS4 = dc_nb11_21[lev]*X[ind_m+1+str_m] + dc_nb11_22[lev]*X[ind_m+str_m] + dc_nb11_23[lev]*X[ind_m-1+str_m] + dc_nb11_24[lev]*X[ind_m-2+str_m]; 	   			    
			   
			    if(j+2*sp==en_iny)
			    {			    
			    	RHS5 = dc_nb12_22[lev]*X[ind_m+2*str_m] + dc_nb12_23[lev]*X[ind_m-1+2*str_m] + dc_nb12_24[lev]*X[ind_m-2+2*str_m];	 //lcr
			    }
			    else
			    {
			    	RHS5 = dc_nb12_21[lev]*X[ind_m+1+2*str_m] + dc_nb12_22[lev]*X[ind_m+2*str_m] + dc_nb12_23[lev]*X[ind_m-1+2*str_m] + dc_nb12_24[lev]*X[ind_m-2+2*str_m];	 
			    }		    			     
				
                            //cout<<dc_nb8_23[lev]<<"\n"<<dc_nb8_22[lev]<<"\n"<<dc_nb8_21[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_23[lev]<<"\n"<<dc_nb9_22[lev]<<"\n"<<dc_nb9_21[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_23[lev]<<"\n"<<dc_nb10_22[lev]<<"\n"<<dc_nb10_21[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_23[lev]<<"\n"<<dc_nb11_22[lev]<<"\n"<<dc_nb11_21[lev]<<"\n"<<dc_nb12_23[lev]<<"\n"<<dc_nb12_22[lev]<<"\n"<<dc_nb12_21[lev]<<"\n"; 	
			    	
			    	//cout<<"----------------------------------------------\n";		    	
			    
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4 + RHS5); 			   
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4 + RHS5); 			    
			  }
			  else //Checked 
			  {    	    	    
			    RHS1 = dc_b8_11[lev]*X[ind_m-2*str_m] + dc_b8_12[lev]*X[ind_m-1-2*str_m] + dc_b8_13[lev]*X[ind_m-2-2*str_m]; 
			    
			    RHS2 = dc_b9_11[lev]*X[ind_m-str_m] + dc_b9_12[lev]*X[ind_m-1-str_m] + dc_b9_13[lev]*X[ind_m-2-str_m]; 
			    
			    RHS3 = dc_b10_11[lev]*X[ind_m] + dc_b10_12[lev]*X[ind_m-1] + dc_b10_13[lev]*X[ind_m-2]; 
			    
			    RHS4 = dc_b11_11[lev]*X[ind_m+str_m] + dc_b11_12[lev]*X[ind_m-1+str_m] + dc_b11_13[lev]*X[ind_m-2+str_m]; 	    	    
			    
			    if(j+2*sp==en_iny)
			    {
			    	RHS5 = dc_b12_12[lev]*X[ind_m-1+2*str_m] + dc_b12_13[lev]*X[ind_m-2+2*str_m];  //lcr
			    }
			    else
			    {
			    	RHS5 = dc_b12_11[lev]*X[ind_m+2*str_m] + dc_b12_12[lev]*X[ind_m-1+2*str_m] + dc_b12_13[lev]*X[ind_m-2+2*str_m];
			    }
			    
			    	//cout<<dc_b8_12[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b8_11[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b9_11[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_11[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_11[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_11[lev]<<"\n";
			    	
			    	//cout<<"---------------------------------------\n";
			    
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4 + RHS5);			   
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4 + RHS5); 			  		  
			  }	  	  
			  
			}
	
			/*****************Case-4(j=ny-2)**********************/ 
	
			if(j==en_iny-sp)
			{
			  if(i==st_inx)  //checked 
			  {	     	    
			    RHS1 = dc_b4_11[lev]*X[ind_m+str_m] + dc_b4_12[lev]*X[ind_m+1+str_m] + dc_b4_13[lev]*X[ind_m+2+str_m]; 

			    RHS2 = dc_b5_11[lev]*X[ind_m] + dc_b5_12[lev]*X[ind_m+1] + dc_b5_13[lev]*X[ind_m+2]; 

			    RHS3 = dc_b6_11[lev]*X[ind_m-str_m] + dc_b6_12[lev]*X[ind_m+1-str_m] + dc_b6_13[lev]*X[ind_m+2-str_m]; 
			    
			    RHS4 = dc_b7_11[lev]*X[ind_m-2*str_m] + dc_b7_12[lev]*X[ind_m+1-2*str_m] + dc_b7_13[lev]*X[ind_m+2-2*str_m] ;			    			    
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4);			     	            
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4);
			    			   			    
			    //cout<<dc_pb7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b4_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"; 
			    
			  }
			  else if(i==st_inx+sp) //Checked
			  {	    	    	    
			    RHS1 = dc_nb4_21[lev]*X[ind_m-1+str_m] + dc_nb4_22[lev]*X[ind_m+str_m] + dc_nb4_23[lev]*X[ind_m+1+str_m] + dc_nb4_24[lev]*X[ind_m+2+str_m]; 
			    
			    RHS2 = dc_nb5_21[lev]*X[ind_m-1] + dc_nb5_22[lev]*X[ind_m] + dc_nb5_23[lev]*X[ind_m+1] + dc_nb5_24[lev]*X[ind_m+2]; 
			    
			    RHS3 = dc_nb6_21[lev]*X[ind_m-1-str_m] + dc_nb6_22[lev]*X[ind_m-str_m] + dc_nb6_23[lev]*X[ind_m+1-str_m] + dc_nb6_24[lev]*X[ind_m+2-str_m]; 
			    
			    RHS4 = dc_nb7_21[lev]*X[ind_m-1-2*str_m] + dc_nb7_22[lev]*X[ind_m-2*str_m] + dc_nb7_23[lev]*X[ind_m+1-2*str_m] + dc_nb7_24[lev]*X[ind_m+2-2*str_m]; 	    
			    
			    //ans(ind_m) = level[lev].coeff[ind]*(RHS1 + RHS2 + RHS3 + RHS4);
			    
			    ans(ind_m) = (RHS1 + RHS2 + RHS3 + RHS4); 
			    
			    //cout<<dc_pb7_21[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb4_21[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"; 			    			     	    	    
			  }
			  else if(i>st_inx+sp && i<en_inx-sp) //Checked 
			  {	    			  		
			  		  			  	
			  	if(i+2*sp==en_inx)
			  	{			  
			  		RHS1 = dc_i4_32[lev]*X[ind_m-2+str_m] + dc_i4_33[lev]*X[ind_m-1+str_m] + dc_i4_34[lev]*X[ind_m+str_m] + dc_i4_35[lev]*X[ind_m+1+str_m] ; //lcr 				  	
			  	}
			  	else
			  	{
			  		RHS1 = dc_i4_32[lev]*X[ind_m-2+str_m] + dc_i4_33[lev]*X[ind_m-1+str_m] + dc_i4_34[lev]*X[ind_m+str_m] + dc_i4_35[lev]*X[ind_m+1+str_m] + dc_i4_36[lev]*X[ind_m+2+str_m];			  	
			  	}		  	
			  		
			  	/*cout<<dc_pb7_33[lev]<<"\n"<<dc_pb7_34[lev]<<"\n"<<dc_pb7_35[lev]<<"\n"<<dc_i7_33[lev]<<"\n"<<dc_i7_34[lev]<<"\n"<<dc_i7_35[lev]<<"\n"<<dc_i6_32[lev]<<"\n"<<dc_i6_33[lev]<<"\n"<<dc_i6_34[lev]<<"\n"<<dc_i6_35[lev]<<"\n"<<dc_i6_36[lev]<<"\n"<<dc_i5_32[lev]<<"\n"<<dc_i5_33[lev]<<"\n"<<dc_i5_34[lev]<<"\n"<<dc_i5_35[lev]<<"\n"<<dc_i5_36[lev]<<"\n"<<dc_i4_32[lev]<<"\n"<<dc_i4_33[lev]<<"\n"<<dc_i4_34[lev]<<"\n"<<dc_i4_35[lev]<<"\n"<<dc_i4_36[lev]<<"\n";		  					  	 	
			  		
			  	cout<<"ind_m= "<<ind_m<<"\n"; 
			  		
			  	cout<<"----------------------------------------------------\n"; */
			  				  	    	    			    
			    
			    	RHS2 = dc_i5_32[lev]*X[ind_m-2] + dc_i5_33[lev]*X[ind_m-1] + dc_i5_34[lev]*X[ind_m] + dc_i5_35[lev]*X[ind_m+1] + dc_i5_36[lev]*X[ind_m+2]; 
			    
			    	RHS3 = dc_i6_32[lev]*X[ind_m-2-str_m] + dc_i6_33[lev]*X[ind_m-1-str_m] + dc_i6_34[lev]*X[ind_m-str_m] + dc_i6_35[lev]*X[ind_m+1-str_m] + dc_i6_36[lev]*X[ind_m+2-str_m]; 
			    
			    	RHS4 = dc_i7_32[lev]*X[ind_m-2-2*str_m] + dc_i7_33[lev]*X[ind_m-1-2*str_m] + dc_i7_34[lev]*X[ind_m-2*str_m] + dc_i7_35[lev]*X[ind_m+1-2*str_m] + dc_i7_36[lev]*X[ind_m+2-2*str_m] ;			    			    
			    
			    	//ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 + RHS4 ) ;
			    	
			    	ans(ind_m) = ( RHS1 + RHS2 + RHS3 + RHS4 ) ;			     	    	 
			  }
			  else if(i==en_inx-sp) //checked 
			  {	    	    	    
			    	RHS1 = dc_nb4_22[lev]*X[ind_m+str_m] + dc_nb4_23[lev]*X[ind_m-1+str_m] + dc_nb4_24[lev]*X[ind_m-2+str_m]; //lcr 
			    
			    	RHS2 = dc_nb5_21[lev]*X[ind_m+1] + dc_nb5_22[lev]*X[ind_m] + dc_nb5_23[lev]*X[ind_m-1] + dc_nb5_24[lev]*X[ind_m-2]; 
			    
			    	RHS3 = dc_nb6_21[lev]*X[ind_m+1-str_m] + dc_nb6_22[lev]*X[ind_m-str_m] + dc_nb6_23[lev]*X[ind_m-1-str_m] + dc_nb6_24[lev]*X[ind_m-2-str_m]; 
			    
			    	RHS4 = dc_nb7_21[lev]*X[ind_m+1-2*str_m] + dc_nb7_22[lev]*X[ind_m-2*str_m] + dc_nb7_23[lev]*X[ind_m-1-2*str_m] +  dc_nb7_24[lev]*X[ind_m-2-2*str_m];	       
			    
			    	//ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 + RHS4 ) ;
			    	
			    	ans(ind_m) = ( RHS1 + RHS2 + RHS3 + RHS4 ) ;	
			    				    
			    	//cout<<"ind_m= "<<ind_m<<"\n";
			    
			    	//cout<<dc_pb7_23[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_21[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_22[lev]<<"\n";   
			  }
			  else //Checked 
			  {	    	    	    
			    RHS1 = dc_b4_12[lev]*X[ind_m-1+str_m] + dc_b4_13[lev]*X[ind_m-2+str_m];   	    //lcr
			    
			    RHS2 = dc_b5_11[lev]*X[ind_m] + dc_b5_12[lev]*X[ind_m-1] + dc_b5_13[lev]*X[ind_m-2]; 
			    
			    RHS3 = dc_b6_11[lev]*X[ind_m-str_m] + dc_b6_12[lev]*X[ind_m-1-str_m] + dc_b6_13[lev]*X[ind_m-2-str_m]; 
			    
			    RHS4 = dc_b7_11[lev]*X[ind_m-2*str_m] + dc_b7_12[lev]*X[ind_m-1-2*str_m] +  dc_b7_13[lev]*X[ind_m-2-2*str_m] ;      				   
			    
			    //ans(ind_m) =  level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 + RHS4 ) ;
			    
			    ans(ind_m) =  ( RHS1 + RHS2 + RHS3 + RHS4 ) ;	
			    
			   // cout<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n";			    		    			    	    
			  }  
			}
		
			/****************************************Case-5(j=en_iny)******************************************/ 
		
			if(j==en_iny)
			{
			  if(i==st_inx) 
			  {			  	    	    	    
			    RHS1 = dc_b1_11[lev]*X[ind_m] + dc_b1_12[lev]*X[ind_m+1] + dc_b1_13[lev]*X[ind_m+2]; 
			    
			    RHS2 = dc_b2_11[lev]*X[ind_m-str_m] + dc_b2_12[lev]*X[ind_m+1-str_m] + dc_b2_13[lev]*X[ind_m+2-str_m]; 
			    
			    RHS3 = dc_b3_11[lev]*X[ind_m-2*str_m] + dc_b3_12[lev]*X[ind_m+1-2*str_m] + dc_b3_13[lev]*X[ind_m+2-2*str_m]; 	        			       	    

			    //ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 );
			    
			    ans(ind_m) = ( RHS1 + RHS2 + RHS3 );	 			    
			    
			    //cout<<dc_ny_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b1_11[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"; 			    
			  }
			  else if(i==st_inx+sp) //Checked 
			  {	    	   	    
			    RHS1 = dc_nb1_21[lev]*X[ind_m-1] + dc_nb1_22[lev]*X[ind_m] + dc_nb1_23[lev]*X[ind_m+1] + dc_nb1_24[lev]*X[ind_m+2] ; 
			    
			    RHS2 = dc_nb2_21[lev]*X[ind_m-1-str_m] + dc_nb2_22[lev]*X[ind_m-str_m] + dc_nb2_23[lev]*X[ind_m+1-str_m] + dc_nb2_24[lev]*X[ind_m+2-str_m]; 
			    
			    RHS3 = dc_nb3_21[lev]*X[ind_m-1-2*str_m] + dc_nb3_22[lev]*X[ind_m-2*str_m] + dc_nb3_23[lev]*X[ind_m+1-2*str_m] + dc_nb3_24[lev]*X[ind_m+2-2*str_m];
		    	       	        	    	   	    
		   	    //ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 ) ;
		   	    
		   	    ans(ind_m) = ( RHS1 + RHS2 + RHS3 ) ;    	    
		   	    
		   	    /*cout<<"ind_m = "<<ind_m<<"\n"; 		   	       	    
		   	    
		   	    cout<<dc_ny_21[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb1_21[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_24[lev]<<"\n";*/		   	    		   	        
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked
			  {	    	    	    			   	     	
			     	if(i+2*sp==en_inx)
			     	{			     				     	
			     		RHS1 = dc_i1_32[lev]*X[ind_m-2] + dc_i1_33[lev]*X[ind_m-1] + dc_i1_34[lev]*X[ind_m] + dc_i1_35[lev]*X[ind_m+1]; //lcr 
			     	}
			     	else
			     	{			     	
			     		RHS1 = dc_i1_32[lev]*X[ind_m-2] + dc_i1_33[lev]*X[ind_m-1] + dc_i1_34[lev]*X[ind_m] + dc_i1_35[lev]*X[ind_m+1] + dc_i1_36[lev]*X[ind_m+2];
			     	}			     	
			     	
			     	/*cout<<"ind_m= "<<ind_m<<"\n"; 
			     		
				cout<<dc_ny_32[lev]<<"\n"<<dc_ny_33[lev]<<"\n"<<dc_ny_34[lev]<<"\n"<<dc_ny_35[lev]<<"\n"<<dc_ny_36[lev]<<"\n"<<dc_ny1_33[lev]<<"\n"<<dc_ny1_34[lev]<<"\n"<<dc_ny1_35[lev]<<"\n"<<dc_i3_33[lev]<<"\n"<<dc_i3_34[lev]<<"\n"<<dc_i3_35[lev]<<"\n"<<dc_i2_32[lev]<<"\n"<<dc_i2_33[lev]<<"\n"<<dc_i2_34[lev]<<"\n"<<dc_i2_35[lev]<<"\n"<<dc_i2_36[lev]<<"\n"<<dc_i1_32[lev]<<"\n"<<dc_i1_33[lev]<<"\n"<<dc_i1_34[lev]<<"\n"<<dc_i1_35[lev]<<"\n"<<dc_i1_36[lev]<<"\n";
				
				cout<<"---------------------------------------------------\n";*/ 
			     						  			     
			     
			     RHS2 = dc_i2_32[lev]*X[ind_m-2-str_m] + dc_i2_33[lev]*X[ind_m-1-str_m] + dc_i2_34[lev]*X[ind_m-str_m] + dc_i2_35[lev]*X[ind_m+1-str_m] + dc_i2_36[lev]*X[ind_m+2-str_m]; 
			     
			     RHS3 =  dc_i3_32[lev]*X[ind_m-2-2*str_m] + dc_i3_33[lev]*X[ind_m-1-2*str_m] + dc_i3_34[lev]*X[ind_m-2*str_m] + dc_i3_35[lev]*X[ind_m+1-2*str_m] + dc_i3_36[lev]*X[ind_m+2-2*str_m];   			        
			     
			     //ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 );			    
			     
			     ans(ind_m) = ( RHS1 + RHS2 + RHS3 );
			  }
			  else if(i==en_inx-sp)  //Checked
			  {    	    	   	    			   	        
			    RHS1 = dc_nb1_22[lev]*X[ind_m] + dc_nb1_23[lev]*X[ind_m-1] + dc_nb1_24[lev]*X[ind_m-2]; 	//lcr 		    			    
			    
			    RHS2 = dc_nb2_21[lev]*X[ind_m+1-str_m] + dc_nb2_22[lev]*X[ind_m-str_m] + dc_nb2_23[lev]*X[ind_m-1-str_m] + dc_nb2_24[lev]*X[ind_m-2-str_m]; 
			    
			    RHS3 = dc_nb3_21[lev]*X[ind_m+1-2*str_m] + dc_nb3_22[lev]*X[ind_m-2*str_m] + dc_nb3_23[lev]*X[ind_m-1-2*str_m] + dc_nb3_24[lev]*X[ind_m-2-2*str_m] ;	        
		   	    
		   	    //ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 );
		   	    
		   	    ans(ind_m) = ( RHS1 + RHS2 + RHS3 );
		   	    
		   	    /*cout<<"ind_m= "<<ind_m<<"\n"; 			    	    
		   	    
		   	    cout<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_21[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"; */
			  }
			  else //Checked 
			  {	  	     	     
			     RHS1 = dc_b1_12[lev]*X[ind_m-1] + dc_b1_13[lev]*X[ind_m-2]; //lcr 
			     
			     RHS2 = dc_b2_11[lev]*X[ind_m-str_m] + dc_b2_12[lev]*X[ind_m-1-str_m] + dc_b2_13[lev]*X[ind_m-2-str_m]; 
			     
			     RHS3 = dc_b3_11[lev]*X[ind_m-2*str_m] + dc_b3_12[lev]*X[ind_m-1-2*str_m] + dc_b3_13[lev]*X[ind_m-2-2*str_m];    	    			     			     	    
			     //ans(ind_m) = level[lev].coeff[ind]*( RHS1 + RHS2 + RHS3 ); //This mulvec part belongs to the pinned point. Hence, it is trimmed out and not involved in computation at all. 
			     
			     ans(ind_m) = ( RHS1 + RHS2 + RHS3 ); 			     			     
			  }	  	
			}											
		
		}	
	}

	subans = ans.submat(0,0,tot_p_sol-2,0); 
	return subans; 
}

/*******************************************************************************/
