#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

void write_restart(fval * fvar,int t)
{
  int ind; 
  ofstream output1,output2; 
        
  output1.open("restart_field.dat");  //Writing the restart field    
  output2.open("restart_count.dat");  //counting the restart
  
  //Writing the restart count file 
  
  output2<<t<<"\n"; 
  output2.close();   
  
  //Writing down U,V,P
  
  //Uface values 
  
  for(int j=0;j<=ny;j++)
  {
    for(int i=0;i<=nx;i++)
    {
      ind = i + j*str_x; 

      output1<<i<<"	"<<j<<"		"<<fvar[ind].u[0]<<"		"<<fvar[ind].u[1]<<"		"<<fvar[ind].u[2]<<"\n";  
    }  
  } 
  
  output1.close();
  output2.close();  
}

/****************************************************************************/

void read_restart(fval * fvar, int & restart_count)
{
  int ind,p,q;      
  ifstream input1, input2;   
  
  input1.open("restart_field.dat");   
  input2.open("restart_count.dat");
  
  input2 >> restart_count; 
  input2.close(); 
  
  //Reading in Uface values  
  for(int j=0;j<=ny;j++)
  {
    for(int i=0;i<=nx;i++)
    {      
      input1 >> p >> q >> fvar[i+j*str_x].u[0]>>fvar[i+j*str_x].u[1]>>fvar[i+j*str_x].u[2]; 
    }  
  }
    
  input1.close();
}

/****************************************************************************/
void mg_restart(fval * fvar, mg_grid * level)
{
	int ind; 
	
	  for(int j=0;j<=ny;j++)
  	  {
	    for(int i=0;i<=nx;i++)
    	    {      
    	    	ind = i + j*str_x; 
      		level[0].phi_s[ind] = fvar[ind].u[2];      		
    	    }  
  	  }	
}
