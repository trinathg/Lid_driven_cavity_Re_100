#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

///////////////////Function overloaded of Writing to a file////////////////////
//Writing data at frequency specified by file_freq in the "init" file

void write_to_file(vertex * node, fval * fvar, int t)
{
	ofstream out_put, vel_put0, vel_put1; //0->x-velocity along y-direction at geometric center and 1->y-velocity along x-direction at geom center  
	ofstream vort1, vort2, flux; 
	
	ofstream outpmtv, outvmtv, outvortmtv; 
	
	int p,p1,p2;
	double udy=0.0, vdx=0.0; 

        tdma1y(fvar,0);
	tdma1x(fvar,1);  

	string file("data_2D");
	
	string filep("pressure_mtv"); 
	string filev("vector_mtv"); 
	string filevort("vorticity_mtv"); 
	
	string file0("x-vel");
	string file1("y-vel");
        string file2("vort-mov");  
  	string file3("vort-cen");
	string file4("mass_flux_gc");

	stringstream tag; 

	tag<<t;

	file = file + "_" + tag.str() + ".vtk";   
	
	filep = filep + "_" + tag.str() + ".dat";   
	filev = filev + "_" + tag.str() + ".dat";   
	filevort = filevort + "_" + tag.str() + ".dat";   
	
	file0 = file0 + "_" + tag.str() + ".dat"; 
	file1 = file1 + "_" + tag.str() + ".dat";
	file2 = file2 + "_" + tag.str() + ".dat";
	file3 = file3 + "_" + tag.str() + ".dat";
	file4 = file4 + "_" + tag.str() + ".dat";

	out_put.open(file.c_str());
	outpmtv.open(filep.c_str());
	outvmtv.open(filev.c_str());
	outvortmtv.open(filevort.c_str());

	vel_put0.open(file0.c_str());
	vel_put1.open(file1.c_str());
        vort1.open(file2.c_str()); 
        vort2.open(file3.c_str());
	flux.open(file4.c_str());

	out_put<<"# vtk DataFile Version 2.0"<<"\n"; 
	out_put<<"2DGrid"<<"\n";
	out_put<<"ASCII"<<"\n";
	out_put<<"DATASET"<<" "<<"STRUCTURED_GRID"<<"\n";
	out_put<<"DIMENSIONS "<<nx+1<<" "<<ny+1<<" 1"<<"\n";
	out_put<<"POINTS"<<" "<<tot_p<<" float"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<node[p].x[0]<<" "<<node[p].x[1]<<" 0.0"<<"\n"; 
		}
	}

	out_put<<"POINT_DATA "<<tot_p<<"\n";
	out_put<<"SCALARS"<<" Pressure"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" Pressure_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[2]<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" divergence"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" div_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].div<<"\n"; 
		}
	}
	out_put<<"SCALARS"<<" residual"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" res_lp"<<"\n";

	for(int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].res<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" F"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" F_lp"<<"\n";

	for(int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].F<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" ps1"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" ps1_lp"<<"\n";

	for(int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[3]<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" ps2"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" ps2_lp"<<"\n";

	for(int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[4]<<"\n"; 
		}
	}

	out_put<<"SCALARS"<<" Vorticity"<<" double"<<"\n";
        out_put<<"LOOKUP_TABLE"<<" Vorticity_lp"<<"\n";

        for(int j=0;j<=ny;j++)
        {
                for (int i=0;i<=nx;i++)
                {
                        p = i + j*str_x;
                        out_put<<fvar[p].uy[0]-fvar[p].ux[1]<<"\n";
                }
        }


	out_put<<"VECTORS"<<" Velocity"<<" double"<<"\n"; 

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[0]<<" "<<fvar[p].u[1]<<" "<<"0.0"<<"\n";
		}
	}

        
	//Writing the x-velocity in the y-direction at the geometric center of the domain 
        
	for(int j=0;j<=ny;j++)
	{
		vel_put0<<node[nx/2 + j*str_x].x[1]<<"	"<<fvar[nx/2 + j*str_x].u[0]<<"\n"; 
	}

	//Writing the y-velocity in the x-direction at the geometric center of the domain 

	for(int i=0;i<=nx;i++)
	{
		vel_put1<<node[i+(ny/2)*str_x].x[0]<<"	"<<fvar[i+(ny/2)*str_x].u[1]<<"\n"; 

	}

	//Writing vorticity along the moving plate

	for(int i=1;i<nx;i++)
	{
		p = i + ny*str_x;                
		vort1<<node[p].x[0]<<"	"<<fvar[p].uy[0] - fvar[p].ux[1]<<"\n";
	}

        //Writing vorticity at the geometric center along a vertical line

	for(int j=0;j<=ny;j++)
        {
		p = nx/2 + j*str_x; 
                vort2<<node[p].x[1]<<"  "<<fvar[p].uy[0] - fvar[p].ux[1]<<"\n";
        }

	//Writing mass fluxes 

	for(int j=1;j<=ny;j++)
	{
		p1 = (nx/2) + (j-1)*str_x; 
		p2 = (nx/2) + j*str_x; 

		udy = udy + (fvar[p1].u[0] + fvar[p2].u[0])*0.5*dy; 
	}

	for(int i=1;i<=nx;i++)
	{
		p = i + (ny/2)*str_x; 		
		
		vdx = vdx + (fvar[p-1].u[1] + fvar[p].u[1])*0.5*dx; 
	}	

	flux<<"udy= "<<udy<<"\n"; 
	flux<<"vdx= "<<vdx<<"\n"; 
	
	/****************************************************************/
	/************Writing data in the PLOTMTV format******************/	

	outpmtv<<"$DATA=CONTOUR\n";
	outpmtv<<"%xmin= 0"<<"\n";
	outpmtv<<"%ymin= 0"<<"\n";
	outpmtv<<"%xmax= "<<l_x<<"\n";
	outpmtv<<"%ymax= "<<l_y<<"\n";
	outpmtv<<"%xlabel=\"x\""<<"\n";
	outpmtv<<"%ylabel=\"y\""<<"\n";
	outpmtv<<"%contstyle=2"<<"\n";
	outpmtv<<"%nx="<<nx+1<<"\n";
	outpmtv<<"%ny="<<ny+1<<"\n";
	outpmtv<<"%cmin=-17000"<<"\n";
	outpmtv<<"%cmax=17000"<<"\n";
	outpmtv<<"%nsteps=50"<<"\n";
	outpmtv<<"%INTERPOLATE=2"<<"\n";
	
	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			outpmtv<<fvar[p].u[2]<<"\n"; 
		}
	}
	
	outvortmtv<<"$DATA=CONTOUR\n";
	outvortmtv<<"%xmin= 0"<<"\n";
	outvortmtv<<"%ymin= 0"<<"\n";
	outvortmtv<<"%xmax= "<<l_x<<"\n";
	outvortmtv<<"%ymax= "<<l_y<<"\n";
	outvortmtv<<"%xlabel=\"x\""<<"\n";
	outvortmtv<<"%ylabel=\"y\""<<"\n";
	outvortmtv<<"%contstyle=2"<<"\n";
	outvortmtv<<"%nx="<<nx+1<<"\n";
	outvortmtv<<"%ny="<<ny+1<<"\n";
	outvortmtv<<"%cmin=-17000"<<"\n";
	outvortmtv<<"%cmax=17000"<<"\n";
	outvortmtv<<"%nsteps=50"<<"\n";
	outvortmtv<<"%INTERPOLATE=2"<<"\n";
	
	for(int j=0;j<=ny;j++)
        {
                for (int i=0;i<=nx;i++)
                {
                        p = i + j*str_x;
                        outvortmtv<<fvar[p].uy[0]-fvar[p].ux[1]<<"\n";
                }
        }
	
	
	outvmtv<<"$DATA=VECTOR\n";	
	outvmtv<<"%xmin=0"<<"	"<<"xmax = "<<l_x<<"\n";
	outvmtv<<"%ymin=0"<<"	"<<"ymax = "<<l_y<<"\n";
	outvmtv<<"%zmin=0"<<"	"<<"zmax = 0"<<"\n";	
	outvmtv<<"%vscale = 0.015"<<"\n";	
	
	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			outvmtv<<node[p].x[0]<<" "<<node[p].x[1]<<" 0.0"<<" "<<fvar[p].u[0]<<" "<<fvar[p].u[1]<<" "<<"0.0"<<"\n";
		}
	}
		       
        vel_put0.close();
	vel_put1.close();
	vort1.close();
	vort2.close(); 
	flux.close();
		
	out_put.close(); 	 
	outpmtv.close(); 
	outvmtv.close();
}

///////////////////Writing to a file////////////////////

void write_to_file(vertex * node, fval * fvar)
{

	ofstream out_put;  
	int p; 

	out_put.open("model_field.vtk");

	out_put<<"# vtk DataFile Version 2.0"<<"\n"; 
	out_put<<"2DGrid"<<"\n";
	out_put<<"ASCII"<<"\n";
	out_put<<"DATASET"<<" "<<"STRUCTURED_GRID"<<"\n";
	out_put<<"DIMENSIONS "<<nx+1<<" "<<ny+1<<" 1"<<"\n";
	out_put<<"POINTS"<<" "<<tot_p<<" float"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<node[p].x[0]<<" "<<node[p].x[1]<<" 0.0"<<"\n"; 
		}
	}

	out_put<<"POINT_DATA "<<tot_p<<"\n";
	out_put<<"SCALARS"<<" Pressure"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" Pressure_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[2]<<"\n"; 			
						
		}
	}

	out_put<<"SCALARS"<<" rhs"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" rhs"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].rhs<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" F"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" F"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].F<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" trunc_error"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" trunc_error"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[2]<<"\n"; 
		}
	}	
	
	
	out_put<<"VECTORS"<<" Velocity"<<" double"<<"\n"; 

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[0]<<" "<<fvar[p].u[1]<<" "<<"0.0"<<"\n";
		}
	}

	out_put.close(); 

}
