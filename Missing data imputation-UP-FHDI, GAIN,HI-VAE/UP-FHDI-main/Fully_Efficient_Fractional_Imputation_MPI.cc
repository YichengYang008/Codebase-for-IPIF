#include <vector>

void Fully_Efficient_Fractional_Imputation(const int ng, const int mg, 
		std::vector<int> v_obsg, std::vector<int> v_mxl, 
		double** y, double** z, const int nrow, const int ncol,
		std::vector<int> v_cn_z_i, double* fwij, const int i_size_v_cn_obsg_ncp,
		double* w, int* id,
		double** fmat)
//Description----------------------
//perform FEFI 
//  Algorithm: impute by using all possible donors
//  final outcome is "fmat" in which each column means that
//  col1: id
//  col2: fid, i.e., id of imputed value
//  col3: sampling weight
//  col4: fractional weights 
//  col5: imputed original data (matrix with column of ncol) 
//  col6: imputed category data (matrix with column of ncol)
//  col7: 1:ng 
//  col8: = col2  (for consistency with FHDI results)
//  col9: = col3  (for consistency with FHDI results)
//
//
// original R code: Dr. Im, J. and Dr. Kim, J. 
// c++ code: 		Dr. Cho, I. 
// All rights reserved
// 
// updated: Nov 14, 2016
//
//IN   : int ng = number of observed donors (=i_size_v_obsg)
//IN   : int mg = Actual row locations that has the same as current missing pattern 
//IN   : std::vector<int> v_obsg = observed donors ordered by half-ascending and -descending manner
//IN   : std::vector<int> v_mxl = a missing row's columns having observed cells  
//IN   : double y(nrow, ncol)= original data matrix with missing cells 
//IN   : double z(nrow, ncol)= categorized matrix of y. 0 for missing cell
//IN   : std::vector<int> v_cn_z_i = actual locations of current missing row in cn
//IN   : double* fwij[i_size_v_cn_obsg_ncp];
//IN   : double* w[nrow]
//IN   : int* id[nrow]
//
//OUT  : double** fmat(ng*mg, 7+2*ncol) see details in the above 
//----------------------
{
	double* wij = new double[ng*mg]; //donors & locations 
	//double** fmat = New_dMatrix(ng*mg, 7+2*ncol); //7columns and two blocks of ncol 
	
	//-----------
	//---local sizes
	//-----------
	const int i_size_v_obsg = (int)v_obsg.size(); 
	const int i_size_v_cn_z_i = (int)v_cn_z_i.size(); 
	
	std::vector<int> v_obsg_times_mg; 
	for(int j=0; j<mg; j++)
		for(int k=0; k<i_size_v_obsg; k++)
		{   v_obsg_times_mg.push_back(v_obsg[k]);   }

	const int i_size_v_obsg_times_mg = (int)v_obsg_times_mg.size(); 

	//----------------------------
	//impute missing cells from original data matrix 
	//----------------------------
	double** d_iy = New_dMatrix(i_size_v_obsg_times_mg, ncol); //imputed original matrix
	double** d_cmat = New_dMatrix(i_size_v_obsg_times_mg, ncol); //imputed category matrix

	const int i_size_v_mxl = (int)v_mxl.size(); 
	int i_mxl = 0; //default of non-missing cell of this row  
	if(i_size_v_mxl >= 1) i_mxl = i_size_v_mxl;  
	
	//default imputed cell is the observed original cells 
	for(int j=0; j<i_size_v_obsg_times_mg; j++)
	{   
		for(int k=0; k<ncol; k++)
		{
			d_iy[j][k] = y[v_obsg_times_mg[j] - 1][k];  //-1 for actual size  
			d_cmat[j][k] = z[v_obsg_times_mg[j] - 1][k];  //-1 for actual loc		
		}
	}  
	
	//impute missing cells using donors  
	if(i_mxl >= 1)
	{
		for(int k=0; k<i_mxl; k++) //column-wise copy
		{   
			int i_temp_row_iy=0; //sequential index in each row 
			for(int j=0; j<i_size_v_cn_z_i; j++) //Note: i_size_v_cn_z_i = length(loc2)
			{
				//Note v_cn_z_i = loc2 //
				double d_temp_iy = y[v_cn_z_i[j] - 1][v_mxl[k]-1];  //-1 for actual 
				
				for( int j2=0; j2< ng; j2++) //repeat 
					d_iy[i_temp_row_iy++][v_mxl[k]-1] = d_temp_iy;   						
			}

		}  
	}
	
	//testout
	////RPrint("in ====== M=FEFI =====");
	////RPrint("mxl: "); //RPrint(v_mxl);
	////RPrint("iy : "); //RPrint(d_iy, i_size_v_obsg_times_mg, ncol);
	////RPrint("cmat : "); //RPrint(d_cmat, i_size_v_obsg_times_mg, ncol);
	
	//------------------------------
	//------------------------------
	//make return matrix fmat[][]
	//------------------------------
	//------------------------------
	for(int j=0; j<mg; j++)
	{
		for(int k=0; k<ng; k++) wij[ng*j + k] = fwij[k]; 
	}
	
	//column 1: id[]// note: v_cn_z_i means loc2 
	for(int j=0; j<i_size_v_cn_z_i; j++)
	{ 
		int i_temp_fmat_col1 = id[v_cn_z_i[j]-1]; //-1 for actual location 
		
		for(int k=0; k<ng; k++) //repeat each id ng times
		{
			fmat[ng*j + k][0] = i_temp_fmat_col1; 
		}
	}		
	
	//column 2: 1:ng repeated mg times // note: v_cn_z_i = loc2 
	for(int j=0; j<mg; j++)
	{ 				
		for(int k=0; k<ng; k++) 
		{
			fmat[ng*j + k][1] = k+1; //store actual number  
		}
	}
				
	//column 3: w[]// note: v_cn_z_i means loc2 
	for(int j=0; j<i_size_v_cn_z_i; j++)
	{ 
		double d_temp_fmat_col3 = w[v_cn_z_i[j]-1]; //-1 for actual location 
		
		for(int k=0; k<ng; k++) 
		{
			fmat[ng*j + k][2] = d_temp_fmat_col3; 
		}
	}
	
	//column 4: wij[]// note: v_cn_z_i means loc2 
	for(int j=0; j<i_size_v_cn_z_i; j++)
	{ 
		for(int k=0; k<ng; k++) 
		{
			fmat[ng*j + k][3] = wij[ng*j + k]; 
		}
	}				
	
	//column set 5: [4, (4+ncol-1)]: d_iy[][ncol]
	int i_begin = 4; //starting point of current column set 
	for(int j=0; j<i_size_v_cn_z_i; j++)
	{ 
		for(int k=0; k<ng; k++) 
		{
			for(int i_col_iy=0; i_col_iy<ncol; i_col_iy++)
			{  fmat[ng*j + k][i_begin+i_col_iy] 
					   = d_iy[ng*j + k][i_col_iy];      }
			 
		}
	}				
	
	//column set 6: [4+ncol, (4+2*ncol-1)]: d_cmat[][ncol]
	i_begin = 4+ncol; //starting point of current column set 
	for(int j=0; j<i_size_v_cn_z_i; j++)
	{ 
		for(int k=0; k<ng; k++) 
		{
			for(int i_col_iy=0; i_col_iy<ncol; i_col_iy++)
			{  fmat[ng*j + k][i_begin+i_col_iy] 
					   = d_cmat[ng*j + k][i_col_iy];      }
			 
		}
	}
	
	//column set 7: at 4+2*ncol. obsg[]
	i_begin = 4+2*ncol; //starting point of current column set 
	for(int j=0; j<i_size_v_cn_z_i; j++) //Note: = mg 
	{ 
		for(int k=0; k<ng; k++) 
		{
			fmat[ng*j + k][i_begin] = v_obsg[k];      
		}
	}			
	
	//column set 8: at 4+2*ncol+1. 1:ng repeated by mg times 
	i_begin = 4+2*ncol+1; //starting point of current column set 
	for(int j=0; j<mg; j++) // 
	{ 
		for(int k=0; k<ng; k++) 
		{
			fmat[ng*j + k][i_begin] = k+1;      
		}
	}				
	
	//column set 9: at 4+2*ncol+2. fwij repeated by mg times 
	i_begin = 4+2*ncol+2; //starting point of current column set 
	for(int j=0; j<mg; j++) // 
	{ 
		for(int k=0; k<ng; k++) 
		{
			fmat[ng*j + k][i_begin] = fwij[k];      
		}
	}			

	//testout
	/*
	//RPrint("in ==== M=FEFI =after making fmat[][]====");
	//RPrint("mxl: "); //RPrint(v_mxl);
	//RPrint("iy : "); //RPrint(d_iy, i_size_v_obsg_times_mg, ncol);
	//RPrint("cmat : "); //RPrint(d_cmat, i_size_v_obsg_times_mg, ncol);
	//RPrint("wij : "); //RPrint(wij, ng*mg);
	//RPrint("fmat : "); //RPrint(fmat, ng*mg, 7+2*ncol);
	*/
	
	//------------
	//Deallocation 
	//------------
	Del_dMatrix(d_iy, i_size_v_obsg_times_mg, ncol);
	Del_dMatrix(d_cmat, i_size_v_obsg_times_mg, ncol);
	
	return; 
}