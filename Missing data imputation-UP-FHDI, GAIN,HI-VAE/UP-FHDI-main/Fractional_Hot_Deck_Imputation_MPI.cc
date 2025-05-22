
void Fractional_Hot_Deck_Imputation(const int i, 
				const int ng, List_FHDI &List_ocsg, const int ncol,
				double** mox, double** y, const int nrow, const int i_M, 
				const int mg, double** z, const int i_mxl, 
				std::vector<int> v_cn_z_i, std::vector<int> v_mxl,
				std::vector<int> v_obsg,
				double* fwij, const int i_size_v_cn_obsg_ncp,
				double* d_obsp, int* i_obsn, 
				const double d_myran, 
				double* w, int* id,
				double** fmat_FHDI)
//Description----------------------
//perform FHDI 
//  Algorithm: impute by using "M" possible donors
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
// updated: March 30, 2017
//
//IN   : int i = the current row 
//IN   : int ng = number of observed donors (=i_size_v_obsg)
//IN   : List_FHDI List_ocsg(nrm); //order records used for variance estimation
//IN   : double** mox(nrow_mox, ncol)  = rows of missing cells 
//IN   : double y(nrow, ncol)= original data matrix with missing cells 
//IN   : int i_M = number of possible donors
//IN   : int mg = Actual row locations that has the same as current missing pattern 
//IN   : double z(nrow, ncol)= categorized matrix of y. 0 for missing cell
//IN   : int i_mxl = number of non-missing column of current row
//IN   : std::vector<int> v_cn_z_i = actual locations of current missing row in cn
//IN   : std::vector<int> v_mxl = a missing row's columns having observed cells  
//IN   : std::vector<int> v_obsg = observed donors ordered by half-ascending and -descending manner
//IN   : double* fwij[i_size_v_cn_obsg_ncp];
//IN   : double* d_obsp = new double[i_size_v_cn_obsg_ncp]; //joint prob of selected donors
//IN   : int*    i_obsn = new int[i_size_v_cn_obsg_ncp]; //counts of the selected donors
//IN   : double d_myran   = random number generated from the uniform distribution
//                          
//IN   : double* w[nrow]
//IN   : int* id[nrow]
//
//OUT  : double** fmat_FHDI(min(ng, i_M) * mg, 7+2*ncol) see details in the above 
//----------------------
{
	//testout
	//RPrint("============== in FHDI  FHDI ==============");
	
	//------------
	//initial data setting
	//------------
	const int i_size_v_obsg = ng; 
	const int i_size_v_cn_z_i = (int)v_cn_z_i.size(); //same as mg and loc2
	
	
	//---------
	//get ocsg, observed donors
	//---------
	double* d_temp_lloc = new double[i_size_v_obsg]; 
	List_ocsg.get_block(i, d_temp_lloc); //get ith stored row from storage 
	//testout
	//RPrint("lloc :"); //RPrint(d_temp_lloc, i_size_v_obsg);

	//--------
	//missing column of current row
	//--------
	std::vector<int> v_rloc; //actual locations of missing
	v_rloc.clear(); //important since this is inside loop
	
	for(int k=0; k<ncol; k++) 
	{ if(fabs(mox[i][k]) < 1e-15) v_rloc.push_back(k+1); }//+1 for actual location
	const int i_size_v_rloc = (int)v_rloc.size(); 
	//testout
	//RPrint("v_rloc :"); //RPrint(v_rloc);
			
	//-------
	//get donors at missing column locations
	//-------
	double** dy = New_dMatrix(i_size_v_obsg, i_size_v_rloc);
	for(int j=0; j<i_size_v_rloc; j++)
	{
		for(int k=0; k<i_size_v_obsg; k++)
		{
			dy[k][j] = y[(int)d_temp_lloc[k]-1][v_rloc[j]-1]; //-1 for actual location 
		}
	}
	//testout
	//RPrint("dy :"); //RPrint(dy, i_size_v_obsg, i_size_v_rloc);
			
	//------------
	//covariance matrix of dy[][]
	//Note: column-wise calculation for the covariance 
	//------------
	double** VM_dy = New_dMatrix(i_size_v_rloc, i_size_v_rloc);
	cov_FHDI(dy, i_size_v_obsg, i_size_v_rloc, VM_dy);
	//testout
	//RPrint("VM_dy :"); //RPrint(VM_dy, i_size_v_rloc, i_size_v_rloc);
			
	//----
	//(a) determine MM depending upon donors and i_M
	//(b) declare memories for FHDI parts
	//----
	int MM = 0; 
	if(i_size_v_obsg <= i_M) MM=i_size_v_obsg; 
	if(i_size_v_obsg >  i_M) MM=i_M; 

	int* i_SN = new int[MM*mg];
	double** d_iy   = New_dMatrix(MM*mg, ncol); //original matrix
	double** d_cmat = New_dMatrix(MM*mg, ncol); //categorized matrix 
	double* wij = new double[ng*mg]; //donors & locations 
	double* d_cs = new double[ng]; 
	double* d_cs_temp = new double[ng]; 
	double* d_Li = new double[ng]; 
	double* d_Ui = new double[ng]; 
	double* d_fefim = new double[i_size_v_rloc];
	int* i_rmg = new int[mg*MM];
	int* i_rM  = new int[mg*MM];
	double* d_SR = new double[mg*MM];
	int* i_obSN = new int[mg*MM];
	//double** fmat_FHDI = New_dMatrix(mg*MM, 7+2*ncol); //7columns and two blocks of ncol 
			
	//------------------------
	//FHDI Case 1: when donors <= M
	//             use all possible donors
	//------------------------
	//int MM=0; 
	if(i_size_v_obsg <= i_M)
	{
		//----
		//index of all donors 
		//----
		for(int j=0; j<mg; j++)
		{
			for(int k=0; k<MM; k++)
				i_SN[j*MM+k] = k+1; //+1 for actual location 
		}
		
		//--------
		//get ready donors from original data
		//--------
		for(int j=0; j<MM*mg; j++)
		{
			for(int k=0; k<ncol; k++)
			{ 
				d_iy[j][k]   = y[v_obsg[i_SN[j]-1]-1][k];  
				d_cmat[j][k] = z[v_obsg[i_SN[j]-1]-1][k];  //-1 for actual loc
			} //-1 for actual location 
		}
		//--------
		//non-missing column consideration of current row
		//--------
		if(i_mxl >= 1)
		{
			for(int k=0; k<i_mxl; k++) //column-wise copy
			{   
				int i_temp_row_iy=0; //sequential index in each row 
				for(int j=0; j<i_size_v_cn_z_i; j++) //Note: i_size_v_cn_z_i = length(loc2)
				{
					//Note v_cn_z_i = loc2 //
					double d_temp_iy = y[v_cn_z_i[j] - 1][v_mxl[k]-1];  //-1 for actual 
		
					for( int j2=0; j2< MM; j2++) //repeat 
						d_iy[i_temp_row_iy++][v_mxl[k]-1] = d_temp_iy;   						
				}

			}  					
		}
		
		//------------
		//weights 
		//------------
		for(int j=0; j<mg; j++)
		{
			for(int k=0; k<ng; k++) wij[ng*j + k] = fwij[k]; 
		}
		
		//testout
		//RPrint("======== in FHDI Case 1 obsg <= M ===========");
		//RPrint("SN :"); //RPrint(i_SN, MM*mg);
		//RPrint("iy :"); //RPrint(d_iy,  MM*mg, ncol);
		//RPrint("d_cmat :"); //RPrint(d_cmat,  MM*mg, ncol);
		//RPrint("wij :"); //RPrint(wij, ng*mg);
	}
	
	//------------------------
	//FHDI Case 2: when donors > M
	//   select donors with probability proportional to size sampling 
	//------------------------
	if(i_size_v_obsg > i_M)
	{
		//-------------
		//cumulative sum of joint probability with pps 
		//-------------
		for(int j=0; j<ng; j++) 
		{    
			d_cs_temp[j] = d_obsp[j]; //default  
			if(i_obsn[j] !=0) d_cs_temp[j] = d_obsp[j]*MM/i_obsn[j];
		} 
		cumsum_FHDI(d_cs_temp, ng, d_cs); 
		
		//----------
		//get ready Li and Ui
		//----------
		d_Li[0] = 0.0; 
		for(int j=1; j<ng; j++)
		{
			d_Li[j] = d_cs[j-1]; //exclude the last entity 
		}
		Copy_dVector(d_cs, ng, d_Ui); 

		//testout
		//RPrint("======== in FHDI Case 2 obsg > M ===========");
		//RPrint("cs :"); //RPrint(d_cs, ng);
		//RPrint("Li :"); //RPrint(d_Li, ng);
		//RPrint("Ui :"); //RPrint(d_Ui, ng);
		
		//-----------------
		//sum(dy*fwij)
		//-----------------
		double d_sum_fefim=0.0;
		for(int j=0; j<i_size_v_rloc; j++)
		{
			d_sum_fefim = 0.0; 
			//---row-wise sum---//
			for(int k=0; k<i_size_v_obsg; k++)
			{
				d_sum_fefim += dy[k][j]*fwij[k]; 
			}
			d_fefim[j] = d_sum_fefim; 
		}
		//testout
		//RPrint("fefim :"); //RPrint(d_fefim, i_size_v_rloc);
		
		//-----------------
		//random location using uniform distribution   
		//using Numerical Recipes of Press et al 2007. 
		//-----------------
		//double d_Rg_myran = myran.doub();
		//cout<<"d_Rg using myran():"<<d_Rg_myran<<endl;
		
		//simple version using standard rand() for CRAN compatibility, March 30, 2017
		double d_Rg = d_myran;
		//cout<<"d_Rg using rand():"<<d_Rg<<endl;
		
		//double d_Rg = 0.0633672; //for debugging !!!
		
		//below codes is not recommended by Press et al. 2007. 
		//below is only available for c++ compiler after 2011 
		//std::default_random_engine generator; 
		//std::uniform_real_distribution<double> distribution(0.0, 1.0); 
		//d_Rg = distribution(generator);
		
		for(int j=0; j<mg; j++) 
		{	for(int k=0; k<MM; k++) i_rmg[j*MM+k]=j+1;  }
	
		for(int j=0; j<mg; j++) 
		{	for(int k=0; k<MM; k++)  i_rM[j*MM+k]=k+1;  }
		
		for(int j=0; j<mg*MM; j++)
		{
			d_SR[j] = (d_Rg+(i_rmg[j]-1))/mg + (i_rM[j]-1); 
		}

		//----------
		//set of SR < Ui
		//----------
		int i_SR_Ui = 0;  
		for(int k=0; k<mg*MM; k++)
		{
			i_SR_Ui = 0; //minimum location where SR <= Ui
			for(int j=0; j<ng; j++) 
			{  if(d_SR[k] <= d_Ui[j]) {i_SR_Ui = j+1; break;} }
			
			i_SN[k] = i_SR_Ui; //Actual location stored 
		}
		//testout
		//RPrint("Rg :"); //RPrint(d_Rg);
		//RPrint("SR :"); //RPrint(d_SR, mg*MM);
		//RPrint("SN :"); //RPrint(i_SN, mg*MM);
		
		//-------------
		// select out M donors from half-asc and -desc observation
		//-------------
		for(int j=0; j<mg*MM; j++)
		{ i_obSN[j] = v_obsg[i_SN[j]-1]; } //-1 for actual location 

		//----------------------------
		//impute missing cells from original data matrix 
		//----------------------------
		for(int j=0; j<mg*MM; j++)
		{   
			for(int k=0; k<ncol; k++)
			{
				d_iy[j][k]   = y[i_obSN[j]-1][k];  //-1 for actual size  
				d_cmat[j][k] = z[i_obSN[j]-1][k];  //-1 for actual loc		
			}
		}

		//impute missing cells using donors  
		if(i_mxl >= 1)
		{
			for(int k=0; k<i_mxl; k++) //column-wise copy
			{   
				int i_temp_row_iy=0; //sequential index in each row 
				for(int j=0; j<mg; j++) //cf. i_size_v_cn_z_i = length(loc2)
				{
					//Note v_cn_z_i = loc2 //
					double d_temp_iy = y[v_cn_z_i[j] - 1][v_mxl[k]-1];  //-1 for actual 
		
					for( int j2=0; j2< MM; j2++) //repeat 
						d_iy[i_temp_row_iy++][v_mxl[k]-1] = d_temp_iy;   						
				}
			}  
		}
		
		//------------
		//diy : not used. So not implemented as of Nov 15, 2016
		//------------
		
		//------------
		//wij
		//------------
		for(int j=0; j<mg*MM; j++) wij[j] = 1.0/MM; 
	}			
	
	//------------------------------
	//------------------------------
	//make return matrix fmat_FHDI[][]
	// Note: row number differs from fmat[][]
	//       must be appended to previous fmat[][] 
	//------------------------------
	//------------------------------
	//column 1: id[]// note: v_cn_z_i means loc2 
	for(int j=0; j<i_size_v_cn_z_i; j++)
	{ 
		int i_temp_fmat_col1 = id[v_cn_z_i[j]-1]; //-1 for actual location 
		
		for(int k=0; k<MM; k++) //repeat each id MM times
		{
			fmat_FHDI[MM*j + k][0] = i_temp_fmat_col1; 
		}
	}		
	
	//column 2: 1:MM repeated mg times // note: v_cn_z_i = loc2 
	for(int j=0; j<mg; j++)
	{ 				
		for(int k=0; k<MM; k++) 
		{
			fmat_FHDI[MM*j + k][1] = k+1; //store actual number  
		}
	}
				
	//column 3: w[]// note: v_cn_z_i means loc2 
	for(int j=0; j<i_size_v_cn_z_i; j++)
	{ 
		double d_temp_fmat_col3 = w[v_cn_z_i[j]-1]; //-1 for actual location 
		
		for(int k=0; k<MM; k++) 
		{
			fmat_FHDI[MM*j + k][2] = d_temp_fmat_col3; 
		}
	}
	
	//column 4: wij[]// note: v_cn_z_i means loc2 
	for(int j=0; j<i_size_v_cn_z_i; j++)
	{ 
		for(int k=0; k<MM; k++) 
		{
			fmat_FHDI[MM*j + k][3] = wij[MM*j + k]; 
		}
	}				
	
	//column set 5: [4, (4+ncol-1)]: d_iy[][ncol]
	int i_begin = 4; //starting point of current column set 
	for(int j=0; j<i_size_v_cn_z_i; j++)
	{ 
		for(int k=0; k<MM; k++) 
		{
			for(int i_col_iy=0; i_col_iy<ncol; i_col_iy++)
			{  fmat_FHDI[MM*j + k][i_begin+i_col_iy] 
					   = d_iy[MM*j + k][i_col_iy];      }
			 
		}
	}				
	
	//column set 6: [4+ncol, (4+2*ncol-1)]: d_cmat[][ncol]
	i_begin = 4+ncol; //starting point of current column set 
	for(int j=0; j<i_size_v_cn_z_i; j++)
	{ 
		for(int k=0; k<MM; k++) 
		{
			for(int i_col_iy=0; i_col_iy<ncol; i_col_iy++)
			{  fmat_FHDI[MM*j + k][i_begin+i_col_iy] 
					   = d_cmat[MM*j + k][i_col_iy];      }
			 
		}
	}
	
	//column set 7: at 4+2*ncol. obsg[]
	i_begin = 4+2*ncol; //starting point of current column set 
	for(int j=0; j<i_size_v_cn_z_i; j++) //Note: = mg 
	{ 
		for(int k=0; k<MM; k++) 
		{
			fmat_FHDI[MM*j + k][i_begin] 
					   = v_obsg[i_SN[MM*j+k]-1];      
		}
	}			
	
	//column set 8: at 4+2*ncol+1. 1:ng repeated by mg times 
	i_begin = 4+2*ncol+1; //starting point of current column set 
	for(int j=0; j<mg; j++) // 
	{ 
		for(int k=0; k<MM; k++) 
		{
			fmat_FHDI[MM*j + k][i_begin] = i_SN[j*MM + k];      
		}
	}				
	
	//column set 9: at 4+2*ncol+2. fwij repeated by mg times 
	i_begin = 4+2*ncol+2; //starting point of current column set 
	for(int j=0; j<mg; j++) // 
	{ 
		for(int k=0; k<MM; k++) 
		{
			fmat_FHDI[MM*j + k][i_begin] = fwij[i_SN[j*MM+k]-1];      
		}
	}			
		
	//testout
	//RPrint("cmat : "); //RPrint(d_cmat, mg*MM, ncol);
	//RPrint("iy : "); //RPrint(d_iy, mg*MM, ncol);
	//RPrint("wij : "); //RPrint(wij, mg*MM);
	//RPrint("fmat_FHDI : "); //RPrint(fmat_FHDI, mg*MM, 7+2*ncol);
	
	
	//-------------------
	//local deallocation
	//-------------------
	delete[] d_temp_lloc; 
	Del_dMatrix(dy, i_size_v_obsg, i_size_v_rloc);	
	Del_dMatrix(VM_dy, i_size_v_rloc, i_size_v_rloc);	

	//----------
	//local deallocation for FHDI parts
	//----------
	delete[] d_cs;
	delete[] d_cs_temp;	
	delete[] d_Li;
	delete[] d_Ui;	
	delete[] d_fefim; 
	delete[] i_rmg; 
	delete[] i_rM; 
	delete[] d_SR; 
	delete[] i_SN; 
	delete[] i_obSN; 
	Del_dMatrix(d_iy,   mg*MM, ncol);
	Del_dMatrix(d_cmat, mg*MM, ncol);
	delete[] wij; 
	//Del_dMatrix(fmat_FHDI, mg*MM, 7+2*ncol);			
		

	return; 
}