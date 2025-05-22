#include <vector>
#include "order_FHDI_binary.cc"
void AGMAT_Extension_cpp(double** mox, const int nrow_mox, 
						 double** uox, const int nrow_uox, 
						 const int ncol, int* id, 
						 std::vector<std::string> v_table_tmvec_row1,
						 std::vector<int> v_table_tmvec_row2,
                         std::string cn[], const int nrow, 
	                     std::vector<std::string> &rst_final)
//Description=========================================
// Augment missing cells in mox using the observed values of uox
//
// Algorithm:  All possible donors will be used to fill in the missing cell 
//             but, if there is no matched donors in uox, this algorithm may fail
//             as of Oct 2016
// for each missing pattern, find all the possible donors
// e.g., 
// (1) a missing row   = 000
// 	   agmat           = all observed rows
// (2) a missing row   = a01
//     agmat           = ac1, af1, a11, ..., az1. 
//
//
// original R code: Dr. Im, J. and Dr. Kim, J. 
// c++ code: 		Dr. Cho, I. 
// All rights reserved
// 
// updated: Oct 26, 2016
//----------------------------------------------------
//IN    : double mox(nrow_mox, ncol)= sorted unique patterns of missing  cells. up to i_count_mox rows are meaningful                           
//IN    : double uox(nrow_uox, ncol)= sorted unique patterns of observed cells. up to i_count_uox rows are meaningful 
//IN    : int id(nrow) = index of row. Default is ACTUAL row number
//IN	: vector<string> v_table_tmvec_row1  = name table of condensed missing patterns
//IN	: vector<int> v_table_tmvec_row2  = counts table of condensed missing patterns
//IN  	: string cn(nrow)		= vector of string to represent each row of z     
//OUT   : rbind_FHDI rst_final(??, ncol) = updated observed rows to be used later 
//                                          number of rows will be determined by this code   
//====================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	const int nr1 = nrow_mox;
	const int nr2 = nrow_uox;
	
	//-------------------
	//initialize rst, the matrix for storage for augmented observations
	//-------------------
	//rbind_FHDI rst(ncol+1); //1+ncol = (row id)+(observed rows used for imputation later) 
	std::vector<std::string> rst;
	std::vector<int> rst_id;
	std::string uox1; //one row of translated string 


	//--------------------
	//Main Loop for all missing rows
	//--------------------
	int* i_temp_x = new int[ncol];
	int i_sum_x = 0;
	std::vector<int> v_cn_z_i; 
	//int* zid = NULL;
	int i_size_zid=0; 
	int i_loc=0;	
	int* i_srst = new int[nr2];
	std::vector<int> loc_srst_nl; 
	
	//-------------
	//LOOP for all missing rows
	//-------------
	for(int i=0; i<nr1; i++)
	{
		//get current row of missing cell 
		for(int j=0; j<ncol; j++) i_temp_x[j] =  (int)mox[i][j]; 
		i_sum_x = sum_FHDI(i_temp_x, ncol);

		std::string s_temp = v_table_tmvec_row1[i]; //name of ith row
		
		v_cn_z_i.clear(); //re-initialize 
		which(cn, nrow, s_temp, v_cn_z_i); //Note: Actual location is returned
		int i_size_v_cn_z_i = (int)v_cn_z_i.size(); //number of locations in cn having s_temp
		
		//----------------------
		//Condition 1: this row's cells are all missing
		//  this case, all observed rows are used for imputation 
		//-----------------------
		if(i_sum_x == 0)
		{
			//-----------------
			//make "zid" which means 
			//the row location of current missing row repeated by number of all observed rows
			//-----------------
			//zid = NULL; //re-initialize; 
			i_size_zid = i_size_v_cn_z_i*nr2;
			int* zid = new int[i_size_zid]; 
			
			for(int j=0; j<i_size_v_cn_z_i; j++)
			{
				for(int k=0; k<nr2; k++) //repeated copy of the id number 
				{
					//NOTE: zid contains ACTUAL id number 
					zid[j*nr2 + k] = id[v_cn_z_i[j]-1]; //-1 for actual location
				}
			}
			
			//-------------------
			//make a matrix that consists of zid & repeated uox for all missing rows
			//-------------------
			//double** rst_temp = New_dMatrix(i_size_zid, ncol+1);
			double* rst_temp = new double[ncol];
			Fill_dVector(rst_temp, ncol, 0.0); //re-initialize

			for(int j=0; j<i_size_v_cn_z_i; j++) //all rows that have the current missing pattern
			{
				for(int k=0; k<nr2; k++) //repeated copy  
				{
					i_loc = j*nr2 + k; //serial number of the entire rows of the matrix
					
					//first column is zid[]
					//rst_temp[i_loc][0] = zid[i_loc];
					rst_id.push_back(zid[i_loc]);
					//from the second through the end columns are occupied with uox
					for(int i_col=0; i_col<ncol; i_col++)
					{
						rst_temp[i_col] = uox[k][i_col];
					}

					Trans1(rst_temp, ncol, uox1);
					rst.push_back(uox1);
				}
			}
			//---
			//Append the entire matrix to rst
			//---
			//rst.bind_blocks(i_size_zid, ncol+1, rst_temp);
			
			//testout
			/*
			RPrint(" ==Condition1. all cells on a row are missing. i: "); RPrint(i);
			RPrint("zid:"); RPrint(zid, i_size_zid);
			RPrint("rst:");
			rst.print_rbind_FHDI(); 
			*/
			
			//---------
			//local deallocation
			//---------
			//Del_dMatrix(rst_temp, i_size_zid, ncol+1);
			delete[] rst_temp;
			delete[] zid;
		}
		
		//----------------------
		//Condition 2: some cells of current row are missing
		//-----------------------
		int nl = 0;
		if(i_sum_x > 0)
		{
			//------
			//number of observed cells on this row
			//------
			nl = 0; 
			for(int j=0; j<ncol; j++) 
			{
				if(mox[i][j]>0) nl++; 
			}
			
			//-------
			//indicator matrix that matches the donors
			//srst: row-wise sum of the indicator matrix 
			//algorithm: 
			// current row   = a01 
			// observed rows = a11, ab1, af1, ... will be selected and stored  
			//-------
			loc_srst_nl.clear(); //re-initialize
			Fill_iVector(i_srst, nr2, 0); //re-initialize 
				
			for(int j=0; j<nr2; j++)
			{
				int i_sum_crst = 0; 
				for(int k=0; k<ncol; k++)
				{
					//Note: in below check, mox is fixed at ith row 
					if(fabs(mox[i][k] - uox[j][k])<1e-3) //part of missing cell = obserbed cell 
					{
						i_sum_crst++; // increment if a cell of missing row = obs. cell 
					}
				}
				//---
				//store how many cells of missing row match those of observed row
				//---
				i_srst[j] = i_sum_crst; 
				if(i_sum_crst==nl) loc_srst_nl.push_back(j+1); //Actual location 				
			}
			//testout
			/*
			RPrint(" ==Condition2. some cells on a row are missing. i: "); RPrint(i);
			RPrint("nl: "); RPrint(nl);
			RPrint("srst: "); RPrint(i_srst, nr2);
			RPrint("loc_srst_nl: "); RPrint(loc_srst_nl);
			*/
			
			//-----
			//total matching rows
			//-----
			const int i_size_loc_srst_nl = (int)loc_srst_nl.size(); 
			if(i_size_loc_srst_nl == 0) //error case
			{cout<<"Error! there is no matched cell!"<<endl; return;}
			
			if(i_size_loc_srst_nl > 0) 
			{
				//-----------------
				//make "zid" which means 
				//the row location of current missing row repeated by number of observed rows
				//-----------------
				//zid = NULL; //re-initialize; 
				i_size_zid = i_size_v_cn_z_i*i_size_loc_srst_nl;
				int* zid = new int[i_size_zid]; 
			
				for(int j=0; j<i_size_v_cn_z_i; j++)
				{
					for(int k=0; k<i_size_loc_srst_nl; k++) //repeated copy of the id number 
					{
						//NOTE: zid contains ACTUAL id number 
						zid[j*i_size_loc_srst_nl + k] = id[v_cn_z_i[j]-1]; //-1 for actual location
					}
				}
			
				//-------------------
				//make a matrix that consists of zid & repeated uox for all missing rows
				//-------------------
				//double** rst_temp2 = New_dMatrix(i_size_zid, ncol+1);
				double* rst_temp2 = new double[ncol];
				Fill_dVector(rst_temp2, ncol, 0.0); //re-initialize
				 
				for(int j=0; j<i_size_v_cn_z_i; j++)
				{
					for(int k=0; k<i_size_loc_srst_nl; k++) //repeated copy of the id number 
					{
						i_loc = j*i_size_loc_srst_nl + k; //serial number of the entire rows of the matrix
						
						//first column is zid[]
						rst_id.push_back(zid[i_loc]);
					
						//from the second through the end columns are occupied with uox
						for(int i_col=0; i_col<ncol; i_col++)
						{
							//rst_temp2[i_loc][i_col+1] = uox[k][i_col]; //cf. condition 1 form
							//rst_temp2[i_loc][i_col+1] = uox[loc_srst_nl[k]-1][i_col]; //-1 for actual location
							rst_temp2[i_col] = uox[loc_srst_nl[k] - 1][i_col]; //-1 for actual location
						}

						Trans1(rst_temp2, ncol, uox1);
						rst.push_back(uox1);
					}
				}
				//---
				//Append the entire matrix to rst
				//---
				//rst.bind_blocks(i_size_zid, ncol+1, rst_temp2);
				
				//testout
				/*
				RPrint("zid:"); RPrint(zid, i_size_zid);
				RPrint("rst:");
				rst.print_rbind_FHDI(); 
				*/
				
				//---------
				//local deallocation
				//---------
				//Del_dMatrix(rst_temp2, i_size_zid, ncol+1);	
				delete[] rst_temp2;
				delete[] zid;
			}
		}

	} //end of LOOP for all missing rows
	
	//----------------
	//re-order rst in terms of id (the first column)
	//----------------
	const int n_row_rst = rst.size(); 
	int* i_rst_id = new int[n_row_rst];
	//for(int i=0; i<n_row_rst; i++) i_rst_id[i] = (int)rst(i,0);
	for (int i = 0; i<n_row_rst; i++) i_rst_id[i] = rst_id[i];
	order_FHDI_binary(i_rst_id, n_row_rst); //returned with the order of rows in ascending magnitude
	//testout
	//RPrint("n_row_rst :"); RPrint(n_row_rst);
	//RPrint("i_rst_id :"); RPrint(i_rst_id, n_row_rst);

    
	//--------------------
	//remove the first column with id
	//store the rst into the final storage
	//--------------------
	//rbind_FHDI rst_final(ncol); //Note: without the first column of id
	//double* d_row_rst 		= new double[ncol+1]; 
	//double* d_row_rst_short = new double[ncol];
	for(int i=0; i<n_row_rst; i++)
	{
		//rst.get_block(i_rst_id[i]-1, d_row_rst); //get a row// -1 for actual loc
		//for(int k=0; k<ncol; k++) d_row_rst_short[k] = d_row_rst[k+1]; //without id  
		//rst_final.append_block(d_row_rst_short);	//append a new row to the final storage 
		rst_final.push_back(rst[i_rst_id[i] - 1]);
	}	

	//testout
	//RPrint("End of AUGMAT =========="); 
	//RPrint("rst_final:"); rst_final.print_rbind_FHDI(); 
	
	//-------
	//local deallocation
	//-------
	delete[] i_temp_x;
	//delete[] zid;	
	delete[] i_srst;
	delete[] i_rst_id;
	//delete[] d_row_rst;
	//delete[] d_row_rst_short;	
	
	return;
}