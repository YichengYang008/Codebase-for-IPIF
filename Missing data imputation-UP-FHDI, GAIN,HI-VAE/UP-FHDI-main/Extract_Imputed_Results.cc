void Extract_Imputed_Results(const int nrow, const int ncol, rbind_FHDI &rbind_ipmat,
                             double* final_full_data)
                             
//Description========================================
//  extract imputed values from the ipmat 
//  with weights and fractional weights
//  
//  Algorithm: 
//  yi = sum( wi*wij*y_ij)/sum(wi*wij)
//  where
//  yi = final vector corresponding to the ith row of original data
//     = {yi1, yi2, ..., yi_ncol}
//  wi = weight of the ith row
//  wij = fractional weight of the jth imputed cell on the ith row
//  y_ij = j_th imputed cell for the missing cell 
//
//  Written by Dr. I. Cho
//  updated Feb, 03, 2017. 
//
//IN   : rbind_FHDI  rbind_ipmat(4+ncol) //column size is 4+ncol (i.e., for R: ID, FID, WGT, FWGT, Variables)
//IN   : int nrow = total rows of original data matrix
//OUT  : double final_full_data(nrow*ncol) = full matrix (in vector form) with original data and imputed values for missing parts
//===================================================  
{
	int i_loc = 0; 
	double* yi = new double[ncol]; 
	Fill_dVector(final_full_data, nrow*ncol, 0.0); 
	
	//--------
	//main loop for all rows of original data matrix
	//--------
	for(int i=0; i<nrow; i++)
	{
		//-----
		//inner loop within the identical ID 
		//-----
		double d_sum_wij = 0.0; 
		Fill_dVector(yi, ncol, 0.0); //initialize ith row data vector 
		for(int j=0; j<nrow; j++) 
		{
			int ID = (int)rbind_ipmat(i_loc, 0) - 1 ; //1st col means ID; "-1" for actual location  
			if(ID == i) //as long as the same ID 
			{
				double wi = rbind_ipmat(i_loc, 2); 
				double wij= rbind_ipmat(i_loc, 3);
				
				//----
				//accumulate fractional weight
				//----
				d_sum_wij += wi*wij; //WGT*FWGT

				//----
				//do weighted summation with all the imputed cells 
				//for current row 
				//----
				for(int i_var=0; i_var<ncol; i_var++)
				{
					yi[i_var] = yi[i_var] + wi*wij * rbind_ipmat(i_loc,4+i_var);	
				}	
		
				//----
				//increment location for next row 
				//----
				i_loc++; 
			}
			
			if(ID > i) {break;} //exit inner loop 
		}
		
		if(fabs(d_sum_wij)== 0.0) 
		{
			cout<<"ERROR! zero sum of fractional weight at the row: "<<i<<endl;
			return;
		}
		
		//-----------
		//calculate the final data matrix
		//-----------
		/*
		for(int i_var=0; i_var<ncol; i_var++)
		{
			final_full_data[i][i_var] = yi[i_var]/d_sum_wij;  
		}
		*/

		//----
		//copy Resp of C++ into return storage
		//----
		for(int i_var=0; i_var< ncol; i_var++) //size of columns of ipmat matrix of C++
		{
			//for(int j=0; j<nrow; j++)
			//{
				double d_temp = yi[i_var]/d_sum_wij; 
				
				//NOTE: R works column-by-column 
				//hence, below sequence is different from C++ ordering 
				final_full_data[i_var*nrow + i] = d_temp;  //note: i=current row
			//}
		}		
		
	} //end of main loop of i = [0, nrow)

	
	//----
	//deallocation
	//----
	delete[] yi; 
	
	return; 
}
