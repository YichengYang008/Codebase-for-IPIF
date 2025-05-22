void Extract_Variance_Results_MPI(const int nrow, const int ncol, 
							 double ** y_bar_i_k_Summary,
							 double* final_variance_data)
                             
//Description========================================
//  Calculate Jackknife Variance Estimator of the mean estimator
//
//  Algorithm: 
//  step1: remove kth data vector (k=1, nrow)
//  step2: calculate kth mean estimate
//  yi^(k) = sum( wi^(k)*wij^(k)*y_ij)/sum(wi^(k)*wij^(k))
//  where
//  yi = final vector corresponding to the ith row of original data
//     = {yi1, yi2, ..., yi_ncol}
//  wi^(k) = weight of the ith row of kth replicate of data
//  wij^(k) = fractional weight of the jth imputed cell on the ith row of kth replicate
//  y_ij = j_th imputed cell for the missing cell 
//
//  Written by Dr. I. Cho
//  updated Feb, 07, 2017. 
//
//IN   : rbind_FHDI  rbind_ipmat(4+ncol) //column size is 4+ncol (i.e., for R: ID, FID, WGT, FWGT, Variables)
//IN   : int nrow = total rows of original data matrix
//IN   : double final_full_data(nrow*ncol) = full matrix (in vector form) with original data and imputed values for missing parts
//IN   : rbind_FHDI  rbind_vrst(nrow) = fractional weights matrix from Jackknife Var Est.
//                                      row number is the same as that of ipmat    
//OUT  : double final_variance_data(ncol) = Jackknife var est vector 
//                                          Variable-wise Variance estimates
//===================================================  
{
	//int i_loc = 0; //sequential index of global total rows 
	//double* yi = new double[ncol]; //column-wise mean 
	//Fill_dVector(final_variance_data, ncol, 0.0); 
	//
	////--------
	////new mean of k_th replicate of data y for Jackknife Var Est
	////--------
	//double** y_bar_i_k = New_dMatrix(nrow, ncol); // for nrow replicates for ncol variables replicates
	//Fill_dMatrix(y_bar_i_k, nrow, ncol, 0.0);
	////cout << "nrow is ----------" << nrow << ", ncol is: " <<ncol << endl;
	//for(int k=0; k<nrow; k++) //Jackknife replicates 
	//{	
	//	//---
	//	//re-initialization! for new jackknifed data
	//	//---
	//	i_loc = 0; 
	//	double d_sum_wij = 0.0; 
	//	Fill_dVector(yi, ncol, 0.0); //initialize vector for column-wise means of all variables  
	//	for(int i=0; i<nrow; i++)
	//	{
	//		//-----
	//		//inner loop within the identical ID 
	//		//-----
	//		for(int j=0; j<nrow; j++) 
	//		{
	//			int ID = (int)rbind_ipmat(i_loc, 0) - 1 ; //1st col means ID; "-1" for actual location  
	//			
	//			//testout
	//			//cout<<"k, i, ID: "<<k<<", "<<i<<" , "<<ID<<endl;
	//			
	//			if(ID == i) //as long as the same ID 
	//			{
	//				double wi = rbind_ipmat(i_loc, 2); 
	//				//double wij= rbind_ipmat(i_loc, 3); //used in mean calculation
	//				
	//				//-------
	//				//NOTE: use the replicated fractional weight in lieu of ipmat
	//				//-------
	//				double wij = rbind_vrst(i_loc, k); //k means current Jackknifed column
	//		        //cout << "k value is: " << k << endl;
	//				//cout << "i_loc is: " << i_loc << endl;
	//				//----
	//				//accumulate fractional weight
	//				//variable-wise weight summation
	//				//----
	//				d_sum_wij += wi*wij; //WGT*FWGT
	//				//----
	//				//do weighted summation with all the imputed cells 
	//				//for current row 
	//				//----
	//				for(int i_var=0; i_var<ncol; i_var++)
	//				{
	//					yi[i_var] = yi[i_var] + wi*wij * rbind_ipmat(i_loc,4+i_var);	
	//				}	
	//		
	//				//----
	//				//increment location for next row 
	//				//----
	//				//cout <<"K: "<<k<<", I_loc inside i: " << i_loc << ", yicheng_counter is: " << yicheng_counter << endl;
	//				//cout << "j: " << j << ", i: " << i <<", ID is: "<<ID<< ", i_loc: " << i_loc << endl;
	//				i_loc++; 
	//			}
	//			
	//			if(ID > i) {break;} //exit inner loop 
	//		}
	//		
	//	} //end of main loop of i = [0, nrow)
	//	//cout << "End inside i loop of i_loc" << endl;
	//	//cout << "I_loc inside k: " << i_loc << endl;
	//	if(fabs(d_sum_wij)== 0.0) 
	//	{
	//		cout<<"ERROR! zero sum of fractional weight at Jackknifed row :"<< k<<endl;
	//		return;
	//	}
	//	
	//		
	//	//----
	//	//store Jackknife mean estimator 
	//	//----
	//	for(int i_var=0; i_var< ncol; i_var++) //size of columns of ipmat matrix of C++
	//	{
	//		//cout << "K: " << k << ", i_loc: " << i_loc << ", d_sum_wij: " << d_sum_wij << endl;
	//		double d_temp = yi[i_var]/d_sum_wij; 
	//			
	//		//-----
	//		//bar{y}_i^(k)
	//		//-----
	//		//NOTE: R works column-by-column 
	//		//hence, below sequence is different from C++ ordering 
	//		//final_full_data[i_var*nrow + i] = d_temp;  //note: i=current row
	//		y_bar_i_k[k][i_var] = d_temp;  //note: k= replicate; i_var = variable
	//	}
	//	
	//}//end of Jackknifed data k = [0, nrow)
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------
	//---------------------
	//Jackknife average of y_bar_i
	//---------------------
	double* y_bar_i = new double[ncol]; 
	for(int i_var=0; i_var<ncol; i_var++)
	{
		double d_temp =0.0; 
		for(int k=0; k<nrow; k++)
		{
			d_temp += y_bar_i_k_Summary[k][i_var];
			//cout << "k value is: " << k << endl;
		}
		y_bar_i[i_var] = d_temp/nrow; 
	}
	//testout
	//cout<<"y_bar_i[0] "<<y_bar_i[0]<<endl;
	
	
	//---------------------
	//Final Jackknife Variance Estimation of the mean 
	//---------------------
	for(int i_var=0; i_var< ncol; i_var++) //size of columns of ipmat matrix of C++
	{
		double d_temp = 0; 
		d_temp = 0.0; 
		for(int k=0; k<nrow; k++) //Jackknife replicate
		{
			d_temp += 
				(y_bar_i_k_Summary[k][i_var] - y_bar_i[i_var])
			   *(y_bar_i_k_Summary[k][i_var] - y_bar_i[i_var]);
		}

		//-----
		//final jackknife variance estimate
		//-----
		final_variance_data[i_var] = d_temp * (nrow-1)/nrow; 
	}
	
	//----
	//deallocation
	//----
	//delete[] yi; 
	//Del_dMatrix(y_bar_i_k_Summary, nrow, ncol);
	delete[] y_bar_i; 
	
	return; 
}
