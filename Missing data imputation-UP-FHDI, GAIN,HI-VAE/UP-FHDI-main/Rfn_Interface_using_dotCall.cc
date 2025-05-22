#define R_NO_REMAP //so that R does not define length(x) which may cause many complie error with fstream

#include <R.h>
#include <Rinternals.h>
#include "Rfn_test.cc" 
#include "Extract_Imputed_Results.cc"
#include "Extract_Variance_Results.cc"


//----------------------------------------------------------------------------
//- This file defines TWO functions. 
// 
// The first one is a C++ function that directly accesses 
// the C++ class. 
// 
// The second one is a C function named ’CWrapper’, that INdirectly accesses 
// the C++ class, via a call to the function.
//----------------------------------------------------------------------------
void Rfn_test_call(double* x, int* r, int * nrow_x, int * ncol_x, 
                   double* k, double* d, int * M, 
				   int * i_option_imputation, int * i_option_variance, 
				   int * id, double* z_UserDefined,
				   rbind_FHDI &rbind_ipmat_FEFI_return,
				   rbind_FHDI &rbind_Resp_FEFI_return, 
				   rbind_FHDI &rbind_irmat_FEFI_return,
				   rbind_FHDI &rbind_ipmat_FHDI_return,
				   rbind_FHDI &rbind_Resp_FHDI_return, 
				   rbind_FHDI &rbind_irmat_FHDI_return,
				   rbind_FHDI &rbind_vrst_FEFI_return, 
				   rbind_FHDI &rbind_vrst_FHDI_return, 
				   
				   rbind_FHDI  &rbind_uox_return,
				   rbind_FHDI  &rbind_mox_return,
				   rbind_FHDI  &rbind_category_return, 
				   
				   std::vector<std::string> &jp_name_return_CellProb,
				   std::vector<double> &jp_prob_return_CellProb,
				   
				   const int i_option_perform, 
				   int* i_option_merge)
{ 

	Rfn_test(x, r, nrow_x, ncol_x, k, d, M, 
	         i_option_imputation, i_option_variance, 
			 id, z_UserDefined,
			 rbind_ipmat_FEFI_return, rbind_Resp_FEFI_return, rbind_irmat_FEFI_return,
			 rbind_ipmat_FHDI_return, rbind_Resp_FHDI_return, rbind_irmat_FHDI_return,
			 rbind_vrst_FEFI_return,  rbind_vrst_FHDI_return,
	         
			 rbind_uox_return, rbind_mox_return, rbind_category_return, 
			 
			 jp_name_return_CellProb, jp_prob_return_CellProb, 
			 
			 i_option_perform,
			 i_option_merge); 
	
	return ;

}
//===============================================================
//===============================================================
//===============================================================
//===============================================================
//===============================================================
extern "C" {
SEXP CWrapper(SEXP x_R, SEXP r_R, SEXP z_R, SEXP i_option_perform_R, 
              SEXP nrow_x_R, SEXP ncol_x_R, 
              SEXP k_R, SEXP d_R, SEXP M_R, 
			  SEXP i_option_imputation_R, SEXP i_option_variance_R, 
			  SEXP id_R, SEXP i_option_merge_R)
//Description -------------------------------------------
//----------------------------------------------------------------------------
//- CWrapper 
//  to perform 
//         (1) Cell Make, (2)Cell Prob, (3) FEFI/FHDI imputation, and 
//         (4) Jackknife Var Est. (if requested)  
// 
// This C function helps invoke the deeper C++ functions. 
// R can access C code, not C++ code directly. 
// This C function provides a C-C++ interface so that R can access the C++
// Note: Rinternals.h is used. In this case, to prevent R's automatic function mapping
//       at the very beginning, #define R_NO_REMAP should be added and also
//       "Rf_" should precede all R functions used herein
// Written by Dr. In Ho Cho
// All rights reserved 
// Jan 17, 2017
//---------------------------------------------------------------------------- 
//IN   : double x_R[nrow_x_R * ncol_x_R] = data matrix in a vector form with missing units
//IN   : int    r_R[nrow_x_R * ncol_x_R] = indices for data missingness 1=observed; 0=missing 			  
//IN   : double z_R[nrow_x_R * ncol_x_R] = user-defined Z matrix (i_option_perform =4) only
//IN   : int    i_option_perform_R = main control option (1: all; 2: CellMake; 3: CellProb
//                                                        4: all by using user-defined z)
//IN   : double k_R[ncol_x_R] = a vector of categories as an initial (can be known for discrete variables) 
//IN   : double d_R[nrow_x_R] = (sampling) weights for units
//IN   : int    M_R	= imputation size (can be generaized by M_i, i denotes missing unit)
//IN   : int 	i_option_imputation_R = 1: Fully Efficient Fractional Imputation
// 									    2: Fractional Hot Deck Imputation
//IN   : int    i_option_variance_R  = 0: skip variance estimation process
//                                      1: perform variance estimation using Jackknife method
//IN   : int    id_R[nrow_x_R] = id numbers of all data 
//IN   : int    i_option_merge = random donor selection in Merge algorithm in Cell Make
//                               0: no 
//                               1: random seed for donor selection
//OUT  : List of 
//       rbind_ipmat(4+ncol) // ID, FID, WGT, FWGT, Imputed Variables, Response indices (i.e., rbind_Resp(ncol+1)) 
//       cured data matrix(nrow, ncol)
//       Mean and Standard Error (2,ncol)
//       rbind_vrst(nrow)    // Jackknife variance estimates  (returned when i_option_variance_R = 1)
//-------------------------------------------------------------------------------
{
	//testout
	//RPrint("Begin CWrapper in Rfn_Interface... ==== "); 

	x_R = PROTECT(Rf_coerceVector(x_R, REALSXP));
	double *x = REAL(x_R);         	//pointer to double vector x[col*row] that contains all data with missing units
	r_R = PROTECT(Rf_coerceVector(r_R, INTSXP));
	int    *r = INTEGER(r_R); 			//pointer to an integer vector r[n_total_x] that contains indices of 0 and 1
	//testout
	//RPrint("in Rfn_Interface... x:  "); 
	//RPrint(x[0]);
	//testout
	//RPrint("in Rfn_Interface... r:  "); 
	//RPrint(r[0]);

	i_option_perform_R = PROTECT(Rf_coerceVector(i_option_perform_R, INTSXP));
	int    *i_option_perform = INTEGER(i_option_perform_R); 			

	z_R = PROTECT(Rf_coerceVector(z_R, REALSXP));
	double *z_UserDefined = REAL(z_R); 
	
	nrow_x_R = PROTECT(Rf_coerceVector(nrow_x_R, INTSXP));
	ncol_x_R = PROTECT(Rf_coerceVector(ncol_x_R, INTSXP));
	int    *nrow_x = INTEGER(nrow_x_R);
    int	   *ncol_x = INTEGER(ncol_x_R); //pointers to integer sizes of x
	//testout
	//RPrint("in Rfn_Interface... nrow_x:  "); 
	//RPrint(nrow_x[0]);
	
	k_R = PROTECT(Rf_coerceVector(k_R, REALSXP)); 
	double *k = REAL(k_R);			//pointer to a double value 
	d_R = PROTECT(Rf_coerceVector(d_R, REALSXP)); 
	double *d = REAL(d_R);  		//pointer to double vector d[nrow_x] 
	M_R = PROTECT(Rf_coerceVector(M_R, INTSXP)); 
	int    *M = INTEGER(M_R); 		//pointer to integer number of donors 

	id_R = PROTECT(Rf_coerceVector(id_R, INTSXP)); 
	int    *id = INTEGER(id_R); 	//pointer to integer number of id of data
	
	i_option_imputation_R = PROTECT(Rf_coerceVector(i_option_imputation_R, INTSXP)); 
	int    *i_option_imputation = INTEGER(i_option_imputation_R); 		
	i_option_variance_R = PROTECT(Rf_coerceVector(i_option_variance_R, INTSXP)); 
	int    *i_option_variance = INTEGER(i_option_variance_R); 		

	i_option_merge_R = PROTECT(Rf_coerceVector(i_option_merge_R, INTSXP)); 
	int    *i_option_merge = INTEGER(i_option_merge_R); 		
	
	//UP TO here 13 "protect" as of Feb 2017
	
	//----------------
	//Basic Error Check
	//as of Feb 3, 2017
	//----------------
	//M
	if(M[0]<1) {RPrint("Error! M is less than 1 "); return(R_NilValue);}
	if(M[0]>nrow_x[0]) {RPrint("Error! M is larger than total rows of data "); return(R_NilValue);}
	//k
	for(int i=0; i<ncol_x[0]; i++)
	{
		if(k[i] < 1){RPrint("Error! some k is less than 1 "); return(R_NilValue);}
		if(k[i] > 35){RPrint("Error! some k is larger than 35 "); return(R_NilValue);}
	}
	
	//--------
	//prep return variables
	//--------
	rbind_FHDI  rbind_ipmat_FEFI_return(4+ncol_x[0]); //column size is 4+ncol
	rbind_FHDI  rbind_Resp_FEFI_return(ncol_x[0]+1);  //separate response matrix 
	rbind_FHDI  rbind_irmat_FEFI_return(5+ncol_x[0]); //column size is 5+ncol    
	rbind_FHDI  rbind_ipmat_FHDI_return(4+ncol_x[0]); //column size is 4+ncol
	rbind_FHDI  rbind_Resp_FHDI_return(ncol_x[0]+1);  //separate response matrix 
	rbind_FHDI  rbind_irmat_FHDI_return(5+ncol_x[0]); //column size is 5+ncol    
	rbind_FHDI  rbind_vrst_FEFI_return(nrow_x[0]);    //variance estimates of FEFI
	rbind_FHDI  rbind_vrst_FHDI_return(nrow_x[0]);    //variance estimates of FHDI

	//below is for output for Cell Make only option
	rbind_FHDI  rbind_uox_return(ncol_x[0]); //unique observed patterns
	rbind_FHDI  rbind_mox_return(ncol_x[0]); //unique observed patterns
	rbind_FHDI  rbind_category_return(ncol_x[0]); //cagetorized matrix 
	
	//below is for output for Cell Prob only option
	std::vector<std::string> jp_name_return_CellProb;   //name of the joint probability table
	std::vector<double> jp_prob_return_CellProb; //the latest joint probability 	
		
	//=====================
	//=====================
	//*********************
	//=====================
	//=====================
	int i_op_p_temp = 1; 
	if(i_option_perform[0] == 4){ i_op_p_temp = 4;} 
	Rfn_test_call(x, r, nrow_x, ncol_x, k, d, M, 
	              i_option_imputation, i_option_variance, id, z_UserDefined,
				  rbind_ipmat_FEFI_return, rbind_Resp_FEFI_return, rbind_irmat_FEFI_return,
				  rbind_ipmat_FHDI_return, rbind_Resp_FHDI_return, rbind_irmat_FHDI_return,
				  rbind_vrst_FEFI_return,  rbind_vrst_FHDI_return, 
				  
				  rbind_uox_return, rbind_mox_return, rbind_category_return,
				  jp_name_return_CellProb, jp_prob_return_CellProb,
				  
				  i_op_p_temp, i_option_merge); //1: perform Entire FEFI/FHDI  
	

	//testout
	//RPrint("in Rfn_Interface... Rfn_test_call has finished  "); 
				  
	//----------
	//FEFI: copy the calculated matrix
	//----------
	if(i_option_imputation[0] == 1) //FEFI
	{
		//--------
		//ipmat return
		//--------
		//const int i_FEFI_ipmat_size_col = rbind_ipmat_FEFI_return.size_col();
		const int i_FEFI_ipmat_size_row = rbind_ipmat_FEFI_return.size_row();
		const int n_ipmat_Resp_total = (4+ncol_x[0]) ; //Feb 7, 2017 		
		SEXP return_ipmat_FEFI = PROTECT(Rf_allocMatrix(REALSXP, i_FEFI_ipmat_size_row, n_ipmat_Resp_total));
		double* d_return_ipmat_FEFI = REAL(return_ipmat_FEFI);
		//----
		//copy ipmat of C++ first into return storage
		//----
		for(int i=0; i< 4+ncol_x[0]; i++) //size of columns of ipmat matrix of C++
		{
			for(int j=0; j<i_FEFI_ipmat_size_row; j++)
			{
				double d_temp = 0.0; 
				d_temp = rbind_ipmat_FEFI_return(j,i); //from ipmat of C++

				//NOTE: R works column-by-column 
				//hence, below sequence is different from C++ ordering 
				d_return_ipmat_FEFI[i*i_FEFI_ipmat_size_row + j] = d_temp; //get all stored values 
			}
		}

		//----------------------
		//----------------------
		//when FEFI's variance estimation is NOT required
		//----------------------
		//----------------------
		if(i_option_variance[0] == 0)
		{
			//----
			//create list for compact return
			//-----
			SEXP list_return;
			PROTECT(list_return = Rf_allocVector(VECSXP,2));

			//-------
			//extract the final full data matrix
			//-------
			SEXP Final_Full_Data = PROTECT(Rf_allocMatrix(REALSXP, nrow_x[0], ncol_x[0]));
			double* d_Final_Full_Data = REAL(Final_Full_Data);		
			Extract_Imputed_Results(nrow_x[0], ncol_x[0], rbind_ipmat_FEFI_return,
									d_Final_Full_Data); 

			//--------
			//Make a final retuan LIST
			//--------
			SET_VECTOR_ELT(list_return, 0, return_ipmat_FEFI);
			SET_VECTOR_ELT(list_return, 1, Final_Full_Data);
			//Rf_setAttrib(list_return, R_NamesSymbol, list_names_return);
			
			
			//----
			//FEFI return
			//----
			UNPROTECT(13+3);
			
			return list_return; 
		}

		//--------------------
		//--------------------
		//when FEFI's variance estimation IS required
		//--------------------
		//--------------------
		if(i_option_variance[0] == 1) 
		{
			const int i_FEFI_vrst_size_row = rbind_vrst_FEFI_return.size_row();
			SEXP return_vrst_FEFI = PROTECT(Rf_allocMatrix(REALSXP, i_FEFI_vrst_size_row, nrow_x[0]));
			double* d_return_vrst_FEFI = REAL(return_vrst_FEFI);
			//----
			//copy vrst of C++ into return storage
			//----
			for(int i=0; i< nrow_x[0]; i++) //size of COLUMNs of vrst matrix of C++
			{
				for(int j=0; j<i_FEFI_vrst_size_row; j++) //rows
				{
					double d_temp = 0.0; 
					d_temp = rbind_vrst_FEFI_return(j,i); //from vrst of C++

					//NOTE: R works column-by-column 
					//hence, below sequence is different from C++ ordering 
					d_return_vrst_FEFI[i*i_FEFI_vrst_size_row + j] = d_temp; //get all stored values 
				}
			}

			//----
			//create list for compact return
			//-----
			SEXP list_return;
			PROTECT(list_return = Rf_allocVector(VECSXP,4));

			//-------
			//extract the final full data matrix
			//-------
			SEXP Final_Full_Data = PROTECT(Rf_allocMatrix(REALSXP, nrow_x[0], ncol_x[0]));
			double* d_Final_Full_Data = REAL(Final_Full_Data);		
			Extract_Imputed_Results(nrow_x[0], ncol_x[0], rbind_ipmat_FEFI_return,
									d_Final_Full_Data); 
									
			//-------
			//extract the final variance data matrix
			//-------
			SEXP Final_Mean_SD_Data = PROTECT(Rf_allocMatrix(REALSXP, 2, ncol_x[0])); //mean & std. dev
			double* d_Final_Mean_SD_Data = REAL(Final_Mean_SD_Data); 
			
			//Variance calculation 
			double* d_Final_Variance_Data = new double[ncol_x[0]];
			Extract_Variance_Results(nrow_x[0], ncol_x[0], 
									rbind_ipmat_FEFI_return,
									d_Final_Full_Data, 
									rbind_vrst_FEFI_return,
									d_Final_Variance_Data); 						

            //Column-wise mean calculation 
			double* d_Final_Column_Mean = new double[ncol_x[0]];
			for(int i_var=0; i_var< ncol_x[0]; i_var++) { //size of columns of ipmat matrix of C++
				double d_temp = 0.0; 
				for(int i=0; i<nrow_x[0]; i++) d_temp += d_Final_Full_Data[i_var*nrow_x[0] + i] ;  //note: i=current row
				d_Final_Column_Mean[i_var] = d_temp/nrow_x[0]; 
			}
			//store mean and standard error for return 
			Fill_dVector(d_Final_Mean_SD_Data, 2*ncol_x[0], 0.0); //initialize 
			for(int i_var=0; i_var<ncol_x[0]; i_var++){
				d_Final_Mean_SD_Data[i_var*2]   = d_Final_Column_Mean[i_var]; //mean of column i_var
				if(d_Final_Variance_Data[i_var]>0.0) d_Final_Mean_SD_Data[i_var*2+1] = sqrt(d_Final_Variance_Data[i_var]);
			}
			//deallocate local array
			delete[] d_Final_Variance_Data;
			delete[] d_Final_Column_Mean; 
			
			//---------------
			//Make a final return LIST 
			//---------------
			SET_VECTOR_ELT(list_return, 0, return_ipmat_FEFI);
			SET_VECTOR_ELT(list_return, 1, Final_Full_Data);
			SET_VECTOR_ELT(list_return, 2, Final_Mean_SD_Data);
			SET_VECTOR_ELT(list_return, 3, return_vrst_FEFI);
			
			//Rf_setAttrib(list_return, R_NamesSymbol, list_names_return);
			
			//----
			//FEFI return
			//----
			UNPROTECT(13+5);			

			return list_return; 
		}
		
	}

	//----------
	//FHDI: copy the calculated matrix
	//----------
	if(i_option_imputation[0] == 2) //FHDI
	{
		//--------
		//ipmat return
		//NOTE: ipmat + Resp for comprehensive matrix in R 
		//--------
		//const int i_FHDI_ipmat_size_col = rbind_ipmat_FHDI_return.size_col();
		const int i_FHDI_ipmat_size_row = rbind_ipmat_FHDI_return.size_row();

        const int n_ipmat_Resp_total = (4+ncol_x[0]); //Feb 7, 2017		
		SEXP return_ipmat_FHDI = PROTECT(Rf_allocMatrix(REALSXP, i_FHDI_ipmat_size_row, n_ipmat_Resp_total));
		double* d_return_ipmat_FHDI = REAL(return_ipmat_FHDI);
		//----
		//copy ipmat of C++ first into return storage
		//----
		for(int i=0; i< 4+ncol_x[0]; i++) //size of columns of ipmat matrix of C++
		{
			for(int j=0; j<i_FHDI_ipmat_size_row; j++)
			{
				double d_temp = 0.0; 
				d_temp = rbind_ipmat_FHDI_return(j,i); //from ipmat of C++

				//NOTE: R works column-by-column 
				//hence, below sequence is different from C++ ordering 
				d_return_ipmat_FHDI[i*i_FHDI_ipmat_size_row + j] = d_temp; //get all stored values 
			}
		}

		//----------------------
		//----------------------
		//when FHDI's variance estimation is NOT required
		//----------------------
		//----------------------
		if(i_option_variance[0] == 0)
		{		
			//----
			//create list for compact return
			//-----
			SEXP list_return;

			PROTECT(list_return = Rf_allocVector(VECSXP,2));
			//-------
			//extract the final full data matrix
			//-------
			SEXP Final_Full_Data = PROTECT(Rf_allocMatrix(REALSXP, nrow_x[0], ncol_x[0]));
			double* d_Final_Full_Data = REAL(Final_Full_Data);		
			Extract_Imputed_Results(nrow_x[0], ncol_x[0], rbind_ipmat_FHDI_return,
									d_Final_Full_Data); 			
			
			SET_VECTOR_ELT(list_return, 0, return_ipmat_FHDI);
			SET_VECTOR_ELT(list_return, 1, Final_Full_Data);
			//Rf_setAttrib(list_return, R_NamesSymbol, list_names_return);			
			
			//----
			//FHDI return
			//----
			UNPROTECT(13+3);			

			return list_return; 
		}

		//--------------------
		//--------------------
		//when FHDI's variance estimation IS required
		//--------------------
		//--------------------
		if(i_option_variance[0] == 1) 
		{
			const int i_FHDI_vrst_size_row = rbind_vrst_FHDI_return.size_row();
			
			SEXP return_vrst_FHDI = PROTECT(Rf_allocMatrix(REALSXP, i_FHDI_vrst_size_row, nrow_x[0]));
			double* d_return_vrst_FHDI = REAL(return_vrst_FHDI);
			//----
			//copy vrst of C++ into return storage
			//----
			for(int i=0; i< nrow_x[0]; i++) //size of COLUMNs of vrst matrix of C++
			{
				for(int j=0; j<i_FHDI_vrst_size_row; j++) //rows
				{
					double d_temp = 0.0; 
					d_temp = rbind_vrst_FHDI_return(j,i); //from vrst of C++

					//NOTE: R works column-by-column 
					//hence, below sequence is different from C++ ordering 
					d_return_vrst_FHDI[i*i_FHDI_vrst_size_row + j] = d_temp; //get all stored values 
				}
			}

			//----
			//create list for compact return
			//-----
			SEXP list_return;

			PROTECT(list_return = Rf_allocVector(VECSXP,4));
			//-------
			//extract the final full data matrix
			//-------
			SEXP Final_Full_Data = PROTECT(Rf_allocMatrix(REALSXP, nrow_x[0], ncol_x[0]));
			double* d_Final_Full_Data = REAL(Final_Full_Data);		
			Extract_Imputed_Results(nrow_x[0], ncol_x[0], rbind_ipmat_FHDI_return,
									d_Final_Full_Data); 			

			//P16----
			//extract the final variance data matrix
			//-------
			SEXP Final_Mean_SD_Data = PROTECT(Rf_allocMatrix(REALSXP, 2, ncol_x[0])); //mean & std. dev
			double* d_Final_Mean_SD_Data = REAL(Final_Mean_SD_Data); 
			
			//Variance calculation 
			double* d_Final_Variance_Data = new double[ncol_x[0]];			
			Extract_Variance_Results(nrow_x[0], ncol_x[0], 
									rbind_ipmat_FHDI_return,
									d_Final_Full_Data, 
									rbind_vrst_FHDI_return,
									d_Final_Variance_Data); 		
									
            //Column-wise mean calculation 
			double* d_Final_Column_Mean = new double[ncol_x[0]];
			for(int i_var=0; i_var< ncol_x[0]; i_var++) { //size of columns of ipmat matrix of C++
				double d_temp = 0.0; 
				for(int i=0; i<nrow_x[0]; i++) d_temp += d_Final_Full_Data[i_var*nrow_x[0] + i] ;  //note: i=current row
				d_Final_Column_Mean[i_var] = d_temp/nrow_x[0]; 
			}
			//store mean and standard error for return 
			Fill_dVector(d_Final_Mean_SD_Data, 2*ncol_x[0], 0.0); //initialize 
			for(int i_var=0; i_var<ncol_x[0]; i_var++){
				d_Final_Mean_SD_Data[i_var*2]   = d_Final_Column_Mean[i_var]; //mean of column i_var
				if(d_Final_Variance_Data[i_var]>0.0) d_Final_Mean_SD_Data[i_var*2+1] = sqrt(d_Final_Variance_Data[i_var]);
			}
			//deallocate local array
			delete[] d_Final_Variance_Data;
			delete[] d_Final_Column_Mean; 									

			//---------------
			//final preparation of return matrices
			//---------------									
			SET_VECTOR_ELT(list_return, 0, return_ipmat_FHDI);
			SET_VECTOR_ELT(list_return, 1, Final_Full_Data );
			SET_VECTOR_ELT(list_return, 2, Final_Mean_SD_Data);
			SET_VECTOR_ELT(list_return, 3, return_vrst_FHDI);
			//Rf_setAttrib(list_return, R_NamesSymbol, list_names_return);
			
			//----
			//FHDI return
			//----
			UNPROTECT(13+5);			
			
			return list_return; 
		}		
	}
	
	//----
	//general return
	//----
	UNPROTECT(13);
	return(R_NilValue); 
}
}

//below is a version for ".C" -----------------------------------------------
/*
extern "C" 
{ 
	void CWrapper(double* x, int* r, int * nrow_x, int * ncol_x, double* k, 
	              double* d, int* M) 
	{ 
		//- Invoke second function which internally can do C++ things. 
		// 
		
		Rfn_test_call(x, r, nrow_x, ncol_x, k, d, M); 

	} 
}
*/

//===============================================================
//===============================================================
//===============================================================
//===============================================================
//===============================================================
extern "C" {
SEXP CWrapper_CellMake(SEXP x_R, SEXP r_R, SEXP nrow_x_R, SEXP ncol_x_R, 
              SEXP k_R, SEXP d_R, SEXP M_R, 
			  SEXP i_option_imputation_R, SEXP i_option_variance_R, 
			  SEXP id_R, SEXP i_option_merge_R)
//Description -------------------------------------------
//----------------------------------------------------------------------------
//- CWrapper 
//  to perform 
//         (1) Cell Make only!
// 
// This C function helps invoke the deeper C++ functions. 
// R can access C code, not C++ code directly. 
// This C function provides a C-C++ interface so that R can access the C++
// Note: Rinternals.h is used. In this case, to prevent R's automatic function mapping
//       at the very beginning, #define R_NO_REMAP should be added and also
//       "Rf_" should precede all R functions used herein
// Written by Dr. In Ho Cho
// All rights reserved 
// Feb 9, 2017
//---------------------------------------------------------------------------- 
//IN   : double x_R[nrow_x_R * ncol_x_R] = data matrix in a vector form with missing units
//IN   : int    r_R[nrow_x_R * ncol_x_R] = indices for data missingness 1=observed; 0=missing 			  
//IN   : double k_R[ncol_x_R] = a vector of categories as an initial (can be known for discrete variables) 
//IN   : double d_R[nrow_x_R] = (sampling) "w" weights for units
//IN   : int    M_R	= imputation size (can be generaized by M_i, i denotes missing unit)
//IN   : int 	i_option_imputation_R = 1: Fully Efficient Fractional Imputation
// 									    2: Fractional Hot Deck Imputation
//IN   : int    i_option_variance_R  = 0: skip variance estimation process
//                                      1: perform variance estimation using Jackknife method
//IN   : int    id_R[nrow_x_R] = id numbers of all data 
//IN   : int    i_option_merge = random donor selection in Merge algorithm in Cell Make
//                               0: no 
//                               1: random seed for donor selection
//OUT:  LIST of 
//     [[1]] ID, WGT, V1, V2, ... 
//     [[2]] categorized V1, V2, ...
//     [[3]] uox: unique patterns of observed rows sorted in the ascending order
//     [[4]] mox: unique patterns of the missing rows sorted in the ascending order
//-------------------------------------------------------------------------------
{
	//testout
	//RPrint("Begin CWrapper_CellMake in Rfn_Interface... ==== "); 

	//P1-2
	x_R = PROTECT(Rf_coerceVector(x_R, REALSXP));
	double *x = REAL(x_R);         	//pointer to double vector x[col*row] that contains all data with missing units
	r_R = PROTECT(Rf_coerceVector(r_R, INTSXP));
	int    *r = INTEGER(r_R); 			//pointer to an integer vector r[n_total_x] that contains indices of 0 and 1
	
	//P3-4
	nrow_x_R = PROTECT(Rf_coerceVector(nrow_x_R, INTSXP));
	ncol_x_R = PROTECT(Rf_coerceVector(ncol_x_R, INTSXP));
	int    *nrow_x = INTEGER(nrow_x_R);
    int	   *ncol_x = INTEGER(ncol_x_R); //pointers to integer sizes of x
	
	//P5-7
	k_R = PROTECT(Rf_coerceVector(k_R, REALSXP)); 
	double *k = REAL(k_R);			//pointer to a double value 
	d_R = PROTECT(Rf_coerceVector(d_R, REALSXP)); 
	double *d = REAL(d_R);  		//pointer to double vector d[nrow_x] 
	M_R = PROTECT(Rf_coerceVector(M_R, INTSXP)); 
	int    *M = INTEGER(M_R); 		//pointer to integer number of donors 

	//P8
	id_R = PROTECT(Rf_coerceVector(id_R, INTSXP)); 
	int    *id = INTEGER(id_R); 	//pointer to integer number of id of data
	
	//options
    //P9-10	
	i_option_imputation_R = PROTECT(Rf_coerceVector(i_option_imputation_R, INTSXP)); 
	int    *i_option_imputation = INTEGER(i_option_imputation_R); 		
	i_option_variance_R = PROTECT(Rf_coerceVector(i_option_variance_R, INTSXP)); 
	int    *i_option_variance = INTEGER(i_option_variance_R); 		

	//P_additional_1
	i_option_merge_R = PROTECT(Rf_coerceVector(i_option_merge_R, INTSXP)); 
	int    *i_option_merge = INTEGER(i_option_merge_R); 		

	
	//--------
	//prep return variables
	//--------
	rbind_FHDI  rbind_ipmat_FEFI_return(4+ncol_x[0]); //column size is 4+ncol
	rbind_FHDI  rbind_Resp_FEFI_return(ncol_x[0]+1);  //separate response matrix 
	rbind_FHDI  rbind_irmat_FEFI_return(5+ncol_x[0]); //column size is 5+ncol    
	rbind_FHDI  rbind_ipmat_FHDI_return(4+ncol_x[0]); //column size is 4+ncol
	rbind_FHDI  rbind_Resp_FHDI_return(ncol_x[0]+1);  //separate response matrix 
	rbind_FHDI  rbind_irmat_FHDI_return(5+ncol_x[0]); //column size is 5+ncol    
	rbind_FHDI  rbind_vrst_FEFI_return(nrow_x[0]);    //variance estimates of FEFI
	rbind_FHDI  rbind_vrst_FHDI_return(nrow_x[0]);    //variance estimates of FHDI

	//below is for output for Cell Make only option
	rbind_FHDI  rbind_uox_return(ncol_x[0]); //unique observed patterns
	rbind_FHDI  rbind_mox_return(ncol_x[0]); //unique observed patterns
	rbind_FHDI  rbind_category_return(ncol_x[0]); //cagetorized matrix 
	
	//below is for output for Cell Prob only option
	std::vector<std::string> jp_name_return_CellProb;   //name of the joint probability table
	std::vector<double> jp_prob_return_CellProb; //the latest joint probability 	

	//user-defined z matrix (i_option_perform =4) only
	double* z_UserDefined = new double[nrow_x[0]*ncol_x[0]];	

	Rfn_test_call(x, r, nrow_x, ncol_x, k, d, M, 
	              i_option_imputation, i_option_variance, id, z_UserDefined,
				  rbind_ipmat_FEFI_return, rbind_Resp_FEFI_return, rbind_irmat_FEFI_return,
				  rbind_ipmat_FHDI_return, rbind_Resp_FHDI_return, rbind_irmat_FHDI_return,
				  rbind_vrst_FEFI_return,  rbind_vrst_FHDI_return, 
				  
				  rbind_uox_return, rbind_mox_return, rbind_category_return, 
				  jp_name_return_CellProb, jp_prob_return_CellProb, 
				  
				  2, i_option_merge); //2: perform only Cell Make  

	delete[] z_UserDefined;
	
	//testout
	//RPrint("in Rfn_Interface... Rfn_test_call has finished  "); 

	//-----------
	//prep ID, WGT, raw data
	//-----------
	//P11
	SEXP return_raw_y = PROTECT(Rf_allocMatrix(REALSXP, nrow_x[0], ncol_x[0]+2));
	double* d_return_raw_y = REAL(return_raw_y);	
	
	//----
	//copy ID, WGT, x of C++ into the return storage
	//----
	for(int i=0; i< ncol_x[0]+2; i++) //size of columns 
	{
		for(int j=0; j<nrow_x[0]; j++)
		{
			double d_temp = 0.0; 
			if(i==0) //ID column
			{
				d_temp = id[j]; 
			}
			if(i==1) //WGT column
			{
				d_temp = d[j]; 
			}
			if(i>1) //raw data matrix
			{
				d_temp = x[(i-2)*nrow_x[0] + j]; //from x of R
			}
			
			//NOTE: R works column-by-column 
			//hence, below sequence is different from C++ ordering 
			d_return_raw_y[i*nrow_x[0] + j] = d_temp; //get all stored values 
		}
	}		

	//-----------
	//prep Categorized matrix 
	//-----------
	//P12
	SEXP return_category = PROTECT(Rf_allocMatrix(REALSXP, nrow_x[0], ncol_x[0]));
	double* d_return_category = REAL(return_category);	
	
	//----
	//copy category matrix from C++ storage
	//----
	for(int i=0; i< ncol_x[0]; i++) //size of columns 
	{
		for(int j=0; j<nrow_x[0]; j++)
		{
			double d_temp = 0.0; 

			d_temp = rbind_category_return(j,i); //from z of C++
			
			//NOTE: R works column-by-column 
			//hence, below sequence is different from C++ ordering 
			d_return_category[i*nrow_x[0] + j] = d_temp; //get all stored values 
		}
	}		
	
	//----------
	//prep uox & mox
	//----------
	const int i_uox_size_col = rbind_uox_return.size_col(); 
	const int i_uox_size_row = rbind_uox_return.size_row(); 
	const int i_mox_size_col = rbind_mox_return.size_col(); 
	const int i_mox_size_row = rbind_mox_return.size_row(); 
	
	//P13
	SEXP return_uox = PROTECT(Rf_allocMatrix(REALSXP, i_uox_size_row, i_uox_size_col));
	double* d_return_uox = REAL(return_uox);
	//P14
	SEXP return_mox = PROTECT(Rf_allocMatrix(REALSXP, i_mox_size_row, i_mox_size_col));
	double* d_return_mox = REAL(return_mox);	
	
	//----
	//copy uox of C++ into the return storage
	//----
	for(int i=0; i< i_uox_size_col; i++) //size of columns 
	{
		for(int j=0; j<i_uox_size_row; j++)
		{
				double d_temp = 0.0; 
				d_temp = rbind_uox_return(j,i); //from uox of C++

				//NOTE: R works column-by-column 
				//hence, below sequence is different from C++ ordering 
				d_return_uox[i*i_uox_size_row + j] = d_temp; //get all stored values 
			}
	}	
	//----
	//copy mox of C++ into the return storage
	//----
	for(int i=0; i< i_mox_size_col; i++) //size of columns 
	{
		for(int j=0; j<i_mox_size_row; j++)
		{
				double d_temp = 0.0; 
				d_temp = rbind_mox_return(j,i); //from mox of C++

				//NOTE: R works column-by-column 
				//hence, below sequence is different from C++ ordering 
				d_return_mox[i*i_mox_size_row + j] = d_temp; //get all stored values 
			}
	}	
	
	//----
	//create list for compact return
	//-----
	SEXP list_return;
	//P15
	PROTECT(list_return = Rf_allocVector(VECSXP,4));

	//--------
	//Make a final retuan LIST
	//--------
	SET_VECTOR_ELT(list_return, 0, return_raw_y);
	SET_VECTOR_ELT(list_return, 1, return_category);
	SET_VECTOR_ELT(list_return, 2, return_uox);
	SET_VECTOR_ELT(list_return, 3, return_mox);
	//Rf_setAttrib(list_return, R_NamesSymbol, list_names_return);
			
	UNPROTECT(15+1);
	
	return list_return; 	
}
}


//===============================================================
//===============================================================
//===============================================================
//===============================================================
//===============================================================
extern "C" {
SEXP CWrapper_CellProb(SEXP x_R, SEXP nrow_x_R, SEXP ncol_x_R, 
                       SEXP d_R,  
			           SEXP id_R)
//Description -------------------------------------------
//----------------------------------------------------------------------------
//- CWrapper 
//  to perform 
//         Cell Prob only!
// 
// This C function helps invoke the deeper C++ functions. 
// R can access C code, not C++ code directly. 
// This C function provides a C-C++ interface so that R can access the C++
// Note: Rinternals.h is used. In this case, to prevent R's automatic function mapping
//       at the very beginning, #define R_NO_REMAP should be added and also
//       "Rf_" should precede all R functions used herein
// Written by Dr. In Ho Cho
// All rights reserved 
// Feb 9, 2017
//---------------------------------------------------------------------------- 
//IN   : double x_R[nrow_x_R * ncol_x_R] = data matrix in a vector form with missing units
//       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//       NOTE! x_R contains categorized values instead of original data values
//       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//IN   : double d_R[nrow_x_R] = (sampling) "w" weights for units
//IN   : int    id_R[nrow_x_R] = id numbers of all data 
//OUT:  LIST of 
//     [[1]] name of joint probability
//     [[2]] joint probability
//-------------------------------------------------------------------------------
{
	//testout
	//RPrint("Begin CWrapper_CellProb in Rfn_Interface... ==== "); 

	//P1
	x_R = PROTECT(Rf_coerceVector(x_R, REALSXP));
	double *x = REAL(x_R);         	//pointer to double vector x[col*row] that contains all data with missing units
	
	//P2-3
	nrow_x_R = PROTECT(Rf_coerceVector(nrow_x_R, INTSXP));
	ncol_x_R = PROTECT(Rf_coerceVector(ncol_x_R, INTSXP));
	int    *nrow_x = INTEGER(nrow_x_R);
    int	   *ncol_x = INTEGER(ncol_x_R); //pointers to integer sizes of x
	
	//P4
	d_R = PROTECT(Rf_coerceVector(d_R, REALSXP)); 
	double *d = REAL(d_R);  		//pointer to double vector d[nrow_x] 

	//P5
	id_R = PROTECT(Rf_coerceVector(id_R, INTSXP)); 
	int    *id = INTEGER(id_R); 	//pointer to integer number of id of data

	//for this Cell Prob Only option, other variables are nullified 
	int    *r = new int[1]; 			//pointer to an integer vector r[n_total_x] that contains indices of 0 and 1
	double *k = new double[1];			//pointer to a double value 
	int    *M = new int[1]; 		//pointer to integer number of donors 
	int    *i_option_imputation = new int[1]; 		
	int    *i_option_variance  = new int[1]; 		
	int    *i_option_merge      = new int[1]; 
		
	//--------
	//prep return variables
	//--------
	rbind_FHDI  rbind_ipmat_FEFI_return(4+ncol_x[0]); //column size is 4+ncol
	rbind_FHDI  rbind_Resp_FEFI_return(ncol_x[0]+1);  //separate response matrix 
	rbind_FHDI  rbind_irmat_FEFI_return(5+ncol_x[0]); //column size is 5+ncol    
	rbind_FHDI  rbind_ipmat_FHDI_return(4+ncol_x[0]); //column size is 4+ncol
	rbind_FHDI  rbind_Resp_FHDI_return(ncol_x[0]+1);  //separate response matrix 
	rbind_FHDI  rbind_irmat_FHDI_return(5+ncol_x[0]); //column size is 5+ncol    
	rbind_FHDI  rbind_vrst_FEFI_return(nrow_x[0]);    //variance estimates of FEFI
	rbind_FHDI  rbind_vrst_FHDI_return(nrow_x[0]);    //variance estimates of FHDI

	//below is for output for Cell Make only option
	rbind_FHDI  rbind_uox_return(ncol_x[0]); //unique observed patterns
	rbind_FHDI  rbind_mox_return(ncol_x[0]); //unique observed patterns
	rbind_FHDI  rbind_category_return(ncol_x[0]); //cagetorized matrix 
	
	//below is for output for Cell Prob only option
	std::vector<std::string> jp_name_return_CellProb;   //name of the joint probability table
	std::vector<double> jp_prob_return_CellProb; //the latest joint probability 	

	//user-defined z matrix (i_option_perform =4) only
	double* z_UserDefined = new double[nrow_x[0]*ncol_x[0]];	

	
	Rfn_test_call(x, r, nrow_x, ncol_x, k, d, M, 
	              i_option_imputation, i_option_variance, id, z_UserDefined,
				  rbind_ipmat_FEFI_return, rbind_Resp_FEFI_return, rbind_irmat_FEFI_return,
				  rbind_ipmat_FHDI_return, rbind_Resp_FHDI_return, rbind_irmat_FHDI_return,
				  rbind_vrst_FEFI_return,  rbind_vrst_FHDI_return, 
				  
				  rbind_uox_return, rbind_mox_return, rbind_category_return, 
				  jp_name_return_CellProb, jp_prob_return_CellProb, 
				  
				  3, i_option_merge); //3: perform only Cell Prob; x has category values   

	delete[] z_UserDefined;
	
	//-----------
    // prep return list
    // (1) names 
    // (2) joint probability
    //-----------
	const int i_size_jp_Final = (int)jp_name_return_CellProb.size();  
    if(i_size_jp_Final <= 0) 
	{cout<<"Error! zero size of the table of joint probability table"<<endl;}

	//P6
	SEXP jp_name_return_Final;
	PROTECT(jp_name_return_Final = Rf_allocVector(STRSXP, i_size_jp_Final)); 
	for(int i=0; i<i_size_jp_Final; i++) 
	{	
		char* ch_temp = (char*)jp_name_return_CellProb[i].c_str(); 
		SET_STRING_ELT(jp_name_return_Final, i, Rf_mkChar(ch_temp));
	}
	
	//P7
	SEXP jp_prob_return_Final;
	PROTECT(jp_prob_return_Final = Rf_allocVector(REALSXP,i_size_jp_Final)); 
	double* d_jp_prob_return_Final = REAL(jp_prob_return_Final);
	for(int i=0; i<i_size_jp_Final; i++) 
	{	
		d_jp_prob_return_Final[i] = jp_prob_return_CellProb[i]; 
	}
	
	//----
	//create list for compact return
	//-----
	SEXP list_return;
	//P8
	PROTECT(list_return = Rf_allocVector(VECSXP,2));

	
	//--------
	//Make a final retuan LIST
	//--------
	SET_VECTOR_ELT(list_return, 0, jp_name_return_Final);
	SET_VECTOR_ELT(list_return, 1, jp_prob_return_Final);
	//Rf_setAttrib(list_return, R_NamesSymbol, list_names_return);
			
	
	UNPROTECT(8);
	
	return list_return; 	
}
}