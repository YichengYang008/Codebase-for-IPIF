
bool ErrorCheck(const int i_option_read_data, const int i_option_ultra, const int i_option_perform, const int i_option_imputation, const int i_option_variance,
	const int i_option_merge, const int nrow_x, const int ncol_x, const int M, const int i_user_defined_datz, const int i_option_collapsing, 
	const int i_option_SIS_type, const int top, const int i_option_cellmake, const int i_option_var_type, const int memory,
	double* k, int* NonCollapsible_categorical, double* d, ofstream& TestOut) {
	//Description=========================================

	// Basic Error Check

	// original R code: Dr. Im, J. and Dr. Kim, J. 

	// c++ code: 		Yicheng Yang 

	// All rights reserved

	// 

	// updated: July 5, 2021

	//=====================================================

	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//---------------------
	//Invalid option check
	//---------------------

	//(1) i_option_read_data 
	if (i_option_read_data != 0 && i_option_read_data != 1) {
		TestOut << "Error! i_option_read_data must be either 0 or 1 for different IO methods" << endl;
		return 0;
	}

	//(2) i_option_ultra
	if (i_option_ultra != 0 && i_option_ultra != 1) {
		TestOut << "Error! i_option_ultra must be either 0 for P-FHDI or 1 for UP-FHDI" << endl;
		return 0;
	}

	//(3) i_option_perform
	if (i_option_perform != 1 && i_option_perform != 2 && i_option_perform != 3 && i_option_perform != 4) {
		TestOut << "ERROR! i_option_perform must be an integer in [1,4] to perform different stages of FHDI based on users’ purpose" << endl;
		return 0;
	}

	//(4) i_option_imputation
	if (i_option_imputation != 1 && i_option_imputation != 2) {
		TestOut << "ERROR! i_option_imputation must be either 1 for FEFI or 2 for FHDI" << endl;
		return 0;
	}

	//(5) i_option_variance
	if (i_option_variance != 0 && i_option_variance != 1) {
		TestOut << "ERROR! i_option_variance must be either 0 to skip variance estimation or 1 to perform variance estimation" << endl;
		return 0;
	}

	//(6) i_option_merge
	if (i_option_merge != 0 && i_option_merge != 1) {
		TestOut << "ERROR! i_option_merge must be either 0 to fix seed or 1 to enable random number generator" << endl;
		return 0;
	}

	//(7) i_donor
	if (M < 1 || M > nrow_x) {
		TestOut << "ERROR! i_donor must be an integer in [1,nrow] where nrow = " << nrow_x << endl;
		return 0;
	}

	//(8) i_user_defined_datz
	if (i_user_defined_datz != 0 && i_user_defined_datz != 1) {
		TestOut << "ERROR! i_user_defined_datz must be either 0 for automatic categorization or 1 for user-defined categorized matrix" << endl;
		return 0;
	}

	//(9) i_option_collapsing
	if (i_option_collapsing < 0 || i_option_collapsing > ncol_x) {
		TestOut << "ERROR! i_option_collapsing must be 0 to disable variable reduction or an integer in [1,ncol] to enable sure independent screening" << endl;
		return 0;
	}

	//(10) i_option_SIS_type
	if (i_option_collapsing > 0 && i_option_collapsing < (ncol_x + 1)) {
		if (i_option_SIS_type != 1 && i_option_SIS_type != 2 && i_option_SIS_type != 3) {
			TestOut << "ERROR! i_option_SIS_type must be 1 for intersection or 2 for union or 3 for the global ranking of simple correlations if i_option_collapsing is an integer in [1,ncol]" << endl;
			return 0;
		}
	}

	//(11) top
	if (top < 1) {
		TestOut << "ERROR! top_correlation must be an integer greater than 0" << endl;
		return 0;
	}

	//(12) i_option_cellmake
	if (i_option_cellmake != 1 && i_option_cellmake != 2) {
		TestOut << "ERROR! i_option_cellmake must be either 1 for cell construction with cell collapsing or 2 for cell construction with k-nearest neighbors" << endl;
		return 0;
	}

	//(13) i_option_var_type
	if (i_option_var_type != 1 && i_option_var_type != 2) {
		TestOut << "ERROR! i_option_var_type must be either 1 for Jackknife variance estimation or 2 for linearized variance estimation" << endl;
		return 0;
	}

	//(14) memory
	if (i_option_ultra == 1) {
		if (memory < 1) {
			TestOut << "ERROR! memory must be an integer greater than 0 for UP-FHDI" << endl;
			return 0;
		}
	}


	//(15) max_overlap_size
	//if (i_option_ultra == 1) {
	//	if (max_overlap_size < 1) {
	//		TestOut << "ERROR! max_pattern must be an integer greater than 0 for UP-FHDI" << endl;
	//		return 0;
	//	}
	//}

	//if (i_option_ultra == 1) {
	//	if (max_overlap_size > nrow_x) {
	//		TestOut << "CAUTION! Since max_overlap_size is greater than nrow_x for UP-FHDI, we suggest reducing max_overlap_size to nrow_x to avoid wasting memory" << endl;
	//	}
	//}

	//(16) max_donor_size
	//if (i_option_ultra == 1) {
	//	if (max_donor_size < 1) {
	//		TestOut << "ERROR! max_donor must be an integer greater than 0 for UP-FHDI" << endl;
	//		return 0;
	//	}
	//}

	//if (i_option_ultra == 1) {
	//	if (max_donor_size > nrow_x) {
	//		TestOut << "CAUTION! Since max_donor_size is greater than nrow_x for UP-FHDI, we suggest reducing max_donor_size to nrow_x to avoid wasting memory" << endl;
	//	}
	//}

	//(17) category k
	for (int i = 0; i < ncol_x; i++) {
		if (k[i] < 1 || k[i] > 35) {
			TestOut << "ERROR! Category must be an integer in [1,35], but k[" << i << "] is " << k[i] << endl;
			return 0;
		}
	}

	//(18) NonCollapsible_categorical
	for (int i = 0; i < ncol_x; i++) {
		if (NonCollapsible_categorical[i] != 0 && NonCollapsible_categorical[i] != 1) {
			TestOut << "ERROR! NonCollapsible_categorical must be either 0 for a continuous variable or 1 for a non-collapsible categorical variable, but NonCollapsible_categorical["<<i<<"] is "<< NonCollapsible_categorical[i]<< endl;
			return 0;
		}
	}

	//(19) sampling weight d
	for (int i = 0; i < nrow_x; i++) {
		if (d[i] <= 0) {
			TestOut << "ERROR! Sampling weight must be greater than 0, but d["<<i<<"] is "<<d[i]<< endl;
			return 0;
		}
	}


	//(20) i_user_defined_datz
	if (i_option_perform == 3) //CellProb only
	{
		if (i_user_defined_datz != 1)
		{
			TestOut << "ERROR! If i_option_perform = 3 (perform CellProb only), i_user_defined_datz must be 1 such that users should provide a user-defined categorized matrix" << endl;
			return 0;
		}
	}

	//(21) i_user_defined_datz
	if (i_option_perform == 4) //do all using a defined datz
	{
		if (i_user_defined_datz != 1)
		{
			TestOut << "ERROR! If i_option_perform = 4 (perform all stages using a user-defined categorized matrix), i_user_defined_datz must be 1 such that users should provide a user-defined categorized matrix" << endl;
			return 0;
		}
	}


	//----------------------------------------
	//Difference between P-FHDI and UP-FHDI
	//----------------------------------------

	//(22) i_option_perform
	if (i_option_ultra == 1) {
		if (i_option_perform != 1) {
			TestOut << "ERROR! UP-FHDI only allows to perform all stages of UP-FHDI such that i_option_perform must be 1" << endl;
			return 0;
		}
	}

	//(23) i_option_imputation
	if (i_option_ultra == 1) {
		if (i_option_imputation != 2) {
			TestOut << "ERROR! UP-FHDI only allows FHDI such that i_option_imputation must be 2" << endl;
			return 0;
		}
	}

	//(24) i_user_defined_datz
	if (i_option_ultra == 1) {
		if (i_user_defined_datz != 0) {
			TestOut << "ERROR! UP-FHDI only allows automatic categorization such that i_user_defined_datz must be 0" << endl;
			return 0;
		}
	}

	//(25) i_option_cellmake
	if (i_option_ultra == 1) {
		if (i_option_cellmake != 2) {
			TestOut << "ERROR! UP-FHDI only allows cell construction with k-nearest neighbors such that i_option_cellmake must be 2" << endl;
			return 0;
		}
	}

	//(26) i_option_var_type
	if (i_option_ultra == 0) {
		if (i_option_var_type != 1) {
			TestOut << "ERROR! P-FHDI only allows the Jackknife variance estimation such that i_option_var_type must be 1" << endl;
			return 0;
		}
	}

	//(27) i_option_ultra
	if (i_option_ultra == 1) {
		if ( (totalnodes - 1) > ncol_x ) {
			TestOut << "ERROR! When total number of slave processors is greater than ncol, UP-FHDI is not applicable for this dataset" << endl;
			return 0;
		}
	}

	//(28) i_option_ultra
	if (ncol_x > 50 && i_option_ultra == 0) {
		TestOut << "CAUTION! Empirically, we suggest UP-FHDI when ncol is greater than 50 " << endl;
	}


	//---------------------
	//Out of memory check
	//---------------------

	//(29) input data volume for P-FHDI
	if (i_option_ultra == 0) {
		double input_volume = 0.0;
		input_volume = (double)(nrow_x*ncol_x * 8 * 6) / pow(10.0, 9.0);//(x_temp + r_temp + z_temp + x + r + z) on all processors
		if (input_volume > memory) {
			TestOut << "ERROR! When P-FHDI is activated, the maximum input data volume is strictly restricted by the available memory of a processor. Please request fewer tasks per node to give each processor more memory" << endl;
			return 0;
		}
	}

	//(30) input data volume for UP-FHDI POSIX IO 
	if (i_option_ultra == 1 && i_option_read_data == 0) {
		double input_volume = 0.0;
		input_volume = (double)(nrow_x*ncol_x * 8 * 2) / pow(10.0, 9.0);//(x + r) on the master processors
		if (input_volume > memory) {
			TestOut << "ERROR! When UP-FHDI is activated, the maximum input data volume of POSIX IO is strictly restricted by the available memory of the master processor. Please request fewer tasks per node to give the master processor more memory" << endl;
			return 0;
		}
	}

	//(31) List_nU
	//if (i_option_ultra == 1 && ncol_x > 20) {
	//	double List_nU_volume = 0.0;
	//	List_nU_volume = (double)(0.7*nrow_x*max_donor_size * 8) / pow(10.0, 9.0);//List_nU on all processors. Note assume 80% recipients
	//	if (List_nU_volume > memory) {
	//		TestOut << "CAUTION! When UP-FHDI is activated, the size of donor list for all recipients are restricted by the available memory of a processor. Please decrease the value of max_donor_size" << endl;
	//	}
	//}

	return 1;
}