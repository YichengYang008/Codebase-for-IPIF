#include "tune_points_var.cc"

void Extract_Imputed_Results_Ultra(const int nrow, const int ncol, int i_nrow_imputation, 
	            
	double** simp_fmat_FHDI, MPI_File fh_final_binary_daty,

	ofstream& TestOut_Slave3)
                             
//Description========================================
//  extract imputed values from the fh_fmat_FHDI 
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
//  Written by Yicheng Yang
//  updated Jan, 06, 2021. 
//
//IN   : const int nrow           = number of rows of daty
//IN   : const int ncol           = number of columns of daty
//IN   : int i_nrow_imputation    = number of rows of imputed values
//IN   : double** simp_fmat_FHDI  = matrix of 4 columns: global id + sampling weight + wij + fwij
//IN   : MPI_File fh_fmat_FHDI    = file string to read imputed matrix in size of (4+ncol) columns: global id + sampling weight + wij + fwij + imputed daty


//OUT  : MPI_File fh_final_binary_daty = file string to write final matrix in size of ncol columns
//===================================================  
{
	//-- MPI variables
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//open file string to read imputed values
	int success = 0;
	MPI_File fh_fmat_FHDI_reading; // File to write imputed values

	success = MPI_File_open(MPI_COMM_WORLD, "./fmat_FHDI_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_fmat_FHDI_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of mox fail to open the file!" << endl;


	//--------------------------------------------
	//Sort simp_fmat_FHDI in ascending order
	//--------------------------------------------

	int* rbind_ipmat = new int[i_nrow_imputation]; // copy of global id of imputed values 
	for (int kk = 0; kk < i_nrow_imputation; kk++) {
		rbind_ipmat[kk] = (int)simp_fmat_FHDI[kk][0];//1st column
	}

	int* order_rbind_ipmat = new int[i_nrow_imputation];
	for (int kk2 = 0; kk2 < i_nrow_imputation; kk2++) {
		order_rbind_ipmat[kk2] = (int)simp_fmat_FHDI[kk2][0];//1st column
	}

	order_FHDI_binary(order_rbind_ipmat, i_nrow_imputation); //mapping from original global id to sorted global id (actual locations)

	//if (mynode == 0) {
	//	for (int k = 0; k < i_nrow_imputation; k++) {
	//		cout<<"order_rbind_ipmat["<<k<<"]: "<< order_rbind_ipmat[k]<<endl;
	//	}
	//}

	double** simp_fmat_FHDI_sorted = New_dMatrix(i_nrow_imputation, 4);// reorder of simp_fmat_FHDI in ascending order
	for (int t = 0; t < i_nrow_imputation; t++) {
		for (int j = 0; j < 4; j++) {
			simp_fmat_FHDI_sorted[t][j] = simp_fmat_FHDI[order_rbind_ipmat[t]-1][j];//actual -1
		}
	}

	//if (mynode == 1) {
	//		cout << "simp_fmat_FHDI_sorted at node " << mynode<< endl;
	//		for (int p = 0; p < i_nrow_imputation; p++) {
	//			for (int k = 0; k < 4; k++) {
	//				cout << setw(20) << simp_fmat_FHDI_sorted[p][k];
	//			}
	//			cout << endl;
	//		}
	//}

	//---------------------------
	//Job assignment
	//---------------------------
	const int L = i_nrow_imputation; //size of d_rw 
	int numWorkPerProc = (int)floor(1.0*L / (1.0*totalnodes - 1));
	int numWorkLocalLast = L - numWorkPerProc * (totalnodes - 2);

	//===============================
	// specify startpoint and end point of imputed values for slave processors
	//===============================
	int startpoint = 0;
	int endpoint = 0;

	if (mynode != 0 && mynode != (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = mynode*numWorkPerProc;
	}

	if (mynode == (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = (mynode - 1)*numWorkPerProc + numWorkLocalLast;
	}

	// adjust boundary indices accordingly to avoid conflict
	int L_temp = 0;// adjusted number of rows of distributed imputed values on each slave processor

	if (mynode != 0) {
		//cout << "Variance startpoint is " << startpoint << " and endpoint is " << endpoint << " at node " << mynode << endl;

		tune_points_var(startpoint, endpoint, rbind_ipmat, order_rbind_ipmat, mynode, totalnodes, i_nrow_imputation);

		//cout << "Variance After tune startpoint is " << startpoint << " and endpoint is " << endpoint << " at node " << mynode << endl;

		L_temp = endpoint - startpoint;
		//cout << "Variance L_temp is " << L_temp << " at node " << mynode << endl;
	}

	//-------------------------------
	//Initialization
	//--------------------------------
	int ncol_fmat = ncol + 4; // 4 columns + ncol : global id + sampling weight + wij + fwij + imputed values
	double* d_1_imputed = new double[ncol_fmat]; //buffer to read imputed values
	
	double** imputed_values = NULL;
	double** final_full_data = NULL;
	int startpoint_final = 0;
	int endpoint_final = 0;
	int nrow_temp = 0;

	//------------------------------------------------
	//Find corresponding boundary indices in daty 
	//------------------------------------------------

	if (mynode != 0) {
		imputed_values = New_dMatrix(L_temp, ncol);//distributed imputed values

		startpoint_final = rbind_ipmat[order_rbind_ipmat[startpoint] - 1] -1; //actual loc
		endpoint_final = rbind_ipmat[order_rbind_ipmat[endpoint-1] - 1] ;

		nrow_temp = endpoint_final - startpoint_final;
		//int nrow_temp = rbind_ipmat[order_rbind_ipmat[endpoint - 1] - 1] - rbind_ipmat[order_rbind_ipmat[startpoint] - 1] + 1;

		//cout << "nrow_temp is " << nrow_temp << "; startpoint_final is " << startpoint_final << "; endpoint_final is " << endpoint_final << " at node " << mynode << endl;

		final_full_data = New_dMatrix(nrow_temp, ncol);//distributed final daty
	}

	//---------------------------------
	//Broadcast number of dustributed daty to all processors
	//-------------------------------
	// To write final distributed final daty to HD
	// one has to know number of rows distributed on each slave processor in advance
	if (mynode != 0) {
		MPI_Send(&nrow_temp, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
	}

	std::vector<int> nrow_size_total;//number of rows of distributed final daty on all slave processor
	if (mynode == 0) {
		int nrow_size_total_recv = 0;
		nrow_size_total.push_back(nrow_size_total_recv);//initialize with 0

		for (int j = 1; j < totalnodes; j = j + 1) {
			MPI_Recv(&nrow_size_total_recv, 1, MPI_INT, j, 2, MPI_COMM_WORLD, &status);
			nrow_size_total.push_back(nrow_size_total_recv);
		}
	}

	if (mynode != 0) {
		nrow_size_total.resize(totalnodes);
	}

	MPI_Bcast(&nrow_size_total[0], totalnodes, MPI_INT, 0, MPI_COMM_WORLD);

	int i_temp = 0;
	for (int t = 0; t < totalnodes; t++) {
		i_temp = i_temp + nrow_size_total[t];
	}

	if (i_temp != nrow) cout<<"ERROR!!!!! The number of distributed final daty is wrong in Extract_Imputed_Results_Ultra!!!!"<<endl;
	//if (mynode == 1) {
	//	for (int t = 0; t < nrow_size_total.size(); t++) {
	//		cout<<"nrow_size_total["<<t<<"]: "<< nrow_size_total[t]<<endl;
	//	}
	//}

	//----------------------------------
	//Read distributed imputed values
	//----------------------------------

	int counter = 0;
	for (int i = startpoint; i < endpoint; i++) {
		int i_temp = order_rbind_ipmat[i] -1;

		MPI_In_uox_mox(ncol_fmat, i_temp, fh_fmat_FHDI_reading, d_1_imputed);//only for validation purpos

		for (int j = 0; j < ncol; j++) {
			imputed_values[counter][j] = d_1_imputed[j + 4];
		}
		counter++;
	}

	//if (mynode ==3 ) {
	//	cout << "imputed_values at node " << mynode<< endl;
	//	for (int p = 0; p < L_temp; p++) {
	//		for (int k = 0; k < ncol; k++) {
	//			cout << setw(20) << imputed_values[p][k];
	//		}
	//		cout << endl;
	//	}
	//}

	//--------
	//main loop for all rows of original data matrix
	//--------
	double* yi = new double[ncol];
	int i_loc = startpoint;
	int i_loc2 = 0;
	int counter2 = 0;

	for (int i = startpoint_final; i < endpoint_final; i++) {
		//-----
		//inner loop within the identical ID 
		//-----
		double d_sum_wij = 0.0;
		Fill_dVector(yi, ncol, 0.0); //initialize ith row data vector 

		for (int j = 0; j < L_temp; j++) {
			//if(mynode==1) cout<<"i_loc is "<< i_loc <<" and i_loc2 is  "<< i_loc2 <<" at i = "<<i<<endl;

			int ID = simp_fmat_FHDI_sorted[i_loc][0] - 1; //1st col means ID; "-1" for actual location 
			//if (mynode == 3) cout << "ID is " << ID <<" at i = "<<i<< endl;
			if (ID == i) {
				double wi = simp_fmat_FHDI_sorted[i_loc][1];
				double wij = simp_fmat_FHDI_sorted[i_loc][2];

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
					yi[i_var] = yi[i_var] + wi*wij * imputed_values[i_loc2][i_var];
				}	
		
				//----
				//increment location for next row 
				//----
				//if(mynode==1) cout<<"i_loc is "<< i_loc <<" and i_loc2 is  "<< i_loc2 <<" at i = "<<i<<endl;
				i_loc++; 
				i_loc2++;
			}

			if ( (ID > i) || (i_loc == endpoint) ) { break; } //exit inner loop
		}
		
		//if (mynode == 1) cout<<"d_sum_wij is "<< d_sum_wij <<" at i = "<<i<<endl;

		if (fabs(d_sum_wij - 1.0) > 1e-15)
		{
			cout << "ERROR! sum of fractional weight at the row: " << i << " is not one, maybe sampling weights are not all 1s !!!!" <<endl;
			//return;
		}

		for (int i_var = 0; i_var< ncol; i_var++) //size of columns of ipmat matrix of C++
		{
			double d_temp = yi[i_var] / d_sum_wij;

			//NOTE: R works column-by-column 
			//hence, below sequence is different from C++ ordering 
			final_full_data[counter2][i_var] = d_temp;  //note: i=current row
		}

		counter2++;

	}


	//======================================
	//Extract mean values of final daty
	//=======================================

	if (mynode != 0) {
		std::vector<double> final_sum_temp(ncol); // column sum on each slave processor
		for (int k = 0; k < ncol; k++) {
			for (int i = 0; i < nrow_temp; i++) {
				final_sum_temp[k] = final_sum_temp[k] + final_full_data[i][k];
			}
		}

		MPI_Send(&final_sum_temp[0], ncol, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}


	if (mynode == 0) {
		std::vector<double> final_sum_recv(ncol);//recv buffer
		std::vector<double> final_sum(ncol);//accumulate column sum over all slave processors
		std::vector<double> colMean(ncol);//final column mean


		for (int j = 1; j < totalnodes; j++) {
			MPI_Recv(&final_sum_recv[0], ncol, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);

			for (int k = 0; k < ncol; k++) {
				final_sum[k] = final_sum[k] + final_sum_recv[k];
			}
		}

		for (int t = 0; t < ncol; t++) colMean[t] = final_sum[t] / nrow;

		//-- Testout column mean
		TestOut_Slave3 << "Mean estimates:" << endl;
		for (int i = 0; i < ncol; i++) {
			//TestOut << "bbfore" << endl;
			//TestOut_Slave3 << setw(20) << colMean[i];
			TestOut_Slave3 << colMean[i] << endl;
			//TestOut << "aafter" << endl;
		}
		//TestOut_Slave3 << endl;
	}


	//Write distributed final daty matrix to local hard drive concurrently
	MPI_Out_uox_mox(nrow_temp, ncol, nrow_size_total, fh_final_binary_daty, final_full_data);
	MPI_Barrier(MPI_COMM_WORLD);

	//close file string to read imputed values
	success = MPI_File_close(&fh_fmat_FHDI_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of imputed matrix fail to close the file!" << endl;


	//if (mynode ==3 ) {
	//		cout << "final_full_data at node " << mynode<< endl;
	//		for (int p = 0; p < nrow_temp; p++) {
	//			for (int k = 0; k < ncol; k++) {
	//				cout << setw(20) << final_full_data[p][k];
	//			}
	//			cout << endl;
	//		}
	//}

	//Deallocation
	delete[] order_rbind_ipmat;
	delete[] rbind_ipmat;
	delete[] yi;
	delete[] d_1_imputed;

	Del_dMatrix(simp_fmat_FHDI_sorted, i_nrow_imputation, 4);

	if (mynode != 0) {
		Del_dMatrix(imputed_values, L_temp, ncol);
		Del_dMatrix(final_full_data, nrow_temp, ncol);
	}

	return; 
}
