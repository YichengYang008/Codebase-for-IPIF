
void Variance_Est_FHDI_Linerization_cpp(MPI_File fh_binary_daty, MPI_File fh_binary_datr, const int nrow, const int ncol, const int memory, double* k,
	List_FHDI &mox_infor, int nrow_mox, double** mox_Imputation_group, int i_size_ml,
	ofstream& TestOut_Slave3) {
	//Description----------------------
	//estimate variance for FHDI using Jackknife method 
	//  Algorithm: 
	//
	// original R code: Dr. Im, J. and Dr. Kim, J. 
	// c++ code: 		Dr. Cho, I. and Yicheng Yang
	// All rights reserved
	// 
	// updated: August 10, 2021
	//
	//----------------------------------------------------
	//IN   : MPI_File fh_binary_daty  = raw daty written column-wisely
	//IN   : MPI_File fh_binary_datr  = raw datr written column-wisely
	//IN   : const int nrow           = number of rows of daty
	//IN   : const int ncol           = number of columns of daty
	//IN   : const int memory         = memory of a core on applied platform
	//IN   : List_FHDI mox_infor      = Actual index list of mox in z
	//IN   : int nrow_mox             = number of rows of mox
	//IN   : double** mox_Imputation_group = (gloabl index of all missing units in z matrix) + (corresponding uox id with the highest conditional probability)
	//IN   : int i_size_ml            = total number of missing units

	//OUT  : TestOut_Slave3
	//-------------------------------------------------------

   	//===========================
	// MPI variables
	//===========================

	//double startup = MPI_Wtime();
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;


	//===========================
	//Read required data
	//===========================

	double Linear1 = MPI_Wtime();

	int success = 0;
	MPI_File fh_datz_reading;//file string to read datz
	MPI_File fh_imputed_daty_reading; //file string to read imputed daty
	MPI_File fh_uox_reading;

	success = MPI_File_open(MPI_COMM_WORLD, "./datz_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_datz_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of z matrix fail to open the file!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./final_daty_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_imputed_daty_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of imputed daty fail to open the file!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./uox_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_uox_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of uox fail to open the file!" << endl;

	const int L = ncol;
	int numWorkPerProc = (int)floor(1.0*L / (1.0*totalnodes - 1));
	int numWorkLocalLast = L - numWorkPerProc * (totalnodes - 2);

	int L_temp = 0;
	if (mynode != (totalnodes - 1) && mynode != 0) L_temp = numWorkPerProc;
	if (mynode == (totalnodes - 1)) L_temp = numWorkLocalLast;

	int startpoint = 0;
	int endpoint = 0;

	if (mynode != 0 && mynode != (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = mynode*numWorkPerProc;

		if (endpoint - startpoint != L_temp) {
			cout << "Linerization boundary ERROR!!!" << endl;
			return;
		}
	}

	if (mynode == (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = (mynode - 1)*numWorkPerProc + numWorkLocalLast;

		if (endpoint - startpoint != L_temp) {
			cout << "Linerization boundary ERROR!!!" << endl;
			return;
		}
	}

	//Memory check: If required memory > 60% available memory ------> Causion!!!
	//Note we have 5 big distributed matrices

	if (mynode != 0) {
		double required_memory = 0.0;
		required_memory = (nrow* L_temp * 8 * 5) / pow(10.0, 9.0);
		if ((0.6*memory) < required_memory) {
			cout << "ERROR!!! The reuiqred memory in Linerization is " << required_memory << ", which is greater than 50% of avaiable memory such that may out of memory!!!!" <<" at node "<<mynode<< endl;
		}
	}


	//---------------------------
	//Read distributed data
	//---------------------------

	double** distr_daty = NULL;
	int** distr_datr = NULL;
	double** distr_datz = NULL;
	double** distr_imputed = NULL;

	if (mynode != 0) {
		double* array_temp_daty = new double[nrow];//daty buffer
		int* array_temp_datr = new int[nrow];//datr buffer
		double* array_temp_datz = new double[nrow];//datz buffer
		double* array_temp_imputed = new double[nrow];//imputed daty buffer

		distr_daty = New_dMatrix(nrow, L_temp);// distributed daty
		distr_datr = New_iMatrix(nrow, L_temp);// distributed datr
		distr_datz = New_dMatrix(nrow, L_temp);// distributed datz
		distr_imputed = New_dMatrix(nrow, L_temp);// distributed daty

		int counter = 0;
		for (int t = startpoint; t < endpoint; t++) {
			double MAMA1 = MPI_Wtime();
			//Note that daty and datr are written column-wisely
			MPI_In_raw(nrow, t, fh_binary_daty, array_temp_daty);
			MPI_In_raw(nrow, t, fh_binary_datr, array_temp_datr);
			//if (t == 1) cout << "YYC Running time of MAMA1 at node " << mynode << " = " << MPI_Wtime() - MAMA1 << endl;

			double MAMA2 = MPI_Wtime();
			//Note that datz and imputed daty are written row-wisely
			MPI_In_datz_column(nrow, ncol, t, fh_datz_reading, array_temp_datz);
			MPI_In_datz_column(nrow, ncol, t, fh_imputed_daty_reading, array_temp_imputed);
			//if (t == 1) cout << "YYC Running time of MAMA2 at node " << mynode << " = " << MPI_Wtime() - MAMA2 << endl;

			for (int m = 0; m < nrow; m++) {
				distr_daty[m][counter] = array_temp_daty[m];
				distr_datr[m][counter] = array_temp_datr[m];
				distr_datz[m][counter] = array_temp_datz[m];
				distr_imputed[m][counter] = array_temp_imputed[m];
			}

			counter++;
		}

		//Deallocation
		delete[] array_temp_daty;
		delete[] array_temp_datr;
		delete[] array_temp_datz;
		delete[] array_temp_imputed;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/*if (mynode == 1) {
		cout << "daty at node " << mynode << endl;
		for (int kk2 = 0; kk2 < nrow; kk2++) {
			for (int kk3 = 0; kk3 < L_temp; kk3++) {
				cout << setw(20) << distr_daty[kk2][kk3];
			}
			cout << endl;
		}

		cout << "datr at node " << mynode << endl;
		for (int kk2 = 0; kk2 < nrow; kk2++) {
			for (int kk3 = 0; kk3 < L_temp; kk3++) {
				cout << setw(20) << distr_datr[kk2][kk3];
			}
			cout << endl;
		}

		cout << "datz at node " << mynode << endl;
		for (int kk2 = 0; kk2 < nrow; kk2++) {
			for (int kk3 = 0; kk3 < L_temp; kk3++) {
				cout << setw(20) << distr_datz[kk2][kk3];
			}
			cout << endl;
		}

		cout << "imputed daty at node " << mynode << endl;
		for (int kk2 = 0; kk2 < nrow; kk2++) {
			for (int kk3 = 0; kk3 < L_temp; kk3++) {
				cout << setw(20) << distr_imputed[kk2][kk3];
			}
			cout << endl;
		}

	}*/


	//if (mynode == 1) {
	//	//cout << "YYC Running time of Cell_make at node " << mynode << " = " << MPI_Wtime() - cell_make_begin << endl;
	//	printf("Yicheng Running time of Linear1 = %f seconds\n", MPI_Wtime() - Linear1);
	//}

	//=======================================================
	//Compute bar{y}_{gl}, r_gl, and n_gl
	//=======================================================

	double Linear2 = MPI_Wtime();

	//if (mynode == 1) {
	//	cout << "mox_Imputation_group: " << endl;
	//	for (int j = 0; j < i_size_ml; j++) {
	//		for (int k = 0; k < 2; k++) {
	//			cout << setw(20) << mox_Imputation_group[j][k];
	//		}
	//		cout << endl;
	//	}
	//}

	//-----------------------------------------
	//Fill in missing cells in distributed datz 
	//with its imputation group
	//-----------------------------------------

	double** distr_datz_full = NULL;

	if (mynode != 0) {

		distr_datz_full = New_dMatrix(nrow, L_temp);// distributed datr
		double* d_read_out = new double[ncol];

		std::vector<int> ml_location; // Actual location of missing patterns in z matrix (not sorted)
		std::vector<double> ml_location_temp; // Actual location of missing patterns in z matrix (not sorted)
		mox_infor.unlist(ml_location_temp);

		int i_v_ol = 0;
		for (int t = 0; t < ml_location_temp.size(); t++) {
			i_v_ol = 0;
			i_v_ol = (int)ml_location_temp[t];
			ml_location.push_back(i_v_ol);
		}


		if (i_size_ml != ml_location.size()) cout<<"ERROR!!!! i_size_ml != ml_location.size() in variance using Linerization!!!!"<<endl;

		//Sort actual locations of missing units and fully observed units
		sort(ml_location.begin(), ml_location.end());


		//Copy distributed datz
		for (int l = 0; l < nrow; l++) {
			for (int j = 0; j < L_temp;j++) {
				distr_datz_full[l][j] = distr_datz[l][j];
			}
		}


		int row_temp = 0;
		int uox_id = 0;

		for (int l = 0; l < i_size_ml; l++) {
			row_temp = 0;
			row_temp = ml_location[l];// Actual location
			uox_id = 0;
			for (int t = 0; t < i_size_ml; t++) {
				if (fabs(row_temp - mox_Imputation_group[t][0]) <1e-15) {
					uox_id = mox_Imputation_group[t][1];// Actual location
					break;
				}
			}

			if (uox_id == 0) cout<<"ERROR!!!! Did not find matched uox in mox_Imputation_group for missing unit "<<l<<endl;

			MPI_In_uox_mox(ncol, uox_id - 1, fh_uox_reading, d_read_out);// read in uox row


			for (int j = 0; j < L_temp;j++) {
				if (distr_datz_full[row_temp - 1][j] < 1e-15) {
					distr_datz_full[row_temp - 1][j] = d_read_out[startpoint + j];
				}
			}

		}


		//if (mynode == 1) {
		//	cout << "distr_datz_full at node " << mynode << endl;
		//	for (int kk2 = 0; kk2 < nrow; kk2++) {
		//		for (int kk3 = 0; kk3 < L_temp; kk3++) {
		//			cout << setw(20) << distr_datz_full[kk2][kk3];
		//		}
		//		cout << endl;
		//	}
		//}

		delete[] d_read_out;

	}// end of mynode



	//Get maximum number of group in k
	//In case some variables have different number of categories

	int kmax = 0;
    
	for (int j = 0; j < ncol; j++) {
		if (kmax < k[j]) {
			kmax = (int)k[j];
		}
	}

    //--------------------------------------
	//Compute bar{y}_{gl}, r_gl, and n_gl
	//--------------------------------------

	double** y_gl = NULL;
	int** r_gl = NULL;
	int** n_gl = NULL;

	int r_gl_temp = 0;
	int n_gl_temp = 0;
	int k_temp = 0;

	if (mynode != 0) {
		y_gl = New_dMatrix(kmax, L_temp);
		r_gl = New_iMatrix(kmax, L_temp);
		n_gl = New_iMatrix(kmax, L_temp);

		//-------------------------
		//Compute bar{y}_{gl}, r_gl
		//-------------------------

		double y_gl_sum = 0.0;
		
		for (int j = 0; j < L_temp; j++) {
			k_temp = 0;
			k_temp = (int)k[startpoint + j];

			for (int g = 0; g < k_temp; g++) {
				y_gl_sum = 0.0;
				r_gl_temp = 0;

				for (int i = 0; i < nrow; i++) {

					if ( fabs(distr_datz[i][j] - (g + 1)) < 1e-15 ) {
						r_gl_temp++;
						y_gl_sum = y_gl_sum + distr_daty[i][j];
					}

				}

				r_gl[g][j] = r_gl_temp;
				y_gl[g][j] = y_gl_sum / r_gl_temp;//Caution!! r_gl_temp can be 0 such that y_gl = -nan

			}

		}//end of main loop

		//---------------
		// Compute n_gl
		//----------------

		for (int j = 0; j < L_temp; j++) {
			k_temp = 0;
			k_temp = (int)k[startpoint + j];

			for (int g = 0; g < k_temp; g++) {
				n_gl_temp = 0;

				for (int i = 0; i < nrow; i++) {

					if ( fabs(distr_datz_full[i][j] - (g + 1)) < 1e-15 ) {
						n_gl_temp++;
					}

				}

				n_gl[g][j] = n_gl_temp;

			}

		}//end of main loop

		int n_gl_sum = 0;
		for (int m = 0; m < L_temp; m++) {

			n_gl_sum = 0;

			for (int f = 0; f < kmax; f++) {
				n_gl_sum = n_gl_sum + n_gl[f][m];
			}

			if (n_gl_sum != nrow) cout<<"ERROR!!! n_gl != nrow in variance linerization !!!!"<<endl;
		}

		//Error check for linearized variance estimation made on July 15, 2021
		//for (int b = 0; b < L_temp; b++) {
		//	for (int t = 0; t < kmax; t++) {
		//		if (r_gl[t][b] == 0) {
		//			cout << "CAUTION!!! The number of imputation cells (r_gl) in category " << t+1 <<" at the variable "<< startpoint + b << " is zero!!!" << endl;
		//		}
		//	}
		//}


		//if (mynode == 1) {
		//	cout << "r_gl at node " << mynode << endl;
		//	for (int kk2 = 0; kk2 < kmax; kk2++) {
		//		for (int kk3 = 0; kk3 < L_temp; kk3++) {
		//			cout << setw(20) << r_gl[kk2][kk3];
		//		}
		//		cout << endl;
		//	}

		//	cout << "y_gl at node " << mynode << endl;
		//	for (int kk2 = 0; kk2 < kmax; kk2++) {
		//		for (int kk3 = 0; kk3 < L_temp; kk3++) {
		//			cout << setw(20) << y_gl[kk2][kk3];
		//		}
		//		cout << endl;
		//	}

		//	cout << "n_gl at node " << mynode << endl;
		//	for (int kk2 = 0; kk2 < kmax; kk2++) {
		//		for (int kk3 = 0; kk3 < L_temp; kk3++) {
		//			cout << setw(20) << n_gl[kk2][kk3];
		//		}
		//		cout << endl;
		//	}
		//}

	}//end of mynode

	//if (mynode == 1) {
	//	//cout << "YYC Running time of Cell_make at node " << mynode << " = " << MPI_Wtime() - cell_make_begin << endl;
	//	printf("Yicheng Running time of Linear2 = %f seconds\n", MPI_Wtime() - Linear2);
	//}

	 //=======================================
	 //Compute variance
	 //=======================================

	double Linear3 = MPI_Wtime();

	if (mynode != 0) {
		double** eta_hat = New_dMatrix(nrow, L_temp);// distributed daty
		double* eta_bar = new double[L_temp];//daty buffer

		//Compute eta_hat
		int group_temp = 0;
		for (int j = 0; j < L_temp; j++) {
			for (int i = 0; i < nrow; i++) {
				if (distr_datr[i][j] == 0) {
					eta_hat[i][j] = distr_imputed[i][j];
				}

				if (distr_datr[i][j] != 0) {
					group_temp = 0;
					group_temp = (int)distr_datz_full[i][j] - 1;
					if (r_gl[group_temp][j] == 0) cout << "ERROR! The denominator (r_gl) in equation to compute linearized variance is zero such that variance can be -nan !!" << endl;
					eta_hat[i][j] = y_gl[group_temp][j] + (double)n_gl[group_temp][j] * (distr_daty[i][j] - y_gl[group_temp][j]) / r_gl[group_temp][j];
				}

			}
		}

		//if (mynode == 1) {
		//	cout << "eta_hat at node " << mynode << endl;
		//	for (int kk2 = 0; kk2 < nrow; kk2++) {
		//		for (int kk3 = 0; kk3 < L_temp; kk3++) {
		//			cout << setw(20) << eta_hat[kk2][kk3];
		//		}
		//		cout << endl;
		//	}
		//}

		//Compute eta_bar
		double sum_temp = 0.0;
		for (int j = 0; j < L_temp; j++) {
			sum_temp = 0.0;
			for (int t = 0; t < nrow; t++) {
				sum_temp = sum_temp + eta_hat[t][j];
			}
			eta_bar[j] = sum_temp / nrow;
		}

		//if (mynode == 1) {
		//	for (int j = 0;j < L_temp; j++) {
		//		cout<<"eta_bar["<<j<<"]: "<< eta_bar[j]<<endl;
		//	}
		//}

		//Compute variance
		std::vector<double> variance_temp;

		double v_temp = 0.0;
		//double para =  1.0 / (nrow*(nrow - 1));//This is incorrect way to compute 
		double mag = (double)nrow*(nrow - 1);//Important change on July 15, 2021!!!
		double para = 1.0 / mag;
		if (para < 1e-15) { cout << "ERROR! Computation of coefficient in linearized variance estimation is incorrect!" << endl; }

		//cout<<"para is "<< para <<" at node "<<mynode<<endl;

		for (int j = 0; j < L_temp; j++) {
			sum_temp = 0.0;

			for (int i = 0; i < nrow; i++) {
				sum_temp = sum_temp + (eta_hat[i][j] - eta_bar[j])*(eta_hat[i][j] - eta_bar[j]);
			}

			v_temp = 0.0;
			v_temp = para * sum_temp;
			variance_temp.push_back(v_temp);
		}

		MPI_Send(&variance_temp[0], L_temp, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

		//if (mynode == 1) {
		//	for (unsigned int t = 0; t < variance_temp.size(); t++) {
		//		cout<<"variance_temp["<<t<<"]: "<< variance_temp[t]<<endl;
		//	}
		//}

		//Deallocation
		Del_dMatrix(eta_hat, nrow, L_temp);
		delete[] eta_bar;

	}//end of mynode

	//-------------------------------------
	//ALL_Gather variance on master node
	//-------------------------------------

	if (mynode == 0) {

		std::vector<double> final_variance_data;

		for (int j = 1; j < totalnodes; j++) {

			if (j != (totalnodes - 1)) {
				std::vector<double> variance_data_recv(numWorkPerProc);
				MPI_Recv(&variance_data_recv[0], numWorkPerProc, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);

				for (int l = 0; l < numWorkPerProc; l++) {
					final_variance_data.push_back(variance_data_recv[l]);
				}
			}

			if (j == (totalnodes - 1)) {
				std::vector<double> variance_data_recv(numWorkLocalLast);
				MPI_Recv(&variance_data_recv[0], numWorkLocalLast, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);

				for (int l = 0; l < numWorkLocalLast; l++) {
					final_variance_data.push_back(variance_data_recv[l]);
				}
			}


		}

		int variance_size = (int)final_variance_data.size();
		if (variance_size != ncol) cout<<"ERROR!!!! variance_size != ncol in variance using linerization!!! "<<endl;

		TestOut_Slave3 << "Linerization Variance Results:" << endl;
		for (int i = 0; i < ncol; i++) {
			TestOut_Slave3 << final_variance_data[i] << endl;
		}

	}//end of mynode

	//if (mynode == 0 || mynode == 1) {
	//	cout << "YYC Running time of Linear3 at node " << mynode << " = " << MPI_Wtime() - Linear3 << endl;
	//}


	//====================
	//Close files 
	//=====================

	success = MPI_File_close(&fh_datz_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of z matrix fail to close the file!" << endl;

	success = MPI_File_close(&fh_imputed_daty_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of imputed daty fail to close the file!" << endl;

	success = MPI_File_close(&fh_uox_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of uox fail to close the file!" << endl;


	//==========================
	//Deallocation
	//===========================

	if (mynode != 0) {

		Del_dMatrix(distr_daty, nrow, L_temp);
		Del_iMatrix(distr_datr, nrow, L_temp);
		Del_dMatrix(distr_datz, nrow, L_temp);
		Del_dMatrix(distr_imputed, nrow, L_temp);

		Del_dMatrix(distr_datz_full, nrow, L_temp);

		Del_dMatrix(y_gl, kmax, L_temp);
		Del_iMatrix(r_gl, kmax, L_temp);
		Del_iMatrix(n_gl, kmax, L_temp);


	}



	return;
}

