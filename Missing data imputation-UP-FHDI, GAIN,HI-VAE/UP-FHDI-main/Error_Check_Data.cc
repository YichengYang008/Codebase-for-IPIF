
bool ErrorCheck_Data(const int nrow, const int ncol, const int L_temp, int startpoint, int numWorkPerProc, int numWorkLocalLast, double** x_temp, int** r_temp, ofstream& TestOut) {
	//Description=========================================

	// Basic Error Check of raw data

	// Check case 1: If there exists identical variables, UP-FHDI will directly terminate
	// Check case 2: If there exists all-missing rows, UP-FHDI will directly terminate
	// Check case 3: If there are less than two fully observed rows, UP-FHDI will directly terminate

	// original R code: Dr. Im, J. and Dr. Kim, J. 

	// c++ code: 		Yicheng Yang 

	// All rights reserved

	// 

	// updated: July 30, 2021

	//----------------------------------------------------

	//IN	: double** x_temp		= distributed daty (nrow,L_temp)

	//IN    : int** r_temp          = distributed datr (nrow,L_temp)

	//IN    : const int L_temp      = distributed number of columns on each processor  

	//=====================================================

	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//--------------------------
	//Check case 1: Same column of daty check
	//-------------------------

	if (mynode != 0) {
		std::vector<int> column_indicator;// 1 if not same; 0 if all same column

		for (int l = 0; l < L_temp; l++) {
			column_indicator.push_back(1);
		}

		std::vector<double> buffer;// hold one column of observed values of daty
		double min = 0.0;
		double max = 0.0;

		for (int j = 0; j < L_temp; j++) {
			buffer.clear();
			min = 0.0;
			max = 0.0;

			for (int i = 0; i < nrow; i++) {
				if (r_temp[i][j] != 0) {
					buffer.push_back(x_temp[i][j]);
				}
			}

			min = min_FHDI(buffer);
			max = max_FHDI(buffer);

			//if (mynode == 1) {
			//	cout<<"buffer at column "<<j<<endl;
			//	for (int t = 0; t < buffer.size(); t++) {
			//		cout << setw(20) << buffer[t] <<endl;
			//	}
			//	cout<<"min is "<< min <<" and max is "<< max <<endl;
			//}

			//if(mynode==1) cout << "min is " << min << " and max is " << max <<" at j = "<<j<< endl;

			if (fabs(max - min) < 1e-15) {
				//cout << "ERROR! The observed values of the " << startpoint + j << "th variable are constant! Please remove this variable to improve data quality! " << endl;
				column_indicator[j] = 0;
			}
		}

		//for (int l = 0; l < column_indicator.size(); l++) {
		//	cout << "column_indicator[" << l << "]: " << column_indicator[l]<<" at node "<<mynode<< endl;
		//}

		if (column_indicator.size() != L_temp) {
			cout<<"ERROR! The updating of column_indicator in Error_Check_Data function is incorrect! "<<endl;
		}

		MPI_Send(&column_indicator[0], L_temp, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}//end of slave processors


	std::vector<int> final_column_indicator;

	if (mynode == 0) {

		for (int j = 1; j < totalnodes; j++) {
			if (j != (totalnodes - 1)) {
				std::vector<int> column_indicator_recv(numWorkPerProc);
				MPI_Recv(&column_indicator_recv[0], numWorkPerProc, MPI_INT, j, 1, MPI_COMM_WORLD, &status);

				for (unsigned int k = 0; k < column_indicator_recv.size(); k++) {
					final_column_indicator.push_back(column_indicator_recv[k]);
				}
			}

			if (j == (totalnodes - 1)) {
				std::vector<int> column_indicator_recv(numWorkLocalLast);
				MPI_Recv(&column_indicator_recv[0], numWorkLocalLast, MPI_INT, j, 1, MPI_COMM_WORLD, &status);

				for (unsigned int k = 0; k < column_indicator_recv.size(); k++) {
					final_column_indicator.push_back(column_indicator_recv[k]);
				}
			}
		}//end of mynode != 0

		if (final_column_indicator.size() != ncol) cout<<"ERROR! Communication incorrect in Error_Check_Data function! "<<endl;

	}//end of master node

	if (mynode != 0) final_column_indicator.resize(ncol);
	MPI_Bcast(&final_column_indicator[0], ncol, MPI_INT, 0, MPI_COMM_WORLD);

	//for (int l = 0; l < final_column_indicator.size(); l++) {
	//	cout << "final_column_indicator[" << l << "]: " << final_column_indicator[l] <<" at node "<<mynode<< endl;
	//}

	for (int t = 0; t < ncol; t++) {
		//if (mynode == (totalnodes - 1)) cout << "final_column_indicator[" << t << "]: " << final_column_indicator[t] << endl;
		if (final_column_indicator[t] == 0) {
			TestOut << "ERROR! The " << t + 1 << "th column of daty is identical! Please remove this variable! " << endl;
			return 0;
		}
	}



	//--------------------------
	//Check case 2: Empty row of daty check
	//-------------------------

	if (mynode != 0) {
		std::vector<int> row_indicator; //sum of each row of datr
		int sum_temp = 0;

		for (int i = 0; i < nrow; i++) {
			sum_temp = 0;
			for (int j = 0; j < L_temp; j++) {
				sum_temp = sum_temp + r_temp[i][j];
			}
			row_indicator.push_back(sum_temp);
		}

		MPI_Send(&row_indicator[0], nrow, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}

	std::vector<int> final_row_indicator(nrow, 0);

	if (mynode == 0) {

		for (int j = 1; j < totalnodes; j++) {
			std::vector<int> row_indicator_recv(nrow);
			MPI_Recv(&row_indicator_recv[0], nrow, MPI_INT, j, 1, MPI_COMM_WORLD, &status);

			for (int i = 0; i < nrow; i++) {
				final_row_indicator[i] = final_row_indicator[i] + row_indicator_recv[i];
			}

		}//end of mynode!=0

	}//end of master node

	MPI_Bcast(&final_row_indicator[0], nrow, MPI_INT, 0, MPI_COMM_WORLD);

	for (int t = 0; t < nrow; t++) {
		//if(mynode==(totalnodes-1)) cout << "final_row_indicator[" << t << "]: " << final_row_indicator[t] << endl;
		if (final_row_indicator[t] == 0) {
			TestOut << "ERROR! The " << t + 1 << "th row of daty is all missing! Please remove this row! " << endl;
			return 0;
		}
	}

	//-----------------------------
	//Check case 3: Few fully observed rows check
	//-----------------------------

	int observed_row = 0;
	for (int t = 0; t < nrow; t++) {
		if (final_row_indicator[t] == ncol) {
			observed_row++;
		}
	}

	//cout<<"observed_row is "<< observed_row <<" at node "<<mynode<<endl;

	if (observed_row < 2) {
		TestOut << "ERROR! The raw data has less than two fully observed rows such that UP-FHDI can not find sufficient donors for each recipient! Please improve the data quality! " << endl;
		return 0;
	}

	return 1;
}