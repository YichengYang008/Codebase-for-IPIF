using namespace std;

#include <iostream>
#include <iomanip>
#include <string>
#include <iostream> // For cout, cerr, endl, and flush 
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>
#include<mpi.h>
#include <stdio.h>
#include <limits.h>

#include "matrix_utility_FHDI.cc" //for local matrix utilities
#include "MPI_IO.cpp"

int main(int argc, char *argv[]) {

	MPI_Init(&argc, &argv);

	//-- MPI variables
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	double d_begin_MPI = MPI_Wtime();

	MPI_File fh_daty_full;// File to read daty column-wisely
	MPI_File fh_daty_column;
	MPI_File fh_daty_row;
	MPI_File fh_datr;

	int success = 0;

	success = MPI_File_open(MPI_COMM_WORLD, "./full_daty_column_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_daty_full);
	if (success != MPI_SUCCESS) cout << "MPI I/O fail to open the file!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./daty_column_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_daty_column);
	if (success != MPI_SUCCESS) cout << "MPI I/O of daty matrix fail to open the file!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./daty_row_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_daty_row);
	if (success != MPI_SUCCESS) cout << "MPI I/O of daty matrix fail to open the file!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./datr_column_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_datr);
	if (success != MPI_SUCCESS) cout << "MPI I/O of datr matrix fail to open the file!" << endl;

	//==================
	//Setup
	//==================

	const int nrow = 24016;
	const int ncol = 2401;
	double rho = 0.5;

	int numWorkPerProc = (int)floor(1.0*ncol / (1.0*totalnodes - 1));
	int numWorkLocalLast = ncol - numWorkPerProc * (totalnodes - 2);

	int L_temp = 0;
	if (mynode != (totalnodes - 1) && mynode != 0) L_temp = numWorkPerProc;
	if (mynode == (totalnodes - 1)) L_temp = numWorkLocalLast;


	int startpoint = 0;
	int endpoint = 0;

	if (mynode != 0 && mynode != (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = mynode*numWorkPerProc;
		if (endpoint - startpoint != L_temp) {
			return 0;
		}
	}

	if (mynode == (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = (mynode - 1)*numWorkPerProc + numWorkLocalLast;
		if (endpoint - startpoint != L_temp) {
			return 0;
		}
	}

	cout << "L_temp is " << L_temp <<", startpoint is "<< startpoint <<", endpoint is "<< endpoint << " at node " << mynode << endl;

	//random number generator
	//std::default_random_engine generator;
	std::random_device rd;
	std::mt19937 generator{ rd() }; // or std::default_random_engine e{rd()};


	//---------------------------------------
	//Read distributed complete daty matrix
	//---------------------------------------

	double **daty = New_dMatrix(nrow, L_temp);
	double* array_temp = new double[nrow];	// Final full matrix

	int counter = 0;
	for (int k = startpoint; k < endpoint; k++) {

		MPI_In_raw(nrow, k, fh_daty_full, array_temp);

		for (int m = 0; m < nrow; m++) {
			daty[m][counter] = array_temp[m];
		}
		counter++;
	}

	//if (mynode == 1) {
	//	cout<<"daty at node "<<mynode<<endl;
	//	for (int i = 0; i < nrow; i++) {
	//		for (int j = 0; j < L_temp; j++) {
	//			cout<< setw(20) << daty[i][j];
	//		}
	//		cout << endl;
	//	}
	//}

	delete[] array_temp;

	//----------------------------------------------
	//Generate datr matrix for first 70% of datr
	//----------------------------------------------
	int nrow_sample = 0;
	nrow_sample = floor(0.8*nrow);
	//cout <<"nrow_sample is "<< nrow_sample <<endl;
	int **datr1 = New_iMatrix(nrow_sample, L_temp);// First 70% of datr
	int bernoulli_number = 0;

	double p = 0.0; // probability of 1s using bernoulli distribution

	for (int t = 0; t < L_temp; t++) {
		if ( (t+1) % 4 == 1) {
			for (int i = 0; i < nrow_sample; i++) {
				p = 1 / (1 + exp(1 + 0.0005*abs(daty[i][t])));
				if ((p - 1.0) > 1e-15) cout << "ERROR in computing probability p at t = " << t << endl;
				if (t==0 && mynode==1) cout<<"p is "<< p <<" at t = "<<t<<endl;

				std::bernoulli_distribution b(p);//bernoulli distribution
				bernoulli_number = b(generator);
				datr1[i][t] = bernoulli_number;
			}
		}
		if ((t + 1) % 4 == 2) {
			for (int i = 0; i < nrow_sample; i++) {
				p = 1 / (1 + exp(1 + 0.0005*abs(daty[i][t])));
				if ((p - 1.0) > 1e-15) cout << "ERROR in computing probability p at t = " << t << endl;
				//if (k == 0) cout << "p is " << p << " at t = " << t << endl;

				std::bernoulli_distribution b(p);//bernoulli distribution
				bernoulli_number = b(generator);
				datr1[i][t] = bernoulli_number;
			}
		}
		if ((t + 1) % 4 == 3) {
			for (int i = 0; i < nrow_sample; i++) {
				p = 1 / (1 + exp(1 + 0.0005*abs(daty[i][t])));
				if ((p - 1.0) > 1e-15) cout << "ERROR in computing probability p at t = " << t << endl;
				//if (k == 0) cout << "p is " << p << " at t = " << t << endl;

				std::bernoulli_distribution b(p);//bernoulli distribution
				bernoulli_number = b(generator);
				datr1[i][t] = bernoulli_number;
			}
		}
		if ((t + 1) % 4 == 0) {
			for (int i = 0; i < nrow_sample; i++) {
				p = 1 / (1 + exp(1 + 0.0005*abs(daty[i][t])));
				if ((p - 1.0) > 1e-15) cout << "ERROR in computing probability p at t = " << t << endl;
				//if (k == 0) cout << "p is " << p << " at t = " << t << endl;

				std::bernoulli_distribution b(p);//bernoulli distribution
				bernoulli_number = b(generator);
				datr1[i][t] = bernoulli_number;
			}
		}
	}


	if (mynode == totalnodes - 1) {
		for (int i = 0; i < nrow_sample; i++) {
			datr1[i][L_temp - 1] = 1;
		}
	}
	//------------------------------
	// Check empty rows
	//-----------------------------

	if (mynode != 0) {
		std::vector<int> row_indicator; //sum of each row of datr
		int check = 0;

		for (int i = 0; i < nrow_sample; i++) {
			check = 0;
			for (int l = 0; l < L_temp; l++) {
				check = check + datr1[i][l];
			}
			row_indicator.push_back(check);
		}

		MPI_Send(&row_indicator[0], nrow_sample, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}

	std::vector<int> final_row_indicator(nrow_sample, 0);

	if (mynode == 0) {

		for (int j = 1; j < totalnodes; j++) {
			std::vector<int> row_indicator_recv(nrow_sample);
			MPI_Recv(&row_indicator_recv[0], nrow_sample, MPI_INT, j, 1, MPI_COMM_WORLD, &status);

			for (int i = 0; i < nrow_sample; i++) {
				final_row_indicator[i] = final_row_indicator[i] + row_indicator_recv[i];
			}

		}//end of mynode!=0

	}//end of master node

	MPI_Bcast(&final_row_indicator[0], nrow_sample, MPI_INT, 0, MPI_COMM_WORLD);

	//if (mynode == 2) {
	//	for (int t = 0; t < final_row_indicator.size(); t++) {
	//		cout<<"final_row_indicator["<<t<<"]: "<< final_row_indicator[t] <<endl;
	//	}
	//}

	MPI_Barrier(MPI_COMM_WORLD);

	//Change all-missing rows to fully observed
	if (mynode != 0) {
		for (int i = 0; i < nrow_sample; i++) {
			if (final_row_indicator[i] == 0) {
				if (mynode == 1) { cout << "CAUTION! The " << i + 1 << "th row of datr1 is all zero!" << endl; }
				for (int l = 0; l < L_temp; l++) {
					datr1[i][l] = 1;
				}
			}
		}
	}


	//TestOut
	//if (mynode == 1) {
	//	cout << "datr1 at node " << mynode << endl;
	//	for (int i = 0; i < nrow_sample; i++) {
	//		for (int j = 0; j < L_temp; j++) {
	//			cout << setw(20) << datr1[i][j];
	//		}
	//		cout << endl;
	//	}
	//}


	//Combine to final datr
	int **datr = New_iMatrix(nrow, L_temp);
	Fill_iMatrix(datr, nrow, L_temp, 1);

	//Replace first 70% with datr1
	for (int i = 0; i < nrow_sample; i++) {
		for (int j = 0; j < L_temp; j++) {
			datr[i][j] = datr1[i][j];
		}
	}

	//TestOut
	// if (mynode == 1) {
	// 	cout << "Final datr at node " << mynode << endl;
	// 	for (int i = 0; i < nrow; i++) {
	// 		for (int j = 0; j < L_temp; j++) {
	// 			cout << setw(20) << datr[i][j];
	// 		}
	// 		cout << endl;
	// 	}
	// }



	//----------------------------------------------`
	//Reflect datr to daty
	//----------------------------------------------
	for (int i = 0; i < nrow_sample; i++) {
		for (int j = 0; j < L_temp; j++) {
			if (datr[i][j] == 0) {
				daty[i][j] = 0.0;
			}
		}
	}

	// if (mynode == 1) {
	// 	cout << "Final daty at node " << mynode << endl;
	// 	for (int i = 0; i < nrow; i++) {
	// 		for (int j = 0; j < L_temp; j++) {
	// 			cout << setw(20) << daty[i][j];
	// 		}
	// 		cout << endl;
	// 	}
	// }

	MPI_Out_daty(nrow, L_temp, numWorkPerProc, fh_daty_column, daty);
	MPI_Out_datz_rowise(nrow, L_temp, ncol, numWorkPerProc, fh_daty_row, daty);

	MPI_Out_datr(nrow, L_temp, numWorkPerProc, fh_datr, datr);

	success = MPI_File_close(&fh_daty_full);
	if (success != MPI_SUCCESS) cout << "MPI I/O of full daty matrix fail to close the file!" << endl;
	success = MPI_File_close(&fh_daty_column);
	if (success != MPI_SUCCESS) cout << "MPI I/O of daty matrix fail to close the file!" << endl;
	success = MPI_File_close(&fh_daty_row);
	if (success != MPI_SUCCESS) cout << "MPI I/O of daty matrix fail to close the file!" << endl;
	success = MPI_File_close(&fh_datr);
	if (success != MPI_SUCCESS) cout << "MPI I/O of datr matrix fail to close the file!" << endl;

	cout << "Generating ultra data in size of " << nrow << " and " << ncol << " has successfully finished at node " << mynode << endl;
	double d_end_MPI = MPI_Wtime();
	if (mynode == 0) cout << "YYC Total running time = " << d_end_MPI - d_begin_MPI;
	//Deallocation
	Del_iMatrix(datr1, nrow_sample, L_temp);
	Del_dMatrix(daty, nrow, L_temp);
	Del_iMatrix(datr, nrow, L_temp);

	//cin.get();
	MPI_Finalize();
	return 0;

}