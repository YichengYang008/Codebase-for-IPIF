//------------------------------------------------
//Parallel Mutual Information (P-MI)
//Last release date: June 10, 2021
//Developers: Yicheng Yang (Iowa State University)
//------------------------------------------------

using namespace std;

#include <iostream>
#include <iomanip>
#include <string>
#include <iostream> // For cout, cerr, endl, and flush 
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>    // std::sort
#include <math.h> 
#include<mpi.h>
#include <stdio.h>
#include <strstream>

#include "matrix_utility_FHDI.cc" //for local matrix utilities
#include "MPI_IO.cpp"
#include "ReadInput_FHDI_MPI.cc"


int main(int argc, char *argv[]) {

	MPI_Init(&argc, &argv);

	//-- MPI variables
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	double d_begin_MPI = MPI_Wtime();

	//--------------------
	// Input basic setup
	//---------------------

	ifstream* pFileInput = new ifstream
	("./input.txt");

	if (pFileInput->bad()) //cf. opposite member fn= good()
	{
		cout << "Can't open input file" << endl;
		return 0;
	}

	int nrow = 0;//user-defined
	int ncol = 0;//user-defined
	int K = 0;// Number of category

	int i_col_target_temp = 0; // target variable from 1 (1-indexed)
	int i_SIS_temp = 0; // Number of selected features

	ReadInput_FHDI_MPI(pFileInput, nrow, ncol, K, i_col_target_temp, i_SIS_temp);

	int i_col_target = 0; // Actual target variable from 0 (0-indexed)
	int i_SIS = 0; // Number of selected predictors
	i_col_target = i_col_target_temp - 1;
	i_SIS = i_SIS_temp - 1;

	//cout<<"nrow is "<< nrow <<" at node "<<mynode<<endl;
	//cout << "ncol is " << ncol << " at node " << mynode << endl;
	//cout << "K is " << K << " at node " << mynode << endl;
	//cout << "i_col_target is " << i_col_target << " at node " << mynode << endl;
	//cout << "i_SIS is " << i_SIS << " at node " << mynode << endl;
	//--------------------
	// Input file string
	//---------------------

	MPI_File fh_daty;//file string to read daty

	int success = 0;
	 
	success = MPI_File_open(MPI_COMM_WORLD, "./daty_column_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_daty);
	if (success != MPI_SUCCESS) cout << "MPI I/O of daty matrix fail to open the file!" << endl;

	//-------------------
	//TestOut files
	//-------------------

	ofstream TestOut_daty;//file string to write daty with selected features
	ofstream TestOut_id;//file string to write id with selected features

	if (mynode == 0) {
		TestOut_daty.open("./selected_daty.txt", ios::out);
		TestOut_id.open("./selected_id.txt", ios::out);
	}

	//==================
	//Setup
	//==================
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
			cout << "naive boundary ERROR!!!" << endl;
			return 0;
		}
	}

	if (mynode == (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = (mynode - 1)*numWorkPerProc + numWorkLocalLast;

		if (endpoint - startpoint != L_temp) {
			cout << "naive boundary ERROR!!!" << endl;
			return 0;
		}
	}

	//cout<<"startpoint is "<< startpoint <<" and endpoint is "<< endpoint <<" and L_temp is "<< L_temp <<" at node "<<mynode<< endl;

	if (i_col_target >= startpoint && i_col_target < endpoint) {
		L_temp = L_temp - 1;
		//cout<<"New L_temp is "<< L_temp <<" where startpoint is "<< startpoint <<" and endpoint is "<< endpoint <<" at node "<<mynode<<endl;
	}


	//================================================
	//Read distributed data and target variable
	//=================================================
	double* array_temp1 = NULL; 
	double* array_target = NULL;
	double** daty = NULL;  


	if (mynode != 0) {
		array_temp1 = new double[nrow];	// daty buffer
		array_target = new double[nrow];	// target variable buffer
		daty = New_dMatrix(nrow, L_temp);// distributed daty 

		// Read target variable on all slave processors
		MPI_In_raw(nrow, i_col_target, fh_daty, array_target);

		// Read distributed daty on all slave processors. Note that target variable exclusive
		int counter = 0;
		for (int k = startpoint; k < endpoint; k++) {

			if (k != i_col_target) {
				MPI_In_raw(nrow, k, fh_daty, array_temp1);

				for (int m = 0; m < nrow; m++) {
					daty[m][counter] = array_temp1[m];
				}
				counter++;
			}
		}

	}//end of mynode
	
	MPI_Barrier(MPI_COMM_WORLD);

	//if (mynode == 1) {
	//		cout << "daty at node " << mynode << endl;
	//		for (int kk2 = 0; kk2 < nrow; kk2++) {
	//			for (int kk3 = 0; kk3 < L_temp; kk3++) {
	//				cout << setw(20) << daty[kk2][kk3];
	//			}
	//			cout << endl;
	//		}

	//		cout << "array_target at node " << mynode << endl;
	//		for (int kk2 = 0; kk2 < nrow; kk2++) {
	//			cout<<"array_target["<<kk2<<"]: "<< array_target[kk2]<<endl;
	//		}
	//}

	//===============================
	//Discretization 
	//================================
	double** datz = NULL; 
	double* targetz = NULL;

	if (mynode != 0) {

		datz = New_dMatrix(nrow, L_temp);// Categorized matrix
		targetz = new double[nrow];	// Categorized target


		double* x_one_column = new double[nrow]; Fill_dVector(x_one_column, nrow, 0.0);
		double* x_one_column_temp = new double[nrow]; Fill_dVector(x_one_column_temp, nrow, 0.0);

		//---------------------------------
		//Categorize distributed daty
		//--------------------------------

		for (int i_col = 0; i_col < L_temp; i_col++) {
		
			for (int i = 0; i<nrow; i++) x_one_column[i] = daty[i][i_col]; //get one column
			for (int i = 0; i<nrow; i++) x_one_column_temp[i] = daty[i][i_col]; //get one column


			double* perc = new double[K - 1]; Fill_dVector(perc, (K - 1), 0.0);


			for (int i = 0; i<(K - 1); i++)

			{

				perc[i] = (i + 1)*(1.0 / K);

			}

			std::sort(x_one_column_temp, x_one_column_temp + nrow);

			//Note: the last quantile (i.e. 100%) is not included, and thus (k_one_column-1) is used

			double* x_quantile = new double[K - 1]; Fill_dVector(x_quantile, (K - 1), 0.0);

			for (int i = 0; i<(K - 1); i++)

			{

				double d_h = (nrow - 1)*perc[i]; //+1 is removed for c++ code 

				x_quantile[i] = x_one_column_temp[int(floor(d_h))]

					+ (d_h - floor(d_h))*(x_one_column_temp[int(floor(d_h) + 1)]

						- x_one_column_temp[int(floor(d_h))]);

				//if (mynode != 0) cout<< "d_h: "<< d_h <<" and n_observed: "<< n_observed <<" at column "<< i_col <<endl;

			}

			for (int i = 0; i < nrow; i++) {
				//---------

				//default category of non-NA unit is 1 as of 0124_2017

				//---------

				datz[i][i_col] = 1; //default 



				//----------

				//consider each quantile

				//----------

				if (x_one_column[i] < x_quantile[0]) { datz[i][i_col] = 1; } //1st category

				if (x_one_column[i] > x_quantile[K - 2]) { datz[i][i_col] = K; } //last category



				for (int j = 1; j<(K - 1); j++)

				{

					if (x_quantile[j - 1] < x_one_column[i] && x_one_column[i] <= x_quantile[j])

					{

						datz[i][i_col] = j + 1; //(j+1)th quantile. Note: j =[0,k_one_column) 

						break;

					}

				}
			}

			delete[] perc;
			delete[] x_quantile;

		}//end of main loop 



		//------------------------------
		//Categorize target variable
		//---------------------------------

		double* array_target_temp = new double[nrow]; Fill_dVector(array_target_temp, nrow, 0.0);
		for (int i = 0; i<nrow; i++) array_target_temp[i] = array_target[i]; //get one column

		double* perc = new double[K - 1]; Fill_dVector(perc, (K - 1), 0.0);


		for (int i = 0; i<(K - 1); i++)

		{

			perc[i] = (i + 1)*(1.0 / K);

		}

		std::sort(array_target_temp, array_target_temp + nrow);

		//Note: the last quantile (i.e. 100%) is not included, and thus (k_one_column-1) is used

		double* x_quantile = new double[K - 1]; Fill_dVector(x_quantile, (K - 1), 0.0);

		for (int i = 0; i<(K - 1); i++)

		{

			double d_h = (nrow - 1)*perc[i]; //+1 is removed for c++ code 

			x_quantile[i] = array_target_temp[int(floor(d_h))]

				+ (d_h - floor(d_h))*(array_target_temp[int(floor(d_h) + 1)]

					- array_target_temp[int(floor(d_h))]);

			//if (mynode != 0) cout<< "d_h: "<< d_h <<" and n_observed: "<< n_observed <<" at column "<< i_col <<endl;

		}

		for (int i = 0; i < nrow; i++) {
			//---------

			//default category of non-NA unit is 1 as of 0124_2017

			//---------

			targetz[i]= 1; //default 



			//----------

			//consider each quantile

			//----------

			if (array_target[i] < x_quantile[0]) { targetz[i] = 1; } //1st category

			if (array_target[i] > x_quantile[K - 2]) { targetz[i] = K; } //last category



			for (int j = 1; j<(K - 1); j++)

			{

				if (x_quantile[j - 1] < array_target[i] && array_target[i] <= x_quantile[j])

				{

					targetz[i] = j + 1; //(j+1)th quantile. Note: j =[0,k_one_column) 

					break;

				}

			}
		}

		delete[] perc;
		delete[] x_quantile;
		delete[] array_target_temp;


		delete[] x_one_column;
		delete[] x_one_column_temp;

	}// end of mynode

	//if (mynode == 3) {
	//	cout << "datz at node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < nrow; kk2++) {
	//		for (int kk3 = 0; kk3 < L_temp; kk3++) {
	//			cout << setw(20) << datz[kk2][kk3];
	//		}
	//		cout << endl;
	//	}

	//	cout << "targetz at node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < nrow; kk2++) {
	//		cout << "targetz[" << kk2 << "]: " << targetz[kk2] << endl;
	//	}
	//}

	//==============================
    //Compute MIs
	//==============================
	std::vector<double> MI_send_temp; //MIs on each slave processor

	if (mynode != 0) {

		double* datz_one_column = new double[nrow]; Fill_dVector(datz_one_column, nrow, 0.0);

		double MI_temp = 0.0;
		double P_x = 0.0; //marginal probability of x
		double P_y = 0.0; //marginal probability of target
		double P_xy = 0.0; //joint probability of x and target

		for (int l = 0; l < L_temp; l++) {

			MI_temp = 0.0;

			for (int i = 0; i<nrow; i++) datz_one_column[i] = datz[i][l]; //get one column

			for (int k1 = 0; k1 < K; k1++) {

				//Compute P_x
				P_x = 0.0;
				int P_x_count = 0; 
				for (int j = 0; j < nrow; j++) {
					if (fabs(datz_one_column[j] - (k1+1)) < 1e-15) {
						P_x_count++;
					}
				}



				P_x = (double)P_x_count / nrow;
				//cout << "P_x_count is " << P_x_count <<" and nrow is "<< nrow <<" and P_x is "<< P_x << " at k1 = " << k1 << " l = " << l << endl;

				for (int k2 = 0; k2 < K; k2++) {
				   
					//Compute P_y
					P_y = 0.0;
					int P_y_count = 0;
					for (int j = 0; j < nrow; j++) {
						if (fabs(targetz[j] - (k2+1)) < 1e-15) {
							P_y_count++;
						}
					}

					P_y = (double)P_y_count / nrow;

					//Compute P_xy
					P_xy = 0.0;
					int P_xy_count = 0;
					for (int t = 0; t < nrow; t++) {
						if ( (fabs(datz_one_column[t] - (k1 + 1)) < 1e-15) && (fabs(targetz[t] - (k2 + 1)) < 1e-15) ) {
							P_xy_count++;
						}
					}

					P_xy = (double)P_xy_count / nrow;
					//if (mynode == 1) cout << "P_x is " << P_x << " and P_y is " << P_y<<" and P_xy is " << P_xy <<" at k1 = "<<k1<<" k2 = "<<k2<<" l = "<<l<< endl;

					//Compute MI
					MI_temp = MI_temp + (double)P_xy*log2(P_xy/(P_x * P_y)); //Note that we should use log with a base of 2

				}
			}

			//if (mynode == 1) cout<<"MI_temp at l = "<<l<<" is "<< MI_temp <<endl;

			MI_send_temp.push_back(MI_temp);

		}//main loop

		//Note that we should let the master know the number of MIs from each slave processor
		//Because one of slave processors has a size of L_temp-1 (target variable exclusive) but we don't know its id. 

		int MI_send_size = 0;
		MI_send_size = (int)MI_send_temp.size();

		MPI_Send(&MI_send_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);

		delete[] datz_one_column;
	}//end of mynode

	std::vector<int> MI_size; // number of MIs of all slave processors
	if (mynode == 0) {
	
		int MI_send_size_recv = 0.0;

		for (int j = 1; j < totalnodes; j++) {
			MPI_Recv(&MI_send_size_recv, 1, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
			MI_size.push_back(MI_send_size_recv);
		}
	
		int MI_size_total = 0;

		for (unsigned k = 0; k < MI_size.size(); k++) {

			MI_size_total = MI_size_total + MI_size[k];
			//cout<<"MI_size["<<k<<"]: "<< MI_size[k]<<" at node "<<k+1<<endl;
		}

		if (MI_size_total != (ncol - 1)) cout<<"ERROR!!! The MI_size_total is incorrect!!!"<<endl;

	}

	MPI_Barrier(MPI_COMM_WORLD);


	//Send MIs on slave processors to the master
	//Now we should know who excludes the target variable

	if (mynode != 0) {
		MPI_Send(&MI_send_temp[0], L_temp, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
	}

	if (mynode == 0) {

		//---------------
		//All_gather MIs 
		//---------------
		std::vector<double> MIs; // final MIs 

		for (int j = 1; j < totalnodes; j++) {

			int MI_size_temp = 0;
			MI_size_temp = MI_size[j - 1];

			std::vector<double> MI_Recv_temp(MI_size_temp); //MI buufer

			//cout << "MI_size_temp is " << MI_size_temp << " at node " << j << endl;

			MPI_Recv(&MI_Recv_temp[0], MI_size_temp, MPI_DOUBLE, j, 2, MPI_COMM_WORLD, &status);

			for (int t = 0; t < MI_size_temp; t++) {
				MIs.push_back(MI_Recv_temp[t]);
			}

		}

		if (MIs.size() != (ncol - 1)) cout << "ERROOR in computing MIs!!!!" << endl;

		for (unsigned l = 0; l < MIs.size(); l++) {
			if (MIs[l] < 1e-15) cout << "ERROR !!! MI[" << l << " is less than 0!!!" << endl;
			//cout << "MIs[" << l << "]: " << MIs[l] << endl;
		}


		//-------------------------------------------------
		// Select i_SIS features who have the largest MIs
		//--------------------------------------------------

		std::vector<int> Predictor; // index of all predictors excluding the target variables from 0
		std::vector<int> selected_Predictor;// i_SIS selected features

		for (int k = 0; k < ncol; k++) {
			if (k != i_col_target) {
				Predictor.push_back(k);
			}
		}

		//for (int k = 0; k < Predictor.size(); k++) {
		//	cout << "Predictor[" << k << "]: " << Predictor[k]<< endl;
		//}

		//---------------------------------------------
		//Find the largest MIs and corresponding index
		//------------------------------------------------

		for (int t = 0; t < i_SIS; t++) {
			double max_temp = 0.0;
			int max_index = 0;
			for (int j = 0; j < (ncol - 1); j++) {
				if (max_temp < MIs[j]) {
					max_temp = MIs[j];
					max_index = j;
				}
			}

			//if (t == 3) {
			//	for (unsigned l = 0; l < MIs.size(); l++) {
			//		cout << "MIs[" << l << "]: " << MIs[l] << endl;
			//	}
			//	cout<<"min_index is "<< min_index <<endl;
			//}


			MIs[max_index] = 0.0;

			selected_Predictor.push_back(Predictor[max_index]);
		}

		//for (unsigned j = 0; j < selected_Predictor.size(); j++) {
		//	cout << "selected_Predictor[" << j << "]: " << selected_Predictor[j] << endl;
		//}

		// Sort selected features for easier reading later
		sort(selected_Predictor.begin(), selected_Predictor.end());

		//Always add the target variables as the last column
		selected_Predictor.push_back(i_col_target);

		//cout<<"Indices of selected features: "<<endl;
		//for (unsigned j = 0; j < selected_Predictor.size(); j++) {
		//	cout << setw(20) <<selected_Predictor[j];
		//}
		//cout << endl;

		if (selected_Predictor.size() != (i_SIS + 1)) cout<<"ERROR!!!! You may not add the target variable!!!"<<endl;
		
		//Read selected features
		double* array_temp2 = new double[nrow];	// daty buffer
		double** selected_daty = New_dMatrix(nrow, i_SIS + 1);// distributed daty 

		int counter = 0;
		for (int k = 0; k < (i_SIS + 1); k++) {

			int i_temp = selected_Predictor[k];

			MPI_In_raw(nrow, i_temp, fh_daty, array_temp2);

			for (int m = 0; m < nrow; m++) {
				selected_daty[m][counter] = array_temp2[m];
			}
			counter++;
		}

		//cout << "selected_daty at node " << mynode << endl;
		//for (int kk2 = 0; kk2 < nrow; kk2++) {
		//	for (int kk3 = 0; kk3 < (i_SIS + 1); kk3++) {
		//		cout << setw(20) << selected_daty[kk2][kk3];
		//	}
		//	cout << endl;
		//}

		//Write selected features (predictors + target) to hard drive colum-wisely
		//MPI_Out_raw(nrow, i_SIS + 1, fh_selected_daty, selected_daty);
		for (int i = 0; i < nrow; i++) {
			for (int j = 0; j < (i_SIS + 1); j++) {
				TestOut_daty << setw(20) << selected_daty[i][j];
			}
			TestOut_daty << endl;
		}


		//Write selected ids to hard drive
		//MPI_Out_id(i_SIS + 1, fh_selected_id, selected_Predictor);
		for (int j = 0; j < (i_SIS + 1); j++) {
			TestOut_id << setw(20) << selected_Predictor[j] + 1; //Print out 1-indexed IDs
		}
		TestOut_id << endl;

		delete[] array_temp2;
		Del_dMatrix(selected_daty, nrow, i_SIS + 1);

	}


	//--------------------
	// Clsoe file strings
	//-------------------
	success = MPI_File_close(&fh_daty);
	if (success != MPI_SUCCESS) cout << "MPI I/O of daty matrix fail to close the file!" << endl;
	
	TestOut_daty.close();
	TestOut_id.close();

	double d_end_MPI = MPI_Wtime();
	//if (mynode == 0) cout << "YYC Total running time = " << d_end_MPI - d_begin_MPI;
	if (mynode == 0) cout <<"Parallel MI has successfully finished !!!! "<<endl;
	//-----------------
	//Deallocation
	//--------------------

	if (mynode != 0) {

		delete[] array_temp1;
		delete[] array_target;
		delete[] targetz;

		Del_dMatrix(daty, nrow, L_temp);
		Del_dMatrix(datz, nrow, L_temp);
	}



	//cin.get();
	MPI_Finalize();
	return 0;

}