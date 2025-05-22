//#include "ran_FHDI.h" //for uniform distribution 
//#include <cstdlib> //for rand() and srand()
//#include <time.h>  //for time()

void Validation_cell_make_knn(const int nrow_uox, const int nrow_mox, const int nrow, const int ncol, const int max_overlap_size,

	MPI_File fh_datz, MPI_File fh_uox, MPI_File fh_mox, int** uox_info_final, int** mox_info_final,

	std::vector<int> v_nD, List_FHDI &List_nU, ofstream& TestOut)
	//Description=========================================
	// Find deficient donors for the recipient who has less than 2 donors by the Euclidean distance
	//
	// Algorithm: 											
	//   For a given missing row at i_reci                 e.g., {12, NA, 4}           
	//   Step 1: compute Euclidean distance from all unique observed patterns to the recipient  
	//   Step 2: Case 0: if the recipient has a donor, then randomly select another one from list
	//   Step 3: Case 1: if the recipient has no donor and the max occurnace of the candidate (who has minimum distance) in the list is >=2, then randomly select another one from list
	//   Step 4: Case 2: if the recipient has no donor and the max occurnace of the candidate (who has minimum distance) in the list is <2, and the size of the list is >=2, 
	//                   then randomly select two donors from list
	//   Step 5: Case 3: if the recipient has no donor and the max occurnace of the candidate (who has minimum distance) in the list is <2, and the size of the list is 1, 
	//                   then select another donor from the list of candicates who has the second minimum distance
	//
	//
	// original R code: Dr. Im, J. and Dr. Kim, J. 
	// c++ code: 		Yicheng Yang. 
	// All rights reserved
	// 
	// updated: July 16, 2020
	//----------------------------------------------------
	//IN    : int i_reci = location of row with missing cell that has the least number of donors

	//====================================================
{
	//-- MPI variables
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//-----------------------------
	//VALIDATION 
	//------------------------------

	if (mynode == 0) {
		//------------
		//uox and mox
		//-------------
		double** mox_vd = New_dMatrix(nrow_mox, ncol);
		double** uox_vd = New_dMatrix(nrow_uox, ncol);

		double*  d_read_out = new double[ncol];

		for (int k = 0; k < nrow_mox; k++) {

			MPI_In_uox_mox(ncol, k, fh_mox, d_read_out);

			for (int m = 0; m < ncol; m++) {
				mox_vd[k][m] = d_read_out[m];
			}
		}

		for (int j = 0; j < nrow_uox; j++) {

			MPI_In_uox_mox(ncol, j, fh_uox, d_read_out);

			for (int m = 0; m < ncol; m++) {
				uox_vd[j][m] = d_read_out[m];
			}
		}

		std::string *cn_mox = new std::string[nrow_mox]; //declaration of concatenated string vector of z
		std::string *cn_uox = new std::string[nrow_uox]; //declaration of concatenated string vector of z

		Trans(mox_vd, nrow_mox, ncol, cn_mox);
		Trans(uox_vd, nrow_uox, ncol, cn_uox);

		std::string *s_mox = new std::string[nrow_mox]; //string vector of observed patterns only
		std::string *s_uox = new std::string[nrow_uox]; //string vector of missing patterns only

		for (int l = 0; l < nrow_mox; l++) s_mox[l] = cn_mox[l];
		for (int l = 0; l < nrow_uox; l++) s_uox[l] = cn_uox[l];

		std::sort(s_mox, s_mox + nrow_mox); //knowing that s_ol[] has i_ol_temp entities
		std::sort(s_uox, s_uox + nrow_uox); //knowing that s_ml[] has i_ml_temp entities

		std::vector<int> order_mox;
		std::vector<int> order_uox;

		std::string s_temp;
		for (int t = 0; t < nrow_mox; t++) {
			s_temp = s_mox[t];
			for (int k = 0; k < nrow_mox; k++) {
				if (s_temp.compare(cn_mox[k]) == 0) {
					order_mox.push_back(k);
				}
			}
		}
		for (int t = 0; t < nrow_uox; t++) {
			s_temp = s_uox[t];
			for (int k = 0; k < nrow_uox; k++) {
				if (s_temp.compare(cn_uox[k]) == 0) {
					order_uox.push_back(k);
				}
			}
		}

		double** mox_sorted = New_dMatrix(nrow_mox, ncol);
		double** uox_sorted = New_dMatrix(nrow_uox, ncol);

		for (int t = 0; t < nrow_mox; t++) {
			for (int l = 0; l < ncol; l++) {
				mox_sorted[t][l] = mox_vd[order_mox[t]][l];
			}
		}

		for (int t = 0; t < nrow_uox; t++) {
			for (int l = 0; l < ncol; l++) {
				uox_sorted[t][l] = uox_vd[order_uox[t]][l];
			}
		}
		//-----------
		//z matrix
		//---------
		double** z_matrix_vd = New_dMatrix(nrow, ncol);
		double*   d_read_out2 = new double[ncol];

		for (int k = 0; k < nrow; k++) {

			MPI_In_datz(nrow, ncol, k, fh_datz, d_read_out2);// read distributed z matrix row by row

			for (int m = 0; m < ncol; m++) {
				z_matrix_vd[k][m] = d_read_out2[m];
			}
		}
		//-----------
		//List_nU
		//---------
		//const int max_donor_size = 10;
		//int** List_nU_temp_total = New_iMatrix(nrow_mox, max_donor_size);
		std::vector<int> List_row;
		cout << "List_nU at node " << mynode << endl;
		//for (int p = 0; p < nrow_mox; p++) {
		//	List_row.clear();
		//	List_nU.get_block_yicheng(order_mox[p], List_row);

		//	for (int j = 0; j < List_row.size(); j++) {
		//		cout << setw(20) << List_row[j];
		//	}
		//	cout << endl;
		//}

		for (int p = 0; p < nrow_mox; p++) {
			List_row.clear();
			List_nU.get_block_yicheng(order_mox[p], List_row);

			cout<<"order_mox["<<p<<"]: "<< order_mox[p] <<endl;
			for (int j = 0; j < List_row.size(); j++) {
				for (int k = 0; k < ncol; k++) {
					cout << setw(20) << uox_vd[List_row[j]-1][k];
				}
				cout << endl;
			}
		}
		//-----------
		//v_nD
		//---------
		cout << "final v_nD at node " << mynode << endl;

		for (int p = 0; p < nrow_mox; p++) {
			cout << "v_nD[" << p << "]: " << v_nD[order_mox[p]] << endl;
		}
		//cout << "final z_matrix_vd matrix from node " << mynode << endl;
		//for (int kk2 = 0; kk2 < nrow; kk2++) {
		//	for (int kk3 = 0; kk3 < ncol; kk3++) {
		//		cout << setw(20) << z_matrix_vd[kk2][kk3];
		//	}
		//	cout << endl;
		//}

		//cout << "final mox_sorted Reading matrix from node " << mynode << endl;
		//for (int kk2 = 0; kk2 < nrow_mox; kk2++) {
		//	for (int kk3 = 0; kk3 < ncol; kk3++) {
		//		cout << setw(20) << mox_sorted[kk2][kk3];
		//	}
		//	cout << endl;
		//}

		//cout << "final uox_sorted Reading matrix from node " << mynode << endl;
		//for (int kk2 = 0; kk2 < nrow_uox; kk2++) {
		//	for (int kk3 = 0; kk3 < ncol; kk3++) {
		//		cout << setw(20) << uox_sorted[kk2][kk3];
		//	}
		//	cout << endl;
		//}


		delete[] d_read_out;
		Del_dMatrix(mox_vd, nrow_mox, ncol);
		Del_dMatrix(uox_vd, nrow_uox, ncol);

		delete[] cn_mox;
		delete[] cn_uox;
		delete[] s_mox;
		delete[] s_uox;

		Del_dMatrix(mox_sorted, nrow_mox, ncol);
		Del_dMatrix(uox_sorted, nrow_uox, ncol);
		Del_dMatrix(z_matrix_vd, nrow, ncol);

		delete[] d_read_out2;
	}

	return; //temporary ending 
}
