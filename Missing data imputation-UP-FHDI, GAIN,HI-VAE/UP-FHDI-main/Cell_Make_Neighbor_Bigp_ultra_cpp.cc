//Fn===========================================================================

//Cell_Make_Neighbor_cpp.cc-----------------------------------------------------------------------------

//Fn===========================================================================
#include "nDAU_Bigp_ultra_cpp_MPI.cc"
#include "Validation_cell_make_bigp_knn.cc"

bool Cell_Make_Neighbor_Bigp_ultra_cpp(const int nrow, const int ncol, const int memory, double* d_k, int* NonCollapsible_categorical,

	MPI_File fh_binary_daty, MPI_File fh_binary_datr,

	List_FHDI &uox_infor, List_FHDI &mox_infor, int &nrow_uox, int &nrow_mox, int** codes,

	List_FHDI &List_nU, const int i_merge, const int i_collapsing, const int i_SIS_type, int top,

	ofstream& TestOut)

	//Description=========================================

	// make cells with the raw data matrix x 

	// categorization takes place 

	// according to the given number of categories stored in d_k(ncol)

	//

	// Algorithm:  for categorization

	// perc: percentiles used to get quantiles, determined by k

	// quan: quantiles if k=4, we quan=(Q1,Q2,Q3) have Q1(=1/4), Q2 (=Median) and Q3(=3/4)

	// 

	// Note: as of Oct 2016, NA values (missing data) is marked by a long integer at the parent "r" code

	//

	// original R code: Dr. Im, J. and Dr. Kim, J. 

	// c++ code: 		Yicheng Yang and Dr. Cho, I. 

	// All rights reserved

	// 

	// updated: August 7, 2021

	//----------------------------------------------------

	//IN	: double d_k(ncol)		= a vector of categories of each column of xalloc

	//IN    : int NonCollapsible_categorical(nrol) = {0,0, .., 1,.. 0} 
	//				index for non-collapsible categorical variables. 
	//				when at least one column has "1" skip cell-collapse procedure
	//				this may casue a potential error of lack of enough donor! 
	//				(2018, 04 21) 
	//											  

	//IN    : List_FHDI uox_infor = Actual index list of uox in z

	//IN    : List_FHDI mox_infor = Actual index list of mox in z

	//IN    : int nrow_uox      = number of rows of uox

	//IN    : int nrow_mox      = number of rows of mox

	//IN    : int i_collapsing  = number of selected variables

	//IN    : int i_SIS_type    = SIS methods: 1. intersection; 2. union; 3. global ranking

	//IN    : int top           = number of top correlation and corrlation ranking for all variables

	//OUT	: MPI_File fh_uox   = compact storage of uox, unique observed rows written in local storage

	//OUT	: MPI_File fh_mox   = compact storage of mox, unique missing rows written in local storage

	//OUT   : List_FHDI List_nU = actual locations of donor in uox of all mox

	//OUT   : codes             = selected variables of all mox

	//IN    : int i_merge = random donor selection in Merge algorithm in Cell Make

	//                   0= no random seed number setting

	//			         1= random seed number setting 

	//====================================================

{
	//-- MPI variables
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	double category_begin = MPI_Wtime();

	//open file string to write z matrix
	int success = 0;
	MPI_File fh_datz;

	success = MPI_File_open(MPI_COMM_WORLD, "./datz_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_datz);
	if (success != MPI_SUCCESS) cout << "MPI I/O of z matrix fail to open the file!" << endl;

	bool b_success_categorize = 0;//0 is False and 1 is True

	b_success_categorize = categorize_ultra_cpp(fh_binary_daty, fh_binary_datr, nrow, ncol, d_k,
		NonCollapsible_categorical, fh_datz, TestOut);

	//close file string to write z matrix
	success = MPI_File_close(&fh_datz);
	if (success != MPI_SUCCESS) cout << "MPI I/O of z matrix fail to close the file!" << endl;

	if (!b_success_categorize)
	{
		//early deallocation 
		//delete[] d_k_Collapsible;

		return 0;
	}

	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_make at node " << mynode << " = " << MPI_Wtime() - cell_make_begin << endl;
	//	printf("Yicheng Running time of category = %f seconds\n", MPI_Wtime() - category_begin);
	//}

	double Zmat_begin = MPI_Wtime();

	//open file string to read z matrix
	MPI_File fh_datz_reading;
	//open file string to write uox and mox
	MPI_File fh_uox;
	MPI_File fh_mox;

	success = MPI_File_open(MPI_COMM_WORLD, "./datz_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_datz_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of z matrix fail to open the file!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./uox_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_uox);
	if (success != MPI_SUCCESS) cout << "MPI I/O of uox fail to open the file!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./mox_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_mox);
	if (success != MPI_SUCCESS) cout << "MPI I/O of mox fail to open the file!" << endl;

	bool b_success_Zmat = 0;//0 is False and 1 is True

	b_success_Zmat = Zmat_Extension_ultra_cpp(nrow, ncol, memory,
		uox_infor, mox_infor, nrow_uox, nrow_mox,
		fh_datz_reading, fh_uox, fh_mox,
		TestOut);

	//close file string to read z matrix
	success = MPI_File_close(&fh_datz_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of z matrix fail to close the file!" << endl;

	//clsoe file string to write uox and mox
	success = MPI_File_close(&fh_uox);
	if (success != MPI_SUCCESS) cout << "MPI I/O of uox fail to close the file!" << endl;

	success = MPI_File_close(&fh_mox);
	if (success != MPI_SUCCESS) cout << "MPI I/O of mox fail to close the file!" << endl;

	if (!b_success_Zmat) {
		return 0;
	}

	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_make at node " << mynode << " = " << MPI_Wtime() - cell_make_begin << endl;
	//	printf("Yicheng Running time of Zmat = %f seconds\n", MPI_Wtime() - Zmat_begin);
	//}

	//cout<<"nrow_uox is "<< nrow_uox <<" and nrow_mox is "<< nrow_mox <<" at node "<<mynode<<endl;

	//---------------------------------------------------
	//Prepare actual locations of fully observed y matrix
	//from the actual locations stored in uox_info_final
	//-------------------------------------------------
	std::vector<int> v_ol;
	std::vector<double> v_ol_temp;

	uox_infor.unlist(v_ol_temp);

	int i_v_ol = 0;
	for (int k = 0; k < v_ol_temp.size(); k++) {
		i_v_ol = 0;
		i_v_ol = (int)v_ol_temp[k];
		if (i_v_ol < 1 || i_v_ol > nrow) { cout<<"ERROR! i_v_ol is incorrect in cell make neighbor! "<<endl; }

		v_ol.push_back(i_v_ol);
	}

	if (v_ol.size() != v_ol_temp.size()) { cout<<"ERROR! v_ol is incorrect in cell make neighbor! "<<endl; }

	//for (int k = 0; k < nrow_uox; k++) {
	//	for (int j = 2; j < (2 + max_overlap_size); j++) {
	//		if (uox_info_final[k][j] != 0) {
	//			v_ol.push_back(uox_info_final[k][j]);
	//		}
	//	}
	//}

	sort(v_ol.begin(), v_ol.end());
	//if (mynode == 1) {
	//	for (int i = 0; i < v_ol.size();i++) cout << "v_ol[" << i << "]: " << v_ol[i] << " at node " << mynode << endl;
	//}


	//int top = 100; // Top "top" rankings exclude itself. Note that  i_option_collapsing < top <= ncol-1
	if ((ncol - 1) < top) {
		top = ncol - 1;
	}

	//TestOut << "top_correlation in cell make bigp with KNN is " << top << endl;

	int** correlation_ranking_top = New_iMatrix(ncol, top);
	double** correlation_top = New_dMatrix(ncol, top);

	double ranking_m_time = MPI_Wtime();
	Ranking_top_ultra(nrow, ncol, memory, top, v_ol, correlation_top, correlation_ranking_top, TestOut);

	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_make at node " << mynode << " = " << MPI_Wtime() - cell_make_begin << endl;
	//	printf("Yicheng Running time of Ranking_top_ultra for ultra data = %f seconds\n", MPI_Wtime() - ranking_m_time);
	//}

	//if (mynode == 2) {
	//	cout << "final correlation_top matrix at node " << mynode << " where top is " << top << endl;
	//	for (int kk2 = 0; kk2 < ncol; kk2++) {
	//		for (int kk3 = 0; kk3 < top; kk3++) {
	//			cout << setw(20) << correlation_top[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//	cout << "final correlation_ranking_top matrix at node " << mynode << " where top is " << top << endl;
	//	for (int kk2 = 0; kk2 < ncol; kk2++) {
	//		for (int kk3 = 0; kk3 < top; kk3++) {
	//			cout << setw(20) << correlation_ranking_top[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//}

	//open file string to read uox and mox
	MPI_File fh_uox_reading;
	MPI_File fh_mox_reading;

	success = MPI_File_open(MPI_COMM_WORLD, "./uox_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_uox_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of uox fail to open the file!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./mox_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_mox_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of mox fail to open the file!" << endl;


	List_nU.initialize(nrow_mox);
	std::vector<int> v_nD;
	int i_cellmake = 2; // Inactiavte the b_success_nDAU because cell make with KNN always have enough donors
	bool b_success_nDAU = 0;//0 is False and 1 is True

	b_success_nDAU = nDAU_Bigp_ultra_cpp_MPI(memory, nrow_uox, nrow_mox, nrow, ncol, i_collapsing, i_SIS_type,

		top, correlation_ranking_top, correlation_top, i_cellmake,

		fh_uox_reading, fh_mox_reading, uox_infor,

		i_merge, d_k,

		v_nD, List_nU, codes, TestOut);

	success = MPI_File_close(&fh_uox_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of uox fail to close the file!" << endl;

	success = MPI_File_close(&fh_mox_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of mox fail to close the file!" << endl;


	if (!b_success_nDAU)
	{
		cout << "Error! nDAU Failed! Change k, check data quality, further break down categorical variables, or so. It may help " << endl;

		//early deallocaiton -----------------
		Del_iMatrix(correlation_ranking_top, ncol, top);
		Del_dMatrix(correlation_top, ncol, top);

		return 0; //abnormal ending 	

		//exit(0);
	}

	if (mynode == 0) {
		double v_nD_max = 0;
		for (int t = 0; t < v_nD.size(); t++) {
			v_nD_max = v_nD_max + (double)v_nD[t];
		}

		double required_memory = 0.0;
		required_memory = (v_nD_max * 8) / pow(10.0, 9.0);

		if ((0.5*memory) < required_memory) {
			cout << "CAUTION! The required memory to store list of all possible donors of all mox is " << required_memory << ", such that there may be memory issue in estimation of cell probability! " << endl;
		}
	}

	//if (mynode == 1) {
	//	for (int j = 0; j < nrow_mox; j++) {
	//		cout << "final v_nD[" << j << "]: " << v_nD[j] << " at node " << mynode << endl;
	//	}

	//	cout << "final List_nU at node " << mynode << endl;
	//	List_nU.print_List_FHDI_yicheng();

	//}

	//-----------------------------
	//VALIDATION 
	//------------------------------

	//Validation_cell_make_bigp_knn(nrow_uox, nrow_mox, nrow, ncol, max_overlap_size, i_collapsing,

	//	fh_datz, fh_uox, fh_mox, uox_info_final, mox_info_final, codes,

	//	v_nD, List_nU, TestOut);


	//if (mynode == 0) {
	//	cout<<"Final z matrix at node "<<mynode<<endl;

	//	for (int kk2 = 0; kk2 < 6; kk2++) {
	//		for (int kk3 = 0; kk3 < ncol; kk3++) {
	//			cout << setw(20) << matrix_temp[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//}
	////
	//////RPrint(matrix_temp, nrow, 4, TestOut);

	//Del_dMatrix(matrix_temp, 6, ncol);

	//delete[] array_temp;

	//Del_iMatrix(codes, nrow, i_collapsing);

	Del_iMatrix(correlation_ranking_top, ncol, top);

	Del_dMatrix(correlation_top, ncol, top);

	return 1;

}
