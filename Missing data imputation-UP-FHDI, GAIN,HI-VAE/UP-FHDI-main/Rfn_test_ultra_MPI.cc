//----------------------------------------
//---------------------------------------- 
#include "Cell_Make_Neighbor_ultra_cpp.cc"
#include "Cell_Make_Neighbor_Bigp_ultra_cpp.cc"
#include "Cell_Prob_Neighbor_ultra_new_cpp_MPI.cc"
#include "Cell_Prob_Neighbor_ultra_cpp_MPI.cc"
#include "FHDI_Neighbor_ultra_cpp_MPI.cc"
#include "Rep_CellP_Neighbor_ultra_MPI.cc"
#include "Variance_Est_FHDI_ultra_cpp.cc"
#include "Extract_Imputed_Results_Ultra.cc"
#include "Variance_Est_FHDI_Linerization_cpp.cc"

using std::cout;
using std::cerr;
using std::endl;


void Rfn_test_ultra_MPI(const int nrow_x, const int ncol_x, double* k,
	double* d_w, const int i_M,
	const int i_option_imputation, const int i_option_variance, const int i_option_collapsing, const int i_option_SIS_type, int top, const int i_option_cellmake, const int i_option_var_type,
	int* id,
	int* NonCollapsible_categorical,
	const int i_option_perform,
	const int i_option_merge, const int memory,
	ofstream& TestOut, ofstream& TestOut_Slave1, ofstream& TestOut_Slave2, ofstream& TestOut_Slave3, double d_begin_MPI)

	//Description================================================================
	// perform FEFI and FHDI based on the theory of Drs. Im, Kim and Fuller
	//
	// Last modified date: September 16, 2021
	// R code written by Dr. Im, J. H. and Dr. Kim, J. G. 
	// C++ code by Yicheng Yang and Dr. I. Cho
	// All rights reserved
	//---
	//IN	: double x(nrow_x, ncol_x) = {y1, y2, ...} raw data containing missing values
	//        Note: as of Oct 5, 2016. missing value is marked by a long number not NA 
	//				to avoid interface error between C++ and R
	//        Note: as of Feb 10, 2017
	//              Cell Prob Only (option =3) uses x as the categorized value matrix 
	// 
	//IN    : double x_raw(nrow_x, ncol_x) = original data matrix containing missing cells
	//IN    : int    r_raw(nrow_x, ncol_x) = index for missing (0) or observed (1)
	//IN	: double k(ncol_x)	= number of total categories per column of x
	//IN	: double d_w(nrow_x) 	= weight of row (default = 1.0)
	//IN    : int    M_donor = number of donors for FHDI
	//IN    : int    i_option_imputation = 1: FEFE; 2:FHDI
	//IN    : int    i_option_variance  = 0: skip variance estimation; 1: perform var. est. 
	//IN    : int    id(nrow_x) = ID of raw data 
	//IN    : double z(nrow_x, ncol_x) = user-defined category matrix (i_option_perform=4 only)
	//
	//IN    : int NonCollapsible_categorical(ncol) = {0,0, .., 1,.. 0} 
	//				index for non-collapsible categorical variables. 
	//				when at least one column has "1" skip cell-collapse procedure
	//				this may casue a potential error of lack of enough donor! 
	//				(2018, 04 21) 		
	//
	//IN    : int i_option_perform = main performance option
	//                            1=perform entire FEFI/FHDI (cellmake, jp, impute, var)
	//                            2=only Cell Make  
	//                            3=cell prob only
	//                            4=perform all FEFI/FHDI but by using User-Defined Z matrix 
	//IN    : int i_option_merge = random donor selection in Merge algorithm in Cell Make
	//                          0= no random seed number setting
	//						    1= random seed number setting 
	//IN    : int i_option_collapsing = number of selected most correlated variables
	//                               0= activate big-n algorithm
	//                              !0= activate big-p algorithm
	//IN    : int i_option_var_type = type of variance estimation
	//                                  1= parallel Jackknife variance estimation
	//                                  2= parallel variance estimation using linerization
	//IN    : int memory            = memory for each core on your platform. Will be used for out of memory check
	//============================================================================

	//Input file notes===============================================
	//daty: daty_column_binary.bin in column-oriented distribution
	//      daty_row_binary.bin in row-oriented distribution. Note that this file is only used in imputation function 
	//      to avoid expensive non-contiguous reading of daty
	//
	//datr: datr_column_binary.bin in column-oriented distribution
	//===============================================================

	//Output file notes===========================================================================================
	//Note that UP-FHDI has non-contiguously writing for z matrix in category function.
	//
	//(1) z matrix: datz_binary.bin in row-oriented distribution. Write non-contiguously in category and read contiguously in ZMAT
	//
	//(2) uox and mox: uox_binary.bin and mox_binary.bin in row-oriented distribution
	//
	//(3) imputed values with fractional weights: fmat_FHDI_binary.bin in row-oriented distribution
	//
	//(4) final imputed values: final_daty_binary.bin in row-oriented distribution
	//=============================================================================================================
{
	//-- MPI variables
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;


	const int nrow = nrow_x;
	const int ncol = ncol_x;
	const int i_imputation = i_option_imputation;
	const int i_variance = i_option_variance;
	const int i_merge = i_option_merge;
	const int i_collapsing = i_option_collapsing;
	const int i_SIS_type = i_option_SIS_type;
	const int i_cellmake = i_option_cellmake;
	const int i_var_type = i_option_var_type;


	//One need recursively read data in zmat and nDAU functions
	//To determine the number of rows to read recursively by empirical equation, 
	//one has to set up memory capacity of a processor in the used platform (e.g., 8GB on Condo)
	//const int memory = 8;

	//-----------------------------------------------
	//open file strings to read daty and datr
	//open file string to write final imputaed values
	//------------------------------------------------

	MPI_File fh_binary_daty;// File to read daty column-wisely
	MPI_File fh_binary_daty_row;// File to read daty row-wisely

	MPI_File fh_binary_datr;// File to read datr
	MPI_File fh_final_binary_daty;// File to write final daty

	int success = 0;// inidicator of accessing file system

	success = MPI_File_open(MPI_COMM_WORLD, "./daty_column_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_binary_daty);
	if (success != MPI_SUCCESS) cout << "MPI I/O fail to open the file!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./daty_row_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_binary_daty_row);
	if (success != MPI_SUCCESS) cout << "MPI I/O fail to open the file!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./datr_column_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_binary_datr);
	if (success != MPI_SUCCESS) cout << "MPI I/O fail to open the file!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./final_daty_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_final_binary_daty);
	if (success != MPI_SUCCESS) cout << "MPI I/O of uox fail to open the file!" << endl;

	//====================================================
	//====================================================
	//Cell_Make(): make imputation cells and donor list
	//======================================================
	//=====================================================
	double cell_make_begin = MPI_Wtime();

	bool b_Cell_Make_Neighbor_ultra = 0;//0 is False and 1 is True
	int nrow_uox = 0;//number of uox
	int nrow_mox = 0;//number of mox
	List_FHDI uox_infor(nrow); // Actual index list of uox in z. size will be updated to nrow_uox inside
	List_FHDI mox_infor(nrow); // Actual index list of mox in z. size will be updated to nrow_mox inside
	List_FHDI List_nU(nrow); //default for the size of nrow, but will be updated in the main loop

	if ((i_option_perform != 4) && (i_cellmake == 2) && (i_collapsing == 0)) {


		b_Cell_Make_Neighbor_ultra = Cell_Make_Neighbor_ultra_cpp(nrow, ncol, memory, k, NonCollapsible_categorical,

			fh_binary_daty, fh_binary_datr,

			uox_infor, mox_infor, nrow_uox, nrow_mox,

			List_nU, i_merge, TestOut);

		//if (mynode == 0) {
		//	cout << " ========= Cell_Make_Neighbor for ultra data.. has successfully finished!" << endl;
		//}

	}

	int** codes = New_iMatrix(nrow, i_collapsing);// The record of i_option_collapsing most correlated variables of each mox

	if ((i_option_perform != 4) && (i_cellmake == 2) && (i_collapsing != 0)) {

		b_Cell_Make_Neighbor_ultra = Cell_Make_Neighbor_Bigp_ultra_cpp(nrow, ncol, memory, k, NonCollapsible_categorical,

			fh_binary_daty, fh_binary_datr,

			uox_infor, mox_infor, nrow_uox, nrow_mox, codes,

			List_nU, i_merge, i_collapsing, i_SIS_type, top,

			TestOut);

		//if (mynode == 0) {
		//	cout << " ========= Cell_Make_Neighbor_Bigp for ultra data.. has successfully finished!" << endl;
		//}

	}

	Del_iMatrix(codes, nrow, i_collapsing);//Keep codes in case we need to check selected variables for all mox

	if (!b_Cell_Make_Neighbor_ultra) {

		success = MPI_File_close(&fh_binary_daty);
		if (success != MPI_SUCCESS) cout << "MPI I/O of daty fail to close the file!" << endl;

		success = MPI_File_close(&fh_binary_daty_row);
		if (success != MPI_SUCCESS) cout << "MPI I/O of daty fail to close the file!" << endl;

		success = MPI_File_close(&fh_binary_datr);
		if (success != MPI_SUCCESS) cout << "MPI I/O of datr fail to close the file!" << endl;

		success = MPI_File_close(&fh_final_binary_daty);
		if (success != MPI_SUCCESS) cout << "MPI I/O of imputed matrix fail to close the file!" << endl;

		return;
	}

	//cout<<"nrow_uox is "<< nrow_uox <<" and nrow_mox is "<< nrow_mox <<" at node "<<mynode<<endl;

	MPI_Barrier(MPI_COMM_WORLD);

	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_make at node " << mynode << " = " << MPI_Wtime() - cell_make_begin << endl;
	//	printf("Yicheng Running time of Cell_make for ultra data = %f seconds\n", MPI_Wtime() - cell_make_begin);
	//}

	//if (mynode == 1) {
	//	cout << "nrow_uox is " << nrow_uox << " and nrow_mox is " << nrow_mox << " at node " << mynode << endl;
	//	
	//	cout << "Final uox_infor from node " << mynode << endl;
	//	uox_infor.print_List_FHDI_yicheng();

	//	cout << "Final mox_infor from node " << mynode << endl;
	//	mox_infor.print_List_FHDI_yicheng();

	//	cout << "Final List_nU at node " << mynode << endl;
	//	List_nU.print_List_FHDI_yicheng();
	//}


	//===================================
	//===================================
	//Cell_prob(): calculate the joint probability of cells
	//===================================
	//===================================
	double cell_prob_begin = MPI_Wtime();

	std::vector<double> jp_prob_return; //the latest joint probability 

	int ml_temp = 0; 
	int i_size_ml = 0;// number of missing units in z matrix
	for (int j = 0; j < nrow_mox; j++) {
		ml_temp = 0;
		mox_infor.get_a_row_size(j, ml_temp);
		i_size_ml = i_size_ml + ml_temp;
	}

	//Note that this information will only be used in variance estimation using Lineriaztion
	//One should know imputation groups for each missing unit with this information
	double** mox_Imputation_group = New_dMatrix(i_size_ml, 2);//(gloabl index of mox in z matrix starting from 1) + (imputation group) 

	if (i_cellmake == 2) {

		Cell_Prob_Neighbor_ultra_new_cpp_MPI(nrow, ncol, List_nU, uox_infor, mox_infor,

			nrow_uox, nrow_mox, d_w, jp_prob_return, mox_Imputation_group, i_size_ml,

			TestOut);

		//if (mynode == 0) {
		//	cout << " ========= Cell_Prob_Neighbor for ultra data.. has successfully finished!" << endl;
		//}
	}


	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_Prob at node " << mynode << " = " << MPI_Wtime() - cell_prob_begin << endl;
	//	printf("Yicheng Running time of Cell_Prob for ultra data = %f seconds\n", MPI_Wtime() - cell_prob_begin);
	//}


	//=====================================
	//=====================================
	//FHDI: Fractional Hot Deck Imputation 
	//=====================================
	//=====================================
	double FHDI_imputation_begin = MPI_Wtime();

	//--------------------------------------------------------------------------------------
	//  To pass return matrix of imputed results (simp_fmat_FHDI) for later use, one have to 
	//   know its dimensionality first, especially number of rows
	//   note that it has 4 columns: global id + sampling weight + wij + fwij
	//---------------------------------------------------------------------------------

	int i_row_fmat = 0;

	//-----------------------------------------
	//Compute nrow of simp_fmat_FHDI of all mox in advance

	std::vector<int> loc_srst_nl;
	int uox_sum = 0;
	int temp_size = 0;

	int i_ol_size = 0;
	std::vector<int> ol; // actual locations of all fully observed rows in daty
	std::vector<double> ol_temp;// Actual location of observed patterns in z matrix (not sorted)

	uox_infor.unlist(ol_temp);
	int i_v_ol = 0;
	for (int k = 0; k < ol_temp.size(); k++) {
		i_v_ol = 0;
		i_v_ol = (int)ol_temp[k];
		ol.push_back(i_v_ol);
	}

	i_ol_size = ol.size(); // number of fully observed rows in daty

	//Accumulate number of rows of simp_fmat_FHDI
	int i_temp2 = 0;
	for (int k = 0; k < nrow_mox; k++) {
		loc_srst_nl.clear();
		uox_sum = 0;
		List_nU.get_block_yicheng(k, loc_srst_nl);
		temp_size = loc_srst_nl.size();
		for (int t = 0; t < temp_size; t++) {
			i_temp2 = 0;
			uox_infor.get_a_row_size(loc_srst_nl[t] - 1, i_temp2);
			uox_sum = uox_sum + i_temp2;
		}

		if (uox_sum > i_M) uox_sum = i_M;

		i_temp2 = 0;
		mox_infor.get_a_row_size(k, i_temp2);
		i_row_fmat = i_row_fmat + i_temp2 * uox_sum;
	}

	//Note that one must add fully observed rows to return matrix as well
	i_row_fmat = i_row_fmat + i_ol_size;

	//if (mynode == 1) { cout << "i_row_fmat is " << i_row_fmat << " and i_ol_size is " << i_ol_size << endl; }

	double** simp_fmat_FHDI = New_dMatrix(i_row_fmat, 4);// global id + sampling weight + wij + fwij

	std::string s_M;
	if ((i_imputation == 2) && (i_cellmake == 2)) {
		s_M = "FHDI";
		FHDI_Neighbor_ultra_cpp(fh_binary_daty_row, nrow, ncol, i_merge,
			nrow_uox, nrow_mox,
			jp_prob_return,
			s_M, i_M, List_nU, ol,
			uox_infor, mox_infor, simp_fmat_FHDI,
			TestOut_Slave1);

	}

	//if (mynode == 1) {
	//	cout << "Final simp_fmat_FHDI at node " << mynode << endl;

	//	for (int kk2 = 0; kk2 < i_row_fmat; kk2++) {
	//		for (int kk3 = 0; kk3 < 4; kk3++) {
	//			cout << setw(20) << simp_fmat_FHDI[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//}

	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_Prob at node " << mynode << " = " << MPI_Wtime() - cell_prob_begin << endl;
	//	printf("Yicheng Running time of FHDI Imputation for ultra data = %f seconds\n", MPI_Wtime() - FHDI_imputation_begin);
	//}

	//-----------------------------------------------
	//-- Extract imputed result
	//-------------------------------------------------
	double extract_begin = MPI_Wtime();

	Extract_Imputed_Results_Ultra(nrow, ncol, i_row_fmat,

		simp_fmat_FHDI, fh_final_binary_daty,

		TestOut_Slave3);

	//================================================
	//Close files
	//Note need to close earlier here becasue we need 
	//read imputed daty in Linerization
	//================================================
	success = MPI_File_close(&fh_final_binary_daty);
	if (success != MPI_SUCCESS) cout << "MPI I/O of imputed matrix fail to close the file!" << endl;


	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_Prob at node " << mynode << " = " << MPI_Wtime() - cell_prob_begin << endl;
	//	printf("Yicheng Running time of Extract imputed values for ultra data = %f seconds\n", MPI_Wtime() - extract_begin);
	//}


	//=========================================================
	//=========================================================
	//Variance: compute joint probability of uox based on 
	//          different replicate weight
	//========================================================
	//========================================================
	double Rep_CellP_begin = MPI_Wtime();

	//compute replicate weights
	RepWeight_FHDI d_rw(nrow);
	std::vector<int> List_rst_prob_size; // very important used in FHDI, the memory of List_rst_prob and List_rst_name will be tremendous if input data is big large. 

	//Rep_cellP is commented in UP-FHDI for now since it will not affect any results

	//if ((i_variance == 1) && (i_imputation == 2) && (i_cellmake == 2) && (i_var_type == 1)) {

	//	Rep_CellP_Neighbor_ultra(nrow, ncol,

	//		uox_infor, mox_infor, nrow_uox, nrow_mox,

	//		d_rw, List_nU, List_rst_prob_size,

	//		TestOut);

	//	//To simplify problem, we don't import it to variance estimation part
	//	//If any size of joint probability is 0, then there is some potential problem

	//	if (List_rst_prob_size.size() == 0) {
	//		cout << "ERROR in Rep_CellP for ulra data, no size of joint probability is added!!!" << endl;
	//		return;
	//	}

	//	for (int k = 0; k < List_rst_prob_size.size(); k++) {
	//		if (List_rst_prob_size[k] == 0) {
	//			cout << "ERROR in Rep_CellP for ulra data, the size of joint probability based on " << k << "th replcate weights are 0!!!" << endl;
	//		}
	//	}

	//}

	//MPI_Barrier(MPI_COMM_WORLD);

	//if ( (mynode == 1) && (i_var_type == 1) ) {
	//	int i_size = 0;
	//	i_size = (int)List_rst_prob_size.size();
	//	cout << "Final List_rst_prob_size is " << i_size << " at node " << mynode << endl;
	//	for (int k = 0; k < i_size; k++) {
	//		cout << "List_rst_prob_size[" << k << "]: " << List_rst_prob_size[k] << endl;
	//	}
	//}

	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_Prob at node " << mynode << " = " << MPI_Wtime() - cell_prob_begin << endl;
	//	printf("Yicheng Running time of Rep_CellP for ultra data = %f seconds\n", MPI_Wtime() - Rep_CellP_begin);
	//}

	//======================================================
	//======================================================
	//compute Jackknife variance 
	//======================================================
	//======================================================
	double variance_begin = MPI_Wtime();

	if ((i_variance == 1) && (i_imputation == 2) && (i_cellmake == 2) && (i_var_type == 1)) {

		Variance_Est_FHDI_ultra_cpp(fh_binary_daty, nrow, ncol,
			d_rw, id, ol,
			simp_fmat_FHDI, i_row_fmat,
			List_nU, uox_infor, mox_infor, nrow_mox,
			TestOut_Slave3);

	}

	//======================================================
	//======================================================
	//compute Linerization variance 
	//======================================================
	//======================================================

	if ((i_variance == 1) && (i_imputation == 2) && (i_cellmake == 2) && (i_var_type == 2)) {

		Variance_Est_FHDI_Linerization_cpp(fh_binary_daty, fh_binary_datr, nrow, ncol, memory, k,
			mox_infor, nrow_mox, mox_Imputation_group, i_size_ml,
			TestOut_Slave3);

	}

	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_Prob at node " << mynode << " = " << MPI_Wtime() - cell_prob_begin << endl;
	//	printf("Yicheng Running time of Variance for ultra data = %f seconds\n", MPI_Wtime() - variance_begin);
	//}



	//====================
	//Final Deallocation
	//=====================
	Del_dMatrix(simp_fmat_FHDI, i_row_fmat, 4);
	Del_dMatrix(mox_Imputation_group, i_size_ml, 2);

	//================================================
	//Close all files
	//================================================
	success = MPI_File_close(&fh_binary_daty);
	if (success != MPI_SUCCESS) cout << "MPI I/O of daty fail to close the file!" << endl;

	success = MPI_File_close(&fh_binary_daty_row);
	if (success != MPI_SUCCESS) cout << "MPI I/O of daty fail to close the file!" << endl;

	success = MPI_File_close(&fh_binary_datr);
	if (success != MPI_SUCCESS) cout << "MPI I/O of datr fail to close the file!" << endl;



	return; //test return of double* 

}



