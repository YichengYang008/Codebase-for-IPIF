//---------------------------------------
//%Development note for MPI version
//%2018 1220
//
// new version of 
//
// categorize_cpp.cc
// Cell_Make_Extension_cpp.cc
// Rfn_test_MPI.cc
// main_MPI.cpp
// nDAU_cpp.cc
//----------------------------------------
//---------------------------------------- 
#include <limits> // For NaN 
#include <iostream> // For TestOut, cerr, endl, and flush 
#include <assert.h> // For assert
#include <algorithm> // For sort 
#include <string>		//For string array
#include <vector>
#include <iomanip>	// setw
//#include <map>

//#include "base_FHDI.cc"
//#include "List_FHDI.cc"
//#include "List_string_FHDI.cc"
//#include "rbind_FHDI.cc"

//#include "categorize_cpp.cc"
#include "categorize_cpp_MPI.cc"
#include "Zmat_Extension_cpp.cc"  
//#include "Zmat_Extension_cpp_MPI.cc"  
//#include "nDAU_cpp.cc"
#include "nDAU_cpp_MPI.cc"
#include "nDAU_Bigp_cpp_MPI.cc"

#include "Cell_Make_Extension_cpp.cc"
#include "Cell_Make_Extension_Bigp_cpp.cc"
#include "Cell_Make_Neighbor_cpp.cc"
#include "Cell_Make_Neighbor_Bigp_cpp.cc"

#include "Cell_Prob_Extension_cpp_MPI.cc"
#include "Cell_Prob_Neighbor_cpp_MPI.cc"
#include "Cell_Prob_Extension_Bigp_cpp_MPI.cc"

#include "FHDI_Extension_cpp_MPI.cc"
#include "FHDI_Neighbor_cpp_MPI.cc"
#include "FHDI_Extension_Bigp_cpp_MPI.cc"

#include "Extract_Imputed_Results.cc"
#include "matrix_utility_FHDI.cc"
#include "Variance_Est_FEFI_Extension_cpp.cc"
#include "Extract_Variance_Results_MPI.cc"

#include "Rep_CellP_MPI_FEFI.cc"
#include "Rep_CellP_MPI_FHDI.cc"
#include "Rep_CellP_Bigp_MPI_FEFI.cc"
#include "Rep_CellP_Bigp_MPI_FHDI.cc"

#include "Rep_CellP_Neighbor_MPI_FEFI.cc"
#include "Rep_CellP_Neighbor_MPI_FHDI.cc"

#include "Variance_Est_FHDI_Extension_cpp.cc"
//#include "Variance_Est_FHDI_Extension_Bigp_cpp.cc"


#define i_DEBUGGING 0 //0=no printout; 1=printout for debugging 

using std::cout; 
using std::cerr; 
using std::endl; 


void Rfn_test_MPI(double** x_raw, int** r_raw, const int nrow_x, const int ncol_x, double* k, 
                   double* d_w, const int i_M,
				   const int i_option_imputation, const int i_option_variance, const int i_option_collapsing, const int i_option_SIS_type, int top, const int i_option_cellmake,
				   int* id, double** z, 
				   int* NonCollapsible_categorical,
				   rbind_FHDI &rbind_ipmat_FEFI,
				   rbind_FHDI &rbind_Resp_FEFI, 
				   rbind_FHDI &rbind_irmat_FEFI,
				   rbind_FHDI &rbind_ipmat_FHDI,
				   rbind_FHDI &rbind_Resp_FHDI, 
				   rbind_FHDI &rbind_irmat_FHDI,
				   rbind_FHDI &rbind_vrst_FEFI, 
				   rbind_FHDI &rbind_vrst_FHDI,

				   rbind_FHDI  &rbind_uox_CellMake,
				   rbind_FHDI  &rbind_mox_CellMake,	
				   rbind_FHDI  &rbind_category_CellMake,	
				   
				   std::vector<std::string> &jp_name_return_CellProb,				   
				   std::vector<double> &jp_prob_return_CellProb,				
				   
				   const int i_option_perform, 
				   const int i_option_merge,
				   ofstream& TestOut, ofstream& TestOut_Slave1, ofstream& TestOut_Slave2, ofstream& TestOut_Slave3, double d_begin_MPI)
				   
//Description================================================================
// perform FEFI and FHDI based on the theory of Drs. Im, Kim and Fuller
//
// Cell_make
// Cell_Prob
//
// October 5, 2016
// R code written by Dr. Im, J. H. and Dr. Kim, J. G. 
// C++ code by Dr. I. Cho
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
//
//OUT   : rbind_FHDI  rbind_ipmat_FEFI(4+ncol) //column size is 4+ncol (i.e., for R: ID, FID, WGT, FWGT, Variables)
//OUT   : rbind_FHDI  rbind_Resp_FEFI(ncol+1)  //separate response matrix  (i.e. for R: unit responses and Resp0)
//OUT   : rbind_FHDI  rbind_irmat_FEFI(5+ncol) //column size is 5+ncol (i.e. for R:ID, FID, OID, ORDER, FEFIW, CELL )
//OUT   : rbind_FHDI  rbind_ipmat_FHDI(4+ncol) //column size is 4+ncol (i.e., for R: ID, FID, WGT, FWGT, Variables)
//OUT   : rbind_FHDI  rbind_Resp_FHDI(ncol+1)  //separate response matrix  (i.e. for R: unit responses and Resp0)
//OUT   : rbind_FHDI  rbind_irmat_FHDI(5+ncol) //column size is 5+ncol (i.e. for R:ID, FID, OID, ORDER, FEFIW, CELL )
//OUT   : rbind_FHDI  rbind_vrst_FEFI(nrow)    //variance estimates of FEFI
//OUT   : rbind_FHDI  rbind_vrst_FHDI(nrow)    //variance estimates of FHDI
//
//OUT   : rbind_FHDI  rbind_uox_CellMake(ncol) //unique patterns of observed rows
//OUT   : rbind_FHDI  rbind_mox_CellMake(ncol) //unique patterns of missing rows
//OUT   : rbind_FHDI  rbind_category_CellMake(ncol) //matrix of categorized values
//
//OUT   : std::vector<double> &jp_prob_return_CellProb  //joint prob for Cell Prob option
//OUT   : std::vector<std::string> &jp_name_return_CellProb //name of jp for Cell Prob option 			
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
//============================================================================
{
	//-- MPI variables
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	const int nrow = nrow_x;
	const int ncol = ncol_x; 
	//int *i_M = new int[ncol]; Copy_iVector(M_donor, ncol, i_M); //defined at main()
	const int i_imputation = i_option_imputation;
	const int i_variance  = i_option_variance;
	const int i_merge = i_option_merge; 
	const int i_collapsing = i_option_collapsing;
	const int i_SIS_type = i_option_SIS_type;
	const int i_cellmake = i_option_cellmake;
	//RPrint("Imputation type (1=FEFI; 2=FHDI):  "); RPrint(i_imputation);
	//TestOut<<"i_variance:  "<<i_variance<<endl;
	//RPrint("variance estimation (0=skip; 1=Jackknife):  "); RPrint(i_variance);
	
	//if(i_DEBUGGING)
	//{
	//Very impotant initial TestOut 
	//TestOut << "======== Begin Rfn_test =====================================" << endl;
	//system("PAUSE"); 

		//TestOut<<"nrow:  "<<nrow<<endl;
		//TestOut<<"ncol:  "<<ncol<<endl;
		//TestOut << "M = " <<i_M<< endl;

		//TestOut<<"x_raw[]"<<endl;
		//for(int i=0; i<nrow; i++)
		//{
		//	for(int j=0; j<ncol; j++){TestOut<<setw(20)<<x_raw[i][j];}
		//	TestOut<<endl;
		//}
		//TestOut<<"r_raw[]"<<endl;
		//for(int i=0; i<nrow; i++)
		//{
		//	for(int j=0; j<ncol; j++){TestOut<<setw(10)<<r_raw[i][j];}
		//	TestOut<<endl;
		//}		
		//TestOut<<"k[]"<<endl;
		//for(int i=0; i<ncol; i++){TestOut<<setw(10)<<k[i];}
		//TestOut << endl;
		//TestOut<<"d_w[]"<<endl;
		//for(int i=0; i<nrow; i++){TestOut<<setw(10)<<d_w[i];}
		//TestOut << endl;
		//TestOut<<"id[]"<<endl;
		//for(int i=0; i<nrow; i++){TestOut<<setw(10)<<id[i];}	
		//TestOut << endl;
	//}	
 	
    

	//defined at main()
	//----------------------------
	//get matrix format of x -> x_raw containing all data and missing cells  
	//based on row-first rule of R 
	//----------------------------
	//double** x_raw = New_dMatrix(nrow,ncol); ////defined at main()
	//for(int i=0; i<ncol; i++){ for(int j=0; j<nrow; j++) x_raw[j][i] = x[j+i*nrow]; } 

	
	
	//-------------
	//Special Individual Running for 
	//Cell Prob ONLY!
	//-------------
	if(i_option_perform ==3) 
	{
		//z is already stored in x_raw in this special option 
		//double** z_CellProb = New_dMatrix(nrow, ncol); //initialized by 0.0

		//===================================
		//ONLY Cell_prob(): calculate the joint probability of cells
		//===================================
		//defined outside 
		//std::vector<std::string> jp_name_return_CellProb;   //name of the joint probability table
		//std::vector<double> jp_prob_return_CellProb; //the latest joint probability 
	
		Cell_Prob_Extension_cpp_MPI(z, nrow, ncol, 
		                          jp_prob_return_CellProb, jp_name_return_CellProb, 
	                              d_w, id, TestOut);
		//testout 
		//RPrint("jp_prob_return_CellProb: ", TestOut); 
		//RPrint(jp_prob_return_CellProb, TestOut);
		//RPrint("jp_name_return_CellProb: ", TestOut);
		//RPrint(jp_name_return_CellProb, TestOut);
			
		//testout
		//RPrint(" ========= Cell_Prob_ONLY.. has successfully finished!", TestOut);

		//=========================
		//Deallocate memories
		//=========================
		//Del_dMatrix(x_raw, nrow, ncol); //defined at main()
		
		return; 		
	}
	
	//defined at main()
	//----------------------------
	//get matrix format of r -> r_raw containing index for missing/observed
	//based on row-first rule of R 
	//----------------------------
	//int** r_raw = New_iMatrix(nrow, ncol); //initialized by 0
	//for(int i=0; i<ncol; i++){ for(int j=0; j<nrow; j++) r_raw[j][i] = r[j+i*nrow]; } 


	//=====================================
	//=====================================
	//Cell_Make task
	//=====================================
	//=====================================
	//defined at main
	//Storages to be returned
	//double** z = New_dMatrix(nrow, ncol); //initialized by 0.0

	//--------------------------------------------------
	//--------------------------------------------------
	//CELL MAKE-----------------------------------------
	//--------------------------------------------------
	//--------------------------------------------------

	int** codes = New_iMatrix(nrow, i_collapsing);// The record of i_option_collapsing most correlated variables of each mox

	//if i_option_perform = 4, user-defined z option, skip Cell_Make
	//April 4, 2017 =======================
	double cell_make_begin = MPI_Wtime();

	if ( (i_option_perform != 4) && (i_cellmake == 1) && (i_collapsing == 0) )
	{
		//cout<<"Begin of Cell_make at node "<<mynode<<endl;
		bool b_Cell_Make = Cell_Make_Extension_cpp(x_raw, nrow, ncol, k,  
								NonCollapsible_categorical,
								z,
								rbind_uox_CellMake, rbind_mox_CellMake, 
								i_merge, TestOut, d_begin_MPI);
		//cout << "Running time after Cell_Make_Extension_cpp at node " << mynode << " = " << MPI_Wtime() - d_begin_MPI << endl;
		
		//if(mynode==0){
		//	cout << "YYC Running time of Cell_make at node " << mynode << " = " << MPI_Wtime() - cell_make_begin << endl;
		//	printf("Yicheng Running time of Cell_make = %f seconds\n", MPI_Wtime() - cell_make_begin);
		//}
		
		//testout
		//if (mynode == 0) {
		//	TestOut << " z matrix from Cell Make with cell collapsing" << endl;
		//	RPrint(z, nrow, ncol, TestOut);
		//}
		//if (mynode == 0) {
		//	cout << " ========= Cell_Make_Extension.. has successfully finished!" << endl;
		//}

	}

	if ( (i_option_perform != 4) && (i_cellmake == 1) && (i_collapsing != 0) ) {

		//double cell_make_begin = MPI_Wtime();
		bool b_Cell_Make = Cell_Make_Extension_Bigp_cpp(x_raw, r_raw, nrow, ncol, k,
			NonCollapsible_categorical,
			z, codes,
			rbind_uox_CellMake, rbind_mox_CellMake,
			i_merge, i_collapsing, i_SIS_type, top, TestOut, d_begin_MPI);

		//testout
		//if (mynode == 0) {
		//	TestOut << " z matrix from Cell Make" << endl;
		//	//RPrint(z, nrow, ncol, TestOut);
		//}
		//if (mynode == 0) {
		//	cout << " ========= Cell_Make_Extension_Bigp.. has successfully finished!" << endl;
		//}
	}

	List_FHDI List_nU(nrow); //default for the size of nrow, but will be updated in the main loop

	if ((i_option_perform != 4) && (i_cellmake == 2) && (i_collapsing == 0))
	{

		bool b_Cell_Make = Cell_Make_Neighbor_cpp(x_raw, nrow, ncol, k,
			NonCollapsible_categorical,
			z,
			rbind_uox_CellMake, rbind_mox_CellMake, List_nU,
			i_merge, TestOut);
		//testout

		//TestOut << "List_nU after cell make neighbor" << endl;
		//List_nU.print_List_FHDI_yicheng(TestOut);

		//TestOut << " z matrix from Cell Make with KNN" << endl;
		//RPrint(z, nrow, ncol, TestOut);

		//cout << " ========= Cell_Make.. has successfully finished!" << endl;
	}

	if ((i_option_perform != 4) && (i_cellmake == 2) && (i_collapsing != 0))
	{

		bool b_Cell_Make = Cell_Make_Neighbor_Bigp_cpp(x_raw, r_raw, nrow, ncol, k,
			NonCollapsible_categorical,
			z, codes,
			rbind_uox_CellMake, rbind_mox_CellMake, List_nU,
			i_merge, i_collapsing, i_SIS_type, top, TestOut);
		//testout

		//TestOut << "List_nU after cell make bigp neighbor" << endl;
		//List_nU.print_List_FHDI_yicheng(TestOut);

		//TestOut << " z matrix from Cell Make Bigp with KNN" << endl;
		//RPrint(z, nrow, ncol, TestOut);

		//cout << " ========= Cell_Make.. has successfully finished!" << endl;
	}

	//if(mynode==0){
	//	//cout << "YYC Running time of Cell_make at node " << mynode << " = " << MPI_Wtime() - cell_make_begin << endl;
	//	printf("Yicheng Running time of Cell_make = %f seconds\n", MPI_Wtime() - cell_make_begin);
	//}

	//TestOut << "Final Codes" << endl;

	//for (int kk2 = 0; kk2 < nrow; kk2++) {
	//	for (int kk3 = 0; kk3 < i_collapsing; kk3++) {
	//		TestOut << setw(20) << codes[kk2][kk3];
	//	}
	//	TestOut << endl;
	//}

	//--------------------------------
	//User-Defined z matrix: Override z of Cell Make 
	//--------------------------------
	if(i_option_perform == 4) //full FEFI/FHDI using user-defined z 
	{
		TestOut<<"Note: Category z matrix is using the user-defined matrix"<<endl;
		//defined at main()
		//for(int i=0; i<ncol; i++){ for(int j=0; j<nrow; j++) z[j][i] = z_UserDefined[j+i*nrow]; }
		//testout
		//RPrint(z, nrow, ncol); 
	}		

	
	if(i_option_perform == 2) //Cell Make only option; prep z to be out
	{
		double* d_temp_z = new double[ncol];
		for(int i=0; i<nrow; i++)
		{
			for(int j=0; j<ncol; j++) d_temp_z[j] = z[i][j]; 
			
			rbind_category_CellMake.append_block(d_temp_z); 
		}

		//=========================
		//Deallocate memories
		//=========================
		delete[] d_temp_z; 

		//Del_dMatrix(z, nrow, ncol); //defined at main
		//Del_dMatrix(x_raw, nrow, ncol); //defined at main()
		//Del_iMatrix(r_raw, nrow, ncol); //defined at main()
		
		return; 
	}
	

	//===================================
	//===================================
	//Cell_prob(): calculate the joint probability of cells
	//===================================
	//===================================
	std::vector<std::string> jp_name_return;   //name of the joint probability table
	std::vector<double> jp_prob_return; //the latest joint probability 
	double cell_prob_begin = MPI_Wtime();

	if ((i_cellmake == 1) && (i_collapsing == 0)) {

		Cell_Prob_Extension_cpp_MPI(z, nrow, ncol, jp_prob_return, jp_name_return,
			d_w, id, TestOut);

		//TestOut<<"jp_name_return:"<<endl;
		//for (int t = 0; t < jp_name_return.size();t++) {
		//	TestOut<<"jp_name_return["<<t<<"]: "<< jp_name_return[t]<<endl;
		//}
		//TestOut << "jp_prob_return:" << endl;
		//for (int t = 0; t < jp_prob_return.size();t++) {
		//	TestOut << "jp_prob_return[" << t << "]: " << jp_prob_return[t] << endl;
		//}

		//if (mynode == 0) {
		//	cout << " ========= Cell_Prob.. has successfully finished!" << endl;
		//}
	}

	if (i_cellmake == 2) {

		Cell_Prob_Neighbor_cpp_MPI(z, nrow, ncol, List_nU, jp_prob_return, jp_name_return,
			d_w, id, TestOut);
	}


	if ((i_cellmake == 1) && (i_collapsing != 0)) {
		Cell_Prob_Extension_Bigp_cpp_MPI(z, nrow, ncol, i_collapsing, jp_prob_return, jp_name_return,
			d_w, id, codes, TestOut);

		//if (mynode == 0) {
		//	cout << " ========= Cell_Prob_Bigp.. has successfully finished!" << endl;
		//}
	}

		//testout 

		//RPrint("after Cell_Prob..() =========Rand() is OFF and n: ");
		//RPrint(nrow);
		
		//RPrint("jp_prob_return: ", TestOut);
		//RPrint(jp_prob_return, TestOut);
		//RPrint("jp_name_return: ", TestOut);
		//RPrint(jp_name_return, TestOut);
		
	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_Prob at node " << mynode << " = " << MPI_Wtime() - cell_prob_begin << endl;
	//	printf("Yicheng Running time of Cell_Prob = %f seconds\n", MPI_Wtime() - cell_prob_begin);
	//}





	//=====================================
	//=====================================
	//FEFI: Fully Efficient Fractional Imputation
	//=====================================
	//Note: as of Nov 15, 2016. FHDI is following the FEFI without returning
	//       					this is computationally efficient
	//=====================================
	//FEFI---------returns----------------- FEFI //
	//get ipmat, Resp (separately), irmat from FEFI results 
	//rbind_FHDI  rbind_ipmat_FEFI(4+ncol); //column size is 4+ncol    //defined at Interface.cc 0117_2017
	//rbind_FHDI  rbind_Resp_FEFI(ncol+1);  //separate response matrix //defined at Interface.cc 0117_2017 
	//rbind_FHDI  rbind_irmat_FEFI(5+ncol); //column size is 5+ncol    //defined at Interface.cc 0117_2017
	//--------
	//FHDI --------returns----------------- FHDI //
	//--------
	//get ipmat, Resp (separately), irmat from FHDI results 
	//rbind_FHDI  rbind_ipmat_FHDI(4+ncol); //column size is 4+ncol	//defined at Interface.cc 0117_2017
	//rbind_FHDI  rbind_Resp_FHDI(ncol+1);  //separate response matrix//defined at Interface.cc 0117_2017  
	//rbind_FHDI  rbind_irmat_FHDI(5+ncol); //column size is 5+ncol	//defined at Interface.cc 0117_2017
	//---------
	//---------
	//other matrices
	//---------
	rbind_FHDI  rbind_uox(ncol); //observed unique categorized matrix 
	rbind_FHDI  rbind_mox(ncol); //missing  unique categorized matrix
	//Note: below Lists contain meaningful items up to i_count_mox rows  
	List_FHDI 	List_ord(nrow); //order records used for variance estimation
	List_FHDI 	List_ocsg(nrow); //order records used for variance estimation

	
	std::string s_M;
	//if (mynode == 0) {
	double FHDI_begin = MPI_Wtime();

	if ((i_imputation == 1) && (i_cellmake == 1) && (i_collapsing == 0)) //FEFI
		{

			s_M = "FEFI";
			FHDI_Extension_cpp(x_raw, z, r_raw,
				nrow, ncol,
				jp_name_return,
				jp_prob_return,
				s_M, i_M,

				rbind_ipmat_FEFI, rbind_Resp_FEFI, rbind_irmat_FEFI,
				rbind_ipmat_FHDI, rbind_Resp_FHDI, rbind_irmat_FHDI,
				rbind_uox, rbind_mox,
				List_ord, List_ocsg, TestOut_Slave1);

			//cout << " ========= FEFI_Extension.. has successfully finished!" << endl;
			//testout
			/*if(mynode==0){
				cout << "YYC Running time of FEFI at node " << mynode << " = " << MPI_Wtime() - FHDI_begin << endl;
				printf("Yicheng Running time of FEFI = %f seconds\n", MPI_Wtime() - FHDI_begin);

			}
			TestOut << "Yicheng Debug---- row of rbind_ipmat_FEFI is " << rbind_ipmat_FEFI.size_row() << endl;
			TestOut << "Yicheng Debug---- row of rbind_Resp_FEFI is " <<  rbind_Resp_FEFI.size_row() << endl;
			TestOut << "Yicheng Debug---- row of rbind_irmat_FEFI is " << rbind_irmat_FEFI.size_row() << endl;

			TestOut << "Yicheng Debug---- row of rbind_ipmat_FHDI is " << rbind_ipmat_FHDI.size_row() << endl;
			TestOut << "Yicheng Debug---- row of rbind_Resp_FHDI is " << rbind_Resp_FHDI.size_row() << endl;
			TestOut << "Yicheng Debug---- row of rbind_irmat_FHDI is " << rbind_irmat_FHDI.size_row() << endl;

			TestOut << "Yicheng Debug---- row of rbind_uox is " << rbind_uox.size_row() << endl;
			TestOut << "Yicheng Debug---- row of rbind_mox is " << rbind_mox.size_row() << endl;
			TestOut << "Yicheng Debug---- row of List_ord is " << List_ord.size_row() << endl;
			TestOut << "Yicheng Debug---- row of List_ocsg is " << List_ocsg.size_row() << endl;*/
			//TestOut << " ========= Main FEFI task has successfully finished!" << endl;
			//if (mynode == 0) {
			//	RPrint("Print out the rbind_ipmat_FEFI---------\n");
			//	rbind_ipmat_FEFI.print_rbind_FHDI();
			//}
		}

	if ((i_imputation == 1) && (i_cellmake == 1) && (i_collapsing != 0)) //FEFI
	{
		//double FHDI_begin = MPI_Wtime();
		s_M = "FEFI";
		FHDI_Extension_Bigp_cpp(x_raw, z, r_raw,
			nrow, ncol, i_collapsing,
			jp_name_return,
			jp_prob_return,
			s_M, i_M, codes,

			rbind_ipmat_FEFI, rbind_Resp_FEFI, rbind_irmat_FEFI,
			rbind_ipmat_FHDI, rbind_Resp_FHDI, rbind_irmat_FHDI,
			rbind_uox, rbind_mox,
			List_ord, List_ocsg, TestOut_Slave1);
		//if (mynode == 0) {
		//	cout << " ========= FEFI_Extension_Bigp.. has successfully finished!" << endl;
		//}
		//testout
		/*if(mynode==0){
		cout << "YYC Running time of FEFI at node " << mynode << " = " << MPI_Wtime() - FHDI_begin << endl;
		printf("Yicheng Running time of FEFI = %f seconds\n", MPI_Wtime() - FHDI_begin);

		}
		TestOut << "Yicheng Debug---- row of rbind_ipmat_FEFI is " << rbind_ipmat_FEFI.size_row() << endl;
		TestOut << "Yicheng Debug---- row of rbind_Resp_FEFI is " <<  rbind_Resp_FEFI.size_row() << endl;
		TestOut << "Yicheng Debug---- row of rbind_irmat_FEFI is " << rbind_irmat_FEFI.size_row() << endl;

		TestOut << "Yicheng Debug---- row of rbind_ipmat_FHDI is " << rbind_ipmat_FHDI.size_row() << endl;
		TestOut << "Yicheng Debug---- row of rbind_Resp_FHDI is " << rbind_Resp_FHDI.size_row() << endl;
		TestOut << "Yicheng Debug---- row of rbind_irmat_FHDI is " << rbind_irmat_FHDI.size_row() << endl;

		TestOut << "Yicheng Debug---- row of rbind_uox is " << rbind_uox.size_row() << endl;
		TestOut << "Yicheng Debug---- row of rbind_mox is " << rbind_mox.size_row() << endl;
		TestOut << "Yicheng Debug---- row of List_ord is " << List_ord.size_row() << endl;
		TestOut << "Yicheng Debug---- row of List_ocsg is " << List_ocsg.size_row() << endl;*/
		//TestOut << " ========= Main FEFI task has successfully finished!" << endl;
		//if (mynode == 0) {
		//	RPrint("Print out the rbind_ipmat_FEFI---------\n");
		//	rbind_ipmat_FEFI.print_rbind_FHDI();
		//}
	}

	if ((i_imputation == 2) && (i_cellmake == 1) && (i_collapsing == 0)) //FHDI
		{
			//double FHDI_begin = MPI_Wtime();
			s_M = "FHDI";
			FHDI_Extension_cpp(x_raw, z, r_raw,
				nrow, ncol,
				jp_name_return,
				jp_prob_return,
				s_M, i_M,

				rbind_ipmat_FEFI, rbind_Resp_FEFI, rbind_irmat_FEFI,
				rbind_ipmat_FHDI, rbind_Resp_FHDI, rbind_irmat_FHDI,
				rbind_uox, rbind_mox,
				List_ord, List_ocsg, TestOut_Slave1);
			//testout
			//if (mynode == 0) {
			//	cout << "YYC Running time of FHDI at node " << mynode << " = " << MPI_Wtime() - FHDI_begin << endl;
			//	printf("Yicheng Running time of FHDI = %f seconds\n", MPI_Wtime() - FHDI_begin);

			//}
			//if (mynode == 0) {
			//	cout << " ========= FHDI_Extension.. has successfully finished!" << endl;
			//}
		}

	if ((i_imputation == 2) && (i_cellmake == 1) && (i_collapsing != 0)) //FHDI
	{
		//double FHDI_begin = MPI_Wtime();
		s_M = "FHDI";
		FHDI_Extension_Bigp_cpp(x_raw, z, r_raw,
			nrow, ncol, i_collapsing,
			jp_name_return,
			jp_prob_return,
			s_M, i_M, codes,

			rbind_ipmat_FEFI, rbind_Resp_FEFI, rbind_irmat_FEFI,
			rbind_ipmat_FHDI, rbind_Resp_FHDI, rbind_irmat_FHDI,
			rbind_uox, rbind_mox,
			List_ord, List_ocsg, TestOut_Slave1);
		//testout
		//if (mynode == 0) {
		//	cout << "YYC Running time of FHDI at node " << mynode << " = " << MPI_Wtime() - FHDI_begin << endl;
		//	printf("Yicheng Running time of FHDI = %f seconds\n", MPI_Wtime() - FHDI_begin);

		//}
		//if (mynode == 0) {
		//	cout << " ========= FHDI_Extension_Bigp.. has successfully finished!" << endl;
		//}
	}
	//} // end of if (mynode == 0)

	

	//--------------------------------------------------------------------------------------
	//FEFI for KNN
	//--------------------------------------------------------------------------------------
	if ((i_imputation == 1) && (i_cellmake == 2)) //FEFI
	{
		s_M = "FEFI";
		    FHDI_Neighbor_cpp(x_raw, z, r_raw,
			nrow, ncol,
			jp_name_return,
			jp_prob_return,
			s_M, i_M, List_nU,

			rbind_ipmat_FEFI, rbind_Resp_FEFI, rbind_irmat_FEFI,
			rbind_ipmat_FHDI, rbind_Resp_FHDI, rbind_irmat_FHDI,
			rbind_uox, rbind_mox,
			List_ord, List_ocsg, TestOut);

		//testout
		//cout << "rbind_ipmat_FEFI with KNN:\n" << endl;
		//rbind_ipmat_FEFI.print_rbind_FHDI_Yicheng(TestOut);

		//cout << " ========= Main FEFI task with KNN (big-n or big-p) has successfully finished!" << endl;
		//cout << " ========= Main FEFI task has successfully finished!" << endl;
	}
	//--------------------------------------------------------------------------------------
	//FHDI for KNN
	//--------------------------------------------------------------------------------------

	if ((i_imputation == 2) && (i_cellmake == 2)) //FHDI
	{
		s_M = "FHDI";
		    FHDI_Neighbor_cpp(x_raw, z, r_raw,
			nrow, ncol,
			jp_name_return,
			jp_prob_return,
			s_M, i_M, List_nU,

			rbind_ipmat_FEFI, rbind_Resp_FEFI, rbind_irmat_FEFI,
			rbind_ipmat_FHDI, rbind_Resp_FHDI, rbind_irmat_FHDI,
			rbind_uox, rbind_mox,
			List_ord, List_ocsg, TestOut);
		//testout
		//cout << "rbind_ipmat_FHDI with KNN:\n" << endl;
		//rbind_ipmat_FHDI.print_rbind_FHDI_Yicheng(TestOut);

		//cout << " ========= Main FHDI task with KNN (big-n or big-p) has successfully finished!" << endl;
		//cout << " ========= Main FHDI task has successfully finished!" << endl;
	}

	//if (mynode == 0) {
	//	//cout << "YYC Running time of FHDI at node " << mynode << " = " << MPI_Wtime() - FHDI_begin << endl;
	//	printf("Yicheng Running time of FHDI = %f seconds\n", MPI_Wtime() - FHDI_begin);

	//}


		//TestOut_Slave1 << "Mynode at " << mynode << endl;
		//TestOut_Slave2 << "Mynode at " << mynode << endl;
		//TestOut_Slave3 << "Mynode at " << mynode << endl;
		//cout <<"Mynode at "<<mynode<<endl;
		if (i_imputation == 1) //FEFI
		{
			double* d_temp = new double[rbind_ipmat_FEFI.size_col()];
			TestOut_Slave2 << setw(20) << "ID" << setw(20) << "FID" << setw(20) << "WT" << setw(20) << "FWT";

			char numstr[21];
			for (int j = 0; j < ncol; j++) {
				std::string name = "y";
				std::string result;

				sprintf(numstr, "%d", j + 1);
				result = name + numstr;

				TestOut_Slave2 << setw(20) << result;
			}
			TestOut_Slave2 << endl;

			for (int i = 0; i < rbind_ipmat_FEFI.size_row(); i++) {
				rbind_ipmat_FEFI.get_block(i, d_temp);
				for (int j = 0; j < rbind_ipmat_FEFI.size_col(); j++) {
					TestOut_Slave2 << setw(20) << d_temp[j];
				}
				TestOut_Slave2 << endl;
			}

			delete[] d_temp;
		}

		if (i_imputation == 2) //FHDI
		{
			double* d_temp = new double[rbind_ipmat_FHDI.size_col()];
			TestOut_Slave2 << setw(20) << "ID" << setw(20) << "FID" << setw(20) << "WT" << setw(20) << "FWT";

			char numstr[21];
			for (int j = 0; j < ncol; j++) {
				std::string name = "y";
				std::string result;

				sprintf(numstr, "%d", j + 1);
				result = name + numstr;

				TestOut_Slave2 << setw(20) << result;
			}
			TestOut_Slave2 << endl;

			for (int i = 0; i < rbind_ipmat_FHDI.size_row(); i++) {
				rbind_ipmat_FHDI.get_block(i, d_temp);
				for (int j = 0; j < rbind_ipmat_FHDI.size_col(); j++) {
					TestOut_Slave2 << setw(20) << d_temp[j];
				}
				TestOut_Slave2 << endl;
			}

			delete[] d_temp;
		}
	//testout
	/*
	RPrint("after Results_... rbind_ipmat_FEFI after binding :");
	rbind_ipmat_FEFI.print_rbind_FHDI(); 
	RPrint("after Results_... rbind_Resp_FEFI after binding :");
	rbind_Resp_FEFI.print_rbind_FHDI(); 						
	RPrint("after Results_... rbind_irmat_FEFI after binding :");
	rbind_irmat_FEFI.print_rbind_FHDI(); 
	//testout
	RPrint("after Results_... rbind_ipmat_FHDI after binding :");
	rbind_ipmat_FHDI.print_rbind_FHDI(); 
	RPrint("after Results_... rbind_Resp_FHDI after binding :");
	rbind_Resp_FHDI.print_rbind_FHDI(); 						
	RPrint("after Results_... rbind_irmat_FHDI after binding :");
	rbind_irmat_FHDI.print_rbind_FHDI(); 	
	//testout of other matrices
	RPrint("after Results_... rbind_uox after binding :");
	rbind_uox.print_rbind_FHDI(); 	
	RPrint("after Results_... rbind_mox after binding :");
	rbind_mox.print_rbind_FHDI(); 	
    RPrint("ord :");   List_ord.print_List_FHDI(); 
	RPrint("ocsg :");  List_ocsg.print_List_FHDI(); 
	*/
	//------Yicheng
	//---------
	//Jackknife weights for variance estimation
	//---------
	
	//double** d_rw = New_dMatrix(nrow, nrow); 
	double Rep_CellP_begin = MPI_Wtime();

	RepWeight_FHDI d_rw(nrow);
	//RepWeight(nrow, d_rw);
	//cout << "Mynode: "<<"Rep_Debug111"<<endl;
	//RPrint("d_rw: "); RPrint(d_rw, nrow, nrow);
	List_FHDI         List_rst_prob(nrow); //only i_nc rows are meaningful
	List_string_FHDI  List_rst_name(nrow); //only i_nc rows are meaningful
	std::vector<int> List_rst_prob_size; // very important used in FHDI, the memory of List_rst_prob and List_rst_name will be tremendous if input data is big large. 
	// however, variance estimation actually does not need List_rst_prob and List_rst_name. It only need make sure size of each List_rst_prob[i] >0
	// thus, List_rst_prob_size recording size of each List_rst_prob[i] will be input in variance estimation to avoid memory issue

	std::vector<std::string> s_ncx;

	if ((i_collapsing == 0) && (i_variance == 1) && (i_imputation == 1) && (i_cellmake ==1)) {
		double** d_cx = New_dMatrix(nrow, ncol);
		Copy_dMatrix(z, nrow, ncol, d_cx);

		//double Rep_CellP_begin = MPI_Wtime();

		Rep_CellP_FEFI(d_cx, nrow, ncol, d_rw, id,
			List_rst_prob,
			List_rst_name,
			s_ncx, TestOut);
		//if (mynode == 0) {
		//	cout << " ========= Rep_CellP for big-n P-FEFI.. has successfully finished!" << endl;
		//}

		Del_dMatrix(d_cx, nrow, ncol);
	}

	if ((i_collapsing != 0) && (i_variance == 1) && (i_imputation == 1) && (i_cellmake == 1)) {
		double** d_cx = New_dMatrix(nrow, ncol);
		Copy_dMatrix(z, nrow, ncol, d_cx);

		//double Rep_CellP_begin = MPI_Wtime();

		Rep_CellP_Bigp_FEFI(d_cx, nrow, ncol, i_collapsing, d_rw, id,
			List_rst_prob,
			List_rst_name,
			s_ncx, codes, TestOut);
		//if (mynode == 0) {
		//	cout << " ========= Rep_CellP for big-p P-FEFI.. has successfully finished!" << endl;
		//}

		Del_dMatrix(d_cx, nrow, ncol);
	}


	if ((i_collapsing == 0) && (i_variance == 1) && (i_imputation == 2) && (i_cellmake == 1)) {
		double** d_cx = New_dMatrix(nrow, ncol);
		Copy_dMatrix(z, nrow, ncol, d_cx);

		//double Rep_CellP_begin = MPI_Wtime();

		Rep_CellP_FHDI(d_cx, nrow, ncol, d_rw, id,
			//List_rst_prob,
			//List_rst_name,
			List_rst_prob_size,
			s_ncx, TestOut);
		//if (mynode == 0) {
		//	cout << " ========= Rep_CellP for big-n P-FHDI has successfully finished!" << endl;
		//}

		Del_dMatrix(d_cx, nrow, ncol);
	}

	if ((i_collapsing != 0) && (i_variance == 1) && (i_imputation == 2) && (i_cellmake == 1)) {
		double** d_cx = New_dMatrix(nrow, ncol);
		Copy_dMatrix(z, nrow, ncol, d_cx);

		//double Rep_CellP_begin = MPI_Wtime();

		Rep_CellP_Bigp_FHDI(d_cx, nrow, ncol, i_collapsing, d_rw, id,
			//List_rst_prob,
			//List_rst_name,
			List_rst_prob_size,
			s_ncx, codes,TestOut);
		//if (mynode == 0) {
		//	cout << " ========= Rep_CellP for big-p P-FHDI.. has successfully finished!" << endl;
		//}

		Del_dMatrix(d_cx, nrow, ncol);
	}

	//--------------------------------------
	//---------------------------------------

	if ((i_variance == 1) && (i_imputation == 1) && (i_cellmake == 2)) {
		double** d_cx = New_dMatrix(nrow, ncol);
		Copy_dMatrix(z, nrow, ncol, d_cx);

		//double Rep_CellP_begin = MPI_Wtime();

		Rep_CellP_Neighbor_FEFI(d_cx, nrow, ncol, d_rw, id, List_nU,
			List_rst_prob,
			List_rst_name,
			s_ncx, TestOut);
		//if (mynode == 0) {
		//	cout << " ========= Rep_CellP with KNN for big-n or big-p P-FEFI.. has successfully finished!" << endl;
		//}

		Del_dMatrix(d_cx, nrow, ncol);
	}

	if ((i_variance == 1) && (i_imputation == 2) && (i_cellmake == 2)) {
		double** d_cx = New_dMatrix(nrow, ncol);
		Copy_dMatrix(z, nrow, ncol, d_cx);

		//double Rep_CellP_begin = MPI_Wtime();

		Rep_CellP_Neighbor_FHDI(d_cx, nrow, ncol, d_rw, id, List_nU,
			//List_rst_prob,
			//List_rst_name,
			List_rst_prob_size,
			s_ncx, TestOut);
		//if (mynode == 0) {
		//	cout << " ========= Rep_CellP with KNN for big-n or big-p P-FHDI.. has successfully finished!" << endl;
		//}

		Del_dMatrix(d_cx, nrow, ncol);
	}

	 //cout << "nrow is " << nrow << ", and i_collapsing is " << i_collapsing << " at node " << mynode << endl;
    // Del_iMatrix(codes, nrow, i_collapsing);
	//if ((i_collapsing != 0) && (i_variance == 1)) {

	//	Rep_CellP_Bigp(d_cx, nrow, ncol, i_collapsing, d_rw, id,
	//		List_rst_prob,
	//		List_rst_name,
	//		s_ncx, codes, TestOut);

	//	if (mynode == 0) {
	//		cout << " ========= Rep_CellP_Bigp.. has successfully finished!" << endl;
	//	}
	//}

	////Del_dMatrix(d_cx, nrow, ncol);

	//if (mynode == 0) {
	//	printf("Yicheng Running time of Rep_CellP = %f seconds\n", MPI_Wtime() - Rep_CellP_begin);
	//}
	//--------------------------Testout
	//cout <<"Mynode: "<<mynode <<", List_rst_name:" << endl;;
	//List_rst_name.print_List_string_FHDI();
	//cout << "Mynode: "<<mynode<<", List_rst_prob:" << endl;
	//List_rst_prob.print_List_FHDI();
	//if (mynode == 0) {
	//	int temp_size = 0;
	//	List_rst_prob.get_a_row_size(3, temp_size);
	//	double *temp = new double[temp_size];
	//	List_rst_prob.get_block(3, temp);
	//	cout << "Lst_test:" << endl;
	//	RPrint(temp, temp_size);
	//}


	//---------------------------
	//---------------------------
	//Variable Estimation using FEFI results
	//---------------------------
	//---------------------------
	s_M = "FEFI"; 
	int nrow_dat2_FEFI = 0; //on ALL Proc

	nrow_dat2_FEFI = rbind_ipmat_FEFI.size_row();

	double** y_bar_i_k_Summary = New_dMatrix(nrow, ncol);

	double variance_begin = MPI_Wtime();

	if ((i_variance == 1) && (i_imputation == 1)) //when FEFI's variance estimation is required
    { 
		//testout 
		//TestOut<<"begin Var Est FEFI"<<endl;
		//TestOut<<"nrow_dat2_FEFI: "<<nrow_dat2_FEFI<<endl;	
		//TestOut<<"nrow: "<<nrow<<endl;			
		//double variance_begin = MPI_Wtime();
		//TestOut << "DEBUG YICHENG" << endl;

		if(nrow_dat2_FEFI<=0) {TestOut<<"ERROR! dimension of ipmat is less than 0!"<<endl;}
	
		 //nrow_dat2_FEFI = rows of w1
		//TestOut << "DEBUG YICHENG2" << endl;

		//TestOut << "nrow_dat2_FEFI is " << nrow_dat2_FEFI << endl;
		//TestOut << "nrow is " << nrow << endl;
		Variance_Est_FEFI_Extension_cpp(x_raw, z, nrow, ncol, 
				d_rw, d_w, id, 
				rbind_ipmat_FEFI, rbind_Resp_FEFI, rbind_irmat_FEFI,
				rbind_ipmat_FHDI, rbind_Resp_FHDI, rbind_irmat_FHDI,
				rbind_uox, rbind_mox, 
				List_ord, List_ocsg, List_rst_prob, List_rst_name, s_ncx,
				s_M,
			y_bar_i_k_Summary, TestOut, TestOut_Slave1, TestOut_Slave2);
		
		//if(mynode==0){
		//	cout << "YYC Running time of Variance_FEFI at node " << mynode << " = " << MPI_Wtime() - variance_begin << endl;
		//	printf("Yicheng Running time of Variance_FEFI = %f seconds\n", MPI_Wtime() - variance_begin);
		//}
		//Del_dMatrix(wmat_FEFI, nrow_dat2_FEFI, nrow);
		//-------
		//prep return
		//-------
		//if (mynode == 0) {
		//	RPrint("YichengYang---------------------------------------------wmat_FEFI");
		//	RPrint(wmat_FEFI, nrow_dat2_FEFI, nrow);
		//}
		//TestOut << "wmat_FEFI" << endl;
		//for (int i = 0; i<nrow_dat2_FEFI; i++)
		//{
		//	for (int j = 0; j<nrow; j++) { TestOut << setw(20) << wmat_FEFI[i][j]; }
		//	TestOut << endl;
		//}
		//-----------------------------------------------------------------------------
		
		//testout
		//std::ofstream TestOut_FEFI("./debug_FEFI.txt", std::ios::out);
		//RPrint("wmat_FEFI");
		//RPrint(wmat_FEFI, nrow_dat2_FEFI, nrow, TestOut_FEFI);

		//deallocation 
		//if (mynode == 0) {
		//	RPrint("y_bar_i_k_Summary--------------------------");
		//	RPrint(y_bar_i_k_Summary, nrow, ncol);
		//}


		//Del_dMatrix(y_bar_i_k_Summary, nrow, ncol);
		//delete[] d_vrst_temp; 
		//------------------------------------------------------------------------
	}
	//---------------------------
	//---------------------------
	//Variable Estimation using FEFI results
	//---------------------------
	//---------------------------
	s_M = "FHDI"; 
	
	const int nrow_dat2_FHDI = rbind_ipmat_FHDI.size_row();

	if ( (i_variance == 1) && (i_imputation == 2)) //when FHDI's variance estimation is required
	{
		//testout 
		//TestOut<<"begin Var Est FHDI"<<endl;		
		double variance_FHDI_begin = MPI_Wtime();
		if (nrow_dat2_FHDI <= 0) { TestOut << "ERROR! dimension of ipmat is less than 0!" << endl; }

		//double** wmat_FHDI = New_dMatrix(nrow_dat2_FHDI, nrow); //nrow_dat2_FHDI = rows of w1

		Variance_Est_FHDI_Extension_cpp(x_raw, z, nrow, ncol,
			d_rw, d_w, id,
			rbind_ipmat_FEFI, rbind_Resp_FEFI, rbind_irmat_FEFI,
			rbind_ipmat_FHDI, rbind_Resp_FHDI, rbind_irmat_FHDI,
			rbind_uox, rbind_mox,
			List_ord, List_ocsg, List_rst_prob_size, s_ncx,
			s_M,
			y_bar_i_k_Summary, TestOut);

		//if (mynode == 0) {
		//	cout << "YYC Running time of Variance_FHDI at node " << mynode << " = " << MPI_Wtime() - variance_FHDI_begin << endl;
		//	printf("Yicheng Running time of Variance_FHDI = %f seconds\n", MPI_Wtime() - variance_FHDI_begin);
		//}
		//-------
		//prep return
		//-------
		//double* d_vrst_temp = new double[nrow]; 
		//for(int i=0; i<nrow_dat2_FHDI; i++)
		//{
		//	for(int j=0; j<nrow; j++) d_vrst_temp[j] = wmat_FHDI[i][j]; 
		//	rbind_vrst_FHDI.append_block(d_vrst_temp); //append row-by-row
		//}

		//testout
		//std::ofstream TestOut_FHDI("C:\\Users\\ihcho_000\\Documents\\Cho_FHDI_extended_Oct2016\\debug_FHDI.txt", std::ios::out);
		//RPrint("wmat_FHDI");
		//RPrint(wmat_FHDI, nrow_dat2_FHDI, nrow, TestOut_FHDI);

		//deallocation		
		//Del_dMatrix(wmat_FHDI, nrow_dat2_FHDI, nrow);
		//delete[] d_vrst_temp; 
		//Del_dMatrix(wmat_FHDI, nrow_dat2_FHDI, nrow);
	}


	//if (mynode == 0) {
	//	//cout << "YYC Running time of Variance_FHDI at node " << mynode << " = " << MPI_Wtime() - variance_FHDI_begin << endl;
	//	printf("Yicheng Running time of Variance_FHDI = %f seconds\n", MPI_Wtime() - variance_begin);
	//}
	//if (mynode == 0) {
	//	cout << " ========= Variance_Estimation.. has successfully finished!" << endl;
	//}
	//Del_dMatrix(y_bar_i_k_Summary,nrow,ncol);
	//-----------------------------------------------
	// Make final matrix including imputed values
	//-----------------------------------------------



	double* d_Final_Full_Data = new double[nrow*ncol];	// Final full matrix
	double* colMean = new double[ncol];					// Column mean of final matrix
	double* final_variance_data = new double[ncol];
	for (int i = 0; i < ncol; i++) {
		colMean[i] = 0.0;
	}

	if (mynode == 0) {
		//-- Extract imputed result
		if (i_imputation == 1) { //FEFI
			Extract_Imputed_Results(nrow, ncol, rbind_ipmat_FEFI, d_Final_Full_Data);
		}
		if (i_imputation == 2) { //FHDI
			Extract_Imputed_Results(nrow, ncol, rbind_ipmat_FHDI, d_Final_Full_Data);
		}

		//-- Testout final matrix
		//TestOut << "Final matrix:" << endl;
		for (int i = 0; i < nrow; i++) {
			for (int j = 0; j < ncol; j++) {
				TestOut_Slave1 << setw(20) << d_Final_Full_Data[i + j*nrow];
				colMean[j] = colMean[j] + d_Final_Full_Data[i + j*nrow] / double(nrow);	// Calculate column mean of the final matrix
				if (j % ncol == ncol-1) TestOut_Slave1 << endl;
			}
			/*if (i == nrow - 1) {
				cout << "col mean: ";
				for (int j = 0; j < ncol; j++) {
					cout << setw(20) << colMean[j];
				}
			}*/
		}
		TestOut_Slave1 << endl;

		//-- Testout column mean
		TestOut_Slave3 << "Mean estimates:" << endl;
		for (int i = 0; i < ncol; i++) {
			//TestOut << "bbfore" << endl;
			TestOut_Slave3 << setw(20) << colMean[i];
			//TestOut << "aafter" << endl;
		}
		TestOut_Slave3 << endl;
		//----------------------Extract variance---------------
		if (i_variance == 1) {
			 Extract_Variance_Results_MPI(nrow, ncol, y_bar_i_k_Summary, final_variance_data);
		}
		//RPrint("Yicheng-------------------------------Print out final variance matrix\n");
		//RPrint(final_variance_data, ncol);

		TestOut_Slave3 << "Variance Results:" << endl;
		for (int i = 0; i < ncol; i++) {
			//TestOut << "bbfore" << endl;
			TestOut_Slave3 << setw(20) << final_variance_data[i];
			//TestOut << "aafter" << endl;
		}
		TestOut_Slave3 << endl;

		TestOut_Slave3 <<"s_op_imputation: "<<endl;
		TestOut_Slave3 << i_option_imputation <<endl;

		TestOut_Slave3 <<"i_option_merge: "<<endl;
		TestOut_Slave3 << i_option_merge << endl;

		TestOut_Slave3<<"M: "<<endl;
		TestOut_Slave3<<i_M<<endl;

	} //end of if (mynode == 0)

	//=========================
	//Deallocate memories
	//=========================
	//Del_dMatrix(z, nrow, ncol); //defined at main
	//Del_dMatrix(x_raw, nrow, ncol); //defined at main()
	//Del_iMatrix(r_raw, nrow, ncol); //defined at main()
	delete[] d_Final_Full_Data;
	delete[] colMean;
	delete[] final_variance_data;
	Del_dMatrix(y_bar_i_k_Summary, nrow, ncol);
	//Del_dMatrix(d_rw, nrow, nrow); 

	//final_variance_data = NULL;
	//y_bar_i_k_Summary = NULL;
	//codes = NULL;

	//cout << "nrow is " << nrow << ", and i_collapsing is " << i_collapsing << " at node " << mynode << endl;
	//Del_iMatrix(codes, nrow, i_collapsing);


	return ; //test return of double* 

}



