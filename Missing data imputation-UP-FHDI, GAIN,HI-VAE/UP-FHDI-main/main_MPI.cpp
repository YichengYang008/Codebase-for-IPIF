//-----------------------------------------------------------------
//Ultra Data-Oriented Parallel Fractional Hot Deck Imputation (UP-FHDI)
//
//A general-purpose, assumption-free imputation software for curing any incomplete data sets
//
//Version: 2.0
//Last release date: September 16, 2021
//Developers: Yicheng Yang, Dr.Jae-Kwang Kim, Dr.In-Ho Cho (Iowa State University)
//Contact: icho@iastate.edu
//License: GPL >= 2
//Depend: Intel MPI module ver.15.0.2 through 19.4
//Repository: https://sites.google.com/site/ichoddcse2017/home/type-of-trainings
//References:
//Serial version R package FHDI, https://cran.r-project.org/packages=FHDI.
//Jongho Im, In Ho Cho, and Jaekwang Kim, 2018, The R Journal, Vol.10(1), 140-154. [https://journal.r-project.org/archive/2018/RJ-2018-020/index.html].
//-----------------------------------------------------------------

using namespace std;    //this is very important for cluster sys.

#include <mpi.h>
#include <iostream>
#include <limits> // For NaN 
#include <iostream> // For cout, cerr, endl, and flush 
#include <assert.h> // For assert
#include <algorithm> // For sort 
#include <string> //For string array
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <strstream>
#include <iomanip>

//Fn===========================================================================
//matrix_utility_FHDI.cc-------------------------------------------------------
//Fn===========================================================================
//below local functions for avoiding error for other compilers 
int    fabs_FHDI(int x)    { if(x>=0)   {return x;} if(x<0)   {return x*-1;}   return x;}
double fabs_FHDI(double x) { if(x>=0.0) {return x;} if(x<0.0) {return x*-1.0;} return x;}

bool isnan_FHDI(double x) { return x!=x; } //added to avoid error regarding std::isnan

#include "MPI_IO.cpp"
#include "ReadWrite_matrix.cpp"
#include "matrix_utility_FHDI.cc" //for local matrix utilities

#include "StringStream_utility.cpp"
#include "ReadInput_FHDI_MPI.cc"
#include "CheckInputdata.cc"
#include "base_FHDI_MPI.cc"
#include "List_FHDI_MPI.cc"
#include "List_string_FHDI_MPI.cc"
#include "rbind_FHDI_MPI.cc"
#include "Rfn_test_MPI.cc"
#include "Rfn_test_ultra_MPI.cc"

#include "Error_Check.cc"


int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	//==============
	//MPI initial setting
	//==============
	int totalnodes, mynode;
	double d_begin_MPI, d_end_MPI, d_begin_MPI_temp; //checking time

	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	
	MPI_Barrier(MPI_COMM_WORLD);
	d_begin_MPI = MPI_Wtime();
	if(mynode == 0) d_begin_MPI_temp = MPI_Wtime();

	//==============================
	//data input
	//==============================

	ifstream* pFileInput = new ifstream
		("./input.txt");

	//============
	//for Window
	//============   
	//ofstream TestOut("./debug.txt", ios::out);
	//ofstream TestOut1("./debug_FEFI.txt", ios::out);
	ofstream TestOut;
	ofstream TestOut_Slave1;
	ofstream TestOut_Slave2;
	ofstream TestOut_Slave3;


	if (mynode == 0) 
	{
		TestOut.open("./debug.txt", ios::out);

		TestOut_Slave1.open("./simp.data.txt", ios::out);

		TestOut_Slave2.open("./fimp.data.txt", ios::out);

		TestOut_Slave3.open("./summary.data.txt", ios::out);
	}

	if (pFileInput->bad()) //cf. opposite member fn= good()
	{
		TestOut << "Can't open input file" << endl;
		return 0;
	}



    //Input parameters required for all cases (16)
	int    i_option_read_data = 0;
	int    i_option_ultra = 0;
	int    i_option_perform = 0;
	int    i_option_imputation = 0;
	int    i_option_variance = 0;
	int    i_option_merge = 0;
	int    nrow_x = 0;
	int	   ncol_x = 0; //pointers to integer sizes of x
	int    M = 0; //number of donors for FHDI
	int    i_user_defined_datz = 0; 
	int    i_option_collapsing = 0;//number of selected most correlated values
	int    i_option_SIS_type = 0;//methods of sure independence screening
	int    top;//top correlation ranking
	int    i_option_cellmake = 0;//methods of cell construction
	int    i_option_var_type = 0;//methods of cell construction
	int    memory = 0; //memory capacity of a processor


	ReadInput_FHDI_MPI(pFileInput, i_option_read_data, i_option_ultra,
		i_option_perform, nrow_x,
		ncol_x, i_option_imputation,
		i_option_variance,
		i_option_merge,
		M,
		i_user_defined_datz, i_option_collapsing, i_option_SIS_type, top, i_option_cellmake, i_option_var_type,
		memory,
		TestOut);

	//cout << "Running time after ReadInput_FHDI_MPI at node " << mynode << " = " << MPI_Wtime() - d_begin_MPI << endl;

	double **x = NULL;
	int **r = NULL;
	double **z = NULL;

	double **x_temp = NULL;
	int **r_temp = NULL;
	double **z_temp = NULL;

	double *k_temp = NULL;
	double *d_temp = NULL;
	int    *id_temp = NULL;
	int    *NonCollapsible_categorical_temp = NULL;

	if (i_option_ultra == 0) {

		x = New_dMatrix(nrow_x, ncol_x);          	//pointer to double vector x[col*row] that contains all data with missing units
		r = New_iMatrix(nrow_x, ncol_x); 			//pointer to an integer vector r[n_total_x] that contains indices of 0 and 1
		z = New_dMatrix(nrow_x, ncol_x);			//not used as of April, 2017

		x_temp = New_dMatrix(nrow_x, ncol_x);          	//pointer to double vector x[col*row] that contains all data with missing units
		r_temp = New_iMatrix(nrow_x, ncol_x); 			//pointer to an integer vector r[n_total_x] that contains indices of 0 and 1
		z_temp = New_dMatrix(nrow_x, ncol_x);			//not used as of April, 2017

		k_temp = new double[ncol_x];			//category for each column 
		d_temp = new double[nrow_x];  		//weight for each row 
		id_temp = new int[nrow_x]; 			//id of each row 
		NonCollapsible_categorical_temp = new int[ncol_x];
	}

	double *k = new double[ncol_x];			//category for each column 
	double *d = new double[nrow_x];  		//weight for each row 
	int    *id = new int[nrow_x]; 			//id of each row 
	int    *NonCollapsible_categorical = new int[ncol_x]; 


	if ( (mynode == 0) && (i_option_read_data == 0) && (i_option_ultra ==0)) {
		//READ daty, datr, category, NonCollapsible_categorical, weight, datz in a single file (input.txt)
		//READ id
		InputData(pFileInput,
			i_option_perform, nrow_x,
			ncol_x, i_option_imputation,
			i_option_variance,
			i_option_merge,
			i_user_defined_datz,

			x_temp, r_temp, z_temp, k_temp, d_temp, NonCollapsible_categorical_temp,
			TestOut);

		for (int i = 0; i < nrow_x; i++) {
			id_temp[i] = i + 1;
		}

	}

	if ((mynode == 0) && (i_option_read_data == 1) && (i_option_ultra == 0)) {

		ifstream* pFileInput_daty = new ifstream
		("./daty.txt");
		ifstream* pFileInput_datr = new ifstream
		("./datr.txt");
		//ifstream* pFileInput_datz = new ifstream
		//("./datz.txt");

		if (pFileInput_daty->bad()) //cf. opposite member fn= good()
		{
			TestOut << "Can't open daty.txt file" << endl;
			return 0;
		}

		if (pFileInput_datr->bad()) //cf. opposite member fn= good()
		{
			TestOut << "Can't open datr.txt file" << endl;
			return 0;
		}

		//READ daty from a single file (daty.txt)
		//READ datr from a single file (datr.txt)
		//READ category, NonCollapsible_categorical, weight, datz from a single file (input.txt)
		//READ id
		InputData_Seperate(pFileInput, pFileInput_daty, pFileInput_datr,
			i_option_perform, nrow_x,
			ncol_x, i_option_imputation,
			i_option_variance,
			i_option_merge,
			i_user_defined_datz,

			x_temp, r_temp, z_temp, k_temp, d_temp, NonCollapsible_categorical_temp,
			TestOut);

		for (int i = 0; i < nrow_x; i++) {
			id_temp[i] = i + 1;
		}
	}

	if ((i_option_read_data == 1) && (i_option_ultra == 1)) {

		//Note there will be no file lock if input.txt file is not
		//generated within the execution of UP-FHDI. Thus,
		//all processors can read input.txt concurrently

		//READ category, NonCollapsible_categorical, weight
		//READ id

		InputData_Seperate_ultra(pFileInput, 
			nrow_x,
			ncol_x, 
			k, 
			d, 
			NonCollapsible_categorical,
			TestOut);

		for (int i = 0; i < nrow_x; i++) {
			id[i] = i + 1;
		}
	}

	if ((i_option_read_data == 0) && (i_option_ultra == 1)) {

		ifstream* pFileInput_daty = new ifstream
		("./daty.txt");
		ifstream* pFileInput_datr = new ifstream
		("./datr.txt");
		//ifstream* pFileInput_datz = new ifstream
		//("./datz.txt");

		if (pFileInput_daty->bad()) //cf. opposite member fn= good()
		{
			TestOut << "Can't open daty.txt file" << endl;
			return 0;
		}

		if (pFileInput_datr->bad()) //cf. opposite member fn= good()
		{
			TestOut << "Can't open datr.txt file" << endl;
			return 0;
		}

		//READ daty from a single file (daty.txt)
		//READ datr from a single file (datr.txt)
		//READ category, NonCollapsible_categorical, weight, datz from a single file (input.txt)
		//READ id
		InputData_Seperate_ultra_txt(pFileInput, pFileInput_daty, pFileInput_datr,
			nrow_x,
			ncol_x,
			k,
			d,
			NonCollapsible_categorical,
			TestOut);

		for (int i = 0; i < nrow_x; i++) {
			id[i] = i + 1;
		}
	}

	//cout<<"NonCollapsible_categorical at node "<<mynode<<endl;
	//for (int t = 0; t < ncol_x; t++) {
	//	cout << setw(20) << NonCollapsible_categorical[t];
	//}
	//cout << endl;

	//cout << "k at node " << mynode << endl;
	//for (int t = 0; t < ncol_x; t++) {
	//	cout << setw(20) << k[t];
	//}
	//cout << endl;

	//cout<<"Finish input file reading at node "<<mynode<<endl;

	if (i_option_ultra == 0) {
		MPI_Bcast(x_temp[0], nrow_x*ncol_x, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(r_temp[0], nrow_x*ncol_x, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(z_temp[0], nrow_x*ncol_x, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		Copy_dMatrix(x_temp, nrow_x, ncol_x, x);
		Copy_iMatrix(r_temp, nrow_x, ncol_x, r);
		Copy_dMatrix(z_temp, nrow_x, ncol_x, z);

		Del_dMatrix(x_temp, nrow_x, ncol_x);          	//pointer to double vector x[col*row] that contains all data with missing units
		Del_iMatrix(r_temp, nrow_x, ncol_x); 			//pointer to an integer vector r[n_total_x] that contains indices of 0 and 1
		Del_dMatrix(z_temp, nrow_x, ncol_x);			//not used as of April, 2017

		MPI_Bcast(k_temp, ncol_x, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(d_temp, nrow_x, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(id_temp, nrow_x, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(NonCollapsible_categorical_temp, ncol_x, MPI_INT, 0, MPI_COMM_WORLD);

		Copy_dVector(k_temp, ncol_x, k);
		Copy_dVector(d_temp, nrow_x, d);
		Copy_iVector(id_temp, nrow_x, id);
		Copy_iVector(NonCollapsible_categorical_temp, ncol_x, NonCollapsible_categorical);

		delete[] k_temp;
		delete[] d_temp;
		delete[] id_temp;
		delete[] NonCollapsible_categorical_temp;

	}

	//cout << "Finish bcast at node " << mynode << endl;
	//cout << "Running time after InputData at node " << mynode << " = " << MPI_Wtime() - d_begin_MPI << endl;

	//testout
	//TestOut << " Input reading has finished successfully ==================" << endl;
	
	//===============================
	//================================
	//Basic Error Check
	//as of July 5, 2021
	//===============================
	//================================

	bool b_ERROR = 0;//0 is False and 1 is True

	b_ERROR = ErrorCheck(i_option_read_data, i_option_ultra, i_option_perform, i_option_imputation, i_option_variance,
		i_option_merge, nrow_x, ncol_x, M, i_user_defined_datz, i_option_collapsing,
		i_option_SIS_type, top, i_option_cellmake, i_option_var_type, memory,
		k, NonCollapsible_categorical, d, TestOut);

	if (!b_ERROR) {
		delete[] k;
		delete[] d;
		delete[] id;
		delete[] NonCollapsible_categorical;

		if (i_option_ultra == 0) {
			Del_dMatrix(x, nrow_x, ncol_x);
			Del_iMatrix(r, nrow_x, ncol_x);
			Del_dMatrix(z, nrow_x, ncol_x);
		}

		MPI_Finalize();
		TestOut.close();
		TestOut_Slave1.close();
		TestOut_Slave2.close();
		TestOut_Slave3.close();

		remove("./simp.data.txt");
		remove("./fimp.data.txt");
		remove("./summary.data.txt");

		return(0);
	}

	if (i_option_ultra == 0) {
		CheckInputdata(x, r, nrow_x, ncol_x, TestOut);

		//-------------
		//NA -> long number in x, i.e., daty
		//-------------
		for (int i = 0; i < nrow_x; i++)
			for (int j = 0; j < ncol_x; j++)
			{
				if (r[i][j] == 0) x[i][j] = 1234567899;
			}
	}

	//================================
	//End of Basic Error Check
	//===============================


	//-------------
	//Non-collapsible case consideration
	//-------------
	//int* NonCollapsible_categorical = new int[ncol_x]; 
	//Fill_iVector(NonCollapsible_categorical, ncol_x, 0); //default setting all zero. To be activated later
	//cout<<"Begin of Rfn_test-----------------------------"<<endl;
	if (i_option_ultra == 0) {

		//--------
		//prep return variables
		//--------
		rbind_FHDI  rbind_ipmat_FEFI_return(4 + ncol_x); //column size is 4+ncol
		rbind_FHDI  rbind_Resp_FEFI_return(ncol_x + 1);  //separate response matrix 
		rbind_FHDI  rbind_irmat_FEFI_return(5 + ncol_x); //column size is 5+ncol    
		rbind_FHDI  rbind_ipmat_FHDI_return(4 + ncol_x); //column size is 4+ncol
		rbind_FHDI  rbind_Resp_FHDI_return(ncol_x + 1);  //separate response matrix 
		rbind_FHDI  rbind_irmat_FHDI_return(5 + ncol_x); //column size is 5+ncol    
		rbind_FHDI  rbind_vrst_FEFI_return(nrow_x);    //variance estimates of FEFI
		rbind_FHDI  rbind_vrst_FHDI_return(nrow_x);    //variance estimates of FHDI

													   //below is for output for Cell Make only option
		rbind_FHDI  rbind_uox_return(ncol_x); //unique observed patterns
		rbind_FHDI  rbind_mox_return(ncol_x); //unique observed patterns
		rbind_FHDI  rbind_category_return(ncol_x); //cagetorized matrix 

												   //below is for output for Cell Prob only option
		std::vector<std::string> jp_name_return_CellProb;   //name of the joint probability table
		std::vector<double> jp_prob_return_CellProb; //the latest joint probability 	


		Rfn_test_MPI(x, r, nrow_x, ncol_x, k, d, M,
			i_option_imputation, i_option_variance, i_option_collapsing, i_option_SIS_type, top, i_option_cellmake,
			id, z,
			NonCollapsible_categorical,
			rbind_ipmat_FEFI_return, rbind_Resp_FEFI_return, rbind_irmat_FEFI_return,
			rbind_ipmat_FHDI_return, rbind_Resp_FHDI_return, rbind_irmat_FHDI_return,
			rbind_vrst_FEFI_return, rbind_vrst_FHDI_return,

			rbind_uox_return, rbind_mox_return, rbind_category_return,

			jp_name_return_CellProb, jp_prob_return_CellProb,

			i_option_perform,
			i_option_merge,
			TestOut, TestOut_Slave1, TestOut_Slave2, TestOut_Slave3, d_begin_MPI);
	}

	if (i_option_ultra == 1) {

		Rfn_test_ultra_MPI(nrow_x, ncol_x, k, d, M,
			i_option_imputation, i_option_variance, i_option_collapsing, i_option_SIS_type, top, i_option_cellmake, i_option_var_type,
			id,
			NonCollapsible_categorical,

			i_option_perform,
			i_option_merge, memory,
			TestOut, TestOut_Slave1, TestOut_Slave2, TestOut_Slave3, d_begin_MPI);
	}

	//cout<<"Rfn_test has successfully finished at node "<<mynode<<endl;

	//-----------------------
	//Deallocation
	//-----------------------
	if (i_option_ultra == 0) {
		Del_dMatrix(x, nrow_x, ncol_x);          	//pointer to double vector x[col*row] that contains all data with missing units
		Del_iMatrix(r, nrow_x, ncol_x); 			//pointer to an integer vector r[n_total_x] that contains indices of 0 and 1
		Del_dMatrix(z, nrow_x, ncol_x);			//not used as of April, 2017
	}

	delete[] k;
	delete[] d;
	delete[] id;
	delete[] NonCollapsible_categorical;

	//system("PAUSE");
	//MPI_Barrier(MPI_COMM_WORLD);
	d_end_MPI = MPI_Wtime();
	//if (mynode == 0) cout << "YYC Total running time = " << d_end_MPI - d_begin_MPI;
	//if (mynode == 0) cout<<"Total Time is "<< d_end_MPI - d_begin_MPI <<endl;
	//if (mynode == 0) printf(" YYC Total Runtime =%f senconds\n", d_end_MPI - d_begin_MPI);

	MPI_Finalize();
	TestOut.close();
	TestOut_Slave1.close();
	TestOut_Slave2.close();
	TestOut_Slave3.close();

	//=============================
	//Delete unnecessary txt files
	//=============================
	if (mynode == 0 && i_option_ultra == 1) {
		remove("./simp.data.txt");
		remove("./fimp.data.txt");
	}

	return 0; 
}																																												 

