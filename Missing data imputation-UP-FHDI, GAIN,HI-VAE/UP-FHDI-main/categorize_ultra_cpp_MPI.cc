
#include "Error_Check_Data.cc"

//Fn===========================================================================

//categorize_cpp.cc-----------------------------------------------------------------------------

//Fn===========================================================================

//namespace FHDI {

bool categorize_ultra_cpp(MPI_File fh_binary_daty, MPI_File fh_binary_datr, const int nrow, const int ncol, double* k,

	int* NonCollapsible_categorical, MPI_File fh_datz, ofstream& TestOut)

	//Description=========================================

	// categorize the data matrix x 

	// according to the given number of categories stored in k(ncol)

	//

	// [Algorithm I] for continuous column (or variable):  
	//
	// perc: percentiles used to get quantiles, determined by k

	// quan: quantiles if k=4, we quan=(Q1,Q2,Q3) have Q1(=1/4), Q2 (=Median) and Q3(=3/4)

	// [Algorithm II] for categorical column (or variable):

	// if a column consists of all integer values, automatically changed to categorical variable
	// and also adjust the k[] 

	// Note: as of Dec 2016, NA values (missing data) is marked by a long number at the parent "r" code

	//                       the long number is 1234567899

	// original R code: Dr. Im, J. and Dr. Kim, J. 

	// c++ code: 		Dr. Cho, I. 

	// All rights reserved

	// 

	// updated: April 9, 2018

	//----------------------------------------------------

	//INOUT	: double k(ncol)		= a vector of categories of each column of xalloc

	//OUT   : MPI_File fh_datz      = catorized matrix corresponding to original matrix x written in local storage

	//IN    : int NonCollapsible_categorical[ncol] = index vector of 0: collapsible; 1: Non-Collapsible
	//====================================================

{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	const int L = ncol; //size of d_rw 
	int numWorkPerProc = (int)floor(1.0*L / (1.0*totalnodes - 1));
	int numWorkLocalLast = L - numWorkPerProc * (totalnodes - 2);

	int L_temp = 0;
	if ((mynode != (totalnodes - 1)) && (mynode != 0)) L_temp = numWorkPerProc; // New Update
	if (mynode == (totalnodes - 1)) L_temp = numWorkLocalLast;

	//cout << "L_temp is " << L_temp << " at node " << mynode << endl;

	int startpoint = 0;
	int endpoint = 0;

	if (mynode != 0 && mynode != (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = mynode*numWorkPerProc;

		//RPrint("(Variance)Strating point and ending point on node ");RPrint(mynode);
		//RPrint("are: \n");
		//RPrint(startpoint); RPrint(endpoint);
		if (endpoint - startpoint != L_temp) {
			TestOut << "category boundary ERROR!!!" << endl;
			return 0;
		}
	}

	if (mynode == (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = (mynode - 1)*numWorkPerProc + numWorkLocalLast;

		//RPrint("(Variance)Strating point and ending point on node ");RPrint(mynode);
		//RPrint("are: \n");
		//RPrint(startpoint); RPrint(endpoint);
		if (endpoint - startpoint != L_temp) {
			TestOut << "category boundary ERROR!!!" << endl;
			return 0;
		}
	}
	//cout << "Category startpoint: " << startpoint << "; endpoint: " << endpoint << " at node " << mynode << endl;

	double* array_temp = new double[nrow];	// Final full matrix
	int* array_temp2 = new int[nrow];	// Final full matrix

	double** x_temp = New_dMatrix(nrow, L_temp);
	int** r_temp = New_iMatrix(nrow, L_temp);
	double** z_temp = New_dMatrix(nrow, L_temp);

	//ifstream ReadIn_binary_daty;
	//ifstream ReadIn_binary_datr;

	//MPI_File fh_binary_daty;
	//MPI_File fh_binary_datr;

	//ReadIn_binary_daty.open("./daty_binary.bin", ios::in | ios::binary); //read general data
	//ReadIn_binary_daty.seekg(0, ios::beg); //set the get pointer to the beginning of the file

	//ReadIn_binary_datr.open("./datr_binary.bin", ios::in | ios::binary); //read general data
	//ReadIn_binary_datr.seekg(0, ios::beg); //set the get pointer to the beginning of the file

	double read_begin = MPI_Wtime();

	//int success = 0;

	//success = MPI_File_open(MPI_COMM_WORLD, "./daty_column_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_binary_daty);
	//if (success != MPI_SUCCESS) cout << "MPI I/O fail to open the file!" << endl;

	//success = MPI_File_open(MPI_COMM_WORLD, "./datr_column_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_binary_datr);
	//if (success != MPI_SUCCESS) cout << "MPI I/O fail to open the file!" << endl;

	int counter = 0;
	for (int k = startpoint; k < endpoint; k++) {

		MPI_In_raw(nrow, k, fh_binary_daty, array_temp);
		MPI_In_raw(nrow, k, fh_binary_datr, array_temp2);

		//ReadWrite_matrix_column(nrow, ncol, ReadIn_binary_daty, k, array_temp);
		//ReadWrite_matrix_column(nrow, ncol, ReadIn_binary_datr, k, array_temp2);

		//for (int j = 0; j < nrow; j++) {
		//	x_temp[j][counter] = array_temp[j];
		//	r_temp[j][counter] = array_temp2[j];
		//}
		//counter++;

		for (int m = 0; m < nrow; m++) {
			x_temp[m][counter] = array_temp[m];
			r_temp[m][counter] = array_temp2[m];
		}
		counter++;

	}

	//success = MPI_File_close(&fh_binary_daty);
	//if (success != MPI_SUCCESS) cout << "MPI I/O fail to close the file!" << endl;
	//success = MPI_File_close(&fh_binary_datr);
	//if (success != MPI_SUCCESS) cout << "MPI I/O fail to close the file!" << endl;

	delete[] array_temp;
	delete[] array_temp2;

	//if (mynode == (totalnodes - 1)) {
	//	cout << "daty Reading test from node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < nrow; kk2++) {
	//		for (int kk3 = 0; kk3 < L_temp; kk3++) {
	//			cout << setw(20) << x_temp[kk2][kk3];
	//		}
	//		cout << endl;
	//	}


	//	cout << "datr Reading test from node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < nrow; kk2++) {
	//		for (int kk3 = 0; kk3 < L_temp; kk3++) {
	//			cout << setw(20) << r_temp[kk2][kk3];
	//		}
	//		cout << endl;
	//	}

	//}

	//----------------------------------------
	//Check input data quality
	//We did not check in main function
	//Since here is where we firstly read data
	//---------------------------------------

	bool b_Error_Data = 0;//0 is False and 1 is True

	b_Error_Data = ErrorCheck_Data(nrow, ncol, L_temp, startpoint, numWorkPerProc, numWorkLocalLast, x_temp, r_temp, TestOut);

	if (!b_Error_Data) {

		Del_dMatrix(x_temp, nrow, L_temp);

		Del_iMatrix(r_temp, nrow, L_temp);

		Del_dMatrix(z_temp, nrow, L_temp);

		return 0;
	}


	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_make at node " << mynode << " = " << MPI_Wtime() - cell_make_begin << endl;
	//	printf("Yicheng Running time of reading in category = %f seconds\n", MPI_Wtime() - read_begin);
	//}

	//if (mynode == 0) {
	//	cout << " daty matrix read test at node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < nrow; kk2++) {
	//		for (int kk3 = 0; kk3 < L_temp; kk3++) {
	//			cout << setw(20) << x_temp[kk2][kk3];
	//		}
	//		cout << endl;
	//	}

	//	cout << " datr matrix read test at node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < nrow; kk2++) {
	//		for (int kk3 = 0; kk3 < L_temp; kk3++) {
	//			cout << setw(20) << r_temp[kk2][kk3];
	//		}
	//		cout << endl;
	//	}

	//}

	for (int kk4 = 0;kk4 < nrow; kk4++) {
		for (int kk5 = 0; kk5<L_temp; kk5++) {
			if (r_temp[kk4][kk5] == 0) x_temp[kk4][kk5] = 1234567899;
		}
	}

	//---------------------------------------
	//Automatically Identify Categorical Columns (Variables)
	//---------------------------------------
	//---------------------------------------
	//global storage for Categorical variable
	//original list of category values
	//Note: assuming the largest category is 35 as of April 9, 2018
	//---------------------------------------
	int** i_category_list_original = New_iMatrix(L_temp, 35);
	Fill_iMatrix(i_category_list_original, L_temp, 35, 0);

	//Index of each column. 0: continuous; [1, 35]: categorical 
	int* i_k_categorical = new int[L_temp];
	for (int i = 0; i < L_temp; i++) i_k_categorical[i] = 0;

	//----------------------
	//Loop for each column
	//----------------------
	double* d_x_one_column = new double[nrow]; //one column of [x]
	for (int i_col = 0; i_col < L_temp; i_col++)
	{
		int i_integer_in_row = 0; //total integer counts of current column 
		bool b_categorical = 0;

		//----
		//get this column
		//----
		for (int i = 0; i < nrow; i++) d_x_one_column[i] = x_temp[i][i_col];

		//----
		//total observed cells in this column
		//----
		int i_total_observed_cells_this_column = 0;
		for (int i = 0; i<nrow; i++)
		{
			//only for the meaningful cell value of current column
			if (fabs_FHDI(d_x_one_column[i] - 1234567899) > 1e-5)
			{
				i_total_observed_cells_this_column++;
			}
		}

		//----
		//check all values in this column are integer
		//----
		for (int i_row = 0; i_row < nrow; i_row++)
		{
			double d_x_one = d_x_one_column[i_row];
			double d_round = (double)round(d_x_one);

			//when current cell value is integer & observed cell only  
			if (fabs_FHDI(d_x_one - d_round) < 1E-10 &&
				fabs_FHDI(d_x_one - 1234567899) > 1e-5)
			{
				i_integer_in_row++;
			}
		}

		//---------
		//when all values are integer
		//---------
		if (i_integer_in_row == i_total_observed_cells_this_column)
		{
			b_categorical = 1; //this column may be categorical
		}

		//testout
		//cout << "i_col:" << i_col << " i_integer_in_row:" << i_integer_in_row << " b_categorical:" << b_categorical << endl;

		//--------
		//find how many categories are
		//--------
		if (b_categorical) //when this column is categorical
		{
			std::vector<double> v_table_value;
			std::vector<int> v_table_count;
			table_cpp(d_x_one_column, nrow,
				v_table_value, v_table_count);

			int n_size_v_table = v_table_value.size(); //how many different categories

													   //---
													   //exception consideration when the missing cell is counted 
													   //as a category in the table
													   //---
			bool b_missing_cell_included = 0;
			if (n_size_v_table>1)
			{   //when the last category turns out to be the NA
				for (int j = 0; j<n_size_v_table; j++)
				{
					if (fabs_FHDI(v_table_value[j] - 1234567899) < 1e-5)
					{
						b_missing_cell_included = 1;
					}
				}
			}
			if (b_missing_cell_included) n_size_v_table = n_size_v_table - 1;

			//if the categories are less than 35 ---------
			if (n_size_v_table >= 1 && n_size_v_table < 36)
			{
				i_k_categorical[i_col] = n_size_v_table; //how many categories 
			}
			//if the categories are larger than 35 ---------
			//considered as continuous as of April 9th, 2018
			if (n_size_v_table > 35)
			{
				i_k_categorical[i_col] = 0; //0 means continuous 
				n_size_v_table = 0; //reset to zero 
			}

			//testout<<
			/*
			cout << "n_size_v_table :" << n_size_v_table << endl;
			cout << "v_table_value[] :" <<  endl;
			for (int i = 0; i < n_size_v_table; i++) cout << v_table_value[i] << " ";
			cout << endl;
			cout << "v_table_count[] :" << endl;
			for (int i = 0; i < n_size_v_table; i++) cout << v_table_count[i] << " ";
			cout << endl;
			*/

			//-------------
			//store the current column's original category values 
			//which may be not consecutive 
			//-------------
			if (i_k_categorical[i_col] >= 1 && i_k_categorical[i_col] < 36)
			{
				for (int j = 0; j < n_size_v_table; j++)
				{
					i_category_list_original[i_col][j] = static_cast<int>(v_table_value[j]);
				}
			}
			//clear vector container
			v_table_value.clear();
			v_table_count.clear();

		}

	}

	//delete local array
	delete[] d_x_one_column;

	//testout
	//RPrint("  i_k_categorical[]: ");
	//RPrint(i_k_categorical, ncol); 
	//RPrint("  i_category_list_original[][]: ");
	//RPrint(i_category_list_original, ncol, 10); 

	/*
	cout << "column,   i_k_categorical[]" << endl;
	for (int i = 0; i < ncol; i++)
	{
	cout << i + 1 << "  ,  " << i_k_categorical[i] << endl;
	}
	cout << endl;

	//testout
	cout << "column, i_category_list_original[i_col][1:35]" << endl;
	for(int i = 0; i < ncol; i++)
	{
	cout << i + 1<<"  :";
	for (int j = 0; j < 10; j++) cout << i_category_list_original[i][j] << "  ";
	cout << endl;
	}
	cout << endl;
	*/


	//-------------------------
	//Override original k[] when there are categorical columns
	//-------------------------
	int counter2 = 0;
	for (int i_col = startpoint; i_col<endpoint; i_col++)
	{
		//-----------------------
		//When this column is NON-collapsible
		//----------------------
		if (NonCollapsible_categorical[i_col] == 1)
		{
			if (i_k_categorical[counter2] >= 1 && i_k_categorical[counter2] <36)
			{
				k[i_col] = static_cast<double>(i_k_categorical[counter2]);
				if (mynode == 0) { RPrint("Note: Non-collapsible categorical variables are identified, and their {k} may be replaced with actual total categories \n"); }
			}
		}
		//-----------------------
		//When this column is COLLAPSIBLE
		//----------------------
		if (NonCollapsible_categorical[i_col] == 0)
		{
			i_k_categorical[counter2] = 0; //nullify the categorical type 
										   //and do not touch user-defined k[]
		}

		counter2++;

	}

	//if (mynode != 0) {
	//	cout << " i_k_categorical at node " << mynode << endl;
	//	for (int kk3 = 0; kk3 < L_temp; kk3++) {
	//		cout << setw(20) << i_k_categorical[kk3];
	//	}
	//	cout << endl;
	//}

	double* x_one_column = new double[nrow]; Fill_dVector(x_one_column, nrow, 0.0);

	double* x_one_column_temp = new double[nrow]; Fill_dVector(x_one_column_temp, nrow, 0.0);




	for (int i_col = 0; i_col<L_temp; i_col++)
	{

		for (int i = 0; i<nrow; i++) x_one_column[i] = x_temp[i][i_col]; //get one column



																		 //------------------------------
																		 //Algorithm II: Categorical Variable (column)
																		 //------------------------------
		if (i_k_categorical[i_col] >= 1 && i_k_categorical[i_col] <36) //if this column is categorical
		{
			for (int i = 0; i<nrow; i++)
			{
				bool b_update_z = 1;

				//only for the meaningful cell value of current column
				if (fabs_FHDI(x_one_column[i] - 1234567899) > 1e-5)
				{
					for (int i_m = 0; i_m < i_k_categorical[i_col]; i_m++) // as of April 2018, maximum categories = 35
					{
						if (fabs_FHDI(x_one_column[i] - i_category_list_original[i_col][i_m]) < 1e-5)
						{
							if (b_update_z) z_temp[i][i_col] = (i_m + 1)*1.0; //Actual Category Number. Stored as double 
							b_update_z = 0; //move to next row 

											//testout
											//RPrint(" x_one_column[i]:"); RPrint(x_one_column[i]); 
											//RPrint(" i_category_list_original[i_col][i_m]:"); 
											//RPrint(i_category_list_original[i_col][i_m]); 
						}
					}
				}
			}
		}

		//------------------------------
		//Algorithm I: Continuous Variable (column)
		//------------------------------
		if (i_k_categorical[i_col] == 0) //if this column is Continuous
		{


			//----------------

			// omit Not Available (NA) values in each column of x

			//----------------

			int i_temp = 0;

			for (int i = 0; i<nrow; i++)

			{

				if (fabs_FHDI(x_one_column[i] - 1234567899) > 1e-5)

					//if(   !std::isnan(x_one_column[i])   ) //non-NA value only	

				{

					x_one_column_temp[i_temp] = x_one_column[i];

					i_temp++;

				}

			}



			//-----------------

			//make percentile except for 1.0

			//-----------------

			int k_one_column = (int)k[startpoint + i_col];

			//cout<<"k_one_column is "<< k_one_column <<" at "<< startpoint + i_col <<"th column"<<endl;

			if (fabs_FHDI(k_one_column) <= 1.0)

			{
				RPrint("Error! in categorize_cpp, k_one_column is <=1.0! \n"); return 0;
			} //error check

			double* perc = new double[k_one_column - 1]; Fill_dVector(perc, (k_one_column - 1), 0.0);



			for (int i = 0; i<(k_one_column - 1); i++)

			{

				perc[i] = (i + 1)*(1.0 / k_one_column);

			}



			//------------------

			//quantile generation

			//the same as Type 7 (default in R)

			//------------------

			int n_observed = i_temp; //actual size of non-NA data in current column of x

			if (n_observed <= nrow)

				//{ std::sort(&x_one_column_temp[0], &x_one_column_temp[n_observed]); }

			{
				std::sort(x_one_column_temp, x_one_column_temp + n_observed);
			}

			//Note: sort happens in [begin, end)

			//Note: use <algorithm> of c++ library. formation: sort(*begin, *end)

			if (n_observed > nrow)  //error case 

			{
				RPrint("Error! n_observed > nrow in categorize()"); return 0;
			}







			//Note: the last quantile (i.e. 100%) is not included, and thus (k_one_column-1) is used

			double* x_quantile = new double[k_one_column - 1]; Fill_dVector(x_quantile, (k_one_column - 1), 0.0);



			for (int i = 0; i<(k_one_column - 1); i++)

			{

				double d_h = (n_observed - 1)*perc[i]; //+1 is removed for c++ code 

				x_quantile[i] = x_one_column_temp[int(floor(d_h))]

					+ (d_h - floor(d_h))*(x_one_column_temp[int(floor(d_h) + 1)]

						- x_one_column_temp[int(floor(d_h))]);

				//if (mynode != 0) cout<< "d_h: "<< d_h <<" and n_observed: "<< n_observed <<" at column "<< i_col <<endl;

			}

		/*	if (mynode != 0) {

				cout << " perc at column " << i_col << " at node " << mynode << endl;
				for (int kk3 = 0; kk3 < (k_one_column - 1); kk3++) {
					cout << setw(20) << perc[kk3];
				}
				cout << endl;

				cout << " x_quantile at column "<< i_col <<" at node " << mynode << endl;
				for (int kk3 = 0; kk3 < (k_one_column - 1); kk3++) {
					cout << setw(20) << x_quantile[kk3];
				}
				cout << endl;
			}*/

			//---------------

			//assign z with category values

			// Note: categories = {1, 2, ...} 

			//---------------

			for (int i = 0; i<nrow; i++)

			{

				if (fabs_FHDI(x_one_column[i] - 1234567899) > 1e-5) //non-NA value only

																	//if(   !std::isnan(x_one_column[i])   ) //non-NA value only

				{

					//---------

					//default category of non-NA unit is 1 as of 0124_2017

					//---------

					z_temp[i][i_col] = 1; //default 



										  //----------

										  //consider each quantile

										  //----------

					if (x_one_column[i] < x_quantile[0]) { z_temp[i][i_col] = 1; } //1st category

					if (x_one_column[i] > x_quantile[k_one_column - 2]) { z_temp[i][i_col] = k_one_column; } //last category



					for (int j = 1; j<(k_one_column - 1); j++)

					{

						if (x_quantile[j - 1] < x_one_column[i] && x_one_column[i] <= x_quantile[j])

						{

							z_temp[i][i_col] = j + 1; //(j+1)th quantile. Note: j =[0,k_one_column) 

							break;

						}

					}

				}
			}



			//--------------

			//local Deallocation

			//--------------

			delete[] perc;

			delete[] x_quantile;
		} //end of continuous variable (column)			

	}

	//if (mynode == 1) {
	//	cout << " z temp matrix at node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < nrow; kk2++) {
	//		for (int kk3 = 0; kk3 < L_temp; kk3++) {
	//			cout << setw(20) << z_temp[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//}

	//******************
	//******************
	//******************
	//Write columnly distributed z matrix to the same file concorrently
	//MPI_Out_datz(nrow, L_temp, numWorkPerProc, fh_datz, z_temp);
	MPI_Out_datz_rowise(nrow, L_temp, ncol, numWorkPerProc, fh_datz, z_temp);

	MPI_Barrier(MPI_COMM_WORLD);
    //====================================================================
	//=====================================================================
	// Validation purpose
	//To recover, just delete this validation part and activate writing of 
	//z matrix in line 661
	//====================================================================
	//=====================================================================

	//if (mynode != 0) {
	//	MPI_Send(z_temp[0], nrow*L_temp, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	//	MPI_Send(x_temp[0], nrow*L_temp, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
	//	MPI_Send(r_temp[0], nrow*L_temp, MPI_INT, 0, 3, MPI_COMM_WORLD);
	//}

	//double** x_sorted = New_dMatrix(nrow, ncol);
	//int** r_sorted = New_iMatrix(nrow, ncol);

	//if (mynode == 0) {

	//	double** z_matrix = New_dMatrix(nrow, ncol);
	//	double** z_temp_recv = New_dMatrix(nrow, numWorkPerProc);
	//	double** z_temp_recv_last = New_dMatrix(nrow, numWorkLocalLast);

	//	double** x_matrix = New_dMatrix(nrow, ncol);
	//	double** x_temp_recv = New_dMatrix(nrow, numWorkPerProc);
	//	double** x_temp_recv_last = New_dMatrix(nrow, numWorkLocalLast);

	//	int** r_matrix = New_iMatrix(nrow, ncol);
	//	int** r_temp_recv = New_iMatrix(nrow, numWorkPerProc);
	//	int** r_temp_recv_last = New_iMatrix(nrow, numWorkLocalLast);

	//	int startpoint_temp = 0;
	//	int end_temp = 0;

	//	for (int j = 1; j < totalnodes; j++) {

	//		if (j != (totalnodes - 1)) {
	//			MPI_Recv(z_temp_recv[0], nrow*numWorkPerProc, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
	//			MPI_Recv(x_temp_recv[0], nrow*numWorkPerProc, MPI_DOUBLE, j, 2, MPI_COMM_WORLD, &status);
	//			MPI_Recv(r_temp_recv[0], nrow*numWorkPerProc, MPI_INT, j, 3, MPI_COMM_WORLD, &status);

	//			startpoint_temp = (j - 1)*numWorkPerProc;
	//			end_temp = j*numWorkPerProc;

	//			cout << "startpoint_temp is " << startpoint_temp << "; end_temp is "<< end_temp <<" at j = "<<j<< endl;

	//			for (int k = 0; k < nrow; k++) {
	//				counter = 0;
	//				for (int l = startpoint_temp; l < end_temp; l++) {
	//					z_matrix[k][l] = z_temp_recv[k][counter];
	//					x_matrix[k][l] = x_temp_recv[k][counter];
	//					r_matrix[k][l] = r_temp_recv[k][counter];
	//					counter++;
	//				}
	//			}
	//		}

	//		if (j == (totalnodes - 1)) {
	//			MPI_Recv(z_temp_recv_last[0], nrow*numWorkLocalLast, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
	//			MPI_Recv(x_temp_recv_last[0], nrow*numWorkLocalLast, MPI_DOUBLE, j, 2, MPI_COMM_WORLD, &status);
	//			MPI_Recv(r_temp_recv_last[0], nrow*numWorkLocalLast, MPI_INT, j, 3, MPI_COMM_WORLD, &status);

	//			startpoint_temp = (j - 1)*numWorkPerProc;
	//			end_temp = (j - 1)*numWorkPerProc + numWorkLocalLast;

	//			cout << "startpoint_temp is " << startpoint_temp << "; end_temp is " << end_temp << " at j = " << j << endl;

	//			for (int k = 0; k < nrow; k++) {
	//				counter = 0;
	//				for (int l = startpoint_temp; l < end_temp; l++) {
	//					z_matrix[k][l] = z_temp_recv_last[k][counter];
	//					x_matrix[k][l] = x_temp_recv_last[k][counter];
	//					r_matrix[k][l] = r_temp_recv_last[k][counter];
	//					counter++;
	//				}
	//			}
	//		}

	//	}//end 

	//	for (int kk4 = 0;kk4 < nrow; kk4++) {
	//		for (int kk5 = 0; kk5<ncol; kk5++) {
	//			if (r_matrix[kk4][kk5] == 0) x_matrix[kk4][kk5] = 0.0;
	//		}
	//	}

	//	/*cout<<"Complete Z matrix is "<<endl;
	//	for (int p = 0; p < nrow; p++) {
	//		for (int m = 0; m < ncol;m++) {
	//			cout << setw(20) << z_matrix[p][m];
	//		}
	//		cout << endl;
	//	}

	//	cout << "Complete X matrix is " << endl;
	//	for (int p = 0; p < nrow; p++) {
	//		for (int m = 0; m < ncol;m++) {
	//			cout << setw(20) << x_matrix[p][m];
	//		}
	//		cout << endl;
	//	}

	//	cout << "Complete R matrix is " << endl;
	//	for (int p = 0; p < nrow; p++) {
	//		for (int m = 0; m < ncol;m++) {
	//			cout << setw(20) << r_matrix[p][m];
	//		}
	//		cout << endl;
	//	}*/


	//	std::string *cn_z_matrix = new std::string[nrow]; //declaration of concatenated string vector of z
	//	std::string *s_z_matrix = new std::string[nrow]; //string vector of observed patterns only
	//	
	//	Trans(z_matrix, nrow, ncol, cn_z_matrix);
	//	std::vector<int> order_z_matrix;

	//	for (int l = 0; l < nrow; l++) s_z_matrix[l] = cn_z_matrix[l];

	//	std::sort(s_z_matrix, s_z_matrix + nrow); //knowing that s_ol[] has i_ol_temp entities

	//	std::string s_temp;
	//	for (int t = 0; t < nrow; t++) {
	//		s_temp = s_z_matrix[t];
	//		for (int k = 0; k < nrow; k++) {
	//			if (s_temp.compare(cn_z_matrix[k]) == 0) {
	//				order_z_matrix.push_back(k);
	//				cn_z_matrix[k] = "99999999";
	//			}
	//		}
	//	}

	//	cout<<"order_z_matrix: "<<endl;
	//	for (int t = 0; t < nrow; t++) cout<<"order_z_matrix["<<t<<"]: "<< order_z_matrix[t] <<endl;

	//	for (int t = 0; t < nrow; t++) {
	//		for (int l = 0; l < ncol; l++) {
	//			x_sorted[t][l] = x_matrix[order_z_matrix[t]][l];
	//		}
	//	}

	//	for (int t = 0; t < nrow; t++) {
	//		for (int l = 0; l < ncol; l++) {
	//			r_sorted[t][l] = r_matrix[order_z_matrix[t]][l];
	//		}
	//	}

	//	//for (int t = 0; t < nrow; t++) {
	//	//	for (int l = 0; l < ncol; l++) {
	//	//		z_matrix_sorted[t][l] = z_matrix[order_z_matrix[t]][l];
	//	//	}
	//	//}

	//	cout << "Complete sorted X matrix is " << endl;
	//	for (int p = 0; p < nrow; p++) {
	//		for (int m = 0; m < ncol;m++) {
	//			cout << setw(20) << x_sorted[p][m];
	//		}
	//		cout << endl;
	//	}

	//	cout << "Complete sorted R matrix is " << endl;
	//	for (int p = 0; p < nrow; p++) {
	//		for (int m = 0; m < ncol;m++) {
	//			cout << setw(20) << r_sorted[p][m];
	//		}
	//		cout << endl;
	//	}

	//	Del_dMatrix(z_temp_recv, nrow, numWorkPerProc);
	//	Del_dMatrix(z_temp_recv_last, nrow, numWorkLocalLast);
	//	Del_dMatrix(z_matrix, nrow, ncol);

	//	Del_dMatrix(x_temp_recv, nrow, numWorkPerProc);
	//	Del_dMatrix(x_temp_recv_last, nrow, numWorkLocalLast);
	//	Del_dMatrix(x_matrix, nrow, ncol);

	//	Del_iMatrix(r_temp_recv, nrow, numWorkPerProc);
	//	Del_iMatrix(r_temp_recv_last, nrow, numWorkLocalLast);
	//	Del_iMatrix(r_matrix, nrow, ncol);
	//	//Del_dMatrix(z_matrix_sorted, nrow, ncol);
	//	delete[] cn_z_matrix;
	//	delete[] s_z_matrix;
	//}

	////MPI_Out_datz(nrow, ncol, fh_datz, z_matrix_sorted);
	//int success = 0;
	//MPI_File fh_sorted_daty;
	//MPI_File fh_sorted_datr;

	//success = MPI_File_open(MPI_COMM_WORLD, "./sorted_daty_column_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_sorted_daty);
	//if (success != MPI_SUCCESS) cout << "MPI I/O of z matrix fail to open the file!" << endl;

	//success = MPI_File_open(MPI_COMM_WORLD, "./sorted_datr_column_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_sorted_datr);
	//if (success != MPI_SUCCESS) cout << "MPI I/O of z matrix fail to open the file!" << endl;

	//MPI_Out_sorted_daty(nrow, ncol, fh_sorted_daty, x_sorted);
	//MPI_Out_sorted_datr(nrow, ncol, fh_sorted_datr, r_sorted);

	//success = MPI_File_close(&fh_sorted_daty);
	//if (success != MPI_SUCCESS) cout << "MPI I/O of z matrix fail to close the file!" << endl;

	//success = MPI_File_close(&fh_sorted_datr);
	//if (success != MPI_SUCCESS) cout << "MPI I/O of z matrix fail to close the file!" << endl;

	//Del_dMatrix(x_sorted, nrow, ncol);
	//Del_iMatrix(r_sorted, nrow, ncol);

	
	//--------------------------------------------------------------------
	//end of validation
	//--------------------------------------------------------------------

	//=====================================================================

	//-------------------

	//Close files

	//--------------------

	//ReadIn_binary_daty.close();
	//ReadIn_binary_datr.close();

	//--------------------

	//Deallocation

	//--------------------

	Del_iMatrix(i_category_list_original, ncol, 35);

	Del_dMatrix(x_temp, nrow, L_temp);

	Del_iMatrix(r_temp, nrow, L_temp);

	Del_dMatrix(z_temp, nrow, L_temp);

	delete[] i_k_categorical;

	delete[] x_one_column;

	delete[] x_one_column_temp;



	return 1;

}






