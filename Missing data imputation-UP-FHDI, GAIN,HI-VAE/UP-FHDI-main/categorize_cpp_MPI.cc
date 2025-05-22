//Fn===========================================================================

//categorize_cpp.cc-----------------------------------------------------------------------------

//Fn===========================================================================

//namespace FHDI {

bool categorize_cpp(double** x, const int nrow, const int ncol, double* k, double** z,
	int* NonCollapsible_categorical)

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

	//IN	: double x(nrow, ncol) 	= {y1, y2, ... } total data containing missing values

	//INOUT	: double k(ncol)		= a vector of categories of each column of xalloc

	//OUT   : double z(nrow, ncol)  = catorized matrix corresponding to original matrix x

	//                                initialized with 0.0 

	//IN    : int NonCollapsible_categorical[ncol] = index vector of 0: collapsible; 1: Non-Collapsible
	//====================================================

{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//---------------------------------------
	//Automatically Identify Categorical Columns (Variables)
	//---------------------------------------
	//---------------------------------------
	//global storage for Categorical variable
	//original list of category values
	//Note: assuming the largest category is 35 as of April 9, 2018
	//---------------------------------------
	int** i_category_list_original = New_iMatrix(ncol, 35);
	Fill_iMatrix(i_category_list_original, ncol, 35, 0);

	//Index of each column. 0: continuous; [1, 35]: categorical 
	int* i_k_categorical = new int[ncol];
	for (int i = 0; i < ncol; i++) i_k_categorical[i] = 0;

	//----------------------
	//Loop for each column
	//----------------------
	double* d_x_one_column = new double[nrow]; //one column of [x]
	for (int i_col = 0; i_col < ncol; i_col++)
	{
		int i_integer_in_row = 0; //total integer counts of current column 
		bool b_categorical = 0;

		//----
		//get this column
		//----
		for (int i = 0; i < nrow; i++) d_x_one_column[i] = x[i][i_col];

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
	for (int i_col = 0; i_col<ncol; i_col++)
	{
		//-----------------------
		//When this column is NON-collapsible
		//----------------------
		if (NonCollapsible_categorical[i_col] == 1)
		{
			if (i_k_categorical[i_col] >= 1 && i_k_categorical[i_col] <36)
			{
				k[i_col] = static_cast<double>(i_k_categorical[i_col]);
				if (mynode == 0) { RPrint("Note: Non-collapsible categorical variables are identified, and their {k} may be replaced with actual total categories \n"); }
			}
		}
		//-----------------------
		//When this column is COLLAPSIBLE
		//----------------------
		if (NonCollapsible_categorical[i_col] == 0)
		{
			i_k_categorical[i_col] = 0; //nullify the categorical type 
										//and do not touch user-defined k[]
		}

	}


	double* x_one_column = new double[nrow]; Fill_dVector(x_one_column, nrow, 0.0);

	double* x_one_column_temp = new double[nrow]; Fill_dVector(x_one_column_temp, nrow, 0.0);




	for (int i_col = 0; i_col<ncol; i_col++)
	{

		for (int i = 0; i<nrow; i++) x_one_column[i] = x[i][i_col]; //get one column



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
							if (b_update_z) z[i][i_col] = (i_m + 1)*1.0; //Actual Category Number. Stored as double 
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

			int k_one_column = (int)k[i_col];

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

			}



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

					z[i][i_col] = 1; //default 



									 //----------

									 //consider each quantile

									 //----------

					if (x_one_column[i] < x_quantile[0]) { z[i][i_col] = 1; } //1st category

					if (x_one_column[i] > x_quantile[k_one_column - 2]) { z[i][i_col] = k_one_column; } //last category



					for (int j = 1; j<(k_one_column - 1); j++)

					{

						if (x_quantile[j - 1] < x_one_column[i] && x_one_column[i] <= x_quantile[j])

						{

							z[i][i_col] = j + 1; //(j+1)th quantile. Note: j =[0,k_one_column) 

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



	//--------------------

	//Deallocation

	//--------------------
	Del_iMatrix(i_category_list_original, ncol, 35);

	delete[] i_k_categorical;

	delete[] x_one_column;

	delete[] x_one_column_temp;



	return 1;

}




//=========================================================

//=========================================================

//=========================================================

//=========================================================

//=========================================================

bool categorize_cpp(double* x, const int nrow, double &k, double* z,
	const int NonCollapsible_categorical_1)

	//Description=========================================

	// categorize a data ARRAY x 

	// according to the given number of category stored in k

	//

	// Algorithm:  

	// perc: percentiles used to get quantiles, determined by k

	// quan: quantiles if k=4, we quan=(Q1,Q2,Q3) have Q1(=1/4), Q2 (=Median) and Q3(=3/4)

	// 

	// Note: as of Oct 2016, NA values (missing data) is marked by 1234567899 at the parent "r" code

	//

	// original R code: Dr. Im, J. and Dr. Kim, J. 

	// c++ code: 		Dr. Cho, I. 

	// All rights reserved

	// 

	// updated: Oct 6, 2016

	//----------------------------------------------------

	//IN	: double x(nrow) 	= a column data containing missing values

	//INOUT	: double k   		= a number of category of the column. May be automatically changed for categorical variable

	//OUT   : double z(nrow)    = catorized array corresponding to original array x

	//                                initialized with 0.0 

	//IN    : int NonCollapsible_categorical_1 = 0: when this column is collapsible; 1: Non-collapsible 
	//====================================================

{

	//testout
	//RPrint("previoius k :"); RPrint(k); 

	//---------------------------------------
	//Automatically Identify Categorical Columns (Variables)
	//---------------------------------------
	//---------------------------------------
	//global storage for Categorical variable
	//original list of category values
	//Note: assuming the largest category is 35 as of April 9, 2018
	//---------------------------------------
	//const int ncol = 1; //for this function only 
	//int** i_category_list_original = New_iMatrix(ncol, 35);
	//Fill_iMatrix(i_category_list_original, ncol, 35, 0);
	int* i_category_list_original = new int[35];
	Fill_iVector(i_category_list_original, 35, 0);

	//Index of each column. 0: continuous; [1, 35]: categorical 
	//int* i_k_categorical = new int[ncol]; 
	//for (int i = 0; i < ncol; i++) i_k_categorical[i] = 0;
	int i_k_categorical = 0;

	//----------------------
	//Loop for each column
	//----------------------
	double* d_x_one_column = new double[nrow]; //one column of [x]
											   //for (int i_col = 0; i_col < ncol; i_col++)
											   //{
	int i_integer_in_row = 0; //total integer counts of current column 
	bool b_categorical = 0;

	//----
	//get this column
	//----
	//for (int i = 0; i < nrow; i++) d_x_one_column[i] = x[i][i_col];
	for (int i = 0; i < nrow; i++) d_x_one_column[i] = x[i]; //for this function only 

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
		if (n_size_v_table > 0 && n_size_v_table < 36)
		{
			//i_k_categorical[i_col] = n_size_v_table; //how many categories 
			i_k_categorical = n_size_v_table; //how many categories 
		}
		//if the categories are larger than 35 ---------
		//considered as continuous as of April 9th, 2018
		if (n_size_v_table > 35)
		{
			//i_k_categorical[i_col] = 0; //0 means continuous 
			i_k_categorical = 0; //0 means continuous 
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
		//if (i_k_categorical[i_col] >= 1 && i_k_categorical[i_col] < 36)
		if (i_k_categorical >= 1 && i_k_categorical < 36)
		{
			for (int j = 0; j < n_size_v_table; j++)
			{
				//i_category_list_original[i_col][j] = static_cast<int>(v_table_value[j]); 
				i_category_list_original[j] = static_cast<int>(v_table_value[j]);

			}
		}
		//clear vector container
		v_table_value.clear();
		v_table_count.clear();

	}

	//} //loop for columns inactivated for this function only 

	//delete local array
	delete[] d_x_one_column;

	//testout
	//RPrint("i_k_categorical :"); RPrint(i_k_categorical); 
	//RPrint("i_category_list_original []:"); RPrint(i_category_list_original, 35); 

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
	//const int i_col = 0 ;
	//if(i_k_categorical[i_col] >= 1 && i_k_categorical[i_col] <36)

	//-------
	//when this column is NON-collapsible
	//-------
	if (NonCollapsible_categorical_1 == 1)
	{
		if (i_k_categorical >= 1 && i_k_categorical <36)
		{
			//k[i_col] = static_cast<double>(i_k_categorical[i_col]); 
			k = static_cast<double>(i_k_categorical);
			//RPrint("Note! Some categorical columns are automatically identified and {k} may be replaced! \n");  
		}
	}
	//-------
	//when this column is Collapsible
	//-------	
	if (NonCollapsible_categorical_1 == 0)
	{
		i_k_categorical = 0; //nullify the categorical type
							 //and do not touch user-defined k
	}



	//below is for matrix version 
	/*
	for(int i_col=0; i_col<ncol; i_col++)
	{
	if(i_k_categorical[i_col] >= 1 && i_k_categorical[i_col] <36)
	{
	k[i_col] = static_cast<double>(i_k_categorical[i_col]);
	RPrint("Caution! some categorical columns are automatically identified and {k} may be replaced!");
	}
	}
	*/

	//testout
	//RPrint("maybe new k :"); RPrint(k); 


	double* x_one_column = new double[nrow]; Fill_dVector(x_one_column, nrow, 0.0);

	double* x_one_column_temp = new double[nrow]; Fill_dVector(x_one_column_temp, nrow, 0.0);



	for (int i = 0; i<nrow; i++) x_one_column[i] = x[i]; //get the one column

														 //testout
														 //RPrint("previous x[]:"); RPrint(x_one_column, nrow); 


														 //------------------------------
														 //Algorithm II: Categorical Variable (column)
														 //------------------------------
														 //if(i_k_categorical[i_col] >= 1 && i_k_categorical[i_col] <36) //if this column is categorical
	if (i_k_categorical >= 1 && i_k_categorical <36) //if this column is categorical
	{
		for (int i = 0; i<nrow; i++)
		{
			bool b_update_z = 1;
			//only for the meaningful cell value of current column
			if (fabs_FHDI(x_one_column[i] - 1234567899) > 1e-5)
			{
				//for(int i_m=0; i_m < i_k_categorical[i_col]; i_m++) // as of April 2018, maximum categories = 35
				for (int i_m = 0; i_m < i_k_categorical; i_m++) // as of April 2018, maximum categories = 35
				{
					//if(fabs_FHDI(x_one_column[i] - i_category_list_original[i_col][i_m]) <1e-5)					
					if (fabs_FHDI(x_one_column[i] - i_category_list_original[i_m]) <1e-5)
					{
						//if(b_update_z) z[i][i_col] = (i_m + 1); //Actual Category Number 
						if (b_update_z) z[i] = static_cast<double>(i_m + 1); //Actual Category Number //for this function only
						b_update_z = 0; //move to next row 
					}
				}
			}
		}
	}

	//------------------------------
	//Algorithm I: Continuous Variable (column)
	//------------------------------
	//if(i_k_categorical[i_col] == 0 ) //if this column is Continuous	
	if (i_k_categorical == 0) //if this column is Continuous	
	{

		//----------------

		// omit Not Available (NA) values in each column of x

		//----------------

		int i_temp = 0;

		for (int i = 0; i<nrow; i++)

		{

			if (fabs_FHDI(x_one_column[i] - 1234567899) > 1e-5)

				//if(  !std::isnan(x_one_column[i])  ) 	

			{

				x_one_column_temp[i_temp] = x_one_column[i];

				i_temp++;

			}

		}



		//-----------------

		//make percentile except for 1.0

		//-----------------

		int k_one_column = (int)k;

		if (fabs_FHDI(k_one_column) <= 1.0) { RPrint("Error! in categorize_cpp, k_one_column is <=1.0!"); return 0; } //error check



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

		}



		//---------------

		//assign z with category values

		// Note: categories = {1, 2, ...} 

		//---------------

		for (int i = 0; i<nrow; i++)

		{
			//----------
			//Avoid error by updating NA z with non-zero value during cell collapse 
			//----------
			z[i] = 0; //for general default  


			if (fabs_FHDI(x_one_column[i] - 1234567899) > 1e-5) //non-NA value only

																//if(   !std::isnan(x_one_column[i])   ) //non-NA value only

			{


				//default category of non-NA unit is 1 as of 0124_2017

				//---------

				z[i] = 1; //default 



				if (x_one_column[i] < x_quantile[0]) { z[i] = 1; } //1st category

				if (x_one_column[i] > x_quantile[k_one_column - 2]) { z[i] = k_one_column; } //last category



				for (int j = 1; j<(k_one_column - 1); j++)

				{

					if (x_quantile[j - 1] < x_one_column[i] && x_one_column[i] <= x_quantile[j])

					{

						z[i] = j + 1; //(j+1)th quantile. Note: j =[0,k_one_column) 

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
	} //end of continuous variable 


	  //--------------------

	  //Deallocation

	  //--------------------
	  //Del_iMatrix(i_category_list_original, ncol, 35);
	  //delete[] i_k_categorical;  
	delete[] i_category_list_original;

	delete[] x_one_column;

	delete[] x_one_column_temp;



	return 1;

}




