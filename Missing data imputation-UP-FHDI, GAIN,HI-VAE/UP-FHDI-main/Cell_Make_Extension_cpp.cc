//Fn===========================================================================

//Cell_Make_Extension_cpp.cc-----------------------------------------------------------------------------

//Fn===========================================================================

//namespace FHDI{

#include "Merge_Extension_cpp.cc"


bool Cell_Make_Extension_cpp(double** x, const int nrow, const int ncol, double* d_k, 

							 int* NonCollapsible_categorical,
							 
							 double** z, 

						 	rbind_FHDI &rbind_uox_CellMake, 

							rbind_FHDI &rbind_mox_CellMake, 

							const int i_merge, 
							
							ofstream& TestOut, double d_begin_MPI) 

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

// c++ code: 		Dr. Cho, I. 

// All rights reserved

// 

// updated: March 28, 2017

//----------------------------------------------------

//IN	: double x(nrow, ncol) 	= {y1, y2, ... } total data containing missing values

//IN	: double d_k(ncol)		= a vector of categories of each column of xalloc

//IN    : int NonCollapsible_categorical(ncol) = {0,0, .., 1,.. 0} 
//				index for non-collapsible categorical variables. 
//				when at least one column has "1" skip cell-collapse procedure
//				this may casue a potential error of lack of enough donor! 
//				(2018, 04 21) 
//											  

//OUT   : double z(nrow, ncol)  = catorized matrix corresponding to original matrix x

//                                initialized with 0.0 

//OUT	: rbind_FHDI rbind_uox_CellMake(ncol); //compact storage of uox, unique observed rows sorted in the ascending order

//OUT	: rbind_FHDI rbind_mox_CellMake(ncol); //compact storage of mox, unique observed rows sorted in the ascending order

//

//IN    : int i_merge = random donor selection in Merge algorithm in Cell Make

//                            0= no random seed number setting

//						      1= random seed number setting 

//====================================================

{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;
	double cell_make_check1 = MPI_Wtime();
	const int n_max_iteration = nrow*2; //maximum number of iterations 

	//const int n_max_iteration = 2;  //temporary

	//-------------------------------------
	//Determine if there is Non-Collapsible Categorical variable
	//-------------------------------------
	double* d_k_Collapsible = new double[ncol]; //k for collapsible variables only 
	Copy_dVector(d_k, ncol, d_k_Collapsible); 
	
	int i_NonCollapsible_categ_total = 0; 
	for(int i=0; i<ncol; i++)
	{
		i_NonCollapsible_categ_total += NonCollapsible_categorical[i];
		
		if(NonCollapsible_categorical[i] == 0) d_k_Collapsible[i] = d_k[i]; //use user-defined k 
		if(NonCollapsible_categorical[i] == 1) d_k_Collapsible[i] = 1; //will be overwritten by actual total categories
	}

	//testout
	/*RPrint("initial d_k[]"); RPrint(d_k, ncol); 
	RPrint("\n"); 
	RPrint("d_k_Collapsible[]"); RPrint(d_k_Collapsible, ncol); 
	RPrint("\n"); 
	RPrint("NonCollapsible_categorical[]"); RPrint(NonCollapsible_categorical, ncol); 
	RPrint("\n"); 
	*/
	

	//-------------------------------------

	//Categorize raw data

	//-------------------------------------
	//Note: when there is non-collapsible variable, 
	//      its associated d_k
	//      is replaced with actual total categories 
	bool b_success_categorize = categorize_cpp(x, nrow, ncol, d_k, z, 
											   NonCollapsible_categorical);


	if(!b_success_categorize) 
	{
		//early deallocation 
		delete[] d_k_Collapsible; 
		
		return 0;
	}

	//testout
//   	RPrint("After initial categorize()  \n"); 

	//RPrint("d_k: \n"); RPrint(d_k,  ncol);

	//RPrint("z: \n"); RPrint(z, nrow, ncol);


	//----------
	//clear category matrix for possible garbage
	//Note: z has only positive integer as category #
	//----------
	/*for(int i=0; i<nrow; i++)
	{
		for(int j=0; j<ncol; j++)
		{
			if(fabs_FHDI(z[i][j] < 1e-3)) z[i][j] = 0.0; 
		}
	}
	*/

	
	//-------------------------------------

	//make a copy of z

	//-------------------------------------

	double** zbase = New_dMatrix(nrow, ncol);

	Copy_dMatrix(z, nrow, ncol, zbase);

		

	//-------------------------------------

	//sort in the order of high missing rate

	//at the end, i_orn has the "actual" column numbers from highest missing rate

	//            to the lowest missing rate

	//-------------------------------------

	int* i_orn = new int[ncol]; 	 	Fill_iVector(i_orn, ncol, 0);

	int* i_orn_temp = new int[ncol]; 	Fill_iVector(i_orn_temp, ncol, 0);

	int* i_orn_temp2 = new int[ncol]; 	Fill_iVector(i_orn_temp2, ncol, 0);



	int i_temp = 0; 

	for(int i_col=0; i_col<ncol; i_col++)

	{

		i_temp = 0;

		for(int i_row=0; i_row<nrow; i_row++)

		{

			if(fabs_FHDI(z[i_row][i_col]) < 1e-5) //count only  "0"

			{ i_temp++; }

		}

		i_orn_temp[i_col] = i_temp; //store how many "0" in this column

	}

	Copy_iVector(i_orn_temp, ncol, i_orn_temp2); //store before sorting 

	//std::sort(&i_orn_temp[0], &i_orn_temp[ncol-1]); //this works well, but not recommended

	std::sort(i_orn_temp, i_orn_temp+ncol);  	

	

	for(int i=0; i<ncol; i++)

	{

		i_temp = i_orn_temp[ncol-1-i]; //reversed searching since the "sort" occurred in ascending order

		for(int j=0; j<ncol; j++)

		{

			if(i_temp == i_orn_temp2[j])

			{

				i_orn[i] = j+1; //store column number (actual number, 1, 2, ...) 

				i_orn_temp2[j] = -1; //not to be found again 

				break;

			}

		}

	}

	

	//-------------------------------------

	//create concatenated vector of z

	//-------------------------------------

	//Note: after Zmat..() all of the below variables are updated 

	//-------------------------------------

	//std::string cn[nrow]; //declaration of concatenated string vector of z

	std::string *cn = new std::string[nrow]; //declaration of concatenated string vector of z

	int* ml = new int[nrow];

	int* ol = new int[nrow];

	double** uox = New_dMatrix(nrow, ncol);

	double** mox = New_dMatrix(nrow, ncol);

	int i_count_ol;

	int i_count_ml; 

	int i_count_uox;

	int i_count_mox; 



	std::vector<int> v_nD; 

	List_FHDI List_nU(nrow); //default for the size of nrow, but will be updated in the main loop

	int* tnU = new int[nrow]; Fill_iVector(tnU, nrow, 0); //this default size will be udpated in the main loop

	int i_cellmake = 1; // actiavte the b_success_nDAU because cell make may not have enough donors

	//cout << " Cell_Make_check1 at node " << mynode << " = " << MPI_Wtime() - cell_make_check1 << endl;

	//============================================

	//============================================

	//Main Loop to update z by merging algorithm

	//============================================

	//============================================

	int i_loop = 0; 
	double cell_make_check_main = MPI_Wtime();
	for(i_loop = 0; i_loop<n_max_iteration; i_loop++)

	{
		double cell_make_check_iteration = MPI_Wtime();
		//testout

		//RPrint("==============================");

		//RPrint("Main loop of Cell_Make.. i+1: "); RPrint(i_loop+1);

		
		double cell_make_check2 = MPI_Wtime();
		bool b_DEBUG_Zmat = false; 

		if(i_loop == -5) b_DEBUG_Zmat = true;


		//----------
		//clear category matrix for possible garbage
		//Note: z has only positive integer as category #
		//----------
		/*for(int i=0; i<nrow; i++)
		{
			for(int j=0; j<ncol; j++)
			{
				if(fabs_FHDI(z[i][j] < 1e-3)) z[i][j] = 0.0; 
			}
		}
		*/
		

		//-----------

		//var; ref of Drs Im, Kim, and Fuller

		//ml = A_M, actual numbers of rows containing missing units

		//ol = A_R, actual numbers of rows having the fully observed units 

		//-----------

		//uox = sorted unique categorized patterns: e.g., a1b, a2c, c1f, d44, d45, etc.

		//mox = sorted unique categorized patterns: e.g., a00, a01, b10, c00, d01, etc.  

		//-----------

		Zmat_Extension_cpp(z, nrow, ncol, cn, 

								ml, ol, i_count_ol, i_count_ml,  

								uox, mox, i_count_uox, i_count_mox,

								b_DEBUG_Zmat,TestOut);

		//if (i_loop == 1) {
		//	cout << " Cell_Make_check2 at node " << mynode << " at iteration " << i_loop << " = " << MPI_Wtime() - cell_make_check2 << endl;
		//}

		if(i_count_ml <= 0 || i_count_ol <= 0)

		{ 
			//FHDI::RPrint("ERROR! i_count_ml or _ol is zero! Change k, check data quality, further break down categorical variables, or so. It may help \n"); 
			cout<<"ERROR! i_count_ml or _ol is zero! Change k, check data quality, further break down categorical variables, or so. It may help "<<endl; 

			//early deallocaiton -----------------
			delete[] d_k_Collapsible; 
			delete[] cn; 
			delete[] ml;
			delete[] ol;
			delete[] tnU;
			Del_dMatrix(zbase, nrow, ncol);
			Del_dMatrix(uox, nrow, ncol);
			Del_dMatrix(mox, nrow, ncol);
			delete[] i_orn;
			delete[] i_orn_temp;
			delete[] i_orn_temp2;
			
			return 0;
		}

		double cell_make_check3 = MPI_Wtime();
		//----------
		//clear category matrix for possible garbage
		//Note: z has only positive integer as category #
		//----------
		/*for(int i=0; i<nrow; i++)
		{
			for(int j=0; j<ncol; j++)
			{
				if(fabs_FHDI(z[i][j] < 1e-3)) z[i][j] = 0.0; 
			}
		}
		*/


	

		//------------------------------------------

		//generate number of donors nD[]

		// List of observed cells serving as donors List_nU: row numbers of donors per each missing row

		// Table of nU, tnU: total number of all donors for each missing row 

		//------------------------------------------

		//re-initialize tnU and List_nU and v_nD

		List_nU.initialize(i_count_mox);

		v_nD = std::vector<int>(); 

	    tnU = NULL; tnU = new int[i_count_uox]; Fill_iVector(tnU, i_count_uox, 0);

		

		bool b_DEBUG_nDAU = 0; 

		//if(i_loop>9) b_DEBUG_nDAU = 1;

		
		double cell_make_check3_nDAU = MPI_Wtime();
		bool b_success_nDAU = nDAU_cpp_MPI(uox, mox, i_count_uox, i_count_mox, ncol,

				 cn, ol, i_count_ol, i_cellmake,

				 v_nD, List_nU, tnU, b_DEBUG_nDAU); 

		//testout

		//RPrint("nDAU_... has been done");

		if(!b_success_nDAU)
		{			
			//FHDI::RPrint("Error! nDAU Failed! Change k, check data quality, further break down categorical variables, or so. It may help \n");
			cout<<"Error! nDAU Failed! Change k, check data quality, further break down categorical variables, or so. It may help "<<endl;

			//early deallocaiton -----------------
			delete[] d_k_Collapsible; 
			delete[] cn; 
			delete[] ml;
			delete[] ol;
			delete[] tnU;
			Del_dMatrix(zbase, nrow, ncol);
			Del_dMatrix(uox, nrow, ncol);
			Del_dMatrix(mox, nrow, ncol);
			delete[] i_orn;
			delete[] i_orn_temp;
			delete[] i_orn_temp2;
			
			//return 0; //abnormal ending 

			exit(0);
		}


		//================================
		//When there are at least one Non-Collapsible Categorical Variable exists
		//Skip Cell-Collapse procedure 
		//as of 2018, 04 21
		//INACTIVATED 2018, 04 25
		//================================
		/*
		if(i_NonCollapsible_categ_total >= 1)
		{
			RPrint("There is at least one Non-Collapsible categorical variable! \n");
			RPrint("So, the automatic cell-collapse won't take place. \n");
			RPrint("This may affect the following FHDI/FEFI imputation's convergence. \n");
			RPrint("If converged imputation is of top interest, (categorical=NULL) may help. \n");
			break;
			
		}
		*/

		
		
		
		
		//================================

		//Check whether there exist enough fully observed units for current 

		//"k" categories 

		//If not, reduce k to (k-1) categories at a specific position

		//================================

		int i_nD_sum = 0; //summation of all nD[] 

		for(unsigned i=0; i<v_nD.size(); i++) {if(v_nD[i] < 2) i_nD_sum += v_nD[i]; }

		//when too small number of donors
		double cell_make_check3_special = MPI_Wtime();
		if(i_count_uox == 2 && i_nD_sum >0)

		{

			//testout
			//FHDI::RPrint(" Special case for small donors!  \n");
			cout<<" Special case for small donors!  "<<endl;

			//-----

			//max of k[ncol]: maximum category number

			//-----

			//double max_k = 0.0; for(int i=0; i<ncol; i++) {if(max_k <k[i]) max_k = k[i];}

			//double max_k = max_FHDI(d_k, ncol);
			double max_k = 0.0;
			double max_k_Collapsible = 0.0; //max of k among only collapsible variables
			//-------------------
			//when all variables are collapsible 
			//-------------------
			if(i_NonCollapsible_categ_total == 0)
			{
				max_k = max_FHDI(d_k, ncol);
			}
			//---------
			//when there is non-collapsible categorical variable
			//---------
			if(i_NonCollapsible_categ_total >= 1)
			{
				max_k_Collapsible = max_FHDI(d_k_Collapsible, ncol);
			}
			
			//testout
			/*RPrint("within Special case \n"); 
			RPrint("i_NonCollapsible_categ_total: "); RPrint(i_NonCollapsible_categ_total); 
			RPrint("max_k: "); RPrint(max_k); 
			RPrint("d_k[]: "); RPrint(d_k, ncol); 
			RPrint("max_k_Collapsible: "); RPrint(max_k_Collapsible); 
			RPrint("d_k_Collapsible[]: "); RPrint(d_k_Collapsible, ncol); 
			*/
			
			std::vector<int> v_maxk; 

			//which(d_k, ncol, max_k, v_maxk); //ACTUAL locations of columns that have the max category k

			//-------------------
			//when all variables are collapsible 
			//-------------------
			if(i_NonCollapsible_categ_total == 0)
			{
				which(d_k, ncol, max_k, v_maxk); //ACTUAL locations of columns that have the max category k
			}
			//---------
			//when there is non-collapsible categorical variable
			//---------
			if(i_NonCollapsible_categ_total >= 1)
			{
				which(d_k_Collapsible, ncol, max_k_Collapsible, v_maxk); //ACTUAL locations of columns that have the max category k
			}

			
			//-----

			//get some of orn of which location is the same as the maxk 

			//-----

			int n_orm = (int)v_maxk.size(); 

			int* i_orm = new int[n_orm]; 

			for(int j=0; j<n_orm; j++) i_orm[j] = i_orn[v_maxk[j]-1]; //Note: actual loc of column

			

			//-----

			//get the first column that has the min(i_orm)

			//-----

			int min_orm = min_FHDI(i_orm, n_orm); 

			//for(int j=0; j<n_orm; j++) {if(min_orm>i_orm[j]) min_orm = i_orm[j];} 

			std::vector<int> v_orm; 

			which(i_orm, n_orm, min_orm, v_orm); //ACTUAL locations of columns that have the min orm

			int i_mc = v_orm[0]; //the first cell that has the min orn

			
			//testout
			//RPrint("i_mc: "); RPrint(i_mc); 

			//------

			//re-categorize with the reduced category (k-1)

			//------

			double* d_x_temp = new double[nrow]; //a column corresponding to min_orm

			double* d_z_temp = new double[nrow]; //new column with categorized values

			for(int i=0; i<nrow; i++) d_x_temp[i] = x[i][i_mc-1]; //Note: actual loc in i_mc

			double d_k_one_dummy = d_k[i_mc -1];
			
			//b_success_categorize = categorize_cpp(d_x_temp, nrow, d_k[i_mc -1], d_z_temp); 
			
			int NonCollapsible_categorical_1 = NonCollapsible_categorical[i_mc-1];
			b_success_categorize = categorize_cpp(d_x_temp, nrow, d_k_one_dummy, d_z_temp,
			                                      NonCollapsible_categorical_1); 

			//testout
			//RPrint("within Special case, after single column categorize() \n"); 
			//RPrint("NonCollapsible_categorical_1: "); RPrint(NonCollapsible_categorical_1); 
			//RPrint("d_k_one_dummy: "); RPrint(d_k_one_dummy); 
												  
			if(!b_success_categorize) 
			{
				//early deallocaiton -----------------
				delete[] d_k_Collapsible; 
				delete[] cn; 
				delete[] ml;
				delete[] ol;
				delete[] tnU;
				Del_dMatrix(zbase, nrow, ncol);
				Del_dMatrix(uox, nrow, ncol);
				Del_dMatrix(mox, nrow, ncol);
				delete[] i_orn;
				delete[] i_orn_temp;
				delete[] i_orn_temp2;
				
				delete[] i_orm; 
				delete[] d_x_temp; 
				delete[] d_z_temp;
				
				return 0; 
			}

			//if this column is collapsible 
			if(NonCollapsible_categorical[i_mc-1] == 0)
			{
				d_k[i_mc -1] = d_k_one_dummy ; //this k value may have been updated for automatic categorical var.
				d_k_Collapsible[i_mc-1] = d_k_one_dummy; 
			}
			//----------
			//clear category matrix for possible garbage
			//Note: d_z_temp has only positive integer as category #
			//----------
			for(int i=0; i<nrow; i++)
			{ if(fabs_FHDI(d_z_temp[i] < 1e-3)) d_z_temp[i] = 0.0; }
			
			
			
			//-----

			//update zbase's one column with the reduced category

			//-----
			//if this column is collapsible 
			if(NonCollapsible_categorical[i_mc-1] == 0)
			{
				for(int i=0; i<nrow; i++) zbase[i][i_mc-1] = d_z_temp[i]; //Note: actual loc in i_mc
			}

			//testout
			//RPrint("In special case, i_mc : "); RPrint(i_mc); 
			//RPrint("d_z_temp[] : "); RPrint(d_z_temp, nrow); 
			
			//-----
			//check too small category number error (April 2018)
			//-----
			//if this column is collapsible 
			if(NonCollapsible_categorical[i_mc-1] == 0)
			{			
				if( fabs_FHDI(d_k[i_mc-1] - 1) < 1.0)
				{
					{
						//FHDI::RPrint("Error! There is not enough observed units or categories. Change k or break down category; it may help  \n "); 
						cout<<"Error! There is not enough observed units or categories. Change k or break down category; it may help  "<<endl; 

						return 0;
					}
				}
			}
			
			//-----
			//reduce the previous category number
			//for the ease of category condensation
			//-----
			//if this column is collapsible 
			if(NonCollapsible_categorical[i_mc-1] == 0)
			{			
				d_k[i_mc-1] = d_k[i_mc-1] - 1;  
				d_k_Collapsible[i_mc-1] = d_k_Collapsible[i_mc-1] - 1; 
			}

			

			//------

			//check Abort condition by looking at min of k

			//------

			int min_k_new = min_FHDI(d_k, ncol); 

			//for(int j=0; j<ncol; j++) {if(min_k_new>k[j]) min_k_new = k[j];}

			if(min_k_new < 2) 

			{
				//FHDI::RPrint("There is not enough observed units in the original data. Thus, automatic cell-collapse has been done!   \n"); 
			 	cout<<"There is not enough observed units in the original data. Thus, automatic cell-collapse has been done!  "<<endl; 

				break;
			}

			//MUST ACTIVATE BREAK after adding a LOOP !!!!!

			
			//----------
			//clear category matrix for possible garbage
			//Note: zbase has only positive integer as category #
			//----------
			for(int i=0; i<nrow; i++)
			{
				for(int j=0; j<ncol; j++)
				{
					if(fabs_FHDI(zbase[i][j] < 1e-3)) zbase[i][j] = 0.0; 
				}
			}

			

			//-----

			//update with new reduced data

			//-----
			

			Copy_dMatrix(zbase, nrow, ncol, z);

			

			//-----

			//re-initialize before calling Zmat_...()

			//-----

			Fill_iVector(ml, nrow, 0);

			Fill_iVector(ol, nrow, 0);

			Fill_dMatrix(uox, nrow, ncol, 0.0);

			Fill_dMatrix(mox, nrow, ncol, 0.0);

			i_count_ol = 0;

			i_count_ml = 0; 

			i_count_uox = 0;

			i_count_mox = 0;

			v_nD = std::vector<int>(); 

			//----------
			//clear category matrix for possible garbage
			//Note: z has only positive integer as category #
			//----------
			/*for(int i=0; i<nrow; i++)
			{
				for(int j=0; j<ncol; j++)
				{
					if(fabs_FHDI(z[i][j] < 1e-3)) z[i][j] = 0.0; 
				}
			}*/
			

			Zmat_Extension_cpp(z, nrow, ncol, cn, 

								ml, ol, i_count_ol, i_count_ml,  

								uox, mox, i_count_uox, i_count_mox,

								b_DEBUG_Zmat,TestOut);

		

			if(i_count_ml <= 0 || i_count_ol <= 0)
			{ 
				//FHDI::RPrint("ERROR! i_count_ml or _ol is zero!   \n");
				cout<<"ERROR! i_count_ml or _ol is zero!   "<<endl;

				//FHDI::RPrint("Change k, further break down categorical variables, or check data quality \n");		

				//early deallocaiton -----------------
				delete[] d_k_Collapsible; 
				delete[] cn; 
				delete[] ml;
				delete[] ol;
				delete[] tnU;
				Del_dMatrix(zbase, nrow, ncol);
				Del_dMatrix(uox, nrow, ncol);
				Del_dMatrix(mox, nrow, ncol);
				delete[] i_orn;
				delete[] i_orn_temp;
				delete[] i_orn_temp2;
				
				delete[] i_orm; 
				delete[] d_x_temp; 
				delete[] d_z_temp;
				
				return 0;
			}

			

			//------------------------------------------

			//generate number of donors nD[]

			// List of observed cells serving as donors List_nU

			// Table of nU tnU

			//NOTE: for large data with many variables, 

			//      The case of no possible donors may arise 

			//      If so, reduction of k begins, say with (k-1)

			//------------------------------------------

			//re-initialize tnU and List_nU and v_nD

			List_nU.initialize(i_count_mox);

			v_nD = std::vector<int>(); 

			tnU = NULL; tnU = new int[i_count_uox]; Fill_iVector(tnU, i_count_uox, 0);

			

			b_success_nDAU =  nDAU_cpp_MPI(uox, mox, i_count_uox, i_count_mox, ncol,

					 cn, ol, i_count_ol, i_cellmake,

					 v_nD, List_nU, tnU, b_DEBUG_nDAU); 			

			if(!b_success_nDAU)
			{			
				//FHDI::RPrint("Error! nDAU Failed! Change k, check data quality, further break down categorical variables, or so. It may help \n");
                cout<<"Error! nDAU Failed! Change k, check data quality, further break down categorical variables, or so. It may help "<<endl;
				//early deallocaiton -----------------
				delete[] d_k_Collapsible; 
				delete[] cn; 
				delete[] ml;
				delete[] ol;
				delete[] tnU;
				Del_dMatrix(zbase, nrow, ncol);
				Del_dMatrix(uox, nrow, ncol);
				Del_dMatrix(mox, nrow, ncol);
				delete[] i_orn;
				delete[] i_orn_temp;
				delete[] i_orn_temp2;
				
				delete[] i_orm; 
				delete[] d_x_temp; 
				delete[] d_z_temp;
				
				//return 0; //abnormal ending 	

				exit(0);
			}
			

			//-----

			//local deallocation

			//-----

			v_maxk = std::vector<int>();

			v_orm  = std::vector<int>();

			delete[] i_orm; 

			delete[] d_x_temp; 

			delete[] d_z_temp; 

		} //end of Special case of too small donors 
		//if (i_loop == 1) {
		//	cout << " Cell_make_check3_special at iteration " << i_loop << " at node " << mynode << " = " << MPI_Wtime() - cell_make_check3_special << endl;
		//}

		//-----

		//find the first location having min of nD[], i.e. minimum donors

		//-----

		int i_min_nD = min_FHDI(v_nD); 

		int i_reci = 0;

		for(int i=0; i<(int)v_nD.size();i++) 

		{if(v_nD[i] == i_min_nD){ i_reci = i; break;}} 

		
		//if (i_loop == 1) {
		//	cout << " Cell_Make_check3 at node " << mynode << " at iteration " << i_loop << " = " << MPI_Wtime() - cell_make_check3 << endl;
		//}

		

		



		//===========================

		//===========================

		//Merge z 

		//===========================

		//===========================

		//testout 

		bool b_DEBUG_Merge = false; 

		if(i_loop ==-3) b_DEBUG_Merge = true; 

		double cell_make_check4 = MPI_Wtime();
		//----------
		//clear category matrix for possible garbage
		//Note: z has only positive integer as category #
		//----------
		/*for(int i=0; i<nrow; i++)
		{
			for(int j=0; j<ncol; j++)
			{
				if(fabs_FHDI(z[i][j] < 1e-3)) z[i][j] = 0.0; 
			}
		}*/
		//------------------
		//when there is at least one non-collapsible categorical variable
		//skip the merge procedure 
		//2018, 0426
		//------------------
		if(i_NonCollapsible_categ_total>=1 && v_nD[i_reci] < 2)
		{
			//FHDI::RPrint("The current data set does not have enough donors while there is at least one non-collapsible categorical variable! \n");
			//FHDI::RPrint("Thus, auto merging procedure won't take place! \n"); 
			
			cout<<"The current data set does not have enough donors while there is at least one non-collapsible categorical variable! "<<endl;
			cout<<"Thus, auto merging procedure won't take place! "<<endl;; 
			
			break; 
		}
		
		if(v_nD[i_reci] < 2) //if donors are less than 2, do MERGE

		{

			Merge_Extension_cpp(i_reci, uox, i_count_uox, 

							mox, i_count_mox, tnU, 

							cn,  ol, i_count_ol,

							z, nrow, ncol, 

							i_merge,

							b_DEBUG_Merge);

		}

		if(v_nD[i_reci] >=2) //if more than 2 donors exist, exit

		{ break; } //finish main loop 

		//cout << " Cell_Make_check4 at node " << mynode << " at iteration " << i_loop << " = " << MPI_Wtime() - cell_make_check4 << endl;


		

		//----

		//unconverged ending

		//-----

		if(i_loop == n_max_iteration-1) 

		{ 
			//FHDI::RPrint(" reached n_max_iteration after step ");

			//FHDI::RPrint("%d ", n_max_iteration);
			
			//FHDI::RPrint(" Change k, check data quality, further break down categorical variables, or so. It may help ");

			cout<<" reached n_max_iteration after step "<<endl;

			cout<< n_max_iteration<<endl;
			
			cout<<" Change k, check data quality, further break down categorical variables, or so. It may help "<<endl;

			//early deallocaiton -----------------
			delete[] d_k_Collapsible; 
			delete[] cn; 
			delete[] ml;
			delete[] ol;
			delete[] tnU;
			Del_dMatrix(zbase, nrow, ncol);
			Del_dMatrix(uox, nrow, ncol);
			Del_dMatrix(mox, nrow, ncol);
			delete[] i_orn;
			delete[] i_orn_temp;
			delete[] i_orn_temp2;
				

			return 0; 

		}


	} //end of main loop

	//testout
	//cout<<"i_count_ol: "<< i_count_ol <<endl;
	//cout <<"i_count_ml: " << i_count_ml << endl;
	//cout << "i_count_uox: " << i_count_uox << endl;
	//cout << "i_count_mox: " << i_count_mox << endl;
	//FHDI::RPrint("converged in Cell_Make after iterations: "); FHDI::RPrint("%d ", i_loop+1);
	//cout << " Cell_Make_check_main at node " << mynode << " = " << MPI_Wtime() - cell_make_check_main << endl;


	//----------
	//clear uox and mox matrix for possible garbage
	//Note: Must have only positive integer as category #
	//----------
	for(int i=0; i<i_count_uox; i++)
	{
		for(int j=0; j<ncol; j++)
		{
			if(fabs_FHDI(uox[i][j] < 1e-3)) uox[i][j] = 0.0; 
		}
	}
	for(int i=0; i<i_count_mox; i++)
	{
		for(int j=0; j<ncol; j++)
		{
			if(fabs_FHDI(mox[i][j] < 1e-3)) mox[i][j] = 0.0; 
		}
	}




	//--------------

	//prepare separate output of Cell Make

	//--------------

	double* d_temp_um = new double[ncol]; 

	for(int i=0; i<i_count_uox; i++)

	{

		for(int j=0; j<ncol; j++) d_temp_um[j] = uox[i][j]; 

		

		rbind_uox_CellMake.append_block(d_temp_um); 

	}

	for(int i=0; i<i_count_mox; i++)

	{

		for(int j=0; j<ncol; j++) d_temp_um[j] = mox[i][j]; 

		

		rbind_mox_CellMake.append_block(d_temp_um); 

	}	

	delete[] d_temp_um; 

	

	

	//testout

	//RPrint(" ========= Cell_Make_Extension.. has successfully finished!");

	//FHDI::RPrint(" ========= FHDI_CellMake has successfully finished!\n");
	//cout<<" ========= FHDI_CellMake has successfully finished!"<<endl;

	

	//-------------------------------------

	//Deallocation

	//-------------------------------------
	delete[] d_k_Collapsible; 
	
	delete[] cn; 

	delete[] i_orn;

	delete[] i_orn_temp;

	delete[] i_orn_temp2;

	delete[] ml;

	delete[] ol;

	delete[] tnU;



	Del_dMatrix(zbase, nrow, ncol);

	Del_dMatrix(uox, nrow, ncol);

	Del_dMatrix(mox, nrow, ncol);

	

	

	return 1;

}



//} //end of namespace
