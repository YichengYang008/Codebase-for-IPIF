//Fn===========================================================================

//Cell_Make_Neighbor_cpp.cc-----------------------------------------------------------------------------

//Fn===========================================================================
#include "KNN.cc"
	bool Cell_Make_Neighbor_cpp(double** x, const int nrow, const int ncol, double* d_k,

		int* NonCollapsible_categorical,

		double** z,

		rbind_FHDI &rbind_uox_CellMake,

		rbind_FHDI &rbind_mox_CellMake,

		List_FHDI &List_nU,

		const int i_merge, ofstream& TestOut)

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

		//IN    : int NonCollapsible_categorical(nrol) = {0,0, .., 1,.. 0} 
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

		const int n_max_iteration = nrow * 2; //maximum number of iterations 

											  //const int n_max_iteration = 2;  //temporary

											  //-------------------------------------
											  //Determine if there is Non-Collapsible Categorical variable
											  //-------------------------------------
		double* d_k_Collapsible = new double[ncol]; //k for collapsible variables only 
		Copy_dVector(d_k, ncol, d_k_Collapsible);

		int i_NonCollapsible_categ_total = 0;
		for (int i = 0; i<ncol; i++)
		{
			i_NonCollapsible_categ_total += NonCollapsible_categorical[i];

			if (NonCollapsible_categorical[i] == 0) d_k_Collapsible[i] = d_k[i]; //use user-defined k 
			if (NonCollapsible_categorical[i] == 1) d_k_Collapsible[i] = 1; //will be overwritten by actual total categories
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

		//TestOut << " z matrix from categorize" << endl;
		//RPrint(z, nrow, ncol, TestOut);

		if (!b_success_categorize)
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

		//double** zbase = New_dMatrix(nrow, ncol);

		//Copy_dMatrix(z, nrow, ncol, zbase);



		//-------------------------------------

		//sort in the order of high missing rate

		//at the end, i_orn has the "actual" column numbers from highest missing rate

		//            to the lowest missing rate

		//-------------------------------------

		int* i_orn = new int[ncol]; 	 	Fill_iVector(i_orn, ncol, 0);

		int* i_orn_temp = new int[ncol]; 	Fill_iVector(i_orn_temp, ncol, 0);

		int* i_orn_temp2 = new int[ncol]; 	Fill_iVector(i_orn_temp2, ncol, 0);



		int i_temp = 0;

		for (int i_col = 0; i_col<ncol; i_col++)

		{

			i_temp = 0;

			for (int i_row = 0; i_row<nrow; i_row++)

			{

				if (fabs_FHDI(z[i_row][i_col]) < 1e-5) //count only  "0"

				{
					i_temp++;
				}

			}

			i_orn_temp[i_col] = i_temp; //store how many "0" in this column

		}

		Copy_iVector(i_orn_temp, ncol, i_orn_temp2); //store before sorting 

													 //std::sort(&i_orn_temp[0], &i_orn_temp[ncol-1]); //this works well, but not recommended

		std::sort(i_orn_temp, i_orn_temp + ncol);



		for (int i = 0; i<ncol; i++)

		{

			i_temp = i_orn_temp[ncol - 1 - i]; //reversed searching since the "sort" occurred in ascending order

			for (int j = 0; j<ncol; j++)

			{

				if (i_temp == i_orn_temp2[j])

				{

					i_orn[i] = j + 1; //store column number (actual number, 1, 2, ...) 

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

		//List_FHDI List_nU(nrow); //default for the size of nrow, but will be updated in the main loop

		int* tnU = new int[nrow]; Fill_iVector(tnU, nrow, 0); //this default size will be udpated in the main loop

		int i_cellmake = 2; // Inactiavte the b_success_nDAU because cell make with KNN always have enough donors

		//============================================

		//============================================

		//Main part to update z by K-nearest-neighbor

		//============================================

		//============================================

		bool b_DEBUG_Zmat = false;

		//if (i_loop == -5) b_DEBUG_Zmat = true;

		Zmat_Extension_cpp(z, nrow, ncol, cn,

			ml, ol, i_count_ol, i_count_ml,

			uox, mox, i_count_uox, i_count_mox,

			b_DEBUG_Zmat, TestOut);


		//cout << "Yicheng Check cell make i_count_uox: " << i_count_uox << ", i_count_mox: " << i_count_mox << endl;
		if (i_count_ml <= 0 || i_count_ol <= 0)

		{
			cout << "ERROR! i_count_ml or _ol is zero! Change k, check data quality, further break down categorical variables, or so. It may help " << endl;

			//early deallocaiton -----------------
			delete[] d_k_Collapsible;
			delete[] cn;
			delete[] ml;
			delete[] ol;
			delete[] tnU;
			//Del_dMatrix(zbase, nrow, ncol);
			Del_dMatrix(uox, nrow, ncol);
			Del_dMatrix(mox, nrow, ncol);
			delete[] i_orn;
			delete[] i_orn_temp;
			delete[] i_orn_temp2;

			return 0;
		}


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



		bool b_success_nDAU = nDAU_cpp_MPI(uox, mox, i_count_uox, i_count_mox, ncol,

			cn, ol, i_count_ol, i_cellmake,

			v_nD, List_nU, tnU, b_DEBUG_nDAU);

		//testout

		//RPrint("nDAU_... has been done");

		if (!b_success_nDAU)
		{
			cout << "Error! nDAU Failed! Change k, check data quality, further break down categorical variables, or so. It may help " << endl;

			//early deallocaiton -----------------
			delete[] d_k_Collapsible;
			delete[] cn;
			delete[] ml;
			delete[] ol;
			delete[] tnU;
			//Del_dMatrix(zbase, nrow, ncol);
			Del_dMatrix(uox, nrow, ncol);
			Del_dMatrix(mox, nrow, ncol);
			delete[] i_orn;
			delete[] i_orn_temp;
			delete[] i_orn_temp2;

			//return 0; //abnormal ending 	

			exit(0);
		}
		
		/*TestOut<<"v_nD before KNN"<<endl;
		for (int j1 = 0; j1 < i_count_mox; j1++) {
			TestOut<<"v_nD["<<j1<<"]: "<< v_nD[j1]<<endl;
		}

		TestOut << "mox in neighbor with i_count_mox = " << i_count_mox << endl;

		for (int kk2 = 0; kk2 < i_count_mox; kk2++) {
			for (int kk3 = 0; kk3 < ncol; kk3++) {
				TestOut << setw(20) << mox[kk2][kk3];
			}
			TestOut << endl;
		}
		TestOut << "uox in neighbor with i_count_uox = " << i_count_uox << endl;

		for (int kk2 = 0; kk2 < i_count_uox; kk2++) {
			for (int kk3 = 0; kk3 < ncol; kk3++) {
				TestOut << setw(20) << uox[kk2][kk3];
			}
			TestOut << endl;
		}*/

		//cout<<"List_nU before cell make neighbor"<<endl;
		//List_nU.print_List_FHDI_yicheng(TestOut);

		for (int i_loop = 0; i_loop < i_count_mox; i_loop++) {

			if (v_nD[i_loop] < 2) {
			
				    KNN(i_loop, uox, i_count_uox,
					mox, i_count_mox, d_k,
					cn, ol, i_count_ol, 
					nrow, ncol, i_merge,
					v_nD, List_nU, TestOut);

			}
		}

		//TestOut << "List_nU after cell make neighbor" << endl;
		//List_nU.print_List_FHDI_yicheng(TestOut);

		//TestOut << "v_nD after KNN" << endl;
		//for (int j1 = 0; j1 < i_count_mox; j1++) {
		//	TestOut << "v_nD[" << j1 << "]: " << v_nD[j1] << endl;
		//}


		//TestOut << "List_nU before sort" << endl;
		//List_nU.print_List_FHDI_yicheng(TestOut);


		//Important!!! Note that List_nU must include actual locations of donors in uox in ascending orders
		//Or there will be mismatch problem in FHDI_Neighbor to compute fractional weights

		std::vector<int> List_temp;

		for (int j2 = 0; j2 < i_count_mox;j2++) {

			List_temp.clear();

			List_nU.get_block_yicheng(j2, List_temp);

			sort(List_temp.begin(),List_temp.end());

			//cout<<"List_temp at "<<j2<<endl;
			//for (int j3 = 0; j3 < List_temp.size();j3++) {
			//	cout<<"List_temp["<<j3<<"]: "<< List_temp[j3]<<endl;
			//}

			List_nU.put_block(j2, List_temp);

		}
		//int i_min_nD = min_FHDI(v_nD);
		//TestOut << "List_nU after KNN" << endl;
		//List_nU.print_List_FHDI_yicheng(TestOut);

		for(int kk=0; kk< v_nD.size();kk++){
		if (v_nD[kk] < 2) TestOut<<"ERROR! The "<<kk<< "th data set after k-nearest-neighbor does not gurantee at least two donors for all unique missing patterns"<<endl;
		}

		//----------
		//clear uox and mox matrix for possible garbage
		//Note: Must have only positive integer as category #
		//----------
		for (int i = 0; i<i_count_uox; i++)
		{
			for (int j = 0; j<ncol; j++)
			{
				if (fabs_FHDI(uox[i][j] < 1e-3)) uox[i][j] = 0.0;
			}
		}
		for (int i = 0; i<i_count_mox; i++)
		{
			for (int j = 0; j<ncol; j++)
			{
				if (fabs_FHDI(mox[i][j] < 1e-3)) mox[i][j] = 0.0;
			}
		}




		//--------------

		//prepare separate output of Cell Make

		//--------------

		double* d_temp_um = new double[ncol];

		for (int i = 0; i<i_count_uox; i++)

		{

			for (int j = 0; j<ncol; j++) d_temp_um[j] = uox[i][j];



			rbind_uox_CellMake.append_block(d_temp_um);

		}

		for (int i = 0; i<i_count_mox; i++)

		{

			for (int j = 0; j<ncol; j++) d_temp_um[j] = mox[i][j];



			rbind_mox_CellMake.append_block(d_temp_um);

		}

		delete[] d_temp_um;





		//testout

		//RPrint(" ========= Cell_Make_Extension.. has successfully finished!");

		//RPrint(" ========= Cell_Make with KNN has successfully finished!\n", TestOut);



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



		//Del_dMatrix(zbase, nrow, ncol);

		Del_dMatrix(uox, nrow, ncol);

		Del_dMatrix(mox, nrow, ncol);





		return 1;

	}
