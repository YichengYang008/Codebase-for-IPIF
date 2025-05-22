//Fn===========================================================================

//nDAU_cpp.cc-----------------------------------------------------------------------------

//Fn===========================================================================

//namespace FHDI{

bool nDAU_Bigp_ultra_cpp_MPI(const int memory, const int nrow_uox, const int nrow_mox, const int nrow, const int ncol, const int i_collapsing, const int i_SIS_type,

	const int top, int** correlation_ranking_top, double** correlation_top, int i_cellmake,

	MPI_File fh_uox, MPI_File fh_mox, List_FHDI &uox_infor,

	const int i_merge, double* d_k,

	std::vector<int> &v_nD, List_FHDI &List_nU, int** codes, ofstream& TestOut)

	//Description=========================================

	// identify information of the missing cells and observed cells

	//

	// Algorithm:  

	// 

	//

	// original R code: Dr. Im, J. and Dr. Kim, J. 

	// c++ code: 		Dr. Cho, I. 

	// All rights reserved

	// 

	// updated: March 28, 2017

	//----------------------------------------------------

	//IN    : int nrow_uox          = number of rows of uox
	//IN    : int nrow_mox          = number of rows of mox
	//IN    : int nrow              = number of rows of raw data
	//IN    : int ncol              = number of columns of raw data  
	//IN    : int i_option_collapsing = choice of big-p algorithm 
	//                              0 = no big-p algorithms
	//                             !0 = perform big-p algorithms
	//IN    : int i_SIS_type            = type of SIS method
	//IN    : int top                   = number of top correlation or correlation ranking
	//IN    : int** correlation_ranking_top = top correlation ranking of all variables
	//IN    : int** correlation_top = top correlation of all variables
	//IN    : i_cellmake            = cell construction method: 1: cell collapsing; 2: cell makw with KNN
	//IN	: MPI_File fh_uox	    = file string to read uox
	//IN	: MPI_File fh_mox	    = file string to read mox
	//IN    : List_FHDI uox_infor   = Actual index list of uox in z
	//IN    : List_FHDI mox_infor   = Actual index list of mox in z
	//IN    : int i_merge           = random number controller
	//IN    : int* d_k              = category of all variables
	//INOUT : std::vector v_nD		= total number of donnors of each missing pattern
	//INOUT : List_FHDI List_nU     = list of observed cells to serve as donors 
	//INOUT : int** codes(nrow)     = storage to record most correlated variables of all mox

	//IMPORTANT SUMMARY OF INDEX 
	//Matrix           starting offset
	//v_oloc                 1
	//List_nU(entity)        1
	//====================================================

{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//if (mynode == 2) {
	//	cout << " correlation_top matrix in nDAU_bigp_ultra at node " << mynode << " where top is " << top << endl;
	//	for (int kk2 = 0; kk2 < ncol; kk2++) {
	//		for (int kk3 = 0; kk3 < top; kk3++) {
	//			cout << setw(20) << correlation_top[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//	cout << " correlation_ranking_top matrix in nDAU_bigp_ultra at node " << mynode << " where top is " << top << endl;
	//	for (int kk2 = 0; kk2 < ncol; kk2++) {
	//		for (int kk3 = 0; kk3 < top; kk3++) {
	//			cout << setw(20) << correlation_ranking_top[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//}
	//cout << "nDAU nrow_uox is " << nrow_uox << " and nrow_mox is " << nrow_mox << " at node " << mynode << endl;
	//============================================================================================
	//============================================================================================
	//============================================================================================

	//nDAU: Find existing donors for all mox and add to List_nu_temp and v_nD_temp on slave processors   

	//============================================================================================
	//============================================================================================
	//============================================================================================
	double nDAU_begin = MPI_Wtime();

	int number_before_split = 0;
	int number_after_split = 0;
	int startpoint = 0;
	int endpoint = 0;
	int splitnode = 0;
	int L_temp = 0;
	//cout<<"nrow_mox at mynode "<<mynode<<" is "<< nrow_mox <<endl;
	//cout << "nrow_uox at mynode " << mynode << " is " << nrow_uox << endl;

	if (nrow_mox % (totalnodes - 1) != 0) { splitnode = nrow_mox % (totalnodes - 1); }
	if (nrow_mox % (totalnodes - 1) == 0) { splitnode = 1; }
	//if (mynode == 0) cout << "Split: " << splitnode << endl;
	number_after_split = floor(1.0*nrow_mox / (1.0*totalnodes - 1));
	number_before_split = 1.0*(nrow_mox - floor(1.0*nrow_mox / (1.0*totalnodes - 1)) *(totalnodes - splitnode - 1)) / splitnode;

	if (mynode >= 1 && mynode <= splitnode) {
		startpoint = (mynode - 1) * number_before_split;
		endpoint = mynode *number_before_split;
		L_temp = number_before_split;
	}
	if (mynode > splitnode) {
		startpoint = splitnode * number_before_split + (mynode - splitnode - 1)*number_after_split;
		endpoint = splitnode * number_before_split + (mynode - splitnode)*number_after_split;
		L_temp = number_after_split;
	}
	if ((number_before_split*splitnode + number_after_split* (totalnodes - splitnode - 1)) != nrow_mox) {
		cout << "Work Assignment Error!!!!" << endl;
		return 0;
	}

	//cout << "Splitnode is " << splitnode << "; number_before_split is " << number_before_split << " and number_after_split is " << number_after_split << "; Mynode " << mynode << ", " << startpoint << "<= x < " << endpoint << endl;
	//cout << "L_temp distributed is " << L_temp << " at node " << mynode << endl;
	//========================================================
	//Read distributed mox matrix on all slave processors
	//=========================================================

	double*  d_read_out = NULL;
	double** mox_temp = NULL;

	if (mynode != 0) {
		d_read_out = new double[ncol];
		mox_temp = New_dMatrix(L_temp, ncol);
	}

	int counter = 0;
	for (int k = startpoint; k < endpoint; k++) {

		MPI_In_uox_mox(ncol, k, fh_mox, d_read_out);

		for (int m = 0; m < ncol; m++) {
			mox_temp[counter][m] = d_read_out[m];
		}
		counter++;
	}

	if (mynode != 0) { delete[] d_read_out; }

	//if (mynode == 3) {
	//	cout << "mox Reading matrix from node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < L_temp; kk2++) {
	//		for (int kk3 = 0; kk3 < ncol; kk3++) {
	//			cout << setw(20) << mox_temp[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//}

	//=====================================================================
	//Find selected variables for distributed moxs on all slave processors
	//=====================================================================
	std::vector<int> v_mxl_temp; // 1. temporary vector for the locaiton of non zeros in mox[i]
	int* mox_buffer = new int[ncol]; //temporary integer array to store mox[i]
	int** codes_Send = NULL;// The codes at each node 

	if(mynode!=0) codes_Send = New_iMatrix(L_temp, i_collapsing);

	for (int t = 0; t < L_temp; t++) {

		int oc = 0; // Get number of observed values in current mox

		for (int k = 0; k < ncol; k++) {

			mox_buffer[k] = (int)mox_temp[t][k];

			if (mox_buffer[k] > 0) {

				oc++;
			}
		}

		//cout << "nDAU mox[" << i << "] at node " << mynode << endl;
		//for (int m = 0; m < ncol; m++) {
		//	cout << setw(20) << mox_buffer[m];
		//}
		//cout << endl;

		v_mxl_temp.clear();

		if (oc < (i_collapsing + 1)) { // oc <= i_option_collapsing

			whichINV(mox_buffer, ncol, 0, v_mxl_temp); //get the location of Non-zero in mox 
		}

		if (oc > i_collapsing) {

			if (i_SIS_type == 1) {

				//correlated_variable_intersection2(ncol, i_collapsing, i, ia_temp, correlation_yicheng, correlation_ranking, v_mxl, TestOut);
				correlated_variable_intersection_ultra(ncol, i_collapsing, top, t, mox_buffer, correlation_top, correlation_ranking_top, v_mxl_temp, TestOut);

			}

			if (i_SIS_type == 2) {

				//correlated_variable_union(ncol, i_collapsing, i, ia_temp, correlation_yicheng, correlation_ranking, v_mxl, TestOut);
				correlated_variable_union_ultra(ncol, i_collapsing, top, t, mox_buffer, correlation_top, correlation_ranking_top, v_mxl_temp, TestOut);
			}

			if (i_SIS_type == 3) {

				//correlated_variable_global(ncol, i_collapsing, ia_temp, correlation_yicheng, v_mxl, TestOut);
				correlated_variable_global_ultra(ncol, i_collapsing, top, mox_buffer, correlation_top, correlation_ranking_top, v_mxl_temp, TestOut);

			}

		}

		for (int b3 = 0; b3 < v_mxl_temp.size(); b3++) {

			codes_Send[t][b3] = v_mxl_temp[b3];
			//cout << "v_mxl_Top[" << b3 << "]: " << v_mxl[b3];
			//cout<<"codes["<<i<<"]["<<b3<<"]: "<< v_mxl[b3] <<" at node "<<mynode<<endl;
		}

	}



	//if (mynode == 3) {
	//	cout<<"Codes at node "<<mynode<<endl;
	//	for (int t = 0; t < L_temp;t++) {
	//		for (int h = 0;h < i_collapsing;h++) {
	//			cout << setw(20) << codes_Send[t][h];
	//		}
	//		cout << endl;
	//	}
	//}

	delete[] mox_buffer;



	//Initialize 
	std::vector<int> v_nD_temp(L_temp);
	List_FHDI List_nU_temp(L_temp); // list of donors in uox for distributed recipients


	std::string s_cn0; //temporary string for a row of mox
	int* ia_temp = new int[ncol]; //temporary integer array  

	//Initilization
	int recursion_size = 0;
	recursion_size = (int)floor((memory * 0.35 * 1e+9) / (ncol * 8));
	if (nrow_uox < recursion_size) recursion_size = nrow_uox;

	int recursion_time = 0; // total number of times to import partial uox matrix

	if (nrow_uox % recursion_size == 0) recursion_time = (1.0*nrow_uox) / recursion_size;
	if (nrow_uox % recursion_size != 0) recursion_time = (int)floor(1.0*nrow_uox / recursion_size) + 1;
	int recursion_size_temp = 0;

	int recursion_size_last = nrow_uox - (recursion_time - 1)*recursion_size;
	if (recursion_size_last < 10) cout<<"ERROR!!! recursion_size_last is less than 10 in nDAU_Bigp_ultra!!!"<<endl;

	if (mynode == 0) recursion_time = 0;
	//if (mynode == 0 || mynode == 1) cout << "nDAU Bigp recursion_time is " << recursion_time << " and recursion_size_last is " << recursion_size_last << " and recursion_size is " << recursion_size << " at node " << mynode << endl;

	//============================================
	//Read recursion_size rows of uox recursively 
	//============================================

	for (int i_recur = 0; i_recur < recursion_time; i_recur++) {
		if (i_recur != (recursion_time - 1)) recursion_size_temp = recursion_size;
		if (i_recur == (recursion_time - 1)) recursion_size_temp = recursion_size_last;
		//cout << "recursion_size_temp is " << recursion_size_temp << " at i = " << i << endl;

		double** uox_recursion = NULL;
		double* array_temp = NULL;

		//boundary index of recursively read z matrix
		startpoint = 0;
		endpoint = 0;

		if (mynode != 0) {

			uox_recursion = New_dMatrix(recursion_size_temp, ncol);
			array_temp = new double[ncol];

			if (i_recur != (recursion_time - 1)) {
				startpoint = i_recur*recursion_size;
				endpoint = (i_recur + 1)*recursion_size;
			}

			if (i_recur == (recursion_time - 1)) {
				startpoint = i_recur*recursion_size;
				endpoint = i_recur*recursion_size + recursion_size_last;
			}
		}

		//cout << "recusive startpoint: " << startpoint << "; recusive endpoint: " << endpoint << " at recursion " << i_recur <<" at node "<<mynode<<endl;

		counter = 0;
		for (int k = startpoint; k < endpoint; k++) {

			MPI_In_uox_mox(ncol, k, fh_uox, array_temp);

			for (int m = 0; m < ncol; m++) {
				uox_recursion[counter][m] = array_temp[m];
			}
			counter++;
		}

		//if (mynode == 1) {
		//	cout << "uox from recusion " << i_recur << " at node " << mynode << endl;
		//	for (int kk2 = 0; kk2 < recursion_size_temp; kk2++) {
		//		for (int kk3 = 0; kk3 < ncol; kk3++) {
		//			cout << setw(20) << uox_recursion[kk2][kk3];
		//		}
		//		cout << endl;
		//	}
		//}

		//===============================================================================
		//Compare resursively read uoxs with distributed mox and update List_nU and v_nD
		//===============================================================================
		for (int i = 0; i < L_temp; i++) {

			//-----

			//find non zero cells of current missing row mox

			//-----

			for (int k = 0; k < ncol; k++) ia_temp[k] = (int)mox_temp[i][k];

			//Note: below will contain ACTUAL location of cells with non-zero observed data

			std::vector<int> v_mxl; //temporary vector for the locaiton of non zeros in mox

			v_mxl.clear();

			//get the location of selected variables in mox 
			for (int t = 0; t < i_collapsing; t++) {

				if (codes_Send[i][t] != 0) {

					v_mxl.push_back(codes_Send[i][t]);

				}

			}

			//-----
			//total number of non zeros in mox
			//-----

			const int nxl = v_mxl.size(); //total number of observed units at current row (>0)

			double* d_rcn0_temp = new double[nxl]; //temporary

			for (int k = 0; k < nxl; k++) d_rcn0_temp[k] = mox_temp[i][v_mxl[k] - 1]; //"-1" is for actual location


			//-----
			//condense non-zero category names in current row of mox
			//-----

			std::string s_rcn0;

			Trans1(d_rcn0_temp, nxl, s_rcn0); //condense one row  


											  //deallocate used array or matrix (2018_0416)
			delete[] d_rcn0_temp;

			//--------------

			//make a list of possible donors 

			//--------------

			std::vector<int> v_oloc;

			//----

			//get all non-zero cells from all observed rows

			//----

			double * d_t1 = new double[v_mxl.size()];

			std::vector<std::string> v_cand; //vector of found string with condensed non-zero observed data

											 //-----

											 //search all observed rows

											 //of which the same columns are non-zero as the current missing row  

											 //-----

			for (int m = 0; m < recursion_size_temp; m++)

			{	//Note: "-1" in v_mxl is from the ACTUAL location info in it

				for (unsigned k = 0; k < v_mxl.size(); k++) d_t1[k] = uox_recursion[m][v_mxl[k] - 1];

				//------

				//condense the found rows with non-zero observed cell only  

				//------

				std::string s_cand_1;

				Trans1(d_t1, v_mxl.size(), s_cand_1);

				v_cand.push_back(s_cand_1); //add more string to the string vector 

			}
			//cout << "Debug_nDAu 4 at node " << mynode << endl;


			//--------------

			//Find the rows of v_cand that match the current non-zero missing pattern

			//Note: below will contain ACTUAL locations of the found rows 

			//--------------

			//Note that it must return global actual location in uox, not local ones
			which(v_cand, s_rcn0, i_recur, recursion_size, v_oloc); //get the actual global local locations of observed cells containing s_rcn0 


																	//------------

																	//local deallocation

																	//------------

			delete[] d_t1;

			//----------------

			//Store oloc into LIST named nU 

			// ith row of nU corresponds to ith row of the List

			//----------------

			int i_oloc_temp = (int)v_oloc.size(); //size of current oloc
												  //cout<<"i_oloc_temp["<<i<<"] at node "<<mynode<<" is "<< i_oloc_temp <<endl;

			if (i_oloc_temp > 0) //only for meaningful oloc

			{
				std::vector<double> d_oloc_temp(v_oloc.begin(), v_oloc.end());//put_block_yicheng only accepts double vector
				List_nU_temp.put_block_yicheng(i, i_oloc_temp, d_oloc_temp);//add new list
			}

			//---------------------

			//number of donors; this case only the matched observed rows become possible donors

			//---------------------

			int i_temp_tocn_sum = 0;
			int uox_size_temp = 0;

			for (int k = 0; k < i_oloc_temp; k++) //accumulate all possible donors 

			{
				uox_size_temp = 0;
				uox_infor.get_a_row_size(v_oloc[k] - 1, uox_size_temp);

				i_temp_tocn_sum += uox_size_temp;
				
				//i_temp_tocn_sum += uox_info_final[v_oloc[k] - 1][1];//second column of uox_info_final has count of all uoxs
			} //-1 for actual loc

			  //v_nD_temp.push_back(i_temp_tocn_sum); //store into integer vector to return 
			int value_temp = 0;
			value_temp = v_nD_temp[i] + i_temp_tocn_sum;

			v_nD_temp[i] = value_temp;


		}//end of main loop

		if (mynode != 0) {
			delete[] array_temp;
			Del_dMatrix(uox_recursion, recursion_size_temp, ncol);
		}

	}//end of recursion

	 //Deallocation
	delete[] ia_temp;

	MPI_Barrier(MPI_COMM_WORLD);

	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_make at node " << mynode << " = " << MPI_Wtime() - cell_make_begin << endl;
	//	printf("Yicheng Running time of nDAU ultra Bigp = %f seconds\n", MPI_Wtime() - nDAU_begin);
	//}

	//if (mynode == 1) {
	//	for (int j = 0; j < v_nD_temp.size(); j++) {
	//		cout << "v_nD_temp[" << j << "]: " << v_nD_temp[j] <<" at node "<<mynode<< endl;
	//	}

	//	cout << "Before List_nU at node " << mynode << endl;
	//	List_nU_temp.print_List_FHDI_yicheng();
	//}

	//============================================================================================
	//============================================================================================
	//============================================================================================

	//KNN: Find deficient donors for mox who has less than two donors and update List_nu_temp and 
	//     v_nD_temp on slave processors

	//============================================================================================
	//============================================================================================
	//============================================================================================
	double KNN_begin = MPI_Wtime();


	//Initialization
	std::vector<int> f_nD_temp; // locations of v_nD_temp who has less than two donors

	if (mynode != 0) {
		for (int m = 0; m < L_temp; m++) {
			if (v_nD_temp[m] < 2) {
				f_nD_temp.push_back(m);
			}
		}
	}

	int L_temp_size = 0;
	if (mynode != 0) L_temp_size = f_nD_temp.size();
	//cout << "L_temp_size is " << L_temp_size << " at node " << mynode << endl;

	List_FHDI min_donor_temp(L_temp_size); //list of uox who has the min distance to moxs
	List_FHDI second_min_donor_temp(L_temp_size); //list of uox who has the second min distance to moxs
	std::vector<double> d_min_temp(L_temp_size); //min dist from uoxs to moxs who have less than two donors
	std::vector<double> d_second_min_temp(L_temp_size); //second min dist from uoxs to moxs who have less than two donors

	if (mynode != 0) {
		for (int k2 = 0; k2 < L_temp_size; k2++) {
			d_min_temp[k2] = 123456789.0;
			d_second_min_temp[k2] = 123456789.0;
		}
	}

	double temp_min = 0.0;//to fetch temp min dist
	double temp_second_min = 0.0;//to fetch temp second min dist


	for (int j_recur = 0; j_recur < recursion_time; j_recur++) {

		if (j_recur != (recursion_time - 1)) recursion_size_temp = recursion_size;
		if (j_recur == (recursion_time - 1)) recursion_size_temp = recursion_size_last;
		//cout << "recursion_size_temp is " << recursion_size_temp << " at j_recur = " << j_recur << endl;

		double** uox_recursion2 = NULL;
		double* array_temp2 = NULL;

		//boundary index of recursively read z matrix
		startpoint = 0;
		endpoint = 0;

		if (mynode != 0) {

			uox_recursion2 = New_dMatrix(recursion_size_temp, ncol);
			array_temp2 = new double[ncol];

			if (j_recur != (recursion_time - 1)) {
				startpoint = j_recur*recursion_size;
				endpoint = (j_recur + 1)*recursion_size;
			}

			if (j_recur == (recursion_time - 1)) {
				startpoint = j_recur*recursion_size;
				endpoint = j_recur*recursion_size + recursion_size_last;
			}
		}

		//cout << "KNN recusive startpoint: " << startpoint << "; recusive endpoint: " << endpoint << " at recursion " << j_recur <<" at node "<<mynode<<endl;

		counter = 0;
		for (int k = startpoint; k < endpoint; k++) {

			MPI_In_uox_mox(ncol, k, fh_uox, array_temp2);

			for (int m = 0; m < ncol; m++) {
				uox_recursion2[counter][m] = array_temp2[m];
			}
			counter++;
		}

		//===============================================================================
		//Compare resursively read uoxs with distributed mox who has less than two donors 
		//and compute minimum and second minimum Euclidean distance for these specified mox
		//===============================================================================

		for (int l = 0; l < L_temp_size; l++) {

			int i_reci = 0;
			i_reci = f_nD_temp[l];

			if (i_merge == 1) std::srand(time(NULL)); //turn on random seed using C++ standard rand() fn 
													  //this will generate purely random numbers and 
													  //will be platform-dependent, i.e., different on Win and Linux 
			if (i_merge == 0) std::srand(123);	//turn on the fixed seed //This is still platform-dependent 
												//maybe, use Numerical Recipe for platform-independent 
												//-----------------
												//Which columns are NOT missing in mox[i_reci][]
												//-----------------
			double* d_mox_row = new double[ncol]; //temporary array
			for (int i = 0; i < ncol; i++) d_mox_row[i] = mox_temp[i_reci][i];
			std::vector<int> v_mxl; //ACTUAL location of non-missing column of mox[i_reci][]

			for (int t = 0; t < i_collapsing; t++) {

				if (codes_Send[i_reci][t] != 0) {

					v_mxl.push_back(codes_Send[i_reci][t]);

				}

			}

			const int i_nxl = (int)v_mxl.size(); //number of non-missing cell on this row
			delete[] d_mox_row;

			//-----------------
			//Find the nearest potential donor cells using "fdis"
			//NOTE: below two matrix and array has nrow_uox rows since it is 
			//related to observed cells uox
			//-----------------`
			double ** d_cand = New_dMatrix(recursion_size_temp, ncol); //NOTE: the column may be flexible for below cases 
			double *  d_fdist = new double[recursion_size_temp];        //distance between entities 
			Fill_dVector(d_fdist, recursion_size_temp, 0.0);

			if (i_nxl == 1) //when the current missing row has only ONE observed cell   
			{
				//------------
				//make a copy of all rows of the one column 
				// that corresponds to the column where the observed cell of current missing row
				// is located 		
				//------------
				for (int i = 0; i < recursion_size_temp; i++)
				{
					d_cand[i][0] = (uox_recursion2[i][v_mxl[0] - 1]) / (d_k[v_mxl[0] - 1]);
				} //-1 for ACTUAL location 

				  //calculate distance using |a-b|^2
				const double d_mox_mxl = (mox_temp[i_reci][v_mxl[0] - 1]) / (d_k[v_mxl[0] - 1]);
				distance2(d_cand, recursion_size_temp, i_nxl, d_mox_mxl,
					d_fdist);
			}//end of i_nxl

			if (i_nxl > 1) //when current missing row has more than one column that has observed cells 
			{
				//------------
				//make a copy of all rows of all columns that correspond to the observed cells 
				//------------
				for (int i = 0; i < recursion_size_temp; i++)
				{
					for (int j = 0; j < i_nxl; j++) //note: i_nxl is the length of v_mxl
					{
						d_cand[i][j] = (uox_recursion2[i][v_mxl[j] - 1]) / (d_k[v_mxl[j] - 1]); // Normalized it by k
					} //-1 for ACTUAL location 
				}
				//-------------
				//calculate distance = sum(|a-b|^2) per row where mox[i][mxl] is the origin
				//-------------	
				double d_sum_dist = 0.0;
				for (int i = 0; i < recursion_size_temp; i++)
				{
					d_sum_dist = 0.0; //re-initialize
					for (int j = 0; j < i_nxl; j++)
					{
						double d_mox_temp = (mox_temp[i_reci][v_mxl[j] - 1]) / (d_k[v_mxl[j] - 1]);// Normalized it by k
						double d_temp1 = d_cand[i][j];
						d_sum_dist += (d_mox_temp - d_temp1)*(d_mox_temp - d_temp1);
					}
					d_fdist[i] = d_sum_dist;
				}

			}//end of i_nxl

			 //-------------
			 //set the distance from obatined donors in uox to mox[i_reci] as 1234567 instead of 0s
			 //if the distnace is 0, that means this uox is already a donor in the list
			 //-------------
			for (int k1 = 0; k1 < recursion_size_temp; k1++) {
				if (d_fdist[k1] == 0) {
					//if(mynode==1) cout<<"Distance at "<<k1<<" is uox at j_recur =" << j_recur << " and l = " << l << " at node " << mynode << endl;
					d_fdist[k1] = 1234567.0;
				}
			}

			//for (int p = 0; p < recursion_size_temp; p++) cout << "d_fdist[" << p << "]: " << d_fdist[p] << "at j_recur = " << j_recur << " and l = " << l << " at node " << mynode << endl;

			//------------
			//find the minimum distance
			//------------
			std::vector<int> v_floc; //ACTUAL location of donors in uox who has the minimum distance
			double d_min_fdist = 0.0;
			if (i_nxl >= 1)
			{
				d_min_fdist = min_FHDI(d_fdist, recursion_size_temp);
				which(d_fdist, recursion_size_temp, d_min_fdist, j_recur, recursion_size, v_floc);
			}
			//if (mynode == 1) cout<<"d_min_fdist is "<< d_min_fdist <<" at j_recur="<< j_recur <<" and l="<< l <<" at node "<<mynode<<endl;
			//------------------------------------------------------------------
			//update min and second min dist and uox list for 
			//disrtibuted mox regarding recursively imported uoxs
			//--------------------------------------------------------------

			//--------------------------------------------
			//Update min dist and corresponding uox list
			//--------------------------------------------
			temp_min = d_min_temp[l];
			//if (mynode == 1)  cout << "d_min_temp[" << l << "]: " << temp_min << "; d_min_fdist["<<l<<"]: "<< d_min_fdist << " at j_recur=" << j_recur << " and l=" << l << endl;
			//Case1: If new min dist > former min dist of resursive uox, skip
			//if ((d_min_fdist - temp_min) > 1e-15) {
				//if (mynode == 1) cout << "Case1 d_min_fdist=" << d_min_fdist << "; temp_min=" << temp_min << " at j_recur = " << j_recur << " and l = " << l << " at node " << mynode << endl;
				//continue;
			//}

			int v_floc_size = v_floc.size();
			//Important!!! Note that if the current min dist has been replaced, this min dist (i.e., temp_min) may be used as second min dist
			//Thus, one has to save temp_min and its corresponding locations in v_floc_temp for later use.
			//Example for Case2:
			//1 recursion: d_min = 0.22; d_second_min = 0.33
			//2 recursion: d_min = 0.11; d_second_min = 0.25
			//2 final    : d_min = 0.11; d_second_min = 0.22
			std::vector<int> v_floc_temp;

			//Case2: If new min dist < former min dist of resursive uox, change min dist and remove former list and add new list
			if ((d_min_fdist - temp_min) < -1e-15) {
				min_donor_temp.get_block_yicheng(l, v_floc_temp);//save former locations for later use for d_second_min
				//if (mynode == 1) cout << "Case2 d_min_fdist=" << d_min_fdist << "; temp_min=" << temp_min << " at j_recur = " << j_recur << " and l = " << l << " at node " << mynode << endl;
				d_min_temp[l] = d_min_fdist;//change min dist
				//if (mynode == 1) cout << "d_min_temp[" << l << "]: " << d_min_temp[l] << "; d_min_fdist: " << d_min_fdist << " at j_recur=" << j_recur << " and l=" << l << " at node " << mynode << endl;
				min_donor_temp.remove_block_yicheng(l);//remove former list
				std::vector<double> v_floc_double(v_floc.begin(), v_floc.end());//put_block_yicheng only accepts double vector
				min_donor_temp.put_block_yicheng(l, v_floc_size, v_floc_double);//add new list
			}

			//Case3: If new min dist = former min dist of resursive uox, keepon adding to former list
			if (fabs(d_min_fdist - temp_min) < 1e-15) {
				//if (mynode == 1) cout << "Case3 d_min_fdist=" << d_min_fdist << "; temp_min=" << temp_min << " at j_recur = " << j_recur << " and l = " << l << " at node " << mynode << endl;
				std::vector<double> v_floc_double(v_floc.begin(), v_floc.end());//put_block_yicheng only accepts double vector
				min_donor_temp.put_block_yicheng(l, v_floc_size, v_floc_double);//add new list
			}


			//----------------------------------------------
			//Find the global second min dist until current recusion
			//---------------------------------------------

			double d_second_min_fdist = 0.0; //global second min dist
			std::vector<int> v_floc_second; //ACTUAL location of donors in uox who has the second minimum distance

			temp_second_min = d_second_min_temp[l];//former global second min dist

			//In the first recursion read, d_second_min_temp is initilized with 1234567.0
			//No need to compare current d_second_min_fdist with former 1234567.0
			//if (j_recur == 0)
			//{
			//	d_second_min_fdist = second_min_FHDI(d_fdist, recursion_size_temp);
			//	which(d_fdist, recursion_size_temp, d_second_min_fdist, j_recur, recursion_size, v_floc_second);
			//}

			//Format
			//1 recuriosn:  temp_min      tem_second_min
			//2 recursion:  d_min_fdist   d_second_min_fdist

			//-----------------------------------------------------------------------------------------------------------------------------------
			//If min dist replacement does not happen and temp_min < d_min, one only need compare current d_min_fdist with former temp_second_min
			//-----------------------------------------------------------------------------------------------------------------------------------
			//Example 1:
			//1 recursion: d_min(temp_min) = 0.22; d_second_min (temp_second_min) = 0.32
			//2 recursion: d_min = 0.33; d_second_min = 0.45
			//2 final    : d_min = 0.22; d_second_min = 0.32

			//Example 2:
			//1 recursion: d_min(temp_min) = 0.22; d_second_min (temp_second_min) = 0.44
			//2 recursion: d_min = 0.33; d_second_min = 0.36
			//2 final    : d_min = 0.22; d_second_min = 0.33

			//Example 3:
			//1 recursion: d_min(temp_min) = 0.22; d_second_min (temp_second_min) = 0.33
			//2 recursion: d_min = 0.33; d_second_min = 0.36
			//2 final    : d_min = 0.22; d_second_min = 0.33 (temp_second_min)

			if ((d_min_fdist - temp_min) > 1e-15)
			{
				//d_second_min_fdist = second_min_FHDI(d_fdist, recursion_size_temp);

				//if (temp_second_min - d_second_min_fdist < -1e-15) {
				//	d_second_min_fdist = temp_second_min;
				//}
				if ((d_min_fdist - temp_second_min) <  -1e-15) d_second_min_fdist = d_min_fdist;
				if ((d_min_fdist - temp_second_min) >= -1e-15) d_second_min_fdist = temp_second_min;

				which(d_fdist, recursion_size_temp, d_second_min_fdist, j_recur, recursion_size, v_floc_second);
			}

			//If min dist replacement does not happen and temp_min = d_min, one only need compare current d_min_fdist with former temp_second_min
			//Example 4:
			//1 recursion: d_min(temp_min) = 0.22; d_second_min (temp_second_min) = 0.37
			//2 recursion: d_min = 0.22; d_second_min = 0.36
			//2 final    : d_min = 0.22; d_second_min = 0.36

			if (fabs(d_min_fdist - temp_min) < 1e-15)
			{
				d_second_min_fdist = second_min_FHDI(d_fdist, recursion_size_temp);

				if (temp_second_min - d_second_min_fdist <= 1e-15) {
					d_second_min_fdist = temp_second_min;
				}

				which(d_fdist, recursion_size_temp, d_second_min_fdist, j_recur, recursion_size, v_floc_second);
			}

			//------------------------------------------------------------------------------------------------------
			//If min dist replacement happens, one only need compare current d_second_min_fdist with former temp_min
			//------------------------------------------------------------------------------------------------------
			//Example 1:
			//1 recursion: d_min(temp_min) = 0.22; d_second_min = 0.33
			//2 recursion: d_min = 0.11; d_second_min = 0.25
			//2 final    : d_min = 0.11; d_second_min = 0.22

			//However, sometimes d_second_min_fdist may not take its place
			//Example 2:
			//1 recursion: d_min(temp_min) = 0.55; d_second_min = 0.66
			//2 recursion: d_min = 0.22; d_second_min = 0.33
			//2 final    : d_min = 0.22; d_second_min = 0.33

			//Example 3:
			//1 recursion: d_min(temp_min) = 0.55; d_second_min = 0.66
			//2 recursion: d_min = 0.22; d_second_min = 0.55
			//2 final    : d_min = 0.22; d_second_min = 0.55 (temp_min)

			//Hence, one need to check
			if ((d_min_fdist - temp_min) < -1e-15)
			{
				d_second_min_fdist = second_min_FHDI(d_fdist, recursion_size_temp);
				if (temp_min - d_second_min_fdist <= 1e-15) {
					d_second_min_fdist = temp_min;
				}
				which(d_fdist, recursion_size_temp, d_second_min_fdist, j_recur, recursion_size, v_floc_second);
			}

			//Note that global second min dist may be same with current d_min_fdist. Min dist replacement must not happen
			//And temp_second_min < current d_second_min
			//Example:
			//1 recursion: d_min(temp_min) = 0.11; d_second_min = 0.22
			//2 recursion: d_min = 0.22; d_second_min = 0.33
			//2 final    : d_min = 0.11; d_second_min = 0.22
			//if (d_second_min_fdist == d_min_fdist) { cout << "Error in KNN ultra of cell make for getting the second minimum distance!!!" << endl;}

			//if (mynode == 1 && l == 0) {
			//	cout << "YEYE d_min_fdist is " << d_min_fdist << " and d_second_min_fdist is " << d_second_min_fdist << endl;
			//	for (int m = 0; m < recursion_size_temp; m++) cout << "YEYE d_fdist[" << m << "]: " << d_fdist[m] << endl;
			//}
			//cout<<"d_min_fdist is "<< d_min_fdist <<" and d_second_min_fdist is "<< d_second_min_fdist <<" at j_recur=" << j_recur << " and l=" << l << endl;

			//--------------------------------------------
			//Update second min dist and corresponding uox list
			//--------------------------------------------

			//if (mynode == 1)  cout << "d_min_temp[" << l << "]: " << temp_min << "; d_min_fdist["<<l<<"]: "<< d_min_fdist << " at j_recur=" << j_recur << " and l=" << l << endl;
			//Case1: If new min dist > former min dist of resursive uox, then it is an error. Impossible!!!
			if ((d_second_min_fdist - temp_second_min) > 1e-15) {
				cout << "ERROR in KNN ultra function to find second min distance!!!!" << endl;
				//continue;
			}


			//If temp_min < d_min_fdist (replacement not happen)
			//1. d_second_min_fdist = d_min_fdist                      (Case 5)
			//2. d_second_min_fdist = temp_second_min                  (Case 4)

			//If temp_min = d_min_fdist (replacement not happen)
			//3. d_second_min_fdist = second_min_FHDI()                (Case 6)
			//4. d_second_min_fdist = temp_second_min                  (Case 4)

			//If temp_min > d_min_fdist (replacement happens)
			//5. d_second_min_fdist = second_min_FHDI()                (Case 3)
			//6. d_second_min_fdist = temp_min                         (Case 2)


			int v_floc_second_size = v_floc_second.size();
			int v_floc_temp_size = v_floc_temp.size();

			//Case2: d_second_min_fdist == temp_min
			//change second min dist and remove former list and add new locations in two steps
			if (fabs(d_second_min_fdist - temp_min)< 1e-15) {
				//if (mynode == 1) cout << "Case2 d_second_min_fdist=" << d_second_min_fdist << "; temp_second_min=" << temp_second_min << " at j_recur = " << j_recur << " and l = " << l << " at node " << mynode << endl;
				d_second_min_temp[l] = d_second_min_fdist;//change min dist
														  //if (mynode == 1) cout << "d_min_temp[" << l << "]: " << d_min_temp[l] << "; d_min_fdist: " << d_min_fdist << " at j_recur=" << j_recur << " and l=" << l << " at node " << mynode << endl;
				second_min_donor_temp.remove_block_yicheng(l);//remove former list

				//Firstly add accumulated locations from min dist
				std::vector<double> v_floc_temp_double(v_floc_temp.begin(), v_floc_temp.end());//put_block_yicheng only accepts double vector
				second_min_donor_temp.put_block_yicheng(l, v_floc_temp_size, v_floc_temp_double);

				//Secondly add locations from second min dist
				std::vector<double> v_floc_second_double(v_floc_second.begin(), v_floc_second.end());//put_block_yicheng only accepts double vector
				second_min_donor_temp.put_block_yicheng(l, v_floc_second_size, v_floc_second_double);//add new list
			}

			//Case3: If Min dist replacement happens and temp_min > d_second_min_fdist
			//change second min dist and remove former list and only add locations in second min dist
			if ((d_min_fdist - temp_min < -1e-15) && (fabs(d_second_min_fdist - temp_min) > 1e-15)) {
				//if (mynode == 1) cout << "Case2 d_second_min_fdist=" << d_second_min_fdist << "; temp_second_min=" << temp_second_min << " at j_recur = " << j_recur << " and l = " << l << " at node " << mynode << endl;
				d_second_min_temp[l] = d_second_min_fdist;//change min dist
														  //if (mynode == 1) cout << "d_min_temp[" << l << "]: " << d_min_temp[l] << "; d_min_fdist: " << d_min_fdist << " at j_recur=" << j_recur << " and l=" << l << " at node " << mynode << endl;
				second_min_donor_temp.remove_block_yicheng(l);//remove former list

				std::vector<double> v_floc_second_double(v_floc_second.begin(), v_floc_second.end());//put_block_yicheng only accepts double vector
				second_min_donor_temp.put_block_yicheng(l, v_floc_second_size, v_floc_second_double);//add new list
			}

			//Case4: d_second_min_fdist == temp_second_min, keep on adding to former list
			if (fabs(d_second_min_fdist - temp_second_min) < 1e-15) {
				//if (mynode == 1) cout << "Case3 d_second_min_fdist=" << d_second_min_fdist << "; temp_second_min=" << temp_second_min << " at j_recur = " << j_recur << " and l = " << l << " at node " << mynode << endl;
				std::vector<double> v_floc_second_double(v_floc_second.begin(), v_floc_second.end());//put_block_yicheng only accepts double vector
				second_min_donor_temp.put_block_yicheng(l, v_floc_second_size, v_floc_second_double);//add new list
			}

			//Case5: d_second_min_fdist == d_min_fdist, remove and add new ones
			if ((d_min_fdist - temp_min > 1e-15) && (fabs(d_second_min_fdist - temp_second_min) > 1e-15)) {
				d_second_min_temp[l] = d_second_min_fdist;//change min dist

				second_min_donor_temp.remove_block_yicheng(l);//remove former list

				std::vector<double> v_floc_second_double(v_floc_second.begin(), v_floc_second.end());//put_block_yicheng only accepts double vector
				second_min_donor_temp.put_block_yicheng(l, v_floc_second_size, v_floc_second_double);//add new list
			}

			//Case6: If Min dist replacement does not happen and d_min_fdist=temp_min, remove and add new ones
			if ((fabs(d_min_fdist - temp_min) < 1e-15) && (fabs(d_second_min_fdist - temp_second_min) > 1e-15)) {
				d_second_min_temp[l] = d_second_min_fdist;//change min dist
														  //if (mynode == 1) cout << "d_min_temp[" << l << "]: " << d_min_temp[l] << "; d_min_fdist: " << d_min_fdist << " at j_recur=" << j_recur << " and l=" << l << " at node " << mynode << endl;
				second_min_donor_temp.remove_block_yicheng(l);//remove former list

				std::vector<double> v_floc_second_double(v_floc_second.begin(), v_floc_second.end());//put_block_yicheng only accepts double vector
				second_min_donor_temp.put_block_yicheng(l, v_floc_second_size, v_floc_second_double);//add new list
			}

			//Deallocation
			Del_dMatrix(d_cand, recursion_size_temp, ncol);
			delete[] d_fdist;

		}//end of main loop

		//if (mynode == 1) {
		//	//for (int j1 = 0; j1 < L_temp_size; j1++) cout << "d_second_min_temp[" << j1 << "]: " << d_second_min_temp[j1] << " at j_recur " << j_recur << endl;
		//	cout << "d_min_temp[12]: " << d_min_temp[12] << " at j_recur " << j_recur << endl;
		//	cout << "d_second_min_temp[12]: " << d_second_min_temp[12] << " at j_recur " << j_recur << endl;
		//	cout << "second_min_donor_temp at j_recur " << j_recur << endl;
		//	second_min_donor_temp.print_List_FHDI_yicheng();
		//}

		if (mynode != 0) {
			delete[] array_temp2;
			Del_dMatrix(uox_recursion2, recursion_size_temp, ncol);
		}

	}//end of recuriosn

	 //TestOut
	//if (mynode == 3) {
	//	//for (int j1 = 0; j1 < L_temp_size; j1++) cout << "final d_min_temp[" << j1 << "]: " << d_min_temp[j1] <<" at node "<<mynode<< endl;
	//	//for (int j1 = 0; j1 < L_temp_size; j1++) cout << "final d_second_min_temp[" << j1 << "]: " << d_second_min_temp[j1] << " at node " << mynode << endl;
	//	cout << "final min_donor_temp at node " << mynode << endl;
	//	min_donor_temp.print_List_FHDI_yicheng();
	//	cout << "final second_min_donor_temp at node " << mynode << endl;
	//	second_min_donor_temp.print_List_FHDI_yicheng();
	//}

	//------------------------------------------------------------------------------------------
	//Update List_nU_temp and v_nD_temp of mox who has less than two donors on slave processors
	//regarding minimum and second minimum Euclidean distance
	//------------------------------------------------------------------------------------------
	int i_size_floc = 0;
	for (int t = 0; t < L_temp_size; t++) {
		int i_mox = 0;
		i_mox = f_nD_temp[t];

		//select out a table of the location information of the minimum distance cells
		//------------
		std::vector<int> v_floc; //ACTUAL location of donors in uox who has the minimum distance
		min_donor_temp.get_block_yicheng(t, v_floc);
		i_size_floc = (int)v_floc.size();
		//if (mynode == 1) {
		//	for (int p = 0; p < i_size_floc; p++) cout << "v_floc[" << p << "]: "<< v_floc[p]<<endl;
		//}

		if (i_size_floc <= 0) { RPrint("Error! floc size is 0!"); return 0; }

		int* i_nf = new int[i_size_floc]; // occurance of all donors in uox for mox[i]
		
		int uox_size_temp2 = 0;
		for (int i = 0; i < i_size_floc; i++) {
			uox_size_temp2 = 0;
			uox_infor.get_a_row_size(v_floc[i] - 1, uox_size_temp2);
			i_nf[i] = uox_size_temp2;
		}
		//for (int i = 0; i<i_size_floc; i++) i_nf[i] = uox_info_final[v_floc[i] - 1][1];//the second column of uox_info_final is occurance
		
		
		const int max_nf = max_FHDI(i_nf, i_size_floc); //highest occueance of all donors

		//-------------
		//find rows that have max nf
		//-------------
		std::vector<int> v_nf_max;
		which(i_nf, i_size_floc, max_nf, v_nf_max); //Actual locations which have max of nf
		const int i_size_nf_max = (int)v_nf_max.size();

		//-------------
		//locations having the minimum distance between missing and observed cells
		//-------------
		std::vector<int> v_xloc;// actual locations of donors in uox who has minimum distance and highest occurance
		for (int i = 0; i < i_size_nf_max; i++) v_xloc.push_back(v_floc[v_nf_max[i] - 1]); //-1 for actual loc
		const int i_size_xloc = (int)v_xloc.size();

		// Example:
		// v_floc = {1,3,7,9} -> occur = {3,3,1,1}
		// v_xloc ={1,3} -> occur = {3,3}
		// i_size_xloc = i_size_nf_max =2

		//---------------------------------------
		//Case 0: if mox[i] has one donor and it only needs one more donor from uox who has the smallest distance

		//a) i_size_xloc >=2, need randomly select 1
		//{3,3,2,1} i_size_xloc=2
		//{1,1} i_size_xloc=2

		//b) i_size_xloc == 1, no random selection
		//{3,2,1,1} i_size_xloc=1
		//{2} ..
		//{1} .. 

		//----------------------------------------

		if (v_nD_temp[i_mox] == 1) {
			//cout << "CASE0: KNN with v_nD_temp[" << i_mox << "] with 1 donor and max_nf =" << max_nf << endl;
			//-------------
			//random number within [1, i_size_xloc]
			//Note: this is ACTUAL location
			//-------------
			int i_loc_rand_temp0 = 1;
			if ((i_merge == 1) && (i_size_xloc >= 2)) i_loc_rand_temp0 = std::rand() % i_size_xloc + 1; //purely random 


			const int i_loc_rand_xloc = v_xloc[i_loc_rand_temp0 - 1]; //-1 for actual loc

			v_nD_temp[i_mox] = 1 + max_nf;

			std::vector<double> d_index_temp;//put_block_yicheng only accepts double vector
			d_index_temp.push_back((double)i_loc_rand_xloc);
			List_nU_temp.put_block_yicheng(i_mox, 1, d_index_temp);//add new list

		}//end of v_nD_temp

		if (v_nD_temp[i_mox] == 0) {
			//-----------------------------
			//Case 1: if mox[i] has no donors and the highest occurance (i.e., max_nf) of donors with the smallest distance is >= 2
			//        then it needs only one more donor from uox who has the smallest distance

			//a) max_nf >= 2
			// {3,3,2,1} i_size_xloc >=2, need randomly select 1

			// {3,2,1,1} i_size_xloc ==1,no random selection
			// {2} ..

			if (max_nf >= 2) {
				//cout << "CASE1: KNN with v_nD_temp[" << i_mox << "] with no donor and max_nf = " << max_nf << endl;
				//-------------
				//random number within [1, i_size_xloc]
				//Note: this is ACTUAL location
				//-------------
				int i_loc_rand_temp0 = 1;
				if ((i_merge == 1) && (i_size_xloc >= 2)) i_loc_rand_temp0 = std::rand() % i_size_xloc + 1; //purely random 


				const int i_loc_rand_xloc = v_xloc[i_loc_rand_temp0 - 1]; //-1 for actual loc

				v_nD_temp[i_mox] = max_nf;

				std::vector<double> d_index_temp;//put_block_yicheng only accepts double vector
				d_index_temp.push_back((double)i_loc_rand_xloc);
				List_nU_temp.put_block_yicheng(i_mox, 1, d_index_temp);//add new list

			}


			//Case 2: if mox[i] has no donors and the highest occurance (i.e., max_nf) of donors with the smallest distance is < 2 and 
			//        v_floc.size() is >= 2,
			//        then it needs two donors from uox who has the smallest distance

			// {1,1,1,1} need randomly select 2

			// {1,1} no random selection 

			if ((max_nf < 2) && (v_floc.size() >= 2)) {
				//cout << "CASE2: KNN with v_nD_temp[" << i_mox << "] with no donor and max_nf = " << max_nf << ", and v_floc size is " << v_floc.size() << endl;
				//-------------
				//random number within [1, i_size_xloc]
				//Note: this is ACTUAL location
				//-------------
				int i_loc_rand_temp0 = 1;
				int i_loc_rand_temp1 = 2;

				if ((i_merge == 1) && (i_size_xloc > 2)) {
					i_loc_rand_temp0 = 0;// Make sure it will go into the while loop
					i_loc_rand_temp1 = 0;
					while (i_loc_rand_temp0 == i_loc_rand_temp1) {
						i_loc_rand_temp0 = std::rand() % i_size_xloc + 1;
						i_loc_rand_temp1 = std::rand() % i_size_xloc + 1;
					}
				}


				if (i_loc_rand_temp0 == i_loc_rand_temp1) {
					cout << "Error in KNN of cell make for random selection!!! i_loc_rand_temp0 = " << i_loc_rand_temp0 << "; i_loc_rand_temp1 = " << i_loc_rand_temp1 << endl;
					return 0;
				}
				//if(b_random) i_loc_rand_temp0 = 1 ; //for debugging // 
				const int i_loc_rand_xloc = v_xloc[i_loc_rand_temp0 - 1]; //-1 for actual loc
				const int i_loc_rand_xloc1 = v_xloc[i_loc_rand_temp1 - 1]; //-1 for actual loc

				v_nD_temp[i_mox] = 2;

				std::vector<double> d_index_temp;//put_block_yicheng only accepts double vector
				d_index_temp.push_back((double)i_loc_rand_xloc);
				d_index_temp.push_back((double)i_loc_rand_xloc1);
				List_nU_temp.put_block_yicheng(i_mox, 2, d_index_temp);//add new list

			}



			if ((max_nf < 2) && (v_floc.size() == 1)) {
				
				std::vector<int> v_floc_second; //ACTUAL location of donors in uox who has the second minimum distance
				second_min_donor_temp.get_block_yicheng(t, v_floc_second);

				const int i_size_floc_second = (int)v_floc_second.size();

				//cout << "CASE3: KNN with v_nD_temp[" << i_mox << "] and t = "<<t<< " with no donor and max_nf = " << max_nf << ", and v_floc size is " << v_floc.size() << " and v_floc_second size is " << i_size_floc_second << endl;

				int* i_nf_second = new int[i_size_floc_second]; // occurance of all donors in uox for mox[i]
				
				int uox_size_temp3 = 0;
				for (int i = 0; i < i_size_floc_second; i++) {
					uox_size_temp3 = 0;
					uox_infor.get_a_row_size(v_floc_second[i] - 1, uox_size_temp3);
					i_nf_second[i] = uox_size_temp3;
				}
				//for (int i = 0; i < i_size_floc_second; i++) i_nf_second[i] = uox_info_final[v_floc_second[i] - 1][1]; //-1 for actual loc
				
				const int max_nf_second = max_FHDI(i_nf_second, i_size_floc_second); //highest occueance of all donors

				//-------------
				//find rows that have max_nf_second
				//-------------
				std::vector<int> v_nf_max_second;
				which(i_nf_second, i_size_floc_second, max_nf_second, v_nf_max_second); //Actual locations which have max of nf
				const int i_size_nf_max_second = (int)v_nf_max_second.size();

				//-------------
				//locations having the second minimum distance between missing and observed cells
				//-------------
				std::vector<int> v_xloc_second;// actual locations of donors in uox who has the seond minimum distance and highest occurance
				for (int i = 0; i < i_size_nf_max_second; i++) v_xloc_second.push_back(v_floc_second[v_nf_max_second[i] - 1]); //-1 for actual loc
				const int i_size_xloc_second = (int)v_xloc_second.size();

				//for (int m = 0;m < i_size_floc_second;m++) cout << "v_floc_second[" << m << "]: " << v_floc_second[m] << endl;
				//for (int m = 0;m < i_size_xloc_second;m++) cout << "v_xloc_second[" << m << "]: " << v_xloc_second[m] << endl;

				v_nD_temp[i_mox] = 1 + max_nf_second;

				int i_loc_rand_temp0 = 1;
				if ((i_merge == 1) && (i_size_xloc_second >= 2)) i_loc_rand_temp0 = std::rand() % i_size_xloc + 1; //purely random 
				const int i_loc_rand_xloc = v_xloc_second[i_loc_rand_temp0 - 1]; //-1 for actual loc
				//cout<<"CASE3: i_loc_rand_xloc is "<< i_loc_rand_xloc <<endl;

				std::vector<double> d_index_temp;//put_block_yicheng only accepts double vector
				d_index_temp.push_back((double)v_floc[0]);
				d_index_temp.push_back((double)i_loc_rand_xloc);
				List_nU_temp.put_block_yicheng(i_mox, 2, d_index_temp);//add new list

				delete[] i_nf_second;
			}



		}//end of v_nD_temp

		//Deallocation

		delete[] i_nf;


		}//end of mox who has less than two donors


	 //------------------------------------------------------------------
	 //Send List_nU_temp and v_nD_temp and codes_Send from slave processors to master
	 //-----------------------------------------------------------------

	 //Send size of List_nU_temp firstly 
	if (mynode != 0) {
		std::vector<int> List_nU_temp_size;

		int i_size_temp = 0;
		for (int k = 0; k < L_temp; k++) {
			i_size_temp = 0;
			List_nU_temp.get_a_row_size(k, i_size_temp);
			List_nU_temp_size.push_back(i_size_temp);
		}

		if (List_nU_temp_size.size() != L_temp) {
			cout << "ERROR! List_nU_temp_size is incorrect!" << endl;
		}

		MPI_Send(&List_nU_temp_size[0], List_nU_temp_size.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&v_nD_temp[0], v_nD_temp.size(), MPI_INT, 0, 2, MPI_COMM_WORLD);
		MPI_Send(codes_Send[0], (L_temp*i_collapsing), MPI_INT, 0, 3, MPI_COMM_WORLD);
	}

	std::vector<int> List_nU_size_total;

	//------------------------------------------------------------------------
	//All_gather size of List_nU_temp and v_nD_temp from slave processors on master 
	//------------------------------------------------------------------------
	if (mynode == 0) {
		//Initialization
		//List_nU
		std::vector<int> List_nU_size_before(number_before_split);
		std::vector<int> List_nU_size_after(number_after_split);

		//v_nD
		std::vector<int> v_nD_temp_recv_before(number_before_split);
		std::vector<int> v_nD_temp_recv_after(number_after_split);

		//codes
		int** codes_Recv_before = New_iMatrix(number_before_split, i_collapsing);
		int** codes_Recv_after = New_iMatrix(number_after_split, i_collapsing);

		int counter6 = 0;
		for (int j = 1; j < totalnodes; j = j + 1) {

			if (j >= 1 && j <= splitnode) {
				//List_nU
				MPI_Recv(&List_nU_size_before[0], number_before_split, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				for (int k = 0; k < number_before_split; k++) {
					List_nU_size_total.push_back(List_nU_size_before[k]);
				}

				//v_nD
				MPI_Recv(&v_nD_temp_recv_before[0], number_before_split, MPI_INT, j, 2, MPI_COMM_WORLD, &status);
				for (int k = 0; k < number_before_split;k++) {
					v_nD.push_back(v_nD_temp_recv_before[k]);
				}

				//codes
				MPI_Recv(codes_Recv_before[0], (number_before_split*i_collapsing), MPI_INT, j, 3, MPI_COMM_WORLD, &status);
				for (int k3 = 0; k3 < number_before_split; k3++) {
					for (int k4 = 0; k4 < i_collapsing; k4++) {
						codes[counter6][k4] = codes_Recv_before[k3][k4];
					}
					counter6++;
				}

			}
			//--------------------------------------
			if (j > splitnode) {
				//List_nU
				MPI_Recv(&List_nU_size_after[0], number_after_split, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				for (int k = 0; k < number_after_split; k++) {
					List_nU_size_total.push_back(List_nU_size_after[k]);
				}

				//v_nD
				MPI_Recv(&v_nD_temp_recv_after[0], number_after_split, MPI_INT, j, 2, MPI_COMM_WORLD, &status);
				for (int k = 0; k < number_after_split;k++) {
					v_nD.push_back(v_nD_temp_recv_after[k]);
				}

				//codes
				MPI_Recv(codes_Recv_after[0], (number_after_split*i_collapsing), MPI_INT, j, 3, MPI_COMM_WORLD, &status);
				for (int k10 = 0; k10 < number_after_split; k10++) {
					for (int k11=0; k11 < i_collapsing; k11++) {
						codes[counter6][k11] = codes_Recv_after[k10][k11];
					}
					counter6++;
				}
			}
		}//end of slaves

		 //Deallocation
		Del_iMatrix(codes_Recv_before, number_before_split, i_collapsing);
		Del_iMatrix(codes_Recv_after, number_after_split, i_collapsing);

	}//end of node 0

	 //------------------------------------
	 //Broadcast to all slaves
	 //-------------------------------------

	//List_nU
	if (mynode != 0) List_nU_size_total.resize(nrow_mox);
	MPI_Bcast(&List_nU_size_total[0], nrow_mox, MPI_INT, 0, MPI_COMM_WORLD);
	
	//v_nD
	if (mynode != 0) v_nD.resize(nrow_mox);
	MPI_Bcast(&v_nD[0], nrow_mox, MPI_INT, 0, MPI_COMM_WORLD);

	//codes
	MPI_Bcast(codes[0], nrow*i_collapsing, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);//In case all processors have List_nU_size_total


    //-----------------------------------------------------------
	//Send cells of List_nU_temp to master processor from slaves
    //-----------------------------------------------------------

	if (mynode != 0) {
		std::vector<double> List_nU_temp_cell;
		List_nU_temp.unlist(List_nU_temp_cell); //unlist all cells of List_nU_temp
		MPI_Send(&List_nU_temp_cell[0], List_nU_temp_cell.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}

	//-----------------------------------------
	//Aggregate cells of List_nU_temp on master
	//-----------------------------------------

	std::vector<double> List_nU_cell_total; // all cells of List_nU 
	if (mynode == 0) {

		int i_size_cell = 0; //number of cells in List_nU_temp on different processors
		int startindex = 0;
		int endindex = 0;

		for (int j = 1; j < totalnodes; j = j + 1) {
			if (j >= 1 && j <= splitnode) {

				i_size_cell = 0;
				startindex = 0;
				endindex = 0;
				startindex = (j - 1) * number_before_split;
				endindex = j *number_before_split;

				//cout << "startindex is " << startindex << " and endindex is " << endindex << " at node " << j << endl;

				for (int t = startindex; t < endindex; t++) {
					i_size_cell = i_size_cell + List_nU_size_total[t];
				}

				std::vector<double> List_nU_cell_before(i_size_cell);
				MPI_Recv(&List_nU_cell_before[0], i_size_cell, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);

				for (int k = 0; k < i_size_cell; k++) {
					List_nU_cell_total.push_back(List_nU_cell_before[k]);
				}
			}

			//-----------------------------

			if (j > splitnode) {

				i_size_cell = 0;
				startindex = 0;
				endindex = 0;
				startindex = splitnode * number_before_split + (j - splitnode - 1)*number_after_split;
				endindex = splitnode * number_before_split + (j - splitnode)*number_after_split;

				//cout << "startindex is " << startindex << " and endindex is " << endindex << " at node " << j << endl;

				for (int t = startindex; t < endindex; t++) {
					i_size_cell = i_size_cell + List_nU_size_total[t];
				}

				std::vector<double> List_nU_cell_after(i_size_cell);
				MPI_Recv(&List_nU_cell_after[0], i_size_cell, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);

				for (int k = 0; k < i_size_cell; k++) {
					List_nU_cell_total.push_back(List_nU_cell_after[k]);
				}
			}
		}//end of mynode
	}//end of node 0

	int i_size_cell_total = 0;

	for (int t = 0; t < nrow_mox; t++) {
		i_size_cell_total = i_size_cell_total + List_nU_size_total[t];
	}

	if (mynode == 0) {

		if (List_nU_cell_total.size() != i_size_cell_total) {
			cout << "ERROR! List_nU_cell_total is incorrect in nDAU function! " << endl;
		}

		//cout<<"max_size is "<< List_nU_cell_total.max_size() <<" and i_size_cell_total is "<< i_size_cell_total <<endl;
		if (i_size_cell_total > List_nU_cell_total.max_size()) {
			cout<<"CAUTION! There is potential memory issue for List_nu in nDAU function since exceeding max length of vector!"<<endl;
		}

	}

	//Broadcast all cells of List_nU to all processors
	if (mynode != 0) List_nU_cell_total.resize(i_size_cell_total);
	MPI_Bcast(&List_nU_cell_total[0], i_size_cell_total, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	//----------------------------------
	//Assemble List_nU on all processors
	//----------------------------------

	std::vector<double> List_nU_buffer;
	int counter8 = 0;
	int i_size_temp2 = 0;

	for (int k = 0; k < nrow_mox; k++) {
		List_nU_buffer.clear();
		i_size_temp2 = List_nU_size_total[k];
		for (int j = 0; j < i_size_temp2; j++) {
			List_nU_buffer.push_back(List_nU_cell_total[counter8]);
			counter8++;
		}

		List_nU.put_block_yicheng(k, List_nU_buffer.size(), List_nU_buffer);
	}



	//if (mynode == 1) {
		//for (int j = 0; j < nrow_mox; j++) {
		//	cout << "v_nD[" << j << "]: " << v_nD[j] << " at node " << mynode << endl;
		//}

		//cout << "List_nU at node " << mynode << endl;
		//List_nU.print_List_FHDI_yicheng();
		//for (int i = 0;i < List_nU.size_row();i++) {
		//	cout << "List_nU[" << i << "]: " << endl;
		//	List_nU.print_one_List_FHDI(i);
		//}
		//cout<<"nDAU codes at node "<<mynode<<endl;
		//for (int i = 0; i < nrow;i++) {
		//	for (int j = 0;j < i_collapsing; j++) {
		//		cout << setw(20) << codes[i][j];
		//	}
		//	cout << endl;
		//}

	//}

	//----------------------------

	//make a table-like information of nU

	//example:

	// List_nU

	// 0:  3,5,7

	// 1:  1,10

	// 2:  23,1,3, 5

	// then,

	// d_v_nU_unlist_temp: 3,5,7,1,10, 23,1,3,5,...

	// v_table_item_List_nU ; v_table_count_List_nU

	// 1                      2 

	// 3                      2 

	// 5                      2 

	// 7                      1 

	// 10                     1 

	// 23                     1  	

	//finally, 

	//tnU

	// 2, 2, 2, 1, 1

	//----------------------------

	std::vector<double> v_nU_unlist;

	List_nU.unlist(v_nU_unlist); //get the list of all entities of List_nU

	int i_size_v_nU_unlist = (int)v_nU_unlist.size();

	//----

	//Error check

	//-----

	if ((i_size_v_nU_unlist <= 0) && (i_cellmake == 2))

	{

		//Rprint("No possible donors with current k. Retry with reduced k \n");

		cout << "Causion!!! No possible donors with current k in nDAU_cpp. Retry with reduced k" << endl;
		//return 0;

	}


	//if (mynode == 0) {
	//	//cout << "YYC Running time of Cell_make at node " << mynode << " = " << MPI_Wtime() - cell_make_begin << endl;
	//	printf("Yicheng Running time of KNN ultra Bigp = %f seconds\n", MPI_Wtime() - KNN_begin);
	//}


	//Deallocation
	if (mynode != 0) {
		Del_dMatrix(mox_temp, L_temp, ncol);
		Del_iMatrix(codes_Send, L_temp, i_collapsing);
	}


	return 1;

}



//} //end of namespace

