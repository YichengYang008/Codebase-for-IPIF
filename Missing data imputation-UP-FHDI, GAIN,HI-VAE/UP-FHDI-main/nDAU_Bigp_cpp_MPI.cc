//Fn===========================================================================

//nDAU_cpp.cc-----------------------------------------------------------------------------

//Fn===========================================================================

//namespace FHDI{
#include "correlated_variable_intersection.cc"
#include "correlated_variable_union.cc"
#include "correlated_variable_global.cc"

bool nDAU_Bigp_cpp_MPI(double** uox, double** mox, const int nrow_uox, const int nrow_mox, const int nrow, const int ncol, const int i_collapsing, const int i_SIS_type,

	std::string cn[], int* ol, const int nrow_ol, const int top, int i_cellmake,

	std::vector<int> &v_nD, List_FHDI &List_nU, int* tnU, int** codes, int** correlation_ranking_top, double** ol_matrix,

	bool b_DEBUG, ofstream& TestOut, int i_loop)

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

	//IN    : double uox(nrow_uox, ncol)= sorted unique patterns of observed cells. up to i_count_uox rows are meaningful 

	//IN    : double mox(nrow_mox, ncol)= sorted unique patterns of missing  cells. up to i_count_mox rows are meaningful                           

	//IN 	: string cn(nrow)		= vector of string to represent each row of z          

	//IN	: int ol(nrow_ol)		= actual location of rows containing ONLY observed cells    

	//IN    : int i_option_collapsing = choice of big-p algorithm 

	//                               0= no big-p algorithms

	//                              !0= perform big-p algorithms

	//INOUT : std::vector v_nD		= total number of donnors of each missing pattern

	//OUT   : List_FHDI List_nU     = list of observed cells to serve as donors 

	//OUT   : int tnU[nrow_uox]		= table format of the total numbers of donors for each missing rows

	//OUT   : int codes(nrow, i_option_collapsing); // storage to record most correlated variables of mox

	//====================================================

{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;


	//initialize 	

	//double* snr2 = new double[nrow_mox];

	//for (int i = 0; i < nrow_mox; i++) { snr2[i] = i + 1; } //actual row number of mox



	//------------

	//make a table of unique strings in cn

	//------------

	std::vector<std::string> v_table_item_cn;// uox

	std::vector<int>	     v_table_count_cn;// i_count_uox



	//observed patterns only

	//std::string s_cn_ol[nrow_ol]; 

	std::string *s_cn_ol = new std::string[nrow_ol];

	for (int i = 0; i < nrow_ol; i++) { s_cn_ol[i] = cn[ol[i] - 1]; } //"-1" for actual location number of row

    std::sort(s_cn_ol, s_cn_ol + nrow_ol); //sort the observed patterns 



	table_cpp(s_cn_ol, nrow_ol,

		v_table_item_cn, v_table_count_cn);


	//std::vector<std::string>().swap(v_table_item_cn);
	delete[] s_cn_ol;
	//for (int i = 0;i < v_table_item_cn.size();i++) {
	//	cout << "v_table_item_cn[" << i << "]: " << v_table_item_cn[i] << endl;
	//}

	//for (int i = 0;i < v_table_count_cn.size();i++) {
	//	cout << "v_table_count_cn[" << i << "]: " << v_table_count_cn[i] << endl;
	//}

	//-----------

	//null string of width of ncol

	//-----------

	std::string s_zval; for (int i = 0; i < ncol; i++) s_zval.append("0");



	//-------------------

	//Loop for mox (missing cells patterns) rows

	//-------------------

	double* d_temp = new double[ncol]; // store mox[i]

	int* ia_temp = new int[ncol]; //temporary integer array to store mox[i]



	std::string s_cn0; //temporary string for a row of mox

	//-------------------------------------

	//example

	// mox[i,] 		13 0 1

	// mxl 			1    3		//location of observed cells in mox

	// nxl			2         	//number of observed cells in mox

	// rcn0			d1			//condensed string of the observed cells in mox (from 35 letters, 1-9 and a-z)

	// cand 		"11" "29" "d1" "b3" "d1" ... 	//condensed strings of uox corresponding to the observed cell columns

	// oloc			3    5     	//location of cand that has the same as rcn0

	//-------------------------------------

	//initialize 
	std::vector<int> v_nD_temp;
	std::vector<int> List_nU_value;
	std::vector<int> List_nU_size;
	int** codes_Send = New_iMatrix(nrow, i_collapsing);// The codes at each node 

	int number_before_split = 0;
	int number_after_split = 0;
	int startpoint = 0;
	int endpoint = 0;
	int splitnode = 0;
	//cout<<"nrow_mox at mynode "<<mynode<<" is "<< nrow_mox <<endl;
	//cout << "nrow_uox at mynode " << mynode << " is " << nrow_uox << endl;
	if (nrow_mox < ( totalnodes-1 )) {
		TestOut <<"Error! The number of mox is smaller than the total number of avaiable nodes in nDAU_Bigp. Please reduce the total number of nodes"<< endl;
	}

	if (nrow_mox % (totalnodes - 1) != 0) { splitnode = nrow_mox % (totalnodes - 1); }
	if (nrow_mox % (totalnodes - 1) == 0) { splitnode = 1; }
	//if (mynode == 0) cout << "Split: " << splitnode << endl;
	number_after_split = floor(1.0*nrow_mox / (1.0*totalnodes - 1));
	number_before_split = 1.0*(nrow_mox - floor(1.0*nrow_mox / (1.0*totalnodes - 1)) *(totalnodes - splitnode - 1)) / splitnode;

	if (mynode >= 1 && mynode <= splitnode) {
		startpoint = (mynode - 1) * number_before_split;
		endpoint = mynode *number_before_split;
	}
	if (mynode > splitnode) {
		startpoint = splitnode * number_before_split + (mynode - splitnode - 1)*number_after_split;
		endpoint = splitnode * number_before_split + (mynode - splitnode)*number_after_split;
	}
	if ((number_before_split*splitnode + number_after_split* (totalnodes - splitnode - 1)) != nrow_mox) {
		TestOut << "Error !!! Work Assignment Error in nDAU_BigP!!!!" << endl;
		return 0;
	}
	//if (mynode == 0) {
	//	cout << "Mynode: " << mynode << ", number_before_split: " << number_before_split << ", number_after_split:" << number_after_split << endl;
	//}
	//cout << "Mynode " << mynode << ", " << startpoint << "<= x < " << endpoint << endl;

	std::vector<int> v_mxl; // 1. temporary vector for the locaiton of non zeros in mox[i]
	                        // 2. temporary vector for the locaitons of selected non zero variable for mox[i]
	std::vector<int> v_oloc;

	std::vector<std::string> v_cand; //vector of found string with condensed non-zero observed data, string format of corresponding mox regarding select variables of mox[i]

	//if ((mynode == 0)||(mynode == 2)) cout<<"nDAU_bigp_1 at "<< i_loop <<" at node "<<mynode<<endl;
	for (int i = startpoint; i < endpoint; i++)

	{

		for (int j = 0; j < ncol; j++) d_temp[j] = mox[i][j];  //ith row of mox


		//cout<<"Debug_nDAu 1 at node "<<mynode<<endl;
		Trans1(d_temp, ncol, s_cn0); //condense a row to string 





		//---------

		//when current missing row is null string, i.e. "   "

		//---------

		if (s_cn0.compare(s_zval) == 0) //0: equal string

		{
			//cout << "Debug_nDAu 2 at node " << mynode << endl;
			//---------------------

			//number of donors; this case all observed cells are possible donors

			//---------------------

			int i_nD_sum = 0;

			for (unsigned k = 0; k < v_table_count_cn.size(); k++) i_nD_sum += v_table_count_cn[k];

			v_nD_temp.push_back(i_nD_sum); //store the number of possible donors into the integer vector to return



			//------

			//store a row of nU into the List storage

			//for this null string row, all observed rows become possible donors

			//------	

			//double* d_nU_temp = new double[nrow_uox]; 
			List_nU_size.push_back(nrow_uox);
			//cout<<"List_nU_size["<<i<<"]:"<< nrow_uox <<endl;

			for (int k = 0; k < nrow_uox; k++) {
				//d_nU_temp[k] = k + 1; 
				List_nU_value.push_back(k + 1);
			}

			//List_nU.put_block(i, nrow_uox, d_nU_temp);

			//cout << "Debug_nDAu 3 at node " << mynode << endl;
			//delete[] d_nU_temp; 

		}



		//----------

		//for general cases for missing units, other than null string

		//----------

		if (s_cn0.compare(s_zval) != 0) //0: equal string

		{
			double nDAU1 = MPI_Wtime();

			//-----

			//find non zero cells of current missing row mox

			//-----

			int oc = 0; // Get number of observed values in current mox

			for (int k = 0; k < ncol; k++) {

				ia_temp[k] = (int)mox[i][k];

				if (ia_temp[k] > 0) {

					oc++;
				}

			}


			//Note: below will contain ACTUAL location of cells with non-zero observed data

			//std::vector<int> v_mxl; //temporary vector for the locaiton of non zeros in mox
			v_mxl.clear();

			if (oc < (i_collapsing + 1)) { // oc <= i_option_collapsing

				whichINV(ia_temp, ncol, 0, v_mxl); //get the location of Non-zero in mox 

			}

			if (oc > i_collapsing) {

				if (i_SIS_type == 1) {

					//correlated_variable_intersection2(ncol, i_collapsing, i, ia_temp, correlation_yicheng, correlation_ranking, v_mxl, TestOut);
					correlated_variable_intersection2(ncol, i_collapsing, top, i, nrow_ol, ia_temp, ol_matrix, correlation_ranking_top, v_mxl, TestOut);

				}

				if (i_SIS_type == 2) {

					//correlated_variable_union(ncol, i_collapsing, i, ia_temp, correlation_yicheng, correlation_ranking, v_mxl, TestOut);
					correlated_variable_union2(ncol, i_collapsing, top, i, nrow_ol, ia_temp, ol_matrix, correlation_ranking_top, v_mxl, TestOut);
				}

				if (i_SIS_type == 3) {

					//correlated_variable_global(ncol, i_collapsing, ia_temp, correlation_yicheng, v_mxl, TestOut);
					correlated_variable_global2(ncol, i_collapsing, nrow_ol, ia_temp, ol_matrix, v_mxl, TestOut);

				}

			}


			for (int b3 = 0; b3 < v_mxl.size(); b3++) {

				codes_Send[i][b3] = v_mxl[b3];
				//cout << "v_mxl_Top[" << b3 << "]: " << v_mxl[b3];
				//cout<<"codes["<<i<<"]["<<b3<<"]: "<< v_mxl[b3] <<" at node "<<mynode<<endl;
			}
			//cout << endl;
			//-----

			//if (mynode == 1 && i_loop == 1) {
			//	cout << "YYC Running time of nDAU1 at iteration " << i_loop <<" at i = "<<i<< " = " << MPI_Wtime() - nDAU1 << endl;
			//}

			double nDAU2 = MPI_Wtime();

			//total number of non zeros in mox

			//-----

			const int nxl = v_mxl.size(); //1. total number of observed units at current row mox[i] (>0)
			                              //2. total number of selected non zero variables in mox[i]

			double* d_rcn0_temp = new double[nxl]; //temporary

			for (int k = 0; k < nxl; k++) d_rcn0_temp[k] = mox[i][v_mxl[k] - 1]; //"-1" is for actual location



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

			//std::vector<int> v_oloc;
			v_oloc.clear();

			if (nxl >= 1) //unlike R code, below algorithm suffices for all cases 

			{

				//----

				//get all non-zero cells from all observed rows

				//----

				double * d_t1 = new double[v_mxl.size()];

				//std::vector<std::string> v_cand; //vector of found string with condensed non-zero observed data, string format of corresponding mox regarding select variables of mox[i]
				v_cand.clear();
				//-----

				//search all observed rows

				//of which the same columns are non-zero as the current missing row  

				//-----

				for (int m = 0; m < nrow_uox; m++)

				{	//Note: "-1" in v_mxl is from the ACTUAL location info in it

					for (unsigned k = 0; k < v_mxl.size(); k++) d_t1[k] = uox[m][v_mxl[k] - 1];



					//------

					//condense the found rows with non-zero observed cell only  

					//------

					std::string s_cand_1;

					Trans1(d_t1, v_mxl.size(), s_cand_1);

					v_cand.push_back(s_cand_1); //add more string to the string vector 

				}
				//cout << "v_cand size : " << v_cand.size() << endl;


				//--------------

				//Find the rows of v_cand that match the current non-zero missing pattern

				//Note: below will contain ACTUAL locations of the found rows 

				//--------------

				which(v_cand, s_rcn0, v_oloc); //get the locations of observed cells containing s_rcn0 



				//------------

				//local deallocation

				//------------

				delete[] d_t1;

			}



			//----------------

			//Store oloc into LIST named nU 

			// ith row of nU corresponds to ith row of the List

			//----------------

			int i_oloc_temp = (int)v_oloc.size(); //size of current oloc
			//cout<<"i_oloc_temp["<<i<<"] at node "<<mynode<<" is "<< i_oloc_temp <<endl;
			
			// Only user for parallel nDAU, which add 0 size to List_nU cutoff points. 
			if (i_oloc_temp == 0) { 
				List_nU_size.push_back(0); 
				//cout << "List_nU_size[" << i << "]:" << i_oloc_temp << endl;
			}
			
			if (i_oloc_temp > 0) //only for meaningful oloc

			{

				//double* d_oloc_temp = new double[i_oloc_temp]; 
				List_nU_size.push_back(i_oloc_temp);
				//cout << "List_nU_size[" << i << "]:" << i_oloc_temp << endl;
				for (int k = 0; k < i_oloc_temp; k++) {
					//d_oloc_temp[k] = v_oloc[k]; 
					List_nU_value.push_back(v_oloc[k]);
				}

				//List_nU.put_block(i, i_oloc_temp, d_oloc_temp); //store ith missing row's donors list into the block  

				//cout << "Debug_nDAu 5 at node " << mynode << endl;
				//delete[] d_oloc_temp; 

			}

			//---------------------

			//number of donors; this case only the matched observed rows become possible donors

			//---------------------

			int i_temp_tocn_sum = 0;

			for (int k = 0; k < i_oloc_temp; k++) //accumulate all possible donors 

			{
				i_temp_tocn_sum += v_table_count_cn[v_oloc[k] - 1];
			} //-1 for actual loc

			v_nD_temp.push_back(i_temp_tocn_sum); //store into integer vector to return 

			//if (mynode == 1 && i_loop == 1) {
			//	cout << "YYC Running time of nDAU2 at iteration " << i_loop <<" at i = "<<i<< " = " << MPI_Wtime() - nDAU2 << endl;
			//}



		} //end of general missing case, other than null string row case

	} //end of the main loop for i of all missing patterns  

	  delete[] d_temp;

	  delete[] ia_temp;
	//cout<<"List_nU_size.size() at mynode "<<mynode<<" is "<< List_nU_size.size() <<endl;
	//if (mynode != 0) {
	//	TestOut<<"codes at node "<<mynode<<endl;
	//	for (int i = 0; i < nrow_mox; i++) {
	//		for (int j = 0; j < i_collapsing; j++) {
	//			TestOut << setw(20) << codes_Send[i][j];
	//		}
	//		TestOut << endl;
	//	}
	//}
	//cout << "nDAU main loops finished at node " << mynode << endl;
	//Send and Recv
	//int List_nU_length_send = List_nU_value.size();
	//int List_nU_length_recv_before = 0;
	//int List_nU_length_recv_after = 0;

	//if (mynode == splitnode) {
	//	MPI_Send(&List_nU_length_send, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	//}

	//if (mynode == (splitnode + 1)) {
	//	MPI_Send(&List_nU_length_send, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	//}
	//cout << "List_nU_length_send at node " << mynode <<" : "<< List_nU_value.size() <<endl;


	//if (mynode == 0) {
	//	MPI_Recv(&List_nU_length_recv_before, 1, MPI_INT, splitnode, 1, MPI_COMM_WORLD, &status);

	//	cout<<"First pass"<<endl;
	//	MPI_Recv(&List_nU_length_recv_after, 1, MPI_INT, (splitnode + 1), 1, MPI_COMM_WORLD, &status);

	//	cout<<"List_nU_length_recv_before: "<< List_nU_length_recv_before <<endl;
	//	cout << "List_nU_length_recv_after: " << List_nU_length_recv_after << endl;
	//}
	//MPI_Barrier(MPI_COMM_WORLD);

	//-----------------------
	//cout << "Debug_nDAu 7 at node " << mynode << endl;
	//if ((mynode == 0) || (mynode == 2)) cout << "nDAU_bigp_2 at " << i_loop << " at node "<<mynode<<endl;


	//List_nU

	std::vector<int> List_nU_vector;


    //v_nD
	std::vector<int> v_nD_temp_recv_before;
	std::vector<int> v_nD_temp_recv_after;
	if (mynode == 0) {
		v_nD_temp_recv_before.resize(number_before_split);
		v_nD_temp_recv_after.resize(number_after_split);
	}

	//List_nU size cutoff points
	std::vector<int> List_nU_size_total;
	std::vector<int> List_nU_size_before;
	std::vector<int> List_nU_size_after;

	if (mynode == 0) {
		List_nU_size_before.resize(number_before_split);
		List_nU_size_after.resize(number_after_split);
	}

	//codes
	int** codes_Recv = NULL;
	
	if (mynode == 0) { 
		codes_Recv = New_iMatrix(nrow, i_collapsing);
	}// The record of i_option_collapsing most correlated variables of each mox
		

	if (mynode != 0) {
		MPI_Send(&List_nU_size[0], List_nU_size.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&List_nU_value[0], List_nU_value.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&v_nD_temp[0], v_nD_temp.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(codes_Send[0], (nrow*i_collapsing), MPI_INT, 0, 1, MPI_COMM_WORLD);

	}

	//if (mynode == 2) {
	//	cout << "List_nU_size: " << List_nU_size.size() << ", and List_nU_value: " << List_nU_value.size() << ", and v_nD_temp_size: " << v_nD_temp.size()<<" at " << i_loop << endl;
	//}
	//if (mynode == 0) cout << "nDAU_bigp_3_1 at " << i_loop << endl;

	if (mynode == 0) {
		for (int j = 1; j < totalnodes; j = j + 1) {
			if (j >= 1 && j <= splitnode) {
				//List_nU_size
				MPI_Recv(&List_nU_size_before[0], number_before_split, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				for (int k = 0; k < number_before_split;k++) {
					List_nU_size_total.push_back(List_nU_size_before[k]);
				}
				//List_nU_vector
				int List_nU_length_recv_before = 0;
				for (int m = 0;m < number_before_split;m++) {
					List_nU_length_recv_before = List_nU_length_recv_before + List_nU_size_before[m];
				}
				//cout << "List_nU_length_recv_before at node " << j << " : " << List_nU_length_recv_before << endl;

				std::vector<int> List_nU_value_recv_before(List_nU_length_recv_before);
				MPI_Recv(&List_nU_value_recv_before[0], List_nU_length_recv_before, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				for (int p = 0; p < List_nU_length_recv_before;p++) {
					List_nU_vector.push_back(List_nU_value_recv_before[p]);
				}
				//std::vector<int>().swap(List_nU_value_recv_before);
				//v_nD
				MPI_Recv(&v_nD_temp_recv_before[0], number_before_split, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				for (int k = 0; k < number_before_split;k++) {
					v_nD.push_back(v_nD_temp_recv_before[k]);
				}

				//codes: Note continous index is distributed in all slaves, not from 0 for all slaves.
				MPI_Recv(codes_Recv[0], (nrow*i_collapsing), MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				for (int l = (j - 1)*number_before_split; l < j*number_before_split;l = l + 1) {
					for (int k1 = 0; k1 < ncol; k1 = k1 + 1) {
						codes[l][k1] = codes_Recv[l][k1];
						//cout << "mynode: " << j << ", l: " << l << ", k1: " << k1 << ", counter: " << counter << endl;
					}
				}
				//cout<<"J: "<< j <<endl;
			}//
			if (j > splitnode) {
				//List_nU_size
				MPI_Recv(&List_nU_size_after[0], number_after_split, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				for (int k = 0; k < number_after_split;k++) {
					List_nU_size_total.push_back(List_nU_size_after[k]);
				}
				//List_nU_vector
				int List_nU_length_recv_after = 0;
				for (int m = 0;m < number_after_split;m++) {
					List_nU_length_recv_after = List_nU_length_recv_after + List_nU_size_after[m];
				}
				//cout << "List_nU_length_recv_after at node " << j << " : " << List_nU_length_recv_after << endl;

				std::vector<int> List_nU_value_recv_after(List_nU_length_recv_after);
				MPI_Recv(&List_nU_value_recv_after[0], List_nU_length_recv_after, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				for (int p = 0; p < List_nU_length_recv_after;p++) {
					List_nU_vector.push_back(List_nU_value_recv_after[p]);
				}
				//std::vector<int>().swap(List_nU_value_recv_after);

				//v_nD
				MPI_Recv(&v_nD_temp_recv_after[0], number_after_split, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				for (int k = 0; k < number_after_split;k++) {
					v_nD.push_back(v_nD_temp_recv_after[k]);
				}

				//codes
				MPI_Recv(codes_Recv[0], (nrow*i_collapsing), MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				for (int l = (splitnode*number_before_split + (j-splitnode-1)*number_after_split); l < (splitnode*number_before_split + (j- splitnode)*number_after_split);l = l + 1) {
					for (int k1 = 0; k1 < ncol; k1 = k1 + 1) {
						codes[l][k1] = codes_Recv[l][k1];
						//cout << "mynode: " << j << ", l: " << l << ", k1: " << k1 << ", counter: " << counter << endl;
					}
				}
				//cout << "J: " << j << endl;

			}//

		}//end of j

		//for (int o = 0;o < List_nU_vector.size();o++) {
		//	cout<<"Master List_nU_vector["<<o<<"]: "<< List_nU_vector[o]<<endl;
		//}

		//for (int o = 0;o < List_nU_size_total.size();o++) {
		//	cout << "Master List_nU_size_total[" << o << "]: " << List_nU_size_total[o] << endl;
		//}

	}//end of mynode==0


	//cout << "Debug_nDAu 8 at node " << mynode << endl;
	if (mynode != 0) v_nD.resize(nrow_mox);
	MPI_Bcast(&v_nD[0], nrow_mox, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Bcast(codes[0], nrow*i_collapsing, MPI_INT, 0, MPI_COMM_WORLD);

	int List_nU_vector_size = 0;
	if (mynode == 0) { List_nU_vector_size = List_nU_vector.size(); }

	MPI_Bcast(&List_nU_vector_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//cout << "List_nU_vector_size at node " << mynode <<" is "<< List_nU_vector_size <<endl;

	if (mynode != 0) List_nU_vector.resize(List_nU_vector_size);
	MPI_Bcast(&List_nU_vector[0], List_nU_vector_size, MPI_INT, 0, MPI_COMM_WORLD);

	if (mynode != 0) List_nU_size_total.resize(nrow_mox);
	MPI_Bcast(&List_nU_size_total[0], nrow_mox, MPI_INT, 0, MPI_COMM_WORLD);
	
	int counter0 = 0;

	for (int v = 0;v < nrow_mox;v++) {
		double* List_nU_buffer = new double[List_nU_size_total[v]];
		for (int k = 0; k < List_nU_size_total[v];k++) {
			List_nU_buffer[k] = List_nU_vector[counter0];
			//cout << "mynode " << mynode << "counter0: " << counter0 << endl;
			counter0++;
		}
		List_nU.put_block(v, List_nU_size_total[v], List_nU_buffer);
		delete[] List_nU_buffer;
	}

	//if ((mynode == 0) || (mynode == 2)) cout << "nDAU_bigp_3 at " << i_loop <<" at node "<<mynode<<endl;
	//Debug Testout
	//cout << "The size of List_nU is " << List_nU.size_row() <<" at node "<<mynode<<endl;
	//cout << "The size of nrow_uox is " << nrow_uox << " at node " << mynode << endl;
	//cout << "The size of nrow_mox is " << nrow_mox << " at node " << mynode << endl;
	//cout << "The size of v_nD is " << v_nD.size() << " at node " << mynode << endl;

	//for (int i = 0;i < nrow_uox;i++) {
	//	cout << "tnU[" << i << "]:" << tnU[i] << endl;
	//}

	//for (int i = 0; i < v_nD.size();i++) {
	//	cout << "v_nD[" << i << "]:" << v_nD[i] << endl;
	//}

	//for (int i = 0;i < List_nU.size_row();i++) {
	//	cout << "List_nU[" << i << "]: " << endl;
	//	List_nU.print_one_List_FHDI(i);

	//}
	//RPrint("List_nU:");List_nU.print_List_FHDI();


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

	//int i_size_v_nU_unlist = (int)List_nU_vector.size();



	//----

	//Error check

	//-----


	if ((i_size_v_nU_unlist <= 0) && (i_cellmake == 2))

	{

		//Rprint("No possible donors with current k. Retry with reduced k \n");

		cout << "Causion!!! No possible donors with current k in nDAU_cpp. Retry with reduced k" << endl;
		//return 0;

	}

	if ((i_size_v_nU_unlist <= 0) && (i_cellmake == 1))

	{

		//Rprint("No possible donors with current k. Retry with reduced k \n");

		cout << "ERROR!!! No possible donors with current k in nDAU_cpp bigp. Retry with reduced k" << endl;

		//exit(0);

		return 0;

	}



	double* d_v_nU_unlist_temp = new double[i_size_v_nU_unlist];

	//for (int k = 0; k < i_size_v_nU_unlist; k++) d_v_nU_unlist_temp[k] = (double)List_nU_vector[k]; //a copy of all donors (row numbers) 

	for (int k = 0; k < i_size_v_nU_unlist; k++) d_v_nU_unlist_temp[k] = v_nU_unlist[k]; //a copy of all donors (row numbers) 

	std::vector<double> v_table_item_List_nU; //names of List_nU

	std::vector<int>	v_table_count_List_nU;//counts of List_nU



	table_cpp(d_v_nU_unlist_temp, i_size_v_nU_unlist,

		v_table_item_List_nU, v_table_count_List_nU);


	delete[] d_v_nU_unlist_temp;
	//testout

	if (b_DEBUG)

	{
		RPrint("table_cpp has been done ");
	}



	Fill_iVector(tnU, nrow_uox, 0);

	for (int i = 0; i < nrow_uox; i++)

	{

		//-----

		//error check

		//sometimes tnU size is less than v_table_item_List_nU 

		//-----

		if (i >= (int)v_table_item_List_nU.size()) { break; }



		//-----

		//search meaningful locations to be stored into tnU 

		//-----

		double d_t2 = v_table_item_List_nU[i]; //Note: Actual number is stored!

		//if(std::isnan(d_t2)==1) {break;} //Exit at the end of the table list. only for meaningful number 
		if (isnan_FHDI(d_t2) == 1) { break; } //Exit at the end of the table list. only for meaningful number 


		for (int j = 1; j < nrow_uox + 1; j++) //"+1" is needed for Actual # stored

		{

			if (fabs_FHDI(d_t2 - j) < 1e-15) { tnU[i] = v_table_count_List_nU[i]; break; }

		}

	}

	//if ((mynode == 0) || (mynode == 2)) cout << "nDAU_bigp_4 at " << i_loop << " at node "<<mynode<<endl;


	//Deallocation

	//-----------------

	//delete[] s_cn_ol;

	//delete[] snr2;

	//delete[] d_temp;

	//delete[] ia_temp;

	//delete[] d_v_nU_unlist_temp;

	//delete[] tnU;
	Del_iMatrix(codes_Send, nrow, i_collapsing);

	if (mynode == 0) {
		Del_iMatrix(codes_Recv, nrow, i_collapsing);
	}

	return 1;

}



//} //end of namespace

