#include <vector>
//#include "order_FHDI_binary.cc"
void Cal_W_Extension_Bigp_cpp(double** mox, const int nrow_mox,
	double** uox, const int nrow_uox, const int i_collapsing,
	const int ncol, int* id, int** codes,
	std::vector<std::string> v_table_tmvec_row1,
	std::vector<int> v_table_tmvec_row2,
	std::vector<double> jp_prob,
	double** d_mx, const int i_size_ml,
	double* w, std::string cn[], const int nrow,
	std::vector<double> &v_rst_final)
	//Description=========================================
	// update weight and joint probability
	//
	// Algorithm:  All possible donors will be used to fill in the missing cell 
	//             but, if there is no matched donors in uox, this algorithm may fail
	//             as of Oct 2016
	//
	// original R code: Dr. Im, J. and Dr. Kim, J. 
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: Oct 26, 2016
	//----------------------------------------------------
	//IN    : double mox(nrow_mox, ncol)= sorted unique patterns of missing  cells. up to i_count_mox rows are meaningful                           
	//IN    : double uox(nrow_uox, ncol)= sorted unique patterns of observed cells. up to i_count_uox rows are meaningful 
	//IN    : int id(nrow) = index of row. Default is ACTUAL row number
	//IN	: vector<string> v_table_tmvec_row1  = name of table of condensed missing patterns
	//IN	: vector<int> v_table_tmvec_row2  = counts of table of condensed missing patterns
	//IN   	: vector<double> jp_prob 	= weighted joint probability of all condensed observed DONORS
	//IN    : double d_mx(i_size_ml, ncol) = copy of all the missing cells 
	//IN    : double w[ml] 				= weights corresponding to missing rows
	//IN  	: string cn(nrow)			= vector of string to represent each row of z    
	//IN    : int i_option_collapsing = choice of big-p algorithm 
	//                               0= no big-p algorithms
	//                              !0= perform big-p algorithms
	//IN   :  int codes(nrow, i_option_collapsing); // storage to record most correlated variables of mox
	//OUT   : std::vector<double> &v_rst_final  = new weights 
	//====================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	double cal_W_check1 = MPI_Wtime();
	const int nr1 = nrow_mox;
	const int nr2 = nrow_uox;
	const int nx = i_size_ml; //rows of mx

							  //-------------------
							  //sum of joint probability
							  //-------------------
	double sum_jp = 0.0;
	const int i_size_jp = (int)jp_prob.size();
	for (int i = 0; i<i_size_jp; i++) sum_jp += jp_prob[i];

	//-------------------
	//initialize rst, the matrix for storage for augmented observations
	//-------------------
	rbind_FHDI rst(2); //number of columns 

					   //--------------------
					   //Main Loop for all missing rows
					   //--------------------
	int* i_temp_x = new int[ncol];
	int i_sum_x = 0;
	std::vector<int> v_cn_z_i;
	//int* zid = NULL;
	int i_size_zid = 0;
	int i_loc = 0;
	int* i_srst = new int[nx];		//for Condition 1&2
	int* i_srst1 = new int[nr2]; 	//for Condition 2
	std::vector<int> loc_srst_ncol;
	std::vector<int> loc_srst_ncol1;
	std::vector<int> v_mxl; // hold the most correlated variables of mox[i]

	//double* w_srst_ncol = NULL;
	//double* jp_zi = NULL;

	//----------------
	//LOOP for all missing rows
	//----------------
	for (int i = 0; i<nr1; i++)
	{
		//---------------------
		// generate sum of rows that indicate the matched rows of mx and mox
		//---------------------
		//indicator matrix that matches the donors
		//srst: row-wise sum of the indicator matrix 
		//-------
		loc_srst_ncol.clear(); //re-initialize
		Fill_iVector(i_srst, nx, 0); //re-initialize 

		for (int j = 0; j<nx; j++)  //Loop for i_size_ml, all the missing rows
		{
			int i_sum_crst = 0;
			for (int k = 0; k<ncol; k++)
			{
				//Note: in below check, mox is fixed at ith row 
				if (fabs(mox[i][k] - d_mx[j][k])<1e-3) //part of missing cell = obserbed cell 
				{
					i_sum_crst++; // increment if a cell of missing row = obs. cell 
				}
			}
			//---
			//store how many cells of the current missing row match those of all missing rows
			//---
			i_srst[j] = i_sum_crst;

			//---
			//store numbers of missing rows that exactly match the current missing row
			//i.e., target rows to be imputed later 
			//---
			if (i_sum_crst == ncol) loc_srst_ncol.push_back(j + 1); //Actual location 				
		}
		//-----
		//how many missing rows have the same missing pattern as the current missing row
		//-----
		const int i_size_loc_srst_ncol = (int)loc_srst_ncol.size();

		//---------------------------------
		//get current row of missing cell 
		//---------------------------------
		for (int j = 0; j<ncol; j++) i_temp_x[j] = mox[i][j];
		i_sum_x = sum_FHDI(i_temp_x, ncol); //how many non-zeros in current missing row

		std::string s_temp = v_table_tmvec_row1[i]; //string name of ith missing row

		v_cn_z_i.clear(); //re-initialize 
		which(cn, nrow, s_temp, v_cn_z_i); //Note: Actual location is returned
		int i_size_v_cn_z_i = (int)v_cn_z_i.size(); //number of locations in cn having s_temp

													//----------------------
													//----------------------
													//Condition 1: this row's cells are all missing
													//----------------------
													//----------------------
		if (i_sum_x == 0)
		{
			//-----------------
			//make "zid" which means 
			//the row location of current missing row repeated by number of observed rows
			//-----------------
			//zid = NULL; //re-initialize; 

						//-----
						// "i_size_v_cn_z_i" means all the missing rows that have the same pattern as the current missing row
						// so, below "i_size_zid" means that all the observed rows (nr2) will fill the missing rows 
						//-----
			i_size_zid = i_size_v_cn_z_i*nr2;
			int* zid = new int[i_size_zid];

			for (int j = 0; j<i_size_v_cn_z_i; j++)
			{
				for (int k = 0; k<nr2; k++) //nr2 times repeated copy with the id number 
				{
					//NOTE: zid contains ACTUAL id number 
					//Meaning id's of the missing rows that have the identical pattern as the current missing rows 
					zid[j*nr2 + k] = id[v_cn_z_i[j] - 1]; //-1 for actual location
				}
			}

			//-------------------
			//get ready w[] at srst = ncol
			//-------------------
			//w_srst_ncol = NULL; //re-initialize
			double* w_srst_ncol = new double[i_size_loc_srst_ncol*nr2];
			for (int j = 0; j<i_size_loc_srst_ncol; j++)  //repeat each entity by nr2 times
			{
				for (int k = 0; k<nr2; k++)
				{
					w_srst_ncol[j*nr2 + k] = w[loc_srst_ncol[j] - 1]; //-1 for actual location
				}
			}

			//-------------------
			//get ready second column
			//-------------------
			const int z_i_now = v_table_tmvec_row2[i];
			//jp_zi = NULL; //re-initialize 
			double* jp_zi = new double[i_size_jp * z_i_now];
			for (int j = 0; j<z_i_now; j++)  //repeat entire jp..[] by z_i_now times 
			{
				for (int k = 0; k<i_size_jp; k++)
					jp_zi[j*i_size_jp + k] = jp_prob[k] / sum_jp;
			}

			//-------------------
			//make a matrix that consists of zid & repeated weights for all missing rows
			//-------------------
			double** rst_temp = New_dMatrix(i_size_zid, 2);

			for (int j = 0; j<i_size_v_cn_z_i; j++)
			{
				for (int k = 0; k<nr2; k++) //repeated copy of the id number 
				{
					i_loc = j*nr2 + k; //serial number of the entire rows of the matrix

									   //first column is zid[]
					rst_temp[i_loc][0] = zid[i_loc];

					//second column joint prob * weight 
					rst_temp[i_loc][1] = jp_zi[i_loc] * w_srst_ncol[i_loc];
				}
			}
			//---
			//Append the entire matrix to rst
			//---
			rst.bind_blocks(i_size_zid, 2, rst_temp);

			//testout
			/*
			RPrint(" == in Cal_W Condition 1 ==== i: "); RPrint(i);
			RPrint("zid:"); RPrint(zid, i_size_zid);
			RPrint("rst:");
			rst.print_rbind_FHDI();
			*/

			//---------
			//local deallocation
			//---------
			Del_dMatrix(rst_temp, i_size_zid, 2);
			delete[] zid;
			delete[] w_srst_ncol;
			delete[] jp_zi;
		}

		//----------------------
		//----------------------
		//Condition 2: some cells of current row are not missing
		//----------------------
		//----------------------
		int nl = 0;
		if (i_sum_x > 0)
		{
			//------
			//number of observed cells on this row
			//------
			nl = 0;
			for (int j = 0; j<ncol; j++)
			{
				if (mox[i][j]>0) nl++; //number of the observed 
			}

			if (nl > i_collapsing) {
				nl = i_collapsing;
			}

			//-------
			//indicator matrix that matches the donors
			//srst: row-wise sum of the indicator matrix 
			//-------
			loc_srst_ncol1.clear(); //re-initialize //Note: this is different from loc_srst_ncol
			Fill_iVector(i_srst1, nr2, 0); //re-initialize 
			v_mxl.clear();

			//inherents the most correlated variables of mox[i]
			for (int k = 0; k < i_collapsing; k++) {
				if (codes[i][k] != 0) {
					v_mxl.push_back(codes[i][k]);
				}
			}

			int v_mxl_size = v_mxl.size();

			for (int j = 0; j<nr2; j++)
			{
				int i_sum_crst = 0;
				for (int k = 0; k<v_mxl_size; k++)
				{
					//Note: in below check, mox is fixed at ith row 
					if (fabs(mox[i][v_mxl[k] - 1] - uox[j][v_mxl[k] - 1])<1e-3) //part of missing cell = obserbed cell 
					{
						i_sum_crst++; // increment if a cell of the current missing row = obs. cell 
					}
				}
				//---
				//store how many cells of the current missing row match those of the observed row
				//---
				i_srst1[j] = i_sum_crst;

				//---
				//store row number of the observed that matches the current missing row
				//---
				if (i_sum_crst == nl) loc_srst_ncol1.push_back(j + 1); //Actual location 				
			}
			//testout
			/*
			RPrint(" == in Cal_W Condition2. ====== i: "); RPrint(i);
			RPrint("nl: "); RPrint(nl);
			RPrint("srst: "); RPrint(i_srst, nr2);
			RPrint("loc_srst_ncol: "); RPrint(loc_srst_ncol);
			*/

			//-----
			//total number of the observed rows that matches the current missing row
			//-----
			const int i_size_loc_srst_ncol1 = (int)loc_srst_ncol1.size();
			if (i_size_loc_srst_ncol1 == 0) //error case
			{
				cout << "Error! there is no matched cell!" << endl; return;
			}

			if (i_size_loc_srst_ncol1 > 0)
			{
				//-----------------
				//make "zid" which means 
				//the row location of current missing row repeated by number of observed rows
				//-----------------
				//zid = NULL; //re-initialize; 
				i_size_zid = i_size_v_cn_z_i * i_size_loc_srst_ncol1; //Note: .._ncol1 is used NOT .._ncol
				int* zid = new int[i_size_zid];

				for (int j = 0; j<i_size_v_cn_z_i; j++)
				{
					for (int k = 0; k<i_size_loc_srst_ncol1; k++) //repeated copy of the id number 
					{
						//NOTE: zid contains ACTUAL id number 
						zid[j*i_size_loc_srst_ncol1 + k] = id[v_cn_z_i[j] - 1]; //-1 for actual location
					}
				}

				//-------------------
				//get ready w[] at srst = ncol
				//-------------------
				//w_srst_ncol = NULL; //re-initialize
				double* w_srst_ncol = new double[i_size_loc_srst_ncol * i_size_loc_srst_ncol1];
				for (int j = 0; j<i_size_loc_srst_ncol; j++)  //repeat each entity 
				{
					for (int k = 0; k<i_size_loc_srst_ncol1; k++)
					{
						//------
						//Note: the weights are pulled out from loc_srst_loc NOT .._loc1 
						//      below "loc_srst_ncol" means the target rows to be imputed later
						//------
						w_srst_ncol[j*i_size_loc_srst_ncol1 + k] = w[loc_srst_ncol[j] - 1]; //-1 for actual location
					}
				}

				//-------------------
				//get ready second column
				//below "loc_srst_ncol1" contains the obs. row numbers that will serve as donor
				//-------------------
				double sum_jp_loc = 0.0; //sum of jp only at location where srst = ncol
				for (int j = 0; j<i_size_loc_srst_ncol1; j++)
					sum_jp_loc += jp_prob[loc_srst_ncol1[j] - 1];

				const int z_i_now = v_table_tmvec_row2[i];
				double* jp_zi = new double[i_size_loc_srst_ncol1 * z_i_now];
				for (int j = 0; j<z_i_now; j++)  //repeat entire jp..[] by z_i_now times 
				{
					for (int k = 0; k<i_size_loc_srst_ncol1; k++)
						jp_zi[j*i_size_loc_srst_ncol1 + k]
						= jp_prob[loc_srst_ncol1[k] - 1] / sum_jp_loc;
				}


				//-------------------
				//make a matrix that consists of zid & repeated uox for all missing rows
				//-------------------
				double** rst_temp2 = New_dMatrix(i_size_zid, 2);

				for (int j = 0; j<i_size_v_cn_z_i; j++)
				{
					for (int k = 0; k<i_size_loc_srst_ncol1; k++) //repeated copy of the id number 
					{
						i_loc = j*i_size_loc_srst_ncol1 + k; //serial number of the entire rows of the matrix

															 //first column is zid[]
						rst_temp2[i_loc][0] = zid[i_loc];

						//second column. joint prob * weight 
						rst_temp2[i_loc][1] = jp_zi[i_loc] * w_srst_ncol[i_loc]; //-1 for actual location
					}
				}
				//---
				//Append the entire matrix to rst
				//---
				rst.bind_blocks(i_size_zid, 2, rst_temp2);

				//testout
				/*
				RPrint("zid:"); RPrint(zid, i_size_zid);
				RPrint("rst:");
				rst.print_rbind_FHDI();
				*/

				//---------
				//local deallocation
				//---------
				Del_dMatrix(rst_temp2, i_size_zid, ncol + 1);
				delete[] zid;
				delete[] w_srst_ncol;
				delete[] jp_zi;
			}
		}

	} //end of LOOP for all missing rows
	//cout << "cal_W_check1 at node " << mynode << " = " << MPI_Wtime() - cal_W_check1 << endl;
	//----------------
	//re-order rst in terms of id (the first column)
	//----------------
	double cal_W_check2 = MPI_Wtime();
	const int n_row_rst = rst.size_row();
	int* i_rst_id = new int[n_row_rst];
	for (int i = 0; i<n_row_rst; i++) i_rst_id[i] = (int)rst(i, 0);
	order_FHDI_binary(i_rst_id, n_row_rst); //returned with the order of rows in ascending magnitude
									 //testout
									 //RPrint("n_row_rst :"); RPrint(n_row_rst);
									 //RPrint("i_rst_id :"); RPrint(i_rst_id, n_row_rst);


									 //--------------------
									 //remove the first column with id
									 //store the rst into the final storage
									 //--------------------
	double* d_row_rst = new double[2];
	double  d_row_rst_short = 0.0;
	for (int i = 0; i<n_row_rst; i++)
	{
		rst.get_block(i_rst_id[i] - 1, d_row_rst); //get a row// -1 for actual loc
		d_row_rst_short = d_row_rst[1]; //without id  
		v_rst_final.push_back(d_row_rst_short);	//append a new row to the final storage 
	}

	//testout
	//RPrint("End of Cal_W =========="); 
	//RPrint("v_rst_final:"); RPrint(v_rst_final); 
	//cout << "cal_W_check2 at node " << mynode << " = " << MPI_Wtime() - cal_W_check2 << endl;

	//-------
	//local deallocation
	//-------
	delete[] i_temp_x;
	//delete[] zid;
	delete[] i_srst;
	delete[] i_srst1;
	//delete[] w_srst_ncol;
	//delete[] jp_zi;
	delete[] i_rst_id;
	delete[] d_row_rst;

	return;
}