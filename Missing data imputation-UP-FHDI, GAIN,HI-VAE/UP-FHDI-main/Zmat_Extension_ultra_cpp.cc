
bool Zmat_Extension_ultra_cpp(const int nrow, const int ncol, const int memory,
	List_FHDI &uox_infor, List_FHDI &mox_infor, int& uox_size, int& mox_size,
	MPI_File fh_datz, MPI_File fh_uox, MPI_File fh_mox,
	ofstream& TestOut)
	//Description=========================================
	// make the condensed expression of z
	//
	// Algorithm:  each row of z will be concatenated as a single string consisting of 35 characters
	// 
	// Note: as of Oct 2016, NA values (missing data) is marked by a long integer at the parent "r" code
	//
	// original R code: Dr. Im, J. and Dr. Kim, J. 
	// c++ code: 		Yicheng Yang and Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: August 4, 2021
	//----------------------------------------------------
	//IN    : nrow                = number of row of daty
	//IN    : ncol                = number of column of daty
	//OUT   : List_FHDI uox_infor = Actual index list of uox in z
	//OUT   : List_FHDI mox_infor = Actual index list of mox in z
	//OUT   : uox_size            = number of uox
	//OUT   : mox_size            = number of mox
	//OUT   : MPI_File fh_datz    = file string to read categorized matrix corresponding to original matrix x
	//OUT   : MPI_File fh_uox     = file string to write unique observed patterns (not sorted)
	//OUT   : MPI_File fh_mox     = file string to write unique missing patterns (not sorted)                      

	//IMPORTANT SUMMARY OF INDEX 
	//Matrix         starting offset
	//summary        1
	//remove_list    0
	//uox_info       1
	//mox_infor      1
	//ol             0
	//ml             0
	//====================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	double read_begin = MPI_Wtime();
	const int L = nrow; //size of d_rw 
	int numWorkPerProc = (int)floor(1.0*L / (1.0*totalnodes - 1));
	int numWorkLocalLast = L - numWorkPerProc * (totalnodes - 2);


	int startpoint = 0;
	int endpoint = 0;

	if (mynode != 0 && mynode != (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = mynode*numWorkPerProc;
	}

	if (mynode == (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = (mynode - 1)*numWorkPerProc + numWorkLocalLast;
	}

	int L_temp = 0;
	if ( (mynode!=0) && (mynode != (totalnodes - 1)) ) L_temp = numWorkPerProc;
	if (mynode == (totalnodes - 1)) L_temp = numWorkLocalLast;

	double** z_temp = New_dMatrix(L_temp, ncol);

	//cout<<"zmat L_temp is "<< L_temp <<" where startpoint is "<< startpoint <<" and endpoint is "<< endpoint <<" at node "<< mynode <<endl;

	//========================================================
	//Read distributed z matrix on all slave processors
	//=========================================================

	MPI_In_datz(L_temp, ncol, numWorkPerProc, fh_datz, z_temp);

	//if (mynode == 1 || mynode==0) {
	//	cout << "Yicheng Running time of reading in Zmat = "<< MPI_Wtime() - read_begin <<" at node "<<mynode<< endl;
	//	printf("Yicheng Running time of reading in Zmat = %f seconds\n", MPI_Wtime() - read_begin);
	//}
	//------------------------
	//TestOut
	//------------------------

	//if (mynode == 1) {
	//	cout << "z MPI Reading matrix from node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < L_temp; kk2++) {
	//		for (int kk3 = 0; kk3 < ncol; kk3++) {
	//			cout << setw(20) << z_temp[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//}


	//=======================================================================================================================
	//Main loop {
	//1. All slave processors recursively read recursion_size rows of complete z matrix 
	//2. On each slave processors, compare distributed z matrix with recursively read z matrix and update counts and overlap
	//}
	//========================================================================================================================

	//////////////////
	//Initilization
	int recursion_size = 0;
	int recursion_time = 0; // total number of times to import partial z matrix
	int recursion_size_temp = 0;
	int recursion_size_last = 0;

	if (mynode != 0) {

		//Assume one can only use 35% of memory 
		recursion_size = (int)floor((memory * 0.35 * 1e+9) / (ncol * 8));
		if (nrow < recursion_size) recursion_size = nrow;

		if (nrow % recursion_size == 0) recursion_time = (1.0*nrow) / recursion_size;
		if (nrow % recursion_size != 0) recursion_time = (int)floor(1.0*nrow / recursion_size) + 1;

		recursion_size_last = nrow - (recursion_time - 1)*recursion_size;
		if (recursion_size_last < 10) cout<<"CAUTION!!! recursion_size_last in Zmat_ultra is less than 10!!!"<<endl;
	}

	//if (mynode == 0 || mynode == 1) cout << "Zmat recursion_time is " << recursion_time << " and recursion_size_last is " << recursion_size_last << " and recursion_size is " << recursion_size << " at node " << mynode << endl;
	////////////////

	List_FHDI summary(L_temp); // list for each pattern (both observed and missing)

	std::vector<int> remove_list;// list to remove deplicate patterns

	std::string *z_temp_cn = NULL; // string format of distributed z matrix


	if (mynode != 0 && mynode != (totalnodes - 1)) {

		z_temp_cn = new string[numWorkPerProc];
		Trans(z_temp, numWorkPerProc, ncol, z_temp_cn);

		//for (int t = 0; t < numWorkPerProc; t++) {
		//	cout<<"z_temp["<<t<<"]: "<< z_temp_cn[t]<<endl;
		//}

	}
	if (mynode == (totalnodes - 1)) {

		z_temp_cn = new string[numWorkLocalLast];
		Trans(z_temp, numWorkLocalLast, ncol, z_temp_cn);

		//cout << "summary matrix from node " << mynode << endl;
		//for (int kk2 = 0; kk2 < numWorkLocalLast; kk2++) {
		//	for (int kk3 = 0; kk3 < (2 + max_overlap_size); kk3++) {
		//		cout << setw(20) << summary[kk2][kk3];
		//	}
		//	cout << endl;
		//}
		//for (int t = 0; t < numWorkLocalLast; t++) {
		//	cout << "z_temp_last[" << t << "]: " << z_temp_cn[t] << endl;
		//}
	}

	// ==================== Main loop of recursive read ===============================================
	// The rule is that we only keep the unique pattern who occurs firstly!!!!!

	double recursion_begin = MPI_Wtime();
	int counter = 0;
	int startpoint_recur = 0;
	int endpoint_recur = 0;

	for (int i = 0; i < recursion_time; i++) {
		//cout<<"Recursion at i = "<<i<<" at node "<<mynode<<endl;

		if (i != (recursion_time - 1)) recursion_size_temp = recursion_size;
		if (i == (recursion_time - 1)) recursion_size_temp = recursion_size_last;
		//cout << "recursion_size_temp is " << recursion_size_temp << " at i = " << i << endl;

		double** z_recursion = NULL;
		double* array_temp = NULL;

		//boundary index of recursively read z matrix
		startpoint_recur = 0;
		endpoint_recur = 0;

		if (mynode != 0) {

			z_recursion = New_dMatrix(recursion_size_temp, ncol);
			array_temp = new double[ncol];

			if (i != (recursion_time - 1)) {
				startpoint_recur = i*recursion_size;
				endpoint_recur = (i + 1)*recursion_size;
			}

			if (i == (recursion_time - 1)) {
				startpoint_recur = i*recursion_size;
				endpoint_recur = i*recursion_size + recursion_size_last;
			}
		}

		//cout << "Recusive startpoint_recur: " << startpoint_recur << "; recusive endpoint_recur: " << endpoint_recur << " at recursion " << i <<" at node "<<mynode<<endl;
		double recursion1 = MPI_Wtime();
		counter = 0;
		for (int k = startpoint_recur; k < endpoint_recur; k++) {

			//******************
			//******************
			//******************
			//MPI_In_datz(nrow, ncol, k, fh_datz_read, array_temp);// read distributed z matrix row by row
			MPI_In_uox_mox(ncol, k, fh_datz, array_temp);//read distributed z matrix row by row. Note that fh_datz must be written column-wisely

			for (int m = 0; m < ncol; m++) {
				z_recursion[counter][m] = array_temp[m];
			}
			counter++;
		}

		//if (i ==0 || i==1) {
		//	cout << "YYC Running time of rucusion1 = " << MPI_Wtime() - recursion1 <<" at i="<<i<< " at node " << mynode << endl;
		//	printf("Yicheng Running time of main recusion in Zmat = %f seconds\n", MPI_Wtime() - recursion_begin);
		//}
		//if (mynode == 1) {
		//	cout << "z matrix from recusion " << i << " at node " << mynode << endl;
		//	for (int kk2 = 0; kk2 < recursion_size_temp; kk2++) {
		//		for (int kk3 = 0; kk3 < ncol; kk3++) {
		//			cout << setw(20) << z_recursion[kk2][kk3];
		//		}
		//		cout << endl;
		//	}
		//}
		double recursion2 = MPI_Wtime();
        if (mynode != 0) {
			std::string *z_recursion_cn = new std::string[recursion_size_temp]; //declaration of concatenated string vector of z
			Trans(z_recursion, recursion_size_temp, ncol, z_recursion_cn);

			//---------------------
			//Comparison process
			//---------------------
			std::string s_temp;
			std::vector<double> summary_buffer;
			double index_temp = 0.0;

			for (int k = 0; k < L_temp; k++) {
				summary_buffer.clear();
				s_temp = z_temp_cn[k]; //get a string 

				for (int l = 0; l < recursion_size_temp; l++) {

					if (s_temp.compare(z_recursion_cn[l]) == 0) {
						// If the global index of resursive z matrix < global index of disrtibuted z matrix
						// It means it is not the first unique one, we need to remove it afterwards
						if ( (k + startpoint + 1) > (l + startpoint_recur + 1) ) {
							remove_list.push_back(k); //Note that one should add local index here for removal later. 
						}

						index_temp = (double)(l + startpoint_recur + 1);
						summary_buffer.push_back(index_temp);

					}

				}

				int i_oloc_temp = (int)summary_buffer.size();
				summary.put_block_yicheng(k, i_oloc_temp, summary_buffer);

			}//end of compare

			delete[] z_recursion_cn;
		}// end of slave processors

		//if (i == 0 || i == 1) {
		//	cout << "YYC Running time of rucusion2 = " << MPI_Wtime() - recursion2 << " at i=" << i << " at node " << mynode << endl;
		//	printf("Yicheng Running time of main recusion in Zmat = %f seconds\n", MPI_Wtime() - recursion_begin);
		//}


		if (mynode != 0) {
			delete[] array_temp;
			Del_dMatrix(z_recursion, recursion_size_temp, ncol);
		}
	}//end of main
	
	//=================================================================================================


	//if (mynode == 2) {
	//	cout << "remove_list from node " << mynode << endl;
	//	for (int j = 0; j < remove_list.size(); j++) {
	//		cout<<"remove_list["<<j<<"]: "<< remove_list[j]<<endl;
	//	}
	//	cout << "summary matrix from node " << mynode << endl;
	//	summary.print_List_FHDI_yicheng();
	//}

	//=====================================================================================
	//Remove overlapped z matrix on slave processors 
	//from distributed z matrix and summary matrix (i.e., counts and overlapped locations)
	//regrading remove_list which contains global index
	//Note that z matrix will be globally unique after this step
	//======================================================================================
	int remove_size = remove_list.size();
	//cout<<"remove_size is "<< remove_size <<" at node "<<mynode<<endl;
	if (remove_size > 0) {
		for (int k = 0; k < remove_size; k++) {
			//remove overlapped ones from disrtibuted matrix
			for (int j = 0; j < ncol;j++) {
				z_temp[remove_list[k]][j] = 0.0;
			}
		}
	}

	//if (mynode == 2) {
	//	cout << "z_temp remove matrix from node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < L_temp; kk2++) {
	//		for (int kk3 = 0; kk3 < ncol; kk3++) {
	//			cout << setw(20) << z_temp[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//}

	//if (mynode == 1 || mynode == 0) {
	//	cout << "Yicheng Running time of main recusion in Zmat = " << MPI_Wtime() - recursion_begin <<" at node "<<mynode<< endl;
	//	printf("Yicheng Running time of main recusion in Zmat = %f seconds\n", MPI_Wtime() - recursion_begin);
	//}
	//MPI_Barrier(MPI_COMM_WORLD);
	//=====================================================================================
	//Split uox and mox from unique z matrix on slave processors
	//======================================================================================
	double split_begin = MPI_Wtime();
	int i_ol_temp = 0; //number of distributed uox
	int i_ml_temp = 0; //number of distributed mox
	double** uox_temp = NULL;
	double** mox_temp = NULL;

	std::vector<int> ml;// local locations of missing patterns in z; It is also local locations of mox in z
	std::vector<int> ol;// local locations of observed patterns in z; It is also local locations of uox in z

	if (mynode != 0) {

		double d_temp = 0.0;
		int zero_count = 0;

		for (int i_row = 0; i_row < L_temp; i_row++) {
			d_temp = 1.0;
			zero_count = 0;

			for (int i_col = 0; i_col < ncol; i_col++) {
				if ( fabs(z_temp[i_row][i_col]) < 1e-15) { d_temp = 0.0; zero_count++; } //found zero, i.e. missing cell
			}

			if (zero_count == ncol) {
				//cout<<"The "<< i_row << "th row of unique distributed z matrix is skipped because it is all 0s (overlapped ones) at node "<<mynode<<endl;
				continue;
			}

			if (fabs_FHDI(d_temp) > 1e-15) //this row has no missing cells

			{
				ol.push_back(i_row);
			} //actual number of the row having no missing cells



			if (fabs_FHDI(d_temp) < 1e-15) //this row has AT LEAST one missing cells

			{
				ml.push_back(i_row);
			}

		}

		i_ol_temp = ol.size();
		i_ml_temp = ml.size();

		uox_temp = New_dMatrix(i_ol_temp, ncol);
		mox_temp = New_dMatrix(i_ml_temp, ncol);

		for (int k = 0; k < i_ol_temp; k++) {
			for (int j = 0; j < ncol;j++) {
				uox_temp[k][j] = z_temp[ol[k]][j];
			}
		}

		for (int k = 0; k < i_ml_temp; k++) {
			for (int j = 0; j < ncol;j++) {
				mox_temp[k][j] = z_temp[ml[k]][j];
			}
		}

		MPI_Send(&i_ol_temp, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&i_ml_temp, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);

		//if (mynode == 1) {
		//	cout<<"ol at node "<<mynode<<endl;
		//	for (int m = 0; m < ol.size(); m++) {
		//		cout << "ol[" << m << "]: " << ol[m] << endl;
		//	}
		//	cout << "ml at node " << mynode << endl;
		//	for (int m = 0; m < ml.size(); m++) {
		//		cout << "ml[" << m << "]: " << ml[m] << endl;
		//	}
		//}

		//if (mynode == 3) {
		//	cout << "uox_temp matrix from node " << mynode << endl;
		//	for (int kk2 = 0; kk2 < i_ol_temp; kk2++) {
		//		for (int kk3 = 0; kk3 < ncol; kk3++) {
		//			cout << setw(20) << uox_temp[kk2][kk3];
		//		}
		//		cout << endl;
		//	}
		//	cout << "mox_temp matrix from node " << mynode << endl;
		//	for (int kk2 = 0; kk2 < i_ml_temp; kk2++) {
		//		for (int kk3 = 0; kk3 < ncol; kk3++) {
		//			cout << setw(20) << mox_temp[kk2][kk3];
		//		}
		//		cout << endl;
		//	}
		//}

	}//end of slave processors
	//if (mynode == 1 || mynode == 0) {
	//	cout << "Yicheng Running time of split1 in Zmat = "<< MPI_Wtime() - split_begin <<" at node "<<mynode<< endl;
	//	printf("Yicheng Running time of split1 in Zmat = %f seconds\n", MPI_Wtime() - split_begin);
	//}
	 //-----------------------------------------------------------------------
	 //Note that to write uox and mox on slave processors to local storage,
	 //we have to know the dimensions of uox and mox on each processor to compute
	 //offset value for MPI_I/O
	 //-----------------------------------------------------------------------
	double split_begin2 = MPI_Wtime();
	std::vector<int> ol_size;//number of uox on all slave processors
	std::vector<int> ml_size;//number of mox on all slave processors

	if (mynode == 0) {
		int i_ol_temp_recv = 0;
		int i_ml_temp_recv = 0;
		ol_size.push_back(i_ol_temp_recv);//Initialization
		ml_size.push_back(i_ml_temp_recv);//Initialization

		for (int j = 1; j < totalnodes; j++) {
			MPI_Recv(&i_ol_temp_recv, 1, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
			ol_size.push_back(i_ol_temp_recv);

			MPI_Recv(&i_ml_temp_recv, 1, MPI_INT, j, 2, MPI_COMM_WORLD, &status);
			ml_size.push_back(i_ml_temp_recv);
			//cout<<"i_ol_temp_recv is "<< i_ol_temp_recv <<" and i_ml_temp_recv is "<< i_ml_temp_recv <<" at node "<<j<<endl;
		}

	}

	//if (mynode == 1 || mynode == 0) {
	//	cout << "Yicheng Running time of split2 in Zmat is " << MPI_Wtime() - split_begin2 << " at node " << mynode << endl;
	//}


	double split_begin20 = MPI_Wtime();
	if (mynode != 0) {
		ol_size.resize(totalnodes);
		ml_size.resize(totalnodes);
	}


	MPI_Bcast(&ol_size[0], totalnodes, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ml_size[0], totalnodes, MPI_INT, 0, MPI_COMM_WORLD);

	//if (mynode == 0) {
	//	for (int t = 0; t < totalnodes; t++) {
	//		cout << "ol_size[" << t << "] : "<< ol_size[t] << endl;
	//	}

	//	for (int t = 0; t < totalnodes; t++) {
	//		cout << "ml_size[" << t << "] : " << ml_size[t] << endl;
	//	}
	//}

	for (int k = 1; k < totalnodes; k++) {
		uox_size = uox_size + ol_size[k];
		mox_size = mox_size + ml_size[k];
	}
	
	//cout << "uox_size is " << uox_size << " and mox_size is " << mox_size << " at node " << mynode << endl;

	//cout << "i_ol_temp is " << i_ol_temp << " and i_ml_temp is " << i_ml_temp << " at node " << mynode << endl;

	//Write uox and mox row-wisely on slave processors to the same files concurrently
	MPI_Out_uox_mox(i_ol_temp, ncol, ol_size, fh_uox, uox_temp);
	MPI_Out_uox_mox(i_ml_temp, ncol, ml_size, fh_mox, mox_temp);
	MPI_Barrier(MPI_COMM_WORLD);
	//if (mynode == 1 || mynode == 0) {
	//	cout << "Yicheng Running time of wrting uox and mox to hard drive in Zmat is " << MPI_Wtime() - split_begin20 << " at node " << mynode << endl;
	//}



	//-----------------------------------------------------------------------
	//Split summary to information of uox and mox
	//-----------------------------------------------------------------------

	//1. All processors must have size of uox_infor and mox_infor firstly
	double split_begin3 = MPI_Wtime();
	if (mynode != 0) {

		std::vector<int> uox_size_temp;//number of occurance for distributed uox on slave prcessors
		int i_size_temp = 0;
		for (int k = 0; k < i_ol_temp;k++) {
			i_size_temp = 0;
			summary.get_a_row_size(ol[k], i_size_temp);
			uox_size_temp.push_back(i_size_temp);
		}

		std::vector<int> mox_size_temp;//number of occurance for distributed mox on slave prcessors
		for (int k = 0; k < i_ml_temp;k++) {
			i_size_temp = 0;
			summary.get_a_row_size(ml[k], i_size_temp);
			mox_size_temp.push_back(i_size_temp);
		}

		MPI_Send(&uox_size_temp[0], uox_size_temp.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&mox_size_temp[0], mox_size_temp.size(), MPI_INT, 0, 2, MPI_COMM_WORLD);
	}//end of slave

	std::vector<int> uox_size_total; //number of occurance for all uox
	std::vector<int> mox_size_total; //number of occurance for all mox

	if (mynode == 0) {

		for (int j = 1; j < totalnodes; j++) {
			
			int i_ol_size = ol_size[j]; //number of distributed uox on jth processor. ol_size[0] = 0 for master processor
			int i_ml_size = ml_size[j]; //number of distributed mox on jth processor. ml_size[0] = 0 for master processor

			std::vector<int> uox_size_recv(i_ol_size);
			std::vector<int> mox_size_recv(i_ml_size);

			MPI_Recv(&uox_size_recv[0], i_ol_size, MPI_INT, j, 1, MPI_COMM_WORLD, &status);

			for (int k = 0; k < i_ol_size; k++) {
				uox_size_total.push_back(uox_size_recv[k]);
			}

			MPI_Recv(&mox_size_recv[0], i_ml_size, MPI_INT, j, 2, MPI_COMM_WORLD, &status);

			for (int k = 0; k < i_ml_size; k++) {
				mox_size_total.push_back(mox_size_recv[k]);
			}

		}//end of mynode

		//cout<<"Node 0: "<<endl;
		//for (int i = 0; i < uox_size_total.size(); i++) {
		//	cout << "uox_size_total[" << i << "]: " << uox_size_total[i] << endl;
		//}
		//for (int i = 0; i < mox_size_total.size(); i++) {
		//	cout << "mox_size_total[" << i << "]: " << mox_size_total[i] << endl;
		//}
	}

	if (mynode != 0) {
		uox_size_total.resize(uox_size);
		mox_size_total.resize(mox_size);
	}

	MPI_Bcast(&uox_size_total[0], uox_size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mox_size_total[0], mox_size, MPI_INT, 0, MPI_COMM_WORLD);

	if (uox_size_total.size() != uox_size) { cout << "uox_size_total is incorrect in Zmat!" << endl; }
	if (mox_size_total.size() != mox_size) { cout << "mox_size_total is incorrect in Zmat!" << endl; }

	int i_ol_size = 0; // total number of observed patterns. Note i_ol_size = number of uox if remove list is 0 for high-dimensional data
	for (int t = 0; t < uox_size; t++) {
		i_ol_size = i_ol_size + uox_size_total[t];
	}

	int i_ml_size = 0; // total number of missing patterns. Note i_ml_size = number of mox if remove list is 0 for high-dimensional data
	for (int t = 0; t < mox_size; t++) {
		i_ml_size = i_ml_size + mox_size_total[t];
	}

	if ((i_ol_size + i_ml_size) != nrow) { cout<<"ERROR! i_ol_size + i_ml_size ! nrow in Zmat function!"<<endl; }

	//cout<<"i_ol_size is "<< i_ol_size <<" and i_ml_size is "<< i_ml_size <<" at node "<<mynode<<endl;

	//if (mynode == 1) {
	//	for (int i = 0; i < uox_size; i++) {
	//		cout << "uox_size_total[" << i << "]: " << uox_size_total[i] << endl;
	//	}
	//	for (int i = 0; i < mox_size; i++) {
	//		cout << "mox_size_total[" << i << "]: " << mox_size_total[i] << endl;
	//	}
	//}

	MPI_Barrier(MPI_COMM_WORLD);


	//Send index of distributed patterns in z  
	if (mynode != 0) {
		std::vector<double> uox_infor_temp; // index of all patterns in z of distributed uox
		std::vector<double> uox_infor_buffer; // index of all patterns in z of one uox

		for (int k = 0; k < i_ol_temp; k++) {
			uox_infor_buffer.clear();
			summary.get_block_yicheng(ol[k], uox_infor_buffer);
			for (int t = 0; t < uox_infor_buffer.size(); t++) {
				uox_infor_temp.push_back(uox_infor_buffer[t]);
			}
		}

		std::vector<double> mox_infor_temp; // index of all patterns in z of distributed mox
		std::vector<double> mox_infor_buffer; // index of all patterns in z of one mox

		for (int k = 0; k < i_ml_temp; k++) {
			mox_infor_buffer.clear();
			summary.get_block_yicheng(ml[k], mox_infor_buffer);
			for (int t = 0; t < mox_infor_buffer.size(); t++) {
				mox_infor_temp.push_back(mox_infor_buffer[t]);
			}
		}

		//cout<<"uox_infor_temp size is "<< uox_infor_temp.size() <<" at node "<<mynode<<endl;
		//cout <<"mox_infor_temp size is " << mox_infor_temp.size() << " at node " << mynode << endl;

		MPI_Send(&uox_infor_temp[0], uox_infor_temp.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&mox_infor_temp[0], mox_infor_temp.size(), MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);

	}


	//if (mynode == 1 || mynode == 0) {
	//	cout << "Yicheng Running time of split3 in Zmat = " << MPI_Wtime() - split_begin3 <<" at node "<<mynode<< endl;
	//	printf("Yicheng Running time of split3 in Zmat = %f seconds\n", MPI_Wtime() - split_begin3);
	//}

	//-----------------------------------------------------------------------
	//All_gather information of mox and uox on master and broadcast to slaves
	//-----------------------------------------------------------------------
	double gather_begin = MPI_Wtime();

	if ((i_ol_size + i_ml_size) != nrow) { cout << "ERROR! i_ol_size + i_ml_size ! nrow in Zmat function!" << endl; }

	std::vector<double> uox_infor_total; //index of all uox in z
	std::vector<double> mox_infor_total; //index of all mox in z

	if (mynode == 0) {
		int i_size_cell = 0;
		int start_point = 0; //start for uox_size_total on each slave processor
		int end_point = 0; // end for uox_size_total on each slave processor

		for (int j = 1; j < totalnodes; j++) {

			//Receive cells of uox
			i_size_cell = 0;
			start_point = 0;
			end_point = 0;

			for (int k = 0; k < j; k++) {
				start_point = start_point + ol_size[k];
			}

			end_point = start_point + ol_size[j];

			//cout << "uox start_point is " << start_point <<" and end_point is "<< end_point << " at node " << j << endl;

			for (int t = start_point; t < end_point; t++) {
				i_size_cell = i_size_cell + uox_size_total[t];
			}

			//cout<<"uox i_size_cell is "<< i_size_cell <<" at node "<<j<<endl;
			std::vector<double> uox_infor_recv(i_size_cell);
			MPI_Recv(&uox_infor_recv[0], i_size_cell, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);

			for (int k = 0; k < i_size_cell; k++) {
				uox_infor_total.push_back(uox_infor_recv[k]);
			}

			//Receive cells of mox
			i_size_cell = 0;
			start_point = 0;
			end_point = 0;

			for (int k = 0; k < j; k++) {
				start_point = start_point + ml_size[k];
			}

			end_point = start_point + ml_size[j];

			//cout << "mox start_point is " << start_point << " and end_point is " << end_point << " at node " << j << endl;

			for (int t = start_point; t < end_point; t++) {
				i_size_cell = i_size_cell + mox_size_total[t];
			}
			
			//cout << "mox i_size_cell is " << i_size_cell << " at node " << j << endl;

			std::vector<double> mox_infor_recv(i_size_cell);
			MPI_Recv(&mox_infor_recv[0], i_size_cell, MPI_DOUBLE, j, 2, MPI_COMM_WORLD, &status);

			for (int k = 0; k < i_size_cell; k++) {
				mox_infor_total.push_back(mox_infor_recv[k]);
			}

		}
		
	}//end of master

	if (mynode != 0) { 
		uox_infor_total.resize(i_ol_size);
		mox_infor_total.resize(i_ml_size);
	}

	MPI_Bcast(&uox_infor_total[0], i_ol_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mox_infor_total[0], i_ml_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	//if (mynode == 1) {
	//	for (int i = 0; i < i_ol_size; i++) {
	//		cout << "uox_infor_total[" << i << "]: " << uox_infor_total[i] << endl;
	//	}

	//	for (int i = 0; i < i_ml_size; i++) {
	//		cout << "mox_infor_total[" << i << "]: " << mox_infor_total[i] << endl;
	//	}
	//}

	//if (mynode == 1 || mynode == 0) {
	//	cout << "Yicheng Running time of gathering in Zmat = " << MPI_Wtime() - gather_begin <<" at node " <<mynode<< endl;
	//	printf("Yicheng Running time of gathering in Zmat = %f seconds\n", MPI_Wtime() - gather_begin);
	//}

	//------------------------------------------------------------
	//Transfer uox infor and mox infor to List_FHDI data structure
	//------------------------------------------------------------

	uox_infor.initialize(uox_size); // index list of uox in z 
	mox_infor.initialize(mox_size); // index list of mox in z 

	std::vector<double> infor_buffer; // hold one row of infor
	int counter2 = 0;
	int i_size_temp2 = 0; 
	for (int k = 0; k < uox_size; k++) {
		i_size_temp2 = 0;
		infor_buffer.clear();
		i_size_temp2 = uox_size_total[k];

		for (int j = 0; j < i_size_temp2; j++) {
			infor_buffer.push_back(uox_infor_total[counter2]);
			counter2++;
		}
		uox_infor.put_block_yicheng(k, i_size_temp2, infor_buffer);
	}

	if (counter2 != i_ol_size) { cout<<"ERROR! counter2 in Zmat is incorrect!"<<endl; }

	counter2 = 0;
	i_size_temp2 = 0;
	for (int k = 0; k < mox_size; k++) {
		i_size_temp2 = 0;
		infor_buffer.clear();
		i_size_temp2 = mox_size_total[k];

		for (int j = 0; j < i_size_temp2; j++) {
			infor_buffer.push_back(mox_infor_total[counter2]);
			counter2++;
		}
		mox_infor.put_block_yicheng(k, i_size_temp2, infor_buffer);
	}

	if (counter2 != i_ml_size) { cout << "ERROR! counter2 in Zmat is incorrect!" << endl; }

	//if (mynode == 1) {
	//	cout << "uox_infor from node " << mynode << endl;
	//	uox_infor.print_List_FHDI_yicheng();

	//	cout << "mox_infor from node " << mynode << endl;
	//	mox_infor.print_List_FHDI_yicheng();
	//}

	//if (mynode == 1) {
	//		cout << "uox_info matrix from node " << mynode << endl;
	//		for (int kk2 = 0; kk2 < nrow; kk2++) {
	//			for (int kk3 = 0; kk3 < (2 + max_overlap_size); kk3++) {
	//				cout << setw(20) << uox_info[kk2][kk3];
	//			}
	//			cout << endl;
	//		}
	//		cout << "mox_info matrix from node " << mynode << endl;
	//		for (int kk2 = 0; kk2 < nrow; kk2++) {
	//			for (int kk3 = 0; kk3 < (2 + max_overlap_size); kk3++) {
	//				cout << setw(20) << mox_info[kk2][kk3];
	//			}
	//			cout << endl;
	//		}
	//}


	//-----------------------
	//Deallocation
	//------------------------
	if (mynode!=0) delete[] z_temp_cn;

	Del_dMatrix(z_temp, L_temp, ncol);


	if (mynode != 0) { 
		Del_dMatrix(uox_temp, i_ol_temp, ncol); 
		Del_dMatrix(mox_temp, i_ml_temp, ncol);
	}


	return 1;

}
