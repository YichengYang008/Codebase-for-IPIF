#include "Fractional_Hot_Deck_Imputation_ultra_MPI.cc"


void FHDI_Neighbor_ultra_cpp(
	MPI_File fh_binary_daty_row, const int nrow, const int ncol, const int i_merge,
	int nrow_uox, int nrow_mox,
	std::vector<double> 	 jp_prob,
	std::string s_M, const int i_M, List_FHDI &List_nU, std::vector<int> ol,
	List_FHDI &uox_infor, List_FHDI &mox_infor, double** simp_fmat_FHDI,
	ofstream& TestOut_Slave1)
	//Description=========================================
	// perform
	// Fractional Hot Deck Imputation Only
	// 
	// Algorithm: FEFI of Dr Jae Kwang. Kim and FHDI of Dr Jong Ho. Im
	//
	// original R code: Dr. Im, J. and Dr. Kim, J. 
	// c++ code: 		Yicheng Yang and Dr. Cho, In-Ho 
	// All rights reserved
	// 
	// updated: August 16, 2021
	//----------------------------------------------------
	//IN   : MPI_File fh_binary_daty_row  = raw daty written row-wisely
	//IN   : const int nrow           = number of rows of daty
	//IN   : const int ncol           = number of columns of daty
	//IN   : const int i_merge        = control of random number generator
	//IN   : int nrow_uox             = number of rows of uox
	//IN   : int nrow_mox             = number of rows of mox
	//IN   : MPI_File fh_uox	      = file string to read uox
	//IN   : MPI_File fh_mox	      = file string to read mox
	//IN   : vector<double> jp_prob   = joint probability of all uox
	//IN   : string s_M               = "FHDI"
	//IN   : const int i_M            = user-defined number of donors
	//IN   : List_FHDI List_nU        = donor list
	//IN   : vector<int> ol           = actual locations of fully observed rows in z matrix
	//IN   : List_FHDI uox_infor      = Actual index list of uox in z
	//IN   : List_FHDI mox_infor      = Actual index list of mox in z

	//OUT  : double** simp_fmat_FHDI  = matrix of 4 columns: global id + sampling weight + wij + fwij
	//OUT  : MPI_File fh_fmat_FHDI    = file string to write imputed matrix in size of (4+ncol) columns: global id + sampling weight + wij + fwij + imputed daty



	//====================================================
{
	//-----------------
	//random location using uniform distribution   
	//using Numerical Recipes of Press et al 2007. 
	//-----------------
	//Ran_FHDI myran(1); 	//not used for CRAN Compatibility
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	double FHDI_begin = MPI_Wtime();

	//open file string to read mox
	int success = 0;
	MPI_File fh_mox_reading;
	//open file to write imputed values
	MPI_File fh_fmat_FHDI; 

	success = MPI_File_open(MPI_COMM_WORLD, "./mox_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_mox_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of mox fail to open the file!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./fmat_FHDI_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_fmat_FHDI);
	if (success != MPI_SUCCESS) cout << "MPI I/O of mox fail to open the file!" << endl;

	//-------------
	//sample weight (default is 1)
	//id array (default is row number)
	//-------------
	double* w = new double[nrow];
	int* id = new int[nrow];
	for (int i = 0; i<nrow; i++)
	{
		w[i] = 1.0;
		id[i] = i + 1; //ACTUAL id
	}

	//--------------------
	//Distribute nrow_mox
	//--------------------
	const int L = nrow_mox;
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
	if (mynode != (totalnodes - 1)) L_temp = numWorkPerProc;
	if (mynode == (totalnodes - 1)) L_temp = numWorkLocalLast;

	//cout << "FHDI start is " << startpoint << " and endpoint is " << endpoint << " at node " << mynode << endl;
	//============================

	//if (mynode == 0) {
	//	cout<<"ol: "<<endl;
	//	for (int t = 0; t < i_ol_size; t++) {
	//		cout<<"ol["<<t<<"]: "<<ol[t]<<endl;
	//	}
	//}

	//-----------------------------------------
	//Distribute rows of fully observed rows
	//------------------------------------------

	int i_ol_size = 0;
	i_ol_size = ol.size();

	int numWorkPerProc_ol = (int)floor(1.0*i_ol_size / (1.0*totalnodes - 1));
	int numWorkLocalLast_ol = i_ol_size - numWorkPerProc_ol * (totalnodes - 2);


	int startpoint_ol = 0;
	int endpoint_ol = 0;

	if (mynode != 0 && mynode != (totalnodes - 1)) {
		startpoint_ol = (mynode - 1)*numWorkPerProc_ol;
		endpoint_ol = mynode*numWorkPerProc_ol;
	}

	if (mynode == (totalnodes - 1)) {
		startpoint_ol = (mynode - 1)*numWorkPerProc_ol;
		endpoint_ol = (mynode - 1)*numWorkPerProc_ol + numWorkLocalLast_ol;
	}

	int ol_temp = 0;
	if (mynode != (totalnodes - 1)) ol_temp = numWorkPerProc_ol;
	if (mynode == (totalnodes - 1)) ol_temp = numWorkLocalLast_ol;


	//===============================================================
	//Compute size of distributed result matrix on slave processors 
	//
	//Note that it should be distributed ol matrix + distributed imputed matrix
	//because one slave processor can not hold all ol matrix
	//===============================================================
	double FHDI_main = MPI_Wtime();

	//Initilization
	std::vector<int> loc_srst_nl; // buffer to hold donors of each mox
	int i_row_fmat_FHDI_total = 0;// total number of rows of distributed result matrix on slave processors
	int uox_sum = 0;
	int temp_size = 0;

	double** fmat_FHDI_temp = NULL; // result matrix of (4 + ncol columns): global id + sampling weight + wij + fwij + imputed y of distributed mox
	double ** simp_fmat_FHDI_temp = NULL;// matrix of 4 columns: global id + sampling weight + wij + fwij of distributed mox


	int ncol_fmat = 4; // 4 columns: global id + sampling weight + wij + fwij 
	if (mynode != 0) {

		//Compute size of distributed result matrix on slave processors 
		int i_temp2 = 0;
		for (int k = startpoint; k < endpoint; k++) {
			loc_srst_nl.clear();
			uox_sum = 0;
			List_nU.get_block_yicheng(k, loc_srst_nl);// get donors of distributed mox
			temp_size = loc_srst_nl.size();
			for (int t = 0; t < temp_size; t++) {
				i_temp2 = 0;
				uox_infor.get_a_row_size(loc_srst_nl[t] - 1, i_temp2);
				uox_sum = uox_sum + i_temp2;// second column of uox_info_final is occurance 
			}

			if (uox_sum > i_M) uox_sum = i_M; // if larger than user-defined donor number i_M, use i_M

			i_temp2 = 0;
			mox_infor.get_a_row_size(k, i_temp2);

			i_row_fmat_FHDI_total = i_row_fmat_FHDI_total + i_temp2 * uox_sum;
		}

		//cout << "first i_row_fmat_FHDI_total is " << i_row_fmat_FHDI_total << " at node " << mynode << endl;

		i_row_fmat_FHDI_total = i_row_fmat_FHDI_total + ol_temp; // note that one should add distributed ol matrix

																 //cout << "FHDI i_row_fmat_FHDI_total is " << i_row_fmat_FHDI_total << " at node " << mynode << endl;

		fmat_FHDI_temp = New_dMatrix(i_row_fmat_FHDI_total, ncol_fmat + ncol); //return matrix from FHDI
		simp_fmat_FHDI_temp = New_dMatrix(i_row_fmat_FHDI_total, ncol_fmat); //return matrix from FHDI
		Fill_dMatrix(fmat_FHDI_temp, i_row_fmat_FHDI_total, ncol_fmat + ncol, 0.0);
		Fill_dMatrix(simp_fmat_FHDI_temp, i_row_fmat_FHDI_total, ncol_fmat, 0.0);

	}

	//----------------------------------------------
	// Global counter for distributed result matrix
	int counter_fmat = 0;
	//--------------------------------------------

	//------------------------------------------------------------------------------
	// Add distributed fully observed matrixs to distributed result matrix firstly
	// One can't add to single processor to avoid memory issue
	//-----------------------------------------------------------------------------
	double* array_temp = new double[ncol];// temp buffer
	int i_ol_temp = 0;
	for (int i = startpoint_ol; i < endpoint_ol; i++) {
		i_ol_temp = ol[i] - 1;
		//MPI_In_daty_row(nrow, ncol, i_ol_temp, fh_binary_daty, array_temp);
		MPI_In_uox_mox(ncol, i_ol_temp, fh_binary_daty_row, array_temp);// read distributed z matrix row by row

		simp_fmat_FHDI_temp[counter_fmat][0] = ol[i]; // global id
		simp_fmat_FHDI_temp[counter_fmat][1] = 1.0; // sampling weight
		simp_fmat_FHDI_temp[counter_fmat][2] = 1.0; // wij
		simp_fmat_FHDI_temp[counter_fmat][3] = 1.0; // fwij

		fmat_FHDI_temp[counter_fmat][0] = ol[i]; // global id
		fmat_FHDI_temp[counter_fmat][1] = 1.0; // sampling weight
		fmat_FHDI_temp[counter_fmat][2] = 1.0; // wij
		fmat_FHDI_temp[counter_fmat][3] = 1.0; // fwij

		for (int k = 0; k < ncol; k++) {
			fmat_FHDI_temp[counter_fmat][k + 4] = array_temp[k]; // original fully observed daty
		}

		//cout<<"counter_fmat = "<< counter_fmat <<endl;
		counter_fmat++;
	}

	//Deallocation
	delete[] array_temp;

	//cout << "counter_fmat is " << counter_fmat << " at node " << mynode << endl;

	//if (mynode == 1) {
	//	cout << "simp_fmat_FHDI_temp at node " << mynode << endl;
	//	for (int j = 0; j < i_row_fmat_FHDI_total; j++) {
	//		for (int t = 0; t < ncol_fmat; t++) {
	//			cout << setw(20) << simp_fmat_FHDI_temp[j][t];
	//		}
	//		cout << endl;
	//	}

	//	cout << "fmat_FHDI_temp at node " << mynode << endl;
	//	for (int j = 0; j < i_row_fmat_FHDI_total; j++) {
	//		for (int t = 0; t < (ncol_fmat + ncol); t++) {
	//			cout << setw(20) << fmat_FHDI_temp[j][t];
	//		}
	//		cout << endl;
	//	}
	//}


	//============================================================
	//Main iteration to add imputed matrix of distributed mox
	//============================================================

	//Initilization
	std::vector<int> v_mxl; //a row's columns having observed cells  
	std::vector<int> v_cn_z_i; // actual locations of a mox in z matrix
	std::vector<int> v_obsg0; // buffer to hold all possible donors of a mox in z matrix (original)
	std::vector<int> v_obsg; // buffer to hold all possible donors of a mox in z matrix (random shuffled)
	int counter = 0;
	double* mox_temp = new double[ncol];// buffer to read a mox

	std::vector<double> v_cp; //normalzied joint probability of donors

	std::vector<double> uox_buffer;
	std::vector<double> mox_buffer;
	for (int i = startpoint; i < endpoint; i++) {
		//cout<<"Imputation at i = "<<i<<" at node "<<mynode<<endl;
		double FHDI_iteration1 = MPI_Wtime();
		//------
		//number of observed cells on this row
		//------
		loc_srst_nl.clear();
		v_mxl.clear();
		v_obsg0.clear();
		v_obsg.clear();
		v_cn_z_i.clear();

		// Read ith row of mox
		MPI_In_uox_mox(ncol, i, fh_mox_reading, mox_temp);

		//if (mynode == 1) {
		//	cout << "mox at i = " << i << endl;
		//	for (int t = 0; t < ncol; t++) {
		//		cout << setw(20) << mox_temp[t];
		//	}
		//	cout << endl;
		//}


		for (int j = 0; j<ncol; j++)
		{
			if (mox_temp[j]>0)
			{
				v_mxl.push_back(j + 1); //Actual non-missing cell location 
			}
		}

		// actual locations of ith mox in z matrix
		mox_buffer.clear();
		mox_infor.get_block_yicheng(i, mox_buffer);
		for (int k = 0; k < mox_buffer.size(); k++) {
			v_cn_z_i.push_back((int)mox_buffer[k]);
		}


		//for (int k = 2; k < (2 + max_overlap_size); k++) {
		//	if (mox_info_final[i][k] != 0) {
		//		v_cn_z_i.push_back(mox_info_final[i][k]);
		//	}
		//}

		//if (mynode == 1) {
		//	cout << "v_cn_z_i at mox = " << i <<" at node "<<mynode<< endl;
		//	for (int m = 0; m < v_cn_z_i.size(); m++) {
		//		cout << setw(20) << v_cn_z_i[m];
		//	}
		//	cout << endl;
		//}

		// Get donors of ith mox
		List_nU.get_block_yicheng(i, loc_srst_nl);

		int i_size_loc_srst_nl = (int)loc_srst_nl.size();

		if (i_size_loc_srst_nl == 0) {
			cout << "Error! there is no matched cell!" << endl; return;
		}

		for (int k = 0; k < i_size_loc_srst_nl; k++) {
			uox_buffer.clear();
			uox_infor.get_block_yicheng(loc_srst_nl[k] - 1, uox_buffer);

			for (int t = 0; t < uox_buffer.size(); t++) {
				v_obsg0.push_back( (int)uox_buffer[t] );
			}
		}

		int i_size_v_obsg = v_obsg0.size();


		//--------------------------------------------------------
		//Replace sorting half-ascending and half-descending order
		//with random shuffle
		//--------------------------------------------------------

		for (int t = 0; t < i_size_v_obsg; t++) {
			v_obsg.push_back(v_obsg0[t]);
		}

		if (i_merge == 1) std::srand(time(NULL)); //turn on random seed using C++ standard rand() fn 
		if (i_merge == 0) std::srand(123);	//turn on the fixed seed //This is still platform-dependent 
		random_shuffle(v_obsg.begin(), v_obsg.end());

		//---------------------------------------------------------
		//mapping from original donor indices to shuffled indices
		//---------------------------------------------------------
		std::vector<int> v_rbsg; v_rbsg.clear();
		match_FHDI(v_obsg, v_obsg0, v_rbsg); //get loc stored in v_rbsg

		//------------------------------
		//------------------------------
		//Compute Fractional Weights (fwij)
		//Fractional weights for FEFI representing sampling w
		//------------------------------
		//------------------------------

		//------
		//calculate fractional weights
		//------
		std::vector<double> fwij(i_size_v_obsg);
		std::vector<double> d_obsp(i_size_v_obsg);
		std::vector<int> i_obsn(i_size_v_obsg);

		//Note that one have to normalize joint probability of 
		//donors firstly
		v_cp.clear(); //normalzied joint probability of donors
		for (int p = 0; p < i_size_loc_srst_nl; p++) {
			v_cp.push_back(jp_prob[loc_srst_nl[p] - 1]);
		}
		int i_size_v_cp = (int)v_cp.size();

		if (i_size_v_cp != i_size_loc_srst_nl) cout << "ERROR in computing fwij in imputation !!!!" << endl;

		//Normalization of jp
		double d_sum_v_cp = 0.0;
		for (int j = 0; j<i_size_v_cp; j++) d_sum_v_cp += v_cp[j];
		if (d_sum_v_cp != 0) {
			for (int j = 0; j < i_size_v_cp; j++) v_cp[j] = v_cp[j] / d_sum_v_cp;
		}

		//compute fwij
		counter = 0;
		int occurnace = 0;
		for (int k = 0; k < i_size_loc_srst_nl; k++) {
			occurnace = 0;
			uox_infor.get_a_row_size(loc_srst_nl[k] - 1, occurnace);
			for (int l = 0; l < occurnace; l++) {
				d_obsp[counter] = v_cp[k];// -1 for actual location 
				i_obsn[counter] = occurnace;// -1 for actual location  
				counter++;
			}
		}

		//cout << "i_size_v_cp is " << i_size_v_cp << " and i_size_loc_srst_nl is " << i_size_loc_srst_nl << " at i=" << i << " at node " << mynode << endl;
		if (counter != i_size_v_obsg) {
			cout << "ERROR in FHDI_Neighbor_ultra at i = " << i << " at node " << mynode << " where counter = " << counter << " and i_size_v_obsg = " << i_size_v_obsg << endl;
			return;
		}

		// Note d_obsp match original donor indices v_obsg0
		// We need to match v_obsg here using mapping

		for (int k = 0; k < i_size_v_obsg; k++) {
			fwij[k] = 1.0; //default for error case  
			if (i_obsn[k] != 0) fwij[k] = d_obsp[v_rbsg[k]-1] / i_obsn[v_rbsg[k] - 1];
			if (i_obsn[k] == 0) cout << "Error! zero count in obsn at k = " << k << endl;
		}




		/*if (mynode == 1) {

			cout << "v_cn_z_i at mox = " << i << endl;
			for (int m = 0; m < v_cn_z_i.size(); m++) {
				cout << setw(20) << v_cn_z_i[m];
			}
			cout << endl;

			cout << "v_obsg at mox = " << i << endl;
			for (int m = 0; m < i_size_v_obsg; m++) {
				cout << setw(20) << v_obsg[m];
			}
			cout << endl;

			cout << "fwij at mox = " << i << endl;
			for (int m = 0; m < i_size_v_obsg; m++) {
				cout << setw(20) << fwij[m];
			}
			cout << endl;

			cout << "i_obsn at mox = " << i << endl;
			for (int m = 0; m < i_size_v_obsg; m++) {
				cout << setw(20) << i_obsn[m];
			}
			cout << endl;

			cout << "d_obsp at mox = " << i << endl;
			for (int m = 0; m < i_size_v_obsg; m++) {
				cout << setw(20) << d_obsp[m];
			}
			cout << endl;
		}*/

		//if (mynode == (totalnodes - 1)) {
		//	cout << "FHDI_iteration1 time at i = " << i << " = " << MPI_Wtime() - FHDI_iteration1 << endl;
		//}

		double FHDI_iteration3 = MPI_Wtime();

		int i_mxl = (int)v_mxl.size();// number of observed variables in current mox


		if (s_M.compare("FHDI") == 0) {

			double d_myran = 0.0;
			d_myran = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);

			Fractional_Hot_Deck_Imputation_ultra(i, fh_binary_daty_row,
				ncol, nrow,
				mox_temp, i_M, i_merge,
				i_mxl,
				v_cn_z_i, v_mxl,
				v_obsg, fwij,
				d_myran, w, id, counter_fmat,
				simp_fmat_FHDI_temp, fmat_FHDI_temp);
		}

		//if (mynode == (totalnodes - 1)) {
		//	cout << "FHDI_iteration3 time at i = " << i << " = " << MPI_Wtime() - FHDI_iteration3 << endl;
		//}


		//if (mynode == (totalnodes - 1)) {
		//	cout << "FHDI new time at i = " << i << " = " << MPI_Wtime() - FHDI_iteration1 << endl;
		//}

	}//end of main loop

	 //cout << "After counter_fmat is " << counter_fmat << " at node " << mynode << endl;
	if (counter_fmat != i_row_fmat_FHDI_total) {
		cout << "ERROR!!! Global counter for fmat is incorrect in FHDI_neighbor_ultra" << endl;
		return;
	}

	//MPI_Barrier(MPI_COMM_WORLD);

	//if (mynode == 1) {
	//	cout << "simp_fmat_FHDI_temp at node " << mynode << endl;
	//	for (int j = 0; j < i_row_fmat_fhdi_total; j++) {
	//		for (int t = 0; t < ncol_fmat; t++) {
	//			cout << setw(20) << simp_fmat_FHDI_temp[j][t];
	//		}
	//		cout << endl;
	//	}

	//	cout << "fmat_FHDI_temp at node " << mynode << endl;
	//	for (int j = 0; j < i_row_fmat_FHDI_total; j++) {
	//		for (int t = 0; t < (ncol_fmat + ncol); t++) {
	//			cout << setw(20) << fmat_FHDI_temp[j][t];
	//		}
	//		cout << endl;
	//	}
	//}

	// Send distributed return matrix to master node
	if (mynode != 0) {
		//Note to write return matrix to local binary file
		//it is critical to be aware of size in advance
		MPI_Send(&i_row_fmat_FHDI_total, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
	}

	std::vector<int> i_row_fmat_FHDI_total_size;//number of rows of uox on all slave processors


	// All_gather i_row_fmat_FHDI_total
	if (mynode == 0) {
		int i_row_fmat_FHDI_total_recv = 0;
		i_row_fmat_FHDI_total_size.push_back(i_row_fmat_FHDI_total_recv);//First disp is initilized as 0

		for (int j = 1; j < totalnodes; j = j + 1) {
			if (j != (totalnodes - 1)) {
				MPI_Recv(&i_row_fmat_FHDI_total_recv, 1, MPI_INT, j, 2, MPI_COMM_WORLD, &status);

				i_row_fmat_FHDI_total_size.push_back(i_row_fmat_FHDI_total_recv);
			}

			if (j == (totalnodes - 1)) {
				MPI_Recv(&i_row_fmat_FHDI_total_recv, 1, MPI_INT, j, 2, MPI_COMM_WORLD, &status);

				i_row_fmat_FHDI_total_size.push_back(i_row_fmat_FHDI_total_recv);
			}

		}

	}

	//All processors have to know size in advance 
	if (mynode != 0) {
		i_row_fmat_FHDI_total_size.resize(totalnodes);
	}
	MPI_Bcast(&i_row_fmat_FHDI_total_size[0], totalnodes, MPI_INT, 0, MPI_COMM_WORLD);

	//if (mynode == 1) {
	//	for (int t = 0; t < i_row_fmat_FHDI_total_size.size(); t++) {
	//		cout<<"i_row_fmat_FHDI_total_size["<<t<<"]: "<< i_row_fmat_FHDI_total_size[t]<<endl;
	//	}
	//}

	//Write distributed return matrix to local hard drive concurrently
	MPI_Out_uox_mox(i_row_fmat_FHDI_total, ncol_fmat + ncol, i_row_fmat_FHDI_total_size, fh_fmat_FHDI, fmat_FHDI_temp);
	MPI_Barrier(MPI_COMM_WORLD);

	// Send distributed return matrix to master node
	if (mynode != 0) {
		MPI_Send(simp_fmat_FHDI_temp[0], i_row_fmat_FHDI_total*ncol_fmat, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}


	//All_gather distributed return matrix
	if (mynode == 0) {
		int counter8 = 0;
		for (int j = 1; j < totalnodes; j++) { // skip first element 0
			int i_row_fmat_size = i_row_fmat_FHDI_total_size[j];
			double** simp_fmat_FHDI_temp_recv = New_dMatrix(i_row_fmat_size, ncol_fmat);

			MPI_Recv(simp_fmat_FHDI_temp_recv[0], i_row_fmat_size*ncol_fmat, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);

			for (int l = 0; l < i_row_fmat_size; l++) {
				for (int m = 0; m < ncol_fmat; m++) {
					simp_fmat_FHDI[counter8][m] = simp_fmat_FHDI_temp_recv[l][m];
				}
				counter8++;
			}

			Del_dMatrix(simp_fmat_FHDI_temp_recv, i_row_fmat_size, ncol_fmat);
		}
	}

	//-------------------------------------
	//Broadcast simp_fmat_FHDI
	//-------------------------------------
	int i_row_fmat = 0;
	for (int t = 0; t < totalnodes; t++) i_row_fmat = i_row_fmat + i_row_fmat_FHDI_total_size[t];

	//cout<<"i_row_fmat is "<< i_row_fmat <<" and ncol_fmat + ncol = "<< ncol_fmat + ncol <<" at node "<<mynode<<endl;

	MPI_Bcast(simp_fmat_FHDI[0], i_row_fmat*ncol_fmat, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	//close file string to read mox
	success = MPI_File_close(&fh_mox_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of mox fail to close the file!" << endl;

	//close file string to write imputed values
	success = MPI_File_close(&fh_fmat_FHDI);
	if (success != MPI_SUCCESS) cout << "MPI I/O of imputed matrix fail to close the file!" << endl;

	//=================================
	//int success = 0;

	//MPI_File fh_binary_daty;
	//success = MPI_File_open(MPI_COMM_WORLD, "./daty_column_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_binary_daty);
	//if (success != MPI_SUCCESS) cout << "MPI I/O fail to open the file!" << endl;

	//double* array_temp = new double[ncol];;
	//int temp = 0;

	//if (mynode == 1) {
	//	MPI_In_daty_row(nrow, ncol, 2, fh_binary_daty, array_temp);// read distributed z matrix row by row
	//	cout << "array_temp: " << endl;
	//	for (int t = 0; t < ncol; t++) {
	//		cout << setw(20) << array_temp[t];
	//	}
	//	cout << endl;
	//}

	//success = MPI_File_close(&fh_binary_daty);
	//if (success != MPI_SUCCESS) cout << "MPI I/O of z matrix fail to close the file!" << endl;
	//===========================

	//Deallocation
	if (mynode != 0) {
		Del_dMatrix(fmat_FHDI_temp, i_row_fmat_FHDI_total, ncol_fmat + ncol);
		Del_dMatrix(simp_fmat_FHDI_temp, i_row_fmat_FHDI_total, ncol_fmat);
	}

	delete[] w;
	delete[] id;
	delete[] mox_temp;


	return;

}
