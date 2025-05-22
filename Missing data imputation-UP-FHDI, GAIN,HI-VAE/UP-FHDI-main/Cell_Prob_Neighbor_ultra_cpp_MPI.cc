//#include "AGMAT_Neighbor_ultra_cpp_MPI.cc"
//#include "Cal_W_Neighbor_ultra_cpp_MPI.cc"

void Cell_Prob_Neighbor_ultra_cpp_MPI(const int nrow, const int ncol, List_FHDI &List_nU, List_FHDI &uox_infor, List_FHDI &mox_infor,
	int nrow_uox, int nrow_mox, double* w, std::vector<double> &jp_prob_return,
	ofstream& TestOut)

	//Description=========================================
	// make joint probability of cells with the categorized matrix z 
	// where 0 means missing data
	//
	//
	// original R code: Dr. Im, J. and Dr. Kim, J. 
	// c++ code: 		Yicheng Yang and Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: August 10, 2021
	//----------------------------------------------------

	//IN    : const int nrow = number of rows of raw data
	//IN    : const int ncol = number of columns of raw data
	//IN    : List_FHDI List_nU = donors list of all mox
	//IN    : List_FHDI uox_infor = Actual index list of uox in z
	//IN    : List_FHDI mox_infor = Actual index list of mox in z
	//IN    : int nrow_uox = number of row of uox
	//IN    : int nrow_mox = number of row of mox
	//IN    : double w(nrow) = weights for rows (default = 1.0)
	//
	//OUT   : std::vector<double>      jp_prob_return  = final updated joint probability of sorted uox 
	//====================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//-------------
	//maximum number of iterations for updating weights
	//-------------
	const int n_maximum_iteration = nrow * 100; //set by user!


	std::vector<int> ol_location;// Actual location of observed patterns in z matrix (not sorted)
	std::vector<double> ol_location_temp;// Actual location of observed patterns in z matrix (not sorted)

	uox_infor.unlist(ol_location_temp);
	int i_v_ol = 0;
	for (int k = 0; k < ol_location_temp.size(); k++) {
		i_v_ol = 0;
		i_v_ol = (int)ol_location_temp[k];
		ol_location.push_back(i_v_ol);
	}

	int i_size_ol = ol_location.size();// number of fully observed patterns in z matrix

	//std::vector<int> ml_location;

	//for (int j = 0; j < nrow_mox; j++) {
	//	for (int k = 2; k < (max_overlap_size + 2); k++) {
	//		if (mox_info_final[j][k] != 0) {
	//			ml_location.push_back(mox_info_final[j][k]);
	//		}
	//	}
	//}

	//int i_size_ml = ml_location.size();
	
	//-------------------
	//Augment observed cells for missing patterns
	//algorithm: 
	//for each missing pattern, find all the possible donors
	//e.g.,
	//mox[1] has three occurance in z matrix and two donors 
	//Actual locations of two donors in z matrix are {3,10}
	//then agmat is: {3,10, 3,10}
	//-------------------

	std::vector<int> agmat;// indices of augmented patterns (i.e., donors) of all mox 

	//Note that agmat is not sorted regarding asecending order of mox in z
	//Thus, updated weights of agmat in cal_w function should not be sorted as well to match up
	AGMAT_Neighbor_ultra_cpp_MPI(mox_infor, nrow_mox, nrow, ncol, List_nU, agmat);

	const int n_row_agmat = agmat.size(); //total number of augumented donors

	//---------
	//Translate 1. existing ox (rows with full observations) & 2. agmat
	//without appending the augmented rows onto the previous ox, i.e. the existing rows with observations
	//---------

	const int i_total_ox_agmat = i_size_ol + n_row_agmat;

	std::vector<int> s_fcd;// indices of all observed patterns in z matrix + indices of augmented patterns 
	
	//Add indices of each uox by its occurnace times
	//This is equvalent to add all observed patterns in z matrix
	//Note that updating of joint probability is based on occurnace of uox only
	//Thus, we can play the trick here
	int counts = 0;
	for (int i = 0; i < nrow_uox; i++) {
		counts = 0;
		uox_infor.get_a_row_size(i, counts);
		for (int j = 0; j < counts; j++) {
			s_fcd.push_back(i + 1);// the actual locations of all observed patterns in uox
		}
	}

	//Add indices of augmented patterns
	for (int k = 0;k < n_row_agmat; k++) {
		s_fcd.push_back(agmat[k]);
	}

	//if (mynode == 0) {
	//	cout << "s_fcd: " << endl;
	//	for (int m = 0; m < s_fcd.size(); m++) {
	//		cout << "s_fcd[" << m << "]: " << s_fcd[m] << endl;
	//	}
	//}

	if (s_fcd.size() != i_total_ox_agmat) {
		cout << "ERROR in cell probbaility estimation of ultra data!!!!" << endl;
	}

	//----------------
	//calculate weighted joint probability
	//Note: Initial jp is calculated only
	//with all the OBSERVED condensed strings in s_ocn
	//NOT with Augmented data matix 
	//-----------------

	//-------------
	//sampling weight of the observed unit
	//-------------
	double* w1 = new double[i_size_ol];
	for (int i = 0; i<i_size_ol; i++) w1[i] = w[ol_location[i] - 1]; //-1 for actual location

	//Compute initial joint probability based on fully observed patterns
	std::vector<double>		 jp_prob;// initial joint probability of all uox
	wpct_initial_FHDI(uox_infor, nrow_uox, w, jp_prob);
	const int i_size_jp_prob = (int)jp_prob.size();

	if (i_size_jp_prob != nrow_uox) {
		cout<<"ERROR in computing initial joint probability for ultra data"<<endl;
		return;
	}

	//if (mynode == 0) {
	//	cout << "jp_prob at node " << mynode << endl;
	//	for (int i = 0; i < jp_prob.size(); i++) {
	//		cout << "jp_prob[" << i << "]: " << jp_prob[i] << endl;
	//	}
	//}


	//===================================
	//===================================
	//Cal_W(): update new weights and the joint probability of cells
	//===================================
	//===================================
	std::vector<double> 		w20; //new storage for updated weights  
	std::vector<double>		 	jp_prob_0; //probability backup in the loop
	std::vector<double>		 	jp_prob_new; // updated probability

	for (int j = 0; j<i_size_jp_prob; j++)
	{
		//jp_name_new.push_back(jp_name[j]); //initialize with jp_prob
		jp_prob_new.push_back(jp_prob[j]); //initialize with jp_prob 
	}

	//MAIN ITERATION for Updating weights =======================
	for (int i_loop = 0; i_loop < n_maximum_iteration; i_loop++) {

		//------------
		//intialize with the updated joint probability
		//------------
		jp_prob_0.clear(); //re-initialize 
		for (int j = 0; j<i_size_jp_prob; j++)
		{
			jp_prob_0.push_back(jp_prob_new[j]); //initialize with jp_prob_new 
		}

		//---------------------------------------
		//update weights, w20[] 
		//Note: prob must be the newest one! i.e. jp_prob_new
		//---------------------------------------
		w20.clear();  //re-initialize 

		//Note that agmat is not sorted regarding asecending order of mox in z
		//Thus, updated weights of agmat in cal_w function should not be sorted as well to match up
		Cal_W_Neighbor_ultra_cpp_MPI(mox_infor, nrow_mox,
			nrow, ncol, jp_prob_new, w, List_nU, w20);

		//-----------
		//combine new weights
		//-----------
		double* w12 = new double[i_total_ox_agmat];
		for (int j = 0; j<i_total_ox_agmat; j++)
		{
			if (j<i_size_ol) //weights of existing w1
			{
				w12[j] = w1[j];
			}
			if (j >= i_size_ol) //updated weights 
			{
				w12[j] = w20[j - i_size_ol];
			}
		}

		//if (mynode == 0) {
		//	cout<<"w12 at iteration "<< i_loop <<endl;
		//	for (int k = 0; k < i_total_ox_agmat; k++) {
		//		cout<<"w12["<<k<<"]: "<< w12[k]<<endl;
		//	}
		//}

		jp_prob_new.clear(); //re-initialize
        //Update weights of augumented donors
		wpct_ultra_serial_FHDI(s_fcd, i_total_ox_agmat, nrow_uox, w12, jp_prob_new);

		//if (mynode == 0) {
		//	cout << "jp_prob_new at iteration " << i_loop << endl;
		//	for (int j = 0; j<i_size_jp_prob; j++)
		//	{
		//		cout << "jp_prob_new[" << j << "]: " << jp_prob_new[j] << endl;
		//	}
		//}
		//------------
		//calculate difference in the joint probability
		//------------
		double dif = 0.0;
		for (int j = 0; j < i_size_jp_prob; j++) {
			dif += (jp_prob_0[j] - jp_prob_new[j])*(jp_prob_0[j] - jp_prob_new[j]);
		}

		if (dif < 1e-6)
		{
			//cout << "Cell_Prob_Extension_cpp_MPI.. has successfully finished at iterations " << i_loop << endl;
			//TestOut<<" Cell_Prob... finished after iterations : "<< i_loop+1<<endl;
			break;
		}
		
		//------------
		//check max iterations
		//------------
		if (i_loop == n_maximum_iteration - 1)
		{
			TestOut << "CAUTION!! max iteration reached in Cell_Prob..()" << endl;
			//RPrint("CAUTION!! max iteration reached in Cell_Prob..()");
		}


		//------------
		//local deallocation
		//------------
		delete[] w12;
	}

	//---------------------------
	//prep return 
	//the latest joint probability
	//---------------------------
	jp_prob_return.clear();

	for (int j = 0; j<i_size_jp_prob; j++)
	{
		jp_prob_return.push_back(jp_prob_new[j]); //return with jp_prob_new	
	}

	//if (mynode == 0) {
	//	cout<<"jp_prob_return: "<<endl;
	//	for (int j = 0; j<i_size_jp_prob; j++)
	//	{
	//		cout<<"jp_prob_return["<<j<<"]: "<< jp_prob_return[j]<<endl;
	//	}
	//}


	delete[] w1;
	return;

}