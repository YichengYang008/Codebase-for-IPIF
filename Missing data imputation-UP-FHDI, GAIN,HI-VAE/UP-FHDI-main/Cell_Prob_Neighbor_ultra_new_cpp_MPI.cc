#include "AGMAT_Neighbor_ultra_cpp_MPI.cc"
#include "Cal_W_Neighbor_ultra_cpp_MPI.cc"
#include "wpct_FHDI.cc"

void Cell_Prob_Neighbor_ultra_new_cpp_MPI(const int nrow, const int ncol, List_FHDI &List_nU, List_FHDI &uox_infor, List_FHDI &mox_infor,
	int nrow_uox, int nrow_mox, double* w, std::vector<double> &jp_prob_return, double** mox_Imputation_group, int i_size_ml,
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
	// updated: August 9, 2021
	//----------------------------------------------------

	//IN    : const int nrow = number of rows of raw data
	//IN    : const int ncol = number of columns of raw data
	//IN    : List_FHDI List_nU = donors list of all mox
	//IN    : List_FHDI uox_infor = Actual index list of uox in z
	//IN    : List_FHDI mox_infor = Actual index list of mox in z
	//IN    : int nrow_uox = number of uox
	//IN    : int nrow_mox = number of mox
	//IN    : double w(nrow) = weights for rows (default = 1.0)
	//IN    : int i_size_ml  = total number of missing units in z matrix
	//
	//OUT   : std::vector<double>      jp_prob_return  = final updated joint probability of sorted uox 
	//OUT   : double** mox_Imputation_group            = (gloabl index of all missing units in z matrix) + (corresponding uox id with the highest conditional probability)
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

	if (i_size_ol + i_size_ml != nrow) cout<<"ERROR!!!! i_size_ol + i_size_ml != nrow in cell prob new ultra!!!"<<endl;


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
	//mox[1] has two occurance in z matrix and two donors 
	//Actual locations of two donors in z matrix are {3,10}
	//then agmat is: {3,10, 3,10}
	//-------------------

	std::vector<int> agmat;// indices of augmented patterns (i.e., donors) of all mox 

	//Note that agmat is not sorted regarding asecending order of mox in z
	//Thus, updated weights of agmat in cal_w function should not be sorted as well to line up
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


	//---------------------------------------------------------------------
	//Information of conditional probability of donors of all missing units
	//This information is only required in variance estimation using Linerization
	//---------------------------------------------------------------------

	List_FHDI List_conditional_prob(i_size_ml);
	//e.g.,
	//mox[0] has two occurance in z matrix {11,21} and two donors 
	//Actual locations of two donors in z matrix are {3,10}
	//then List_conditional_prob is:
	//1     11      {0.333,0.6666}
	//1     21      {0.333,0.6666}

	List_FHDI agmat_conditional_prob(i_size_ml); //(donor location in uox) of all missing units
	//then agmat_conditional_prob is:
	//              {3, 10}
	//              {3, 10}

	int counter3 = 0;
	int counter4 = 0;
	int mox_size_temp = 0;
	int donor_size_temp = 0;
	std::vector<double>		 agmat_temp;//donor buffer

	for (int t = 0; t < nrow_mox; t++) {
		mox_size_temp = 0;
		mox_infor.get_a_row_size(t, mox_size_temp);
		//mox_size_temp = mox_info_final[t][1];//occurance

		donor_size_temp = 0;
		List_nU.get_a_row_size(t, donor_size_temp);

		for (int k = 0; k < mox_size_temp; k++) {
			agmat_temp.clear();
			for (int m = 0; m < donor_size_temp; m++) {
				agmat_temp.push_back((double)agmat[counter4]);
				counter4++;
			}

			agmat_conditional_prob.put_block_yicheng(counter3, donor_size_temp, agmat_temp);
			counter3++;
		}

	}



	//----------------------------------------
	//Initialize first two columns of List_conditional_prob
	//----------------------------------------
	int counter = 0; 
	std::vector<double> mox_buffer; // 
	std::vector<double> conditional_buffer; //
	int mox_size_temp2 = 0; 

	for (int t = 0; t < nrow_mox; t++) {
		mox_size_temp2 = 0;
		mox_buffer.clear();
		mox_infor.get_a_row_size(t, mox_size_temp2);
		mox_infor.get_block_yicheng(t, mox_buffer);

		if (mox_buffer.size() != mox_size_temp2) { cout<<"ERROR! mox_buffer is incorrect in cell probability!"<<endl; }

		for (int l = 0; l < mox_size_temp2; l++) {
			conditional_buffer.clear();
			conditional_buffer.push_back(t + 1);
			conditional_buffer.push_back( (int)mox_buffer[l] );

			List_conditional_prob.put_block_yicheng(counter, 2, conditional_buffer);

			counter++;
		}
	}

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

		double cell_prob_each2 = MPI_Wtime();

		jp_prob_new.clear(); //re-initialize
        //Update weights of augumented donors
		wpct_ultra_FHDI(s_fcd, i_total_ox_agmat, nrow_uox, w12, jp_prob_new, i_loop);

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

		//if (mynode == 1) {
		//	cout << "wpct_ultra_FHDI time is " << MPI_Wtime() - cell_prob_each2 << " at i_loop = " << i_loop << endl;
		//}

		if (dif < 1e-6)
		{

			//------------------------------------------
			//Store corresponding conditional probability
			//of all missing units
			//-------------------------------------------
			int counter1 = 0;
			int counter2 = 0;
			int mox_size_temp = 0;
			int donor_size_temp = 0;
			std::vector<double>		 conditional_temp;//donor buffer

			for (int t = 0; t < nrow_mox; t++) {
				mox_size_temp = 0;
				mox_infor.get_a_row_size(t, mox_size_temp);
				//mox_size_temp = mox_info_final[t][1];

				donor_size_temp = 0;
				List_nU.get_a_row_size(t, donor_size_temp);

				for (int k = 0; k < mox_size_temp; k++) {
					conditional_temp.clear();
					for (int m = 0; m < donor_size_temp; m++) {
						conditional_temp.push_back(w20[counter2]);
						counter2++;
					}
					List_conditional_prob.put_block_yicheng(counter1, donor_size_temp, conditional_temp);

					counter1++;
				}

			}

			if (counter2 != w20.size()) cout<<"ERROE!!! counter2 != w20.size() in cell prob new ultra!!!"<<endl;


			//if(mynode==1) cout << "Cell_Prob_Extension_cpp_MPI.. has successfully finished at iterations " << i_loop << endl;
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

	//if (mynode == 1) {

	//	cout << "List_conditional_prob at node " << mynode << endl;
	//	List_conditional_prob.print_List_FHDI_yicheng();


	//	cout << "agmat_conditional_prob at node " << mynode << endl;
	//	agmat_conditional_prob.print_List_FHDI_yicheng();
	//}

	//--------------------------------------------------------------
	//Determination of uox with the highest conditional probability
	//for all missing units in z matrix
	//---------------------------------------------------------------

	double max_temp = 0.0;
	int index_temp = 0;
	std::vector<double>		 agmat_buffer;//donor buffer
	std::vector<double>		 probability_buffer;//conditional probability buffer

	for (int t = 0; t < i_size_ml; t++) {
		probability_buffer.clear();
		List_conditional_prob.get_block_yicheng(t, probability_buffer);
		mox_Imputation_group[t][0] = probability_buffer[1];// Actual locations of all missing units

		max_temp = 0.0;
		index_temp = 0;

		//Get location of uox who has the highest cp
		int i_size_temp = probability_buffer.size();
		for (int k = 2; k < i_size_temp; k++) {
			if (probability_buffer[k] > max_temp) {
				max_temp = probability_buffer[k];
				index_temp = k - 2;//to match agmat_conditional_prob
			}
		}

		agmat_buffer.clear();
		agmat_conditional_prob.get_block_yicheng(t, agmat_buffer);
		mox_Imputation_group[t][1] = agmat_buffer[index_temp];
	}



	//if (mynode == 1) {
	//	cout<<"jp_prob_return: "<<endl;
	//	for (int j = 0; j<i_size_jp_prob; j++)
	//	{
	//		cout<<"jp_prob_return["<<j<<"]: "<< jp_prob_return[j]<<endl;
	//	}

	//	cout << "mox_Imputation_group: " << endl;
	//	for (int j = 0; j < i_size_ml; j++) {
	//		for (int k = 0; k < 2; k++) {
	//			cout << setw(20) << mox_Imputation_group[j][k];
	//		}
	//		cout << endl;
	//	}
	//}

	//Deallocation
	delete[] w1;

	return;

}