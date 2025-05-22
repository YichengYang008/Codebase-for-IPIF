//#include <vector>
//#include "order_FHDI_Yicheng.cc"
void Cal_W_Neighbor_ultra_new_cpp_MPI(int** uox_info_final, const int nrow_uox,
	int** mox_info_final, const int nrow_mox,
	const int nrow, const int ncol, const int max_overlap_size, std::vector<double> jp_prob, double* w,
	List_FHDI &List_nU, List_FHDI &List_conditional_prob, std::vector<double> &v_rst_final)
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
	//IN    : int** uox_info_final = (gloabl index of uox in z matrix) + occurnace + (actual locations of an unique pattern in categorized matrix)                          
	//IN    : int nrow_uox = number of row of uox
	//IN    : int** mox_info_final = (gloabl index of mox in z matrix) + occurnace + (actual locations of an unique pattern in categorized matrix)
	//IN    : int nrow_mox = number of row of mox
	//IN    : const int nrow = number of rows of raw data
	//IN    : const int ncol = number of columns of raw data
	//IN    : const int max_overlap_size = maximum number of actual locations of an unique pattern in categorized matrix
	//IN    : std::vector<double> jp_prob_return = former joint probability of uox
	//IN    : double w(nrow) = weights for rows (default = 1.0)
	//IN    : List_FHDI List_nU = donors list of all mox

	//OUT   : std::vector<double> &v_rst_final  = new weights 
	//====================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	const int nr1 = nrow_mox;
	const int nr2 = nrow_uox;


	std::vector<int> mox_location; // actual location of mox in z matrix
	std::vector<int> loc_srst_ncol1;// actual locations of donors of each mox in uox
	int i_size_zid = 0;

	std::vector<double> w_srst_ncol; // weights for mox
	std::vector<double> jp_zi;// updated joint probability

	double w_srst_ncol_conditional; // weights for mox
	std::vector<double> jp_zi_conditional;// updated joint probability
	std::vector<double> cp_temp;// updated joint probability

	//----------------
	//LOOP for all missing rows
	//----------------
	for (int i = 0; i < nr1; i++) {
		mox_location.clear(); //re-initialize
		loc_srst_ncol1.clear(); //re-initialize //Note: this is different from loc_srst_ncol

		for (int k = 2; k < (max_overlap_size + 2); k++) {
			if (mox_info_final[i][k] != 0) {
				mox_location.push_back(mox_info_final[i][k]);
			}
		}

		int i_size_v_cn_z_i = (int)mox_location.size(); //occurance of this mox in z

		List_nU.get_block_yicheng(i, loc_srst_ncol1);
		//-----
		//total number of the observed rows that matches the current missing row
		//-----
		const int i_size_loc_srst_ncol1 = (int)loc_srst_ncol1.size();
		if (i_size_loc_srst_ncol1 == 0) //error case
		{
			cout << "Error! there is no matched cell!" << endl; return;
		}

		//-------------------
		//get ready w[] at srst = ncol
		//-------------------

		w_srst_ncol.clear();

		for (int j = 0; j<i_size_v_cn_z_i; j++)  //repeat each entity 
		{
			for (int k = 0; k<i_size_loc_srst_ncol1; k++)
			{
				//------
				//Note: the weights are pulled out from loc_srst_loc NOT .._loc1 
				//      below "loc_srst_ncol" means the target rows to be imputed later
				//------
				w_srst_ncol.push_back(w[mox_location[j] - 1]); //-1 for actual location
			}
		}

		w_srst_ncol_conditional = w[mox_location[0] - 1]; //-1 for actual location

		//-------------------
		//get ready second column
		//below "loc_srst_ncol1" contains the obs. row numbers that will serve as donor
		//-------------------
		double sum_jp_loc = 0.0; //sum of jp only at location where srst = ncol
		for (int j = 0; j < i_size_loc_srst_ncol1; j++) {
			sum_jp_loc += jp_prob[loc_srst_ncol1[j] - 1];
		}

		jp_zi.clear();
		jp_zi_conditional.clear();

		for (int j = 0; j < i_size_v_cn_z_i; j++)  //repeat entire jp..[] by z_i_now times 
		{
			for (int k = 0; k < i_size_loc_srst_ncol1; k++) {
				jp_zi.push_back(jp_prob[loc_srst_ncol1[k] - 1] / sum_jp_loc);
			}
		}

		for (int k = 0; k < i_size_loc_srst_ncol1; k++) {
			jp_zi_conditional.push_back(jp_prob[loc_srst_ncol1[k] - 1] / sum_jp_loc);
		}


		i_size_zid = jp_zi.size();
		if (w_srst_ncol.size() - i_size_zid != 0) {
			cout<<"ERROR in cal_w_neighbor for ultra "<<endl;
			return;
		}

		for (int j = 0;j < i_size_zid; j++) {
			v_rst_final.push_back(jp_zi[j] * w_srst_ncol[j]);
		}

		for (int j = 0; j < i_size_loc_srst_ncol1; j++) {
			cp_temp.push_back(jp_zi_conditional[j] * w_srst_ncol_conditional);
		}

		List_conditional_prob.put_block_yicheng_replace(i, i_size_loc_srst_ncol1, cp_temp);
	}


	return;
}