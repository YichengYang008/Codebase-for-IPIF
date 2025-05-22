//#include <vector>
//#include "order_FHDI_binary.cc"
void AGMAT_Neighbor_ultra_cpp_MPI(List_FHDI &mox_infor, const int nrow_mox,
	const int nrow, const int ncol,
	List_FHDI &List_nU, std::vector<int> &rst_final)
	//Description=========================================
	// Augment missing cells in mox using the observed values of uox
	//
	// Algorithm:  All possible donors will be used to fill in the missing cell 
	//             but, if there is no matched donors in uox, this algorithm may fail
	//             as of Oct 2016
	// for each missing pattern, find all the possible donors
	// e.g., 
	// (1) a missing row   = 000
	// 	   agmat           = all observed rows
	// (2) a missing row   = a01
	//     agmat           = ac1, af1, a11, ..., az1. 
	//
	//
	// original R code: Dr. Im, J. and Dr. Kim, J. 
	// c++ code: 		Yicheng Yang and Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: August 9, 2021
	//----------------------------------------------------
	//IN	: List_FHDI mox_infor   = Actual index list of mox in z                      
	//IN    : int nrow_uox          = number of of mox
	//IN    : const int nrow        = number of rows of raw data
	//IN    : const int ncol        = number of columns of raw data
	//IN    : List_FHDI List_nU     = donors list of all mox

	//OUT   : vector<int> rst_final = fully observed patterns + augumented donors where number of rows will be determined by this code   
	//====================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	const int nr1 = nrow_mox;

	std::vector<int> loc_srst_nl; // actual locations of donors of each mox in uox
	int i_size_v_cn_z_i = 0;

	//-------------
	//LOOP for all missing rows
	//-------------
	for (int i = 0; i < nr1; i++)
	{
		loc_srst_nl.clear(); //re-initialize

		List_nU.get_block_yicheng(i, loc_srst_nl);

		//-----
		//total matching rows
		//-----
		const int i_size_loc_srst_nl = (int)loc_srst_nl.size();

		i_size_v_cn_z_i = 0;
		mox_infor.get_a_row_size(i, i_size_v_cn_z_i);
		//int i_size_v_cn_z_i = (int)mox_location.size(); //number of locations in cn having s_temp

		if (i_size_loc_srst_nl == 0) //error case
		{
			cout << "Error! there is no matched cell!" << endl; return;
		}

		//e.g.,
		//If i_size_v_cn_z_i = 3 and loc_srst_nl = {3,5}
		//rst_final = {3,5, 3,5, 3,5}

		for (int j = 0; j < i_size_v_cn_z_i; j++)
		{
			for (int k = 0; k < i_size_loc_srst_nl; k++) //repeated copy of the id number 
			{
				rst_final.push_back(loc_srst_nl[k]);
			}
		}



	}

	//TestOut
	//if (mynode == 0) {
	//	cout << "rst_final: " << endl;
	//	for (int m = 0; m < rst_final.size(); m++) {
	//		cout << "rst_final[" << m << "]: " << rst_final[m] << endl;
	//	}
	//}

	return;
}