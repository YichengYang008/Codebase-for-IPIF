//#include "Cell_Prob_Extension_cpp.cc"

void Rep_CellP_Neighbor_ultra(const int nrow, const int ncol,
	List_FHDI &uox_infor, List_FHDI &mox_infor, int nrow_uox, int nrow_mox,
	RepWeight_FHDI &d_rw, List_FHDI &List_nU, std::vector<int> &List_rst_prob_size,
	ofstream& TestOut)
//Description============================================
// compute joint cell probability using replicate weight rw
// 
// R code: Dr. Im, J., and Dr. Kim, J. 
// C++   : Yicheng Yang and Dr. Cho, I.
// All rights reserved
// Last update: August 10, 2021
//

//IN    : const int nrow = number of rows of raw data
//IN    : const int ncol = number of columns of raw data
//IN    : List_FHDI uox_infor = Actual index list of uox in z
//IN    : List_FHDI mox_infor = Actual index list of mox in z
//IN    : int nrow_uox = number of row of uox
//IN    : int nrow_mox = number of row of mox
//IN    : RepWeight_FHDI d_rw = replicate weights
//IN    : List_FHDI List_nU = donors list of all mox
//
//OUT  : std::vector<int> List_rst_prob_size = size of joint cell probability computed using different replicate weights

//======================================================== 
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//indices of unique patterns who appear firstly: mox infor + uox infor
	//Note that it is not sorted
	std::vector<int> s_ncx;

	//Add mox indices
	std::vector<double> mox_buffer; // index of all patterns in z of distributed uox
	for (int k = 0; k < nrow_mox; k++) {
		mox_buffer.clear();
		mox_infor.get_block_yicheng(k, mox_buffer);
		s_ncx.push_back((int)mox_buffer[0]);
	}

	//Add uox indices
	std::vector<double> uox_buffer; // index of all patterns in z of distributed uox
	for (int t = 0; t < nrow_uox; t++) {
		uox_buffer.clear();
		uox_infor.get_block_yicheng(t, uox_buffer);
		s_ncx.push_back((int)uox_buffer[0]);
	}


	//Now, i_nc means the total number of unique patterns
	const int i_nc = s_ncx.size();

	//Initialization
	std::vector<double> jp_prob_return; 
	std::vector<int> List_rst_prob_size_send;

	double* w_UserDefined = new double[nrow]; 

	//if (mynode == 0) cout << "The value of i_nc in RepCellP is " << i_nc << endl;
	
	if (i_nc < totalnodes) {
		cout << "i_nc is smaller than totalnodes in Rep_cellP, ERROR!!!" << endl;
		return;
	}

	//------------------------------------------------------ Job assignment
	int number_before_split = 0;
	int number_after_split = 0;
	int startpoint = 0;
	int endpoint = 0;
	int splitnode = 0;

	if (i_nc % (totalnodes - 1) != 0) { splitnode = i_nc % (totalnodes - 1); }
	if (i_nc % (totalnodes - 1) == 0) { splitnode = 1; }
	//if (mynode == 0) cout << "Split: " << splitnode << endl;
	number_after_split = floor(1.0*i_nc / (1.0*totalnodes - 1));
	number_before_split = 1.0*(i_nc - floor(1.0*i_nc / (1.0*totalnodes - 1)) *(totalnodes - splitnode - 1)) / splitnode;

	if (mynode >= 1 && mynode <= splitnode) {
		startpoint = (mynode - 1) * number_before_split;
		endpoint = mynode *number_before_split;
	}
	if (mynode > splitnode) {
		startpoint = splitnode * number_before_split + (mynode - splitnode - 1)*number_after_split;
		endpoint = splitnode * number_before_split + (mynode - splitnode)*number_after_split;
	}
	if ((number_before_split*splitnode + number_after_split* (totalnodes - splitnode - 1)) != i_nc) {
		TestOut << "Work Assignment Error!!!!" << endl;
		return;
	}

	int i_loc = 0;
	for (int i = startpoint; i<endpoint; i++)
	{
   
		i_loc = s_ncx[i] - 1;

		jp_prob_return.clear(); 

		for (int j = 0; j<nrow; j++) w_UserDefined[j] = d_rw(j, i_loc);
		
		//cout<< "w_UserDefined at i = " <<i<<" at node "<<mynode<<endl;
		//for (int p = 0; p < nrow; p++) {
		//	cout<< setw(20) << w_UserDefined[p];
		//}
		//cout << endl;

		//-----------------------------
		//calculate joint probability 
		//using the Jackknife replicate weights
		//------------------------------

		Cell_Prob_Neighbor_ultra_cpp_MPI(nrow, ncol, List_nU, uox_infor, mox_infor,

		nrow_uox, nrow_mox, w_UserDefined, jp_prob_return,

		TestOut);		
								
		//-------------
		//prep return
		//------------
		List_rst_prob_size_send.push_back(jp_prob_return.size());
	}
	
	if (mynode != 0) {
		MPI_Send(&List_rst_prob_size_send[0], List_rst_prob_size_send.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);
	}

	std::vector<int> List_rst_prob_size_recv_before(number_before_split);
	std::vector<int> List_rst_prob_size_recv_after(number_after_split);

	//All_gather
	if (mynode == 0) {

		for (int j = 1; j < totalnodes; j = j + 1) {

			if (j >= 1 && j <= splitnode) {
				MPI_Recv(&List_rst_prob_size_recv_before[0], number_before_split, MPI_INT, j, 1, MPI_COMM_WORLD, &status);

				for (int k = 0; k < number_before_split;k++) {
					List_rst_prob_size.push_back(List_rst_prob_size_recv_before[k]);
				}
			}

			if (j > splitnode) {
				MPI_Recv(&List_rst_prob_size_recv_after[0], number_after_split, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				
				for (int k = 0; k < number_after_split;k++) {
					List_rst_prob_size.push_back(List_rst_prob_size_recv_after[k]);
				}

			}//end of j > splitnode
		}//end of j
	}//end of mynode=0

	//resize it before broadcast
	if (mynode != 0) { 
		List_rst_prob_size.resize(i_nc); 
	}

	//Broadcast it to all slave processors
	MPI_Bcast(&List_rst_prob_size[0], i_nc, MPI_INT, 0, MPI_COMM_WORLD);


	//Deallocation
	delete[] w_UserDefined;
	
	return;				   
}
