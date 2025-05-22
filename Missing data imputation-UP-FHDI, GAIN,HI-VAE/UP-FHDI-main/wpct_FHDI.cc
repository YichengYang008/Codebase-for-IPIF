
void wpct_initial_FHDI(List_FHDI &uox_infor, const int nrow_uox, const double* w, std::vector<double> &jp_prob)
//Description=====================================
//  calculate initial weighted probability of all uox
//  using the given weight array in w[]
//
//  written by Dr I. Cho
//  All right reserved
//
//  Algorithm: similar to "R" wpct()
//
//IN   : List_FHDI uox_infor        = Actual index list of mox in z  
//IN   : int nrow_uox               = number of uox
//IN   : double w(nrow) = weights for rows of z matrix (default = 1.0)

//OUT  : std::vector<double>      jp_prob  = weighted joint probability of the table 
//================================================
{

	//---------------
	//find new accumulated weights for each category
	//---------------
	double* d_weight = new double[nrow_uox]; // accumulated weights of all uox
	Fill_dVector(d_weight, nrow_uox, 0.0);

	std::vector<double> uox_buffer; // actual locations of each uox in z matrix

	int i_temp = 0;
	for (int i = 0; i<nrow_uox; i++) //loop for all uox
	{
		uox_buffer.clear();

		uox_infor.get_block_yicheng(i, uox_buffer);
		int uox_buffer_size = uox_buffer.size();

		for (int k = 0; k < uox_buffer_size; k++) {
			i_temp = 0; 
			i_temp = (int)(uox_buffer[k] - 1);
			d_weight[i] = d_weight[i] + w[i_temp];
		}
	}

	//-----------------
	//sum of d_weight 
	//-----------------
	double d_sum_w = 0.0;
	for (int i = 0; i<nrow_uox; i++) d_sum_w += d_weight[i];
	if (d_sum_w == 0.0) { cout << "Error! zero sum of weights in wpct from base_FHDI_MPI.cc" << endl; return; }
	//cout<<"d_sum_w is "<< d_sum_w <<endl;
	//for (int l = 0; l < nrow_uox; l++) cout<<"d_weight["<<l<<"]: "<< d_weight[l]<<endl;
	//------------------
	//prep return
	//------------------
	for (int i = 0; i<nrow_uox; i++)
	{
		jp_prob.push_back(d_weight[i] / d_sum_w);
	}

	//------------------
	//Deallocation
	//------------------
	delete[] d_weight;

}

void wpct_ultra_FHDI(std::vector<int> s_0, const int n, const int nrow_uox,
	const double* w, std::vector<double> &jp_prob, int i_loop)
	//Description=====================================
	//  calculate weighted probability of the string array 
	//  using the given weight array in w[]
	//
	//  written by Dr I. Cho
	//  All right reserved
	//
	//  Algorithm: similar to "R" wpct()
	//
	//IN   : std::vector<int> s_0 = indices of observed patterns + augumented donors 
	//IN   : const int n          = total number of s_0
	//IN   : int nrow_uox         = number of uox
	//IN   : double w[n]  	      = user-defined weight used for proportional weights

	//OUT  : std::vector<double>      jp_prob  = weighted joint probability of uox
	//================================================
{

	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//--------------------
	//Distribute nrow_uox
	//--------------------
	const int L = nrow_uox;
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

	//------------------------------------------------
	//find new accumulated weights for each category
	//----------------------------------------------
	if (mynode != 0) {
		double* d_weight_temp = new double[L_temp];// accumulated weights of distributed uox
		Fill_dVector(d_weight_temp, L_temp, 0.0);

		int counter = 0;

		for (int i = startpoint; i < endpoint; i++) {

			for (int j = 0; j < n; j++) {
				if ((i + 1) == s_0[j]) {
					d_weight_temp[counter] = d_weight_temp[counter] + w[j];  //accumulate the weight of this category
				}
			}
			counter++;
		}

		MPI_Send(d_weight_temp, L_temp, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

		//deallocation
		delete[] d_weight_temp;

	}//end of mynode

	double* d_weight = new double[nrow_uox];// accumulated weights of all uox

	if (mynode == 0) {

		double* d_weight_recv = new double[numWorkPerProc];//
		double* d_weight_recv_last = new double[numWorkLocalLast];// 
		int counter2 = 0;

		for (int j = 1; j < totalnodes; j = j + 1) {

			if (j != (totalnodes - 1)) {
				MPI_Recv(d_weight_recv, numWorkPerProc, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);

				for (int t = 0; t < numWorkPerProc; t++) {
					d_weight[counter2] = d_weight_recv[t];
					counter2++;
				}

			}

			if (j == (totalnodes - 1)) {
				MPI_Recv(d_weight_recv_last, numWorkLocalLast, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);

				for (int t = 0; t <numWorkLocalLast; t++) {
					d_weight[counter2] = d_weight_recv_last[t];
					counter2++;
				}

			}


		}//end of mynode


		 //deallocation
		delete[] d_weight_recv;
		delete[] d_weight_recv_last;

	}//end of master node

	MPI_Bcast(d_weight, nrow_uox, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//if (mynode == 1 && i_loop==2) {
	//	cout << "d_weight after Bcast at node " << mynode <<" at i_loop = "<< i_loop << endl;
	//	for (int t = 0; t < nrow_uox; t++) {
	//		cout << "d_weight[" << t << "]: " << d_weight[t] << endl;
	//	}
	//}

	//-----------------
	//sum of d_weight 
	//-----------------
	double d_sum_w = 0.0;
	for (int i = 0; i<nrow_uox; i++) d_sum_w += d_weight[i];
	if (d_sum_w == 0.0) { cout << "Error! zero sum of weights in wpct from base_FHDI_MPI.cc" << endl; return; }


	//------------------
	//prep return
	//------------------
	for (int i = 0; i<nrow_uox; i++)
	{
		jp_prob.push_back(d_weight[i] / d_sum_w);
	}

	//------------------
	//Deallocation
	//------------------
	delete[] d_weight;

}


void wpct_ultra_serial_FHDI(std::vector<int> s_0, const int n, const int nrow_uox,
	const double* w, std::vector<double> &jp_prob)
	//Description=====================================
	//  calculate weighted probability of the string array 
	//  using the given weight array in w[]
	//
	//  written by Dr I. Cho
	//  All right reserved
	//
	//  Algorithm: similar to "R" wpct()
	//
	//IN   : std::vector<int> s_0 = indices of observed patterns + augumented donors 
	//IN   : const int n          = total number of s_0
	//IN   : int nrow_uox         = number of uox
	//IN   : double w[n]  	      = user-defined weight used for proportional weights

	//OUT  : std::vector<double>      jp_prob  = weighted joint probability of uox
	//================================================
{

	//---------------
	//find new accumulated weights for each category
	//---------------
	double* d_weight = new double[nrow_uox];// accumulated weights of all uox
	Fill_dVector(d_weight, nrow_uox, 0.0);

	for (int i = 0; i < nrow_uox; i++) {
		for (int j = 0; j < n; j++) {
			if ((i + 1) == s_0[j]) {
				d_weight[i] = d_weight[i] + w[j];  //accumulate the weight of this category
			}
		}
	}

	//-----------------
	//sum of d_weight 
	//-----------------
	double d_sum_w = 0.0;
	for (int i = 0; i<nrow_uox; i++) d_sum_w += d_weight[i];
	if (d_sum_w == 0.0) { cout << "Error! zero sum of weights in wpct from base_FHDI_MPI.cc" << endl; return; }


	//------------------
	//prep return
	//------------------
	for (int i = 0; i<nrow_uox; i++)
	{
		jp_prob.push_back(d_weight[i] / d_sum_w);
	}

	//------------------
	//Deallocation
	//------------------
	delete[] d_weight;

}