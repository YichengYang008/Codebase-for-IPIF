
void Ranking_m(const int nrow, const int ncol, double** x_raw, int** r_raw, double ** correlation_yicheng, int** correlation_ranking, ofstream& TestOut)

//Description=========================================
//Make correlation matrix correlation_yicheng
//Make ranking matrix of each variable according to correlation matrix in descending order

//IN	: double x(nrow, ncol) 	= {y1, y2, ... } total data containing missing values
//IN	: double r(nrow, ncol) 	= {y1, y2, ... } total response inndicators containing missing values
//OUT   : double correlation_yicheng(ncol, ncol);
//OUT   : int correlation_ranking(ncol, ncol-1); // Ranking of correlation of each variable in descending order. 
// Note it excludes itself from ranking
//=====================================================
{
	//----------------
	//Prepare fully observed y matrix
	//---------------------
	std::vector<int> ol;
	int d_temp = 0;
	for (int i_row = 0; i_row < nrow; i_row++)
	{
		d_temp = 1.0;
		for (int i_col = 0; i_col < ncol; i_col++)
		{
			if (r_raw[i_row][i_col] == 0) { d_temp = 0.0; break; } //found zero, i.e. missing cell
		}

		if (fabs(d_temp) > 1e-15) //this row has no missing cells
		{
			ol.push_back(i_row);
		} //actual number of the row having no missing cells
	}
	int nrow_ol = ol.size();

	double** ol_matrix = New_dMatrix(nrow_ol, ncol);

	for (int i = 0;i < nrow_ol;i++) {
		for (int j = 0; j < ncol; j++) {
			ol_matrix[i][j] = x_raw[ol[i]][j];
		}
	}
	//TestOut << "ol_matrix[]" << endl;
	//for (int i = 0; i < nrow_ol; i++)
	//{
	//	for (int j = 0; j < ncol; j++) { TestOut << setw(20) << ol_matrix[i][j]; }
	//	TestOut << endl;
	//}

	//-----------------------
	//Compute corrrlation matrix
	//-----------------------

	//double** correlation_yicheng = New_dMatrix(ncol, ncol);
	correlation_FHDI(ol_matrix, nrow_ol, ncol, correlation_yicheng);
	TestOut << "correlation matrix: " << endl;
	for (int i = 0; i < ncol; i++) {
		for (int j = 0; j < ncol; j++) {
			TestOut << setw(20) << correlation_yicheng[i][j];
		}
		TestOut << endl;
	}

	//-----------------------
	//Compute ranking matrix 
	//----------------------

	int** correlation_m_temp = New_iMatrix(ncol, (ncol - 1));

	for (int i = 0; i < ncol; i++) {

		double* d_source_temp = new double[ncol];
		int* i_return = new int[ncol]; //order of score actual loc

		for (int j = 0; j < ncol; j++) {
			d_source_temp[j] = abs(correlation_yicheng[i][j]);// Note the ranking of correlation matrix should be based on absolute value !!!
		}
		order_FHDI(d_source_temp, ncol, i_return);

		// Note i_return_temp must exclude itself priorly in case that i_return have several correlations of 1s
		std::vector<int> i_return_temp;
		for (int j1 = 0; j1 < ncol; j1++) {
			if (i_return[j1] != (i + 1)) {
				i_return_temp.push_back(i_return[j1]);
			}
		}

		if (i_return_temp.size() != (ncol - 1)) TestOut << "Error! The ranking of simple correlation is not correct!" << endl;

		for (int k3 = 0; k3 < (ncol - 1); k3++) {

			correlation_m_temp[i][k3] = i_return_temp[k3];

		}

		delete[] d_source_temp;
		delete[] i_return;
	}

	TestOut << "correlation_m_temp in ascending order" << endl;
	for (int kk2 = 0; kk2 < ncol; kk2++) {
		for (int kk3 = 0; kk3 < (ncol - 1); kk3++) {
			TestOut << setw(20) << correlation_m_temp[kk2][kk3];
		}
		TestOut << endl;
	}

	// Reverse the ranking matrix in the descending order
	for (int i = 0; i < ncol; i++) {
		for (int j = 0; j < (ncol - 1); j++) {
			correlation_ranking[i][j] = correlation_m_temp[i][ncol - 2 - j];
		}
	}

	TestOut << "correlation_ranking in descending order" << endl;
	for (int kk2 = 0; kk2 < ncol; kk2++) {
		for (int kk3 = 0; kk3 < (ncol - 1); kk3++) {
			TestOut << setw(20) << correlation_ranking[kk2][kk3];
		}
		TestOut << endl;
	}

	Del_dMatrix(ol_matrix, nrow_ol, ncol);
	//Del_dMatrix(correlation_yicheng, ncol, ncol);
	Del_iMatrix(correlation_m_temp, ncol, (ncol - 1));

	return;
}



void Ranking_top(const int nrow_ol, const int ncol, const int top, double** ol_matrix, int** correlation_ranking_top, ofstream& TestOut)

//Description=========================================
//Make correlation ranking matrix correlation_ranking
//Make ranking matrix of each variable according to correlation matrix in descending order

//IN	: double x(nrow, ncol) 	= {y1, y2, ... } total data containing missing values
//IN	: double r(nrow, ncol) 	= {y1, y2, ... } total response inndicators containing missing values
//OUT   : double correlation_yicheng(ncol, ncol);
//OUT   : int correlation_ranking(ncol, ncol-1); // Ranking of correlation of each variable in descending order. 
// Note it excludes itself from ranking
//=====================================================
{

	//-----------------------
	//Compute corrrlation matrix
	//-----------------------

	double* x1 = new double[nrow_ol];
	double* x2 = new double[nrow_ol];
	double d_sum = 0.0;
	std::vector<double> cov;

	for (int j = 0; j<ncol; j++) //from the first column to the second last column
	{
		cov.clear();

		for (int j_next = 0; j_next<ncol; j_next++) //next column including itself 
		{
			for (int i = 0; i<nrow_ol; i++)
			{
				x1[i] = ol_matrix[i][j];   //jth column
				x2[i] = ol_matrix[i][j_next];//next column
			}

			//---
			//each column's mean
			//---
			double x1_mean = 0.0; double x2_mean = 0.0;
			for (int i = 0; i<nrow_ol; i++)
			{
				x1_mean += x1[i];   //jth column
				x2_mean += x2[i];   //next column
			}
			x1_mean = x1_mean / nrow_ol;
			x2_mean = x2_mean / nrow_ol;
			//TestOut<<"x1_mean is "<< x1_mean <<", and x2_mean is "<< x2_mean <<" at j_next = "<< j_next <<endl;
			//-----
			//calculate covariance of two columns
			//-----
			d_sum = 0.0;
			for (int i_1 = 0; i_1<nrow_ol; i_1++)
			{
				d_sum += (x1[i_1] - x1_mean)*(x2[i_1] - x2_mean);
			}

			//----------------
			//calculate variance of two columns
			//----------------
			double x1_var = 0.0; double x2_var = 0.0;
			double var_sum = 0.0;
			for (int i_2 = 0; i_2 < nrow_ol;i_2++) {
				var_sum = var_sum + (x1[i_2] - x1_mean)*(x1[i_2] - x1_mean);
			}
			x1_var = var_sum;

			var_sum = 0.0;
			for (int i_3 = 0; i_3 < nrow_ol;i_3++) {
				var_sum = var_sum + (x2[i_3] - x2_mean)*(x2[i_3] - x2_mean);
			}
			x2_var = var_sum;

			//TestOut << "d_sum is " << d_sum  <<" and x1_var is "<< x1_var <<", x2_var is "<< x2_var << " at j_next = " << j_next << endl;
			//---------
			//store covariance using symmetry property
			//---------
			cov.push_back(d_sum / sqrt(x1_var* x2_var));

			//cov[j][j_next] = d_sum / sqrt(x1_var* x2_var);
			//cov[j_next][j] = d_sum / sqrt(x1_var* x2_var);
		}


		std::vector<int> i_return;
		double* d_source_temp = new double[ncol];

		for (int k = 0; k < ncol; k++) {
			d_source_temp[k] = abs(cov[k]);// Note the ranking of correlation matrix should be based on absolute value !!!
										   //TestOut<<"cov["<<k<<"] : "<< abs(cov[k]) <<" at j = "<<j<<endl;
		}

		order_FHDI(d_source_temp, ncol, i_return);

		if (i_return.size() != ncol) cout<<"ERROR!!! Something wrong in the Ranking_m function. "<<endl;

		for (int t = 0; t < top; t++) {
			correlation_ranking_top[j][t] = i_return[ncol - 2 - t]; //exclude the rank of itself
		}


		delete[] d_source_temp;
	}


	//---------
	//Deallocation
	//---------
	delete[] x1;
	delete[] x2;

	//TestOut << "correlation_ranking_top in descending order with top = " <<top<< endl;
	//for (int kk2 = 0; kk2 < ncol; kk2++) {
	//	for (int kk3 = 0; kk3 < top; kk3++) {
	//		TestOut << setw(20) << correlation_ranking_top[kk2][kk3];
	//	}
	//	TestOut << endl;
	//}

	return;

}






void Ranking_top_ultra(const int nrow, const int ncol, const int memory, const int top, std::vector<int> v_ol, double** correlation_top,
	int** correlation_ranking_top, ofstream& TestOut)

//Description=========================================
//Make correlation ranking matrix correlation_ranking
//Make ranking matrix of each variable according to correlation matrix in descending order

//IN	: double x(nrow, ncol) 	= {y1, y2, ... } total data containing missing values
//IN	: double r(nrow, ncol) 	= {y1, y2, ... } total response inndicators containing missing values
//OUT   : double correlation_yicheng(ncol, ncol);

//OUT   : int correlation_top(ncol, ncol-1); // Top correlation of each variable in descending order. 
//OUT   : int correlation_ranking_top(ncol, ncol-1); // Top ranking of correlation of each variable in descending order. 
//        Note it excludes itself from ranking
//=====================================================
{
	//-- MPI variables
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//Note that original top = min(ncol-1, 100) which excludes itself
	//i_top = 1 + min(top-1, 100) which includes itself and the correlation and rank of 
	//itself will be excluded afterwards
	int i_top = 0;
	i_top = top + 1;
	
	//cout<<"i_top in Ranking_top_ultra is "<< i_top <<endl;


	//--------------------------------------------------
	//Distribute all variables to all slave processors
	//---------------------------------------------------
	const int L = ncol;  
	int numWorkPerProc = (int)floor(1.0*L / (1.0*totalnodes - 1));
	int numWorkLocalLast = L - numWorkPerProc * (totalnodes - 2);

	int L_temp = 0;
	if ((mynode != 0)&&(mynode != (totalnodes - 1))) L_temp = numWorkPerProc;
	if (mynode == (totalnodes - 1)) L_temp = numWorkLocalLast;

	//cout << "L_temp is " << L_temp << " at node " << mynode << endl;

	int startpoint = 0;
	int endpoint = 0;

	if (mynode != 0 && mynode != (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = mynode*numWorkPerProc;

		//RPrint("(Variance)Strating point and ending point on node ");RPrint(mynode);
		//RPrint("are: \n");
		//RPrint(startpoint); RPrint(endpoint);
		if (endpoint - startpoint != L_temp) {
			TestOut << "category boundary ERROR!!!" << endl;
			return;
		}
	}

	if (mynode == (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = (mynode - 1)*numWorkPerProc + numWorkLocalLast;

		//RPrint("(Variance)Strating point and ending point on node ");RPrint(mynode);
		//RPrint("are: \n");
		//RPrint(startpoint); RPrint(endpoint);
		if (endpoint - startpoint != L_temp) {
			TestOut << "category boundary ERROR!!!" << endl;
			return;
		}
	}
	//cout << "Ranking_m startpoint: " << startpoint << "; endpoint: " << endpoint << " at node " << mynode << endl;
	int nrow_ol = v_ol.size();
	double** ol_matrix = New_dMatrix(nrow_ol, L_temp);
	double* array_temp = new double[nrow];

	//-----------------------------------------------------
	//Read distributed variables on each slave processors
	//------------------------------------------------------
	MPI_File fh_binary_daty;
	int success = 0;

	success = MPI_File_open(MPI_COMM_WORLD, "./daty_column_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_binary_daty);
	if (success != MPI_SUCCESS) cout << "MPI I/O fail to open the file!" << endl;

	int counter = 0;
	for (int k = startpoint; k < endpoint; k++) {

		MPI_In_raw(nrow, k, fh_binary_daty, array_temp);

		for (int m = 0; m < nrow_ol; m++) {
			ol_matrix[m][counter] = array_temp[v_ol[m] - 1];
		}
		counter++;
	}

	delete[] array_temp;
	//if (mynode == 1) {
	//	cout << " ol matrix read test at node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < nrow_ol; kk2++) {
	//		for (int kk3 = 0; kk3 < L_temp; kk3++) {
	//			cout << setw(20) << ol_matrix[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//}

	//Initilization
	double** correlation_top_temp = NULL;//top correlation on each slave processor
	int** correlation_ranking_top_temp = NULL;//top correlation rank on each slave processor

	if (mynode != 0) {
		correlation_top_temp = New_dMatrix(L_temp, i_top);
		correlation_ranking_top_temp = New_iMatrix(L_temp, i_top);

		Fill_dMatrix(correlation_top_temp, L_temp, i_top, 0.0);
		Fill_iMatrix(correlation_ranking_top_temp, L_temp, i_top, 0);
	}

	//Initilization
	int recursion_size = 0;
	recursion_size = (int)floor((memory * 0.35 * 1e+9) / (nrow * 8));
	if (ncol < recursion_size) recursion_size = ncol;

	int recursion_time = 0; // total number of times to import partial uox matrix

	if (ncol % recursion_size == 0) recursion_time = ncol/recursion_size;
	if (ncol % recursion_size != 0) recursion_time = (int)floor(1.0*ncol / recursion_size) + 1;

	int recursion_size_last = ncol - (recursion_time - 1)*recursion_size;
	if (recursion_size_last < 10) cout<<"Causion!!! recursion_size_last in ranking_m is less than 10!!!"<<endl;

	if (mynode == 0) recursion_time = 0;
	//if (mynode == 0 || mynode == 1) cout << "Ranking_m recursion_time is " << recursion_time << " and recursion_size_last is " << recursion_size_last << " and recursion_size is " << recursion_size << " at node " << mynode << endl;

	int recursion_size_temp = 0;

	//if(mynode==0) cout << "recursion_time is " << recursion_time << " and recursion_size_last is " << recursion_size_last << endl;

	double** ol_recursion = NULL;//resursive read variables within recursion, note that it should be from fully observed rows
	double* array_temp2 = NULL;//buffer to hold one column of daty

	double** correlation_recur = NULL;//temporary tank of (top correlation in former iteration) + (new generated correlation in current iteration)
	int** correlation_ranking_recur = NULL;//temporary tank of (top correlation rank in former iteration) + (new generated correlation rank in current iteration)

	for (int i_recur = 0; i_recur < recursion_time; i_recur++) {
		if (i_recur != (recursion_time - 1)) recursion_size_temp = recursion_size;
		if (i_recur == (recursion_time - 1)) recursion_size_temp = recursion_size_last;

		//boundary index of recursively read z matrix
		startpoint = 0;
		endpoint = 0;

		if (mynode != 0) {

			ol_recursion = New_dMatrix(nrow_ol, recursion_size_temp);
			array_temp2 = new double[nrow];

			if (i_recur != (recursion_time - 1)) {
				startpoint = i_recur*recursion_size;
				endpoint = (i_recur + 1)*recursion_size;
			}

			if (i_recur == (recursion_time - 1)) {
				startpoint = i_recur*recursion_size;
				endpoint = i_recur*recursion_size + recursion_size_last;
			}
		}

		//if(mynode==0) cout << "recusive startpoint: " << startpoint << "; recusive endpoint: " << endpoint << " at recursion " << i_recur <<" at node "<<mynode<<endl;
		//if(mynode==0) cout<<"L_temp is "<< L_temp<<" and recursion_size_temp is " << recursion_size_temp <<" at i_recur "<< i_recur <<" at node "<<mynode<<endl;
		
		//read recursive read variables
		counter = 0;
		for (int k = startpoint; k < endpoint; k++) {

			MPI_In_raw(nrow, k, fh_binary_daty, array_temp2);

			for (int m = 0; m < nrow_ol; m++){
				ol_recursion[m][counter] = array_temp2[v_ol[m] - 1];//Note that we only need fully observed rows
			}
			counter++;
		}

		//if (mynode == 1) {
		//	cout << " ol recusion read test at node " << mynode <<" at i_recur = "<< i_recur <<endl;
		//	for (int kk2 = 0; kk2 < nrow_ol; kk2++) {
		//		for (int kk3 = 0; kk3 < recursion_size_temp; kk3++) {
		//			cout << setw(20) << ol_recursion[kk2][kk3];
		//		}
		//		cout << endl;
		//	}
		//}

		//if ( (mynode != 0) && (i_recur == 0)) {
		//	correlation_recur = New_dMatrix(L_temp, recursion_size_temp);
		//	correlation_ranking_recur = New_iMatrix(L_temp, recursion_size_temp);
		//}


		//If i_recur = 0, it is a temporary tank of new generated correlation in first iteration
		//If i_recur = 0, it is a temporary tank of new generated correlation rank in first iteration
		if ((mynode != 0) && (i_recur == 0)) {
			correlation_recur = New_dMatrix(L_temp, recursion_size_temp);
			correlation_ranking_recur = New_iMatrix(L_temp, recursion_size_temp);
			
			Fill_dMatrix(correlation_recur, L_temp, recursion_size_temp, 0.0);
			Fill_iMatrix(correlation_ranking_recur, L_temp, recursion_size_temp, 0);
		}

		//If i_recur != 0, it is a temporary tank of (top correlation in former iteration) + (new generated correlation in current iteration)
		//If i_recur != 0, it is a temporary tank of (top correlation rank in former iteration) + (new generated correlation rank in current iteration)
		if ( (mynode != 0) && (i_recur != 0)) {
			correlation_recur = New_dMatrix(L_temp, recursion_size_temp + i_top);
			correlation_ranking_recur = New_iMatrix(L_temp, recursion_size_temp + i_top);

			Fill_dMatrix(correlation_recur, L_temp, recursion_size_temp + i_top, 0.0);
			Fill_iMatrix(correlation_ranking_recur, L_temp, recursion_size_temp + i_top, 0);

			//Fill in top correlation or rank firstly
			for (int k = 0; k < L_temp; k++) {
				for (int t = 0; t < i_top; t++) {
					correlation_recur[k][t] = correlation_top_temp[k][t];
					correlation_ranking_recur[k][t] = correlation_ranking_top_temp[k][t];
				}
			}
		}
/*
		if (mynode == 1) {
			if (i_recur == 0) {
				cout << " correlation_recur matrix at  i_recur = " << i_recur << endl;
				for (int kk2 = 0; kk2 < L_temp; kk2++) {
					for (int kk3 = 0; kk3 < recursion_size_temp; kk3++) {
						cout << setw(20) << correlation_recur[kk2][kk3];
					}
					cout << endl;
				}
				cout << " correlation_ranking_recur matrix at  i_recur = " << i_recur << endl;
				for (int kk2 = 0; kk2 < L_temp; kk2++) {
					for (int kk3 = 0; kk3 < recursion_size_temp; kk3++) {
						cout << setw(20) << correlation_ranking_recur[kk2][kk3];
					}
					cout << endl;
				}
			}
			if (i_recur != 0) {
				cout << " correlation_recur matrix at  i_recur = " << i_recur << endl;
				for (int kk2 = 0; kk2 < L_temp; kk2++) {
					for (int kk3 = 0; kk3 < (recursion_size_temp + top); kk3++) {
						cout << setw(20) << correlation_recur[kk2][kk3];
					}
					cout << endl;
				}
				cout << " correlation_ranking_recur matrix at  i_recur = " << i_recur << endl;
				for (int kk2 = 0; kk2 < L_temp; kk2++) {
					for (int kk3 = 0; kk3 < (recursion_size_temp + top); kk3++) {
						cout << setw(20) << correlation_ranking_recur[kk2][kk3];
					}
					cout << endl;
				}
			}
		
		}*/
		//if ( (mynode != 0) && (i_recur == (recursion_time - 1)) ) {
		//	correlation_recur = New_dMatrix(L_temp, recursion_size_temp + top);
		//	correlation_ranking_recur = New_dMatrix(L_temp, recursion_size_temp + top);
		//}


		//-----------------------
		//Compute corrrlation matrix
		//-----------------------

		double* x1 = new double[nrow_ol];
		double* x2 = new double[nrow_ol];
		double d_sum = 0.0;

		for (int j = 0; j<L_temp; j++) //distributed variables
		{

			for (int j_next = 0; j_next<recursion_size_temp; j_next++) //recursively read variables
			{
				for (int i = 0; i<nrow_ol; i++)
				{
					x1[i] = ol_matrix[i][j];   //jth column
					x2[i] = ol_recursion[i][j_next];//next column
				}

				//---
				//each column's mean
				//---
				double x1_mean = 0.0; double x2_mean = 0.0;
				for (int i = 0; i<nrow_ol; i++)
				{
					x1_mean += x1[i];   //jth column
					x2_mean += x2[i];   //next column
				}
				x1_mean = x1_mean / nrow_ol;
				x2_mean = x2_mean / nrow_ol;
				//TestOut<<"x1_mean is "<< x1_mean <<", and x2_mean is "<< x2_mean <<" at j_next = "<< j_next <<endl;
				//-----
				//calculate covariance of two columns
				//-----
				d_sum = 0.0;
				for (int i_1 = 0; i_1<nrow_ol; i_1++)
				{
					d_sum += (x1[i_1] - x1_mean)*(x2[i_1] - x2_mean);
				}

				//----------------
				//calculate variance of two columns
				//----------------
				double x1_var = 0.0; double x2_var = 0.0;
				double var_sum = 0.0;
				for (int i_2 = 0; i_2 < nrow_ol;i_2++) {
					var_sum = var_sum + (x1[i_2] - x1_mean)*(x1[i_2] - x1_mean);
				}
				x1_var = var_sum;

				var_sum = 0.0;
				for (int i_3 = 0; i_3 < nrow_ol;i_3++) {
					var_sum = var_sum + (x2[i_3] - x2_mean)*(x2[i_3] - x2_mean);
				}
				x2_var = var_sum;

				//TestOut << "d_sum is " << d_sum  <<" and x1_var is "<< x1_var <<", x2_var is "<< x2_var << " at j_next = " << j_next << endl;
				//---------
				//store covariance using symmetry property
				//---------
				//if(mynode==1) cout<<"d_sum / sqrt(x1_var* x2_var) is "<< d_sum / sqrt(x1_var* x2_var) <<" at j="<<j<<" at j_next="<< j_next <<" at i_recur="<< i_recur <<endl;

				if(i_recur == 0) correlation_recur[j][j_next] = d_sum / sqrt(x1_var* x2_var);
				if(i_recur != 0) correlation_recur[j][i_top + j_next] = d_sum / sqrt(x1_var* x2_var);

				if (i_recur == 0) correlation_ranking_recur[j][j_next] = startpoint + j_next + 1;//actual location
				if (i_recur != 0) correlation_ranking_recur[j][i_top + j_next] = startpoint + j_next + 1;//actual location
			} //end of j_next

			/*if (mynode == 1) {
				if (i_recur == 0) {
					cout << " correlation_recur matrix at  i_recur = " << i_recur << " at j = " << j << endl;
					for (int kk2 = 0; kk2 < L_temp; kk2++) {
						for (int kk3 = 0; kk3 < recursion_size_temp; kk3++) {
							cout << setw(20) << correlation_recur[kk2][kk3];
						}
						cout << endl;
					}
					cout << " correlation_ranking_recur matrix at  i_recur = " << i_recur << " at j = " << j << endl;
					for (int kk2 = 0; kk2 < L_temp; kk2++) {
						for (int kk3 = 0; kk3 < recursion_size_temp; kk3++) {
							cout << setw(20) << correlation_ranking_recur[kk2][kk3];
						}
						cout << endl;
					}
				}
				if (i_recur != 0) {
					cout << " correlation_recur matrix at  i_recur = " << i_recur << " at j = " << j << endl;
					for (int kk2 = 0; kk2 < L_temp; kk2++) {
						for (int kk3 = 0; kk3 < (recursion_size_temp + top); kk3++) {
							cout << setw(20) << correlation_recur[kk2][kk3];
						}
						cout << endl;
					}
					cout << " correlation_ranking_recur matrix at  i_recur = " << i_recur << " at j = " << j << endl;
					for (int kk2 = 0; kk2 < L_temp; kk2++) {
						for (int kk3 = 0; kk3 < (recursion_size_temp + top); kk3++) {
							cout << setw(20) << correlation_ranking_recur[kk2][kk3];
						}
						cout << endl;
					}
				}
			
			}*/

			std::vector<int> i_return;

			if (i_recur == 0) {

				double* d_source_temp = new double[recursion_size_temp];
				double* d_source = new double[recursion_size_temp];
				int* i_source = new int[recursion_size_temp];

				for (int k = 0; k < recursion_size_temp; k++) {
					d_source_temp[k] = abs(correlation_recur[j][k]);// Note the ranking of correlation matrix should be based on absolute value !!!
					d_source[k] = abs(correlation_recur[j][k]);// Note the ranking of correlation matrix should be based on absolute value !!!
					i_source[k] = correlation_ranking_recur[j][k];
				}


				std::sort(d_source_temp, d_source_temp + recursion_size_temp);

				//cout<<"d_source_temp at j="<<j<<" at i_recur="<< i_recur <<endl;
				//for (int k = 0; k < recursion_size_temp; k++) {
				//	cout << setw(20) << d_source_temp[k];
				//}
				//cout << endl;

				for (int t = 0; t < i_top; t++) {
					correlation_top_temp[j][t] = d_source_temp[recursion_size_temp - 1 - t]; //include the rank of itself
				}

				//Sort in ascending order and return mapping
				//Example: d_source            = [1.1, 0.7]
				//         i_return            = [2, 1]
				//correlation_ranking_top_temp = [2, 1, 0, 0, 0, 0, 0, 0] 
				order_FHDI(d_source, i_source, recursion_size_temp, i_return);

				for (int t = 0; t < i_top; t++) {
					correlation_ranking_top_temp[j][t] = i_return[recursion_size_temp - 1 - t]; //exclude the rank of itself
				}

				//if (mynode == 1) {
				//	cout << " correlation_top_temp matrix at  i_recur = " << i_recur << endl;
				//	for (int kk2 = 0; kk2 < L_temp; kk2++) {
				//		for (int kk3 = 0; kk3 < top; kk3++) {
				//			cout << setw(20) << correlation_top_temp[kk2][kk3];
				//		}
				//		cout << endl;
				//	}
				//	cout << " correlation_ranking_top_temp matrix at  i_recur = " << i_recur << endl;
				//	for (int kk2 = 0; kk2 < L_temp; kk2++) {
				//		for (int kk3 = 0; kk3 < top; kk3++) {
				//			cout << setw(20) << correlation_ranking_top_temp[kk2][kk3];
				//		}
				//		cout << endl;
				//	}

				//}

				delete[] d_source_temp;
				delete[] d_source;
				delete[] i_source;
			} 

			if (i_recur != 0) {
				double* d_source_temp = new double[recursion_size_temp + i_top];
				double* d_source = new double[recursion_size_temp + i_top];
				int* i_source = new int[recursion_size_temp + i_top];

				for (int k = 0; k < (recursion_size_temp + i_top); k++) {
					d_source_temp[k] = abs(correlation_recur[j][k]);// Note the ranking of correlation matrix should be based on absolute value !!!	
					d_source[k] = abs(correlation_recur[j][k]);// Note the ranking of correlation matrix should be based on absolute value !!!
					i_source[k] = correlation_ranking_recur[j][k];//TestOut<<"cov["<<k<<"] : "<< abs(cov[k]) <<" at j = "<<j<<endl;
				}


				std::sort(d_source_temp, d_source_temp + (recursion_size_temp + i_top));

				for (int t = 0; t < i_top; t++) {
					correlation_top_temp[j][t] = d_source_temp[recursion_size_temp + i_top - 1 - t]; //exclude the rank of itself
				}

				//cout<<"i_source at i_recur="<< i_recur<<" at j="<<j<<endl;
				//for (int m = 0;m < (recursion_size_temp + top);m++) {
				//	cout << setw(20) << i_source[m];
				//}
				//cout << endl;
				//cout << "d_source at i_recur=" << i_recur << " at j=" << j << endl;
				//for (int m = 0;m < (recursion_size_temp + top);m++) {
				//	cout << setw(20) << d_source[m];
				//}
				//cout << endl;

				//Sort in ascending order and return mapping
				//ncol = 8, top = 8, recustion_temp = 2
				//Example: d_source            = [1.1, 0.7, 0, 0, 0, 0, 0, 0, 0.9, 0.1] size of 10
				//         i_source            = [2, 1, 0, 0, 0, 0, 0, 0, 3, 4] size of 10
				//sorted d_source_temp         = [0, 0, 0, 0, 0, 0, 0.1, 0.7, 0.9, 1.1] size of 10
				//correlation_top_temp         = [1,1, 0.9, 0.7, 0.1, 0, 0, 0, 0] size of 8

				//         i_return            = [0, 0, 0, 0, 0, 0, 4, 1, 3, 2] size of 10
				//correlation_ranking_top_temp = [2, 3, 4, 1, 0, 0, 0, 0] size of 8

				order_FHDI(d_source, i_source, recursion_size_temp + i_top, i_return);

				//cout<<"i_return at i_recur="<< i_recur<<" at j="<<j<<endl;
				//for (int m = 0;m < i_return.size();m++) {
				//	cout << setw(20) << i_return[m];
				//}
				//cout << endl;

				for (int t = 0; t < i_top; t++) {
					correlation_ranking_top_temp[j][t] = i_return[recursion_size_temp + i_top - 1 - t]; //exclude the rank of itself
				}

				//if (mynode == 1) {
				//	cout << " correlation_top_temp matrix at  i_recur = " << i_recur <<" at j="<<j<< endl;
				//	for (int kk2 = 0; kk2 < L_temp; kk2++) {
				//		for (int kk3 = 0; kk3 < top; kk3++) {
				//			cout << setw(20) << correlation_top_temp[kk2][kk3];
				//		}
				//		cout << endl;
				//	}
				//	cout << " correlation_ranking_top_temp matrix at  i_recur = " << i_recur <<" at j="<<j<< endl;
				//	for (int kk2 = 0; kk2 < L_temp; kk2++) {
				//		for (int kk3 = 0; kk3 < top; kk3++) {
				//			cout << setw(20) << correlation_ranking_top_temp[kk2][kk3];
				//		}
				//		cout << endl;
				//	}

				//}

				delete[] d_source_temp;
				delete[] d_source;
				delete[] i_source;
			}

		}//end of L_temp

		 //---------
		 //Deallocation
		 //---------
		delete[] x1;
		delete[] x2;

		 //if (mynode == 1) {
			// cout << " correlation_top_temp matrix at  i_recur = " << i_recur << endl;
			// for (int kk2 = 0; kk2 < L_temp; kk2++) {
			//	 for (int kk3 = 0; kk3 < top; kk3++) {
			//		 cout << setw(20) << correlation_top_temp[kk2][kk3];
			//	 }
			//	 cout << endl;
			// }
			// cout << " correlation_ranking_top_temp matrix at  i_recur = " << i_recur<< endl;
			// for (int kk2 = 0; kk2 < L_temp; kk2++) {
			//	 for (int kk3 = 0; kk3 < top; kk3++) {
			//		 cout << setw(20) << correlation_ranking_top_temp[kk2][kk3];
			//	 }
			//	 cout << endl;
			// }
		 //
		 //}

		if ((mynode != 0) && (i_recur == 0)) {
			Del_dMatrix(correlation_recur, L_temp, recursion_size_temp);
			Del_iMatrix(correlation_ranking_recur, L_temp, recursion_size_temp);
		}
		if ((mynode != 0) && (i_recur != 0)) {
			Del_dMatrix(correlation_recur, L_temp, recursion_size_temp + i_top);
			Del_iMatrix(correlation_ranking_recur, L_temp, recursion_size_temp + i_top);
		}

		if (mynode != 0) { 
			Del_dMatrix(ol_recursion, nrow_ol, recursion_size_temp); 
			delete[] array_temp2;
		}
	}//end of recursion


	//TestOut
	//if (mynode == 3) {
	//	cout << "final correlation_top_temp matrix at node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < L_temp; kk2++) {
	//		for (int kk3 = 0; kk3 < top; kk3++) {
	//			cout << setw(20) << correlation_top_temp[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//	cout << "final correlation_ranking_top_temp matrix at node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < L_temp; kk2++) {
	//		for (int kk3 = 0; kk3 < top; kk3++) {
	//			cout << setw(20) << correlation_ranking_top_temp[kk2][kk3];
	//		}
	//		cout << endl;
	//	}

	//}

	//Send distributed top correlation and top ranking to the master processor
	if (mynode != 0) {
		MPI_Send(correlation_top_temp[0], L_temp*i_top, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		MPI_Send(correlation_ranking_top_temp[0], L_temp*i_top, MPI_INT, 0, 2, MPI_COMM_WORLD);
	}

	//All_gather distributed top correlation and top ranking on the master processor
	counter = 0;
	if (mynode == 0) {
		double** correlation_top_temp_recv = New_dMatrix(numWorkPerProc, i_top);
		double** correlation_top_temp_recv_last = New_dMatrix(numWorkLocalLast, i_top);

		int** correlation_ranking_top_temp_recv = New_iMatrix(numWorkPerProc, i_top);
		int** correlation_ranking_top_temp_recv_last = New_iMatrix(numWorkLocalLast, i_top);

		for (int j = 1; j < totalnodes; j = j + 1) {
			if (j != (totalnodes - 1)){
				MPI_Recv(correlation_top_temp_recv[0], numWorkPerProc*i_top, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
				MPI_Recv(correlation_ranking_top_temp_recv[0], numWorkPerProc*i_top, MPI_INT, j, 2, MPI_COMM_WORLD, &status);

				//Note that correlation_top and correlation_ranking_top should exclude correlation or ranking of itself
				//correlation_top_temp and correlation_ranking_top_temp has correlation or ranking of itself where top = i_top - 1
				for (int k = 0; k < numWorkPerProc; k++) {
					for (int h = 1; h < i_top; h++) {
						correlation_top[counter][h-1] = correlation_top_temp_recv[k][h];
						correlation_ranking_top[counter][h-1] = correlation_ranking_top_temp_recv[k][h];
					}
					counter++;
				}
			}

			if (j == (totalnodes - 1)){
				MPI_Recv(correlation_top_temp_recv_last[0], numWorkLocalLast*i_top, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
				MPI_Recv(correlation_ranking_top_temp_recv_last[0], numWorkLocalLast*i_top, MPI_INT, j, 2, MPI_COMM_WORLD, &status);

				for (int k = 0; k < numWorkLocalLast; k++) {
					for (int h = 1; h < i_top; h++) {
						correlation_top[counter][h-1] = correlation_top_temp_recv_last[k][h];
						correlation_ranking_top[counter][h-1] = correlation_ranking_top_temp_recv_last[k][h];
					}
					counter++;
				}
			}

		}//end of slace nodes

		 	//cout << "final correlation_top matrix at node " << mynode <<" where top is "<<top<< endl;
		 	//for (int kk2 = 0; kk2 < ncol; kk2++) {
		 	//	for (int kk3 = 0; kk3 < top-1; kk3++) {
		 	//		cout << setw(20) << correlation_top[kk2][kk3];
		 	//	}
		 	//	cout << endl;
		 	//}
		 	//cout << "final correlation_ranking_top matrix at node " << mynode <<" where top is "<<top<< endl;
		 	//for (int kk2 = 0; kk2 < ncol; kk2++) {
		 	//	for (int kk3 = 0; kk3 < top-1; kk3++) {
		 	//		cout << setw(20) << correlation_ranking_top[kk2][kk3];
		 	//	}
		 	//	cout << endl;
		 	//}

		Del_dMatrix(correlation_top_temp_recv, numWorkPerProc, i_top);
		Del_dMatrix(correlation_top_temp_recv_last, numWorkLocalLast, i_top);

		Del_iMatrix(correlation_ranking_top_temp_recv, numWorkPerProc, i_top);
		Del_iMatrix(correlation_ranking_top_temp_recv_last, numWorkLocalLast, i_top);

	}//end of master node

	//Broadcast completed top correlation and top ranking to all slave processors
	MPI_Bcast(correlation_top[0], ncol*(i_top -1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(correlation_ranking_top[0], ncol*(i_top - 1), MPI_INT, 0, MPI_COMM_WORLD);

	//if (mynode == 3) {
	//	cout << "final correlation_top matrix at node " << mynode <<" where i_top is "<< i_top << endl;
	//	for (int kk2 = 0; kk2 < ncol; kk2++) {
	//		for (int kk3 = 0; kk3 < i_top -1; kk3++) {
	//			cout << setw(20) << correlation_top[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//	cout << "final correlation_ranking_top matrix at node " << mynode <<" where i_top is "<< i_top << endl;
	//	for (int kk2 = 0; kk2 < ncol; kk2++) {
	//		for (int kk3 = 0; kk3 < i_top -1; kk3++) {
	//			cout << setw(20) << correlation_ranking_top[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//}

	success = MPI_File_close(&fh_binary_daty);
	if (success != MPI_SUCCESS) cout << "MPI I/O fail to close the file!" << endl;

	Del_dMatrix(ol_matrix, nrow_ol, L_temp);
	if (mynode != 0) {
		Del_dMatrix(correlation_top_temp, L_temp, i_top);
		Del_iMatrix(correlation_ranking_top_temp, L_temp, i_top);
	}

	return;

}