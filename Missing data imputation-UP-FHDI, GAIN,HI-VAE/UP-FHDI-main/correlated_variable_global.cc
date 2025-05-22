
void correlated_variable_global(const int ncol, const int i_option_collapsing, int* ia_temp,
	double **correlation_yicheng, std::vector<int> &v_mxl, ofstream& TestOut)

	//Description=========================================
	//Select the most i_option_collapsing correlated variables from all observed variables of each mox
	//Algorithm:
	//Select the most correlated variables using majority vote. The last one is selected based on the highest correlation

	//IN    : int i_option_collapsing = choice of big-p algorithm 
	//                               0= no big-p algorithms
	//                              !0= perform big-p algorithms
	//IN    : int ia_temp(ncol) = copy of mox[i]
	//IN    : double correlation_yicheng(ncol, ncol);// correlation matrix

	//OUT   : int v_mxl(i_option_collapsing); // the actual location of most correlated variables of mox[i]
	//=====================================================

{

	std::vector<int> v_lm; //temporary vector for the locaiton of missing values in mox

	whichINVNOT(ia_temp, ncol, 0, v_lm); //get the actual location of missing variables of mox[i] 

	int v_lm_size = v_lm.size(); //number of missing variables in mox[i]

	std::vector<int> v_lo; //temporary vector for the actual locaiton of observed values in mox[i]

						   //for (int a1 = 0; a1 < ncol; a1++) {
						   //	for (int a2 = 0; a2 < v_lm_size; a2++) {
						   //		if ((a1 + 1) != v_lm[a2]) {
						   //			continue;
						   //			
						   //		}
						   //		v_lo.push_back(a1 + 1);
						   //	}
						   //}
	whichINV(ia_temp, ncol, 0, v_lo); //get the location of Non-zero in mox 

	int v_lo_size = v_lo.size(); //number of observed variables in mox[i]

								 //for (int a3 = 0; a3 < v_lo_size; a3++) {
								 //	TestOut <<"v_lo["<<a3<<"]: "<< v_lo[a3] << endl;
								 //}

	if (v_lm_size + v_lo_size != ncol) TestOut << "Error in correlated_variable_gloabl!!!!!" << endl;
	//========================================================
	// Pickout the correlation matrix of missing variables only
	//========================================================

	double** correlation_temp = New_dMatrix(v_lm_size, (ncol - v_lm_size));

	//TestOut << "v_lm at mox i=  "<< endl;
	//for (int kk1 = 0; kk1 < v_lm_size; kk1++) {
	//	TestOut << v_lm[kk1] << endl;
	//}

	for (int k6 = 0; k6 < v_lm_size; k6++) {
		for (int k7 = 0; k7 < (ncol - v_lm_size); k7++) {
			correlation_temp[k6][k7] = abs(correlation_yicheng[v_lm[k6] - 1][v_lo[k7] - 1]);
		}
	}


	std::vector<double> r_star; // vector of the highest correlation of each obeserved variable of current mox
	std::vector<double> r_temp; // temporary vector to hold coreelations of between an observed variable and all missing variables

								//=================================
								// Select out the vector of the highest correlation of each obeserved variable of current mox
	for (int k8 = 0; k8 < v_lo_size; k8++) {
		r_temp.clear();
		for (int k9 = 0; k9 < v_lm_size; k9++) {
			r_temp.push_back(correlation_temp[k9][k8]);
		}

		r_star.push_back(max_FHDI(r_temp));

	}

	for (int b1 = 0; b1 < r_star.size(); b1++) {
		TestOut << "r_star[" << b1 << "]: " << r_star[b1] << endl;
	}
	//===================================
	// Add variables whose correlation is among the top of the largest i_option_collapsing

	for (int k10 = 0; k10 < i_option_collapsing; k10++) {
		int max_corr = 0;

		for (int k11 = 0; k11 < v_lo_size; k11++) {
			if (r_star[max_corr] < r_star[k11]) {
				max_corr = k11;
			}
		}

		//TestOut<<"max_corr: "<< max_corr <<endl;
		v_mxl.push_back(v_lo[max_corr]); // the actual location of the variable

		r_star[max_corr] = 0;
	}


	sort(v_mxl.begin(), v_mxl.end());

	if (v_mxl.size() != i_option_collapsing) {
		TestOut << "ERROE! The global ranking of correlation is wrong to get " << i_option_collapsing << " selected variables" << endl;
		exit(0);
	}
	//TestOut << "correlation_temp2 after max_occur at mox i=  " << i << endl;
	//for (int kk2 = 0; kk2 < v_lm_size; kk2++) {
	//	for (int kk3 = 0; kk3 < (ncol - v_lm_size); kk3++) {
	//		TestOut << setw(20) << correlation_temp2[kk2][kk3];
	//	}
	//	TestOut << endl;
	//}

	//----------------------------------
	//Del_iMatrix(correlation_temp, v_lm_size, (ncol - 1));
	Del_dMatrix(correlation_temp, v_lm_size, (ncol - v_lm_size));

	return;
}

void correlated_variable_global2(const int ncol, const int i_option_collapsing, int nrow_ol, int* ia_temp,
	double** ol_matrix, std::vector<int> &v_mxl, ofstream& TestOut)

	//Description=========================================
	//Select the most i_option_collapsing correlated variables from all observed variables of each mox
	//Algorithm:
	//Select the most correlated variables using majority vote. The last one is selected based on the highest correlation

	//IN    : int i_option_collapsing = choice of big-p algorithm 
	//                               0= no big-p algorithms
	//                              !0= perform big-p algorithms
	//IN    : int ia_temp(ncol) = copy of mox[i]
	//IN    : double correlation_yicheng(ncol, ncol);// correlation matrix

	//OUT   : int v_mxl(i_option_collapsing); // the actual location of most correlated variables of mox[i]
	//=====================================================

{

	std::vector<int> v_lm; //temporary vector for the locaiton of missing values in mox

	whichINVNOT(ia_temp, ncol, 0, v_lm); //get the actual location of missing variables of mox[i] 

	int v_lm_size = v_lm.size(); //number of missing variables in mox[i]

	std::vector<int> v_lo; //temporary vector for the actual locaiton of observed values in mox[i]

						   //for (int a1 = 0; a1 < ncol; a1++) {
						   //	for (int a2 = 0; a2 < v_lm_size; a2++) {
						   //		if ((a1 + 1) != v_lm[a2]) {
						   //			continue;
						   //			
						   //		}
						   //		v_lo.push_back(a1 + 1);
						   //	}
						   //}
	whichINV(ia_temp, ncol, 0, v_lo); //get the actual location of Non-zero in mox 

	int v_lo_size = v_lo.size(); //number of observed variables in mox[i]

								 //for (int a3 = 0; a3 < v_lo_size; a3++) {
								 //	TestOut <<"v_lo["<<a3<<"]: "<< v_lo[a3] << endl;
								 //}

	if (v_lm_size + v_lo_size != ncol) TestOut << "Error in correlated_variable_gloabl!!!!!" << endl;
	//========================================================
	// Pickout the correlation matrix of missing variables only
	//========================================================

	//double** correlation_temp = New_dMatrix(v_lm_size, (ncol - v_lm_size));

	//TestOut << "v_lm at mox i=  "<< endl;
	//for (int kk1 = 0; kk1 < v_lm_size; kk1++) {
	//	TestOut << v_lm[kk1] << endl;
	//}

	//for (int k6 = 0; k6 < v_lm_size; k6++) {
	//	for (int k7 = 0; k7 < (ncol - v_lm_size); k7++) {
	//		correlation_temp[k6][k7] = abs(correlation_yicheng[v_lm[k6] - 1][v_lo[k7] - 1]);
	//	}
	//}



	std::vector<double> r_star; // vector of the highest correlation of each obeserved variable of current mox
	std::vector<double> r_temp; // temporary vector to hold coreelations of between an observed variable and all missing variables

								//=================================
								// Select out the vector of the highest correlation of each obeserved variable of current mox

	double* x1 = new double[nrow_ol];
	double* x2 = new double[nrow_ol];
	double d_sum = 0.0;

	for (int k8 = 0; k8 < v_lo_size; k8++) {

		r_temp.clear();

		for (int k9 = 0; k9 < v_lm_size; k9++) {

			for (int k10 = 0; k10<nrow_ol; k10++)
			{
				x1[k10] = ol_matrix[k10][v_lo[k8] - 1];
				x2[k10] = ol_matrix[k10][v_lm[k9] - 1];
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
			//TestOut << "x1_mean is " << x1_mean << ", and x2_mean is " << x2_mean << " at j_next = " << j_next << endl;
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

			r_temp.push_back(abs(d_sum / sqrt(x1_var* x2_var)));
			//r_temp.push_back(correlation_temp[k9][k8]);
		}

		r_star.push_back(max_FHDI(r_temp));

	}


	//for (int b1 = 0; b1 < r_star.size(); b1++) {
	//	TestOut <<"r_star_Top["<<b1<<"]: "<< r_star[b1] << endl;
	//}
	//===================================
	// Add variables whose correlation is among the top of the largest i_option_collapsing

	for (int k10 = 0; k10 < i_option_collapsing; k10++) {
		int max_corr = 0;

		for (int k11 = 0; k11 < v_lo_size; k11++) {
			if (r_star[max_corr] < r_star[k11]) {
				max_corr = k11;
			}
		}

		//TestOut<<"max_corr: "<< max_corr <<endl;
		v_mxl.push_back(v_lo[max_corr]); // the actual location of the variable

		r_star[max_corr] = 0;
	}

	//TestOut<<"v_mxl top size is "<< v_mxl.size()<<endl;

	sort(v_mxl.begin(), v_mxl.end());

	if (v_mxl.size() != i_option_collapsing) {
		TestOut << "ERROE! The global ranking of correlation is wrong to get " << i_option_collapsing << " selected variables" << endl;
		exit(0);
	}
	//TestOut << "correlation_temp2 after max_occur at mox i=  " << i << endl;
	//for (int kk2 = 0; kk2 < v_lm_size; kk2++) {
	//	for (int kk3 = 0; kk3 < (ncol - v_lm_size); kk3++) {
	//		TestOut << setw(20) << correlation_temp2[kk2][kk3];
	//	}
	//	TestOut << endl;
	//}

	//----------------------------------
	//Del_iMatrix(correlation_temp, v_lm_size, (ncol - 1));
	//Del_dMatrix(correlation_temp, v_lm_size, (ncol - v_lm_size));

	//---------
	//Deallocation
	//---------
	delete[] x1;
	delete[] x2;

	return;
}


void correlated_variable_global_ultra(const int ncol, const int i_option_collapsing, const int top, int* ia_temp,
	double** correlation_top, int** correlation_ranking_top, std::vector<int> &v_mxl, ofstream& TestOut)

	//Description=========================================
	//Select the most i_option_collapsing correlated variables from all observed variables of each mox
	//Algorithm:
	//Select the most correlated variables using majority vote. The last one is selected based on the highest correlation

	//IN    : ncol                    = number of variables in raw data
	//IN    : int i_option_collapsing = choice of big-p algorithm 
	//                              0 = no big-p algorithms
	//                             !0 = perform big-p algorithms
	//IN    : int top                 = number of top correlation or correlation ranking
	//IN    : int ia_temp(ncol) = copy of mox[i]
	//IN    : int** correlation_ranking_top = top correlation ranking of all variables
	//IN    : int** correlation_top         = top correlation of all variables
	//INOUT   : int v_mxl(i_option_collapsing)= the actual location of selected variables of mox[i]
	//=====================================================

{

	std::vector<int> v_lm; //temporary vector for the locaiton of missing values in mox

	whichINVNOT(ia_temp, ncol, 0, v_lm); //get the actual location of missing variables of mox[i] 

	int v_lm_size = v_lm.size(); //number of missing variables in mox[i]

	std::vector<int> v_lo; //temporary vector for the actual locaiton of observed values in mox[i]

	whichINV(ia_temp, ncol, 0, v_lo); //get the location of Non-zero in mox 

	int v_lo_size = v_lo.size(); //number of observed variables in mox[i]

	//for (int a3 = 0; a3 < v_lo_size; a3++) {
	//	TestOut <<"v_lo["<<a3<<"]: "<< v_lo[a3] << endl;
	//}

	if (v_lm_size + v_lo_size != ncol) cout << "Error in correlated_variable_gloabl!!!!!" << endl;
	//========================================================
	// Pickout the correlation matrix of missing variables only
	//========================================================


	std::vector<double> r_star; // vector of the highest correlation of each obeserved variable of current mox

	//for (int k = 0; k < v_lo.size(); k++) {
	//	cout<<"v_lo["<<k<<"]: "<< v_lo[k]<<endl;
	//}
	//for (int k = 0; k < v_lm.size(); k++) {
	//	cout << "v_lm[" << k << "]: " << v_lm[k] << endl;
	//}
	//=================================
	// Select out the vector of the highest correlation of each obeserved variable of current mox
	//Example: mox = [2,0,0,0,3,0,3,2]
	//    v2  v3  v4  v6
	//v1  x   x   x   x
	//v5  
	//v7  
	//v8  

	for (int k8 = 0; k8 < v_lo_size; k8++) {
		for (int k9 = 0; k9 < top; k9++) {
			for (int k10 = 0; k10 < v_lm_size; k10++) {
				if (correlation_ranking_top[v_lo[k8] - 1][k9] == v_lm[k10]) {
					r_star.push_back(correlation_top[v_lo[k8] - 1][k9]);
					k10 = v_lm_size;//breaking out nested for loop
					k9 = top;//breaking out nested for loop
				}
				if (k10 == (v_lm_size - 1) && k9 == (top - 1)) { cout << "ERROR!!! One can not find corresponding correlation for qualified variables" << endl; return;}
			}
		}
	}

	//for (int k8 = 0; k8 < v_lo_size; k8++) {
	//	r_temp.clear();
	//	for (int k9 = 0; k9 < v_lm_size; k9++) {
	//		r_temp.push_back(correlation_temp[k9][k8]);
	//	}

	//	r_star.push_back(max_FHDI(r_temp));

	//}

	if (r_star.size() != v_lo_size) { cout << "ERROR in correlated_variable_global_ultra function" << endl; return; }

	//for (int b1 = 0; b1 < r_star.size(); b1++) {
	//	cout << "r_star[" << b1 << "]: " << r_star[b1] << endl;
	//}
	//===================================
	// Add variables whose correlation is among the top of the largest i_option_collapsing

	for (int k10 = 0; k10 < i_option_collapsing; k10++) {
		int max_corr = 0;

		for (int k11 = 0; k11 < v_lo_size; k11++) {
			if (r_star[max_corr] < r_star[k11]) {
				max_corr = k11;
			}
		}

		//TestOut<<"max_corr: "<< max_corr <<endl;
		v_mxl.push_back(v_lo[max_corr]); // the actual location of the variable

		r_star[max_corr] = 0;
	}


	sort(v_mxl.begin(), v_mxl.end());

	if (v_mxl.size() != i_option_collapsing) {
		cout << "ERROE! The global ranking of correlation is wrong to get " << i_option_collapsing << " selected variables" << endl;
		exit(0);
	}
	//TestOut << "correlation_temp2 after max_occur at mox i=  " << i << endl;
	//for (int kk2 = 0; kk2 < v_lm_size; kk2++) {
	//	for (int kk3 = 0; kk3 < (ncol - v_lm_size); kk3++) {
	//		TestOut << setw(20) << correlation_temp2[kk2][kk3];
	//	}
	//	TestOut << endl;
	//}


	return;
}





