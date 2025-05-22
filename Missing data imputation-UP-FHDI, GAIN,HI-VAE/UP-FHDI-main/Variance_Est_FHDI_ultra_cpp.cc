
void Variance_Est_FHDI_ultra_cpp(MPI_File fh_binary_daty, const int nrow, const int ncol,
	RepWeight_FHDI &d_rw, int* id, std::vector<int> ol,
	double** simp_fmat_FHDI, int i_row_fmat,
	List_FHDI &List_nU, List_FHDI &uox_infor, List_FHDI &mox_infor, int nrow_mox,
	ofstream& TestOut_Slave3) {
	//Description----------------------
	//estimate variance for FHDI using Jackknife method 
	//  Algorithm: 
	//
	// original R code: Dr. Im, J. and Dr. Kim, J. 
	// c++ code: 		Dr. Cho, I. and Yicheng Yang
	// All rights reserved
	// 
	// updated: August 16, 2021
	//
	//----------------------------------------------------
	//IN   : MPI_File fh_binary_daty  = raw daty written column-wisely
	//IN   : const int nrow           = number of rows of daty
	//IN   : const int ncol           = number of columns of daty
	//IN   : RepWeight_FHDI &d_rw     = replicate weights
	//IN   : int* id                  = global id of daty
	//IN   : std::vector<int> ol      = indices of fully observed rows in daty
	//IN   : double** simp_fmat_FHDI  = matrix of 4 columns: global id + sampling weight + wij + fwij
	//IN   : MPI_File fh_fmat_FHDI    = file string to write imputed matrix in size of (4+ncol) columns: global id + sampling weight + wij + fwij + imputed daty
	//IN   : int i_row_fmat           = number of rows of imputed values
	//IN   : MPI_File fh_mox	      = file string to read mox
	//IN   : List_FHDI List_nU        = donors list of all mox in uox
	//IN   : List_FHDI uox_infor      = Actual index list of uox in z
	//IN   : List_FHDI mox_infor      = Actual index list of mox in z
	//IN   : int nrow_mox             = number of rows of mox

	//OUT  : TestOut_Slave3
	//-------------------------------------------------------

   	//===========================
	// MPI variables
	//===========================

	//double startup = MPI_Wtime();
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	int success = 0;
	MPI_File fh_mox_reading;
	MPI_File fh_fmat_FHDI_reading; // File to write imputed values

	success = MPI_File_open(MPI_COMM_WORLD, "./mox_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_mox_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of mox fail to open the file mox_binary.bin!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./fmat_FHDI_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_fmat_FHDI_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of mox fail to open the file fmat_FHDI_binary.bin!" << endl;


	//--------------------
	//put 1 into fully observed rows
	//--------------------
	std::vector<double> d_rr0;
	for (int t = 0; t < nrow; t++) {
		d_rr0.push_back(0.0);
	}

	int i_ol = (int)ol.size();;
	for (int k = 0; k < i_ol; k++) {
		d_rr0[ol[k] - 1] = 1.0;
	}

	//if (mynode == 0) {
	//	for (unsigned int t = 0; t < d_rr0.size(); t++) {
	//		TestOut_Slave3 <<"d_rr0["<<t<<"]: "<< d_rr0[t]<<endl;
	//	}
	//}


	//---------------------------
	//Job assignment
	//---------------------------
	const int i_nrow_imputation = i_row_fmat;

	const int L = nrow; //size of d_rw 
	int numWorkPerProc = (int)floor(1.0*L / (1.0*totalnodes - 1));
	int numWorkLocalLast = L - numWorkPerProc * (totalnodes - 2);

	int L_temp = 0;
	if (mynode != (totalnodes - 1)) L_temp = numWorkPerProc;
	if (mynode == (totalnodes - 1)) L_temp = numWorkLocalLast;

	//===============================
	// specify startpoint and end point for slave processors
	//===============================
	int startpoint = 0;
	int endpoint = 0;

	if (mynode != 0 && mynode != (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = mynode*numWorkPerProc;

		//RPrint("(Variance)Strating point and ending point on node ");RPrint(mynode);
		//RPrint("are: \n");
		//RPrint(startpoint); RPrint(endpoint);
		if (endpoint - startpoint != L_temp) {
			cout << "y_bar_i_k boundary ERROR!!!" << endl;
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
			cout << "y_bar_i_k boundary ERROR!!!" << endl;
			return;
		}
	}

	//cout<<"Variance i_nrow_imputation is "<< i_nrow_imputation <<" where startpoint is "<< startpoint <<" and endpoint is "<< endpoint <<" at node "<<mynode<<endl;
	//-----------------------
	//Initialization
	//--------------------------
	const int ncol_simp = 4;// 4 columns: global id + sampling weight + wij + fwij 
	const int ncol_fmat = ncol + 4; // 4 columns + ncol : global id + sampling weight + wij + fwij + imputed values

	//double* rw0 = new double[nrow];//buffer of replicate weights
	//double* Rw = new double[i_nrow_imputation];//replicate weights of imputed values
	//double* wijk = new double[i_nrow_imputation]; //fractional weights

	//double* wmatLocal = new double[i_nrow_imputation];//

	std::vector<double> rw0;
	std::vector<double> Rw;
	std::vector<double> wijk;
	std::vector<double> wmatLocal;

	double* yi = new double[ncol];
	double** y_bar_i_k = NULL;
	if (mynode != 0) {
		y_bar_i_k = New_dMatrix(L_temp, ncol);
		Fill_dMatrix(y_bar_i_k, L_temp, ncol, 0.0);
	}

	//if (mynode == 1) {
	//	cout << "Final simp_fmat_FHDI at node " << mynode << endl;

	//	for (int kk2 = 0; kk2 < i_nrow_imputation; kk2++) {
	//		for (int kk3 = 0; kk3 < 4; kk3++) {
	//			cout << setw(20) << simp_fmat_FHDI[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
	//}

	int* rbind_ipmat = new int[i_nrow_imputation];
	for (int kk = 0; kk < i_nrow_imputation;kk++) {
		rbind_ipmat[kk] = (int)simp_fmat_FHDI[kk][0];//1st column
	}

	//if (mynode == 1) {
	//	for (int t = 0; t < i_nrow_imputation; t++) cout <<"rbind_ipmat["<<t<<"]: "<< rbind_ipmat[t]<<endl;
	//	for (int t = 0; t < nrow; t++) cout << "id[" << t << "]: " << id[t] << endl;
	//}

	//====================================
    //TEstOut
	//double* rw0_test = new double[nrow];//buffer of replicate weights
	//double* Rw_test = new double[i_nrow_imputation];//replicate weights of imputed values

	//for (int i = 0; i<nrow; i++) rw0_test[i] = d_rw(i, 0); //l_th column 
	//for (int i = 0; i < i_nrow_imputation; i++) {
	//	Rw_test[i] = rw0_test[rbind_ipmat[i] - 1];
	//}

	////if (mynode == 1) {

	////	for (int j = 0; j < nrow; j++) cout << "rw0_test[" << j << "]: " << rw0_test[j] << endl;
	////	for (int j = 0; j < i_nrow_imputation; j++) cout << "Rw_test[" << j << "]: " << Rw_test[j] << endl;
	////}

	//delete[] rw0_test;
	//delete[] Rw_test;
	//======================================

	MPI_Barrier(MPI_COMM_WORLD);

	double* d_1_mox = new double[ncol]; //buffer to read mox
	double* d_1_daty = new double[ncol]; //buffer to read daty
	double* d_1_imputed = new double[ncol_fmat]; //buffer to read daty
	std::vector<int> v_lg; //Actual locations of mox whose donors' fractional weights need updated
	std::vector<double> mox_buffer;

	std::vector<int> loc_srst_nl; // buffer to hold donors of each mox
	int i_size_loc_srst_nl = 0;
	std::vector<int> v_obsg0; // buffer to hold all possible donors of a mox in z matrix (original)
	std::vector<int> uox_buffer;
	//------------------------
	//Test normalization
	//if (mynode == 0) {
	//	double** daty_test = New_dMatrix(4, 6);
	//	for (int p = 0; p < 4; p++) {
	//		for (int j = 0; j < 6; j++) {
	//			daty_test[p][j] = (double)p;
	//		}
	//	}

	//	cout << "daty_test at node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < 4; kk2++) {
	//		for (int kk3 = 0; kk3 < 6; kk3++) {
	//			cout << setw(20) << daty_test[kk2][kk3];
	//		}
	//		cout << endl;
	//	}

	//	double** daty_normal = New_dMatrix(4, 6);
	//	normalization(daty_test, daty_normal, 4,6);

	//	cout << "daty_normal at node " << mynode << endl;
	//	for (int kk2 = 0; kk2 < 4; kk2++) {
	//		for (int kk3 = 0; kk3 < 6; kk3++) {
	//			cout << setw(20) << daty_normal[kk2][kk3];
	//		}
	//		cout << endl;
	//	}
 //		//----------------------

	//}
	//-------------------------
	//if (mynode != 0) cout << " Get started MAMA at node " << mynode << endl;
	//------------------------------
	//------------------------------
	//Main loop for L replications
	//------------------------------
	//------------------------------
	for (int l = startpoint; l < endpoint; l++) {
		double FHDI_iteration1 = MPI_Wtime();
		//-------
		//replicate weight from lth column
		//-------
		//Fill_dVector(rw0, ncol, 0.0); 
		rw0.clear();
		for (int i = 0; i < nrow; i++) rw0.push_back((double)d_rw(i, l)); //l_th column 
		//for (int i = 0; i<nrow; i++) rw0[i] = (double)d_rw(i, l); //l_th column 

		//if (mynode !=0 ) cout << "L_Replication at L = " << l << " at node " << mynode << endl;

		//Fill_dVector(Rw, i_nrow_imputation, 0.0);
		Rw.clear();
		for (int i = 0; i < i_nrow_imputation; i++) {
			//Rw[i] = rw0[rbind_ipmat[i] - 1];
			Rw.push_back((double)rw0[rbind_ipmat[i] - 1]);
			//if (l == 0) cout << "Rw[" << i << "]: " << Rw[i] << " at l = " << l << " at node " << mynode << endl;
		}
		
		//Fill_dVector(wijk, i_nrow_imputation, 0.0);
		wijk.clear();
		//if (mynode !=0 ) cout << "variance7 at L = " << l << " at node " << mynode << endl;
		//---------------------------------------
		//Initilization of fractional weights 
		//---------------------------------------
		for (int k = 0; k < i_nrow_imputation; k++) {
			//wijk[k] = (double)simp_fmat_FHDI[k][2];//3rd column
			wijk.push_back((double)simp_fmat_FHDI[k][2]);
			//if(l ==0) cout <<"wijk["<<k<<"]: "<< wijk[k] <<" at l = "<<l<<" at node "<<mynode<<endl;
		}
		//if (mynode!=0) cout << "variance0 at L = " << l <<" at node "<<mynode<< endl;
		//----------------------------------------
		//Updating fractional weights wijk
		//
		//1. if the deleted is missing unit, no further action is taken
		//2. if the deleted is observed unit, then the fractional weights of imputed values of several mox are re-computed 
		//----------------------------------------

		if (fabs(d_rr0[l]) > 0) {
			//---------------------
			//locations of the deleted unit in observed list
			//---------------------
			//std::vector<int> v_lg; //Actual locations of mox whose donors' fractional weights need updated
			v_lg.clear();
			for (int i = 0; i < nrow_mox; i++) //all missing patterns
			{
				loc_srst_nl.clear();
				v_obsg0.clear();
				List_nU.get_block_yicheng(i, loc_srst_nl);

				i_size_loc_srst_nl = (int)loc_srst_nl.size();

				if (i_size_loc_srst_nl == 0) { cout<<"ERROR! Number of donors of a mox is zero in Jackknife variance estimation!"<<endl; }

				for (int b = 0; b < i_size_loc_srst_nl; b++) {
					uox_buffer.clear();
					uox_infor.get_block_yicheng(loc_srst_nl[b] - 1, uox_buffer);

					for (int t = 0; t < uox_buffer.size(); t++) {
						v_obsg0.push_back(uox_buffer[t]);
					}
				}

				for (int k = 0; k < v_obsg0.size(); k++) //donor rows for the jth missing pattern
				{
					if ((l + 1) == (v_obsg0[k])) {
						v_lg.push_back(i+1);
						break;
					}
				}
			}

			const int nlg = (int)v_lg.size();
			//cout<<"nlg is "<< nlg <<" at l="<<l<<" at node "<<mynode<<endl;
			//if (mynode == (totalnodes - 1)) {
			//	for (int p = 0; p < nlg; p++) {
			//		cout<<"v_lg["<<p<<"]: "<< v_lg[p]<<" at l = "<<l<<endl;
			//	}
			//}

			for (int j = 0; j < nlg; j++){

				int i_row_lg = 0;
				i_row_lg = (int)v_lg[j] - 1; // one mox

				//read onw row of mox
				MPI_In_uox_mox(ncol, i_row_lg, fh_mox_reading, d_1_mox);

				//---
				//actual col number of missing cell in current missing row
				//---
				std::vector<int> v_rloc; v_rloc.clear();
				for (int k = 0; k<ncol; k++)
				{
					if (fabs(d_1_mox[k]) < 1e-15) { v_rloc.push_back(k + 1); } //actual col
				}

				int nrloc = (int)v_rloc.size();

				//if (mynode == 1) {
				//	cout << "d_1_mox[" << i_row_lg << "] where nrloc is "<< nrloc << endl;
				//	for (int t = 0; t < ncol; t++) {
				//		cout << setw(20) << d_1_mox[t];
				//	}
				//	cout << endl;
				//}

				//-------
				//location of this mox in z matrix 
				//-------
				std::vector<int> v_mlog; v_mlog.clear();
				mox_buffer.clear();
				mox_infor.get_block_yicheng(i_row_lg, mox_buffer);
				for (int t = 0; t < mox_buffer.size(); t++) {
					v_mlog.push_back((int)mox_buffer[t]);
				}

				//for (int t = 2; t < (2 + max_overlap_size); t++) {
				//	if (mox_info_final[i_row_lg][t] != 0) {
				//		v_mlog.push_back(mox_info_final[i_row_lg][t]);
				//	}
				//}

				const int nmlog = (int)v_mlog.size();//number of a mox in z matrix
				//if (mynode !=0 ) {
				//	for (int t = 0; t < nmlog; t++) cout << "v_mlog[" << t << "]: " << v_mlog[t] <<" at j= "<<j<< endl;
				//}

				if (nmlog == 0) cout<<"ERROR in nmlog at l = "<<l<<endl;
				//-----------------------------------
				//locate these mox in imputed indices
				//------------------------------------
				std::vector<int> v_elog; v_elog.clear();// actual locations of chosen imputed values in rbind

				//binary search can't be used here because rbind_ipmat is not sorted
				for (int k1 = 0; k1 < nmlog; k1++)
				{
					for (int k2 = 0; k2 < i_nrow_imputation; k2++) {
						if (v_mlog[k1] == rbind_ipmat[k2]) {
							v_elog.push_back(k2 + 1);
						}
					}
				}

				const int i_size_v_elog = (int)v_elog.size();//number of chosen imputed values 
				if (i_size_v_elog == 0) cout<<"ERROR in v_elog at l = "<<l<<endl;
				//for (int p = 0; p < i_size_v_elog; p++) {
				//	cout<<"v_elog["<<p<<"]: "<< v_elog[p]<<" at l="<<l<<" at j="<<j<<" at node "<<mynode<<endl;
				//}

				//tank to hold daty[l] + chosen donors
				int i_size_v_elog1 = 0;
				i_size_v_elog1 = i_size_v_elog + 1;

				double** dy_FHDI = New_dMatrix(i_size_v_elog1, nrloc);

				//-----------------------
				// read lth row of daty
				//----------------------
				MPI_In_daty_row(nrow, ncol, l, fh_binary_daty, d_1_daty);
				//if (mynode == 1) {
				//	cout<<"d_1_daty at l = "<<l<<endl;
				//	for (int t = 0; t < ncol; t++) {
				//		cout << setw(20) << d_1_daty[t];
				//	}
				//	cout << endl;
				//}

				for (int t = 0; t < nrloc; t++) dy_FHDI[0][t] = d_1_daty[v_rloc[t]-1];

				//-----------------------
				// read chosen imputed values
				//----------------------

				int counter = 1;
				for (int p = 0; p < i_size_v_elog; p++) {
					int i_temp = 0;
					i_temp = (int)v_elog[p]-1;
					MPI_In_uox_mox(ncol_fmat, i_temp, fh_fmat_FHDI_reading, d_1_imputed);//only for validation purpos

				//if (mynode == 1) {
				//	cout<<"d_1_imputed at l = "<<l<<" and p = "<<p<<endl;
				//	for (int t = 0; t < ncol + 4; t++) {
				//		cout << setw(20) << d_1_imputed[t];
				//	}
				//	cout << endl;
				//}

					for(int t=0; t < nrloc; t++) dy_FHDI[counter][t] = d_1_imputed[v_rloc[t] - 1 + 4];// read from 5th columns
					counter++;
				}

				//if (mynode != 0) {
				//	cout<<"dy_FHDI at l = "<<l<<" and j = "<<j<<endl;
				//	for (int p = 0; p < i_size_v_elog1; p++) {
				//		for (int k = 0; k < nrloc; k++) {
				//			cout << setw(20) << dy_FHDI[p][k];
				//		}
				//		cout << endl;
				//	}
				//}

				//Normalize lth daty + chosen imputed values to (0~1)
				double** dy_FHDI_normalized = New_dMatrix(i_size_v_elog1, nrloc);
				normalization(dy_FHDI, dy_FHDI_normalized, i_size_v_elog1, nrloc);

				//if (mynode != 0) {
				//	cout << "dy_FHDI_normalized at l = " << l << " and j = " << j << endl;
				//	for (int p = 0; p < i_size_v_elog1; p++) {
				//		for (int k = 0; k < nrloc; k++) {
				//			cout << setw(20) << dy_FHDI_normalized[p][k];
				//		}
				//		cout << endl;
				//	}

				//}

				//-----------------------------------------------------------------
				//compute Euclidean distance between lth daty and chosen imputed values
				//-----------------------------------------------------------------

				double* d_score = new double[i_size_v_elog];
				double d_sum_dist = 0.0;
				for (int k = 0; k < i_size_v_elog; k++) {
					d_sum_dist = 0.0;
					for (int m = 0; m < nrloc; m++) {
						double d_daty_temp = dy_FHDI_normalized[0][m];
						double d_imputed_temp = dy_FHDI_normalized[k+1][m];
						d_sum_dist += (d_daty_temp - d_imputed_temp)*(d_daty_temp - d_imputed_temp);
					}
					d_score[k] = d_sum_dist;
				}

				//if (mynode != 0) {
				//	cout << "d_score at l = " << l << endl;
				//	for (int t = 0; t < i_size_v_elog; t++) {
				//		cout << setw(20) << d_score[t];
				//	}
				//	cout << endl;
				//}

				const int MM = i_size_v_elog / nmlog; //number of donors of a mox
				if (MM < 2) cout<<"ERROR in MM!!!!!"<<endl;

				std::vector<int> i_eloc;// actual locations of candidates who have the minimum Euclidean distance
				i_eloc.clear();
				std::vector<double> d_score_temp;// score buffer to hold donors of a mox
				int min = 0; // location of donor who has the minimum Euclidean distance

				for (int k1 = 0; k1 < nmlog; k1++){
					d_score_temp.clear();
					for (int k2 = 0; k2 < MM; k2++) {
						d_score_temp.push_back(d_score[k1*MM + k2]);
					}

					min = 0;
					for (int k3 = 0; k3 < MM; k3++) {
						if (d_score_temp[min] > d_score_temp[k3]) {
							min = k3;
						}
					}

					i_eloc.push_back(k1*MM + min + 1);
				}

				//if (mynode !=0 ) {
				//	cout << "i_eloc at l = " << l << endl;
				//	for (int t = 0; t < i_eloc.size(); t++) {
				//		cout << setw(20) << i_eloc[t];
				//	}
				//	cout << endl;
				//}

				//actual locations in v_elog excluding i_eloc
				std::vector<int> i_without_eloc;
				i_without_eloc.clear();
				for (int m = 0; m < i_size_v_elog; m++) {
					for (int m1 = 0; m1 < nmlog; m1++) {
						if (m != (i_eloc[m1] - 1) ) {
							i_without_eloc.push_back(m + 1);
						}
					}
				}

				//if (mynode != 0) {
				//	cout << "i_without_eloc at l = " << l << endl;
				//	for (int t = 0; t < i_without_eloc.size(); t++) {
				//		cout << setw(20) << i_without_eloc[t];
				//	}
				//	cout << endl;
				//}
				//-----------------------------------------
				//extract weights
				//-----------------------------------------
				//-----------
				//when zero size continue to next iteration
				//-----------

				double* ewijk = new double[i_size_v_elog];
				double* fefiw = new double[i_size_v_elog];

				for (int k1 = 0; k1<i_size_v_elog; k1++)
				{
					ewijk[k1] = wijk[v_elog[k1] - 1]; //-1 for actual loc
													  //5th column is FEFIW of fhdi[[2]] in r version 
					fefiw[k1] = (double)simp_fmat_FHDI[v_elog[k1] - 1][3]; //5th column
				}

				//if (mynode == 1) {
				//	cout << "Before ewijk at l = " << l << endl;
				//	for (int t = 0; t < i_size_v_elog; t++) {
				//		cout << setw(20) << ewijk[t];
				//	}
				//	cout << endl;

				//}

				std::vector<double> d_maxew(nmlog);// w_{ir} - w_{ir,FEFI}
				d_maxew.clear();
				std::vector<double> d_maxval(nmlog);// w_{ir,FEFI}
				d_maxval.clear();

				for (int k1 = 0; k1<nmlog; k1++)
				{
					double d_temp = ewijk[i_eloc[k1] - 1] - fefiw[i_eloc[k1] - 1]; //-1 for actual loc
					d_maxew[k1] = 0.0;
					if (d_temp > 0.0) d_maxew[k1] = d_temp;

					d_maxval[k1] = ewijk[i_eloc[k1] - 1] - d_maxew[k1];
				}

				//------
				//ewijk update
				//------
				for (int k1 = 0; k1<nmlog; k1++)
				{
					ewijk[i_eloc[k1] - 1] = d_maxew[k1];
				}

				//if (mynode != 0) {
				//	cout << "After ewijk at l = " << l << endl;
				//	for (int t = 0; t < i_size_v_elog; t++) {
				//		cout << setw(20) << ewijk[t];
				//	}
				//	cout << endl;
				//}



				std::vector<double> d_maxval_extended(i_size_v_elog - nmlog); // w_{ij,FEFI}/ [sum(w_ij,FEFI)-1]
				d_maxval_extended.clear();
				int i_maxval = 0;
				for (int k1 = 0; k1<nmlog; k1++)
				{
					double d_temp = (double)d_maxval[k1] / (MM - 1);
					for (int k2 = 0; k2<(MM - 1); k2++)
					{
						d_maxval_extended[i_maxval++] = d_temp;
					}
				}

				//if (mynode != 0) {
				//	cout << "d_maxval_extended at l = " << l << endl;
				//	for (int t = 0; t < i_size_v_elog - nmlog; t++) {
				//		cout << setw(20) << d_maxval_extended[t];
				//	}
				//	cout << endl;

				//}
				//---------------
				//update ewijk with extended maxval
				//---------------
				//exclude eloc locations
				//-----------
				for (int k1 = 0; k1<(i_size_v_elog - nmlog); k1++)
				{
					int i_loc_ew = (int)i_without_eloc[k1] - 1; //-1 for actual loc
					ewijk[i_loc_ew] = ewijk[i_loc_ew] + d_maxval_extended[k1];
				}

				//if (mynode != 0) {
				//	cout << "Final ewijk at l = " << l << endl;
				//	for (int t = 0; t < i_size_v_elog; t++) {
				//		cout << setw(20) << ewijk[t];
				//	}
				//	cout << endl;
				//}

				//----------------
				//final update wijk with ewijk
				//----------------
				for (int k1 = 0; k1<i_size_v_elog; k1++)
				{
					wijk[v_elog[k1] - 1] = (double)ewijk[k1];
				}

				//----------
				//local deallocation 
				//----------
				delete[] d_score;
				delete[] ewijk;
				delete[] fefiw;

				Del_dMatrix(dy_FHDI, i_size_v_elog1, nrloc);
				Del_dMatrix(dy_FHDI_normalized, i_size_v_elog1, nrloc);

			}//end of mox

		}//end of updating fractional weights

		//if (mynode != 0) cout << "variance1 at L = " << l <<" at node "<<mynode<< endl;
		wmatLocal.clear();
		for (int k4 = 0; k4 < i_nrow_imputation; k4++)
		{
			//wmatLocal[k4] = Rw[k4] * wijk[k4];
			//double d_temp3 = 0.0;
			//d_temp3 = (double)Rw[k4] * wijk[k4];
			wmatLocal.push_back(Rw[k4] * wijk[k4]);
			//if (l == 0) cout << "wmatLocal[" << k4 << "]: " << wmatLocal[k4] << " at l = " << l << " at node " << mynode << endl;
		}


		double d_sum_wij = 0.0;
		Fill_dVector(yi, ncol, 0.0); //initialize vector for column-wise means of all variables  

		//if (mynode== (totalnodes-1)) {
		//	cout << "FHDI_iteration1 time at l = " << l << " = " << MPI_Wtime() - FHDI_iteration1 << endl;
		//}

		double FHDI_iteration2 = MPI_Wtime();
		for (int i = 0; i < i_nrow_imputation; i++) {
			double wi = 0.0;
			wi = (double)simp_fmat_FHDI[i][1];
			double wij = 0.0;
			wij = (double)wmatLocal[i];

			d_sum_wij = d_sum_wij + wi*wij;

			MPI_In_uox_mox(ncol_fmat, i, fh_fmat_FHDI_reading, d_1_imputed);//only for validation purpos

			for (int i_var = 0; i_var < ncol; i_var++) {
				yi[i_var] = yi[i_var] + wi*wij * d_1_imputed[4 + i_var];
			}

		}

		//if (mynode == (totalnodes - 1)) {
		//	cout << "FHDI_iteration2 time at l = " << l << " = " << MPI_Wtime() - FHDI_iteration2 << endl;
		//}

		if (fabs(d_sum_wij) == 0.0)
		{
			cout << "ERROR! zero sum of fractional weight at Jackknifed row :" << l << endl;
			return;
		}

		for (int i_var = 0; i_var < ncol; i_var++) //size of columns of ipmat matrix of C++
		{
			//cout << "K: " << k << ", i_loc: " << i_loc << ", d_sum_wij: " << d_sum_wij << endl;
			double d_temp = (double)yi[i_var] / d_sum_wij;

			//-----
			//bar{y}_i^(k)
			//-----
			//NOTE: R works column-by-column 
			//hence, below sequence is different from C++ ordering 
			//final_full_data[i_var*nrow + i] = d_temp;  //note: i=current row
			y_bar_i_k[l - (mynode - 1)*numWorkPerProc][i_var] = d_temp;  //note: k= replicate; i_var = variable
																		
		}

		//if (mynode == (totalnodes - 1)) {
		//	cout << "FHDI each time at l = " << l << " = " << MPI_Wtime() - FHDI_iteration1 << endl;
		//}

	
	}//end of main loop

	MPI_Barrier(MPI_COMM_WORLD);

	//======================================
	//======================================
	//Extract Variance
	//======================================

	//---------------------------------
	//Jackknife average of y_bar_i_k
	//--------------------------------

	if (mynode != 0) {
		std::vector<double> y_bar_i;
		for (int i_var = 0; i_var<ncol; i_var++)
		{
			double d_temp = 0.0;
			for (int k = 0; k<L_temp; k++)
			{
				d_temp += y_bar_i_k[k][i_var];
			}

			y_bar_i.push_back(d_temp);
		}

		MPI_Send(&y_bar_i[0], ncol, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}


	std::vector<double> y_bar_i_mean(ncol);// mean of y_bar_i_k

	if (mynode == 0) {
		std::vector<double> y_bar_i_recv(ncol);
		std::vector<double> y_bar_i_mean_temp(ncol);

		for (int j = 1; j < totalnodes; j++) {
			MPI_Recv(&y_bar_i_recv[0], ncol, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
			for (int k = 0; k < ncol; k++) {
				y_bar_i_mean_temp[k] = y_bar_i_mean_temp[k] + y_bar_i_recv[k];
			}
		}

		for (int t = 0; t < ncol; t++) y_bar_i_mean[t] = y_bar_i_mean_temp[t] / nrow;
	}

	//Broadcast mean of y_bar_i_k
	MPI_Bcast(&y_bar_i_mean[0], ncol, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//if (mynode == 1) {
	//	cout<<"Mean: "<<endl;
	//	for (int t = 0; t < ncol; t++) {
	//		cout << setw(20) << y_bar_i_mean[t];
	//	}
	//	cout << endl;
	//}

	//---------------------
	//Final Jackknife Variance Estimation of the mean 
	//---------------------
	if (mynode != 0) {
		std::vector<double> variance_data(ncol); //y_bar_i_k - y_bar_i_mean

		for (int i_var = 0; i_var < ncol; i_var++) {

			double d_temp = 0;
			d_temp = 0.0;
			for (int k = 0; k<L_temp; k++) //Jackknife replicate
			{
				d_temp +=
					(y_bar_i_k[k][i_var] - y_bar_i_mean[i_var])
					*(y_bar_i_k[k][i_var] - y_bar_i_mean[i_var]);
			}
			variance_data[i_var] = d_temp;
		}

		MPI_Send(&variance_data[0], ncol, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

	}

	if (mynode == 0) {
		std::vector<double> variance_data_recv(ncol);
		std::vector<double> variance_sum(ncol);
		std::vector<double> final_variance_data(ncol);
		for (int j = 1; j < totalnodes; j++) {
			MPI_Recv(&variance_data_recv[0], ncol, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
			for (int k = 0; k < ncol; k++) {
				variance_sum[k] = variance_sum[k] + variance_data_recv[k];
			}
		}
		for (int k = 0; k < ncol; k++) {
			final_variance_data[k] = variance_sum[k]* (nrow - 1) / nrow;
		}

		TestOut_Slave3 << "Jackknife Variance Results:" << endl;
		for (int i = 0; i < ncol; i++) {
			//TestOut << "bbfore" << endl;
			//TestOut_Slave3 << setw(20) << final_variance_data[i];
			TestOut_Slave3 << final_variance_data[i] << endl;
			//TestOut << "aafter" << endl;
		}
		//TestOut_Slave3 << endl;
	}

	success = MPI_File_close(&fh_mox_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of mox fail to close the file!" << endl;

	success = MPI_File_close(&fh_fmat_FHDI_reading);
	if (success != MPI_SUCCESS) cout << "MPI I/O of imputed matrix fail to close the file!" << endl;
	//======================================

	//==========================
	//Deallocation
	//===========================
	//delete[] rw0;
	//delete[] Rw;
	//delete[] wijk;
	//delete[] wmatLocal;
	delete[] yi;
	delete[] rbind_ipmat;

	delete[] d_1_mox;
	delete[] d_1_daty;
	delete[] d_1_imputed;

	if (mynode != 0) {
		Del_dMatrix(y_bar_i_k, L_temp, ncol);
	}



	return;
}

