
void Fractional_Hot_Deck_Imputation_ultra(const int i, MPI_File fh_binary_daty_row,
	            const int ncol, const int nrow,
				double* mox_temp, const int i_M, const int i_merge,
				const int i_mxl, 
				std::vector<int> v_cn_z_i, std::vector<int> v_mxl,
				std::vector<int> v_obsg, std::vector<double> fwij,
	            double d_myran, double* w, int* id, int &counter_fmat,
	            double** simp_fmat_FHDI_temp, double** fmat_FHDI_temp)
//Description----------------------
//perform FHDI 
//  Algorithm: impute by using "M" possible donors
//  final outcome is "fmat" in which each column means that
//  col1: global id
//  col2: sampling weight
//  col3: fractional weights (wij)
//  col4: fractional weights (fwij)
//  col5: imputed original data (matrix with column of ncol) 
//
//
// original R code: Dr. Im, J. and Dr. Kim, J. 
// c++ code: 		Dr. Cho, I. 
// All rights reserved
// 
// updated: August 10, 2021
//
//IN   : int i                     = the current row 
//IN   : MPI_File fh_binary_daty   = raw daty written column-wisely
//IN   : double* mox_temp          = ith row of mox
//IN   : int i_M                   = user-defined number of possible donors
//IN   : const int i_merge         = control of random number generator: 0 (fixed); 1 (random)
//IN   : const int ncol            = number of columns of daty
//IN   : const int nrow            = number of rows of daty
//IN   : int i_mxl                 = number of non-missing variables of mox[i]
//IN   : std::vector<int> v_cn_z_i = actual locations of mox[i] in z matrix
//IN   : std::vector<int> v_mxl    = observed variables of mox[i] 
//IN   : std::vector<int> v_obsg   = actual locations of all possible donors of mox[i]
//IN   : double* fwij              = fractional weights of all possible donors of mox[i]                          
//IN   : double* w[nrow]           = sampling weights of daty. Default as 1
//IN   : int* id[nrow]             = global indices of daty

//OUT  : int counter_fmat          = global counter for fmat_FHDI_temp to continue from adding fully observed rows
//OUT  : double** simp_fmat_FHDI_temp = distributed matrix of 4 columns: global id + sampling weight + wij + fwij
//OUT  : double** fmat_FHDI_temp      = distributed matrix of (4+ncol) columns: global id + sampling weight + wij + fwij + imputed daty
//----------------------
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;
	//------------
	//initial data setting
	//------------
	const int i_size_v_obsg = v_obsg.size();
	const int i_size_v_cn_z_i = (int)v_cn_z_i.size(); //same as mg and loc2
	
	//if (mynode == 1) cout<<"d_myran is "<< d_myran <<" at node "<<mynode<<endl;
		
	//----
	//(a) determine MM depending upon donors and i_M
	//----
	int MM = 0; 
	if(i_size_v_obsg <= i_M) MM=i_size_v_obsg; 
	if(i_size_v_obsg >  i_M) MM=i_M; 

    //e.g.,
	//i_size_v_cn_z_i = 2 and MM = 5
	//i_SN = [1,2,3,4,5,1,2,3,4,5]
	int* i_SN = new int[MM*i_size_v_cn_z_i];// actual locations of selected donors in v_obsg
	double** d_iy = New_dMatrix(MM*i_size_v_cn_z_i, ncol); //imputed matrix
	double* wij = new double[MM*i_size_v_cn_z_i]; //fractional weights
	std::vector<double> fwij_selected; //selected fractional weights associated with selected dononors
	fwij_selected.clear();

	//int success = 0;

	//MPI_File fh_daty;
	//success = MPI_File_open(MPI_COMM_WORLD, "./daty_column_binary.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_daty);
	//if (success != MPI_SUCCESS) cout << "MPI I/O fail to open the file!" << endl;
	
	//cout << "MM is " << MM << " at i = " << i << " and i_M is " << i_M << " and i_size_v_obsg is " << i_size_v_obsg<<" and i_size_v_cn_z_i is "<< i_size_v_cn_z_i << endl;
	//------------------------
	//FHDI Case 1: when donors <= M
	//             use all possible donors
	//------------------------
	int counter = 0;
	double* array_temp = new double[ncol]; //buffer to hold one row of daty
	int temp = 0;

	if(i_size_v_obsg <= i_M)
	{
		//----
		//index of all donors 
		//----

		for(int j=0; j<i_size_v_cn_z_i; j++)
		{
			for (int k = 0; k < MM; k++) {
				i_SN[counter] = k + 1; //+1 for actual location 
				counter++;
			}
		}
		
		//--------
		//get ready donors from original data
		//--------
		for(int j=0; j<MM*i_size_v_cn_z_i; j++)
		{
			temp = 0;
			temp = v_obsg[i_SN[j] - 1] - 1;
			//MPI_In_daty_row(nrow, ncol, temp, fh_binary_daty, array_temp);// read distributed z matrix row by row
			MPI_In_uox_mox(ncol, temp, fh_binary_daty_row, array_temp);

			//if (mynode == 1 && i == 0) {
			//	cout<<"array_temp: "<<endl;
			//	for (int t = 0; t < ncol; t++) {
			//		cout << setw(20) << array_temp[t];
			//	}
			//	cout << endl;
			//}

			for(int k=0; k<ncol; k++)
			{ 
				d_iy[j][k] = array_temp[k]; // get donors of mox[i] from daty
			} 
		}

		//if (mynode == 1 && i == 3) {
		//	cout<<"d_iy before at i="<<i<<" at node "<<mynode<<endl;
		//	for (int p = 0; p < MM*i_size_v_cn_z_i; p++) {
		//		for (int m = 0; m < ncol; m++) {
		//			cout << setw(20) << d_iy[p][m];
		//		}
		//		cout << endl;
		//	}
		//}

		//--------
		//non-missing column consideration of current row
		//Replace d_iy at observde variables of mox[i] with original raw data
		//e.g.,
		//i_size_v_cn_z_i = 2 where v_cn_z_i = [2,6]
		//mox[i] : 3 , 0, 2, 0;
 		//daty[2]: 1.1, 0, 2.2, 0;
		//daty[6]: 2.1, 0, 3.2, 0;
		//It has three donors:
		// 0.1, 0.2, 0.3, 0.4
		// 0.3, 0.5, 0.6, 0.7
		// 0.4, 0.1, 0.9, 1.0
		//Thus, d_iy is:
		//
		// 1.1, 0.2, 2.2, 0.4
		// 1.1, 0.5, 2.2, 0.7
		// 1.1, 0.1, 2.2, 1.0

		// 2.1, 0.2, 3.2, 0.4
		// 2.1, 0.5, 3.2, 0.7
		// 2.1, 0.1, 3.2, 1.0
		//--------
		counter = 0;
		for (int j = 0; j < i_size_v_cn_z_i; j++) {
			temp = 0; 
			temp = v_cn_z_i[j] - 1;
			//MPI_In_daty_row(nrow, ncol, temp, fh_binary_daty, array_temp);// read distributed z matrix row by row
			MPI_In_uox_mox(ncol, temp, fh_binary_daty_row, array_temp);
		    
			for (int k = 0; k < MM; k++) {
				for (int t = 0; t < i_mxl; t++) {
					d_iy[counter][v_mxl[t] - 1] = array_temp[v_mxl[t] - 1];
				}
				counter++;
			}
		}
		//if (mynode == 1 && i == 3) {
		//	cout << "d_iy after at i=" << i << " at node " << mynode << endl;
		//	for (int p = 0; p < MM*i_size_v_cn_z_i; p++) {
		//		for (int m = 0; m < ncol; m++) {
		//			cout << setw(20) << d_iy[p][m];
		//		}
		//		cout << endl;
		//	}
		//}
		
		//------------
		//weights 
		//------------

		for (int j = 0; j<i_size_v_cn_z_i; j++)
		{
			for (int t = 0; t < MM; t++) {
				fwij_selected.push_back(fwij[t]);
			}
		}

		if ( fwij_selected.size() != (MM*i_size_v_cn_z_i) ) cout<<"ERROR!!!! Selected fwij is incorrect!!!1"<<endl;

		for (int t = 0; t < i_size_v_cn_z_i*MM; t++)
		{
			wij[t] = fwij_selected[t];
		}
		
	}//end of Case 1
	

	//-----------------------------------------------------
	//FHDI Case 2: when donors > M
	//             use tailed systematic sampling to select M donors
	//-----------------------------------------------------


	if(i_size_v_obsg > i_M){

		//cout << "mox " << i << " has i_size_v_obsg = " << i_size_v_obsg << " donors" << endl;

		//cout << "v_cn_z_i at mox = " << i << endl;
		//for (int t = 0; t < i_size_v_cn_z_i; t++) {
		//	cout << "v_cn_z_i[" << t << "]: " << v_cn_z_i[t] << endl;
		//}

		//--------------------------------------------
		//cumulative sum of joint probability with pps 
		//----------------------------------------------
		double* d_cs_temp = new double[i_size_v_obsg];
		double* d_Ui = new double[i_size_v_obsg];

		for (int j = 0; j < i_size_v_obsg; j++)
		{
			d_cs_temp[j] = fwij[j] * MM;
		}

		//Construct interval for sampling scheme
		cumsum_FHDI(d_cs_temp, i_size_v_obsg, d_Ui);

		//if (mynode == 1) {
		//	cout << "d_Ui at mox " << i << endl;
		//	for (int t = 0; t<i_size_v_obsg; t++)
		//	{
		//		cout << setw(20) << d_Ui[t];
		//	}
		//	cout << endl;
		//}

		//random generator
		double d_Rg = 0.0;
		d_Rg = d_myran;

		//Random numbers shoot in the interval
		int* i_rmg = new int[i_size_v_cn_z_i * MM];
		int* i_rM = new int[i_size_v_cn_z_i * MM];
		double* d_SR = new double[i_size_v_cn_z_i * MM];

		for (int j = 0; j<i_size_v_cn_z_i; j++)
		{
			for (int k = 0; k<MM; k++) i_rmg[j*MM + k] = j + 1;
		}

		for (int j = 0; j<i_size_v_cn_z_i; j++)
		{
			for (int k = 0; k<MM; k++)  i_rM[j*MM + k] = k + 1;
		}

		for (int j = 0; j<i_size_v_cn_z_i * MM; j++)
		{
			d_SR[j] = (d_Rg + (i_rmg[j] - 1)) / i_size_v_cn_z_i + (i_rM[j] - 1);
		}

		//if (mynode == 1) {
		//	cout << "i_rmg at mox " << i << endl;
		//	for (int t = 0; t<i_size_v_cn_z_i*MM; t++)
		//	{
		//		cout << setw(20) << i_rmg[t];
		//	}
		//	cout << endl;

		//	cout << "i_rM at mox " << i << endl;
		//	for (int t = 0; t<i_size_v_cn_z_i*MM; t++)
		//	{
		//		cout << setw(20) << i_rM[t];
		//	}
		//	cout << endl;

		//	cout << "d_SR at mox " << i << endl;
		//	for (int t = 0; t<i_size_v_cn_z_i*MM; t++)
		//	{
		//		cout << setw(20) << d_SR[t];
		//	}

		//	cout << endl;
		//}
		

		//-------------------------------
		//Determine id of selecetd donors
		//set of SR < Ui
		//---------------------------------

		int i_SR_Ui = 0;
		for (int k = 0; k < i_size_v_cn_z_i*MM; k++)
		{
			i_SR_Ui = 0; //minimum location where SR <= Ui
			for (int j = 0; j< i_size_v_obsg; j++)
			{
				if (d_SR[k] <= d_Ui[j]) { i_SR_Ui = j + 1; break; }
			}

			i_SN[k] = i_SR_Ui; //Actual location stored 
		}

		//if (mynode == 1) {
		//	cout << "i_SN is below i_size_v_cn_z_i = " << i_size_v_cn_z_i <<" at mox "<<i<< endl;
		//	for (int j = 0; j < MM*i_size_v_cn_z_i; j++)
		//	{
		//		cout << setw(20) << i_SN[j];
		//	}
		//	cout << endl;
		//}
		//---------------------------------------------------------------


		//-----------------------------------
		//get ready donors from original data
		//------------------------------------

		for (int j = 0; j< MM*i_size_v_cn_z_i; j++)
		{
			temp = 0;
			temp = v_obsg[i_SN[j] - 1] - 1;
			//MPI_In_daty_row(nrow, ncol, temp, fh_binary_daty, array_temp);// read distributed z matrix row by row
			MPI_In_uox_mox(ncol, temp, fh_binary_daty_row, array_temp);


			for (int k = 0; k<ncol; k++)
			{
				d_iy[j][k] = array_temp[k];
			} //-1 for actual location 
		}

		//if (mynode == 1 && i == 3) {
		//	cout << "d_iy before at i=" << i << " at node " << mynode << endl;
		//	for (int p = 0; p < MM*i_size_v_cn_z_i; p++) {
		//		for (int m = 0; m < ncol; m++) {
		//			cout << setw(20) << d_iy[p][m];
		//		}
		//		cout << endl;
		//	}
		//}

		//--------
		//non-missing column consideration of current row
		//--------
		counter = 0;
		for (int j = 0; j < i_size_v_cn_z_i; j++) {
			temp = 0;
			temp = v_cn_z_i[j] - 1;
			//MPI_In_daty_row(nrow, ncol, temp, fh_binary_daty, array_temp);// read distributed z matrix row by row
			MPI_In_uox_mox(ncol, temp, fh_binary_daty_row, array_temp);

			for (int k = 0; k < MM; k++) {
				for (int t = 0; t < i_mxl; t++) {
					d_iy[counter][v_mxl[t] - 1] = array_temp[v_mxl[t] - 1];
				}
				counter++;
			}
		}
		//if (mynode == 1 && i == 3) {
		//	cout << "d_iy after at i=" << i << " at node " << mynode << endl;
		//	for (int p = 0; p < MM*i_size_v_cn_z_i; p++) {
		//		for (int m = 0; m < ncol; m++) {
		//			cout << setw(20) << d_iy[p][m];
		//		}
		//		cout << endl;
		//	}
		//}

		//------------
		//updated wij
		//------------

		//Note that fwij match v_obsg
		//fwij_selected should match selected donors

		for (int t = 0; t < i_size_v_cn_z_i*MM; t++) fwij_selected.push_back( fwij[i_SN[t] - 1] );

		if (fwij_selected.size() != (i_size_v_cn_z_i*MM)) cout << "ERROR!!!! Wrong fractional weights of selecetd donors!!!!" << endl;

		//if (mynode == 1) {
		//	cout << "fwij_selected at mox " << i << endl;
		//	for (int t = 0; t < fwij_selected.size(); t++) {
		//		cout << "fwij_selected[" << t << "]: " << fwij_selected[t] << endl;
		//	}
		//}

		for (int j = 0; j<i_size_v_cn_z_i*MM; j++) wij[j] = 1.0 / MM;


		//Local deallocation
		delete[] d_cs_temp;
		delete[] d_Ui;
		delete[] i_rmg;
		delete[] i_rM;
		delete[] d_SR;
	}			
	
	//------------------------------
	//------------------------------
	//make distributed return matrix fmat_FHDI_temp[][] and simp_fmat_FHDI_temp[][]
	//------------------------------
	//------------------------------
	//column 1: global id[]
	//column 2: w
	//column 3: wij
	//column 4: fwij
	//column 5~ : imputed d_iy

	int counter4 = 0;

	for (int j = 0; j<i_size_v_cn_z_i; j++)
	{
		for (int k = 0; k<MM; k++) //repeat each id MM times
		{
			//if (mynode == 1) cout << " Inside counter_fmat is " << counter_fmat << " at i= " << i << " at node " << mynode << endl;

			fmat_FHDI_temp[counter_fmat][0] = id[v_cn_z_i[j] - 1];
			fmat_FHDI_temp[counter_fmat][1] = w[v_cn_z_i[j] - 1];
			fmat_FHDI_temp[counter_fmat][2] = wij[counter4];
			fmat_FHDI_temp[counter_fmat][3] = fwij_selected[counter4];

			simp_fmat_FHDI_temp[counter_fmat][0] = id[v_cn_z_i[j] - 1];
			simp_fmat_FHDI_temp[counter_fmat][1] = w[v_cn_z_i[j] - 1];
			simp_fmat_FHDI_temp[counter_fmat][2] = wij[counter4];
			simp_fmat_FHDI_temp[counter_fmat][3] = fwij_selected[counter4];

			for (int h = 0; h < ncol; h++) {
				fmat_FHDI_temp[counter_fmat][h + 4] = d_iy[counter4][h];
			}

			counter4++;
			counter_fmat++;

		}
	}



	//success = MPI_File_close(&fh_daty);
	//if (success != MPI_SUCCESS) cout << "MPI I/O of z matrix fail to close the file!" << endl;
	

	//----------
	//local deallocation for FHDI parts
	//----------

	delete[] i_SN; 
	Del_dMatrix(d_iy, MM*i_size_v_cn_z_i, ncol);
	delete[] wij; 
	delete[] array_temp;		
		

	return; 
}