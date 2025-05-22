#include "matrix_utility_FHDI.cc"
#include <vector>
#include <string>
#include <mpi.h>
//#include "binary.cpp"
//below are defined in Variance_Est_FEFI_Extension_cpp.cc

//void RepWeight(const int n, double** d_rw) 

//void Rep_CellP(double** d_cx, const int nrow, const int ncol, double** d_rw, int*  id, 
//			   List_FHDI        &List_rst_prob,
//			   List_string_FHDI &List_rst_name,
//			   std::vector<std::string> &s_ncx)


void Variance_Est_FHDI_Extension_cpp(double** y, double** z, const int nrow, const int ncol,
	RepWeight_FHDI &d_rw, double* w, int* id,
	rbind_FHDI  &rbind_ipmat_FEFI,
	rbind_FHDI  &rbind_Resp_FEFI,
	rbind_FHDI  &rbind_irmat_FEFI,
	rbind_FHDI  &rbind_ipmat_FHDI,
	rbind_FHDI  &rbind_Resp_FHDI,
	rbind_FHDI  &rbind_irmat_FHDI,
	rbind_FHDI  &rbind_uox,
	rbind_FHDI  &rbind_mox,
	List_FHDI 	&List_ord,
	List_FHDI 	&List_ocsg,
	//List_FHDI   &List_rst_prob,
	//List_string_FHDI &List_rst_name,
	std::vector<int> &List_rst_prob_size,
	std::vector<std::string> &s_ncx,
	std::string s_M,
	double** y_bar_i_k_Summary, ofstream& TestOut)

	//Description----------------------
	//estimate variance for FHDI using Jackknife method 
	//  Algorithm: 
	//
	// original R code: Dr. Im, J. and Dr. Kim, J. 
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: June 22, 2020
	//
	//IN   : double y(nrow, ncol)= original data matrix with missing cells 
	//IN   : double z(nrow, ncol)= categorized matrix of y. 0 for missing cell
	//IN   : double d_rw[nrow, nrow] = replicate weights 
	//IN   : double w(nrow) = sampling weight (default = 1.0)
	//IN   : int    id(nrow) = id number of each row (default = 1 to nrow)
	//FEFI --------returns----------------- FEFI //
	//IN   : rbind_FHDI  rbind_ipmat_FEFI(4+ncol); //column size is 4+ncol
	//IN   : rbind_FHDI  rbind_Resp_FEFI(ncol+1);  //separate response matrix  
	//IN   : rbind_FHDI  rbind_irmat_FEFI(5+ncol); //column size is 5+ncol
	//FHDI --------returns----------------- FHDI //
	//IN   : rbind_FHDI  rbind_ipmat_FHDI(4+ncol); //column size is 4+ncol
	//IN   : rbind_FHDI  rbind_Resp_FHDI(ncol+1);  //separate response matrix  
	//IN   : rbind_FHDI  rbind_irmat_FHDI(5+ncol); //column size is 5+ncol
	//other matrices
	//IN   : rbind_FHDI  rbind_uox(ncol); //observed unique categorized matrix 
	//IN   : rbind_FHDI  rbind_mox(ncol); //missing  unique categorized matrix
	//Note: below Lists contain meaningful items up to i_count_mox rows  
	//IN   : List_FHDI 	List_ord(nrow); //order records used for variance estimation
	//IN   : List_FHDI 	List_ocsg(nrow); //order records used for variance estimation
	//IN   : std::string s_M = "FEFI" = fully efficient fractional imputation
	//						   "FHDI" = fractional hot deck imputation
	//OUT  : double** wmat = New_dMatrix(nrow_dat2_FHDI, L=nrow); //nrow_dat2_FHDI = rows of w1
	//
	//    
	//Data Structure Note
	//----------------------
	//dat1 in R version:  id   w  y_matrix  z_matrix
	//     in C++ ver  :  id   w  y         z
	//----------------------
	//dat2 in R version:  id  FID  WGT  FWGT  imputed matrix     Response(0/1)
	//     in C++ ver  :  ----- rbind_ipmat_FEFI------------     rbind_Resp_FEFI 
	//                 :  ----- rbind_ipmat_FHDI------------     rbind_Resp_FHDI
	//----------------------
	//
	//ipmat  = final imputation results
	//     	col1: ID 	= unit index
	//		col2: FID 	= ID of fractionally imputed value
	// 		col3: WGT 	= weight 
	//		col4: FWGT	= Frational weight
	//		col5: Variables 
	//		col6: Responses (separately in  rbind_Resp_...)
	//------------------------
	//
	//irmat  = imputation results related to the categorized matrix 
	//     	col1: ID 	= unit index
	//		col2: FID 	= ID of fractionally imputed value
	//		col3: OID	= original rank of the imputed value
	//		col4: ORDER = SN(selected donor)
	//		col5: FEFIW	= Fefi weights 
	//		col6: CELL	= cells 
	//----------------------
{
	//testout
	//RPrint("=========== Begin Variance Estimation of FHDI ================");

	//below is defined by User 
	/*	//-------------
	//sample weight (default is 1)
	//id array (default is row number)
	//-------------
	double* w = new double[nrow];
	int* id   = new int[nrow];
	for(int i=0; i<nrow; i++)
	{
	w[i] = 1.0;
	id[i] = i+1; //ACTUAL id
	}
	*/
	//===========================
	// MPI variables
	//===========================
	double startup = MPI_Wtime();
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//----------------------------
	//Basic constants declaration
	//----------------------------
	/*
	const int n  = nrow;
	const int nr = nrow;
	const int nc = ncol;
	const int nr1 = rbind_uox.size_row();
	const int nr2 = rbind_mox.size_row();
	const int nrow_uox 		 = rbind_uox.size_row();
	*/
	const int nrow_dat2_FEFI = rbind_ipmat_FEFI.size_row();
	const int nrow_dat2_FHDI = rbind_ipmat_FHDI.size_row();
	const int nrow_mox = rbind_mox.size_row();
	const int L = nrow; //size of d_rw 

	int numWorkPerProc = (int)floor(1.0*L / (1.0*totalnodes - 1));
	int numWorkLocalLast = L - numWorkPerProc * (totalnodes - 2);

	/*if (mynode == 0) {
		RPrint("FHDI The totalnodes is---------------------\n");
		RPrint(totalnodes);
		RPrint("---------------------------\n");

		RPrint("FHDI The numWorkPerProc is---------------------\n");
		RPrint(numWorkPerProc);
		RPrint("---------------------------\n");

		RPrint("FHDI The numWorkLocalLast is---------------------\n");
		RPrint(numWorkLocalLast);
		RPrint("---------------------------\n");

		RPrint("FHDI The L is---------------------\n");
		RPrint(L);
		RPrint("---------------------------\n");
	}*/
	//--------------------
	//get ready id table of FEFI
	//--------------------
	double* d_id_FHDI = new double[nrow_dat2_FHDI];
	for (int i = 0; i<nrow_dat2_FHDI; i++) d_id_FHDI[i] = rbind_ipmat_FHDI(i, 0); //id, 1st col 
	std::vector<double> v_table_name_id_FHDI; //same as "nimp" in R version
	std::vector<int>    v_table_count_id_FHDI;//same as "nimp" in R version 
	table_cpp_Yicheng(d_id_FHDI, nrow_dat2_FHDI, v_table_name_id_FHDI, v_table_count_id_FHDI);

	//std::vector<double> v_table_name_id_FHDI2; //same as "nimp" in R version
	//std::vector<int>    v_table_count_id_FHDI2;//same as "nimp" in R version 
	//table_cpp_Yicheng(d_id_FHDI, nrow_dat2_FHDI, v_table_name_id_FHDI2, v_table_count_id_FHDI2);

	//for (int i = 0;i < nrow_dat2_FHDI;i++) {
	//	cout<<"v_table_name_id_FHDI["<<i<<"]: "<< v_table_name_id_FHDI [i]<<endl;
	//}
	//for (int i = 0;i < nrow_dat2_FHDI;i++) {
	//	cout << "v_table_name_id_FHDI2[" << i << "]: " << v_table_name_id_FHDI2[i] << endl;
	//}
	//for (int i = 0;i < nrow_dat2_FHDI;i++) {
	//	cout << "v_table_count_id_FHDI[" << i << "]: " << v_table_count_id_FHDI[i] << endl;
	//}
	//for (int i = 0;i < nrow_dat2_FHDI;i++) {
	//	cout << "v_table_count_id_FHDI2[" << i << "]: " << v_table_count_id_FHDI2[i] << endl;
	//}
	

	//if (mynode == 1) {
	//	cout << "YYC Running time of startup11 at node " << mynode << " is " << MPI_Wtime() - startup << endl;
	//}
	double startup12 = MPI_Wtime();
	//---------------------
	//imputed real data matrix
	//---------------------
	int nrow_d_iy = nrow_dat2_FHDI; //default
	if (s_M == "FEFI") { nrow_d_iy = nrow_dat2_FEFI; }
	if (s_M == "FHDI") { nrow_d_iy = nrow_dat2_FHDI; }
	double** d_iy = New_dMatrix(nrow_d_iy, ncol); //imputed matrix of real values 
	for (int i = 0; i<ncol; i++)
	{
		for (int j = 0; j<nrow_d_iy; j++)
		{
			if (s_M == "FEFI")
			{
				d_iy[j][i] = rbind_ipmat_FEFI(j, 4 + i);
			} //col5~ncol contains imputed real values 
			if (s_M == "FHDI")
			{
				d_iy[j][i] = rbind_ipmat_FHDI(j, 4 + i);
			} //col5~ncol contains imputed real values 

		}
	}
	//if (mynode == 1) {
	//	cout << "YYC Running time of startup12 at node " << mynode << " is " << MPI_Wtime() - startup12 << endl;
	//}
	double startup13 = MPI_Wtime();
	//----------------------
	//categorized matrix
	//----------------------
	double** d_cx = New_dMatrix(nrow, ncol);
	Copy_dMatrix(z, nrow, ncol, d_cx);

	//---------------------
	//ocg, observed donors for each missing pattern. 
	//---------------------
	//the same as List_ocsg[nrow_mox]
	//---------------------
	int* i_locg = new int[nrow_mox]; //length of donor rows for each missing pattern 
	for (int i = 0; i<nrow_mox; i++)
	{
		int i_temp = 0;
		List_ocsg.get_a_row_size(i, i_temp);
		i_locg[i] = i_temp; //meaning how many rows used as the donor for ith missing pattern 
	}
	//if (mynode == 1) {
	//	cout << "YYC Running time of startup13 at node " << mynode << " is " << MPI_Wtime() - startup13 << endl;
	//}

	//if (mynode == 1) {
	//	cout << "YYC Running time of startup1 at node " << mynode << " is " << MPI_Wtime() - startup << endl;
	//}
	double startup2 = MPI_Wtime();
	//testout
	/*
	RPrint("==== in Variance_Est_Extension_cpp ========");
	RPrint("id: "); RPrint(id, n);
	RPrint("n : "); RPrint(n);
	RPrint("nr: "); RPrint(nr);
	RPrint("nc: "); RPrint(nc);
	RPrint("--------dat1: id (above)and w,  y and z ");
	RPrint("w: "); RPrint(w, n);
	RPrint("y: "); RPrint(y, nrow, ncol);
	RPrint("z: "); RPrint(z, nrow, ncol);
	RPrint("--------dat2: ipmat_FEFI     Resp_FEFI ------- ");
	rbind_ipmat_FEFI.print_rbind_FHDI();
	rbind_Resp_FEFI.print_rbind_FHDI();
	RPrint("--------dat2: ipmat_FHDI     Resp_FHDI ------- ");
	rbind_ipmat_FHDI.print_rbind_FHDI();
	rbind_Resp_FHDI.print_rbind_FHDI();
	RPrint("iy : (imputed real data)"); RPrint(d_iy, nrow_d_iy, ncol);
	RPrint("cx : (categorized matrix)"); RPrint(d_cx, nrow, ncol);
	RPrint("ocg: ");  List_ocsg.print_List_FHDI();
	RPrint("locg: "); RPrint(i_locg, nrow_mox);
	RPrint("nr1: "); RPrint(nr1);
	RPrint("nr2: "); RPrint(nr2);
	*/

	//------------------------
	//cell probability using replicate weight
	//------------------------
	//List_FHDI         List_rst_prob(nrow); //only i_nc rows are meaningful
	//List_string_FHDI  List_rst_name(nrow); //only i_nc rows are meaningful
	//std::vector<std::string> s_ncx;

	//Rep_CellP(d_cx, nrow, ncol, d_rw, id, 
	// List_rst_prob,
	// List_rst_name,
	// s_ncx);
	//const int i_nc = (int)s_ncx.size(); //meaningful rows of List's 

	//testout
	//RPrint("List_rst_name"); List_rst_name.print_List_string_FHDI(); //njp in R ver
	//RPrint("List_rst_prob"); List_rst_prob.print_List_FHDI(); 		 //Tcellp in R version

	//--------------------
	//put 1 into fully observed rows
	//--------------------
	double* d_rr0 = new double[nrow];
	for (int i = 0; i<nrow; i++)
	{
		double d_prod = 1.0;
		for (int j = 0; j<ncol; j++)
		{
			d_prod = d_prod*d_cx[i][j];
		}
		//-----
		//0: at least one missing;  1: all observed
		//-----
		d_rr0[i] = d_prod;
		if (fabs(d_prod) >0.0) d_rr0[i] = 1;
	}

	//--------------
	//calculate w1 = sampling weight
	//--------------
	//std::string cn[nrow]; 
	std::string *cn = new std::string[nrow];
	Trans(d_cx, nrow, ncol, cn);
	double* d_w1 = new double[nrow_dat2_FHDI];
	for (int i = 0; i<nrow_dat2_FHDI; i++)
		d_w1[i] = rbind_ipmat_FHDI(i, 2); //3rd column contains WGT

										  //testout
										  //RPrint("rr0: "); RPrint(d_rr0, nrow);
										  //RPrint("w1: "); RPrint(d_w1, nrow_dat2_FHDI);
	//if (mynode == 1) {
	//	cout << "YYC Running time of startup2 at node " << mynode << " is " << MPI_Wtime() - startup2 << endl;
	//}
	//------------------
	//make Covariance matrix for the subsequent replication process
	//------------------
	double startup3 = MPI_Wtime();
	//int* i_lloc;
	//std::vector<int> v_mox_0;
	//double** d_dy;
	//double** V_var; //covariance matrix of dy
	//List_FHDI List_V(nrow_mox); //storage of covariance matrix
	//							//Note: store cov mat by "row-first" rule 
	////if (mynode == 0) cout << "nrow_mox: " << nrow_mox << endl;
	//for (int i = 0; i<nrow_mox; i++)
	//{
	//	//if (mynode == 0) {
	//	//	cout << "FHDI MOX at " << i << endl;
	//	//}
	//	//cout << "Mynode:" << mynode<<" nrow_mox["<<i<<"]"<<endl;
	//	//---------------------------
	//	//how many donor rows for the ith missing pattern
	//	//------------------------
	//	int i_size_lloc = i_locg[i];
	//	//-----------
	//	//when zero size continue to next iteration
	//	//-----------
	//	if (i_size_lloc <= 0) { continue; }
	//	i_lloc = new int[i_size_lloc]; //vector of actual donor row numbers 
	//	for (int j = 0; j< i_size_lloc; j++) i_lloc[j] = (int)List_ocsg(i, j); //ith row, jth entity

	//																		   //-----
	//																		   //find missing columns in current missing row
	//																		   //------
	//	v_mox_0.clear();
	//	for (int j = 0; j<ncol; j++)
	//	{
	//		if (rbind_mox(i, j) == 0.0) v_mox_0.push_back(j + 1); //ACTUAL zero column id 
	//	}
	//	const int i_size_v_mox_0 = (int)v_mox_0.size();

	//	//----------------------------------
	//	//extract matrix of missing patterns
	//	//----------------------------------
	//	//-----------
	//	//when zero size continue to next iteration
	//	//-----------
	//	if (i_size_lloc <= 0) { continue; }
	//	if (i_size_v_mox_0 <= 0) { continue; }

	//	cout << "i_size_lloc and i_size_v_mox_0 are " << i_size_lloc << ",  " << i_size_v_mox_0 << " at mox " << i << " at mynode " << mynode << endl;
	//	d_dy = New_dMatrix(i_size_lloc, i_size_v_mox_0);
	//	V_var = New_dMatrix(i_size_v_mox_0, i_size_v_mox_0); //column-wise covariance 
	//	for (int j = 0; j< i_size_lloc; j++) //LOOP for donor rows 
	//	{
	//		for (int k = 0; k<i_size_v_mox_0; k++) //LOOP for missing columns 
	//		{
	//			d_dy[j][k] = y[i_lloc[j] - 1][v_mox_0[k] - 1]; //-1 for actual location 
	//		}
	//	}
	//	//----------
	//	//"Estimated covariance" of d_dy by column-to-column method
	//	//----------
	//	cov_FHDI(d_dy, i_size_lloc, i_size_v_mox_0, V_var);


	//	//----------
	//	//store the covariance matrix 
	//	//row-first rule
	//	//----------
	//	//const int i_List_V_col = i_size_v_mox_0*i_size_v_mox_0;  
	//	//double* d_V_temp = new double[i_List_V_col];
	//	//for(int j=0; j<i_size_v_mox_0; j++) 
	//	//{
	//	//	for(int k=0; k<i_size_v_mox_0; k++)
	//	//		d_V_temp[j*i_size_v_mox_0 + k]= V_var[k][j]; 
	//	//}

	//	//List_V.put_block(i, i_List_V_col, d_V_temp); //ith covariance matrix 
	//	List_V.put_block(i, i_size_v_mox_0, i_size_v_mox_0, V_var); //direct matrix saving

	//																//---
	//																//local deallocation
	//																//---
	//	delete[] i_lloc;
	//	Del_dMatrix(d_dy, i_size_lloc, i_size_v_mox_0);
	//	Del_dMatrix(V_var, i_size_v_mox_0, i_size_v_mox_0);
	//	//delete[] d_V_temp; 

	//}
	//if (mynode == 1) {
	//	cout << "YYC Running time of startup3 at node " << mynode << " is " << MPI_Wtime() - startup3 << endl;
	//}
	double startup4 = MPI_Wtime();
	//testout
	//RPrint("  List_V");
	//List_V.print_List_FHDI(); 

	//----------------
	//wmat: Replication Weights
	//----------------
	//double** wmat = New_dMatrix(nrow_dat2_FHDI, L); //nrow_dat2_FHDI = rows of w1

	//------------------------------
	//------------------------------
	//MAIn loop for L replications
	//------------------------------
	//------------------------------
	double* rw0 = new double[nrow];

	int i_sum_Rw = 0;
	for (int i = 0; i<nrow; i++) i_sum_Rw += v_table_count_id_FHDI[i];
	//cout << "i_sum_Rw is " << i_sum_Rw << " at mynode " << mynode << endl;

	double* Rw = new double[i_sum_Rw];

	double* wijk = new double[nrow_dat2_FHDI]; //FWGT from ipmat 

	int i_nrow_imputation = nrow;
	if (s_M == "FEFI") i_nrow_imputation = nrow_dat2_FEFI;
	if (s_M == "FHDI") i_nrow_imputation = nrow_dat2_FHDI;
	//cout << "i_nrow_imputation first: " << i_nrow_imputation << endl;
	//if (mynode == 1) {
	//	cout << "YYC Running time of startup4 at node " << mynode << " is " << MPI_Wtime() - startup4 << endl;
	//}
	//================================
	// parallelization of main loop
	//================================
	double* wmatLocal = new double[i_nrow_imputation];

	int L_temp = 0;
	if (mynode != (totalnodes - 1)) L_temp = numWorkPerProc;
	if (mynode == (totalnodes - 1)) L_temp = numWorkLocalLast;

	//cout << "Mynode: " << mynode << ", L_temp: " << L_temp << endl;
	double** y_bar_i_k = New_dMatrix(L_temp, ncol);
	Fill_dMatrix(y_bar_i_k, L_temp, ncol, 0.0);

	int i_loc = 0;
	double* yi = new double[ncol];
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
			TestOut << "y_bar_i_k boundary ERROR!!!" << endl;
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
			TestOut << "y_bar_i_k boundary ERROR!!!" << endl;
			return;
		}
	}
	//if (mynode == 1) {
	//	cout << "FHDI Mynode " << mynode <<",  range are " << startpoint << "<= x < " << endpoint << endl;
	//}

	//if (mynode == 1) {
	//	cout << "YYC Running time of startup3 at node " << mynode << " is " << MPI_Wtime() - startup3 << endl;
	//	cout << "YYC Running time of startup at node " << mynode << " is " << MPI_Wtime() - startup << endl;
	//}


	int* rbind_ipmat = new int[nrow_dat2_FHDI];
	for (int kk=0; kk < nrow_dat2_FHDI;kk++) {
		rbind_ipmat[kk] = rbind_ipmat_FHDI(kk, 0);
		//if (mynode == 1) {
		//	cout << "rbind[" << kk<< "]: " << rbind_ipmat[kk] << endl;
		//}
		
	}
	//if (mynode == 1) {
	//	cout<<"left brefore: "<<left<<" and right before: "<<right<<endl;
	//}
	//binary(left, right, target, nrow_dat2_FHDI, rbind_ipmat);
	//if (mynode == 1) {
	//	cout << "left after: " << left << " and right after: " << right << endl;
	//}

	double L_begin = MPI_Wtime();
	//int yicheng_counter = 0;
	for (int l = startpoint; l<endpoint; l = l + 1)
	{
		double L_Rep_begin = MPI_Wtime();
		//cout << "Successfully get started of L-replication of FHDI at mynode" << mynode << endl;
		//if (mynode == 1) {
		//	cout << "FHDI L_rep at L =  " << l << " at mynode " << mynode << endl;
		//}
		//cout<<"Mynode: "<<mynode<<", l: "<<l<<endl;
		//-------
		//replicate weight from lth column
		//-------
		for (int i = 0; i<nrow; i++) rw0[i] = d_rw(i, l); //l_th column 
		int i_sum = 0;
		for (int i = 0; i<nrow; i++)
		{
			for (int j = 0; j<v_table_count_id_FHDI[i]; j++) Rw[i_sum++] = rw0[i];
		}

		//---------
		//FWGT of ipmat
		//----------
		for (int i = 0; i<nrow_dat2_FHDI; i++)
			wijk[i] = rbind_ipmat_FHDI(i, 3); //4th column is FWGT 

											  //----------
											  //joint probability associated with current string
											  //-----------
		std::string cn_current = cn[l]; //lth string 
		std::vector<int> v_ncx_cn;
		which(s_ncx, cn_current, v_ncx_cn); //actual location 
											//const int i_size_v_ncx_cn = (int)v_ncx_cn.size(); //MUST BE "1"

		int i_size_cellp = 0;
		//List_rst_prob.get_a_row_size(v_ncx_cn[0] - 1, i_size_cellp); //get a size of the row in the list 
		i_size_cellp = List_rst_prob_size[v_ncx_cn[0] - 1];
																	 //-----------
																	 //when zero size continue to next iteration
																	 //-----------
		if (i_size_cellp <= 0) { continue; }

		//cout << "i_size_cellp is "<< i_size_cellp <<" at L="<<l<<endl;
		//double* d_cellp = new double[i_size_cellp];
		//List_rst_prob.get_block(v_ncx_cn[0] - 1, d_cellp); //-1 for actual row location 

		//testout
		//RPrint(" ======== in Main Loop l+1: "); RPrint(l+1); 
		//RPrint("Rw"); RPrint(Rw, i_sum_Rw); 
		//RPrint("wijk"); RPrint(wijk, nrow_dat2_FHDI); 
		//RPrint("d_cellp"); RPrint(d_cellp, i_size_cellp); 

		//----------------------------------------
		//1. if the deleted is missing unit, no further action is taken
		//2. if the deleted is observed unit, then the fractional weights are re-computed 
		//----------------------------------------
		int* idd = new int[nrow_mox]; //location of the deleted donor in ocg 
		Fill_iVector(idd, nrow_mox, 0);


		//if ((l==1)||(l==2)) {
		//		cout << "L_Rep_1 running time within one loop L " << l << " is " << MPI_Wtime() - L_Rep_begin << endl;
		//}


		double L_Rep_begin_2 = MPI_Wtime();
		if (fabs(d_rr0[l]) > 0)
		{
			//---------------------
			//locations of the deleted unit in observed list
			//---------------------
			std::vector<int> v_lg; //Actual locations 
			v_lg.clear();
			for (int j = 0; j<nrow_mox; j++) //all missing patterns
			{
				for (int k = 0; k<i_locg[j]; k++) //donor rows for the jth missing pattern
				{
					int i_temp_lg = (int)List_ocsg(j, k);
					if (i_temp_lg == (l + 1)) //+1 for actual location  
					{
						v_lg.push_back(j + 1); //actual id of jth missing pattern  
						idd[j] = k + 1; //actual number of donor rows for jth missing pattern 
						break;
					}
				}
			}
			const int nlg = (int)v_lg.size();

			//if ((l == 1) || (l == 2)) {
			//		cout << "L_Rep_2 running time within one loop L " << l << " is " << MPI_Wtime() - L_Rep_begin_2 << endl;
			//}

			//testout
			//RPrint(" in condition rr0[l]!=0 at l+1 ="); RPrint(l+1);
			//RPrint("idd:"); RPrint(idd, nrow_mox);
			//RPrint("lg :"); RPrint(v_lg);
			//RPrint("nlg:"); RPrint(nlg); 

			//--------------------------
			//Adjust fractional weights for all units in lg
			//--------------------------
			double L_Rep_begin_3 = MPI_Wtime();
			if (nlg>0)
			{
				for (int j = 0; j<nlg; j++)
				{
					//if (mynode == 1) {
					//	cout<<"Nlg is "<<j<<" at L "<<l<<endl;
					//}
					int i_row_lg = v_lg[j] - 1; // row number [0,...) 
					double* d_1_mox = new double[ncol];
					for (int k = 0; k<ncol; k++) d_1_mox[k] = rbind_mox(i_row_lg, k);

					//---
					//actual col number of missing cell in current missing row
					//---
					std::vector<int> v_rloc; v_rloc.clear();
					for (int k = 0; k<ncol; k++)
					{
						if (d_1_mox[k] == 0.0) { v_rloc.push_back(k + 1); } //actual col
					}

					//---
					//number of missing columns in current missing row
					//---
					const int nrloc = (int)v_rloc.size();
					std::string cng;
					Trans1(d_1_mox, ncol, cng);

					//-------
					//location of cn which has cng
					//-------
					std::vector<int> v_mlog; v_mlog.clear();
					which(cn, nrow, cng, v_mlog);
					const int nmlog = (int)v_mlog.size();

					//testout
					//RPrint("rloc: "); RPrint(v_rloc);
					//RPrint("mlog: "); RPrint(v_mlog);


					//------------------------
					//------------------------
					//FHDI
					//------------------------
					//------------------------
					std::vector<int> v_elog; v_elog.clear();
					if (s_M.compare("FHDI") == 0) //0=equal 
					{
						//-----
						//find locations of mlog in dat2$ID
						//v_mlog contains the row numbers that have the same string as
						//current missing row 
						//nmlog = n(v_mlog)
						//-----
						v_elog.clear();
						for (int k1 = 0; k1<nmlog; k1++) //loop for mlog
						{
							int left = 0;
							int right = 0;
							int i_temp1 = id[v_mlog[k1] - 1]; //dat1$ID in R version  
							//for (int k2 = 0; k2<nrow_dat2_FHDI; k2++)
							//{
							//	int i_temp2 = rbind_ipmat_FHDI(k2, 0); //1st col is dat2$ID
							//	if (i_temp1 == i_temp2)
							//	{
							//		v_elog.push_back(k2 + 1); //actual location 
							//	}
							//}
							binary(left, right, i_temp1, nrow_dat2_FHDI, rbind_ipmat);
							if (left != (-1) && right != (-1)) {
								for (int k2 = left;k2 < (right + 1);k2++) {
									v_elog.push_back(k2 + 1);
								}
							}
						}
						const int i_size_v_elog = (int)v_elog.size();
						
						//------------------------------------
						//set of donors for missing columns
						//------------------------------------
						//-----------
						//when zero size continue to next iteration
						//-----------
						if (i_size_v_elog <= 0) { continue; }
						if (nrloc <= 0) { continue; }
						double** dy_FHDI = New_dMatrix(i_size_v_elog, nrloc);
						for (int k1 = 0; k1<nrloc; k1++)
						{
							for (int k2 = 0; k2<i_size_v_elog; k2++)
							{
								dy_FHDI[k2][k1] = d_iy[v_elog[k2] - 1][v_rloc[k1] - 1];
							}
						}
						//testout
						//RPrint(" = FHDI l+1: "); RPrint(l+1); 
						//RPrint(" elog: "); RPrint(v_elog); 
						//RPrint(" dy  : "); RPrint(dy_FHDI,i_size_v_elog, nrloc); 

						//---------------------
						// nrloc >= 1: number of missing columns in current missing row
						//---------------------
						if (nrloc >= 1)
						{
							//----
							//make dk matrix
							//filled with l_th original data 
							//at missing column locations  
							//Note: this is the Jackknifed row
							//----------------------------------------
							//-----------
							//when zero size continue to next iteration
							//-----------
							if (i_size_v_elog <= 0) { continue; }
							if (nrloc <= 0) { continue; }


							double** dk = New_dMatrix(i_size_v_elog, nrloc);
							for (int k1 = 0; k1<nrloc; k1++)
							{
								double d_temp_dk = y[l][v_rloc[k1] - 1]; //-1 for actual location  
								for (int k2 = 0; k2<i_size_v_elog; k2++)
								{
									dk[k2][k1] = d_temp_dk;
								}
							}
							//testout
							//RPrint("i_size_v_elog:"); RPrint(i_size_v_elog); 
							//RPrint("nrloc        :"); RPrint(nrloc); 
							//RPrint("dk:"); RPrint(dk, i_size_v_elog, nrloc); 

							//---------
							//difference between   dy         and dk 
							//i.e., difference b/w donor rows and jackknifed row  
							//-------------------------------------------
							//-----------
							//when zero size continue to next iteration
							//-----------
							if (i_size_v_elog <= 0) { continue; }
							if (nrloc <= 0) { continue; }
							double** diff = New_dMatrix(i_size_v_elog, nrloc);
							for (int k1 = 0; k1<nrloc; k1++)
							{
								for (int k2 = 0; k2<i_size_v_elog; k2++)
								{
									diff[k2][k1] = dk[k2][k1] - dy_FHDI[k2][k1];
								}
							}


							//----------
							//get l_th covariance matrix
							//--------------------------------------------
							//-----------
							//when zero size continue to next iteration
							//-----------
							//if (mynode == 1) {
							//cout << "nrloc_2 is " << nrloc << " at L " << l << " at mynode " << mynode <<" at nlg "<<j<< endl;
							//}
							if (nrloc <= 0) { continue; }
							double** V_var_l = New_dMatrix(nrloc, nrloc);
							int i_loc_lg_j = v_lg[j] - 1; //-1 for actual location 
							//List_V.get_block(i_loc_lg_j, nrloc, nrloc, V_var_l); //matrix read by row-first rule
                            //-----------------------------------------------------------------------------------------------------
							//Updated on 6/22/2020. List_V adds covariance matrix of each mox. It leads to out of memory issue if p is very large
							//To avoid excessive requirement of memory, covariance matrix of each mox will not be generated in advance. 
							int i_size_lloc = i_locg[i_loc_lg_j];

							int* i_lloc = new int[i_size_lloc]; //vector of actual donor row numbers 
							for (int a1 = 0; a1 < i_size_lloc; a1++) i_lloc[a1] = (int)List_ocsg(i_loc_lg_j, a1); //ith row, jth entity

							std::vector<int> v_mox_0;
							v_mox_0.clear();
							for (int a2 = 0; a2 < ncol; a2++)
							{
								if (rbind_mox(i_loc_lg_j, a2) == 0.0) v_mox_0.push_back(a2 + 1); //ACTUAL zero column id 
							}
							const int i_size_v_mox_0 = (int)v_mox_0.size();
						

							double** d_dy = New_dMatrix(i_size_lloc, i_size_v_mox_0);

							for (int k3 = 0; k3 < i_size_lloc; k3++) //LOOP for donor rows 
							{
								for (int k4 = 0; k4 < i_size_v_mox_0; k4++) //LOOP for missing columns 
								{
									//cout<<"row: "<< i_lloc[k3] - 1 <<", column: "<< v_mox_0[k4] - 1 <<endl;
									d_dy[k3][k4] = y[i_lloc[k3] - 1][v_mox_0[k4] - 1]; //-1 for actual location 
								}
							}
							//cout << "Debug i_loc_lg_j is " << i_loc_lg_j << " at j = " << j << endl;

							//----------
							//"Estimated covariance" of d_dy by column-to-column method
							//----------
							cov_FHDI(d_dy, i_size_lloc, i_size_v_mox_0, V_var_l);
							//cout<<"Debug i_size_lloc is "<< i_size_lloc <<" at j = "<<j<<endl;

							delete[] i_lloc;
							Del_dMatrix(d_dy, i_size_lloc, i_size_v_mox_0);
							//----------------------------------------------------------------------------------------------

							if (nrloc <= 0) { continue; }
							double** V_inv = New_dMatrix(nrloc, nrloc);

							bool b_success_V_inv = true; //false when abrupt exit due to zero diagonal
							if (nrloc > 1) b_success_V_inv = Inverse_dMatrix_FHDI(V_var_l, nrloc, V_inv);
							if (nrloc == 1) V_inv[0][0] = 1.0 / V_var_l[0][0];
							//if (mynode == 1) {
							//	cout << "FHDI inside L_replication 3 at L = "<< l << " at mynode "<<mynode<<" at nlg "<<j<<endl;
							//}
							if (!b_success_V_inv)
							{
								//testout
								/*cout<<"nrloc: "<<nrloc<<endl;
								cout<<"l: "<<l<<",  L :"<<j<<endl;
								cout<<"j: "<<j<<",  nlg :"<<nlg<<endl;
								cout<<"V_var_l[][]"<<endl;
								RPrint(V_var_l, nrloc, nrloc);
								for(int i_V_var = 0; i_V_var<nrloc; i_V_var++)
								{
								cout<<i_V_var<<":  "<<V_var_l[i_V_var][i_V_var]<<endl;
								}*/

								//----
								//below is an option to abort 
								//zero-diagonal Variance matrix
								//However, if donor vector was zero
								//vector, exceptional consideration
								//is needed
								//Feb 07, 2017
								//----
								//cout<<endl<<"Error! zero diagonal term in V_var_l"<<endl; return; 

								//----
								//special remedy for zero diagonal covariance matrix
								//Feb 7 2017
								//This is enough since we are interested in 
								//relative ordering to find minimum distance later
								//----
								for (int i_inv1 = 0; i_inv1<nrloc; i_inv1++)
								{
									double d_temp_inv = V_var_l[i_inv1][i_inv1];

									//for zero diagonal term: a big number 
									if (fabs(d_temp_inv) <= 1e-13)
									{
										V_inv[i_inv1][i_inv1] = 1E15;
									}

									//for non-zero diagonal term: simple inverse
									if (fabs(d_temp_inv) > 1e-13)
									{
										V_inv[i_inv1][i_inv1] = 1.0 / d_temp_inv;
									}
								}
							}

							//-----------
							//when zero size continue to next iteration
							//-----------
							if (nrloc <= 0) { continue; }
							if (i_size_v_elog <= 0) { continue; }
							double** diff_T = New_dMatrix(nrloc, i_size_v_elog);
							for (int k1 = 0; k1<nrloc; k1++)
							{
								for (int k2 = 0; k2<i_size_v_elog; k2++)
								{
									diff_T[k1][k2] = diff[k2][k1];
								}
							}
							//-----------
							//when zero size continue to next iteration
							//-----------
							if (i_size_v_elog <= 0) { continue; }


							double*diff_V_diffT = new double[i_size_v_elog];
						
							
							dMatrix_Mul_AtBA_Yicheng(diff_T, nrloc, i_size_v_elog,
								V_inv, diff_V_diffT);

							//-----------
							//score = diagonal terms
							//-----------
							double* d_score = new double[i_size_v_elog];
							for (int k1 = 0; k1<i_size_v_elog; k1++)
								d_score[k1] = diff_V_diffT[k1];


							//testout
							//RPrint(" = FHDI l+1: "); RPrint(l+1); 
							//RPrint("score        :"); RPrint(d_score, i_size_v_elog); 

							//------
							//MM calculation
							//------
							if (nmlog == 0)
							{
								TestOut << "Caution! zero mlog!" << endl; continue;
							}

							const int MM = i_size_v_elog / nmlog;
							const int ncol_imt = (int)floor(i_size_v_elog / MM);
			

							//-----------
							//when zero size continue to next iteration
							//-----------
							if (ncol_imt <= 0) { continue; }
							if (MM <= 0) { continue; }
							double** d_imt = New_dMatrix(MM, ncol_imt);
							int i_score = 0;
							for (int k1 = 0; k1<ncol_imt; k1++)
							{
								for (int k2 = 0; k2<MM; k2++) //row-first rule 
								{
									d_imt[k2][k1] = d_score[i_score++];
								}
							}

							//-----------------------------------------
							//extract weights
							//-----------------------------------------
							//-----------
							//when zero size continue to next iteration
							//-----------
							if (i_size_v_elog <= 0) { continue; }

							double* ewijk = new double[i_size_v_elog];
							double* fefiw = new double[i_size_v_elog];

							for (int k1 = 0; k1<i_size_v_elog; k1++)
							{
								ewijk[k1] = wijk[v_elog[k1] - 1]; //-1 for actual loc
																  //5th column is FEFIW of fhdi[[2]] in r version 
								fefiw[k1] = rbind_irmat_FHDI(v_elog[k1] - 1, 4); //5th column
							}

							//testout
							//RPrint("MM :"); RPrint(MM);
							//RPrint("ewijk :"); RPrint(ewijk, i_size_v_elog);
							//RPrint("fefiw :"); RPrint(fefiw, i_size_v_elog);
							//RPrint("d_imt :"); RPrint(d_imt, MM, ncol_imt);

							//-------------------
							//find eloc
							//-------------------
							std::vector<int> LMM; LMM.clear();
							for (int k1 = 0; k1<nmlog; k1++) LMM.push_back(k1*MM);
							//-----------
							//when zero size continue to next iteration
							//-----------
							if (nmlog <= 0) { continue; }
							int* i_eloc = new int[nmlog];
							//column-wise min location (Actual)
							for (int k1 = 0; k1<nmlog; k1++)
							{
								double d_temp = d_imt[0][k1]; //1st row  
								int    i_min = 1;
								for (int k2 = 1; k2<MM; k2++)
								{
									if ((d_temp - d_imt[k2][k1]) > 1e-3)
									{
										d_temp = d_imt[k2][k1];
										i_min = k2 + 1; //actual row location 
									}
								}
								//----
								//store the column-wise min location (Actual)
								//----
								i_eloc[k1] = i_min + LMM[k1];
								//testout
								//RPrint("LMM[k1] :"); RPrint(LMM[k1]);
								//RPrint("i_eloc[k1] :"); RPrint(i_eloc[k1]);
								//RPrint("i_min :"); RPrint(i_min);

							}

							//-------------
							//maximum difference (>0) between ewijk and fefiw
							//-------------------------------------
							//-----------
							//when zero size continue to next iteration
							//-----------
							if (nmlog <= 0) { continue; }
							//if (mynode == 1) {
							//	cout << "FHDI inside L_replication 6 at L = " << l << " at mynode " << mynode << " at nlg "<<j<<endl;
							//}
							double* d_maxew = new double[nmlog]; //nmlog = length of i_eloc
							double* d_maxval = new double[nmlog]; //nmlog = length of i_eloc

							for (int k1 = 0; k1<nmlog; k1++)
							{
								double d_temp = ewijk[i_eloc[k1] - 1] - fefiw[i_eloc[k1] - 1]; //-1 for actual loc
								d_maxew[k1] = 0.0;
								if (d_temp > 0.0) d_maxew[k1] = d_temp;

								d_maxval[k1] = ewijk[i_eloc[k1] - 1] - d_maxew[k1];
							}
							//testout
							//RPrint("d_maxew :"); RPrint(d_maxew, nmlog);
							//RPrint("d_maxval :"); RPrint(d_maxval, nmlog);

							//------
							//ewijk update
							//------
							for (int k1 = 0; k1<nmlog; k1++)
							{
								ewijk[i_eloc[k1] - 1] = d_maxew[k1];
							}

							//-------------
							//extend maxval array
							//-------------
							if (MM == 1) { TestOut << "Error! MM is 1 in Var FHDI" << endl; return; }
							//-----------
							//when zero size continue to next iteration
							//-----------
							if (nmlog <= 0) { continue; }
							

							double* d_maxval_extended = new double[nmlog*(MM - 1)];
							int i_maxval = 0;
							for (int k1 = 0; k1<nmlog; k1++)
							{
								double d_temp = d_maxval[k1] / (MM - 1);
								for (int k2 = 0; k2<(MM - 1); k2++)
								{
									d_maxval_extended[i_maxval++] = d_temp;
								}
							}

							//---------------
							//update ewijk with extended maxval
							//---------------
							//exclude eloc locations
							//-----------
							//when zero size continue to next iteration
							//-----------
							if ((i_size_v_elog - nmlog) <= 0) { continue; }
							int* i_without_eloc = new int[i_size_v_elog - nmlog];
							int i_eloc_temp = 0;
							for (int k1 = 0; k1<i_size_v_elog; k1++)
							{
								bool b_same = false;
								for (int k2 = 0; k2<nmlog; k2++)
								{
									if (i_eloc[k2] - 1 == k1)
									{
										b_same = true;
										break;
									}
								}
								if (!b_same) { i_without_eloc[i_eloc_temp++] = k1; }
							}
							//store values
							
							

							for (int k1 = 0; k1<(i_size_v_elog - nmlog); k1++)
							{
								int i_loc_ew = i_without_eloc[k1];
								ewijk[i_loc_ew] = ewijk[i_loc_ew] + d_maxval_extended[k1];
							}
							//testout
							//RPrint("ewijk :"); RPrint(ewijk, i_size_v_elog);


							//----------------
							//final update wijk with ewijk
							//----------------
							for (int k1 = 0; k1<i_size_v_elog; k1++)
							{
								wijk[v_elog[k1] - 1] = ewijk[k1];

								//testout
								//RPrint("v_elog[k1]-1 :"); RPrint(v_elog[k1]-1);
								//RPrint("wijk[..] :"); RPrint(wijk[v_elog[k1]-1]);

							}

							//if (mynode == 1) {
							//	cout << "FHDI inside L_replication finish at L = " << l << " at mynode " << mynode << " at nlg " << j << endl;
							//}

							//---
							//local deallocation
							//---
							Del_dMatrix(dk, i_size_v_elog, nrloc);
							Del_dMatrix(diff, i_size_v_elog, nrloc);
							Del_dMatrix(V_var_l, nrloc, nrloc);
							Del_dMatrix(V_inv, nrloc, nrloc);
							Del_dMatrix(diff_T, nrloc, i_size_v_elog);
							//Del_dMatrix(diff_V_diffT, i_size_v_elog, i_size_v_elog);
							delete[] diff_V_diffT;
							delete[] d_score;
							Del_dMatrix(d_imt, MM, ncol_imt);
							delete[] ewijk;
							delete[] fefiw;
							delete[] i_eloc;
							delete[] d_maxew;
							delete[] d_maxval;
							delete[] d_maxval_extended;
							delete[] i_without_eloc;

							//NOTE: dy_FHDI has different order from R version
							//as of Nov 27, 2016. It looks fine overall, but 
							//need to check later !!!!
						}


						//----------
						//local deallocation 
						//----------
						Del_dMatrix(dy_FHDI, i_size_v_elog, nrloc);
					}

					//local deallocation 
					delete[] d_1_mox;
				}
			}


			//if ((l == 1) || (l == 2)) {
			//	cout << "L_Rep_3 running time within one loop L " << l << " is " << MPI_Wtime() - L_Rep_begin_3 << endl;
			//}


		}

		//if ((l == 1) || (l == 2)) {
		//		cout << "L_Rep_fabs running time within one loop L " << l << " is " << MPI_Wtime() - L_Rep_begin_2 << endl;
		//}

		//-------------------------------
		//store the updated weights
		//-------------------------------
		double L_Rep_begin_41 = MPI_Wtime();
		int i_nrow_imputation = nrow;
		if (s_M == "FEFI") i_nrow_imputation = nrow_dat2_FEFI;
		if (s_M == "FHDI") i_nrow_imputation = nrow_dat2_FHDI;

		for (int k1 = 0; k1 < i_nrow_imputation; k1++)
		{
			wmatLocal[k1] = Rw[k1] * wijk[k1];
			//yicheng_counter = yicheng_counter + 1;
		}

		//if ((l == 1) || (l == 2)) {
		//	cout << "L_Rep_41 running time within one loop L " << l << " is " << MPI_Wtime() - L_Rep_begin_41 << endl;
		//}

		double L_Rep_begin_4 = MPI_Wtime();
		i_loc = 0;
		double d_sum_wij = 0.0;
		Fill_dVector(yi, ncol, 0.0); //initialize vector for column-wise means of all variables  
		for (int i = 0; i < L; i++)
		{
			//-----
			//inner loop within the identical ID 
			//-----
			for (int j = 0; j < L; j++)
			{
				int ID = (int)rbind_ipmat_FHDI(i_loc, 0) - 1; //1st col means ID; "-1" for actual location  

															  //testout
															  //cout<<"k, i, ID: "<<k<<", "<<i<<" , "<<ID<<endl;

				if (ID == i) //as long as the same ID 
				{
					double wi = rbind_ipmat_FHDI(i_loc, 2);
					//double wij= rbind_ipmat(i_loc, 3); //used in mean calculation

					//-------
					//NOTE: use the replicated fractional weight in lieu of ipmat
					//-------
					double wij;
					wij = wmatLocal[i_loc]; //k means current Jackknifed column

											//cout << "Mynode: " << mynode << " L: " << l << " j: " << j << ", i: " << i << ", ID is: " << ID << ", l: " << l << ", i_loc: " << i_loc<< endl;																					//cout << "k value is: " << k << endl;
											//cout << "i_loc is: " << i_loc << endl;
											//----
											//accumulate fractional weight
											//variable-wise weight summation
											//----
					d_sum_wij += wi*wij; //WGT*FWGT
										 //----
										 //do weighted summation with all the imputed cells 
										 //for current row 
										 //----
					for (int i_var = 0; i_var < ncol; i_var++)
					{
						yi[i_var] = yi[i_var] + wi*wij * rbind_ipmat_FHDI(i_loc, 4 + i_var);
						//cout << "i_ioc::" << i_loc << ", L::" << l << ", i_var:: " << i_var << endl;
					}

					//----
					//increment location for next row 
					//----
					//cout <<"K: "<<k<<", I_loc inside i: " << i_loc << ", yicheng_counter is: " << yicheng_counter << endl;
					//yang_counter = yang_counter + 1;
					i_loc++;

				}

				if (ID > i) { break; } //exit inner loop 
			}

		} //end of main loop of i = [0, nrow)
		  //cout << "End inside i loop of i_loc" << endl;
		  //cout << "I_loc inside k: " << i_loc << endl;

		//if ((l == 1) || (l == 2)) {
		//		cout << "L_Rep_4 running time within one loop L " << l << " is " << MPI_Wtime() - L_Rep_begin_4 << endl;
		//}

		if (fabs(d_sum_wij) == 0.0)
		{
			TestOut << "ERROR! zero sum of fractional weight at Jackknifed row :" << l << endl;
			return;
		}
		//for (int i_var = 0; i_var < ncol; i_var++)
		//{
		//	cout << "L: " << l << ", yi[" << i_var << "]: " << yi[i_var] << endl;

		//}
		//cout << "L " << l << " Dundun1" << endl;
		//----
		//store Jackknife mean estimator 
		//----
		double L_Rep_begin_5 = MPI_Wtime();
		for (int i_var = 0; i_var < ncol; i_var++) //size of columns of ipmat matrix of C++
		{
			//cout << "K: " << k << ", i_loc: " << i_loc << ", d_sum_wij: " << d_sum_wij << endl;
			double d_temp = yi[i_var] / d_sum_wij;

			//-----
			//bar{y}_i^(k)
			//-----
			//NOTE: R works column-by-column 
			//hence, below sequence is different from C++ ordering 
			//final_full_data[i_var*nrow + i] = d_temp;  //note: i=current row
			y_bar_i_k[l - (mynode - 1)*numWorkPerProc][i_var] = d_temp;  //note: k= replicate; i_var = variable
																		 //cout <<"L: " <<l<<", l - (mynode - 1)*numWorkPerProc = " << l - (mynode - 1)*numWorkPerProc <<", i_var = " << i_var <<endl;
		}

		//if ((l == 1) || (l == 2)) {
		//		cout << "L_Rep_5 running time within one loop L " << l << " is " << MPI_Wtime() - L_Rep_begin_5 << endl;
		//}


		//testout
		/*
		double* d_temp_wmat1 = new double[i_nrow_imputation];
		for(int j=0; j<i_nrow_imputation; j++) d_temp_wmat1[j] = wmat[j][l];
		RPrint("wmat[,l]:");
		RPrint(d_temp_wmat1, i_nrow_imputation);
		delete[] d_temp_wmat1;
		*/

		//--------------------
		//local deallocation
		//--------------------
		delete[] idd;
		//delete[] d_cellp;

		//if ((l == 1) || (l == 2)) {
		//		cout << "L_Rep iteration running time within one loop L " << l << " is " << MPI_Wtime() - L_Rep_begin << endl;
		//}
	} //end of main loop for L

	delete[] yi;
	//cout << "YYC Running time of L_Replication at node " << mynode << " is " << MPI_Wtime() - L_begin << endl;
	//cout << "i_nrow_imputation after: " << i_nrow_imputation << endl;

	//--------------------------------------------------------------------------------------------------------------------------------------------
	//double i_bar_begin = MPI_Wtime();

	//double** y_bar_i_k_last = New_dMatrix(numWorkLocalLast, ncol);

	double** y_bar_i_k_Recv = New_dMatrix(numWorkPerProc, ncol);
	double** y_bar_i_k_Recv_last = New_dMatrix(numWorkLocalLast, ncol);

	//-------------------------------------------------------------------------------------------------------------------------------------------

	double communication_begin = MPI_Wtime();
	if (mynode != 0) {
		//cout << "send in mynode " << mynode << ", L and ncol are " << L_temp << "," << ncol << endl;
		MPI_Send(y_bar_i_k[0], (L_temp*ncol), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}

	//cout << "I'm there" << endl;
	if (mynode == 0) {
		for (int j = 1; j < totalnodes; j = j + 1) {

			//--------------------------------------------------------


			if (j != (totalnodes - 1)) {
				//cout << "Recv in mynode " << mynode << "from slave " << j << ", L and ncol are " << numWorkPerProc << "," << ncol << endl;
				MPI_Recv(y_bar_i_k_Recv[0], (numWorkPerProc*ncol), MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
				int counter = 0;
				for (int l = (j - 1)*numWorkPerProc; l < j*numWorkPerProc;l = l + 1) {
					for (int k1 = 0; k1 < ncol; k1 = k1 + 1) {
						y_bar_i_k_Summary[l][k1] = y_bar_i_k_Recv[counter][k1];
						//cout << "mynode: " << j << ", l: " << l << ", k1: " << k1 << ", counter: " << counter << endl;
					}
					counter = counter + 1;
				}
				//cout << "I'm at node " << j << endl;
			}

			if (j == (totalnodes - 1)) {
				//cout << "Recv in mynode " << mynode << "from slave " << j << ", L and ncol are " << L_temp << "," << ncol << endl;
				MPI_Recv(y_bar_i_k_Recv_last[0], (numWorkLocalLast*ncol), MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
				int counter = 0;
				for (int l = (j - 1)*numWorkPerProc; l < ((j - 1)*numWorkPerProc + numWorkLocalLast); l = l + 1) {
					for (int k1 = 0; k1 < ncol; k1 = k1 + 1) {
						y_bar_i_k_Summary[l][k1] = y_bar_i_k_Recv_last[counter][k1];
						//cout << "mynode: " << j << ", l: " << l << ", k1: " << k1 << ", counter: " << counter << endl;
					}
					counter = counter + 1;
				}
				//cout << "I'm at last node " << j << endl;
			}
			//------------------------------------------------------- 

		}//end of if
	}//end of mynode=0

	//if (mynode == 0) cout << "I'm out--------" << endl;
	//----------------------------------------------------------------------------------------------------------------------------------------


	MPI_Barrier(MPI_COMM_WORLD);
	//cout << "YYC Running time of variance_communication at node " << mynode << " is " << MPI_Wtime() - communication_begin << endl;


	delete[] wmatLocal;
	//delete[] wmatLocallast;

	//delete[] wmatRecv;//test
	//delete[] wmatRecvlast;//test
	//Del_dMatrix(wmat_buffer, i_nrow_imputation, L);//test

	//if (mynode == 0) cout << " ========= Variance_Est_FHDI.. has successfully finished!" << endl;

	//-------------
	//deallocation
	//-------------
	//delete[] w; 
	//delete[] id;
	delete[] cn;
	delete[] d_id_FHDI;
	Del_dMatrix(d_iy, nrow_d_iy, ncol);
	Del_dMatrix(d_cx, nrow, ncol);
	delete[] i_locg;
	delete[] d_rr0;
	delete[] d_w1;
	//Del_dMatrix(wmat, nrow_dat2_FHDI, L);
	delete[] rw0;
	delete[] Rw;
	delete[] wijk;
	delete[] rbind_ipmat;

	Del_dMatrix(y_bar_i_k, L_temp, ncol);
	Del_dMatrix(y_bar_i_k_Recv, numWorkPerProc, ncol);
	Del_dMatrix(y_bar_i_k_Recv_last, numWorkLocalLast, ncol);
	return;
}

