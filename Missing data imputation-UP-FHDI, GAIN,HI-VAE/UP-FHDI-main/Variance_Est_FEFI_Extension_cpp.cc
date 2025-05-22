#include "matrix_utility_FHDI.cc"
#include <vector>
#include <string>
#include <mpi.h>
//#include "Rep_CellP_MPI.cc"
//#include "Cell_Prob_Extension_cpp_MPI.cc"

//void RepWeight(const int n, double** d_rw)
////Description -------------------------------------
//// Jackknife replicate weight for simpler random sampling
//// 
//// original R code: Dr. Im, J. and Dr. Kim, J. 
//// c++ code: 		Dr. Cho, I. 
//// All rights reserved
//// 
//// updated: Nov 17, 2016
////IN   : int n = matrix dimension
////OUT  : double d_rw[n,n] = replicate weights 
////--------------------------------------------------
//{
//
//
//	const double d_rw0 = (1.0*n) / (n - 1);
//
//	//--------
//	//initialize
//	//--------
//	Fill_dMatrix(d_rw, n, n, d_rw0);
//
//	//--------
//	//put 0 into diagonal terms
//	//--------
//	for (int i = 0; i<n; i++)
//	{
//		d_rw[i][i] = 0.0;
//	}
//
//	return;
//}

//void Rep_CellP(double** d_cx, const int nrow, const int ncol, double** d_rw, int*  id,
//
//	List_FHDI        &List_rst_prob,
//	List_string_FHDI &List_rst_name,
//	std::vector<std::string> &s_ncx, ofstream& TestOut)
//	//Description============================================
//	// compute cell probability using replicate weight rw
//	// 
//	// R code: Dr. Im, J., and Dr. Kim, J. 
//	// C++   : Dr. Cho, I.
//	// All rights reserved
//	// Last update: March 28, 2017
//	//
//
//	//IN   : double d_cx[nrow, ncol] = categoraized matrix
//	//IN   : double d_rw[nrow, nrow] = replicate weights
//	//IN   : int    id[nrow] = index of rows
//	//
//	//below two lists have meaningful values up to i_nc rows  
//	//OUT  : List_FHDI List_rst_prob(nrow->i_nc); //list of joint probabilities for all missing patterns 
//	//OUT  : List_string_FHDI List_rst_name(nrow->i_nc); //names of joint probabilities for all missing patterns 
//	//OUT  : std::vector<std::string> s_ncx; //uniqe cn0
//	//======================================================== 
//{
//	//--------------
//	//make a condensed expression "cn0" of cx, i.e. z
//	//--------------
//	//std::string cn0[nrow];
//	//std::string cn0_backup[nrow];
//	int mynode, totalnodes;
//	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
//	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
//	MPI_Status status;
//
//	double Rep_CellP_begin = MPI_Wtime();
//	std::string *cn0 = new std::string[nrow];
//	std::string *cn0_backup = new std::string[nrow];
//
//	Trans(d_cx, nrow, ncol, cn0);//make a condensed string expression with a given array
//
//	for (int i = 0; i<nrow; i++) cn0_backup[i] = cn0[i];
//
//	//---------------------
//	//SORT & UNIQUE patterns of cn0
//	//---------------------
//	//std::string s_cn0_temp[nrow]; 
//	std::string *s_cn0_temp = new std::string[nrow];
//	for (int i = 0; i<nrow; i++) s_cn0_temp[i] = cn0[i];
//	std::sort(s_cn0_temp, s_cn0_temp + nrow);
//
//	//------------
//	//memorize observed patterns 
//	//------------
//	//std::vector<std::string> s_ncx; //uniqe cn0
//
//	int i_count_cn0 = 0; //total number of unique cn0 
//	std::string s_temp;
//	for (int i = 0; i<nrow; i++)
//	{
//		s_temp = s_cn0_temp[i]; //get a string from the sorted strings 
//		for (int j = 0; j<nrow; j++) //search all rows 
//		{
//			//----
//			//below condition is needed for finding UNIQUE pattern
//			//----
//			if (s_temp.compare(cn0_backup[j]) == 0) //0: equal string
//			{
//				s_ncx.push_back(cn0_backup[j]);  //store the found observed pattern
//
//												 //----------
//												 //remove all identical string after the current string
//												 //----------
//				for (int k = j; k<nrow; k++)
//				{
//					if (s_temp.compare(cn0_backup[k]) == 0) //0: equal string
//					{
//						cn0_backup[k] = ""; //nullify for the next search
//					}
//				}
//
//				i_count_cn0++;
//				break;
//			}
//
//		}
//	}
//	//Now, i_count_cn0 means the total number of unique sorted strings
//	const int i_nc = i_count_cn0;
//	cout << "Yang Running time of Rep_cell1 at node " << mynode << " is " << MPI_Wtime() - Rep_CellP_begin << endl;
//	//testout 
//
//	//testout
//	//RPrint("=====in Rep_CellP ========");
//	//RPrint("s_ncx"); RPrint(s_ncx);
//	//RPrint("i_nc"); RPrint(i_nc);
//	/*cout<<"=====in Rep_CellP ========"<<endl;
//	cout<<"s_ncx"<<endl;
//	for(int i=0; i<(int)s_ncx.size(); i++) cout<<s_ncx[i]<<" ,  ";
//	cout<<endl;
//	cout<<"i_nc: "<<i_nc<<endl;
//	*/
//
//	//-----------------------------
//	//calculate joint probability and names of all missing patterns
//	//using the Jackknife replicate weights
//	//------------------------------
//	//List_FHDI        List_rst_prob(i_nc); //list of joint probabilities for all missing patterns 
//	//List_string_FHDI List_rst_name(i_nc); //names of joint probabilities for all missing patterns 
//
//	std::vector<double> jp_prob_return;
//	std::vector<std::string> jp_name_return;
//	//std::vector<double> w_UserDefined; 
//	double* w_UserDefined = new double[nrow];
//	double Rep_CellP_begin2 = MPI_Wtime();
//	if (mynode == 0) cout << "The value of i_inc is " << i_nc << endl;
//	for (int i = 0; i<i_nc; i++)
//	{
//		//---
//		//search current missing pattern from all strings
//		//---
//		std::string s_temp = s_ncx[i];
//		int i_loc = 0;
//		for (int j = 0; j<nrow; j++)
//		{
//			if (s_temp.compare(cn0[j]) == 0)
//			{
//				i_loc = j;
//				break;
//			}
//		}
//		//testout
//		/*
//		cout<<"loop i (1:i_nc) :"<<i<<"  found i_loc:"<<i_loc<<endl;
//		if(i==7)
//		{
//		for(int j_temp=0; j_temp<nrow; j_temp++) cout<<d_rw[j_temp][i_loc]<<",  ";
//		}
//		cout<<endl;
//		*/
//		
//		//----
//		//joint probability and names
//		//----
//		
//		jp_prob_return.clear();
//		jp_name_return.clear();
//		//w_UserDefined.clear();
//		//for(int j=0; j<nrow; j++) w_UserDefined.push_back(d_rw[j][i_loc]) ; 
//		for (int j = 0; j<nrow; j++) w_UserDefined[j] = d_rw[j][i_loc];
//
//		Cell_Prob_Extension_cpp(d_cx, nrow, ncol,
//			jp_prob_return,
//			jp_name_return,
//			w_UserDefined, id, TestOut);
//
//		//---
//		//prep return
//		//---
//		List_rst_prob.put_block(i, jp_prob_return); //jth row has joint prob
//		List_rst_name.put_block(i, jp_name_return); //jth row has name of the joint prob
//
//	}
//	cout << "Yang Running time of Rep_cell2 at node " << mynode << " is " << MPI_Wtime() - Rep_CellP_begin2 << endl;
//	//testout
//	RPrint("YYC List_rst_name"); List_rst_name.print_List_string_FHDI();
//	RPrint("YYC List_rst_prob"); List_rst_prob.print_List_FHDI();
//
//
//	delete[] cn0;
//	delete[] cn0_backup;
//	delete[] s_cn0_temp;
//
//	delete[] w_UserDefined;
//	
//	return;
//}


void Variance_Est_FEFI_Extension_cpp(double** y, double** z, const int nrow, const int ncol,
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
	List_FHDI   &List_rst_prob,
	List_string_FHDI &List_rst_name,
	std::vector<std::string> &s_ncx,
	std::string s_M,
	double** y_bar_i_k_Summary, ofstream& TestOut, ofstream& TestOut_Slave1, ofstream& TestOut_Slave2)

	//Description----------------------
	//estimate variance for FEFI using Jackknife method 
	//  Algorithm: 
	//
	// original R code: Dr. Im, J. and Dr. Kim, J. 
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: March 28, 2017
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
	//
	//OUT  : double** wmat = New_dMatrix(nrow_dat2_FEFI, L=nrow); //nrow_dat2_FHDI = rows of w1
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
	//RPrint("=========== Begin Variance Estimation of FEFI ================");
	//cout<<"=========== Begin Variance Estimation of FEFI  ================"<<endl;


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
	int nrow_dat2_FEFI = rbind_ipmat_FEFI.size_row();
	int nrow_dat2_FHDI = rbind_ipmat_FHDI.size_row();
	int nrow_mox = rbind_mox.size_row();
	int L = nrow; //size of d_rw 

	//if (mynode == 0) cout << "nrow_dat2_FEFI is " << nrow_dat2_FEFI << endl;


					   //-- Cacluate number of works assigned to each processor
   //int numWorkPerProc = (int)floor((double)L / (double)(totalnodes - 1));
	int numWorkPerProc = (int)floor(1.0*L / (1.0*totalnodes - 1));
	int numWorkLocalLast = L - numWorkPerProc * (totalnodes - 2);

	//if (mynode == 0) {
	//	RPrint("The totalnodes (variance) is---------------------\n");
	//	RPrint(totalnodes);
	//	RPrint("---------------------------\n");

	//	RPrint("The numWorkPerProc (variance) is---------------------\n");
	//	RPrint(numWorkPerProc);
	//	RPrint("---------------------------\n");

	//	RPrint("The numWorkLocalLast (variance) is---------------------\n");
	//	RPrint(numWorkLocalLast);
	//	RPrint("---------------------------\n");

	//	RPrint("The L (variance) is---------------------\n");
	//	RPrint(L);
	//	RPrint("---------------------------\n");
	//}
	//--------------------
	//get ready id table of FEFI
	//--------------------
	double* d_id_FEFI = new double[nrow_dat2_FEFI];
	for (int i = 0; i < nrow_dat2_FEFI; i++) d_id_FEFI[i] = rbind_ipmat_FEFI(i, 0); //id, 1st col 
	std::vector<double> v_table_name_id_FEFI; //same as "nimp" in R version
	std::vector<int>    v_table_count_id_FEFI;//same as "nimp" in R version 
	table_cpp(d_id_FEFI, nrow_dat2_FEFI, v_table_name_id_FEFI, v_table_count_id_FEFI);

	//testout 
	//cout<<"Main Var Est FEFI: after table_"<<endl;

	//---------------------
	//imputed real data matrix
	//---------------------
	int nrow_d_iy = nrow_dat2_FEFI; //default
	if (s_M == "FEFI") { nrow_d_iy = nrow_dat2_FEFI; }
	if (s_M == "FHDI") { nrow_d_iy = nrow_dat2_FHDI; }
	double** d_iy = New_dMatrix(nrow_d_iy, ncol); //imputed matrix of real values 
	for (int i = 0; i < ncol; i++)
	{
		for (int j = 0; j < nrow_d_iy; j++)
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
	//cout << "YYC Running time of startup1 at node " <<mynode<<" is "<< MPI_Wtime() - startup << endl;
	//testout 
	/*
	cout<<"Main Var Est FEFI: after getting ipmat"<<endl;
	cout<<"d_iy matrix"<<endl;
	for(int i=0; i<nrow_d_iy; i++)
	{
	for(int j=0; j<ncol; j++)
	{
	cout<<d_iy[i][j]<<",  ";
	}
	cout<<endl;
	}
	*/

	//----------------------
	//categorized matrix
	//----------------------
	//double startup2 = MPI_Wtime();
	double** d_cx = New_dMatrix(nrow, ncol);
	Copy_dMatrix(z, nrow, ncol, d_cx);

	//testout 
	/*
	cout<<"Main Var Est FEFI: after getting z = d_cx"<<endl;
	cout<<"d_cx matrix"<<endl;
	for(int i=0; i<nrow; i++)
	{
	for(int j=0; j<ncol; j++)
	{
	cout<<d_cx[i][j]<<",  ";
	}
	cout<<endl;
	}
	*/

	//---------------------
	//ocg, observed donors for each missing pattern. 
	//---------------------
	//the same as List_ocsg[nrow_mox]
	//---------------------
	int* i_locg = new int[nrow_mox]; //length of each list of ocg
	for (int i = 0; i < nrow_mox; i++)
	{
		int i_temp = 0;
		List_ocsg.get_a_row_size(i, i_temp);
		i_locg[i] = i_temp;
	}

	//testout 
	//cout<<"Main Var Est FEFI: after i_locg"<<endl;

	//testout
	//RPrint("==========DEBUG: after i_locg ================");

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
	//	List_rst_prob,
	//	List_rst_name,
	//	s_ncx, TestOut);
	//cout << "YYC Running time of startup2 at node " << mynode << " is " << MPI_Wtime() - startup2 << endl;
	//testout 
	//cout<<"Main Var Est FEFI: after Rep_CellP"<<endl;			  
	//const int i_nc = (int)s_ncx.size(); //meaningful rows of List's 
	//testout 
	/*cout<<"Main Var Est FEFI: after Rep_CellP"<<endl;
	cout<<"d_rw matrix 4th row"<<endl;
	for(int i=0; i<nrow; i++)
	{
	cout<<d_rw[4][i]<<",  ";
	}
	cout<<endl;
	*/

	//testout
	//if (mynode == 0) {
	//	cout << "Mynode: " << mynode << ", List_rst_name:" << endl;;
	//	List_rst_name.print_List_string_FHDI(); //njp in R ver
	//	cout << "Mynode: " << mynode << ", List_rst_prob:" << endl;;
	//	List_rst_prob.print_List_FHDI(); 		 //Tcellp in R version
	//}
	//testout
	//RPrint("==========DEBUG: after Rep_CellP ================");	

	//--------------------
	//put 1 into fully observed rows
	//--------------------
	double startup3 = MPI_Wtime();
	double* d_rr0 = new double[nrow];
	for (int i = 0; i < nrow; i++)
	{
		double d_prod = 1.0;
		for (int j = 0; j < ncol; j++)
		{
			d_prod = d_prod*d_cx[i][j];
		}
		//-----
		//0: at least one missing;  1: all observed
		//-----
		d_rr0[i] = d_prod;
		if (fabs(d_prod) > 0.0) d_rr0[i] = 1;
	}

	//--------------
	//calculate w1 = sampling weight
	//--------------
	//std::string cn[nrow]; 
	std::string *cn = new std::string[nrow];
	Trans(d_cx, nrow, ncol, cn);
	double* d_w1 = new double[nrow_dat2_FEFI];
	for (int i = 0; i < nrow_dat2_FEFI; i++)
		d_w1[i] = rbind_ipmat_FEFI(i, 2); //3rd column contains WGT

										  //testout
										  //RPrint("rr0: "); RPrint(d_rr0, nrow);
										  //RPrint("w1: "); RPrint(d_w1, nrow_dat2_FEFI);
										  //testout 
										  //cout<<"Main Var Est FEFI: after d_w1"<<endl;			  

										  //------------------
										  //make Covariance matrix for following replication process
										  //------------------
	//int* i_lloc;
	//std::vector<int> v_mox_0;
	//double** d_dy;
	//double** V_var; //covariance matrix of dy
	//List_FHDI List_V(nrow_mox); //storage of covariance matrix
	//							//Note: store cov mat by "row-first" rule 
	//for (int i = 0; i<nrow_mox; i++)
	//{
	//	int i_size_lloc = i_locg[i]; //length of the ocg associated with current missing row

	//								 //-----------
	//								 //when zero size continue to next iteration
	//								 //-----------
	//	if (i_size_lloc <= 0) { continue; }
	//	i_lloc = new int[i_size_lloc];
	//	for (int j = 0; j< i_size_lloc; j++) i_lloc[j] = (int)List_ocsg(i, j); //ith row, jth entity

	//																		   //-----
	//																		   //find missing column in current row
	//																		   //------
	//	v_mox_0.clear();
	//	for (int j = 0; j<ncol; j++)
	//	{
	//		if (rbind_mox(i, j) == 0.0) v_mox_0.push_back(j + 1); //ACTUAL zero column id 
	//	}
	//	const int i_size_v_mox_0 = (int)v_mox_0.size();

	//	//-------
	//	//extract matrix of missing patterns
	//	//-------
	//	d_dy = New_dMatrix(i_size_lloc, i_size_v_mox_0);
	//	V_var = New_dMatrix(i_size_v_mox_0, i_size_v_mox_0); //column-wise covariance 
	//	for (int j = 0; j< i_size_lloc; j++)
	//	{
	//		for (int k = 0; k<i_size_v_mox_0; k++)
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

	//testout 
	//cout<<"Main Var Est FEFI: after V_var"<<endl;			  

	//testout
	//RPrint("==========DEBUG: after V_var ================");	

	//testout
	//RPrint("FEFI  List_V");
	//List_V.print_List_FHDI(); 

	//----------------
	//wmat: Replication Weights
	//----------------
	//double** wmat_buffer = New_dMatrix(nrow_dat2_FEFI, L); //nrow_dat2_FEFI = rows of w1

	//------------------------------
	//------------------------------
	//MAIn loop for L replications
	//------------------------------
	//------------------------------
	double* rw0 = new double[nrow];

	int i_sum_Rw = 0;
	for (int i = 0; i < nrow; i++) i_sum_Rw += v_table_count_id_FEFI[i];
	double* Rw = new double[i_sum_Rw];

	double* wijk = new double[nrow_dat2_FEFI]; //FWGT from ipmat 

											   //testout
											   //RPrint("==========DEBUG: begin main loop of l=L ==============");	
											   //testout 
											   //cout<<"Main Var Est FEFI: begin main loop of l=L"<<endl;			  
											   //testout
	//RPrint("==========DEBUG: begin main loop of l=L ==============\n");
	//testout 
	//cout << "Main Var Est FEFI: begin main loop of l=L" << endl;


	int i_nrow_imputation = nrow;
	if (s_M == "FEFI") i_nrow_imputation = nrow_dat2_FEFI;
	if (s_M == "FHDI") i_nrow_imputation = nrow_dat2_FHDI;
	//cout<<"(Variance) i_nrow_imputation is "<< i_nrow_imputation<< endl;

	//================================
	// parallelization of main loop
	//================================
			 //wmat = New_dMatrix(i_nrow_imputation, L);
	//double** wmatRecv = New_dMatrix(i_nrow_imputation, numWorkPerProc);	// wmat received from slave processors
	//double** wmatRecvlast = New_dMatrix(i_nrow_imputation, numWorkLocalLast);// receive from last processor
	//double** wmatLocal = New_dMatrix(i_nrow_imputation, numWorkPerProc);
	//double** wmatLocallast = New_dMatrix(i_nrow_imputation, numWorkLocalLast);

	//double* wmatRecv = new double[i_nrow_imputation*numWorkPerProc];
	//double* wmatRecvlast = new double[i_nrow_imputation*numWorkLocalLast];
	double* wmatLocal = new double[i_nrow_imputation];
	//double* wmatLocallast = new double[i_nrow_imputation*numWorkLocalLast];

	int L_temp = 0;
	if (mynode != (totalnodes - 1)) L_temp = numWorkPerProc;
	if (mynode == (totalnodes - 1)) L_temp = numWorkLocalLast;

	//cout<<"Mynode: "<<mynode<<", L_temp: "<< L_temp <<endl;
	double** y_bar_i_k = New_dMatrix(L_temp, ncol);
	Fill_dMatrix(y_bar_i_k, L_temp, ncol, 0.0);

	int i_loc = 0;
	double* yi = new double[ncol];
	// wmat for slave processors

	//===============================
	// specify startpoint and end point for slave processors
	//===============================
	int startpoint = 0;
	int endpoint = 0;

	if (mynode != 0 && mynode != (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = mynode*numWorkPerProc;

		// RPrint("(Variance)Strating point and ending point on node ");RPrint(mynode);
		 //RPrint("are: \n");
		// RPrint(startpoint); RPrint(endpoint);
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
	   // RPrint(startpoint); RPrint(endpoint);
		if (endpoint - startpoint != L_temp) {
			TestOut << "y_bar_i_k boundary ERROR!!!" << endl;
			return;
		}
	}
	//cout << "Mynode " << mynode << ", " << startpoint << "<= x < " << endpoint << endl;
	//cout << "YYC Running time of startup3 at node " << mynode << " is " << MPI_Wtime() - startup3 << endl;
	//cout << "YYC Running time of startup at node " << mynode << " is " << MPI_Wtime() - startup << endl;

	//if (mynode == 0) cout << "FEFI Get Started 2" << endl;

	double L_begin = MPI_Wtime();
	//int yicheng_counter = 0;
	for (int l = startpoint; l < endpoint; l = l + 1)
		//---------------------------------------------------------------------------------------Same for all processors except for MPI_Send(Yicheng)
	{
		Fill_dVector(wmatLocal, i_nrow_imputation, 0.0);
		//RPrint("Yicheng, mynode is ");RPrint(mynode);RPrint(", L is ");RPrint(l);RPrint("\n");
		//RPrint("Debug Yicheng_1, mynode is ");RPrint(mynode);RPrint("\n"); 
		//wmatLocal = new double[i_nrow_imputation];
		//testout

		//RPrint("==========DEBUG: Inside main loop of l and on Proc: \n"); 
		//RPrint(mynode);RPrint("\n");
		//RPrint(l); 		
		//cout<<"main loop of l:  "<<l<<endl;

		//-------
		//replicate weight from lth column
		//-------
		for (int i = 0; i < nrow; i++) rw0[i] = d_rw(i, l); //l_th column 
		int i_sum = 0;
		for (int i = 0; i < nrow; i++)
		{
			for (int j = 0; j < v_table_count_id_FEFI[i]; j++) Rw[i_sum++] = rw0[i];
		}

		//---------
		//FWGT of ipmat
		//----------
		for (int i = 0; i < nrow_dat2_FEFI; i++)
			wijk[i] = rbind_ipmat_FEFI(i, 3); //4th column is FWGT 

											  //----------
											  //joint probability associated with current string
											  //-----------
		std::string cn_current = cn[l]; //lth string 
		std::vector<int> v_ncx_cn;
		which(s_ncx, cn_current, v_ncx_cn); //actual location 
											//const int i_size_v_ncx_cn = (int)v_ncx_cn.size(); //MUST BE "1"

		int i_size_cellp = 0;
		List_rst_prob.get_a_row_size(v_ncx_cn[0] - 1, i_size_cellp); //get a size of the row in the list 

																	 //-----------
																	 //when zero size continue to next iteration
																	 //-----------
		if (i_size_cellp <= 0) { continue; }
		double* d_cellp = new double[i_size_cellp];
		List_rst_prob.get_block(v_ncx_cn[0] - 1, d_cellp); //-1 for actual row location 

														   //testout
														   //RPrint(" ======== in Main Loop l+1: "); RPrint(l+1); 
														   //RPrint("Rw"); RPrint(Rw, i_sum_Rw); 
														   //RPrint("wijk"); RPrint(wijk, nrow_dat2_FEFI); 
														   //RPrint("d_cellp"); RPrint(d_cellp, i_size_cellp); 

														   //----------------------------------------
														   //1. if the deleted is missing unit, no further action is taken
														   //2. if the deleted is observed unit, then the fractional weights are re-computed 
														   //----------------------------------------
		int* idd = new int[nrow_mox]; //location of the deleted donor in ocg 
		Fill_iVector(idd, nrow_mox, 0);
		//RPrint("Debug Yicheng_2, mynode is ");RPrint(mynode);RPrint(" L is ");RPrint(l);RPrint("\n");
		if (fabs(d_rr0[l]) > 0)
		{
			//---------------------
			//locations of the deleted unit in observed list
			//---------------------
			std::vector<int> v_lg; //Actual locations 
			v_lg.clear();
			for (int j = 0; j < nrow_mox; j++) //list length 
			{
				for (int k = 0; k < i_locg[j]; k++) //a row in the List
				{
					int i_temp_lg = (int)List_ocsg(j, k);
					if (i_temp_lg == (l + 1)) //+1 for actual location  
					{
						v_lg.push_back(j + 1); //actual row location 
						idd[j] = k + 1; //actual location 
						break;
					}
				}
			}
			const int nlg = (int)v_lg.size();

			//testout
		//	if (mynode == 1)
		//	{
		//		RPrint(" in condition rr0[l]!=0 at l+1 ="); RPrint(l + 1);
		//		RPrint("idd:"); RPrint(idd, nrow_mox);
		//		RPrint("lg :"); RPrint(v_lg);
		//		RPrint("nlg:"); RPrint(nlg);
		//	}
			//--------------------------
			//Adjust fractional weights for all units in lg
			//--------------------------
			if (nlg > 0)
			{
				for (int j = 0; j < nlg; j++)
				{
					int i_row_lg = v_lg[j] - 1; // row number [0,...) 
					double* d_1_mox = new double[ncol];
					for (int k = 0; k < ncol; k++) d_1_mox[k] = rbind_mox(i_row_lg, k);

					//---
					//actual col number of missing cell in current missing row
					//---
					std::vector<int> v_rloc; v_rloc.clear();
					for (int k = 0; k < ncol; k++)
					{
						if (d_1_mox[k] == 0.0) { v_rloc.push_back(k + 1); } //actual col
					}
					//const int nrloc = (int)v_rloc.size();
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
					//FEFI
					//------------------------
					//------------------------
					std::vector<int> v_elog;
					if (s_M.compare("FEFI") == 0) //0=equal 
					{
						//-----
						//find locations of mlog in dat2$ID
						//v_mlog contains the row numbers that have the same string as
						//current missing row 
						//nmlog = n(v_mlog)
						//-----
						v_elog.clear();
						for (int k1 = 0; k1 < nmlog; k1++) //loop for mlog
						{
							int i_temp1 = id[v_mlog[k1] - 1]; //dat1$ID in R version  
							for (int k2 = 0; k2 < nrow_dat2_FEFI; k2++)
							{
								int i_temp2 = rbind_ipmat_FEFI(k2, 0); //1st col is dat2$ID
								if (i_temp1 == i_temp2)
								{
									v_elog.push_back(k2 + 1); //actual location 
								}
							}
						}
						const int i_size_v_elog = (int)v_elog.size();
						//testout
						//RPrint(" elog: "); RPrint(v_elog); 

						//----
						//donor id
						//-----
						std::vector<int> v_did; v_did.clear();

						int i_temp_lg_j = v_lg[j] - 1; //-1 for actual location  
						for (int k1 = 0; k1 < i_locg[i_temp_lg_j]; k1++) //a row in the List 
						{
							v_did.push_back((int)List_ocsg(i_temp_lg_j, k1));
						}
						const int i_size_nic = (int)v_did.size(); //length of did 

																  //----
																  //donor string patterns
																  //----
																  //std::string s_icn[i_size_nic];
						std::string *s_icn = new std::string[i_size_nic];
						for (int k1 = 0; k1 < i_size_nic; k1++)
							s_icn[k1] = cn[v_did[k1] - 1]; //-1 for actual location  		

														   //testout 
														   //RPrint("did :"); RPrint(v_did);
														   //RPrint("icn :"); RPrint(s_icn, i_size_nic);

														   //------
														   //unique icn
														   //-------
						std::vector<std::string> v_unique_icn; v_unique_icn.clear();
						//std::string s_icn_backup[i_size_nic];
						std::string *s_icn_backup = new std::string[i_size_nic];
						for (int k1 = 0; k1 < i_size_nic; k1++) s_icn_backup[k1] = s_icn[k1];

						for (int k1 = 0; k1 < i_size_nic; k1++)
						{
							std::string s_uicn_temp = s_icn[k1];
							for (int k2 = 0; k2 < i_size_nic; k2++)
							{
								if (s_uicn_temp.compare(s_icn_backup[k2]) == 0)
								{
									//store the found unique string pattern 
									v_unique_icn.push_back(s_uicn_temp);

									//nullify all the remaining unit that has the same string
									for (int k3 = k2; k3 < i_size_nic; k3++)
									{
										if (s_uicn_temp.compare(s_icn_backup[k3]) == 0)
											s_icn_backup[k3] = ""; //nullify
									}

									break;
								}
							}
						}
						const int i_unic = (int)v_unique_icn.size();
						//RPrint("Debug Yicheng_3, mynode is ");RPrint(mynode);RPrint(" L is ");RPrint(l);RPrint("\n"); 
							//testout 
						//	if (mynode = 1)
						//	{
						//		RPrint("v_unique_icn :"); RPrint(v_unique_icn);
						//		RPrint("i_unic :"); 	  RPrint(i_unic);
						//	}

							//----------
							//sort the unique donors
							//----------
							//std::string s_unique_icn_sorted[i_unic]; //sorted donors
						std::string *s_unique_icn_sorted = new std::string[i_unic]; //sorted donors
						for (int k1 = 0; k1 < i_unic; k1++) s_unique_icn_sorted[k1] = v_unique_icn[k1];
						std::sort(s_unique_icn_sorted, s_unique_icn_sorted + i_unic);

						//---------
						//match the sorted unique donors to njp
						//---------
						//make "njp" 
						int i_temp2 = 0; List_rst_name.get_a_row_size(0, i_temp2);
						//std::string s_njp[i_temp2]; 
						std::string *s_njp = new std::string[i_temp2];
						for (int k1 = 0; k1 < i_temp2; k1++) s_njp[k1] = List_rst_name(0, k1);

						std::vector<int> v_icn_njp; v_icn_njp.clear();
						for (int k1 = 0; k1 < i_unic; k1++)
						{
							std::string s_temp_cp = s_unique_icn_sorted[k1];

							for (int k2 = 0; k2 < i_temp2; k2++)
							{
								if (s_temp_cp.compare(s_njp[k2]) == 0)
								{
									v_icn_njp.push_back(k2 + 1); break;
								} //+1 for actual location 
							}
						}
						const int i_size_v_icn_njp = (int)v_icn_njp.size();
						//RPrint("Debug Yicheng_4, mynode is ");RPrint(mynode);RPrint(" L is ");RPrint(l);RPrint("\n");

						//-------------------------------------
						//get joint probability at the matched location only
						//-------------------------------------
						//-----------
						//when zero size continue to next iteration
						//-----------
						if (i_size_v_icn_njp <= 0) { continue; }
						double* d_cp = new double[i_size_v_icn_njp]; //joint prob
																	 //std::string s_ncp[i_size_v_icn_njp];          //names of the matches
						std::string *s_ncp = new std::string[i_size_v_icn_njp];          //names of the matches
						double d_sum_cp = 0.0;
						for (int k1 = 0; k1 < i_size_v_icn_njp; k1++)
						{
							d_cp[k1] = d_cellp[v_icn_njp[k1] - 1]; //-1 for actual location 
							s_ncp[k1] = s_njp[v_icn_njp[k1] - 1];

							d_sum_cp += d_cp[k1]; //summation of d_cp[]
						}
						for (int k1 = 0; k1 < i_size_v_icn_njp; k1++)
						{
							d_cp[k1] = d_cp[k1] / d_sum_cp;
						}

						//testout
						//RPrint("njp"); RPrint(s_njp, i_temp2); 
						//RPrint("match unique icn and njp :"); RPrint(v_icn_njp); 
						//RPrint("cp"); RPrint(d_cp, i_size_v_icn_njp);
						//RPrint("ncp"); RPrint(s_ncp, i_size_v_icn_njp);						

						//-----------
						//match donors 
						//-----------
						std::vector<int> v_obsp; v_obsp.clear();
						std::vector<double> v_obsp_cp; v_obsp_cp.clear();
						match_FHDI(s_icn, i_size_nic,
							s_ncp, i_size_v_icn_njp,
							v_obsp);
						const int i_size_v_obsp = (int)v_obsp.size();
						for (int k1 = 0; k1 < i_size_v_obsp; k1++)
							v_obsp_cp.push_back(d_cp[v_obsp[k1] - 1]);

						//---------------------------------------
						//Replicated sampling weights for donors
						//---------------------------------------
						//-----------
						//when zero size continue to next iteration
						//-----------
						if (i_size_nic <= 0) { continue; }
						double* d_drw0 = new double[i_size_nic]; //length of did
						double d_sum_drw0 = 0.0;
						for (int k1 = 0; k1 < i_size_nic; k1++)
						{
							d_drw0[k1] = rw0[v_did[k1] - 1];
							d_sum_drw0 += d_drw0[k1];
						}

						//----------
						//update weighted probability
						//----------
						std::vector<std::string> jp_name_icn; jp_name_icn.clear();
						std::vector<double> 	 jp_prob_icn; jp_prob_icn.clear();

						wpct_FHDI(s_icn, i_size_nic, d_drw0,
							jp_name_icn, jp_prob_icn);
						const int i_size_jp_prob_icn = (int)jp_prob_icn.size(); //=i_unic
						for (int k1 = 0; k1 < i_size_jp_prob_icn; k1++)
							jp_prob_icn[k1] = jp_prob_icn[k1] * d_sum_drw0;  //ws.icn

																			 //testout
																			 //RPrint("obsp: "); RPrint(v_obsp);
																			 //RPrint("jp_name_icn :"); RPrint(jp_name_icn);
																			 //RPrint("jp_prob_icn = ws.icn :"); RPrint(jp_prob_icn);



						const int i_current_locg = i_locg[v_lg[j] - 1]; //-1 for actual loc
																		//-----------
																		//when zero size continue to next iteration
																		//-----------
						if (i_size_v_obsp <= 0) { continue; }
						double* d_Fefiw = new double[i_size_v_obsp];
						//---------------
						//weight update I: at deleted donor locations
						//---------------
						if (i_unic == i_current_locg)
						{
							for (int k1 = 0; k1 < i_size_v_obsp; k1++)
								d_Fefiw[k1] = v_obsp_cp[k1];
						}
						//-----------------------------------------------
						//weight update II: at deleted donor locations
						//-----------------------------------------------
						//-----------
						//when zero size continue to next iteration
						//-----------
						if (i_size_nic <= 0) { continue; }
						double* d_drw1 = new double[i_size_nic];
						Fill_dVector(d_drw1, i_size_nic, 0.0); //initialized with 0
						if (i_unic < i_current_locg)
						{
							//-----------
							//match donors 
							//-----------
							std::vector<int> v_temp; v_temp.clear();
							match_FHDI(s_icn, i_size_nic,
								s_unique_icn_sorted, i_unic,
								v_temp);

							for (int k1 = 0; k1 < i_size_nic; k1++)
							{
								if (jp_prob_icn[v_temp[k1] - 1] != 0.0)
									d_drw1[k1] = d_drw0[k1] / jp_prob_icn[v_temp[k1] - 1];
							}

							//------------
							//final updated weights
							//------------
							for (int k1 = 0; k1 < i_size_nic; k1++)
								d_Fefiw[k1] = v_obsp_cp[k1] * d_drw1[k1];

							//testout
						//	if (mynode == 1) {
						//		RPrint("d_drw0 :"); RPrint(d_drw0, i_size_nic);
						//		RPrint("jp_prob_icn[v_temp[k1]-1] :");
						//		for (int k1 = 0; k1 < i_size_nic; k1++)
						//			RPrint(jp_prob_icn[v_temp[k1] - 1]);
						//		RPrint("jp_name_icn[v_temp[k1]-1] :");
						//		std::vector<std::string> v_s_temp; v_s_temp.clear();
						//		for (int k1 = 0; k1 < i_size_nic; k1++)
						//		{
						//			v_s_temp.push_back(jp_name_icn[v_temp[k1] - 1]);
						//		}
						//		RPrint(v_s_temp);
						//		RPrint("v_obsp_cp :"); RPrint(v_obsp_cp);
						//		RPrint("d_drw1 :"); RPrint(d_drw1, i_size_nic);
						//		RPrint("d_Fefiw :"); RPrint(d_Fefiw, i_size_nic);
						//	}

						}

						//------------
						//update wijk: repeated copy, if needed
						//------------
						int i_cycle_wijk = (int)floor(i_size_v_elog*1.0 / i_size_nic*1.0); //expected evenly divisible
						int i_loc_elog = 0; //sequential id

											//testout
											//RPrint("v_elog :"); RPrint(v_elog);
											//RPrint("wijk (updated one only):"); 

						for (int k1 = 0; k1 < i_cycle_wijk;k1++)
						{
							for (int k2 = 0; k2 < i_size_nic; k2++)
							{
								wijk[v_elog[i_loc_elog++] - 1] = d_Fefiw[k2];
								//testout
								//RPrint(d_Fefiw[k2]);
							}
						}

						//RPrint("Yicheng Debug------------------wijk:");RPrint(wijk, i_size_nic);
						//---
						//local deallocation
						//---
						delete[] s_icn;
						delete[] s_icn_backup;
						delete[] d_cp;
						delete[] d_drw0;
						delete[] d_Fefiw;
						delete[] d_drw1;
						delete[] s_unique_icn_sorted;
						delete[] s_njp;
						delete[] s_ncp;
					}

					//local deallocation 
					delete[] d_1_mox;
				}//end of if j<nlg
			} //end of nlg>0

		}// end of if (fabs(d_rr0[l]) > 0) 

		//RPrint("Debug Yicheng_5, mynode is ");RPrint(mynode);RPrint(" L is ");RPrint(l);RPrint("\n");
		//-------------------------------
		//store the updated weights
		//-------------------------------
		int i_nrow_imputation = nrow;
		if (s_M == "FEFI") i_nrow_imputation = nrow_dat2_FEFI;
		if (s_M == "FHDI") i_nrow_imputation = nrow_dat2_FHDI;

		//--------------------------------------------
		// local matrix collection of each processor
		// size of wmatLocal= [i_nrow_imputation,numWorkPerProc];
		// size of wamatLocallast = [i_nrow_imputation,numWorkLocalLast];
		//--------------------------------------------
		//int l_temp = l - (mynode - 1)*numWorkPerProc;


		for (int k1 = 0; k1 < i_nrow_imputation; k1++)
		{
			wmatLocal[k1] = Rw[k1] * wijk[k1];
			//yicheng_counter = yicheng_counter + 1;
		}
		//cout << "wmatlocal finished at L = " << l << endl;
		//--------------------------------------------------------------------------------------------------------------------------
		//---
		//re-initialization! for new jackknifed data
		//---
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
				int ID = (int)rbind_ipmat_FEFI(i_loc, 0) - 1; //1st col means ID; "-1" for actual location  

															  //testout
															  //cout<<"k, i, ID: "<<k<<", "<<i<<" , "<<ID<<endl;

				if (ID == i) //as long as the same ID 
				{
					double wi = rbind_ipmat_FEFI(i_loc, 2);
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
						yi[i_var] = yi[i_var] + wi*wij * rbind_ipmat_FEFI(i_loc, 4 + i_var);
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
		//cout << "L " << l << " Dundun2" << endl;
		//end of Jackknifed data k = [0, nrow)
		 //cout << "y_bar_i_k on node " << mynode << "is--------------------" << endl;
		 //RPrint(y_bar_i_k,L_temp,ncol);
		 //for (int i = 0; i < L_temp; i++)
		 // {

		 // 	for (int j = 0; j < ncol; j++)
		 // 	{
		 // 		TestOut <<i<<" : "<< y_bar_i_k[i][j] << " , ";
		 // 	}
		 // 	TestOut << endl;
		 //
		 //}
		//delete[] yi;
		//--------------------------------------------------------------------------------------------------------------------------
		//cout << "Where is L: " << l << endl;
		//--------------------
		//local deallocation
		//--------------------
		delete[] idd;
		delete[] d_cellp;
		//---------------------------------------------------------------------------------------End of Same for all processors (Yicheng)
	} // End of L replication
	//cout << "YYC Running time of L_Replication at node " << mynode << " is " << MPI_Wtime() - L_begin << endl;
	//cout << "Size of wmatLocal and wmatlocalLast are " << yicheng_counter << " at node " << mynode << " (rbind size is "<< nrow_dat2_FEFI <<" )"<<endl;
	delete[] yi;

	//double i_bar_begin = MPI_Wtime();
	//int L_temp = 0;
	//if (mynode != (totalnodes - 1)) L_temp = numWorkPerProc;
	//if (mynode == (totalnodes - 1)) L_temp = numWorkLocalLast;


	//double** y_bar_i_k = New_dMatrix(L_temp, ncol);
	//Fill_dMatrix(y_bar_i_k, L_temp, ncol, 0.0);

	//double** y_bar_i_k_last = New_dMatrix(numWorkLocalLast, ncol);


	double** y_bar_i_k_Recv = New_dMatrix(numWorkPerProc, ncol);
	double** y_bar_i_k_Recv_last = New_dMatrix(numWorkLocalLast, ncol);



	//if (mynode != 0) {
	//	 // for nrow replicates for ncol variables replicates
	//	int yang_counter = 0;
	//	int i_loc = 0;
	//	double* yi = new double[ncol];
	//	cout << "mynode_lala: " << mynode << ", L_temp is: " << L_temp << endl;
	//	for (int k = 0; k < L_temp; k++) //Jackknife replicates 
	//	{
	//		//---
	//		//re-initialization! for new jackknifed data
	//		//---
	//		i_loc = 0;
	//		double d_sum_wij = 0.0;
	//		Fill_dVector(yi, ncol, 0.0); //initialize vector for column-wise means of all variables  
	//		for (int i = 0; i < L; i++)
	//		{
	//			//-----
	//			//inner loop within the identical ID 
	//			//-----
	//			for (int j = 0; j < L; j++)
	//			{
	//				int ID = (int)rbind_ipmat_FEFI(i_loc, 0) - 1; //1st col means ID; "-1" for actual location  

	//														 //testout
	//														 //cout<<"k, i, ID: "<<k<<", "<<i<<" , "<<ID<<endl;

	//				if (ID == i) //as long as the same ID 
	//				{
	//					double wi = rbind_ipmat_FEFI(i_loc, 2);
	//					//double wij= rbind_ipmat(i_loc, 3); //used in mean calculation

	//					//-------
	//					//NOTE: use the replicated fractional weight in lieu of ipmat
	//					//-------
	//					double wij;
	//					if (mynode != 0 && mynode != (totalnodes - 1)) { wij = wmatLocal[yang_counter]; } //k means current Jackknifed column
	//					if (mynode == (totalnodes - 1)) { wij = wmatLocallast[yang_counter]; }
	//					//cout << "Mynode: " << mynode << " j: " << j << ", i: " << i << ", ID is: " << ID << ", K: " << k << ", i_loc: " << i_loc << ", yang_counter: " << yang_counter << endl;																					//cout << "k value is: " << k << endl;
	//					//cout << "i_loc is: " << i_loc << endl;
	//					//----
	//					//accumulate fractional weight
	//					//variable-wise weight summation
	//					//----
	//					d_sum_wij += wi*wij; //WGT*FWGT
	//										 //----
	//										 //do weighted summation with all the imputed cells 
	//										 //for current row 
	//										 //----
	//					for (int i_var = 0; i_var < ncol; i_var++)
	//					{
	//						yi[i_var] = yi[i_var] + wi*wij * rbind_ipmat_FEFI(i_loc, 4 + i_var);
	//					}

	//					//----
	//					//increment location for next row 
	//					//----
	//					//cout <<"K: "<<k<<", I_loc inside i: " << i_loc << ", yicheng_counter is: " << yicheng_counter << endl;
	//					yang_counter = yang_counter + 1;
	//					i_loc++;
	//					
	//				}

	//				if (ID > i) { break; } //exit inner loop 
	//			}

	//		} //end of main loop of i = [0, nrow)
	//		//cout << "End inside i loop of i_loc" << endl;
	//		//cout << "I_loc inside k: " << i_loc << endl;
	//		if (fabs(d_sum_wij) == 0.0)
	//		{
	//			cout << "ERROR! zero sum of fractional weight at Jackknifed row :" << k << endl;
	//			return;
	//		}


	//		//----
	//		//store Jackknife mean estimator 
	//		//----
	//		for (int i_var = 0; i_var < ncol; i_var++) //size of columns of ipmat matrix of C++
	//		{
	//			//cout << "K: " << k << ", i_loc: " << i_loc << ", d_sum_wij: " << d_sum_wij << endl;
	//			double d_temp = yi[i_var] / d_sum_wij;

	//			//-----
	//			//bar{y}_i^(k)
	//			//-----
	//			//NOTE: R works column-by-column 
	//			//hence, below sequence is different from C++ ordering 
	//			//final_full_data[i_var*nrow + i] = d_temp;  //note: i=current row
	//			y_bar_i_k[k][i_var] = d_temp;  //note: k= replicate; i_var = variable
	//		}

	//	}//end of Jackknifed data k = [0, nrow)
	//	//cout << "y_bar_i_k on node " << mynode << "is--------------------" << endl;
	//	//RPrint(y_bar_i_k,L_temp,ncol);
	//	//for (int i = 0; i < L_temp; i++)
	//	// {

	//	// 	for (int j = 0; j < ncol; j++)
	//	// 	{
	//	// 		TestOut <<i<<" : "<< y_bar_i_k[i][j] << " , ";
	//	// 	}
	//	// 	TestOut << endl;
	//	//
	//	//}
	//	delete[] yi;
	//}//end of mynode!= 0
	//cout << "YYC Running time of y_bar_i_k at node " << mynode << " is " << MPI_Wtime() - i_bar_begin << endl;
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

	//if (mynode==0) cout << "I'm out--------" << endl;
//----------------------------------------------------------------------------------------------------------------------------------------


	MPI_Barrier(MPI_COMM_WORLD);
	//cout << "YYC Running time of variance_communication at node " << mynode << " is " << MPI_Wtime() - communication_begin << endl;

	//MPI_Bcast(wmat_FEFI[0], i_nrow_imputation*L, MPI_DOUBLE,0, MPI_COMM_WORLD);



	//if (mynode == 0) {
	//	//RPrint("yicheng debug------------------congratulation: final wmat in master processor:\n", TestOut);
	//	//RPrint(wmat, i_nrow_imputation, L, TestOut);
	//	TestOut << "yicheng debug------------------congratulation: final wmat in master processor: " << mynode << endl;
	//	TestOut << " i_nrow_imputation: " << i_nrow_imputation << endl;
	//	TestOut << " L: " << L << endl;
	//	for (int i = 0; i < i_nrow_imputation; i++)
	//	{

	//		for (int j = 0; j < L; j++)
	//		{
	//			TestOut <<i<<" : "<< wmat[i][j] << " , ";
	//		}
	//		TestOut << endl;
	//	}
	//	
	//}
	//---------------------------------------------------
	//TestOut << "wmat" << endl;
	//for (int i = 0; i<i_nrow_imputation; i++)
	//{
	//	for (int j = 0; j<L; j++) { TestOut << setw(20) << wmat[i][j]; }
	//	TestOut << endl;
	//}
	//if (mynode == 0) {
	//	RPrint("Yicheng----------------------------------I'm here");
	//}
	//-------------------------------------------------------
	// memory deallocation
	//Del_dMatrix(wmat_buffer, i_nrow_imputation, L);
	//Del_dMatrix(wmatRecv, i_nrow_imputation, numWorkPerProc);
	//Del_dMatrix(wmatRecvlast, i_nrow_imputation, numWorkLocalLast);
	//Del_dMatrix(wmatLocal, i_nrow_imputation, numWorkPerProc);
	//Del_dMatrix(wmatLocallast, i_nrow_imputation, numWorkLocalLast);
	//delete[] wmatRecv;
	//delete[] wmatRecvlast;
	delete[] wmatLocal;
	//delete[] wmatLocallast;

	//delete[] wmatRecv;//test
	//delete[] wmatRecvlast;//test
	//Del_dMatrix(wmat_buffer, i_nrow_imputation, L);//test

	//testout
	//RPrint(" ========= Variance_Est_FEFI.. has successfully finished!");
	//if (mynode == 0) cout << " ========= Variance_Est_FEFI.. has successfully finished!" << endl;
	//testout 
	//cout<<"Main Var Est FEFI: Variance_Est_FEFI.. has successfully finished!"<<endl;			  

	//-------------
	//deallocation
	//-------------
	//delete[] w; 
	//delete[] id; 
	delete[] cn;
	delete[] d_id_FEFI;
	Del_dMatrix(d_iy, nrow_d_iy, ncol);
	Del_dMatrix(d_cx, nrow, ncol);
	delete[] i_locg;
	delete[] d_rr0;
	delete[] d_w1;

	delete[] rw0;
	delete[] Rw;
	delete[] wijk;


	//cout<<"I'm out 1"<<endl;
	//cout << "mynode lalalal" << mynode << ", L_temp is " << L_temp << endl;
	Del_dMatrix(y_bar_i_k, L_temp, ncol);


	Del_dMatrix(y_bar_i_k_Recv, numWorkPerProc, ncol);
	Del_dMatrix(y_bar_i_k_Recv_last, numWorkLocalLast, ncol);

	//cout << "I'm out 3" << endl;
	return;
}

