//#include "ran_FHDI.h" //for uniform distribution //Not used for CRAN Compatibility
#include <cstdlib> //for rand() and srand()
#include <vector>

//#include "Fully_Efficient_Fractional_Imputation_MPI.cc"
//#include "Fractional_Hot_Deck_Imputation_MPI.cc"
//#include "Results_Fully_Efficient_Fractional_Imputation_MPI.cc"
//#include "Results_Fractional_Hot_Deck_Imputation_MPI.cc"

//void yorder(double** y, const int nrow, const int ncol,
//	double* mox_1,
//	std::vector<int> v_loc, int* i_ym_return);

void FHDI_Extension_Bigp_cpp(double** y, double** z, int** r,
	const int nrow, const int ncol, const int i_collapsing,
	std::vector<std::string> jp_name,
	std::vector<double> 	 jp_prob,
	std::string s_M, const int i_M, int** codes,

	rbind_FHDI &rbind_ipmat_FEFI,
	rbind_FHDI &rbind_Resp_FEFI,
	rbind_FHDI &rbind_irmat_FEFI,

	rbind_FHDI &rbind_ipmat_FHDI,
	rbind_FHDI &rbind_Resp_FHDI,
	rbind_FHDI &rbind_irmat_FHDI,

	rbind_FHDI &rbind_uox,
	rbind_FHDI &rbind_mox,
	List_FHDI  &List_ord,
	List_FHDI  &List_ocsg,
	ofstream& TestOut_Slave1)
	//Description=========================================
	// perform
	// Fully Efficient Fractional Imputation OR
	// Fractional Hot Deck Imputation
	// 
	// Algorithm: FEFI of Dr Jae Kwang. Kim and FHDI of Dr Jong Ho. Im
	//
	// original R code: Dr. Im, J. and Dr. Kim, J. 
	// c++ code: 		Dr. Cho, In-Ho 
	// All rights reserved
	// 
	// updated: March 28, 2017
	//----------------------------------------------------
	//IN    : double y(nrow, ncol)= original data matrix with missing cells 
	//IN    : double z(nrow, ncol)= categorized matrix of y. 0 for missing cell
	//IN    : int    r(nrow, ncol) = index matrix of missing unit (0)/observed unit (1)  
	//IN	: vector<string> jp_name  = name of table of joint probability
	//IN	: vector<double> jp_prob  = joint probability 
	//IN  	: string s_M = "FEFI" fully efficient fractional imputation 
	// 					   "FHDI" Fractional Hot Deck Imputation  
	//IN    : int i_M = number of donors used for FHDI
	//IN    : int i_option_collapsing = choice of big-p algorithm 
	//                               0= no big-p algorithms
	//                              !0= perform big-p algorithms
	//IN   : int codes(nrow, i_option_collapsing); // storage to record most correlated variables of mox
	//OUT   : rbind_FHDI  rbind_ipmat_FEFI(4+ncol); //column size is 4+ncol (i.e., for R: ID, FID, WGT, FWGT, Variables)
	//OUT   : rbind_FHDI  rbind_Resp_FEFI(ncol+1);  //separate response matrix  (i.e. for R: unit responses and Resp0)
	//OUT   : rbind_FHDI  rbind_irmat_FEFI(5+ncol); //column size is 5+ncol (i.e. for R:ID, FID, OID, ORDER, FEFIW, CELL )
	//OUT   : rbind_FHDI  rbind_ipmat_FHDI(4+ncol); //column size is 4+ncol
	//OUT   : rbind_FHDI  rbind_Resp_FHDI(ncol+1);  //separate response matrix  
	//OUT   : rbind_FHDI  rbind_irmat_FHDI(5+ncol); //column size is 5+ncol
	//OUT   : rbind_FHDI  rbind_uox (ncol) //observed unique categorized matrix
	//OUT   : rbind_FHDI  rbind_mox (ncol) //missing  unique categorized matrix
	//OUT   : List_FHDI   List_ord(nrow) //but meaningful up to i_count_mox rows
	//OUT   : List_FHDI   List_ocsg(nrow)//but meaningful up to i_count_mox rows
	//====================================================
{
	//-----------------
	//random location using uniform distribution   
	//using Numerical Recipes of Press et al 2007. 
	//-----------------
	//Ran_FHDI myran(1); 	//not used for CRAN Compatibility
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	double FHDI_begin = MPI_Wtime();

	std::srand(123);
	double d_myran = 0.0;


	//-------------
	//column-wise sum of r matrix
	//-------------
	int* i_rn = new int[ncol];
	int i_temp = 0;
	for (int i = 0; i<ncol; i++)
	{
		i_temp = 0;
		for (int j = 0; j<nrow; j++) i_temp += r[j][i];
		i_rn[i] = i_temp;
	}

	//-------------
	//sample weight (default is 1)
	//id array (default is row number)
	//-------------
	double* w = new double[nrow];
	int* id = new int[nrow];
	for (int i = 0; i<nrow; i++)
	{
		w[i] = 1.0;
		id[i] = i + 1; //ACTUAL id
	}

	//--------------
	//locations of missing cells (ml) and observed cells (ol)
	//Note: unlike in Cell_Make..(), std vector is used here
	//--------------
	std::vector<int> ol; //Actual row number 
	std::vector<int> ml; //Actual row number 

	double d_temp = 0.0;
	for (int i_row = 0; i_row<nrow; i_row++)
	{
		d_temp = 1.0;
		for (int i_col = 0; i_col<ncol; i_col++)
		{
			if (z[i_row][i_col] == 0) { d_temp = 0.0; break; } //found zero, i.e. missing cell
		}

		if (fabs(d_temp) > 1e-15) //this row has no missing cells
		{
			ol.push_back(i_row + 1);
		} //actual number of the row having no missing cells

		if (fabs(d_temp) < 1e-15) //this row has AT LEAST one missing cells
		{
			ml.push_back(i_row + 1);
		}  //actual number of the row having missing cells
	}
	const int i_size_ol = (int)ol.size();
	const int i_size_ml = (int)ml.size();
	if (i_size_ol == 0) { cout << "Error! no observed unit in FHDI_Extension.." << endl; return; }
	if (i_size_ml == 0) { cout << "Error! no missing  unit in FHDI_Extension.." << endl; return; }

	//--------------
	//Rows of observed RAW data 
	//Rows of missing  RAW data
	//--------------
	double** d_oy = New_dMatrix(i_size_ol, ncol);
	double** d_my = New_dMatrix(i_size_ml, ncol);
	for (int i = 0; i<i_size_ol; i++)
	{
		for (int j = 0; j<ncol; j++)
		{
			d_oy[i][j] = y[ol[i] - 1][j]; //-1 for Actual loc
		}
	}
	for (int i = 0; i<i_size_ml; i++)
	{
		for (int j = 0; j<ncol; j++)
		{
			d_my[i][j] = y[ml[i] - 1][j];	//-1 for Actual loc
		}
	}

	//--------------
	//Rows of observed data 
	//Rows of missing data
	//--------------
	double** d_ox = New_dMatrix(i_size_ol, ncol);
	double** d_mx = New_dMatrix(i_size_ml, ncol);
	for (int i = 0; i<i_size_ol; i++)
	{
		for (int j = 0; j<ncol; j++)
		{
			d_ox[i][j] = z[ol[i] - 1][j]; //-1 for Actual loc
		}
	}
	for (int i = 0; i<i_size_ml; i++)
	{
		for (int j = 0; j<ncol; j++)
		{
			d_mx[i][j] = z[ml[i] - 1][j];	//-1 for Actual loc
		}
	}

	//----------------
	//weights corresponding to missing/observed rows. 
	//select out weights at missing rows and observed rows
	//----------------
	double* w_ml = new double[i_size_ml]; //same as mw
	double* w_ol = new double[i_size_ol]; //same as ow
	for (int i = 0; i<i_size_ml; i++) w_ml[i] = w[ml[i] - 1]; //-1 for actual loc
	for (int i = 0; i<i_size_ol; i++) w_ol[i] = w[ol[i] - 1]; //-1 for actual loc

															  //----------------
															  //index corresponding to missing/observed rows. 
															  //select out weights at missing rows and observed rows
															  //----------------
	int* id_ml = new int[i_size_ml]; //same as mid
	int* id_ol = new int[i_size_ol]; //same as oid
	for (int i = 0; i<i_size_ml; i++) id_ml[i] = id[ml[i] - 1]; //-1 for actual loc
	for (int i = 0; i<i_size_ol; i++) id_ol[i] = id[ol[i] - 1]; //-1 for actual loc

																//cout << "FHDI_startup Running time is " << MPI_Wtime() - FHDI_begin << " at node " << mynode << endl;
																//-------------------------------
																//Step 1: generate uox and mox
																//-------------------------------
																//make UNIQUE patterns of z by cn
																//--------------
																//transform z into condensed string format
																//--------------
																//std::string cn[nrow]; //declaration of concatenated vector of z
	double FHDI_uoxmox = MPI_Wtime();

	std::string *cn = new std::string[nrow]; //declaration of concatenated vector of z
	Trans(z, nrow, ncol, cn);

	//---------------
	//Rows of Condensed Strings with Observed cells
	//                          with Missing  cells
	//---------------
	//std::string s_ocn[i_size_ol];
	//std::string s_mcn[i_size_ml];
	std::string *s_ocn = new std::string[i_size_ol];
	std::string *s_mcn = new std::string[i_size_ml];
	for (int i = 0; i<i_size_ol; i++) s_ocn[i] = cn[ol[i] - 1]; //-1 for actual row
	for (int i = 0; i<i_size_ml; i++) s_mcn[i] = cn[ml[i] - 1]; //-1 for actual row

																//std::string s_ocn_temp[i_size_ol]; //string vector of observed patterns only
																//std::string s_mcn_temp[i_size_ml]; //string vector of missing patterns only
	std::string *s_ocn_temp = new std::string[i_size_ol]; //string vector of observed patterns only
	std::string *s_mcn_temp = new std::string[i_size_ml]; //string vector of missing patterns only	
	for (int i = 0; i<i_size_ol; i++) { s_ocn_temp[i] = s_ocn[i]; }
	for (int i = 0; i<i_size_ml; i++) { s_mcn_temp[i] = s_mcn[i]; }
	//sort 		
	std::sort(s_ocn_temp, s_ocn_temp + i_size_ol); //knowing that s_ocn_temp[] has i_size_ol entities
	std::sort(s_mcn_temp, s_mcn_temp + i_size_ml); //knowing that s_mcn_temp[] has i_size_ml entities
												   //cout << "0.5M debug1 at node " << mynode << endl;
												   //------------
												   //memorize observed patterns. Only unique patterns are stored  
												   //------------
	double** uox = New_dMatrix(nrow, ncol);
	double** mox = New_dMatrix(nrow, ncol);
	if (i_size_ol == 0 || i_size_ol >= nrow) {
		cout << "Error!!!! " << endl;
	}
	//cout << "total i_size_ol: " << i_size_ol << endl;
	int i_count_uox = 0; //total number of unique uox 
	std::string s_temp;

	for (int i = 0; i < i_size_ol; i++)
	{
		s_temp = s_ocn_temp[i]; //get a string 
		if (i == 0 || (i > 0 && s_temp.compare(s_ocn_temp[i - 1]) != 0)) {
			//if (mynode == 0) {
			//	cout<<"III: "<<i<<endl;
			//}
			for (int j = 0;j < nrow;j++) {
				if (s_temp.compare(cn[j]) == 0) {
					for (int k = 0; k < ncol; k++)
					{
						uox[i_count_uox][k] = z[j][k];
					} //store the found observed pattern
					i_count_uox++;
					break;
				}
			}
		}

	}

	//Now, i_count_uox means the total number of unique observed patterns
	//cout << "0.5M debug2 at node " << mynode << endl;
	//cout<<"Debug i_count_uox = "<< i_count_uox <<endl;

	//if (mynode == 0) {
	//	TestOut_Slave1 << "Debug_UOX" << endl;
	//	for (int k = 0;k < i_count_uox;k++) {
	//		for (int j = 0; j < ncol; j++)
	//		{
	//			TestOut_Slave1 << k << " : " << uox[k][j] << " , ";
	//		}
	//		TestOut_Slave1 << endl;
	//	}
	//}

	//------------
	//memorize missing patterns 
	//------------
	if (i_size_ml == 0 || i_size_ml >= nrow) {
		cout << "Error!!!! " << endl;
	}
	int i_count_mox = 0; //total number of unique mox 
						 //cout<<"total i_size_ml: "<< i_size_ml <<endl;

	for (int i = 0; i < i_size_ml; i++)
	{
		s_temp = s_mcn_temp[i]; //get a string 
		if (i == 0 || (i > 0 && s_temp.compare(s_mcn_temp[i - 1]) != 0)) {
			//if (mynode == 0) {
			//	cout << "BBB: " << i << endl;
			//}
			for (int j = 0;j < nrow;j++) {
				if (s_temp.compare(cn[j]) == 0) {
					for (int k = 0; k < ncol; k++)
					{
						mox[i_count_mox][k] = z[j][k];
					} //store the found observed pattern
					i_count_mox++;
					break;
				}
			}
		}

	}

	//cout << "Debug i_count_mox = " << i_count_mox << endl;
	//if (mynode == 0) {
	//	TestOut_Slave1 << "Debug_MOX" << endl;
	//	for (int k = 0;k < i_count_mox;k++) {
	//		for (int j = 0; j < ncol; j++)
	//		{
	//			TestOut_Slave1 << k << " : " << mox[k][j] << " , ";
	//		}
	//		TestOut_Slave1 << endl;
	//	}
	//}
	//Now, i_count_mox means the total number of unique missing patterns
	//cout << "0.5M debug3 at node " << mynode << endl;
	//----------------
	//additional check for unique observed and missing patterns
	//----------------
	//observed patterns//
	d_temp = 0.0;
	double** uox_final = New_dMatrix(nrow, ncol);
	for (int j = 0; j<ncol; j++) { uox_final[0][j] = uox[0][j]; } //first row initialization
	int i_count_uox_final = 1; //starting from the second row

	for (int i = 1; i<i_count_uox; i++) //starting from the second one
	{
		d_temp = 0.0; //initialize 
		for (int j = 0; j<ncol; j++) { d_temp += fabs(uox[i][j] - uox[i - 1][j]); } //difference of adjacent rows

		if (d_temp > 1e-3) //adjacent rows are NOT the same each other
		{
			for (int j = 0; j<ncol; j++) { uox_final[i_count_uox_final][j] = uox[i][j]; }
			i_count_uox_final++;
		}
	}
	i_count_uox = i_count_uox_final; //replace with the accurate value
									 //store the final matrix 
	for (int i = 0; i<i_count_uox; i++)
	{
		for (int j = 0; j<ncol; j++)
			uox[i][j] = uox_final[i][j];
	}
	Del_dMatrix(uox_final, nrow, ncol);
	//cout << "0.5M debug4 at node " << mynode << endl;
	//--------------------------
	//missing patterns//
	//--------------------------
	double** mox_final = New_dMatrix(nrow, ncol);
	for (int j = 0; j<ncol; j++) { mox_final[0][j] = mox[0][j]; } //first row initialization
	int i_count_mox_final = 1; //starting from the second row

	for (int i = 1; i<i_count_mox; i++) //starting from the second one
	{
		d_temp = 0.0; //initialize
		for (int j = 0; j<ncol; j++) { d_temp += fabs(mox[i][j] - mox[i - 1][j]); } //difference of adjacent rows

		if (d_temp > 1e-3) //adjacent rows are NOT the same each other
		{
			for (int j = 0; j<ncol; j++) { mox_final[i_count_mox_final][j] = mox[i][j]; }
			i_count_mox_final++;
		}
	}
	i_count_mox = i_count_mox_final; //replace with the accurate value

									 //store the final matrix	
	for (int i = 0; i<i_count_mox; i++)
	{
		for (int j = 0; j<ncol; j++)
			mox[i][j] = mox_final[i][j];
	}
	Del_dMatrix(mox_final, nrow, ncol);
	//cout << "0.5M debug5 at node " << mynode << endl;
	//!!!!! now uox and mox have the UNIQUE observed and missing patterns
	//!!!!! i_count_mox and _uox have the final number of meaningful rows of mox and uox, respectively
	const int nrm = i_count_mox;
	const int nru = i_count_uox;

	//-------------------------------------------
	//-------------------------------------------
	//Step 2: Impute missing cells
	//        using all possible donors per missing pattern
	//-------------------------------------------
	//-------------------------------------------
	int* i_temp_x = new int[ncol];
	int i_sum_x = 0;

	std::vector<int> v_mxl; //a row's columns having observed cells  
	rbind_FHDI rbind_icell(ncol); //all possible donor cells 

	std::vector<int> v_cn_z_i;
	//int* zid = NULL;
	//int i_size_zid=0; 
	//int i_loc=0;	
	int* i_srst = new int[nru];
	std::vector<int> loc_srst_nl;
	double* d_temp_cn = new double[ncol];

	//List_FHDI List_ord(nrm); //order records used for variance estimation
	//List_FHDI List_ocsg(nrm); //order records used for variance estimation

	//--------------------------------
	//--------------------------------
	//--------------------------------
	//Main Loop for FEFI and FHDI
	//--------------------------------
	//--------------------------------
	//--------------------------------
	rbind_FHDI rbind_imat_FEFI(7 + 2 * ncol); //large storage that will accumulate fmat from FEFI  
	rbind_FHDI rbind_imat_FHDI(7 + 2 * ncol); //large storage that will accumulate fmat from FHDI
											  //cout << "FHDI_uoxmox Successful Running time is " << MPI_Wtime() - FHDI_uoxmox << " at node " << mynode << endl;

	double FHDI_main = MPI_Wtime();
	for (int i = 0; i<nrm; i++)
	{

		//cout << "NRM: " << i <<" at node "<<mynode<<endl; 

		//get current row of missing cell 
		for (int j = 0; j<ncol; j++) i_temp_x[j] = (int)mox[i][j];
		i_sum_x = sum_FHDI(i_temp_x, ncol);
		//cout<<"i_sum_x is "<< i_sum_x <<" at nrm "<<i<<" in mynode "<<mynode<<endl;

		//-------
		//re-initialization for this missing row 
		//-------
		rbind_icell.initialize(ncol);

		//----------------------
		//Condition 1: this row's cells are all missing
		//-----------------------
		if (i_sum_x == 0)
		{
			v_mxl.clear(); //no missing cells  
			rbind_icell.bind_blocks(i_count_uox, ncol, uox); //fine due to row-based copy
		}

		//cout << "nrmdebug_1111 at nrm of " << i <<" at node "<<mynode<< endl;
		//----------------------
		//Condition 2: this row's cells are partly missing
		//-----------------------
		int nl = 0;
		if (i_sum_x > 0)
		{
			//------
			//number of observed cells on this row
			//------
			nl = 0;
			v_mxl.clear();
			for (int j = 0; j<ncol; j++)
			{
				if (mox[i][j]>0)
				{
					nl++;
					//v_mxl.push_back(j + 1); //Actual non-missing cell location 
				}
			}

			if (nl > i_collapsing) {
				nl = i_collapsing;
			}

			//cout << "nrmdebug_1112 at nrm of " << i << " at node " << mynode << endl;
			//-------
			//indicator matrix that matches the donors
			//srst: row-wise sum of the indicator matrix 
			//-------
			loc_srst_nl.clear(); //re-initialize
			Fill_iVector(i_srst, nru, 0); //re-initialize 

										  //inherents the most correlated variables of mox[i]

			for (int k = 0; k < i_collapsing; k++) {
				if (codes[i][k] != 0) {
					v_mxl.push_back(codes[i][k]);
					//TestOut << "code[" << i << "]: " << codes[i][k] << endl;
				}
			}

			int v_mxl_size = v_mxl.size();

			for (int j = 0; j<nru; j++)
			{
				int i_sum_crst = 0;
				for (int k = 0; k<v_mxl_size; k++)
				{
					//Note: in below check, mox is fixed at ith row 
					if (fabs(mox[i][v_mxl[k] - 1] - uox[j][v_mxl[k] - 1])<1e-3) //part of missing cell = obserbed cell 
					{
						i_sum_crst++; // increment if a cell of missing row = obs. cell 
					}
				}
				//---
				//store how many cells of missing row match those of observed row
				//---
				i_srst[j] = i_sum_crst;
				if (i_sum_crst == nl) loc_srst_nl.push_back(j + 1); //Actual location 				
			}

			//cout << "nrmdebug_1 at nrm of " << i << " at mynode "<<mynode<<endl;
			//-----
			//total matching rows
			//-----
			const int i_size_loc_srst_nl = (int)loc_srst_nl.size();
			if (i_size_loc_srst_nl == 0) //error case
			{
				cout << "Error! there is no matched cell!" << endl; return;
			}

			if (i_size_loc_srst_nl > 0)
			{
				double* d_temp_srst = new double[ncol];
				for (int j = 0; j<i_size_loc_srst_nl; j++)
				{
					for (int k = 0; k<ncol; k++)
					{
						d_temp_srst[k] = uox[loc_srst_nl[j] - 1][k];
					}//-1 for actual loc

					rbind_icell.append_block(d_temp_srst); //ncol is the same 
				}
				delete[] d_temp_srst;
			}
		}


		//cout << "nrmdebug_2 at nrm of " << i <<" at mynode "<<mynode<< endl;

		//testout
		//RPrint(" ==after getting icell at i: "); RPrint(i);
		//RPrint("nl: "); RPrint(nl);
		//RPrint("rbind_icell: "); rbind_icell.print_rbind_FHDI();


		//----------------------------
		//step 3: Assign donors to missing cell
		//----------------------------
		const int nic = rbind_icell.size_row(); //get the number of total rows
												//std::string s_icn[nic];
		std::string * s_icn = new std::string[nic];
		double* d_temp_icell = new double[ncol];
		std::string s_icn_temp;
		for (int j = 0; j<nic; j++)
		{
			for (int k = 0; k<ncol; k++) d_temp_icell[k] = rbind_icell(j, k); //get jth row
			Trans1(d_temp_icell, ncol, s_icn_temp); //transform one row
			s_icn[j] = s_icn_temp;
		}
		delete[] d_temp_icell;

		//------------------
		//search locations where icn = name of jp_name
		//------------------
		const int i_size_jp_name = (int)jp_name.size();
		std::vector<double>      v_cp;  //selected joint probability 
		std::vector<std::string> v_ncp; //names of the selected joint probability
		v_cp.clear();
		v_ncp.clear();

		//cout<<"nic: "<<nic<<" at nrm "<<i<<" at node "<<mynode<<endl;
		//cout << "i_size_jp_name: " << i_size_jp_name << " at nrm " << i << " at node " << mynode << endl;

		for (int j = 0; j<nic; j++)
		{
			s_temp = s_icn[j]; //one donor 
			for (int k = 0; k<i_size_jp_name; k++) //search all names of jp
			{
				if (s_temp.compare(jp_name[k]) == 0) //0 means the same string
				{
					v_cp.push_back(jp_prob[k]); //store the joint probability 
					v_ncp.push_back(jp_name[k]); //store name
					break; //stop searching after finding the first match 
				}
			}
		}
		const int i_size_v_cp = (int)v_cp.size();
		//cout<<"i_size_v_cp: "<< i_size_v_cp <<" in nrm "<<i<<" at mynode "<<mynode<<endl;
		double d_sum_v_cp = 0.0;
		for (int j = 0; j<i_size_v_cp; j++) d_sum_v_cp += v_cp[j];
		if (d_sum_v_cp != 0)
		{
			for (int j = 0; j<i_size_v_cp; j++) v_cp[j] = v_cp[j] / d_sum_v_cp;
		}


		//cout << "nrmdebug_3 at nrm of " << i << " at mynode "<<mynode<<endl;

		//-----------
		//transform current missing row to string for step 3
		//-----------
		for (int j = 0; j<ncol; j++) d_temp_cn[j] = mox[i][j];
		Trans1(d_temp_cn, ncol, s_temp);
		v_cn_z_i.clear(); //re-initialize 
		which(cn, nrow, s_temp, v_cn_z_i); //Note: Actual location is returned
		int i_size_v_cn_z_i = (int)v_cn_z_i.size(); //number of locations in cn having s_temp
		const int mg = i_size_v_cn_z_i; //Note: loc2 =  i_size_v_cn_z_i

										//testout
										//RPrint("======= in Step 3 ====");
										//RPrint("loc2 :"); RPrint(v_cn_z_i);
										//RPrint("mg   :"); RPrint(mg);
										//RPrint("s_icn   :"); RPrint(s_icn, nic);

										//---------------------
										//select out all cells that have s_icn[1:nic] from cn
										//---------------------
		std::vector<int> v_obsg0; v_obsg0.clear();
		for (int j = 0; j<nic; j++)
		{
			s_temp = s_icn[j];
			for (int k = 0; k<nrow; k++)
			{
				if (s_temp.compare(cn[k]) == 0) //0=equal string
				{
					v_obsg0.push_back(k + 1);// Actual location of all donors of ith mox
				}//ACTUAL location stored. No exit 
			}
		}
		const int i_size_v_obsg0 = (int)v_obsg0.size();
		//cout << " i_size_v_obsg0: " << i_size_v_obsg0 <<" at nrm of " << i << " at mynode " << mynode << endl;
		//---------------
		//sort the found cells
		//---------------
		int* i_obsg0_sorted = new int[i_size_v_obsg0];
		for (int k = 0; k<i_size_v_obsg0; k++) i_obsg0_sorted[k] = v_obsg0[k];
		std::sort(i_obsg0_sorted, i_obsg0_sorted + i_size_v_obsg0);
		for (int k = 0; k<i_size_v_obsg0; k++) v_obsg0[k] = i_obsg0_sorted[k];// sort Actual location of all donors of ith mox
		delete[] i_obsg0_sorted;

		//cout << "nrmdebug_4 at nrm of " << i << " at mynode " << mynode << endl;

		//testout
		//RPrint("v_obsg0 :"); RPrint(v_obsg0);
		//RPrint("=========== begin yorder()  with d_temp_cn :"); RPrint(d_temp_cn, ncol);

		//------------------
		//half-ascending and -descending ordering 
		//------------------
		std::vector<int> v_obsg; //half-asc and desc obsg0

		int* i_ym_return = new int[i_size_v_obsg0]; //half-asc and desc 
		yorder(y, nrow, ncol,
			d_temp_cn,
			v_obsg0, i_ym_return);

		v_obsg.clear();
		for (int j = 0; j<i_size_v_obsg0; j++) v_obsg.push_back(i_ym_return[j]);
		//below is temporary for wrong yorder()
		//for(int j=0; j<i_size_v_obsg0; j++) v_obsg.push_back(v_obsg0[j]); 

		const int i_size_v_obsg = (int)v_obsg.size();
		delete[] i_ym_return;

		//cout << "i_size_v_obsg at nrm of " << i_size_v_obsg << " at mynode " << mynode << endl;
		//testout
		//RPrint("============= after yorder() ================");
		//RPrint("v_obsg0 :"); RPrint(v_obsg0);
		//RPrint("v_obsg :");  RPrint(v_obsg);

		//-------------
		// find positions of matches between obsg in obsg0
		// then Store them into List 
		//-------------
		std::vector<int> v_rbsg; v_rbsg.clear();
		match_FHDI(v_obsg, v_obsg0, v_rbsg); //get loc stored in v_rbsg
		const int i_size_v_rbsg = (int)v_rbsg.size();
		//store rbsg
		double* d_temp_rbsg = new double[i_size_v_rbsg];
		for (int k = 0; k<i_size_v_rbsg; k++) d_temp_rbsg[k] = v_rbsg[k];
		List_ord.put_block(i, i_size_v_rbsg, d_temp_rbsg); //put into storage as ith row
		delete[] d_temp_rbsg;
		//store obsg
		double* d_temp_obsg = new double[i_size_v_obsg];
		for (int k = 0; k<i_size_v_obsg; k++) d_temp_obsg[k] = v_obsg[k];
		List_ocsg.put_block(i, i_size_v_obsg, d_temp_obsg); //put into storage as ith row
		delete[] d_temp_obsg;

		const int ng = i_size_v_obsg;
		//testout
		//RPrint("v_rbsg :"); RPrint(v_rbsg);
		//RPrint("ord :");   List_ord.print_one_List_FHDI(i); 
		//RPrint("ocsg :");  List_ocsg.print_one_List_FHDI(i); 
		//RPrint("ng   :"); RPrint(ng);

		//------------------------------
		//------------------------------
		//Compute Fractional Weights (fwij)
		//Fractional weights for FEFI representing sampling w
		//------------------------------
		//------------------------------
		//std::string cn_obsg[i_size_v_obsg]; //cn at locations of obsg
		std::string * cn_obsg = new std::string[i_size_v_obsg]; //cn at locations of obsg
		for (int k = 0; k<i_size_v_obsg; k++)
		{
			cn_obsg[k] = cn[v_obsg[k] - 1];
		}//-1 for actual location 
		 //cout << "nrmdebug_5_1 at nrm of " << i << " at mynode " << mynode << endl;

		 //-----
		 //make a table of cn at obsg locations
		 //------
		std::vector<std::string> v_table_name_cn_obsg;
		std::vector<int>         v_table_count_cn_obsg;

		table_cpp(cn_obsg, i_size_v_obsg,
			v_table_name_cn_obsg, v_table_count_cn_obsg);
		//const int i_size_table_cn_obsg = (int)v_table_count_cn_obsg.size(); 

		//------
		//get joint probability of the selected donor 
		//------
		std::vector<int> v_cn_obsg_ncp; //positions of cn_obsg in ncp 
		match_FHDI(cn_obsg, i_size_v_obsg, v_ncp,
			v_cn_obsg_ncp); //Note: Actual locations are returned 
		const int i_size_v_cn_obsg_ncp = (int)v_cn_obsg_ncp.size();

		//------
		//calculate fractional weights
		//------
		double* fwij = new double[i_size_v_cn_obsg_ncp]; //fractional weights 
		double* d_obsp = new double[i_size_v_cn_obsg_ncp]; //joint prob of selected donors
		int*    i_obsn = new int[i_size_v_cn_obsg_ncp]; //counts of the selected donors
		for (int k = 0; k<i_size_v_cn_obsg_ncp; k++)
		{
			d_obsp[k] = v_cp[v_cn_obsg_ncp[k] - 1];// -1 for actual location 
			i_obsn[k] = v_table_count_cn_obsg[v_cn_obsg_ncp[k] - 1];// -1 for actual location  

			fwij[k] = 1.0; //default for error case  
			if (i_obsn[k] != 0) fwij[k] = d_obsp[k] / i_obsn[k];
			if (i_obsn[k] == 0) cout << "Error! zero count in obsn!" << endl;
		}
		//testout
		//RPrint("fwij[] :"); RPrint(fwij, i_size_v_cn_obsg_ncp);
		//cout << "nrmdebug_6 at nrm of " << i << " at mynode " << mynode << endl;

		//===============================================
		// Note that v_mxl hereafter this function indicates locations of observed variables of mox[i], similar functionality of response indicator, not selected variables for donors. Written by Yicheng Yang
		//===============================================

		v_mxl.clear();
		for (int j = 0; j < ncol; j++)

		{

			if (mox[i][j] > 0)

			{

				v_mxl.push_back(j + 1); //Actual non-missing cell location 

			}

		}

		//----------------------
		//FEFI Imputation
		//  Algorithm: impute by using all possible donors
		//  final outcome is "fmat" in which each column means that
		//  col1: id
		//  col2: fid, i.e., id of imputed value
		//  col3: sampling weight
		//  col4: fractional weights 
		//  col5: imputed original data (matrix with column of ncol) 
		//  col6: imputed category data (matrix with column of ncol)
		//  col7: 1:ng 
		//  col8: = col2  (for consistency with FHDI results)
		//  col9: = col3  (for consistency with FHDI results)
		//----------------------
		std::vector<int> v_obsg_times_mg; v_obsg_times_mg.clear();
		if (s_M.compare("FEFI") == 0) //0=equal string
		{
			double** fmat_FEFI = New_dMatrix(ng*mg, 7 + 2 * ncol); //7columns and two blocks of ncol 

			Fully_Efficient_Fractional_Imputation(ng, mg,
				v_obsg, v_mxl,
				y, z, nrow, ncol,
				v_cn_z_i, fwij, i_size_v_cn_obsg_ncp,
				w, id,
				fmat_FEFI);

			//testout
			//RPrint("in ==== M=FEFI =after making fmat_FEFI[][]====");
			//RPrint("fmat_FEFI : "); RPrint(fmat_FEFI, ng*mg, 7+2*ncol);

			//------------------
			//Append fmat_FEFI onto global storage imat
			//------------------
			rbind_imat_FEFI.bind_blocks(ng*mg, 7 + 2 * ncol, fmat_FEFI);

			//-----------------------
			//local deallocation
			//-----------------------
			Del_dMatrix(fmat_FEFI, ng*mg, 7 + 2 * ncol);
		}


		//------------------------------------
		//FHDI
		//Fractional Hot Deck Imputation
		//------------------------------------
		//if(s_M.compare("FHDI") == 0) //0= equal string

		const int i_mxl = (int)v_mxl.size();

		//prepare return matrix. Note the different row size from fmat of FEFI 
		int i_row_fmat_FHDI = i_M*mg; //default row size of return matrix of FHDI
		if (i_size_v_obsg <= i_M) i_row_fmat_FHDI = i_size_v_obsg*mg;  //if donors are less than i_M

		if (s_M.compare("FHDI") == 0) //0= equal string
		{
			double** fmat_FHDI = New_dMatrix(i_row_fmat_FHDI, 7 + 2 * ncol); //return matrix from FHDI

			d_myran = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);

			Fractional_Hot_Deck_Imputation(i,
				ng, List_ocsg, ncol,
				mox, y, nrow, i_M,
				mg, z, i_mxl,
				v_cn_z_i, v_mxl,
				v_obsg,
				fwij, i_size_v_cn_obsg_ncp,
				d_obsp, i_obsn,
				d_myran,
				w, id,
				fmat_FHDI);

			//------------------
			//Append fmat_FHDI onto global storage imat
			//------------------
			rbind_imat_FHDI.bind_blocks(i_row_fmat_FHDI, 7 + 2 * ncol, fmat_FHDI);

			//-----------------------
			//local deallocation
			//-----------------------
			Del_dMatrix(fmat_FHDI, i_row_fmat_FHDI, 7 + 2 * ncol);
		}
		//if (mynode == 0) {
		//	cout << "nrmdebug_7 at nrm of " << i << endl;
		//}
		//-----------------------
		//local deallocation
		//-----------------------
		delete[] s_icn;
		delete[] cn_obsg;
		delete[] fwij;
		delete[] d_obsp;
		delete[] i_obsn;

	} //end of Main loop for all rows of missing patterns 
	  //if (mynode == 0) cout << "FHDI_main Running time is " << MPI_Wtime() - FHDI_main << endl;
	  //testout
	  //RPrint("========= at the end of FHDI Extension ==============");
	  //RPrint("acculated imat_FEFI:"); rbind_imat_FEFI.print_rbind_FHDI();
	  //RPrint("acculated imat_FHDI:"); rbind_imat_FHDI.print_rbind_FHDI();

	  //---------------
	  //---------------
	  //Step 4: construct output results
	  //---------------
	  //------------------------------------------------------
	  //ipmat  = final imputation results
	  //     	col1: ID 	= unit index
	  //		col2: FID 	= ID of fractionally imputed value
	  // 		col3: WGT 	= weight 
	  //		col4: FWGT	= Frational weight
	  //		col5: Variables 
	  //		col6: Responses
	  //irmat  = imputation results related to the categorized matrix 
	  //     	col1: ID 	= unit index
	  //		col2: FID 	= ID of fractionally imputed value
	  //		col3: OID	= original rank of the imputed value
	  //		col4: ORDER = SN(selected donor)
	  //		col5: FEFIW	= Fefi weights 
	  //		col6: CELL	= cells 
	  //----------------------------------------------------
	  //FEFI                           FEFI //
	  //get ipmat, Resp (separately), irmat from FEFI results 
	  //rbind_FHDI  rbind_ipmat_FEFI(4+ncol); //column size is 4+ncol
	  //rbind_FHDI  rbind_Resp_FEFI(ncol+1);  //separate response matrix  
	  //rbind_FHDI  rbind_irmat_FEFI(5+ncol); //column size is 5+ncol
	double FHDI_FEFI = MPI_Wtime();
	if (s_M.compare("FEFI") == 0) //0= equal string
	{
		Results_Fully_Efficient_Fractional_Imputation(i_size_ol,
			ncol, nrow,
			id_ol, w_ol, d_oy, d_ox,
			rbind_imat_FEFI, r,

			rbind_ipmat_FEFI, rbind_Resp_FEFI, rbind_irmat_FEFI);
		//testout
		/*
		RPrint("after Results_... rbind_ipmat_FEFI after binding :");
		rbind_ipmat_FEFI.print_rbind_FHDI();
		RPrint("after Results_... rbind_Resp_FEFI after binding :");
		rbind_Resp_FEFI.print_rbind_FHDI();
		RPrint("after Results_... rbind_irmat_FEFI after binding :");
		rbind_irmat_FEFI.print_rbind_FHDI();
		*/
	}
	//cout << "0.5M debug6 at node " << mynode << endl;
	//FHDI ------------------------- FHDI //
	//get ipmat, Resp (separately), irmat from FHDI results 
	//rbind_FHDI  rbind_ipmat_FHDI(4+ncol); //column size is 4+ncol
	//rbind_FHDI  rbind_Resp_FHDI(ncol+1);  //separate response matrix  
	//rbind_FHDI  rbind_irmat_FHDI(5+ncol); //column size is 5+ncol

	if (s_M.compare("FHDI") == 0) //0= equal string
	{
		Results_Fractional_Hot_Deck_Imputation(i_size_ol,
			ncol, nrow,
			id_ol, w_ol, d_oy, d_ox,
			rbind_imat_FHDI, r,

			rbind_ipmat_FHDI, rbind_Resp_FHDI, rbind_irmat_FHDI);
		//testout
		/*
		RPrint("after Results_... rbind_ipmat_FHDI after binding :");
		rbind_ipmat_FHDI.print_rbind_FHDI();
		RPrint("after Results_... rbind_Resp_FHDI after binding :");
		rbind_Resp_FHDI.print_rbind_FHDI();
		RPrint("after Results_... rbind_irmat_FHDI after binding :");
		rbind_irmat_FHDI.print_rbind_FHDI();
		*/
	}
	//cout << "0.5M debug7 at node " << mynode << endl;
	//if (mynode == 0) cout << "FHDI_FEFI Running time is " << MPI_Wtime() - FHDI_FEFI << endl;
	//------
	//prep returns of other matrices
	//------
	rbind_uox.bind_blocks(i_count_uox, ncol, uox);
	rbind_mox.bind_blocks(i_count_mox, ncol, mox);

	//testout
	//cout << " ========= FHDI_Extension.. has successfully finished!" << endl;


	//----------------
	//Deallocation
	//----------------
	delete[] cn;
	delete[] s_ocn;
	delete[] s_mcn;
	delete[] s_ocn_temp;
	delete[] s_mcn_temp;

	delete[] i_rn;
	delete[] w;
	delete[] id;
	Del_dMatrix(d_oy, i_size_ol, ncol);
	Del_dMatrix(d_my, i_size_ml, ncol);
	Del_dMatrix(d_ox, i_size_ol, ncol);
	Del_dMatrix(d_mx, i_size_ml, ncol);

	delete[] w_ml;
	delete[] w_ol;
	delete[] id_ml;
	delete[] id_ol;

	Del_dMatrix(uox, nrow, ncol);
	Del_dMatrix(mox, nrow, ncol);

	delete[] i_temp_x;
	delete[] i_srst;
	delete[] d_temp_cn;


	return;

}

//void yorder(double** y, const int nrow, const int ncol,
//	double* mox_1,
//	std::vector<int> v_loc,
//	int* i_ym_return)
//	//Description=========================================
//	// order the donors in  
//	// half-ascending and half-descending manner
//	// 
//	//
//	// original R code: Dr. Im, J. and Dr. Kim, J. 
//	// c++ code: 		Dr. Cho, I. 
//	// All rights reserved
//	// 
//	// updated: Nov 07, 2016
//	//----------------------------------------------------
//	//IN    : double y(nrow, ncol)= original data matrix
//	//IN    : double mox_1(ncol)= ith row of unique missing patterns
//	//IN	: vector<int> v_loc  = ACTUAL locations of donors of ith mox
//	//OUT   : int i_ym_return(i_size_v_loc)  = index in half-ascending and -descending order
//	//====================================================
//{
//	//----------------
//	//get ready original data at columns of zero mox_i
//	//----------------
//	int mynode, totalnodes;
//	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
//	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
//	MPI_Status status;
//
//	//cout << "yorder_1 at mynode " << mynode << endl;
//	std::vector<int> rloc; //ACTUAL location of zero in mox_1
//	for (int i = 0; i<ncol; i++)
//	{
//		if (mox_1[i] == 0) rloc.push_back(i + 1);
//	}
//	const int i_size_rloc = (int)rloc.size();
//
//	double** yv = New_dMatrix(nrow, i_size_rloc);
//
//	for (int i = 0; i<nrow; i++)
//	{
//		for (int j = 0; j< i_size_rloc; j++)
//		{
//			yv[i][j] = y[i][rloc[j] - 1]; //-1 for actual location 	
//		}
//	}
//	//testout
//	//RPrint(" in yorder() yv:"); RPrint(yv, nrow, i_size_rloc);
//	//cout << "yorder_2 at mynode " << mynode << endl;
//	//---------------------
//	//extract donors at the known locations only
//	//---------------------
//	const int i_size_v_loc = (int)v_loc.size();
//	double** yl = New_dMatrix(i_size_v_loc, i_size_rloc);
//	for (int i = 0; i<i_size_v_loc; i++)
//	{
//		for (int j = 0; j<i_size_rloc; j++)
//		{
//			yl[i][j] = yv[v_loc[i] - 1][j]; //-1 for actual location 
//		}
//	}
//	//testout
//	//RPrint(" in yorder() yl:"); RPrint(yl, i_size_v_loc, i_size_rloc);
//
//	//----------------------
//	//column-wise mean calculation
//	//----------------------
//	double d_sum = 0.0;
//	double* my = new double[i_size_rloc]; //mean of each column
//	for (int j = 0; j<i_size_rloc; j++)
//	{
//		d_sum = 0.0;
//		for (int i = 0; i<i_size_v_loc; i++)
//		{
//			d_sum += yl[i][j];
//		}
//		my[j] = d_sum / i_size_v_loc;
//	}
//	//cout << "yorder_3 at mynode " << mynode << endl;
//	//----------------------
//	//(1) combine loc & yl       -> ym
//	//(2) ym(without 1st column) -> dym  
//	//----------------------
//	double** ym = New_dMatrix(i_size_v_loc, i_size_rloc + 1);
//	double** ym_sorted = New_dMatrix(i_size_v_loc, i_size_rloc + 1);
//	for (int i = 0; i<i_size_v_loc; i++)
//	{
//		ym[i][0] = v_loc[i]; //first column is loc (Actual loc)
//		for (int j = 0; j<i_size_rloc; j++) ym[i][j + 1] = yl[i][j];
//	}
//	//testout
//	//RPrint(" in yorder() ym:"); RPrint(ym, i_size_v_loc, i_size_rloc+1);
//	//------------
//	//get ready return 
//	//------------
//	for (int j = 0; j<i_size_v_loc; j++) //default
//	{
//		i_ym_return[j] = (int)ym[j][0]; //default is unsorted locations   
//	}
//
//	//cout << "i_size_v_loc: "<< i_size_v_loc <<" at mynode " << mynode << endl;
//	//-------
//	//covariance matrix of yl
//	//-------
//	double** VM = New_dMatrix(i_size_rloc, i_size_rloc);
//	double** VM_backup = New_dMatrix(i_size_rloc, i_size_rloc); //temp
//	double** dif = New_dMatrix(i_size_v_loc, i_size_rloc); //yl - mean
//	double** VM_inv = New_dMatrix(i_size_rloc, i_size_rloc); //inverse of VM
//	double** dif_T = New_dMatrix(i_size_rloc, i_size_v_loc); //transpose of dif
//	//double** mat_temp = New_dMatrix(i_size_v_loc, i_size_v_loc); //temporary matrix
//	double*  mat_temp = new double[i_size_v_loc];
//	double*  score = new double[i_size_v_loc];
//	if (i_size_v_loc > 1)
//	{
//		//----------
//		//"Estimated covariance" of yl by column-to-column method
//		//----------
//		cov_FHDI(yl, i_size_v_loc, i_size_rloc, VM);
//
//		//----------
//		//column-wise difference between yl[][] - mean(yl)
//		//----------
//		for (int j = 0; j<i_size_rloc; j++) //each column
//		{
//			for (int k = 0; k<i_size_v_loc; k++)
//			{
//				dif[k][j] = yl[k][j] - my[j];
//			}
//		}
//		//cout << "yorder_5 at mynode " << mynode << endl;
//		//----------
//		//matrix multiplication among dif * inverse(VM) * dif^T
//		//the "score" is similar to the Mahalanobis distance (MD) function
//		//  MD measures the number of stdev from a point P to the mean of distribution Deallocation
//		//  along each principal component axis. 
//		//  Unitless and scale-invariant and taking into account the data's correlations.
//		//  MD corresponds to standard Euclidean distance in the transformed space. 
//		//----------
//		Copy_dMatrix(VM, i_size_rloc, i_size_rloc, VM_backup);
//
//		bool b_success_VM = Inverse_dMatrix_FHDI(VM_backup, i_size_rloc, VM_inv); //backup is to avoid pivotting of VM
//		if (!b_success_VM) { cout<<"CAUTION! Inverse matrix of VM maybe incorrect"<<endl; }
//		for (int j = 0; j<i_size_v_loc; j++)
//		{
//			for (int k = 0; k<i_size_rloc; k++) dif_T[k][j] = dif[j][k];
//		}
//		//cout<<"i_size_rloc and i_size_v_loc are "<< i_size_rloc<<" , "<< i_size_v_loc <<endl;
//		dMatrix_Mul_AtBA_Yicheng(dif_T, i_size_rloc, i_size_v_loc,
//			VM_inv, mat_temp);
//		//diagonal terms are score
//		for (int j = 0; j<i_size_v_loc; j++) score[j] = mat_temp[j];
//
//		//if (mynode == 0) {
//		//	for (int j = 0; j < i_size_v_loc; j++) {
//		//		cout << "Score[" << j << "]:" << score[j] << endl;
//		//	}
//		//}
//		//------------
//		//order the score array
//		//------------
//		int* i_return = new int[i_size_v_loc]; //order of score actual loc
//		order_FHDI(score, i_size_v_loc, i_return);
//		//cout << "yorder_6 at mynode " << mynode << endl;
//		//-------------
//		//sorted ym
//		//-------------
//		int i_loc_ym = 0;
//		for (int j = 0; j<i_size_v_loc; j++)
//		{
//			i_loc_ym = i_return[j]; //actual location 
//			for (int k = 0; k<i_size_rloc + 1; k++)
//			{
//				ym_sorted[j][k] = ym[i_loc_ym - 1][k]; //-1 for actual loc
//			}
//		}
//
//		//-------------
//		//half ascending ym's first column
//		//-------------
//		int* i_ym_ascending = new int[i_size_v_loc];
//		int* i_ym_descending = new int[i_size_v_loc];
//		for (int j = 0; j<i_size_v_loc; j++)
//		{
//			i_ym_ascending[j] = (int)ym_sorted[j][0];
//			i_ym_descending[j] = (int)ym_sorted[i_size_v_loc - 1 - j][0];
//		}
//
//		//------------
//		//get ready return 
//		//------------
//		int i_temp = 0;
//		for (int j = 0; j<i_size_v_loc; j = j + 2) //ascending
//		{
//			i_ym_return[j] = i_ym_ascending[i_temp];
//			//below condition is necessary when odd sized vector 
//			if (j <= (i_size_v_loc - 2)) i_ym_return[j + 1] = i_ym_descending[i_temp];
//			i_temp++;
//		}
//		//cout << "yorder_7 at mynode " << mynode << endl;
//		//testout
//		/*
//		RPrint(" mean my          :"); RPrint(my, i_size_rloc);
//		RPrint(" Covariance matrix:"); RPrint(VM, i_size_rloc, i_size_rloc);
//		RPrint(" dif        matrix:"); RPrint(dif, i_size_v_loc, i_size_rloc);
//		RPrint(" score            :"); RPrint(score, i_size_v_loc);
//		//RPrint(" i_return            :"); RPrint(i_return, i_size_v_loc);
//		//RPrint(" i_ym_return            :"); RPrint(i_ym_return, i_size_v_loc);
//		//RPrint(" i_ym_ascending            :"); RPrint(i_ym_ascending, i_size_v_loc);
//		//RPrint(" i_ym_descending            :"); RPrint(i_ym_descending, i_size_v_loc);
//		*/
//
//		//-----
//		//local deallocation
//		//-----
//		delete[] i_return;
//		delete[] i_ym_ascending;
//		delete[] i_ym_descending;
//	}
//
//	//cout << "yorder_8 at mynode " << mynode << endl;
//	//----------------------
//	//Deallocation
//	//----------------------
//	Del_dMatrix(yv, nrow, i_size_rloc);
//	Del_dMatrix(yl, i_size_v_loc, i_size_rloc);
//	delete[] my;
//	Del_dMatrix(ym, i_size_v_loc, i_size_rloc + 1);
//	Del_dMatrix(ym_sorted, i_size_v_loc, i_size_rloc + 1);
//	Del_dMatrix(VM, i_size_rloc, i_size_rloc);
//	Del_dMatrix(VM_backup, i_size_rloc, i_size_rloc);
//	Del_dMatrix(dif, i_size_v_loc, i_size_rloc);
//	Del_dMatrix(VM_inv, i_size_rloc, i_size_rloc); //inverse of VM
//	Del_dMatrix(dif_T, i_size_rloc, i_size_v_loc); //transpose of dif	
//	//Del_dMatrix(mat_temp, i_size_v_loc, i_size_v_loc);
//	delete[] mat_temp;
//	delete[] score;
//
//	return;
//}