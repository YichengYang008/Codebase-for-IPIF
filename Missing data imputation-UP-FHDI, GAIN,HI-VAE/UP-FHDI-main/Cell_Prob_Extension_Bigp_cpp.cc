//#include "AGMAT_Extension_cpp_MPI.cc"
#include "Cal_W_Extension_Bigp_cpp.cc"
#include <vector>

void Cell_Prob_Extension_Bigp_cpp(double** z, const int nrow, const int ncol, const int i_collapsing,
							 std::vector<double> &jp_prob_return,
							 std::vector<std::string> &jp_name_return, 
							 double* w, int* id, int** codes,
						     ofstream& TestOut)

//Description=========================================
// make joint probability of cells with the categorized matrix z 
// where 0 means missing data
//
//
// original R code: Dr. Im, J. and Dr. Kim, J. 
// c++ code: 		Dr. Cho, I. 
// All rights reserved
// 
// updated: March 28, 2017
//----------------------------------------------------
//IN    : double z(nrow, ncol)  = catorized matrix corresponding to original matrix x
//                                initialized with 0.0 
//IN    : double w(nrow) = weights for rows (default = 1.0)
//IN    : int    id(nrow) = row id (default = sequential row numbers)
//IN    : int i_option_collapsing = choice of big-p algorithm 
//                               0= no big-p algorithms
//                              !0= perform big-p algorithms
//IN    : int codes(nrow, i_option_collapsing); // storage to record most correlated variables of mox
//OUT   : std::vector<double>      jp_prob_return  = final updated joint probability 
//OUT   : std::vector<std::string> jp_name_return  = name of final updated joint probability 
//====================================================
{
	//testout
	//TestOut << "==========Begin Cell Prob...===========" << endl;

	//-------------
	//maximum number of iterations for updating weights
	//-------------
	const int n_maximum_iteration = nrow*100; //set by user!
	

	//--------------
	//locations of missing cells (ml) and observed cells (ol)
	//Note: unlike in Cell_Make..(), std vector is used here
	//--------------
	std::vector<int> ol; //Actual row number 
	std::vector<int> ml; //Actual row number 

	double d_temp=0.0; 
	for(int i_row=0; i_row<nrow; i_row++)
	{
		d_temp=1.0; 
		for(int i_col=0; i_col<ncol; i_col++)
		{
			if(z[i_row][i_col] == 0) {d_temp=0.0; break;} //found zero, i.e. missing cell
		}
		
		if(fabs(d_temp) > 1e-15 ) //this row has no missing cells
		{ol.push_back(i_row + 1);} //actual number of the row having no missing cells
		
		if(fabs(d_temp) < 1e-15) //this row has AT LEAST one missing cells
		{ml.push_back(i_row + 1);}  //actual number of the row having missing cells
	}
	const int i_size_ol = (int)ol.size(); 
	const int i_size_ml = (int)ml.size(); 
	if(i_size_ol ==0) {TestOut<<"Error! no observed unit in Cell_Prov.."<<endl; return; }
	if(i_size_ml ==0) {TestOut<<"Error! no missing  unit in Cell_Prov.."<<endl; return; }

	//----------------
	//weights corresponding to missing rows. Will be used for Cal_W..() later
	//select out weights at missing rows
	//----------------
	double* w_ml = new double[i_size_ml];
	for(int i=0; i<i_size_ml; i++) w_ml[i]  = w[ml[i] - 1] ; //-1 for actual loc
	
	//--------------
	//Rows having only observed data (categorized) 
	//Rows having AT LEAST one missing data (categorized)
	//--------------
	double** d_ox = New_dMatrix(i_size_ol, ncol);
	double** d_mx = New_dMatrix(i_size_ml, ncol);
	for(int i=0; i<i_size_ol; i++) 
	{
		for(int j=0; j<ncol; j++) 
		{
			d_ox[i][j] = z[ol[i]-1][j]; //-1 for Actual loc
		}
	}
	for(int i=0; i<i_size_ml; i++) 
	{
		for(int j=0; j<ncol; j++) 
		{
			d_mx[i][j] = z[ml[i]-1][j];	//-1 for Actual loc
		}
	}
	
	//--------------
	//transform z into condensed string format
	//--------------
	//std::string cn[nrow]; //declaration of concatenated vector of z
	//std::string cn0[nrow]; //backup of cn
	std::string *cn = new std::string[nrow]; //declaration of concatenated vector of z
	std::string *cn0 = new std::string[nrow]; //backup of cn	
	Trans(z, nrow, ncol, cn);
	for(int i=0; i<nrow; i++) cn0[i] = cn[i]; 
	
	//---------------
	//Rows of Condensed Strings with Observed cells
	//                          with AT LEAT One Missing  cell
	//---------------
	//std::string s_ocn[i_size_ol];
	//std::string s_mcn[i_size_ml];
	std::string *s_ocn = new std::string[i_size_ol];
	std::string *s_mcn = new std::string[i_size_ml];	
	for(int i=0; i<i_size_ol; i++) s_ocn[i] = cn[ol[i]-1]; //-1 for actual row
	for(int i=0; i<i_size_ml; i++) s_mcn[i] = cn[ml[i]-1]; //-1 for actual row

	//---------------------
	//---------------------
	//make UNIQUE patterns of z by cn
	//i.e., uox and mox
	//---------------------
	//step. Sort the "cn"
	//---------------------
	//std::string s_ocn_temp[i_size_ol]; //string vector of observed patterns only
	//std::string s_mcn_temp[i_size_ml]; //string vector of missing patterns only
	std::string *s_ocn_temp = new std::string[i_size_ol]; //string vector of observed patterns only
	std::string *s_mcn_temp = new std::string[i_size_ml]; //string vector of missing patterns only	
	for(int i=0; i<i_size_ol; i++) {s_ocn_temp[i] = s_ocn[i];} 
	for(int i=0; i<i_size_ml; i++) {s_mcn_temp[i] = s_mcn[i];} 
		
	std::sort(s_ocn_temp, s_ocn_temp+i_size_ol); //knowing that s_ocn_temp[] has i_size_ol entities
	std::sort(s_mcn_temp, s_mcn_temp+i_size_ml); //knowing that s_mcn_temp[] has i_size_ml entities
	
	//------------
	//memorize observed patterns 
	//------------
	double** uox = New_dMatrix(nrow, ncol);
	double** mox = New_dMatrix(nrow, ncol);
	
	int i_count_uox = 0; //total number of unique uox 
	std::string s_temp ; 
	//for(int i=0; i<i_size_ol; i++)
	//{
	//	s_temp = s_ocn_temp[i]; //get a string 
	//	for(int j=0; j<nrow; j++) //search all rows 
	//	{
	//		//----
	//		//below condition is needed for finding UNIQUE pattern
	//		//----
	//		//if(j==0 && s_temp == cn[j]) 
	//		//if(i==0 && s_temp == cn[j]) //with first string, find the same string in cn 
	//		if(i==0 && s_temp.compare(cn[j]) == 0) //0: equal string
	//		{
	//			for(int k=0; k<ncol; k++) 
	//			{uox[i_count_uox][k] = z[j][k]; } //store the found observed pattern
	//			i_count_uox++; 
	//			break; 
	//		}
	//		//if(j>0 && s_temp == cn[j] && s_temp != cn[j-1])
	//		//if(i>0 && s_temp == cn[j] && s_temp != s_ocn_temp[i-1]) //find UNIQUE matching 
	//		if(i>0 && s_temp.compare(cn[j]) == 0 && s_temp.compare(s_ocn_temp[i-1]) != 0) 
	//		{
	//			for(int k=0; k<ncol; k++) 
	//			{uox[i_count_uox][k] = z[j][k]; } //store the found observed pattern				
	//			i_count_uox++; 
	//			break; 
	//		}
	//	}
	//}
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
	delete[] s_ocn_temp;
	//------------
	//memorize missing patterns 
	//------------
	int i_count_mox = 0; //total number of unique mox 
	 
	//for(int i=0; i<i_size_ml; i++)
	//{
	//	s_temp = s_mcn_temp[i]; //get a string 
	//	for(int j=0; j<nrow; j++) //search all rows 
	//	{
	//		//----
	//		//below condition is needed for finding unique pattern
	//		//----
	//		//if(j==0 && s_temp == cn[j]) 
	//		//if(i==0 && s_temp == cn[j]) //with first string, find matching string in cn
	//		if(i==0 && s_temp.compare(cn[j]) == 0 ) //0: equal string 
	//		{
	//			for(int k=0; k<ncol; k++) 
	//			{mox[i_count_mox][k] = z[j][k]; } //store the found missing pattern
	//			i_count_mox++; 
	//			break; 
	//		}
	//		//if(j>0 && s_temp == cn[j] && s_temp != cn[j-1])
	//		//if(i>0 && s_temp == cn[j] && s_temp != s_mcn_temp[i-1]) //find UNIQUE matching string
	//		if(i>0 && s_temp.compare(cn[j]) == 0 && s_temp.compare(s_mcn_temp[i-1]) != 0) //0: equal
	//		{
	//			for(int k=0; k<ncol; k++) 
	//			{mox[i_count_mox][k] = z[j][k]; } //store the found missing pattern				
	//			i_count_mox++; 
	//			break;
	//		}
	//	}
	//}
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
	//Now, i_count_mox means the total number of unique missing patterns
	
	//----------------
	//additional check for unique observed and missing patterns
	//----------------
	//observed patterns//
	d_temp = 0.0; 
	double** uox_final = New_dMatrix(nrow, ncol); 
	for(int j=0; j<ncol; j++) {uox_final[0][j] = uox[0][j]; } //first row initialization
	int i_count_uox_final = 1; //starting from the second row

	for(int i=1; i<i_count_uox; i++) //starting from the second one
	{
		d_temp = 0.0; //initialize 
		for(int j=0; j<ncol; j++) {d_temp += fabs(uox[i][j] - uox[i-1][j]) ;} //difference of adjacent rows
		
		if(d_temp > 1e-3) //adjacent rows are NOT the same each other
		{
			for(int j=0; j<ncol; j++) {uox_final[i_count_uox_final][j] = uox[i][j];} 
			i_count_uox_final++; 
		}
	}
	i_count_uox = i_count_uox_final; //replace with the accurate value
	//store the final matrix 
	for(int i=0; i<i_count_uox; i++) 
	{ 
		for(int j =0; j<ncol; j++) 
			uox[i][j] = uox_final[i][j]; 
	}
	Del_dMatrix(uox_final, nrow, ncol);
	
	//--------------------------
	//missing patterns//
	//--------------------------
	double** mox_final = New_dMatrix(nrow, ncol); 
	for(int j=0; j<ncol; j++) {mox_final[0][j] = mox[0][j]; } //first row initialization
	int i_count_mox_final = 1; //starting from the second row

	for(int i=1; i<i_count_mox; i++) //starting from the second one
	{
		d_temp = 0.0; //initialize
		for(int j=0; j<ncol; j++) {d_temp += fabs(mox[i][j] - mox[i-1][j]) ;} //difference of adjacent rows
		
		if(d_temp > 1e-3) //adjacent rows are NOT the same each other
		{
			for(int j=0; j<ncol; j++) {mox_final[i_count_mox_final][j] = mox[i][j];} 
			i_count_mox_final++; 
		}
	}
	i_count_mox = i_count_mox_final; //replace with the accurate value

	//store the final matrix	
	for(int i=0; i<i_count_mox; i++) 
	{ 
		for(int j =0; j<ncol; j++) 
			mox[i][j] = mox_final[i][j]; 
	}
	Del_dMatrix(mox_final, nrow, ncol);
	
    //!!!!! now uox and mox have the UNIQUE observed and missing patterns
	//!!!!! i_count_mox and _uox have the final number of meaningful rows of mox and uox, respectively
	
	//------------------
	//sort mcn and make a table
	//------------------
	for(int i =0; i<i_size_ml; i++) s_mcn_temp[i] = s_mcn[i]; 
	std::sort(s_mcn_temp, s_mcn_temp+i_size_ml);
	
	std::vector<std::string> v_table_tmvec_row1; //names of the table
	std::vector<int> 		 v_table_tmvec_row2; //counts of the table
	table_cpp(s_mcn_temp, i_size_ml, v_table_tmvec_row1, v_table_tmvec_row2);
	
	delete[] s_mcn_temp;
	//testout
	//RPrint("============ Cell_Prob before AGMAT==============", TestOut);
	/*
	RPrint("z :"); RPrint(z, nrow, ncol);
	RPrint("ol:"); RPrint(ol);
	RPrint("ml:"); RPrint(ml);
	RPrint("ox:"); RPrint(d_ox, i_size_ol, ncol);
	RPrint("mx:"); RPrint(d_mx, i_size_ml, ncol);
	RPrint("ocn:"); RPrint(s_ocn, i_size_ol);
	RPrint("mcn:"); RPrint(s_mcn, i_size_ml);
	RPrint("uox:"); RPrint(uox, i_count_uox, ncol);
	RPrint("mox:"); RPrint(mox, i_count_mox, ncol);
	RPrint("Table of mcn. counts: "); RPrint(v_table_tmvec_row2);
	*/
	
	//-------------------
	//Augment observed cells for missing patterns
	//algorithm: 
	// for each missing pattern, find all the possible donors
	// e.g., 
	// (1) a missing row   = 000
	// 	   agmat           = all observed rows
	// (2) a missing row   = a01
	//     agmat           = ac1, af1, a11, ..., az1. 
	//-------------------
	//rbind_FHDI agmat(ncol); //Note: without the first column of id
	std::vector<std::string> agmat;

	AGMAT_Extension_Bigp_cpp(mox, i_count_mox, 
						uox, i_count_uox, 
						ncol, id, 
						v_table_tmvec_row1,
						v_table_tmvec_row2,
                        cn, nrow, i_collapsing, codes,
						agmat, TestOut); 
	//const int n_row_agmat = agmat.size_row(); //get the number of rows 
	const int n_row_agmat = agmat.size(); //get the number of rows 
	//---------
	//Translate 1. existing ox (rows with full observations) & 2. agmat
	//without appending the augmented rows onto the previous ox, i.e. the existing rows with observations
	//---------
	const int i_total_ox_agmat =  i_size_ol + n_row_agmat;
	double* d_fmat1 = new double[ncol]; //one row of observations
	std::string s_fcd1; //one row of translated string 
	//std::string s_fcd[i_total_ox_agmat]; //total rows of translated ox and augmat
	std::string *s_fcd = new std::string[i_total_ox_agmat]; //total rows of translated ox and augmat
	
	for(int i=0; i<i_total_ox_agmat; i++)
	{
		if(i<i_size_ol) //up to all ox rows 
		{
			for(int j=0; j<ncol; j++) d_fmat1[j] = d_ox[i][j];
			Trans1(d_fmat1, ncol, s_fcd1);
			s_fcd[i] = s_fcd1; 
		}
		if(i>=i_size_ol) // all augmat rows 
		{
			//Fill_dVector(d_fmat1, ncol, 0.0); //re-initialize
			//agmat.get_block( i-i_size_ol, d_fmat1); //row number of agmat = 0...n_row_agmat-1
			//Trans1(d_fmat1, ncol, s_fcd1);
			s_fcd[i] = agmat[i - i_size_ol];
		}		
	}
	
	//-------------
	//sampling weight of the observed unit
	//-------------
	double* w1 = new double[i_size_ol];
	for(int i=0; i<i_size_ol; i++) w1[i] = w[ol[i] - 1]; //-1 for actual location
	
	//testout
	/*
	RPrint("After AUGMAT =========="); 
	RPrint("agmat:"); agmat.print_rbind_FHDI(); 
	RPrint("s_fcd:"); RPrint(s_fcd, i_total_ox_agmat); 
	RPrint("w1:"); RPrint(w1, i_size_ol); 
	*/
	
	//----------------
	//calculate weighted joint probability
	//Note: Initial jp is calculated only
	//with all the OBSERVED condensed strings in s_ocn
	//NOT with Augmented data matix 
	//-----------------
	std::vector<std::string> jp_name; 
	std::vector<double>		 jp_prob;
	wpct_FHDI(s_ocn, i_size_ol, w1, jp_name, jp_prob); 
	const int i_size_jp_prob = (int)jp_prob.size(); 
               
	//testout
	//RPrint("After wpct initial ==========", TestOut);
	/*
	RPrint("After wpct initial =========="); 
	std::string s_temp_jp[i_size_jp_prob];
	for(int i=0; i<i_size_jp_prob; i++) s_temp_jp[i] = jp_name[i]; 
	RPrint("jp_name:"); RPrint(s_temp_jp, i_size_jp_prob); 
	RPrint("jp_prob:"); RPrint(jp_prob); 
	*/
	
	//===================================
	//===================================
	//Cal_W(): update new weights and the joint probability of cells
	//===================================
	//===================================
	std::vector<double> 		w20; //new storage for updated weights  
	std::vector<double>		 	jp_prob_0; //probability backup in the loop
	std::vector<std::string> 	jp_name_new; 
	std::vector<double>		 	jp_prob_new;	

	for(int j=0; j<i_size_jp_prob; j++) 
	{
		jp_name_new.push_back(jp_name[j]); //initialize with jp_prob
		jp_prob_new.push_back(jp_prob[j]); //initialize with jp_prob 
	}

	//MAIN ITERATION for Updating weights =======================
	for(int i_loop=0; i_loop<n_maximum_iteration; i_loop++)
	{
		//------------
		//intialize with the updated joint probability
		//------------
		jp_prob_0.clear(); //re-initialize 
		for(int j=0; j<i_size_jp_prob; j++) 
		{
			jp_prob_0.push_back(jp_prob_new[j]); //initialize with jp_prob_new 
		}

		//---------------------------------------
		//update weights, w20[] 
		//Note: prob must be the newest one! i.e. jp_prob_new
		//---------------------------------------
		w20.clear();  //re-initialize 
		Cal_W_Extension_Bigp_cpp(mox, i_count_mox, 
						uox, i_count_uox, i_collapsing,
						ncol, id, codes,
						v_table_tmvec_row1,
						v_table_tmvec_row2,
						jp_prob_new,
						d_mx, i_size_ml, 
						w_ml, cn0, nrow, 
						w20);
		
		//testout
		
		//RPrint(" ----------------------Cal_W.. finished at i_loop: ", TestOut); 
		//RPrint(i_loop, TestOut);					
		//RPrint("w20:", TestOut); RPrint(w20, TestOut);
		
		
		//-----------
		//combine new weights
		//-----------
		double* w12 = new double[i_total_ox_agmat]; 
		for(int j=0; j<i_total_ox_agmat; j++)
		{
			if(j<i_size_ol) //weights of existing w1
			{ w12[j] = w1[j];}
			if(j>=i_size_ol) //updated weights 
			{ w12[j] = w20[j-i_size_ol]; }
		}
	
		//-----------
		//new joint probability
		//Note: Unlike the initial jp, 
		//new augmented matrix along with the updated weight vector are used for the jp update 
		//-----------
		jp_name_new.clear(); //re-initialize
		jp_prob_new.clear(); //re-initialize
		wpct_FHDI(s_fcd, i_total_ox_agmat, w12, jp_name_new, jp_prob_new); 
		
		//------------
		//calculate difference in the joint probability
		//------------
		double dif = 0.0; 
		for(int j=0; j<i_size_jp_prob; j++) 
			dif += (jp_prob_0[j] - jp_prob_new[j])*(jp_prob_0[j] - jp_prob_new[j]);
		
		//testout
		//RPrint(" in Cal_W ----------- i_loop: "); RPrint(i_loop);
		//RPrint(" dif: "); RPrint(dif);
		//RPrint(" jp_prob_0: "); RPrint(jp_prob_0);
		//RPrint(" jp_prob_new: "); RPrint(jp_prob_new);
		
		if(dif < 1e-6) 
		{
			//TestOut<<" Cell_Prob... finished after iterations : "<< i_loop+1<<endl;
			break; 
		}
	
		//------------
		//check max iterations
		//------------
		if(i_loop == n_maximum_iteration-1)
		{
			TestOut<<"CAUTION!! max iteration reached in Cell_Prob..()"<<endl;
			//RPrint("CAUTION!! max iteration reached in Cell_Prob..()");
		}
		
		//------------
		//local deallocation
		//------------
		delete[] w12; 
	}
	
	
	//---------------------------
	//prep return 
	//the latest joint probability
	//---------------------------
	jp_prob_return.clear(); 
	jp_name_return.clear();
	for(int j=0; j<i_size_jp_prob; j++) 
	{
		jp_prob_return.push_back(jp_prob_new[j]); //return with jp_prob_new
		jp_name_return.push_back(jp_name_new[j]); //return with jp_name_new		
	}

	//testout
	//RPrint(" ========= Cell_Prob_Extension.. has successfully finished!", TestOut);

	//-------------
	//deallocation
	//-------------
	//delete[] w;
	delete[] w_ml; 
	delete[] w1;
	//delete[] id; 
	delete[] d_fmat1; 
	Del_dMatrix(d_ox, i_size_ol, ncol);
	Del_dMatrix(d_mx, i_size_ml, ncol);
	Del_dMatrix(mox, nrow, ncol);
	Del_dMatrix(uox, nrow, ncol);	
	
	delete[] cn; 
	delete[] cn0;
	delete[] s_ocn; 
	delete[] s_mcn; 
	delete[] s_fcd;
	
	return;
	
}