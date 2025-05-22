//#include "ran_FHDI.h" //for uniform distribution 
#include <cstdlib> //for rand() and srand()
#include <time.h>  //for time()

void KNN_Bigp(const int i_reci, double** uox, const int nrow_uox,
					     double** mox, const int nrow_mox, double* d_k, int** codes, const int i_option_collapsing,
						 std::string cn[],  int* ol, const int nrow_ol,
						 const int nrow, const int ncol, const int i_merge,
	                     std::vector<int> &v_nD, List_FHDI &List_nU, ofstream& TestOut)
//Description=========================================
// Find deficient donors for the recipient who has less than 2 donors by the Euclidean distance
//
// Algorithm: 											
//   For a given missing row at i_reci                 e.g., {12, NA, 4}           
//   Step 1: compute Euclidean distance from all unique observed patterns to the recipient  
//   Step 2: Case 0: if the recipient has a donor, then randomly select another one from list
//   Step 3: Case 1: if the recipient has no donor and the max occurnace of the candidate (who has minimum distance) in the list is >=2, then randomly select another one from list
//   Step 4: Case 2: if the recipient has no donor and the max occurnace of the candidate (who has minimum distance) in the list is <2, and the size of the list is >=2, 
//                   then randomly select two donors from list
//   Step 5: Case 3: if the recipient has no donor and the max occurnace of the candidate (who has minimum distance) in the list is <2, and the size of the list is 1, 
//                   then select another donor from the list of candicates who has the second minimum distance
//
//
// original R code: Dr. Im, J. and Dr. Kim, J. 
// c++ code: 		Yicheng Yang. 
// All rights reserved
// 
// updated: July 16, 2020
//----------------------------------------------------
//IN    : int i_reci = location of row with missing cell that has the least number of donors
//
//IN    : double uox(nrow_uox, ncol)= sorted unique patterns of observed cells. up to i_count_uox rows are meaningful 
//IN    : double mox(nrow_mox, ncol)= sorted unique patterns of missing  cells. up to i_count_mox rows are meaningful                           
//IN 	: string cn(nrow)		= vector of string to represent each row of z          
//IN	: int ol(nrow_ol)		= actual location of rows containing ONLY observed cells 
//INOUT : v_nD(nrow_mox)       = total number of donors of all mox 
//INOUT : List_nU(nrow_mox)    = List of actual location of donors of all mox in uox
//====================================================
{
	//below setting is for debugging random sampling
	const bool b_random = 0; //0: use rand(); 1:deterministic for debugging 
	
	//---------------
	//initialize random number generator
	//---------------
	if(i_merge == 1) std::srand(time(NULL)); //turn on random seed using C++ standard rand() fn 
	                                    //this will generate purely random numbers and 
										//will be platform-dependent, i.e., different on Win and Linux 
    if(i_merge == 0) std::srand(123);	//turn on the fixed seed //This is still platform-dependent 
	                                    //maybe, use Numerical Recipe for platform-independent  
	//TestOut<<"i_merge is "<< i_merge <<endl;
	int test_int = 0;
	if ((i_merge == 1)) test_int = std::rand() % 5 + 1; //purely random 
	//TestOut<<"test_int is "<< test_int << " at i_reci = "<< i_reci <<endl;

	//---------------
	//std::string cn_ol[nrow_ol];
	std::string *cn_ol = new std::string[nrow_ol];
	
	for(int i=0; i<nrow_ol; i++) cn_ol[i] = cn[ol[i] - 1]; //Note: -1 for actual loc
	std::vector<std::string> v_table_cn_ol_row1; //names of the table
	std::vector<int> 		 v_table_cn_ol_row2; //counts of the table
	table_cpp(cn_ol, nrow_ol, v_table_cn_ol_row1, v_table_cn_ol_row2);


	const int i_size_v_table_cn_ol_row2 = (int)v_table_cn_ol_row2.size();// same with nrow_uox
	
	//for (int p = 0; p < i_size_v_table_cn_ol_row2;p++) {
	//	TestOut<<"v_table_cn_ol_row1["<<p<<"]: "<< v_table_cn_ol_row1[p]<< setw(20) <<"v_table_cn_ol_row2["<<p<<"]: "<< v_table_cn_ol_row2[p]<<endl;
	//}
	//----------------
	//get a string of the current row having missing cells, which has the least observed donors
	//----------------
	double* d_cn0 = new double[ncol];
	for(int i=0; i<ncol; i++) d_cn0[i] = mox[i_reci][i];
	std::string cn0; 
	Trans1(d_cn0, ncol, cn0);
	
	
	//----------------
	//ACTUAL locations of other missing rows that have the same string as mox[i]
	//----------------
	std::vector<int> v_mloc;
    which(cn, nrow, cn0, v_mloc);	
	const int i_nml = (int)v_mloc.size(); 
	//testout
	/*
	if(b_DEBUG){
	RPrint("========in Merge============");
	RPrint("ol: "); RPrint(ol, nrow_ol);
	RPrint("cn_ol: "); RPrint(cn_ol, nrow_ol);
	RPrint("nrow_uox: "); RPrint(nrow_uox);
	RPrint("nrow_mox: "); RPrint(nrow_mox);
	RPrint("i_reci: "); RPrint(i_reci);
	RPrint("v_table_cn_ol_row2: "); RPrint(v_table_cn_ol_row2);
	RPrint("d_cn0[]: "); RPrint(d_cn0, ncol);
	RPrint("cn0[]: "); RPrint(&cn0,1 );
	RPrint("v_mloc: "); RPrint(v_mloc);
	RPrint("i_nml: "); RPrint(i_nml);
	}
	*/
	
	//TestOut << "List_nU inside KNN at i_recv = " << i_reci << endl;
	//List_nU.print_List_FHDI_yicheng(TestOut);
	
	//-----------------
	//Which columns are NOT missing in mox[i_reci][]
	//-----------------
	double* d_mox_row = new double[ncol]; //temporary array
	for(int i=0; i<ncol; i++) d_mox_row[i] = mox[i_reci][i];

	//TestOut<<"mox["<< i_reci <<"] is"<<endl;
	//for (int k1 = 0; k1 < ncol; k1++) {
	//	TestOut << setw(20) << d_mox_row[k1];
	//}
	//TestOut << endl;

	if (v_nD[i_reci] == 1) {

		double* d_oloc_temp11 = new double[1];

		List_nU.get_block(i_reci, d_oloc_temp11);

		int ll = d_oloc_temp11[0];

		//List_nU.print_one_List_FHDI_yicheng(i_reci, TestOut);
		//TestOut << "existing donor uox[" << ll << "] is " << endl;
		
		//for (int k1 = 0; k1 < ncol; k1++) {
		//	TestOut << setw(20) << uox[ll][k1];
		//}
		//TestOut << endl;

		delete[] d_oloc_temp11;
	}


	std::vector<int> v_mxl; //ACTUAL location of non-missing column of mox[i_reci][]
	//whichINV(d_mox_row, ncol, 0.0, v_mxl);

	for (int k = 0; k < i_option_collapsing; k++) {
		//if (codes[i_reci][i] == 0) { TestOut << "Error! Check the correlated_variables function, it is not correct !!!" << endl; }

		if (codes[i_reci][k] != 0) { // This is for case that the number of observed values in mox[i_recv] is smaller than i_option_collapsing
			v_mxl.push_back(codes[i_reci][k]);
		}
	}


	const int i_nxl = (int)v_mxl.size(); //number of non-missing cell on this row
	delete[] d_mox_row; 
	
	//-----------------
	//Find the nearest potential donor cells using "fdis"
	//NOTE: below two matrix and array has nrow_uox rows since it is 
	//related to observed cells uox
	//-----------------
	double ** d_cand = New_dMatrix(nrow_uox, ncol); //NOTE: the column may be flexible for below cases 
	double *  d_fdist= new double[nrow_uox];        //distance between entities 
	Fill_dVector(d_fdist, nrow_uox, 0.0);
	
	if(i_nxl == 1) //when the current missing row has only ONE observed cell   
	{
		//------------
		//make a copy of all rows of the one column 
		// that corresponds to the column where the observed cell of current missing row
        // is located 		
		//------------
		for (int i = 0; i < nrow_uox; i++)
		{
			d_cand[i][0] = (uox[i][v_mxl[0] - 1]) / (d_k[v_mxl[0] - 1]);
		} //-1 for ACTUAL location 
		
		//calculate distance using |a-b|^2
		const double d_mox_mxl = (mox[i_reci][v_mxl[0] - 1]) / (d_k[v_mxl[0] - 1]);
		distance2(d_cand, nrow_uox, i_nxl, d_mox_mxl, 
                  d_fdist);	 
		//---
		//NOTE: only the 1st value contains meaningful distance	
		//to avoid error in finding the minimum distance, 
		//in below, minimum searching needs due consideration
		//---
		
		//testout
		//if(b_DEBUG){
		//RPrint("successful so far2: i_nxl == 1 ");
		//RPrint("d_fdist[]: "); RPrint(d_fdist, nrow_uox);
		//RPrint("d_cand[][]: "); RPrint(d_cand, nrow_uox, ncol);
		//}
	}
	
	
	if(i_nxl >1 ) //when current missing row has more than one column that has observed cells 
	{
		//------------
		//make a copy of all rows of all columns that correspond to the observed cells 
		//------------
		for (int i = 0; i < nrow_uox; i++)
		{
			for (int j = 0; j < i_nxl; j++) //note: i_nxl is the length of v_mxl
			{
				d_cand[i][j] = (uox[i][v_mxl[j] - 1]) / (d_k[v_mxl[j] - 1]);
			} //-1 for ACTUAL location 
		}
		//-------------
		//calculate distance = sum(|a-b|^2) per row where mox[i][mxl] is the origin
		//-------------	
		double d_sum_dist = 0.0; 
		for(int i=0; i<nrow_uox; i++)
		{
			d_sum_dist = 0.0; //re-initialize
			for(int j=0; j<i_nxl; j++)
			{
				double d_mox_temp = (mox[i_reci][v_mxl[j] - 1]) / (d_k[v_mxl[j] - 1]);
				double d_temp1 = d_cand[i][j]; 
				d_sum_dist +=  (d_mox_temp - d_temp1)*(d_mox_temp - d_temp1);
			}
			d_fdist[i] = d_sum_dist; 
		}
		//testout
		//if(b_DEBUG){
		//RPrint("successful so far3: : i_nxl > 1 ");
		//RPrint("d_fdist[]: "); RPrint(d_fdist, nrow_uox);
		//RPrint("d_cand[][]: "); RPrint(d_cand, nrow_uox, ncol);
		//}
	}
	
	//-------------
	//set the distance from obatined donors in uox to mox[i_reci] as 1234567 instead of 0s
	//if the distnace is 0, that means this uox is already a donor in the list
	//-------------
	for (int k1 = 0; k1 < nrow_uox; k1++) {
		if (d_fdist[k1] == 0) {
			d_fdist[k1] = 1234567.0;
		}
	}

	//TestOut<<"d_fdist is"<<endl;
	//for (int k11 = 0; k11 < nrow_uox; k11++) {
	//	TestOut<<"d_fdist["<<k11<<"]: "<< d_fdist [k11]<<endl;
	//}
	//------------
	//find the minimum distance
	//------------
	std::vector<int> v_floc; //ACTUAL location of donors in uox who has the minimum distance
	double d_min_fdist = 0.0;
	if(i_nxl>=1) 
	{
		d_min_fdist = min_FHDI(d_fdist, nrow_uox);
		which(d_fdist, nrow_uox, d_min_fdist, v_floc); 
	}

	//TestOut<<"d_min_fdist is "<< d_min_fdist <<endl;
	const int i_size_floc = (int)v_floc.size();
	if(i_size_floc <=0) { RPrint("Error! floc size is 0!"); return;}

	//TestOut << "v_floc is" << endl;
	//for (int k15 = 0; k15 < i_size_floc; k15++) {
	//	TestOut << "v_floc[" << k15 << "]: " << v_floc[k15] << endl;
	//}
    //testout
	//if(b_DEBUG){
	//RPrint("under the condition of i_nxl: "); RPrint(i_nxl);
	//RPrint("d_min_fdist: "); RPrint(d_min_fdist);
	//RPrint("v_floc: "); RPrint(v_floc);
	//}
	
	//------------
	//select out a table of the location information of the minimum distance cells
	//------------
	int* i_nf = new int[i_size_floc]; // occurance of all donors in uox for mox[i]
    for(int i=0; i<i_size_floc; i++) i_nf[i] = v_table_cn_ol_row2[v_floc[i]-1]; //-1 for actual loc
	const int max_nf = max_FHDI(i_nf, i_size_floc); //highest occueance of all donors

	//testout
	//if(b_DEBUG){
	//RPrint("max_nf: "); RPrint(max_nf);
	//RPrint("i_nf[]: "); RPrint(i_nf, i_size_floc);
	//}

	//-------------
	//find rows that have max nf
    //-------------
	std::vector<int> v_nf_max;
	which(i_nf, i_size_floc, max_nf, v_nf_max); //Actual locations which have max of nf
	const int i_size_nf_max = (int)v_nf_max.size();

	//-------------
	//locations having the minimum distance between missing and observed cells
	//-------------
	std::vector<int> v_xloc;// actual locations of donors in uox who has minimum distance and highest occurance
	for (int i = 0; i < i_size_nf_max; i++) v_xloc.push_back(v_floc[v_nf_max[i] - 1]); //-1 for actual loc
	const int i_size_xloc = (int)v_xloc.size();

	//TestOut << "v_xloc is" << endl;
	//for (int k12 = 0; k12 < i_size_xloc; k12++) {
	//	TestOut << "v_xloc[" << k12 << "]: " << v_xloc[k12] << endl;
	//}

	// Example:
	// v_floc = {1,3,7,9} -> occur = {3,3,1,1}
	// v_xloc ={1,3} -> occur = {3,3}
	// i_size_xloc = i_size_nf_max =2

	//---------------------------------------
	//Case 0: if mox[i] has one donor and it only needs one more donor from uox who has the smallest distance

	//a) i_size_xloc >=2, need randomly select 1
	//{3,3,2,1} i_size_xloc=2
	//{1,1} i_size_xloc=2

	//b) i_size_xloc == 1, no random selection
	//{3,2,1,1} i_size_xloc=1
	//{2} ..
	//{1} .. 

	//----------------------------------------
	//TestOut<<"List_nU Before:"<<endl;
	//List_nU.print_one_List_FHDI_yicheng(i_reci,TestOut);

	if (v_nD[i_reci] == 1) {
		//TestOut<<"CASE0: KNN with v_nD["<< i_reci <<"] with 1 donor and max_nf ="<< max_nf <<endl;
		//-------------
		//random number within [1, i_size_xloc]
		//Note: this is ACTUAL location
		//-------------
		int i_loc_rand_temp0 = 1;
		if ((i_merge == 1)&&(i_size_xloc >=2)) i_loc_rand_temp0 = std::rand() % i_size_xloc + 1; //purely random 

		//if(b_random) i_loc_rand_temp0 = 1 ; //for debugging // 
		const int i_loc_rand_xloc = v_xloc[i_loc_rand_temp0 - 1]; //-1 for actual loc

		v_nD[i_reci] = 1 + max_nf;

		//TestOut << "i_loc_rand_xloc is " << i_loc_rand_xloc << endl;
		//double* d_oloc_temp1 = new double[1];

		//List_nU.get_block(i_reci, d_oloc_temp1);

		//TestOut<<"Case0 d_oloc_temp1 is "<< d_oloc_temp1 [0]<<endl;

		//double* d_oloc_temp = new double[2];

		//d_oloc_temp[0] = d_oloc_temp1[0]; // Note it should be acutual location

		//TestOut<<"!!!!Actual donor location in uox is "<< d_oloc_temp1[0] <<endl;
		//
		//d_oloc_temp[1] = i_loc_rand_xloc; // Note it should be acutual location

		//TestOut<<"d_oloc_temp is"<<endl;
		//for (int k1 = 0; k1 < 2; k1++) {
		//	TestOut << setw(20) << d_oloc_temp[k1];
		//}
		//TestOut << endl;

		std::vector<double> d_oloc_temp;

		d_oloc_temp.push_back(i_loc_rand_xloc);

		//sort(d_oloc_temp.begin(), d_oloc_temp.end());

		List_nU.put_block_yicheng(i_reci, 1, d_oloc_temp);

		//delete[] d_oloc_temp1;

		//delete[] d_oloc_temp;

	}



	if (v_nD[i_reci] == 0) {
	

    //-----------------------------
	//Case 1: if mox[i] has no donors and the highest occurance (i.e., max_nf) of donors with the smallest distance is >= 2
	//        then it needs only one more donor from uox who has the smallest distance

	//a) max_nf >= 2
	// {3,3,2,1} i_size_xloc >=2, need randomly select 1

    // {3,2,1,1} i_size_xloc ==1,no random selection
	// {2} ..

		if (max_nf >= 2) {
			//TestOut << "CASE1: KNN with v_nD[" << i_reci << "] with no donor and max_nf = "<< max_nf << endl;
			//-------------
			//random number within [1, i_size_xloc]
			//Note: this is ACTUAL location
			//-------------
			int i_loc_rand_temp0 = 1;
			if ((i_merge == 1) && (i_size_xloc >= 2)) i_loc_rand_temp0 = std::rand() % i_size_xloc + 1; //purely random 

			//if(b_random) i_loc_rand_temp0 = 1 ; //for debugging // 
			const int i_loc_rand_xloc = v_xloc[i_loc_rand_temp0 - 1]; //-1 for actual loc

			v_nD[i_reci] = max_nf;

			//double* d_oloc_temp = new double[1];

			//d_oloc_temp[0] = i_loc_rand_xloc; // Note it should be acutual location

			//TestOut << "i_loc_rand_xloc is " << i_loc_rand_xloc << endl;
			//for (int k1 = 0; k1 < 1; k1++) {
			//	TestOut << setw(20) << d_oloc_temp[k1];
			//}
			//TestOut << endl;
			std::vector<double> d_oloc_temp;
			d_oloc_temp.push_back(i_loc_rand_xloc);

			//sort(d_oloc_temp.begin(), d_oloc_temp.end());

			List_nU.put_block_yicheng(i_reci, 1, d_oloc_temp);


			//delete[] d_oloc_temp;

		}
	

    //Case 2: if mox[i] has no donors and the highest occurance (i.e., max_nf) of donors with the smallest distance is < 2 and 
    //        v_floc.size() is >= 2,
	//        then it needs two donors from uox who has the smallest distance

	// {1,1,1,1} need randomly select 2

	// {1,1} no random selection 

		if ((max_nf < 2) && (v_floc.size() >= 2)) {
			//TestOut << "CASE2: KNN with v_nD[" << i_reci << "] with no donor and max_nf = " << max_nf <<", and v_floc size is "<< v_floc.size() << endl;
			//-------------
			//random number within [1, i_size_xloc]
			//Note: this is ACTUAL location
			//-------------
			int i_loc_rand_temp0 = 1;
			int i_loc_rand_temp1 = 2;

			if ((i_merge == 1) && (i_size_xloc > 2)) {
				i_loc_rand_temp0 = 0;// Make sure it will go into the while loop
				i_loc_rand_temp1 = 0;
				while (i_loc_rand_temp0 == i_loc_rand_temp1) {
					i_loc_rand_temp0 = std::rand() % i_size_xloc + 1;
					i_loc_rand_temp1 = std::rand() % i_size_xloc + 1;
				}
			}

			//if ((i_merge == 1) && (i_size_xloc > 2)) i_loc_rand_temp0 = std::rand() % i_size_xloc + 1; //purely random 
			//if ((i_merge == 1) && (i_size_xloc > 2)) i_loc_rand_temp1 = std::rand() % i_size_xloc + 1; //purely random 

			//if (i_loc_rand_temp0 == i_loc_rand_temp1) { RPrint("Error in KNN of cell make for random selection!!!"); return; }

			if (i_loc_rand_temp0 == i_loc_rand_temp1) {
				TestOut<<"Error in KNN of cell make for random selection!!! i_loc_rand_temp0 = "<< i_loc_rand_temp0 <<"; i_loc_rand_temp1 = "<< i_loc_rand_temp1 <<endl;
				return;
			}
			//if(b_random) i_loc_rand_temp0 = 1 ; //for debugging // 
			const int i_loc_rand_xloc = v_xloc[i_loc_rand_temp0 - 1]; //-1 for actual loc
			const int i_loc_rand_xloc1 = v_xloc[i_loc_rand_temp1 - 1]; //-1 for actual loc

			v_nD[i_reci] = 2;

			//double* d_oloc_temp = new double[2];

			//d_oloc_temp[0] = i_loc_rand_xloc; // Note it should be acutual location
			//d_oloc_temp[1] = i_loc_rand_xloc1; // Note it should be acutual location

			//TestOut << "d_oloc_temp is" << endl;
			//for (int k1 = 0; k1 < 2; k1++) {
			//	TestOut << setw(20) << d_oloc_temp[k1];
			//}
			//TestOut << endl;
			std::vector<double> d_oloc_temp;
			d_oloc_temp.push_back(i_loc_rand_xloc);
			d_oloc_temp.push_back(i_loc_rand_xloc1);

			//sort(d_oloc_temp.begin(), d_oloc_temp.end());

			List_nU.put_block_yicheng(i_reci, 2, d_oloc_temp);


			//delete[] d_oloc_temp;
		}

	//Case 3: if mox[i] has no donors and the highest occurance (i.e., max_nf) of donors with the smallest distance is < 2 and 
	//        v_floc.size() is = 1,
	//        then it needs two donors from uox who has the smallest distance

	// {1} no random selection and find in the second minimum set
		if ((max_nf < 2) && (v_floc.size() == 1)) {
			//TestOut << "CASE3: KNN with v_nD[" << i_reci << "] with no donor and max_nf = " << max_nf << ", and v_floc size is " << v_floc.size() << endl;
			std::vector<int> v_floc_second; //ACTUAL location of donors in uox who has the second minimum distance
			double d_second_min_fdist = 0.0;
			if (i_nxl >= 1)
			{
				d_second_min_fdist = second_min_FHDI(d_fdist, nrow_uox);
				which(d_fdist, nrow_uox, d_second_min_fdist, v_floc_second);
			}

			if(d_second_min_fdist == d_min_fdist) { RPrint("Error in KNN of cell make for getting the second minimum distance!!!"); return; }
			const int i_size_floc_second = (int)v_floc_second.size();
         
			int* i_nf_second = new int[i_size_floc_second]; // occurance of all donors in uox for mox[i]
			for (int i = 0; i < i_size_floc_second; i++) i_nf_second[i] = v_table_cn_ol_row2[v_floc_second[i] - 1]; //-1 for actual loc
			const int max_nf_second = max_FHDI(i_nf_second, i_size_floc_second); //highest occueance of all donors

			//-------------
	        //find rows that have max_nf_second
	        //-------------
			std::vector<int> v_nf_max_second;
			which(i_nf_second, i_size_floc_second, max_nf_second, v_nf_max_second); //Actual locations which have max of nf
			const int i_size_nf_max_second = (int)v_nf_max_second.size();

			//-------------
            //locations having the second minimum distance between missing and observed cells
            //-------------
			std::vector<int> v_xloc_second;// actual locations of donors in uox who has the seond minimum distance and highest occurance
			for (int i = 0; i < i_size_nf_max_second; i++) v_xloc_second.push_back(v_floc_second[v_nf_max_second[i] - 1]); //-1 for actual loc
			const int i_size_xloc_second = (int)v_xloc_second.size();

			//for (int m = 0;m < i_size_xloc_second;m++) TestOut<<"v_xloc_second["<<m<<"]: "<< v_xloc_second [m]<< endl;

			v_nD[i_reci] = 1 + max_nf_second;

			int i_loc_rand_temp0 = 1;
			if ((i_merge == 1) && (i_size_xloc_second >= 2)) i_loc_rand_temp0 = std::rand() % i_size_xloc + 1; //purely random 
			const int i_loc_rand_xloc = v_xloc_second[i_loc_rand_temp0 - 1]; //-1 for actual loc

			//double* d_oloc_temp = new double[2];


			//d_oloc_temp[0] = v_floc[0];
			//d_oloc_temp[1] = i_loc_rand_xloc; // Note it should be acutual location

			//TestOut << "d_oloc_temp is" << endl;
			//for (int k1 = 0; k1 < 2; k1++) {
			//	TestOut << setw(20) << d_oloc_temp[k1];
			//}
			//TestOut << endl;

			std::vector<double> d_oloc_temp;
			d_oloc_temp.push_back(v_floc[0]);
			d_oloc_temp.push_back(i_loc_rand_xloc);

			//sort(d_oloc_temp.begin(), d_oloc_temp.end());

			List_nU.put_block_yicheng(i_reci, 2, d_oloc_temp);


			//delete[] d_oloc_temp;

			delete[] i_nf_second;
		}
		


	}
	
	//TestOut << "List_nU After:" << endl;
	//List_nU.print_one_List_FHDI_yicheng(i_reci, TestOut);
	//List_nU.print_List_FHDI_yicheng(TestOut);
	//------------------
	//Deallocation
	//------------------
	delete[] d_cn0; 
	Del_dMatrix(d_cand, nrow_uox, ncol);
    delete[] d_fdist;	
	delete[] i_nf; 
	
	delete[] cn_ol;
	
	return; //temporary ending 
}
