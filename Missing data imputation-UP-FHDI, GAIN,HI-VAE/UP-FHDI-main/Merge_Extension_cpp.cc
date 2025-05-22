//#include "ran_FHDI.h" //for uniform distribution 

//namespace FHDI{
#include <cstdlib> //for rand() and srand()
#include <time.h>  //for time()

void Merge_Extension_cpp(const int i_reci, double** uox, const int nrow_uox, 
					     double** mox, const int nrow_mox, int* tnU, 
						 std::string cn[],  int* ol, const int nrow_ol,
						 double** z, const int nrow, const int ncol, 
						 const int i_merge, 
						 const bool b_DEBUG)
//Description=========================================
// Merge the categorized matrix z 
//
// Algorithm: 											
//   For a given missing row at i_reci                 e.g., {12, NA, 4}           
//   Step 1: find columns having the observed cells:   e.g., 1st and 3rd columns   
//   Step 2: over all rows, search other observed at the 1st and 3rd columns 
//   Step 3: calculate relative distance, sum(|a-b|^2) 
//   Step 4: among rows that have the shortest distance, randomly select donors 
//   Step 5: fill the missing cell with the selected donors 
//
//
// original R code: Dr. Im, J. and Dr. Kim, J. 
// c++ code: 		Dr. Cho, I. 
// All rights reserved
// 
// updated: March 28, 2017
//----------------------------------------------------
//IN    : int i_reci = location of row with missing cell that has the least number of donors
//
//IN    : double uox(nrow_uox, ncol)= sorted unique patterns of observed cells. up to i_count_uox rows are meaningful 
//IN    : double mox(nrow_mox, ncol)= sorted unique patterns of missing  cells. up to i_count_mox rows are meaningful                           
//IN    : int tnU[nrow_uox]		= table format of the number of donors
//IN 	: string cn(nrow)		= vector of string to represent each row of z          
//IN	: int ol(nrow_ol)		= actual location of rows containing ONLY observed cells 
//IN    : int i_merge = random donor selection in Merge algorithm in Cell Make
//                            0= no random seed number setting
//						      1= random seed number setting 
//INOUT : double z(nrow, ncol)  = updated category matrix corresponding to original matrix x   
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
    if(i_merge == 0) std::srand(123);	    //turn on the fixed seed //This is still platform-dependent 
	                                    //maybe, use Numerical Recipe for platform-independent  
	
	//---------------
	//make a table of cn[full observed cells only]
	//---------------
	//std::string cn_ol[nrow_ol];
	std::string *cn_ol = new std::string[nrow_ol];
	
	for(int i=0; i<nrow_ol; i++) cn_ol[i] = cn[ol[i] - 1]; //Note: -1 for actual loc
	std::vector<std::string> v_table_cn_ol_row1; //names of the table
	std::vector<int> 		 v_table_cn_ol_row2; //counts of the table
	table_cpp(cn_ol, nrow_ol, v_table_cn_ol_row1, v_table_cn_ol_row2);
	const int i_size_v_table_cn_ol_row2 = (int)v_table_cn_ol_row2.size();	
	
	//----------------
	//get a string of the current row having missing cells, which has the least observed donors
	//----------------
	double* d_cn0 = new double[ncol];
	for(int i=0; i<ncol; i++) d_cn0[i] = mox[i_reci][i];
	std::string cn0; 
	Trans1(d_cn0, ncol, cn0);
	
	
	//----------------
	//ACTUAL locations of other missing rows that have the same string as cn0
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
	
	
	
	//-----------------
	//Which columns are NOT missing in mox[i_reci][]
	//-----------------
	double* d_mox_row = new double[ncol]; //temporary array
	for(int i=0; i<ncol; i++) d_mox_row[i] = mox[i_reci][i];
	std::vector<int> v_mxl; //ACTUAL location of non-missing column of mox[i_reci][]
	whichINV(d_mox_row, ncol, 0.0, v_mxl);
	const int i_nxl = (int)v_mxl.size(); //number of non-missing cell on this row
	delete[] d_mox_row; 
	
	//testout
	if(b_DEBUG){
	RPrint("v_mxl: "); RPrint(v_mxl);
	RPrint("i_nxl: "); RPrint(i_nxl);
	}
	
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
		for(int i=0; i<nrow_uox; i++) 
		{d_cand[i][0] = uox[i][v_mxl[0]-1]; } //-1 for ACTUAL location 
		
		//calculate distance using |a-b|^2
		const double d_mox_mxl = mox[i_reci][v_mxl[0]-1];
		distance2(d_cand, nrow_uox, i_nxl, d_mox_mxl, 
                  d_fdist);	 
		//---
		//NOTE: only the 1st value contains meaningful distance	
		//to avoid error in finding the minimum distance, 
		//in below, minimum searching needs due consideration
		//---
		
		//testout
		if(b_DEBUG){
		RPrint("successful so far2: i_nxl == 1 ");
		RPrint("d_fdist[]: "); RPrint(d_fdist, nrow_uox);
		RPrint("d_cand[][]: "); RPrint(d_cand, nrow_uox, ncol);
		}
	}
	
	
	if(i_nxl >1 ) //when current missing row has more than one column that has observed cells 
	{
		//------------
		//make a copy of all rows of all columns that correspond to the observed cells 
		//------------
		for(int i=0; i<nrow_uox; i++) 
		{
			for(int j=0; j<i_nxl; j++) //note: i_nxl is the length of v_mxl
			{ d_cand[i][j] = uox[i][v_mxl[j]-1]; } //-1 for ACTUAL location 
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
				double d_mox_temp = mox[i_reci][v_mxl[j]-1];
				double d_temp1 = d_cand[i][j]; 
				d_sum_dist +=  (d_mox_temp - d_temp1)*(d_mox_temp - d_temp1);
			}
			d_fdist[i] = d_sum_dist; 
		}
		//testout
		if(b_DEBUG){
		RPrint("successful so far3: : i_nxl > 1 ");
		RPrint("d_fdist[]: "); RPrint(d_fdist, nrow_uox);
		RPrint("d_cand[][]: "); RPrint(d_cand, nrow_uox, ncol);
		}
	}
	
	
	//------------
	//find the minimum distance
	//------------
	std::vector<int> v_floc; //ACTUAL location of minimum dist. in fdist
	double d_min_fdist = 0.0;
	if(i_nxl>=1) 
	{
		d_min_fdist = min_FHDI(d_fdist, nrow_uox);
		which(d_fdist, nrow_uox, d_min_fdist, v_floc); 
	}
	const int i_size_floc = (int)v_floc.size();
	if(i_size_floc <=0) { RPrint("Error! floc size is 0!"); return;}

    //testout
	if(b_DEBUG){
	RPrint("under the condition of i_nxl: "); RPrint(i_nxl);
	RPrint("d_min_fdist: "); RPrint(d_min_fdist);
	RPrint("v_floc: "); RPrint(v_floc);
	}
	
	//------------
	//select out a table of the location information of the minimum distance cells
	//------------
	int* i_nf = new int[i_size_floc]; 
    for(int i=0; i<i_size_floc; i++) i_nf[i] = v_table_cn_ol_row2[v_floc[i]-1]; //-1 for actual loc
	const int max_nf = max_FHDI(i_nf, i_size_floc); 

	//testout
	if(b_DEBUG){
	RPrint("max_nf: "); RPrint(max_nf);
	RPrint("i_nf[]: "); RPrint(i_nf, i_size_floc);
	}
	
	//-----------------------------
	//Case 1: more than 2 rows that have the smallest distance
	//-----------------------------
	if(max_nf>=2)
	{
		//-------------
		//find rows that have max nf
		//-------------
		std::vector<int> v_nf_max; 
		which(i_nf, i_size_floc, max_nf, v_nf_max); //Actual locations which have max of nf
		const int i_size_nf_max = (int)v_nf_max.size(); 
		
		//-------------
		//locations having the minimum distance between missing and observed cells
		//-------------
		std::vector<int> v_xloc;
		for(int i=0; i<i_size_nf_max; i++) v_xloc.push_back(v_floc[v_nf_max[i]-1]); //-1 for actual loc
		const int i_size_xloc = (int)v_xloc.size(); 
		
		//-------------
		//random number within [1, i_size_xloc]
		//Note: this is ACTUAL location
		//-------------
		int i_loc_rand_temp0 = 1; 
		if(i_merge == 1) i_loc_rand_temp0 = std::rand()%i_size_xloc + 1; //purely random 

		if(b_random) i_loc_rand_temp0 = 1 ; //for debugging // 
		const int i_loc_rand_xloc = v_xloc[i_loc_rand_temp0-1]; //-1 for actual loc
		
		//-------------
		//update z 
		//with the randomly selected row of donor  
		//-------------
		for(int i=0; i<i_nml; i++)//row-wise copy. Note: i_nml is the size of v_mloc[]
		{
			for(int j=0; j<i_nxl; j++)//column-wise copy. Note: i_nxl is the size of v_mxl[]
			{
				z[v_mloc[i]-1][v_mxl[j]-1] = uox[i_loc_rand_xloc-1][v_mxl[j]-1]; //-1 for actual location
			}			
		}
		//testout
		//RPrint("after Case 1");
		//RPrint("i_loc_rand_temp0: "); RPrint(i_loc_rand_temp0);
		//RPrint("i_loc_rand_xloc: ");  RPrint(i_loc_rand_xloc);
		//RPrint("z[][]: "); RPrint(z, nrow, ncol);
		
	}
	//testout
	//RPrint("successful so far5 ");
	
	
	//-----------------------------
	//Case 2: less than 2 rows that have the smallest distance
	//-----------------------------
	if(max_nf<2)
	{
		//-------------
		//find rows that is nf=1
		//-------------
		std::vector<int> v_nf_1; 
		which(i_nf, i_size_floc, 1, v_nf_1); //Actual locations which have 1 in nf
		const int i_size_nf_1 = (int)v_nf_1.size(); 
		
		//-------------
		//reduce floc to have rows with only nf=1
		//-------------
		std::vector<int> v_floc_1; 
		for(int i=0; i<i_size_nf_1; i++) 
		{
			int i_temp=v_nf_1[i]-1 ; //-1 for actual location 
			v_floc_1.push_back(v_floc[i_temp]); 
		}
		const int i_size_v_floc_1 = (int)v_floc_1.size(); 
		//testout
		//RPrint("successful so far5-1 ");
		//RPrint("i_size_nf_1: "); RPrint(i_size_nf_1);
		//RPrint("i_size_v_floc_1: "); RPrint(i_size_v_floc_1);
		
		//-------------
		//random number within [1, i_size_v_floc_1]
		//Note: this is ACTUAL location
		//-------------
		int i_loc_rand_temp = 1; 
		if(i_merge == 1) i_loc_rand_temp = std::rand()%i_size_v_floc_1 + 1; //purely random 

		if(b_random) i_loc_rand_temp = 1 ; //for debugging 
		const int i_loc_rand_floc = v_floc_1[i_loc_rand_temp-1]; //-1 for actual location
		//testout
		//RPrint("i_loc_rand_temp: "); RPrint(i_loc_rand_temp);
		//RPrint("i_loc_rand_floc: "); RPrint(i_loc_rand_floc);
		
		//-------------
		//Make a donor row and update pcell without the donor cell 
		//with the randomly selected row number  
		//-------------
		double* dcell = new double[ncol];
		for(int i=0; i<ncol; i++) dcell[i] = uox[i_loc_rand_floc-1][i]; //-1 for actual location
		//NOTE: pcell has (nrow_uox -1) rows
		double** pcell = New_dMatrix(nrow_uox-1, ncol); //by excluding the dcell
		for(int i=0; i<nrow_uox; i++)//row-wise copy. 
		{
			//---------------
			//below two conditions it to skip the dcell row
			//---------------
			if(i < (i_loc_rand_floc-1) ) 
			{
				for(int j=0; j<ncol; j++)//column-wise copy. 
				{pcell[i][j] = uox[i][j];} 
			}
			if(i > (i_loc_rand_floc-1) ) //to skip the dcell row 
			{
				for(int j=0; j<ncol; j++)//column-wise copy. 
				{pcell[i-1][j] = uox[i][j];} //note the -1 (one row shift)  		
			}
		}
		//testout
		//RPrint("successful so far5-2 ");		
		//RPrint("Inside Case 2");
		//RPrint("dcell[]: "); RPrint(dcell, ncol);
		//RPrint("pcell[][]: "); RPrint(pcell, nrow_uox-1, ncol);

		//-----------------------
		//calculate relative distance between pcell and dcell
		//-----------------------
		const int nrp = nrow_uox-1;  //total number of rows of pcell
		double* d_sdis = new double[nrp]; 
		Fill_dVector(d_sdis, nrp, 0.0);
		
		//-------------
		//distance between pcell and dcell
		//(1) when pcell is an array
		//-------------
		double sdis = 0.0; 
		if(nrp == 1)
		{
			for(int i=0; i<ncol; i++) 
				sdis += (pcell[0][i] - dcell[i])*(pcell[0][i] - dcell[i]);		
		}
		d_sdis[0] = sdis; 
		//testout
		//RPrint("successful so far5-3 ");				
		
		//-------------
		//distance between pcell and dcell
		//(2) when pcell is indeed a Matrix
		//-------------
		if(nrp > 1)
		{
			//-------------
			//calculate distance = sum(|a-b|^2) per row 
			//-------------	
			double d_sum_sdis = 0.0; 
			for(int i=0; i<nrp; i++)
			{
				d_sum_sdis = 0.0; //re-initialize
				for(int j=0; j<ncol; j++)
				{
					d_sum_sdis +=  (pcell[i][j]-dcell[j])*(pcell[i][j]-dcell[j]);
				}
				d_sdis[i] = d_sum_sdis; //NEEDS to check!
			}			
		}
		//testout
		//RPrint("successful so far5-4 ");						
		//RPrint("nrp: "); RPrint(nrp);
		//RPrint("d_sdis[]: "); RPrint(d_sdis, nrp);
		
		//-----
		//find the minimum sdis[]
		//-----
		const double d_min_sdis = min_FHDI(d_sdis, nrp);
		std::vector<int> v_sloc; //ACTUAL locations of min sdis 
		which(d_sdis, nrp, d_min_sdis, v_sloc);
		const int i_size_v_sloc = (int)v_sloc.size(); 
		//testout
		//RPrint("successful so far5-5 ");						
		//RPrint("d_min_sdis: "); RPrint(d_min_sdis);
		//RPrint("v_sloc[]: "); RPrint(v_sloc);
		
		//-----
		//exclude floc row from the table of cn[ol]
		//-----
		int* i_socn = new int[i_size_v_table_cn_ol_row2-1]; //reduced vector of tocn
		for(int i =0; i<i_size_v_table_cn_ol_row2; i++)
		{
			if(i < (i_loc_rand_floc-1)) //except for the floc actual location
			{
				i_socn[i] = v_table_cn_ol_row2[i];
			}
			if(i > (i_loc_rand_floc-1)) //except for the floc actual location
			{
				i_socn[i-1] = v_table_cn_ol_row2[i]; //Note: one index shift 
			}
		}
		//testout
		//RPrint("successful so far5-6 ");								
		//RPrint("i_socn[]: "); RPrint(i_socn, i_size_v_table_cn_ol_row2-1 );
		
		//-----
		//exclude floc row from tnU[]
		//------
		int* i_snU = new int[nrow_uox -1];//reduced array of tnU[]
		for(int i=0; i<nrow_uox; i++)
		{
			if(i < (i_loc_rand_floc-1)) //except for the floc actual location
			{
				i_snU[i] = tnU[i];
			}
			if(i > (i_loc_rand_floc-1)) //except for the floc actual location
			{
				i_snU[i-1] = tnU[i]; //Note: one index shift 
			}
		}
		//testout
		//RPrint("successful so far5-7 ");								
		//RPrint("i_snU[]: "); RPrint(i_snU, nrow_uox -1);
		
		//------
		//select out sloc rows from socn
		//------
		int* i_ns = new int[i_size_v_sloc]; //part of socn at sloc
		for(int i=0; i<i_size_v_sloc; i++) i_ns[i] = i_socn[(int)v_sloc[i]-1]; //-1 for actual location 
		const int max_i_ns = max_FHDI(i_ns, i_size_v_sloc);
		std::vector<int> v_i_ns; //Actual locations of max ns
		which(i_ns, i_size_v_sloc, max_i_ns, v_i_ns);
		const int i_size_v_i_ns = (int)v_i_ns.size();
		int* i_xloc = new int[i_size_v_i_ns];
		for(int i=0; i<i_size_v_i_ns; i++)
		{
			i_xloc[i] = v_sloc[(int)v_i_ns[i]-1]; //-1 for actual location
		}
		
		//------
		//get a random integer between [1, length(x_loc)]
		//------
		int i_loc_rand_temp2 = 1; 
		if(i_merge == 1) i_loc_rand_temp2 = std::rand()%i_size_v_i_ns + 1; //purely random 

		if(b_random) i_loc_rand_temp2 = 1 ; //for debugging
		const int i_loc_rand_xloc  = i_xloc[i_loc_rand_temp2-1]; //-1 for actual location

		//------
		//select out a row at floc from tnU and at xloc from snU
		//------
		const int i_crip1 = tnU[i_loc_rand_floc-1]; //-1 for actual location
		const int i_crip2 = i_snU[i_loc_rand_xloc-1]; //-1 for actual location	

		//testout
		//RPrint("successful so far6 ");
		//RPrint("i_loc_rand_temp2 : "); RPrint(i_loc_rand_temp2);
		//RPrint("i_loc_rand_xloc : "); RPrint(i_loc_rand_xloc);
		//RPrint("i_xloc[]: "); RPrint(i_xloc, i_size_v_i_ns);
		//RPrint("i_crip1: "); RPrint(i_crip1);
		//RPrint("i_crip2: "); RPrint(i_crip2);
		
		//========================
		//sub case 1: max of ns >= 2
		//========================
		if(max_i_ns >= 2)
		{
			//-----------
			//find rows that have the same string as the row at floc
			//-----------
			std::string s_ncn;
			double* d_uox_a_row = new double[ncol]; //temp array of a row of uox
			for(int i=0; i<ncol; i++) d_uox_a_row[i] = uox[i_loc_rand_floc-1][i]; //-1 for actual location
			Trans1(d_uox_a_row, ncol, s_ncn);
			std::vector<int> v_uloc; //actual locations where cn = ncn
			which(cn, nrow, s_ncn, v_uloc); 
			const int i_nul = (int)v_uloc.size(); //size of uloc
			
			//------------
			//replace z at row=mloc and column=mxl
			//------------
			for(int i=0; i<i_nml; i++)//row-wise copy. Note: i_nml is the size of v_mloc[]
			{
				for(int j=0; j<i_nxl; j++)//column-wise copy. Note: i_nxl is the size of v_mxl[]
				{
					z[v_mloc[i]-1][v_mxl[j]-1] 
					= pcell[i_loc_rand_xloc-1][v_mxl[j]-1]; //-1 for actual location
				}			
			}			
			//------------
			//replace z at row=uloc and all columns
			//------------
			for(int i=0; i<i_nul; i++)//row-wise copy. Note: i_nul is the size of v_uloc[]
			{
				for(int j=0; j<ncol; j++)//column-wise copy. Note: all columns
				{
					z[v_uloc[i]-1][j] 
					= pcell[i_loc_rand_xloc-1][j]; //-1 for actual location
				}			
			}			
			delete[] d_uox_a_row; 
			//testout
			//RPrint("===== after sub case 1============ ");
			//RPrint("v_uloc"); RPrint(v_uloc);
			//RPrint("i_nul"); RPrint(i_nul);
			//RPrint("updated z[][]: "); RPrint(z, nrow, ncol);

		}

		//========================
		//sub case 2: max of ns < 2 && crip1>=crip2
		//========================
		if(max_i_ns < 2 && i_crip1 >= i_crip2)
		{
			//-----------
			//find rows that have the same string as the row at floc
			//-----------
			std::string s_ncn2;
			double* d_uox_a_row2 = new double[ncol]; //temp array of a row of uox
			for(int i=0; i<ncol; i++) d_uox_a_row2[i] = pcell[i_loc_rand_xloc-1][i]; //-1 for actual location
			Trans1(d_uox_a_row2, ncol, s_ncn2);
			std::vector<int> v_uloc2; //actual locations where cn = ncn
			which(cn, nrow, s_ncn2, v_uloc2); 
			const int i_nul2 = (int)v_uloc2.size(); //size of uloc
			
			//------------
			//replace z at row=mloc and column=mxl
			//------------
			for(int i=0; i<i_nml; i++)//row-wise copy. Note: i_nml is the size of v_mloc[]
			{
				for(int j=0; j<i_nxl; j++)//column-wise copy. Note: i_nxl is the size of v_mxl[]
				{
					z[v_mloc[i]-1][v_mxl[j]-1] 
					= uox[i_loc_rand_floc-1][v_mxl[j]-1]; //-1 for actual location
				}			
			}			
			//------------
			//replace z at row=uloc and all columns
			//------------
			for(int i=0; i<i_nul2; i++)//row-wise copy. Note: i_nul2 is the size of v_uloc2[]
			{
				for(int j=0; j<ncol; j++)//column-wise copy. Note: all columns
				{
					z[v_uloc2[i]-1][j] 
					= uox[i_loc_rand_floc-1][j]; //-1 for actual location
				}			
			}			
			delete[] d_uox_a_row2; 
			//testout
			//RPrint("===== after sub case 2============ ");
			//RPrint("v_uloc2"); RPrint(v_uloc2);
			//RPrint("i_nul2"); RPrint(i_nul2);
			//RPrint("updated z[][]: "); RPrint(z, nrow, ncol);
			
		}		

		//========================
		//sub case 3: max of ns < 2 && crip1 < crip2
		//========================
		if(max_i_ns < 2 && i_crip1 < i_crip2)
		{
			//-----------
			//find rows that have the same string as the row at floc
			//-----------
			std::string s_ncn3;
			double* d_uox_a_row3 = new double[ncol]; //temp array of a row of uox
			for(int i=0; i<ncol; i++) d_uox_a_row3[i] = uox[i_loc_rand_floc-1][i]; //-1 for actual location
			Trans1(d_uox_a_row3, ncol, s_ncn3);
			std::vector<int> v_uloc3; //actual locations where cn = ncn
			which(cn, nrow, s_ncn3, v_uloc3); 
			const int i_nul3 = (int)v_uloc3.size(); //size of uloc
			
			//------------
			//replace z at row=mloc and column=mxl
			//------------
			for(int i=0; i<i_nml; i++)//row-wise copy. Note: i_nml is the size of v_mloc[]
			{
				for(int j=0; j<i_nxl; j++)//column-wise copy. Note: i_nxl is the size of v_mxl[]
				{
					z[v_mloc[i]-1][v_mxl[j]-1] 
					= pcell[i_loc_rand_xloc-1][v_mxl[j]-1]; //-1 for actual location
				}			
			}			
			//------------
			//replace z at row=uloc and all columns
			//------------
			for(int i=0; i<i_nul3; i++)//row-wise copy. Note: i_nul3 is the size of v_uloc3[]
			{
				for(int j=0; j<ncol; j++)//column-wise copy. Note: all columns
				{
					z[v_uloc3[i]-1][j] 
					= pcell[i_loc_rand_xloc-1][j]; //-1 for actual location
				}			
			}			
			delete[] d_uox_a_row3; 
			//testout
			//RPrint("===== after sub case 3============ ");
			//RPrint("v_uloc3"); RPrint(v_uloc3);
			//RPrint("i_nul3"); RPrint(i_nul3);
			//RPrint("updated z[][]: "); RPrint(z, nrow, ncol);
			
		}		
		//testout
		//RPrint("successful so far9 ");
		
		//-----------------------
		//local deallocation of memory used in Case 2
		//-----------------------
		delete[] dcell; 
		Del_dMatrix(pcell, nrow_uox-1, ncol);
		delete[] d_sdis; 
		delete[] i_socn; 
		delete[] i_snU; 
		delete[] i_xloc;
	} // end of case 2 nf< 2
	
	
	
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
//} //end of namespace