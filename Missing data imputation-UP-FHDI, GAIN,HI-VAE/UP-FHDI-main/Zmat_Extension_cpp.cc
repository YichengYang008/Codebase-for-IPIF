//namespace FHDI
//{
//local fn declaration
//void Trans(double** z, const int nrow, const int ncol, std::string cn[]); //many rows case
//void Trans1(double* z, const int n, std::string &cn); //one row case
#include "tune_points_string.cc"
void Zmat_Extension_cpp(double** z, const int nrow, const int ncol, std::string cn[],
	int* ml, int* ol, int& i_count_ol, int& i_count_ml,
	double** uox, double** mox, int &i_count_uox, int &i_count_mox,
	const bool b_DEBUG, ofstream& TestOut)
	//Description=========================================
	// make the condensed expression of z
	//
	// Algorithm:  each row of z will be concatenated as a single string consisting of 35 characters
	// 
	// Note: as of Oct 2016, NA values (missing data) is marked by a long integer at the parent "r" code
	//
	// original R code: Dr. Im, J. and Dr. Kim, J. 
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: March 28, 2017
	//----------------------------------------------------
	//IN   	: double z(nrow, ncol)  = categorized matrix corresponding to original matrix x
	//OUT	: string cn(nrow)		= vector of string to represent each row of the caterorized z          
	//OUT	: int ml(nrow)			= actual location of rows containing AT LEAST ONE missing cells
	//OUT	: int ol(nrow)			= actual location of rows containing ONLY observed cells  
	//OUT   : int i_count_ol		= total number of ol rows that containing observed cells
	//OUT   : int i_count_ml		= total number of ml rows that containing missing cells 
	//OUT   : double uox(nrow, ncol)= sorted unique categorized patterns of observed cells. up to i_count_uox rows are meaningful 
	//OUT   : double mox(nrow, ncol)= sorted unique categorized patterns of missing  cells. up to i_count_mox rows are meaningful                           
	//OUT   : int i_count_uox		= total number of uox rows that containing meaningful cells
	//OUT   : int i_count_mox		= total number of mox rows that containing meaningful cells
	//====================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;
	//--------------
	//make a condensed expression "cn" of z
	//--------------
	double zmat_check1_begin = MPI_Wtime();
	Trans(z, nrow, ncol, cn);

	//testout
	if (b_DEBUG)
	{
		RPrint("=====in Zmat... ========");
		RPrint("cn"); RPrint(cn, nrow);
		RPrint("z"); RPrint(z, nrow, ncol);
	}

	//--------------
	//locations of missing cells (ml) and observed cells (ol)
	//--------------
	Fill_iVector(ml, nrow, 0); Fill_iVector(ol, nrow, 0); //Initialization
	double d_temp = 0.0;
	int i_ol_temp = 0;
	int i_ml_temp = 0;
	for (int i_row = 0; i_row<nrow; i_row++)
	{
		d_temp = 1.0;
		for (int i_col = 0; i_col<ncol; i_col++)
		{
			if (z[i_row][i_col] == 0) { d_temp = 0.0; break; } //found zero, i.e. missing cell
		}

		if (fabs(d_temp) > 1e-15) //this row has no missing cells
		{
			ol[i_ol_temp] = i_row + 1; i_ol_temp++;
		} //actual number of the row having no missing cells

		if (fabs(d_temp) < 1e-15) //this row has AT LEAST one missing cells
		{
			ml[i_ml_temp] = i_row + 1; i_ml_temp++;
		}  //actual number of the row having missing cells
	}

	if (i_ol_temp == 0) { TestOut << "Error! no observed unit in Zmat_Extension" << endl; return; }
	//if(i_ml_temp ==0) {cout<<"Error! no missing unit"<<endl; return; }
	//cout << "zmat_check1_begin at node " << mynode << " = " << MPI_Wtime() - zmat_check1_begin << endl;

	double zmat_check2_begin = MPI_Wtime();
	i_count_ol = i_ol_temp; //update the actual value 
	i_count_ml = i_ml_temp; //update the actual value

							//---------------------
							//make UNIQUE patterns of z by cn
							//i.e., uox and mox
							//---------------------
							//step . Sort the "cn" in the ascending order 
							//---------------------
							//std::string s_ol[i_ol_temp]; //string vector of observed patterns only
							//std::string s_ml[i_ml_temp]; //string vector of missing patterns only
	std::string *s_ol = new std::string[i_ol_temp]; //string vector of observed patterns only
	std::string *s_ml = new std::string[i_ml_temp]; //string vector of missing patterns only	
	for (int i = 0; i<i_ol_temp; i++) { s_ol[i] = cn[ol[i] - 1]; } //"-1" since ol contains actual row number
	for (int i = 0; i<i_ml_temp; i++) { s_ml[i] = cn[ml[i] - 1]; } //"-1" since ml contains actual row number	

	std::sort(s_ol, s_ol + i_ol_temp); //knowing that s_ol[] has i_ol_temp entities
	std::sort(s_ml, s_ml + i_ml_temp); //knowing that s_ml[] has i_ml_temp entities
	//cout << "zmat_check2_begin at node " << mynode << " = " << MPI_Wtime() - zmat_check2_begin << endl;
	//------------
	//memorize observed patterns 
	//------------
	//--------------------------------------
	//cyclic distribute
	int number_before_split = 0;
	int number_after_split = 0;
	int startpoint = 0;
	int endpoint = 0;
	int splitnode = 0;
	double** uox_temp = New_dMatrix(i_ol_temp, ncol);
	double** uox_recv = New_dMatrix(i_ol_temp, ncol);
	//cout<<"nrow_mox at mynode "<<mynode<<" is "<< nrow_mox <<endl;
	//cout << "nrow_uox at mynode " << mynode << " is " << nrow_uox << endl;
	//cout << "i_ol_temp: " << i_ol_temp << endl;
	//cout << "i_ml_temp: " << i_ml_temp << endl;

	if (i_ol_temp < (totalnodes - 1)) {
		TestOut<<"Error! The number of observed instances is smaller than the total number of available nodes"<<endl;
	}

	if(i_ol_temp > totalnodes) {
		if (i_ol_temp % (totalnodes - 1) != 0) { splitnode = i_ol_temp % (totalnodes - 1); }
		if (i_ol_temp % (totalnodes - 1) == 0) { splitnode = 1; }
		//if (mynode == 0) cout << "Split: " << splitnode << endl;
		number_after_split = floor(1.0*i_ol_temp / (1.0*totalnodes - 1));
		number_before_split = 1.0*(i_ol_temp - floor(1.0*i_ol_temp / (1.0*totalnodes - 1)) *(totalnodes - splitnode - 1)) / splitnode;

		if (mynode >= 1 && mynode <= splitnode) {
			startpoint = (mynode - 1) * number_before_split;
			endpoint = mynode *number_before_split;
		}
		if (mynode > splitnode) {
			startpoint = splitnode * number_before_split + (mynode - splitnode - 1)*number_after_split;
			endpoint = splitnode * number_before_split + (mynode - splitnode)*number_after_split;
		}
		if ((number_before_split*splitnode + number_after_split* (totalnodes - splitnode - 1)) != i_ol_temp) {
			TestOut << "Error! Work Assignment is wrong in Zmat_Extension!!!!" << endl;
			return;
		}

		//for (int o = 0;o < i_ol_temp;o++) {
		//	cout << "s_ol[" << o << "]: " << s_ol[o] << " ---at node " << mynode << endl;
		//}
		//cout << "Before Zmat starting: " << startpoint << " and ending: " << endpoint << " on node " << mynode << endl;
		tune_points_string(startpoint, endpoint, s_ol, mynode, totalnodes, i_ol_temp);
		//cout << "After Zmat starting: " << startpoint << " and ending: " << endpoint << " on node " << mynode << endl;

		//if (mynode == 0) cout << "ZMAT 1" << endl;
		//-------------------------------------
		double zmat_check3_begin = MPI_Wtime();
		int i_count_uox_temp = 0; //total number of unique uox 
		std::string s_temp;
		for (int i = startpoint; i<endpoint; i++)
		{
			//if (mynode == 1) {
			//	cout << "new zmat " << endl;
			//}
			s_temp = s_ol[i]; //get a string 

			if (i == 0 || (i > 0 && s_temp.compare(s_ol[i - 1]) != 0)) {
				//if (mynode == 1) {
				//	cout << "new zmat "<< endl;
				//}
				for (int j = 0;j < nrow;j++) {
					if (s_temp.compare(cn[j]) == 0) {
						for (int k = 0; k < ncol; k++)
						{
							uox_temp[i_count_uox_temp][k] = z[j][k];
						} //store the found observed pattern
						i_count_uox_temp++;
						break;
					}
				}
			}

		}
		//RPrint("uox_temp at node ");RPrint(mynode);
		//RPrint(uox_temp, i_ol_temp, ncol);
		//cout<<"zmat debug1 on node "<<mynode<<endl;
		//cout<<"i_ol_temp: "<< i_ol_temp <<" at node "<<mynode<<endl;
		if (mynode != 0) {
			MPI_Send(uox_temp[0], (i_ol_temp*ncol), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		}
		//cout << "zmat debug2 on node " << mynode << endl;
		if (mynode == 0) {
			i_count_uox = 0;
			for (int j = 1; j < totalnodes; j = j + 1) {
				MPI_Recv(uox_recv[0], (i_ol_temp*ncol), MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
				//cout<<"I'm here zmat at node "<<j<<endl;
				//RPrint("uox_recv: ");
				//RPrint(uox_recv, i_ol_temp, ncol);

				int length = 0;
				for (int l = 0; l < i_ol_temp; l++) {
					//cout << "l: " << l << endl;
					int sum_temp = 0;
					for (int m = 0; m < ncol; m++) {
						sum_temp = sum_temp + uox_recv[l][m];
						if (sum_temp > 1e-3) {
							length = l + 1;
							break;
						}
					}
				}

				//cout<<"length at node "<<j<<" is "<<length<<endl;

				int counter1 = 0;
				int flag = i_count_uox;
				int flag_end = i_count_uox + length;
				//cout<<"i_count_uox + length = "<< i_count_uox + length <<endl;

				for (int k = flag; k < flag_end;k++) {
					for (int m = 0;m < ncol;m++) {
						uox[k][m] = uox_recv[counter1][m];
					}
					//cout <<"counter1: "<< counter1 <<endl;
					//cout << "i_count_uox: " << i_count_uox << endl;
					counter1 = counter1 + 1;
					i_count_uox++;
				}
				//cout << "counter1 and i_count_uox are " << counter1 << " " << i_count_uox << endl;

			}
		}
		//cout << "zmat debug3 on node " << mynode << endl;

		MPI_Bcast(&i_count_uox, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(uox[0], (nrow*ncol), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}// end of parallel ol version
	else {
		i_count_uox = 0; //total number of unique uox 
		std::string s_temp;
		for (int i = 0; i<i_ol_temp; i++)
		{
			s_temp = s_ol[i]; //get a string 

			if (i == 0 || (i > 0 && s_temp.compare(s_ol[i - 1]) != 0)) {
				//if (mynode == 0) {
				//	cout << "III: " << i << endl;
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
	}// end of serial ol version
	
	//cout<<"final i_count_uox: "<< i_count_uox <<endl;
	//RPrint("final uox: ");RPrint("at node ");RPrint(mynode);
	//RPrint(uox, nrow, ncol);

	//cout << "zmat_check3_begin at node " << mynode << " = " << MPI_Wtime() - zmat_check3_begin << endl;
	//Now, i_count_uox means the total number of unique observed patterns
	//if (mynode == 0) cout<<"ZMAT 2"<<endl;

	//------------
	//memorize missing patterns 
	//------------
	 number_before_split = 0;
	 number_after_split = 0;
	 startpoint = 0;
	 endpoint = 0;
	 splitnode = 0;
	 double** mox_temp = New_dMatrix(i_ml_temp, ncol);
	 double** mox_recv = New_dMatrix(i_ml_temp, ncol);
	//cout<<"nrow_mox at mynode "<<mynode<<" is "<< nrow_mox <<endl;
	//cout << "nrow_uox at mynode " << mynode << " is " << nrow_uox << endl;
	 if (i_ml_temp > totalnodes) {
		 if (i_ml_temp % (totalnodes - 1) != 0) { splitnode = i_ml_temp % (totalnodes - 1); }
		 if (i_ml_temp % (totalnodes - 1) == 0) { splitnode = 1; }
		// if (mynode == 0) cout << "Split: " << splitnode << endl;
		 number_after_split = floor(1.0*i_ml_temp / (1.0*totalnodes - 1));
		 number_before_split = 1.0*(i_ml_temp - floor(1.0*i_ml_temp / (1.0*totalnodes - 1)) *(totalnodes - splitnode - 1)) / splitnode;

		 if (mynode >= 1 && mynode <= splitnode) {
			 startpoint = (mynode - 1) * number_before_split;
			 endpoint = mynode *number_before_split;
		 }
		 if (mynode > splitnode) {
			 startpoint = splitnode * number_before_split + (mynode - splitnode - 1)*number_after_split;
			 endpoint = splitnode * number_before_split + (mynode - splitnode)*number_after_split;
		 }
		 if ((number_before_split*splitnode + number_after_split* (totalnodes - splitnode - 1)) != i_ml_temp) {
			 TestOut << "Work Assignment Error!!!!" << endl;
			 return;
		 }
		 //cout << "Before Zmat starting: " << startpoint << " and ending: " << endpoint << " on node " << mynode << endl;
		 tune_points_string(startpoint, endpoint, s_ml, mynode, totalnodes, i_ml_temp);
		 //cout << "After Zmat starting: " << startpoint << " and ending: " << endpoint << " on node " << mynode << endl;



		 double zmat_check4_begin = MPI_Wtime();
		 int i_count_mox_temp = 0; //total number of unique mox 
		 std::string s_temp;
		 for (int i = startpoint; i < endpoint; i++)
		 {
			 s_temp = s_ml[i]; //get a string 

			 if (i == 0 || (i > 0 && s_temp.compare(s_ml[i - 1]) != 0)) {
				 //if (mynode == 0) {
					// cout << "III: " << i << endl;
				 //}
				 for (int j = 0;j < nrow;j++) {
					 if (s_temp.compare(cn[j]) == 0) {
						 for (int k = 0; k < ncol; k++)
						 {
							 mox_temp[i_count_mox_temp][k] = z[j][k];
						 } //store the found observed pattern
						 i_count_mox_temp++;
						 break;
					 }
				 }
			 }

		 }

		 if (mynode != 0) {
			 MPI_Send(mox_temp[0], (i_ml_temp*ncol), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		 }
		 //cout << "zmat debug2 on node " << mynode << endl;
		 if (mynode == 0) {
			 i_count_mox = 0;
			 for (int j = 1; j < totalnodes; j = j + 1) {
				 MPI_Recv(mox_recv[0], (i_ml_temp*ncol), MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
				 //cout << "I'm here zmat at node " << j << endl;
				 //RPrint("uox_recv: ");
				 //RPrint(uox_recv, i_ol_temp, ncol);

				 int length = 0;
				 for (int l = 0; l < i_ml_temp; l++) {
					 //cout << "l: " << l << endl;
					 int sum_temp = 0;
					 for (int m = 0; m < ncol; m++) {
						 sum_temp = sum_temp + mox_recv[l][m];
						 if (sum_temp > 1e-3) {
							 length = l + 1;
							 break;
						 }
					 }
				 }

				 //cout << "length at node " << j << " is " << length << endl;

				 int counter1 = 0;
				 int flag2 = i_count_mox;
				 int flag_end2 = i_count_mox + length;

				 for (int k = flag2; k < flag_end2;k++) {
					 for (int m = 0;m < ncol;m++) {
						 mox[k][m] = mox_recv[counter1][m];
					 }
					 counter1 = counter1 + 1;
					 i_count_mox++;
				 }

			 }
		 }
		 //cout << "zmat debug3 on node " << mynode << endl;

		 MPI_Bcast(&i_count_mox, 1, MPI_INT, 0, MPI_COMM_WORLD);
		 MPI_Bcast(mox[0], (nrow*ncol), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	 } //end of parallel 
	 else {
		 i_count_mox = 0; //total number of unique uox 
		 std::string s_temp;
		 for (int i = 0; i<i_ml_temp; i++)
		 {
			 s_temp = s_ml[i]; //get a string 

			 if (i == 0 || (i > 0 && s_temp.compare(s_ml[i - 1]) != 0)) {
				 //if (mynode == 0) {
					// cout << "III: " << i << endl;
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
	 }//end of serial
	//cout << "final i_count_mox: " << i_count_mox << endl;
	//RPrint("final mox: ");RPrint("at node ");RPrint(mynode);
	//RPrint(mox, nrow, ncol);

	// if (mynode == 0) cout << "ZMAT 3" << endl;
	//Now, i_count_mox means the total number of unique missing patterns
	//cout << "zmat_check4_begin at node " << mynode << " = " << MPI_Wtime() - zmat_check4_begin << endl;
	//----------------
	//additional check for unique observed and missing patterns
	//----------------
	//observed patterns//
	double zmat_check5_begin = MPI_Wtime();
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
	Copy_dMatrix(uox_final, nrow, ncol, uox);
	//cout << "zmat_check5_begin at node " << mynode << " = " << MPI_Wtime() - zmat_check5_begin << endl;

	double zmat_check6_begin = MPI_Wtime();
	//missing patterns//
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
	Copy_dMatrix(mox_final, nrow, ncol, mox);
	//cout << "zmat_check6_begin at node " << mynode << " = " << MPI_Wtime() - zmat_check6_begin << endl;
	//------------------
	//Deallocation
	//------------------
	Del_dMatrix(uox_final, nrow, ncol);
	Del_dMatrix(mox_final, nrow, ncol);
	Del_dMatrix(uox_temp, i_ol_temp, ncol);
	Del_dMatrix(uox_recv, i_ol_temp, ncol);
	Del_dMatrix(mox_temp, i_ml_temp, ncol);
	Del_dMatrix(mox_recv, i_ml_temp, ncol);
	delete[] s_ol;
	delete[] s_ml;

	return;

}

/*
void Trans(double** z, const int nrow, const int ncol, std::string cn[])
//Description=========================================
// make a condensed expression of z
//
// Algorithm:  each row of z will be concatenated as a single string consisting of 35 characters
//
// Note: as of Oct 2016, NA values (missing data) is marked by 125711.131723 at the parent "r" code
//
// original R code: Dr. Im, J. and Dr. Kim, J.
// c++ code: 		Dr. Cho, I.
// All rights reserved
//
// updated: Oct 6, 2016
//----------------------------------------------------
//IN   	: double z(nrow, ncol)  = categorized matrix corresponding to original matrix x
//OUT	: string cn(nrow)		= vector of string to represent each row of z
//====================================================
{
const char ch_db[35] = {'1', '2', '3', '4', '5', '6', '7', '8', '9',
'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',
'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r',
's', 't', 'u', 'v', 'w', 'x', 'y', 'z'};
char ch_temp;
int  i_temp=0;
std::string	s_all;
for(int i_row=0; i_row<nrow; i_row++)
{
s_all = ""; //initialize with empty string
for(int i_col = 0; i_col<ncol; i_col++)
{
i_temp = (int)z[i_row][i_col];

ch_temp = '0'; //default character is zero
if(i_temp>=1 && i_temp<=35)
{
ch_temp = ch_db[i_temp-1];
}
s_all.append(&ch_temp);
}

//---------
//store the condensed string
//---------
cn[i_row] = s_all;
}

return;


}

void Trans1(double* z, const int n, std::string &cn)
//Description=========================================
// make a condensed expression of a double array, z
//
// Algorithm:  z will be concatenated as a single string consisting of 35 characters
//
// original R code: Dr. Im, J. and Dr. Kim, J.
// c++ code: 		Dr. Cho, I.
// All rights reserved
//
// updated: Oct 6, 2016
//----------------------------------------------------
//IN   	: double z(ncol)  =  categorized array corresponding to a row of original matrix x
//OUT	: string cn		  =  a string to represent the given row of z
//====================================================
{
const char ch_db[35] = {'1', '2', '3', '4', '5', '6', '7', '8', '9',
'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',
'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r',
's', 't', 'u', 'v', 'w', 'x', 'y', 'z'};
char ch_temp;
int  i_temp=0;
std::string	s_all;
s_all = ""; //initialize with empty string
for(int i = 0; i<n; i++)
{
i_temp = (int)z[i];

ch_temp = '0'; //default character is zero
if(i_temp>=1 && i_temp<=35)
{
ch_temp = ch_db[i_temp-1];
}
s_all.append(&ch_temp);
}

//---------
//store the condensed string
//---------
cn = s_all;

return;
}
*/

//} //end of namespace