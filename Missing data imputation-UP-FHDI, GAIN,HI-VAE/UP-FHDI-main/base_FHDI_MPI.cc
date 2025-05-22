#include <limits> // For NaN 
#include <iostream> // For cout, cerr, endl, and flush 
#include <assert.h> // For assert
#include <algorithm> // For sort 
#include <string>    //For string array
#include <vector>
#include <map>
#include "findFrequencyUtil.cpp"
//#include <R.h>
//#include <Rmath.h>

#include "base_FHDI_MPI.h"


//------------------------
// Definitions of local base functions
// for FHDI
//------------------------

void table_cpp(std::string cn[], const int nrow,
	std::vector<std::string> &v_table_row1, std::vector<int> &v_table_row2)
	//Description=========================================
	// make a table of given STRING vector 
	//
	// Algorithm: count unique items in the given cn[]  
	// 
	//
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: Oct 6, 2016
	//----------------------------------------------------
	//IN 	: string cn(nrow)		= vector of string to represent each row of z          
	//OUT   : std::vector<string> v_table_row1= 1st row: unique item in string format
	//OUT   : std::vector<int> v_table_row2 = 2nd row: count of the unique item  
	//====================================================
{
	std::string s_temp;
	//std::string cn_temp[nrow]; 
	std::string * cn_temp = new std::string[nrow];
	for (int i = 0; i<nrow; i++) { cn_temp[i] = cn[i]; } //make a copy of original cn[]

														 //-----------
														 //internal sorting of the cn_temp[]
														 //just like "table" of R
														 //-----------
	std::sort(cn_temp, cn_temp + nrow);


	const std::string s_null = ""; //empty string 
	int i_temp = 0;
	for (int i = 0; i<nrow; i++)
	{
		i_temp = 0; //re-initialize
		s_temp = cn_temp[i];
		//-----
		//search s_temp
		//-----
		if (s_temp.compare(s_null) != 0) //NOT an empty cell 
		{
			for (int j = i; j<nrow; j++) //count item including myself
			{
				if (s_temp.compare(cn_temp[j]) == 0) //0: equal string
				{
					i_temp++;    //count the same string in cn 
					if (j>i) cn_temp[j] = s_null; //delete the same string just found 
				}
			}
			//store the found unique string and its count 
			if (i_temp > 0) //there is at least one unique item
			{
				v_table_row1.push_back(s_temp);
				v_table_row2.push_back(i_temp); //actual total number of the unique string 
			}
		}
	}


	delete[] cn_temp;

	return;
}

void table_cpp(double* d_source, const int nrow,
	std::vector<double> &v_table_row1, std::vector<int> &v_table_row2)
	//Description=========================================
	// make a table of given DOUBLE array 
	//
	// Algorithm: count unique items in the given d_source[]  
	// 
	//
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: Oct 6, 2016
	//----------------------------------------------------
	//IN 	: double d_source(nrow)		= double array           
	//OUT   : std::vector<double> v_table_row1= 1st row: unique item in double format
	//OUT   : std::vector<int>    v_table_row2 = 2nd row: count of the unique item  
	//====================================================
{
	double d_temp;
	double* d_source_temp = new double[nrow];
	for (int i = 0; i<nrow; i++) { d_source_temp[i] = d_source[i]; } //make a copy of original array

																	 //-----------
																	 //internal sorting of the cn_temp[]
																	 //just like "table" of R
																	 //-----------
	std::sort(d_source_temp, d_source_temp + nrow);


	int i_temp = 0;

	for (int i = 0; i<nrow; i++)
	{
		i_temp = 0; //re-initialize
		d_temp = d_source_temp[i];
		if (std::isnan(d_temp) != 1) //only meaningful value
		{
			//-----
			//search d_temp
			//-----
			for (int j = i; j<nrow; j++) //count item including myself
			{
				if (fabs(d_temp - d_source_temp[j])<1e-15) //~0: equal value
				{
					i_temp++;    //count the same double in d_soure 
					if (j>i) d_source_temp[j] = nan(""); //delete the same double just found 
				}
			}
			//store the found unique string and its count 
			if (i_temp > 0) //there is at least one unique item
			{
				v_table_row1.push_back(d_temp);
				v_table_row2.push_back(i_temp); //actual total number of the unique double 
			}
		}

	}

	delete[] d_source_temp;
	return;
}

void table_cpp_Yicheng(std::string cn[], const int nrow,
	std::vector<std::string> &v_table_row1, std::vector<int> &v_table_row2)
	//Description=========================================
	// make a table of given STRING vector 
	//
	// Algorithm: count unique items in the given cn[]  
	// 
	//
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: Oct 6, 2016
	//----------------------------------------------------
	//IN 	: string cn(nrow)		= vector of string to represent each row of z          
	//OUT   : std::vector<string> v_table_row1= 1st row: unique item in string format
	//OUT   : std::vector<int> v_table_row2 = 2nd row: count of the unique item  
	//====================================================
{
	//std::string cn_temp[nrow];
	//cout<<"table_cpp_yicheng"<<endl;
	std::string * cn_temp = new std::string[nrow];
	for (int i = 0; i<nrow; i++) { cn_temp[i] = cn[i]; } //make a copy of original cn[]

														 //-----------
														 //internal sorting of the cn_temp[]
														 //just like "table" of R
														 //-----------
	std::sort(cn_temp, cn_temp + nrow);

	int i_temp = 0;
	std::string s_temp;

	for (int i = 0;i<nrow;i++) {
		s_temp = cn_temp[i];
		if (i == 0) {
			v_table_row1.push_back(s_temp);
		}
		if (i > 0 && s_temp.compare(cn_temp[i - 1]) != 0) {
			v_table_row2.push_back(i_temp);
			i_temp = 0;
			v_table_row1.push_back(s_temp);
		}
		i_temp++;
		if (i == (nrow - 1)) {
			v_table_row2.push_back(i_temp);
		}
	}

	delete[] cn_temp;

	return;
}


void table_cpp_Yicheng(double* d_source, const int nrow,
	std::vector<double> &v_table_row1, std::vector<int> &v_table_row2)
	//Description=========================================
	// make a table of given DOUBLE array 
	//
	// Algorithm: count unique items in the given d_source[]  
	// 
	//
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: Oct 6, 2016
	//----------------------------------------------------
	//IN 	: double d_source(nrow)		= double array           
	//OUT   : std::vector<double> v_table_row1= 1st row: unique item in double format
	//OUT   : std::vector<int>    v_table_row2 = 2nd row: count of the unique item  
	//====================================================
{
	double* d_source_temp = new double[nrow];
	for (int i = 0; i<nrow; i++) { d_source_temp[i] = d_source[i]; } //make a copy of original array

																	 //-----------
																	 //internal sorting of the cn_temp[]
																	 //just like "table" of R
																	 //-----------
	std::sort(d_source_temp, d_source_temp + nrow);
	vector<int> freq(d_source_temp[nrow - 1] + 1, 0);

	// Fill the vector with frequency 
	findFrequencyUtil(d_source_temp, 0, nrow - 1, freq);

	// Print the frequencies 
	for (int i = 0; i <= d_source_temp[nrow - 1]; i++) {
		if (freq[i] != 0) {
			//cout << "Element " << i << " occurs " << freq[i] << " times" << endl;
			v_table_row1.push_back(i);
			v_table_row2.push_back(freq[i]);
		}
	}



	delete[] d_source_temp;
	return;
}

//This table function is only used for correlated variables
void table_cpp_yicheng(int* d_source, int nrow,
	std::vector<int> &v_table_row1, std::vector<int> &v_table_row2, ofstream& TestOut)
	//Description=========================================
	// make a table of given DOUBLE array 
	//
	// Algorithm: count unique items in the given d_source[]  
	// 
	//
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: Oct 6, 2016
	//----------------------------------------------------
	//IN 	: double d_source(nrow)		= double array           
	//OUT   : std::vector<double> v_table_row1= 1st row: unique item in double format
	//OUT   : std::vector<int>    v_table_row2 = 2nd row: count of the unique item  
	//====================================================
{
	//TestOut<<"Table_Yicheng"<<endl;
	int d_temp;
	int* d_source_temp = new int[nrow];
	for (int i = 0; i<nrow; i++) { d_source_temp[i] = d_source[i]; } //make a copy of original array

																	 //-----------
																	 //internal sorting of the cn_temp[]
																	 //just like "table" of R
																	 //-----------
	std::sort(d_source_temp, d_source_temp + nrow);


	int i_temp = 0;

	for (int i = 0; i < nrow; i++)
	{
		i_temp = 0; //re-initialize
		d_temp = d_source_temp[i];

		//-----
		//search d_temp
		//-----
		if (d_temp != 1234567890) {
			for (int j = i; j < nrow; j++) //count item including myself
			{
				if (fabs(d_temp - d_source_temp[j]) < 1e-15) //~0: equal value
				{
					i_temp++;    //count the same double in d_soure 
					if (j > i) d_source_temp[j] = 1234567890; //delete the same double just found 
				}
			}
			//store the found unique string and its count 
			if (i_temp > 0) //there is at least one unique item
			{
				v_table_row1.push_back(d_temp);
				v_table_row2.push_back(i_temp); //actual total number of the unique double 
			}
		}

	}
	delete[] d_source_temp;
	return;
}


void Trans(double** z, const int nrow, const int ncol, std::string cn[])
//Description=========================================
// make a condensed expression of z
//
// Algorithm:  each row of z will be concatenated as a single string consisting of 35 characters
// 
// Note: as of Oct 2016, NA values (missing data) is marked by a long integer at the parent "r" code
// Note: as of Apr 2017, the use of combination of char and string appears to cause error
//                       in Ubuntu platform
//                       Hence, a uniform use of string is recommended as below
// original R code: Dr. Im, J. and Dr. Kim, J. 
// c++ code: 		Dr. Cho, I. 
// All rights reserved
// 
// updated: April 4, 2017
//----------------------------------------------------
//IN   	: double z(nrow, ncol)  = categorized matrix corresponding to original matrix x
//OUT	: string cn(nrow)		= vector of string to represent each row of z                                
//====================================================
{
	const std::string ch_db[35] = { "1", "2", "3", "4", "5", "6", "7", "8", "9",
		"a", "b", "c", "d", "e", "f", "g", "h", "i",
		"j", "k", "l", "m", "n", "o", "p", "q", "r",
		"s", "t", "u", "v", "w", "x", "y", "z" };
	std::string ch_temp;
	int  i_temp = 0;

	for (int i_row = 0; i_row<nrow; i_row++)
	{
		std::string	s_all;

		for (int i_col = 0; i_col<ncol; i_col++)
		{
			i_temp = (int)z[i_row][i_col];

			ch_temp = "0"; //default character is zero
			if (i_temp >= 1 && i_temp <= 35)
			{
				ch_temp = ch_db[i_temp - 1];
			}
			s_all.append(ch_temp);
		}

		//---------
		//store the condensed string
		//---------
		cn[i_row] = s_all;
	}

	return;


}

//void Trans_new(double** z, const int nrow, const int ncol, std::string cn[]) {
//
//	int  i_temp = 0;
//	for (int i = 0; i < nrow; i++) {
//		std::string str;
//		for (int j = 0; j < ncol; j++) {
//			i_temp = (int)z[i][j];
//			str = str + std::to_string(i_temp);
//		}
//		cn[i] = str;
//	}
//
//}


void Trans1(double* z, const int n, std::string &cn)
//Description=========================================
// make a condensed expression of a double array, z
//
// Algorithm:  z will be concatenated as a single string consisting of 35 characters
//
// Note: as of Apr 2017, the use of combination of char and string appears to cause error
//                       in Ubuntu platform
//                       Hence, a uniform use of string is recommended as below
// 
// original R code: Dr. Im, J. and Dr. Kim, J. 
// c++ code: 		Dr. Cho, I. 
// All rights reserved
// 
// updated: April 4, 2017
//----------------------------------------------------
//IN   	: double z(n)  =  categorized array corresponding to a row of original matrix x
//OUT	: string cn		  =  a string to represent the given row of z                                
//====================================================
{
	const std::string ch_db[35] = { "1", "2", "3", "4", "5", "6", "7", "8", "9",
		"a", "b", "c", "d", "e", "f", "g", "h", "i",
		"j", "k", "l", "m", "n", "o", "p", "q", "r",
		"s", "t", "u", "v", "w", "x", "y", "z" };
	std::string ch_temp;
	int  i_temp = 0;
	std::string	s_all;

	for (int i = 0; i<n; i++)
	{
		i_temp = (int)z[i];

		ch_temp = "0"; //default character is zero
		if (i_temp >= 1 && i_temp <= 35)
		{
			ch_temp = ch_db[i_temp - 1];
		}
		s_all.append(ch_temp);
	}

	//---------
	//store the condensed string
	//---------
	cn = s_all;

	return;
}

//-------------------------
//local function for "which" of R 
//-------------------------
// return ACTUAL location having the same integer as i_target 
//-------------------------
void which(int* i_vector, const int n, const int i_target, std::vector<int> &v_location)
{
	if (n <= 0) { cout << "Error! n<=0! in which()" << endl; return; }
	for (int i = 0; i<n; i++)
	{
		if (i_vector[i] == i_target) v_location.push_back(i + 1); //actual location 
	}
	return;
}

//-------------------------
// return ACTUAL location having the same double as d_target 
//-------------------------
void which(double* d_vector, const int n, const double d_target, std::vector<int> &v_location)
{
	if (n <= 0) { cout << "Error! n<=0! in which()" << endl; return; }
	for (int i = 0; i<n; i++)
	{
		if (fabs(d_vector[i] - d_target)<1e-15) v_location.push_back(i + 1); //actual location 
	}

	return;
}

//Note that this function is particular for returning ACTUAL GLOBAL location having the same STRING as s_target in nDAU_ultra function
void which(double* d_vector, const int n, const double d_target, int i_recursion, int recursion_size, std::vector<int> &v_location)
{
	if (n <= 0) { cout << "Error! n<=0! in which()" << endl; return; }

	int offset = 0;
	offset = i_recursion*recursion_size;

	for (int i = 0; i<n; i++)
	{
		if (fabs(d_vector[i] - d_target)<1e-15) v_location.push_back(i + offset + 1); //actual location 
	}

	return;
}

//-------------------------
// return ACTUAL location having the same STRING as s_target 
//-------------------------
void which(std::vector<std::string> s_vector, const  std::string s_target,
	std::vector<int> &v_location)
{
	const int n = s_vector.size();
	if (n <= 0) { cout << "Error! n<=0! in which s_vector()" << endl;   return; }
	for (int i = 0; i<n; i++)
	{
		if (s_vector[i].compare(s_target) == 0) //0: equal string
		{
			v_location.push_back(i + 1);
		} //actual location 
	}
	return;
}

//Note that this function is particular for returning ACTUAL GLOBAL location having the same STRING as s_target in nDAU_ultra function
void which(std::vector<std::string> s_vector, const  std::string s_target, int i_recursion, int recursion_size,
	std::vector<int> &v_location)
{
	const int n = s_vector.size();
	if (n <= 0) { cout << "Error! n<=0! in which s_vector()" << endl;   return; }

	int offset = 0;
	offset = i_recursion*recursion_size;

	for (int i = 0; i<n; i++)
	{
		if (s_vector[i].compare(s_target) == 0) //0: equal string
		{
			v_location.push_back(i + offset + 1);
		} //actual location 
	}
	return;
}

//-------------------------
// return ACTUAL location having the same STRING in an ARRAY as s_target 
//-------------------------
void which(std::string s_array[], const int n, const  std::string s_target,
	std::vector<int> &v_location)
{
	if (n <= 0) { cout << "Error! n<=0! in which s_array()" << endl;   return; }
	for (int i = 0; i<n; i++)
	{
		if (s_array[i].compare(s_target) == 0) //0: equal string
		{
			v_location.push_back(i + 1);
		} //actual location 
	}
	return;
}

//-------------------------
// return ACTUAL location having the DIFFERENT integer from i_target 
//-------------------------
void whichINV(int* i_vector, const int n, const int i_target, std::vector<int> &v_location)
{
	if (n <= 0) { cout << "Error! n<=0! in which()" << endl;  return; }
	for (int i = 0; i<n; i++)
	{
		if (i_vector[i] != i_target) v_location.push_back(i + 1); //actual location 
	}
	return;
}

//-------------------------
// return ACTUAL location having the same integer from i_target (Yicheng)
//-------------------------
void whichINVNOT(int* i_vector, const int n, const int i_target, std::vector<int> &v_location)
{
	if (n <= 0) { cout << "Error! n<=0! in which()" << endl;  return; }
	for (int i = 0; i<n; i++)
	{
		if (i_vector[i] == i_target) v_location.push_back(i + 1); //actual location 
	}
	return;
}


//-------------------------
// return ACTUAL location having the DIFFERENT double from d_target 
//-------------------------
void whichINV(double* d_vector, const int n, const double d_target, std::vector<int> &v_location)
{
	if (n <= 0) { cout << "Error! n<=0! in which()" << endl; return; }
	for (int i = 0; i<n; i++)
	{
		if (fabs(d_vector[i] - d_target) > 1e-15) v_location.push_back(i + 1); //actual location 
	}
	return;
}

//------------------------------- 
// Rprint: print out double matrix on R Console
//-------------------------------
void RPrint(double** d_debug, const int nrow, const int ncol)
{
	if (nrow <= 0 || ncol <= 0)
	{
		cout << "Error! nrow or ncol<=0! in printing d_debug[][]" << endl;
	}

	for (int i = 0; i<nrow; i++)
	{
		for (int j = 0; j<ncol; j++)
		{
			cout << d_debug[i][j] << "    ";
		}
		cout << endl;
	}
	return;
}


//------------------------------- 
// Rprint: print out double matrix on debug output file 
//-------------------------------
void RPrint(double** d_debug, const int nrow, const int ncol, ofstream& TestOut)
{
	if (nrow <= 0 || ncol <= 0) { TestOut << "Error! nrow or ncol<=0! in printing d_debug[][]" << endl;  return; }

	for (int i = 0; i<nrow; i++)
	{
		for (int j = 0; j<ncol; j++)
		{
			TestOut << setw(10) << d_debug[i][j];
		}
		TestOut << endl;
	}

	return;
}

void RPrint_Yicheng(double** d_debug, const int nrow, const int ncol, ofstream& TestOut1)
{
	if (nrow <= 0 || ncol <= 0) { TestOut1 << "Error! nrow or ncol<=0! in printing d_debug[][]" << endl;  return; }

	for (int i = 0; i<nrow; i++)
	{
		for (int j = 0; j<ncol; j++)
		{
			TestOut1 << setw(10) << d_debug[i][j];
		}
		TestOut1 << endl;
	}

	return;
}

//------------------------------- 
// Rprint: print out double vector on R console
//-------------------------------
void RPrint(double* d_debug, const int n)
{
	if (n <= 0) { cout << "Error! n<=0! in printing d_debug[]" << endl;  return; }

	for (int i = 0; i<n; i++) { cout << "%g " << d_debug[i]; cout << "      " << endl; }


	return;
}
void RPrint_Yicheng(double* d_debug, const int n)
{
	if (n <= 0) { cout << "Error! n<=0! in printing d_debug[]" << endl;  return; }

	for (int i = 0; i<n; i++) { cout << " " << d_debug[i]; }


	return;
}

//------------------------------- 
// Rprint: print out double vector on Output file
//-------------------------------
void RPrint(double* d_debug, const int n, ofstream& TestOut)
{
	if (n <= 0) { TestOut << "Error! n<=0! in printing d_debug[]" << endl;  return; }

	for (int i = 0; i<n; i++) { TestOut << setw(10) << d_debug[i] << endl; }


	return;
}

//------------------------------- 
// Rprint: print out integer vector on R console
//-------------------------------
void RPrint(int* i_debug, const int n)
{
	if (n <= 0) { cout << "Error! n<=0! in printing i_debug[]" << endl;   return; }
	for (int i = 0; i<n; i++) { cout << "%d " << i_debug[i]; cout << "      " << endl; }


	return;
}

//------------------------------- 
// Rprint: print out integer vector on Output file
//-------------------------------
void RPrint(int* i_debug, const int n, ofstream& TestOut)
{
	if (n <= 0) { TestOut << "Error! n<=0! in printing i_debug[]" << endl;   return; }
	for (int i = 0; i<n; i++) { TestOut << setw(10) << i_debug[i] << endl; }


	return;
}

//------------------------------- 
// Rprint: print out vector of integer on R console
//-------------------------------
void RPrint(std::vector<int> i_debug)
{
	const int n = i_debug.size();
	if (n <= 0) { cout << "Error! n<=0! in vector<int>" << endl;   return; }

	for (int i = 0; i<n; i++) { cout << "%d " << i_debug[i]; cout << "      " << endl; }


	return;
}

//------------------------------- 
// Rprint: print out vector of integer on Output File
//-------------------------------
void RPrint(std::vector<int> i_debug, ofstream& TestOut)
{
	const int n = i_debug.size();
	if (n <= 0) { TestOut << "Error! n<=0! in vector<int>" << endl;   return; }

	for (int i = 0; i<n; i++) { TestOut << setw(10) << i_debug[i] << endl; }


	return;
}

//------------------------------- 
// Rprint: print out vector of double on R console
//-------------------------------
void RPrint(std::vector<double> d_debug)
{
	const int n = d_debug.size();
	if (n <= 0) { cout << "Error! n<=0! in vector<double>" << endl;   return; }

	for (int i = 0; i<n; i++) { cout << "%g " << d_debug[i]; cout << "      " << endl; }

	return;
}

//------------------------------- 
// Rprint: print out vector of double on Output File
//-------------------------------
void RPrint(std::vector<double> d_debug, ofstream& TestOut)
{
	const int n = d_debug.size();
	if (n <= 0) { TestOut << "Error! n<=0! in vector<double>" << endl;   return; }

	for (int i = 0; i<n; i++) { TestOut << setw(10) << d_debug[i] << endl; }

	return;
}

//------------------------------- 
// Rprint: print out double array on Output File
//-------------------------------
void RPrint_Yicheng_Output(double* d_debug, const int n, ofstream& TestOut)
{
	if (n <= 0) { cout << "Error! n<=0! in printing d_debug[]" << endl;  return; }

	for (int i = 0; i<n; i++) { TestOut << setw(20) << d_debug[i]; }
	TestOut << endl;

	return;
}

//------------------------------- 
// Rprint: print out string array on R console
//-------------------------------
void RPrint(std::string s_debug[], const int n)
{
	if (n <= 0) { cout << "Error! n<=0! in string[]" << endl;   return; }

	for (int i = 0; i<n; i++)
	{
		const char * ch_temp = s_debug[i].c_str();
		cout << "%s " << ch_temp; cout << "      " << endl;

	}

	return;
}

//------------------------------- 
// Rprint: print out string array on Output File
//-------------------------------
void RPrint(std::string s_debug[], const int n, ofstream& TestOut)
{
	if (n <= 0) { TestOut << "Error! n<=0! in string[]" << endl;   return; }

	for (int i = 0; i<n; i++)
	{
		const char * ch_temp = s_debug[i].c_str();
		TestOut << setw(10) << ch_temp << endl;

	}

	return;
}

//------------------------------- 
// Rprint: print out string vector on R console
//-------------------------------
void RPrint(std::vector<std::string> v_sdebug)
{
	const int n = (int)v_sdebug.size();
	if (n <= 0) { cout << "Error! n<=0! in string[]" << endl;   return; }

	for (int i = 0; i<n; i++)
	{
		const char * ch_temp = v_sdebug[i].c_str();
		cout << "%s " << ch_temp; cout << "      " << endl;
		//cout<<"%s ", s_debug[i]); 

	}


	return;
}

//------------------------------- 
// Rprint: print out string vector on Output File
//-------------------------------
void RPrint(std::vector<std::string> v_sdebug, ofstream& TestOut)
{
	const int n = (int)v_sdebug.size();
	if (n <= 0) { TestOut << "Error! n<=0! in string[]" << endl;   return; }

	for (int i = 0; i<n; i++)
	{
		const char * ch_temp = v_sdebug[i].c_str();
		TestOut << setw(10) << ch_temp << endl;
		//cout<<"%s ", s_debug[i]); 

	}


	return;
}

//------------------------------- 
// Rprint: print one integer on R console
//-------------------------------
void RPrint(const int i_target)
{
	cout << "%d " << i_target;

	return;
}

//------------------------------- 
// Rprint: print one integer on Output File
//-------------------------------
void RPrint(const int i_target, ofstream& TestOut)
{
	TestOut << setw(10) << i_target << endl;

	return;
}

//------------------------------- 
// Rprint: print one double on R console
//-------------------------------
void RPrint(const double d_target)
{
	cout << "%g " << d_target;

	return;
}

//------------------------------- 
// Rprint: print one double on Output File
//-------------------------------
void RPrint(const double d_target, ofstream& TestOut)
{
	TestOut << setw(10) << d_target << endl;

	return;
}

//------------------------------- 
// Rprint: print out string on R console
//-------------------------------
void RPrint(const char *vString)
{
	cout << "%s" << vString;   return;
}

//------------------------------- 
// Rprint: print out string on Output File
//-------------------------------
void RPrint(const char *vString, ofstream& TestOut)
{
	TestOut << setw(10) << vString << endl;   return;
}

//-------------------------------
//basic tools for vector, array
//-------------------------------
int sum_FHDI(std::vector<int> i_source)
{
	int i_sum = 0;
	int i_n = (int)i_source.size();
	for (int i = 0; i<i_n; i++) i_sum += i_source[i];

	return i_sum;
}
int sum_FHDI(int* i_source, const int n_size)
{
	int i_sum = 0;
	int i_n = n_size;
	for (int i = 0; i<i_n; i++) i_sum += i_source[i];

	return i_sum;
}
//---------------------------------
//max value of int vector
//---------------------------------
int max_FHDI(std::vector<int> i_source)
{
	int max = i_source[0];
	int i_n = (int)i_source.size();
	for (int i = 0; i<i_n; i++) { if (max < i_source[i]) max = i_source[i]; }

	return max;
}
//---------------------------------
//max value of double vector
//---------------------------------
double max_FHDI(std::vector<double> i_source)
{
	double max = i_source[0];
	int i_n = (int)i_source.size();
	for (int i = 0; i<i_n; i++) { if (max < i_source[i]) max = i_source[i]; }

	return max;
}
//---------------------------------
//min value of int vector
//---------------------------------
int min_FHDI(std::vector<int> i_source)
{
	int min = i_source[0];
	int i_n = (int)i_source.size();
	for (int i = 0; i<i_n; i++) { if (min>i_source[i]) min = i_source[i]; }

	return min;
}
//---------------------------------
//min value of double vector
//---------------------------------
double min_FHDI(std::vector<double> i_source)
{
	double min = i_source[0];
	int i_n = (int)i_source.size();
	for (int i = 0; i<i_n; i++) { if (min>i_source[i]) min = i_source[i]; }

	return min;
}
//---------------------------------
//max value of double array
//---------------------------------
double max_FHDI(double* k, const int n)
{
	double max_k = k[0];
	for (int i = 0; i<n; i++)
	{
		if (max_k < k[i]) max_k = k[i];
	}

	return max_k;
}
//---------------------------------
//min value of double array
//---------------------------------
double min_FHDI(double* k, const int n)
{
	double min_k = k[0];
	for (int i = 0; i<n; i++)
	{
		if (min_k > k[i]) min_k = k[i];
	}

	return min_k;
}
//----------------------------------
//second min value of double array written by Yicheng; Note that if arr = {1.1, 2.2, 1.1, 2.3}, it will return 2.2
//---------------------------------
double second_min_FHDI(double arr[], int n) {
	double smallest, secondSmallest;
	if ((arr[0] - arr[1]) < -1e-15) {
		smallest = arr[0];
		secondSmallest = arr[1];
	}

	if (fabs(arr[0] - arr[1]) < 1e-15) {
		smallest = arr[0];
		for (int j = 0; j < n;j++) {
			if (arr[j] != smallest) {
				secondSmallest = arr[j];
			}
		}
	}
	else {
		smallest = arr[1];
		secondSmallest = arr[0];
	}
	for (int i = 0; i < n; i++) {

		if ((smallest - arr[i]) > 1e-15) {
			secondSmallest = smallest;
			smallest = arr[i];
		}

		else if (((arr[i] - secondSmallest) < -1e-15) && ((arr[i] - smallest) > 1e-15)) {
			secondSmallest = arr[i];
		}
	}
	return secondSmallest;
}
//---------------------------------
//max value of integer array
//---------------------------------
int max_FHDI(int* k, const int n)
{
	int max_k = k[0];
	for (int i = 0; i<n; i++)
	{
		if (max_k < k[i]) max_k = k[i];
	}

	return max_k;
}
//---------------------------------
//min value of integer array
//---------------------------------
int min_FHDI(int* k, const int n)
{
	int min_k = k[0];
	for (int i = 0; i<n; i++)
	{
		if (min_k > k[i]) min_k = k[i];
	}

	return min_k;
}

//--------------------------------
//calculate absolute distance^2 between Matrix entities and a double 
//--------------------------------
void distance2(double** d_mat, const int nrow, const int ncol, const double d_origin,
	double* d_distance)
	//Description----------------------------------------
	//calculate the absolute distance^2 between all entities of the given matrix 
	// and a given origin
	//IN   : double d_mat(nrow, ncol)   = source matrix with double values
	//IN   : double d_origin 		 	= origin 
	//OUT  : double d_distance(nrow) = sum(|a - b|^2) per row 
	//----------------------------------------------------
{
	Fill_dVector(d_distance, nrow, 0.0);
	double d_sum = 0.0;
	for (int i = 0; i<nrow; i++)
	{
		d_sum = 0.0; //reinitialization
		for (int j = 0; j<ncol; j++)
		{
			d_sum += (d_mat[i][j] - d_origin)*(d_mat[i][j] - d_origin);
		}
		d_distance[i] = d_sum;
	}
	return;
}


void order_FHDI(int* i_original, const int n)
//Description ================================
// Order the positive integer array in ascending order
//
//INOUT   : int i_original_0(n) returned with the ordered (Actual) cell numbers     
//          i_original > 0
//=============================================
{

	int* i_source = new int[n];
	int* i_order = new int[n];

	for (int i = 0; i<n; i++)
	{
		i_source[i] = i_original[i]; //backup
		i_order[i] = i + 1; //default
	}

	//-----------
	//leverage sorting library
	//-----------
	std::sort(i_source, i_source + n);
	int i_now = 0;

	i_order[0] = 1; //first cell location as default
	for (int i = 0; i<n; i++)
	{
		i_now = i_source[i];
		//----------------
		//comparisons from the first entiry to now 
		//----------------
		for (int j = 0; j<n; j++)
		{
			if (fabs(i_now - i_original[j])<1e-3)
			{
				i_order[i] = j + 1; //Actual location
				i_original[j] = -1; //dummy value
				break;
			}
		}
	}
	//---prep return
	for (int i = 0; i<n; i++)
	{
		i_original[i] = i_order[i]; //backup
	}

	delete[] i_source;
	delete[] i_order;

	return;
}

void order_FHDI(double* d_original_0, const int n, int* i_return)
//Description ================================
// Order the positive double-precision array in ascending order
//
//IN   : double d_original_0(n) = original array of double-precision float numbers
//              d_original > 0.0 
//OUT  : int i_return(n)    = returned with the ordered (Actual) cell numbers..]   
//
//=============================================
{
	//Note: below backup is different from integer version
	double* d_original = new double[n]; //backup
	Copy_dVector(d_original_0, n, d_original);

	double* d_source = new double[n];
	int* i_order = new int[n];

	for (int i = 0; i<n; i++)
	{
		d_source[i] = d_original[i]; //backup
		i_order[i] = i + 1; //default
	}

	//-----------
	//leverage sorting library
	//-----------
	std::sort(d_source, d_source + n);
	double d_now = 0;

	i_order[0] = 1; //first cell location as default
	for (int i = 0; i<n; i++)
	{
		d_now = d_source[i];
		//----------------
		//comparisons from the first entiry to now 
		//----------------
		for (int j = 0; j<n; j++)
		{
			if (fabs(d_now - d_original[j])<1e-15)
			{
				i_order[i] = j + 1; //Actual location
				d_original[j] = -1.0; //dummy value
				break;
			}
		}
	}
	//---prep return
	for (int i = 0; i<n; i++)
	{
		i_return[i] = i_order[i]; //backup
	}

	delete[] d_original;
	delete[] d_source;
	delete[] i_order;

	return;
}

void order_FHDI(double* d_original_0, const int n, std::vector<int> &i_return)
//Description ================================
// Order the positive double-precision array in ascending order
//
//IN   : double d_original_0(n) = original array of double-precision float numbers
//              d_original > 0.0 
//OUT  : int i_return(n)    = returned with the ordered (Actual) cell numbers..]   
//
//=============================================
{
	//Note: below backup is different from integer version
	double* d_original = new double[n]; //backup
	Copy_dVector(d_original_0, n, d_original);

	double* d_source = new double[n];
	int* i_order = new int[n];

	for (int i = 0; i<n; i++)
	{
		d_source[i] = d_original[i]; //backup
		i_order[i] = i + 1; //default
	}

	//-----------
	//leverage sorting library
	//-----------
	std::sort(d_source, d_source + n);
	double d_now = 0;

	i_order[0] = 1; //first cell location as default
	for (int i = 0; i<n; i++)
	{
		d_now = d_source[i];
		//----------------
		//comparisons from the first entiry to now 
		//----------------
		for (int j = 0; j<n; j++)
		{
			if (fabs(d_now - d_original[j])<1e-15)
			{
				i_order[i] = j + 1; //Actual location
				d_original[j] = -1.0; //dummy value
				break;
			}
		}
	}
	//---prep return
	for (int i = 0; i<n; i++)
	{
		i_return.push_back(i_order[i]); //backup
	}

	delete[] d_original;
	delete[] d_source;
	delete[] i_order;

	return;
}

void order_FHDI(double* d_original_0, int* i_original_0, const int n, std::vector<int> &i_return)
//Description ================================
// Order the positive double-precision array in ascending order
// Note that actual index of this array is given in i_original_0

//IN   : double d_original_0(n) = original array of double-precision float numbers
//              d_original > 0.0 
//IN   : double d_original_0(n) = original actual index of double array

//OUT  : int i_return(n)    = returned with the ordered (Actual) cell numbers..]   
//
//=============================================
{
	//Note: below backup is different from integer version
	double* d_original = new double[n]; //backup
	Copy_dVector(d_original_0, n, d_original);

	int* i_original = new int[n]; //backup
	Copy_iVector(i_original_0, n, i_original);

	double* d_source = new double[n];
	int* i_order = new int[n];

	for (int i = 0; i<n; i++)
	{
		d_source[i] = d_original[i]; //backup
		i_order[i] = 0; //default
	}

	//-----------
	//leverage sorting library
	//-----------
	std::sort(d_source, d_source + n);
	double d_now = 0;

	for (int i = 0; i<n; i++)
	{
		d_now = d_source[i];
		//----------------
		//comparisons from the first entiry to now 
		//----------------
		for (int j = 0; j<n; j++)
		{
			if (fabs(d_now - d_original[j])<1e-15)
			{
				i_order[i] = i_original[j]; //Note that i_original has actula locations already, no need to +1
				d_original[j] = -1.0; //dummy value
				break;
			}
		}
	}
	//---prep return
	for (int i = 0; i<n; i++)
	{
		//if (i_order[i] == 0) cout << "ERROR!!!! in order_FHDI in ranking_m for ultra data" << endl;
		i_return.push_back(i_order[i]); //backup
	}

	delete[] d_original;
	delete[] i_original;
	delete[] d_source;
	delete[] i_order;

	return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void wpct_FHDI(std::string s_0[], const int n, const double* w,
	std::vector<std::string> &jp_name, std::vector<double> &jp_prob)
	//Description=====================================
	//  calculate weighted probability of the string array 
	//  using the given weight array in w[]
	//
	//  written by Dr I. Cho
	//  All right reserved
	//
	//  Algorithm: similar to "R" wpct()
	//
	//IN   : string s_0[n] 	= target array of string
	//IN   : double w[n]  	= user-defined weight used for proportional weights
	//OUT  : std::vector<std::string> jp_name  = names of joint probability table
	//OUT  : std::vector<double>      jp_prob  = weighted joint probability of the table 
	//================================================
{

	//---------------
	//make a table of s_0[n]
	//---------------
	std::vector<std::string> v_table_row1; //names of the table
	std::vector<int> 		 v_table_row2; //counts of the table
	table_cpp(s_0, n, v_table_row1, v_table_row2);
	const int i_size_v_table = (int)v_table_row2.size();

	//---------------
	//find new accumulated weights for each category
	//---------------
	double* d_weight = new double[i_size_v_table];
	Fill_dVector(d_weight, i_size_v_table, 0.0);

	std::string s_temp;
	int i_count = 0;
	for (int i = 0; i<i_size_v_table; i++) //loop for table names 
	{
		s_temp = v_table_row1[i];
		i_count = 0; //re-initialize 

					 //-----------
					 //search and get the weight of current string
					 //-----------
		for (int j = 0; j<n; j++)
		{
			if (s_temp.compare(s_0[j]) == 0) //0 means equal string
			{
				d_weight[i] = d_weight[i] + w[j];  //accumulate the weight of this category
				i_count++;
				if (i_count == v_table_row2[i]) { break; }
			}
		}
	}

	//-----------------
	//sum of d_weight 
	//-----------------
	double d_sum_w = 0.0;
	for (int i = 0; i<i_size_v_table; i++) d_sum_w += d_weight[i];
	if (d_sum_w == 0.0) { cout << "Error! zero sum of weights in wpct from base_FHDI_MPI.cc" << endl; return; }


	//------------------
	//prep return
	//------------------
	for (int i = 0; i<i_size_v_table; i++)
	{
		jp_name.push_back(v_table_row1[i]);
		jp_prob.push_back(d_weight[i] / d_sum_w);
	}

	//------------------
	//Deallocation
	//------------------
	delete[] d_weight;

}

void wpct_FHDI_Yicheng(std::string s_0[], const int n, const double* w,
	std::vector<std::string> &jp_name, std::vector<double> &jp_prob)
	//Description=====================================
	//  calculate weighted probability of the string array 
	//  using the given weight array in w[]
	//
	//  written by Dr I. Cho
	//  All right reserved
	//
	//  Algorithm: similar to "R" wpct()
	//
	//IN   : string s_0[n] 	= target array of string
	//IN   : double w[n]  	= user-defined weight used for proportional weights
	//OUT  : std::vector<std::string> jp_name  = names of joint probability table
	//OUT  : std::vector<double>      jp_prob  = weighted joint probability of the table 
	//================================================
{
	//---------------
	//make a table of s_0[n]
	//---------------
	std::vector<std::string> v_table_row1; //names of the table
	std::vector<int> 		 v_table_row2; //counts of the table

	table_cpp_Yicheng(s_0, n, v_table_row1, v_table_row2);
	const int i_size_v_table = (int)v_table_row2.size();

	//---------------
	//find new accumulated weights for each category
	//---------------
	double* d_weight = new double[i_size_v_table];
	Fill_dVector(d_weight, i_size_v_table, 0.0);

	std::string s_temp;
	int i_count = 0;
	for (int i = 0; i<i_size_v_table; i++) //loop for table names 
	{
		s_temp = v_table_row1[i];
		i_count = 0; //re-initialize 

					 //-----------
					 //search and get the weight of current string
					 //-----------
		for (int j = 0; j<n; j++)
		{
			if (s_temp.compare(s_0[j]) == 0) //0 means equal string
			{
				d_weight[i] = d_weight[i] + w[j];  //accumulate the weight of this category
				i_count++;
				if (i_count == v_table_row2[i]) { break; }
			}
		}
	}

	//-----------------
	//sum of d_weight 
	//-----------------
	double d_sum_w = 0.0;
	for (int i = 0; i<i_size_v_table; i++) d_sum_w += d_weight[i];
	if (d_sum_w == 0.0) { cout << "Error! zero sum of weights in wpct from base_FHDI_MPI.cc" << endl; return; }


	//------------------
	//prep return
	//------------------
	for (int i = 0; i<i_size_v_table; i++)
	{
		jp_name.push_back(v_table_row1[i]);
		jp_prob.push_back(d_weight[i] / d_sum_w);
	}

	//------------------
	//Deallocation
	//------------------
	delete[] d_weight;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void cov_FHDI(double** x, const int nrow, const int ncol, double** cov)
//Description=====================================
//  calculate covariance in a column-to-column manner
//  Note: this generate the "estimated covariance" NOT "population covariance"
//        thus, at the end, 1/(n-1) not 1/n
//  The same as var() of "R"
//
//  written by Dr I. Cho
//  All right reserved
//
//
//IN   : double x[nrow, ncol]  =origianl matrix 
//OUT  : double cov[ncol, ncol]= covariance matrix. cov[0][1] means cov of column 0 and col 1
//================================================
{
	double* x1 = new double[nrow];
	double* x2 = new double[nrow];
	double d_sum = 0.0;
	Fill_dMatrix(cov, ncol, ncol, 0.0);

	//----------
	//get ready two columns
	//----------
	for (int j = 0; j<ncol; j++) //from the first column to the second last column
	{
		for (int j_next = j; j_next<ncol; j_next++) //next column including itself 
		{
			for (int i = 0; i<nrow; i++)
			{
				x1[i] = x[i][j];   //jth column
				x2[i] = x[i][j_next];//next column
			}

			//---
			//each column's mean
			//---
			double x1_mean = 0.0; double x2_mean = 0.0;
			for (int i = 0; i<nrow; i++)
			{
				x1_mean += x1[i];   //jth column
				x2_mean += x2[i];   //next column
			}
			x1_mean = x1_mean / nrow;
			x2_mean = x2_mean / nrow;

			//-----
			//calculate covariance of two columns
			//-----
			d_sum = 0.0;
			for (int i_1 = 0; i_1<nrow; i_1++)
			{
				d_sum += (x1[i_1] - x1_mean)*(x2[i_1] - x2_mean);
			}
			d_sum = d_sum / (nrow - 1);

			//---------
			//store covariance using symmetry property
			//---------
			cov[j][j_next] = d_sum;
			cov[j_next][j] = d_sum;
		}
	}


	//---------
	//Deallocation
	//---------
	delete[] x1;
	delete[] x2;

	return;
}

void correlation_FHDI(double** x, const int nrow, const int ncol, double** cov)
//Description=====================================
//  calculate covariance in a column-to-column manner
//  Note: this generate the "estimated covariance" NOT "population covariance"
//        thus, at the end, 1/(n-1) not 1/n
//  The same as var() of "R"
//
//  written by Dr I. Cho
//  All right reserved
//
//
//IN   : double x[nrow, ncol]  =origianl matrix 
//OUT  : double cov[ncol, ncol]= covariance matrix. cov[0][1] means cov of column 0 and col 1
//================================================
{
	double* x1 = new double[nrow];
	double* x2 = new double[nrow];
	double d_sum = 0.0;
	Fill_dMatrix(cov, ncol, ncol, 0.0);

	//----------
	//get ready two columns
	//----------
	for (int j = 0; j<ncol; j++) //from the first column to the second last column
	{
		for (int j_next = j; j_next<ncol; j_next++) //next column including itself 
		{
			for (int i = 0; i<nrow; i++)
			{
				x1[i] = x[i][j];   //jth column
				x2[i] = x[i][j_next];//next column
			}

			//---
			//each column's mean
			//---
			double x1_mean = 0.0; double x2_mean = 0.0;
			for (int i = 0; i<nrow; i++)
			{
				x1_mean += x1[i];   //jth column
				x2_mean += x2[i];   //next column
			}
			x1_mean = x1_mean / nrow;
			x2_mean = x2_mean / nrow;

			//-----
			//calculate covariance of two columns
			//-----
			d_sum = 0.0;
			for (int i_1 = 0; i_1<nrow; i_1++)
			{
				d_sum += (x1[i_1] - x1_mean)*(x2[i_1] - x2_mean);
			}

			//----------------
			//calculate variance of two columns
			//----------------
			double x1_var = 0.0; double x2_var = 0.0;
			double var_sum = 0.0;
			for (int i_2 = 0; i_2 < nrow;i_2++) {
				var_sum = var_sum + (x1[i_2] - x1_mean)*(x1[i_2] - x1_mean);
			}
			x1_var = var_sum;

			var_sum = 0.0;
			for (int i_3 = 0; i_3 < nrow;i_3++) {
				var_sum = var_sum + (x2[i_3] - x2_mean)*(x2[i_3] - x2_mean);
			}
			x2_var = var_sum;
			//---------
			//store covariance using symmetry property
			//---------
			cov[j][j_next] = d_sum / sqrt(x1_var* x2_var);
			cov[j_next][j] = d_sum / sqrt(x1_var* x2_var);
		}

	}


	//---------
	//Deallocation
	//---------
	delete[] x1;
	delete[] x2;

	return;
}

void match_FHDI(std::string cn[], const int nrow,
	std::string cn_large[], const int nrow_large,
	std::vector<int> &v_match)
	//Description=========================================
	// find a vector of the positions of first matches of cn in cn_large 
	//
	// Algorithm: the same as "match() in R"  
	// 
	//
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: Nov 10, 2016
	//----------------------------------------------------
	//IN 	: string cn(nrow)		        = vector of string 
	//IN 	: string cn_large(nrow_large)	= large vector of strings 
	//
	//OUT   : std::vector<int> v_match = ACTUAL positions of the first matches 
	//====================================================
{
	std::string s_temp;

	const std::string s_null = ""; //empty string 
	for (int i = 0; i<nrow; i++)
	{
		s_temp = cn[i];
		//-----
		//search s_temp
		//-----
		if (s_temp.compare(s_null) != 0) //NOT an empty cell 
		{
			for (int j = 0; j<nrow_large; j++) //find the first match in cn_large
			{
				if (s_temp.compare(cn_large[j]) == 0) //0: equal string
				{
					v_match.push_back(j + 1); //+1 for actual location
					break;
				}
			}
		}
	}

	return;
}

void match_FHDI(std::string cn[], const int nrow,
	std::vector<std::string> v_cn_large,
	std::vector<int> &v_match)
	//Description=========================================
	// find a vector of the positions of first matches of cn in cn_large 
	//
	// Algorithm: the same as "match() in R"  
	// 
	//
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: Nov 10, 2016
	//----------------------------------------------------
	//IN 	: string cn(nrow)		        = vector of string 
	//IN 	: std::vector<std::string> v_cn_large(nrow_large)	= large vector of strings 
	//
	//OUT   : std::vector<int> v_match = ACTUAL positions of the first matches 
	//====================================================
{
	std::string s_temp, s_temp_large;
	const int nrow_large = (int)v_cn_large.size();

	const std::string s_null = ""; //empty string 
	for (int i = 0; i<nrow; i++)
	{
		s_temp = cn[i];
		//-----
		//search s_temp
		//-----
		if (s_temp.compare(s_null) != 0) //NOT an empty cell 
		{
			for (int j = 0; j<nrow_large; j++) //find the first match in cn_large
			{
				s_temp_large = v_cn_large[j];
				if (s_temp.compare(s_temp_large) == 0) //0: equal string
				{
					v_match.push_back(j + 1); //+1 for actual location
					break;
				}
			}
		}
	}

	return;
}


void match_FHDI(std::vector<int> v_cn, std::vector<int> v_cn_large,
	std::vector<int> &v_match)
	//Description=========================================
	// find a vector of the positions of first matches of cn in cn_large 
	//
	// Algorithm: the same as "match() in R"  
	// 
	//
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: Nov 10, 2016
	//----------------------------------------------------
	//IN 	: std::vector<int> v_cn(nrow)   = vector of integer 
	//IN 	: std::vector<int> v_cn_large(nrow_large)	= large vector of integer 
	//
	//OUT   : std::vector<int> v_match = ACTUAL positions of the first matches 
	//====================================================
{
	int i_temp;
	const int nrow = (int)v_cn.size();
	const int nrow_large = (int)v_cn_large.size();

	for (int i = 0; i<nrow; i++)
	{
		i_temp = v_cn[i];
		for (int j = 0; j<nrow_large; j++) //find the first match in cn_large
		{
			if (i_temp == v_cn_large[j])
			{
				v_match.push_back(j + 1); //+1 for actual location
				break;
			}
		}
	}

	return;
}

void match_FHDI(int* i_cn, const int nrow, int* i_cn_large, const int nrow_large,
	std::vector<int> &v_match)
	//Description=========================================
	// find a vector of the positions of first matches of cn in cn_large 
	//
	// Algorithm: the same as "match() in R"  
	// 
	//
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: Nov 10, 2016
	//----------------------------------------------------
	//IN 	: int cn(nrow)		        = vector of integer
	//IN 	: int cn_large(nrow_large)	= large vector of integers 
	//
	//OUT   : std::vector<int> v_match = ACTUAL positions of the first matches 
	//====================================================
{
	int i_temp;

	for (int i = 0; i<nrow; i++)
	{
		i_temp = i_cn[i];
		//-----
		//search i_temp
		//-----
		for (int j = 0; j<nrow_large; j++) //find the first match in cn_large
		{
			if (i_temp == i_cn_large[j])
			{
				v_match.push_back(j + 1); //+1 for actual location
				break;
			}
		}
	}

	return;
}

void match_FHDI(double* d_cn, const int nrow, double* d_cn_large, const int nrow_large,
	std::vector<int> &v_match)
	//Description=========================================
	// find a vector of the positions of first matches of cn in cn_large 
	//
	// Algorithm: the same as "match() in R"  
	// 
	//
	// c++ code: 		Dr. Cho, I. 
	// All rights reserved
	// 
	// updated: Nov 10, 2016
	//----------------------------------------------------
	//IN 	: double cn(nrow)	        = vector of double
	//IN 	: double cn_large(nrow_large)	= large vector of doubles 
	//
	//OUT   : std::vector<int> v_match = ACTUAL positions of the first matches 
	//====================================================
{
	double d_temp;

	for (int i = 0; i<nrow; i++)
	{
		d_temp = d_cn[i];
		//-----
		//search d_temp
		//-----
		for (int j = 0; j<nrow_large; j++) //find the first match in cn_large
		{
			if (fabs(d_temp - d_cn_large[j])<1e-15)
			{
				v_match.push_back(j + 1); //+1 for actual location
				break;
			}
		}
	}

	return;
}

void cumsum_FHDI(double* d_original, const int n, double* d_return)
//Description=========================================
// return cumulative sum of the original elements  
//
// Algorithm: the same as "cumsum() in R"  
// 
//
// c++ code: 		Dr. Cho, I. 
// All rights reserved
// 
// updated: Nov 10, 2016
//----------------------------------------------------
//IN 	: double d_original(n)	= original double
//OUT 	: double d_return(n)	= cumlative summation of elements  
//====================================================
{
	double d_sum = 0.0;

	for (int i = 0; i<n; i++)
	{
		d_sum += d_original[i];
		d_return[i] = d_sum;
	}

	return;
}


void normalization(double** daty, double** daty_normalized, int nrow, int ncol)
//Description=========================================
// Normalize the target matrix daty to 0~1 range column-wisely
//
// c++ code: Yicheng Yang
// All rights reserved
// 
// updated: Dec 30, 2020
//----------------------------------------------------
//IN 	: double** daty  	        = input matrix
//IN 	: int nrow	                = number of rows of input matrix
//IN 	: int ncol	                = number of columns of input matrix
//
//OUT   : double** daty_normalized  = normalized matrix 
//====================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	double* min = new double[ncol];
	double* max = new double[ncol];
	std::vector<double> buffer;// hold one column of daty

	for (int j = 0; j < ncol; j++) {
		buffer.clear();
		for (int i = 0; i < nrow; i++) {
			buffer.push_back(daty[i][j]);
		}

		//if (mynode == 1) {
		//	cout<<"buffer at column "<<j<<endl;
		//	for (int t = 0; t < buffer.size(); t++) {
		//		cout << setw(20) << buffer[t];
		//	}
		//	cout << endl;
		//}

		min[j] = min_FHDI(buffer);
		max[j] = max_FHDI(buffer);
	}

	//for (int i = 0; i < ncol; i++) {
	//	if (fabs(max[i] - min[i]) < 1e-15) {
	//		cout << "ERROR in normalization function such that min = max" << endl;
	//	}
	//}
	//if (mynode == 1) {
	//	cout<<"min: "<<endl;
	//	for (int t = 0; t < ncol; t++) {
	//		cout << setw(20) << min[t];
	//	}
	//	cout << endl;

	//	cout << "max: " << endl;
	//	for (int t = 0; t < ncol; t++) {
	//		cout << setw(20) << max[t];
	//	}
	//	cout << endl;

	//}

	//for (int p = 0; p < nrow; p++) {
	//	for (int j = 0; j < ncol; j++) {
	//		daty_normalized[p][j] = (daty[p][j] - min[j]) / (max[j] - min[j]);
	//	}
	//}

	for (int j = 0; j < ncol; j++) {
		for (int i = 0; i < nrow; i++) {
			if (fabs(max[j] - min[j]) >= 1e-15) {
				daty_normalized[i][j] = (daty[i][j] - min[j]) / (max[j] - min[j]);
			}

			if (fabs(max[j] - min[j]) < 1e-15) {
				daty_normalized[i][j] = 0.0;
			}
		}
	}

	delete[] min;
	delete[] max;

	return;
}


//-----------------------------------------------------
//-----------------------------------------------------
//for replacing unneccessarily large matrix d_rw[][]
//Jan 18, 2019
//----------------------------------------------------- 
class RepWeight_FHDI {
public:
	int size_row() const { return _size_row; }

	//----------
	//(i,j) operator overloading
	//----------
	double operator() (int i, int j) const
	{
		double d_element = 1.0*_size_row / (_size_row - 1.0); //n/(n-1)


		if (i == j) { d_element = 0.0; } //zero diagonal

										 //exceptions
		if (i<0 || i >= _size_row || j<0 || j >= _size_row)
		{
			d_element = 0.0;
		}

		return d_element;
	}
public:
	RepWeight_FHDI(int n) : _size_row(n) { }; //constructor
	~RepWeight_FHDI() { }; //desctructor

private:
	int _size_row;

private:
	RepWeight_FHDI(const RepWeight_FHDI &);
	const RepWeight_FHDI & operator = (const RepWeight_FHDI &);

};
