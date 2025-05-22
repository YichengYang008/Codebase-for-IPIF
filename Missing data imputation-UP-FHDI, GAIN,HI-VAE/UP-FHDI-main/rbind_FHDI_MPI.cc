
#include "rbind_FHDI_MPI.h"
//======================
//Compact row-based binding of matrix  
//to replace "R" rbind() function
//
//
//Note: 1. only the given number of elements are stored in the compact storage
//      2. but the access can happen through the normal index (i,j)  like c++, starting from 0
//      3. all operation on an element outside the size of each list entity is null
//
// Last update: Oct 27, 2016
//
// by Dr. Cho, I. 
// All rights reservd
//======================
//================
//implementation of rbind_FHDI class
//================
void rbind_FHDI::initialize(int new_size_col)
{
	_size_col = new_size_col; 
	_v_block.clear(); //return this to size 0
}

void rbind_FHDI::unlist(std::vector<double> & d_value)
//Description==================================
//	get all the stored non-null values from _v_block 
//  to d_value[]
//  like R's "unlist()" 
//
//OUT  : double d_value[n_size]  where n_size must be known before calling this fn. 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	int n_size = (int)_v_block.size();
		
	for(int i=0; i<n_size; i++) {d_value.push_back(_v_block[i]);} 

    return;
}

void rbind_FHDI::put_entire_block(std::vector<double> d_value)
//Description==================================
//	put the entire block into the storage _v_block 
//  from d_value[]
//
//IN   : double d_value[n_size]  where n_size must be known before calling this fn. 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	int n_size = (int)d_value.size();
		
	for(int i=0; i<n_size; i++) {_v_block.push_back(d_value[i]);} 

    return;
}

void rbind_FHDI::append_block(double* d_value)
//Description==================================
//	append the a row into the storage _v_block 
//  from d_value[]
//
//IN   : double d_value[n_col]  where n_col is the size_col must be known before calling this fn. 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
		
	for(int i=0; i<_size_col; i++) {_v_block.push_back(d_value[i]);} 

    return;
}

void rbind_FHDI::get_block(const int i_row, double* d_value)
//Description==================================
//	get stored block at the i_row of the matrix 
//
//IN   : int  i_row    = target row in the matrix [0,...]
//OUT  : double d_value[n_size]  where n_size must be known before calling this fn. 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	
	for(int i=0; i< _size_col; i++) {d_value[i] = _v_block[i_row*_size_col + i];} 

    return;
}


void rbind_FHDI::bind_blocks(const int nrow, const int ncol, double ** d_value)
//Description==================================
//	append the new matrix into block
//
//IN  : double d_value[nrow, ncol=_size_col]
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	if(ncol != _size_col){ cout<<"Error! column does not match!"<<endl; return;}
	
    for(int i=0; i<nrow; i++)
	{
		for(int j=0; j<ncol; j++)
		{
			_v_block.push_back(d_value[i][j]);
		}
	}

    return;
}

void rbind_FHDI::matrix_rbind(const int nrow, const int ncol, double ** d_value)
//Description==================================
//	return the entire matrix from block
//
//OUT  : double d_value[nrow=stored rows, ncol=_size_col]
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	const int i_size_row = (*this).size_row(); //total number of rows
	const int i_size_col = (*this).size_col(); //total number of columns
	if(nrow != i_size_row){ cout<<"Error! total rows do not match!"<<endl; return;}
	if(ncol != i_size_col){ cout<<"Error! total columns do not match!"<<endl; return;}
	
    for(int i=0; i<nrow; i++)
	{
		for(int j=0; j<ncol; j++)
		{
			d_value[i][j] = _v_block[i*i_size_col + j];
		}
	}

    return;
}


//=====================
//print out rbind_FHDI
//=====================
void rbind_FHDI::print_rbind_FHDI()
{
	int n_row = (*this).size_row(); 
	int n_col = (*this).size_col();

	double* d_temp = new double[n_col];
	//RPrint("Inside print function-------\n");
	for(int i=0; i<n_row; i++)
	{
			(*this).get_block(i, d_temp);
			cout<<endl; 
			//RPrint(i);
			RPrint_Yicheng(d_temp, n_col);
	}
	
	delete[] d_temp;
	return; 
}

void rbind_FHDI::print_rbind_FHDI_Yicheng(ofstream& TestOut)
{
	int n_row = (*this).size_row();
	int n_col = (*this).size_col();

	double* d_temp = new double[n_col];
	RPrint("Inside print function-------\n");
	for (int i = 0; i<n_row; i++)
	{
		(*this).get_block(i, d_temp);
		cout << endl;
		//RPrint(i);
		RPrint_Yicheng_Output(d_temp, n_col, TestOut);
	}

	delete[] d_temp;
	return;
}

//void rbind_FHDI::TestOut_rbind_FHDI(ofstream TestOut)//Yicheng
//{
//	int n_row = (*this).size_row();
//	int n_col = (*this).size_col();
//
//	double* d_temp = new double[n_col];
//	//RPrint("Inside print function-------\n");
//	for (int i = 0; i<n_row; i++)
//	{
//		(*this).get_block(i, d_temp);
//		for (int j = 0; j < n_col; j++) {
//			TestOut << setw(20) << d_temp[j];
//		}
//	}
//	TestOut << endl;
//
//	delete[] d_temp;
//	return;
//}


//=====================
//Constructor
//=====================
rbind_FHDI::rbind_FHDI(int size_col)
: _size_col(size_col) 
{
//none 
}

//=====================
//Destructor
//=====================
rbind_FHDI::~rbind_FHDI()
{
	_v_block.clear(); 
}

