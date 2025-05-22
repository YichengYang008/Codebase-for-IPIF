#include "List_string_FHDI_MPI.h"
//======================
//Compact "String" LIST class declaration with the KNOWN ROW numbers
//to replace "R" list class
//
//
//Note: 1. only the given number of elements are stored in the compact storage
//      2. but the access can happen through the normal index (i,j)  like c++, starting from 0
//      3. all operation on an element outside the size of each list entity is null
//
// Last update: Nov 23, 2016
//
// by Dr. Cho, I. 
// All rights reservd
//======================

//================
//implementation of List_string_FHDI class
//================
void List_string_FHDI::initialize(int new_size_row)
{
	_size_row = new_size_row; 
	
	_n_each_row_size = NULL;
	_n_each_row_size = new int[new_size_row];
	for(int i=0; i<new_size_row; i++) _n_each_row_size[i] = 0 ; 
	
	_v_block.clear(); //return this to size 0
}

void List_string_FHDI::unlist(std::vector<std::string> & s_value)
//Description==================================
//	get all the stored non-null strings from _v_block 
//  to s_value[]
//  like R's "unlist()" 
//
//IN   : int  i_row    = target row in the list
//OUT  : std::vector<std::string> s_value[n_size]  where n_size must be known before calling this fn. 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	int n_size = _v_block.size();
		
	for(int i=0; i<n_size; i++) {s_value.push_back(_v_block[i]);} 

    return;
}

void List_string_FHDI::put_entire_block(std::vector<std::string> s_value)
//Description==================================
//	put the entire block into the storage _v_block 
//  from s_value[]
//
//IN   : std::vector<std::string> s_value[n_size]  where n_size must be known before calling this fn. 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	int n_size = (int)s_value.size();
		
	for(int i=0; i<n_size; i++) {_v_block.push_back(s_value[i]);} 

    return;
}

void List_string_FHDI::get_block(const int i_row, std::string s_value[])
//Description==================================
//	get stored block at the i_row of the list 
//
//IN   : int  i_row    = target row in the list
//OUT  : std::string s_value[n_size]  where n_size must be known before calling this fn. 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	int n_size = _n_each_row_size[i_row]; 
	
	//accumulated size of all the previous rows in the list
	int i_sum = 0; for(int k=0; k<i_row; k++) {i_sum += _n_each_row_size[k];}
	
	for(int i=0; i<n_size; i++) {s_value[i] = _v_block[i_sum + i];} 

    return;
}


void List_string_FHDI::put_block(const int i_row, const int n_size, std::string s_value[])
//Description==================================
//	put the new row into block
//  Note: 1. if current row i_row was not stored before, just append it by using push_back
//        2. if this row has been stored before, replacement takes place at the row
//
//IN  : int i_row  = target row number of the list (from 0 like c++ index)
//IN  : int n_size = ACTUAL size of the current row 
//IN  : std::string s_value[n_size] 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	int n_existing_size = _n_each_row_size[i_row]; 
	//----------------
	//first time input
	//----------------
	if(n_existing_size == 0)
	{
		//----------------
		//store the new data into the _v_block 
		//----------------
		for(int i=0; i<n_size; i++) _v_block.push_back(s_value[i]);
		
		//-----------------
		//update the size of current row of list
		//-----------------
		_n_each_row_size[i_row] = n_size; 
			
	}
	//---------------
	//replace existing stored data
	//---------------
	if(n_existing_size > 0)
	{
		//accumulated size of all the previous rows in the list
		int i_sum = 0; for(int k=0; k<i_row; k++) {i_sum += _n_each_row_size[k];}
		
		for(int i=0; i<n_size; i++) {_v_block[i_sum + i] = s_value[i];} 
	}

    return;
}

void List_string_FHDI::put_block(const int i_row, std::vector<std::string> s_value)
//Description==================================
//	put the new row store in vector into block
//  Note: 1. if current row i_row was not stored before, just append it by using push_back
//        2. if this row has been stored before, replacement takes place at the row
//
//IN  : int i_row  = target row number of the list (from 0 like c++ index)
//IN  : std::vector<std::string> s_value[n_size] 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	const int n_size = (int)s_value.size(); 
	
	int n_existing_size = _n_each_row_size[i_row]; 
	//----------------
	//first time input
	//----------------
	if(n_existing_size == 0)
	{
		//----------------
		//store the new data into the _v_block 
		//----------------
		for(int i=0; i<n_size; i++) _v_block.push_back(s_value[i]);
		
		//-----------------
		//update the size of current row of list
		//-----------------
		_n_each_row_size[i_row] = n_size; 
			
	}
	//---------------
	//replace existing stored data
	//---------------
	if(n_existing_size > 0)
	{
		//accumulated size of all the previous rows in the list
		int i_sum = 0; for(int k=0; k<i_row; k++) {i_sum += _n_each_row_size[k];}
		
		for(int i=0; i<n_size; i++) {_v_block[i_sum + i] = s_value[i];} 
	}

    return;
}


void List_string_FHDI::get_a_row_size(const int i_row, int & i_value)
//Description==================================
//	get stored _n_each_row_size of a row at i_row 
//
//IN   : int i_row 	= the row number of the list
//OUT  : int i_value 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	i_value = _n_each_row_size[i_row]; 

    return;
}


void List_string_FHDI::put_a_row_size(const int i_row, int i_value)
//Description==================================
//	put the new size of the list into _n_each_row_size
//  
//
//IN  : int i_row   = the row number of the list
//IN  : int i_value 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	_n_each_row_size[i_row] = i_value; 

    return;
}

//=====================
//print out List_string_FHDI
//=====================
void List_string_FHDI::print_List_string_FHDI()
{
	int n_row = (*this).size_row();
    std::string* s_temp;	
	for(int i=0; i<n_row; i++)
	{
		int i_temp = 0; (*this).get_a_row_size(i, i_temp);
		if(i_temp>0) //only for meaningful row 
		{
			
			s_temp = new std::string[i_temp];
			(*this).get_block(i, s_temp);
			RPrint(i);
			RPrint(s_temp, i_temp);
			
			delete[] s_temp; 
			
		}
	}
	return; 
}

//=====================
//print out ONE Row of List_string_FHDI
//=====================
void List_string_FHDI::print_one_List_string_FHDI(const int i_row)
{
	int n_row = (*this).size_row(); 
	std::string* s_temp; 
	if(i_row < n_row) 
	{
		int i = i_row; //target row number
		
		int i_temp = 0; (*this).get_a_row_size(i, i_temp);
		if(i_temp>0) //only for meaningful row 
		{
			s_temp = new std::string[i_temp];
			(*this).get_block(i, s_temp);
			RPrint(i);
			RPrint(s_temp, i_temp);
			
			delete[] s_temp; 
			 
		}
	}
	return; 
}
//=====================
//Constructor
//=====================
List_string_FHDI::List_string_FHDI(int size_row)
: _size_row(size_row), _n_each_row_size(new int[size_row])
{
	//=======
	//initialize with 0
	//=======
	for(int i=0; i<size_row; i++) _n_each_row_size[i] = 0 ; 

}

//=====================
//Destructor
//=====================
List_string_FHDI::~List_string_FHDI()
{
	delete[] _n_each_row_size ; 
}
