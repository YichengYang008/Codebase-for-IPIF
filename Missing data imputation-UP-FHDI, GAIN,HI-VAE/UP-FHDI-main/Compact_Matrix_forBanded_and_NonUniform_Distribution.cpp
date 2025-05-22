//================
//implementation of Compact_Matrix_forBanded_and_NonUniform_Distribution class
//================

void Compact_Matrix_forBanded_and_NonUniform_Distribution::clear(double value)
//Description==================================
//   a member fn  to initialize matrix with value (default =0.0)
//IN   : double value
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	for(int i=0; i< _size_row*_size_col; i++) _block[i] =value ;

    return;
}


void Compact_Matrix_forBanded_and_NonUniform_Distribution::put_one_element(double value, const int i, const int j)
//Description==================================
//	put items into the matrix at the very first construction step 
//	with global indices(i,j) or [RU]local 
//	to proper local indices of comp[RU]local
//
//   RU_local(i,j)  --> comp_RU_local(i,j)
//   (r, n)             (r, 2*Nband)
//
//IN   : double value
//IN   : int i, j for indices of [RU]local(i       ,    j   )   
//                       ranging from     [0,r-1] and [0,n-1] respectively
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	const double tolerance = 1.e-15; 
 

	//std::cout<<"i,j ="<<i<<'\t'<<j<<" value="<<value<<std::endl;
	//std::cout<<"_first_non_zero ="<<_first_non_zero[j]<<std::endl;

	//=================
	//put the value into the proper location
	//=================
    if(fabs(value) > tolerance ) //only for non zero value
	{
		//=======
		//new first non zero position
		//Shfting of _block is necessary 
		//=======
		if(j < _first_non_zero[i]) //Note _first.. is set with default with n_equation-1 in constructor
		{
			//std::cout<<"new first and shifting begin"<<std::endl;
			//std::cout<<"\n before shifting _block ===="<<std::endl;
			//for(int k=0; k<_size_row; k++) std::cout<<_block[j*_size_row + k]<<'\t';

			int i_first_prev = _first_non_zero[i] ;

			//======
			//new first non zero index
			//======
			_first_non_zero[i] = j ;

			int i_gap = i_first_prev - j ;

			//======
			//prep shifting
			//with _size_col = 2*Nband
			//======
			double* d_buffer = new double[_size_col] ;
			for(int k=0; k<_size_col; k++) d_buffer[k] = _block[i*_size_col + k] ;
			
			//=====
			//shifting with _size_col = 2*Nband
			//by pushing back
			//=====
			for(int k=0; k<_size_col; k++)
			{
				if(k<  i_gap) _block[i*_size_col + k] = 0.0 ; //cleaning
				if(k>= i_gap) _block[i*_size_col + k] = d_buffer[k - i_gap] ;  //copy the shifted values
			}
			delete[] d_buffer ; 


			//======
			//put new first one
			//======
			_block[i*_size_col + 0] = value ;

			//std::cout<<"\n after shifting _block ===="<<std::endl;
			//for(int k=0; k<_size_row; k++) std::cout<<_block[j*_size_row + k]<<'\t';
				
		
		}

		//=======
		//NOT a new first non zero position
		//so just put into proper location
		//=======
		int i_loc = j - _first_non_zero[i] ;

		if(i_loc >= _size_col) { //neglect one outside 2*Nband to avoid non zero value outside bandwidth, which is not realistic case
		}
		else{
			_block[i*_size_col + i_loc] = value ;
		}


	}

	//std::cout<<"\n at the end of putting ===="<<std::endl;
	//for(int k=0; k<_size_row; k++) std::cout<<_block[j*_size_row + k]<<'\t';
	//std::cout<<'\n';

    return;
}


void Compact_Matrix_forBanded_and_NonUniform_Distribution::putBYadd_one_element(double value, const int i, const int j)
//Description==================================
//	put by addition of items into the matrix at the very first construction step 
//	with global indices(i,j) or [RU]local 
//	to proper local indices of comp[RU]local
//
//   RU_local(i,j)  - by + -> comp_RU_local(i,j)
//   (r, n)             (r, 2*Nband)
//
//IN   : double value
//IN   : int i, j for indices of [RU]local(i       ,    j   )   
//                       ranging from     [0,r-1] and [0,n-1] respectively
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	const double tolerance = 1.e-15; 
 

	//std::cout<<"i,j ="<<i<<'\t'<<j<<" value="<<value<<std::endl;
	//std::cout<<"_first_non_zero ="<<_first_non_zero[j]<<std::endl;

	//=================
	//put  the value by + into the proper location
	//=================
    if(fabs(value) > tolerance ) //only for non zero value
	{
		//=======
		//new first non zero position
		//Shfting of _block is necessary 
		//Since new first item is detected, no addition is necessary and only put is carried out
		//=======
		if(j < _first_non_zero[i]) //Note _first.. is set with default with n_equation-1 in constructor
		{
			//std::cout<<"new first and shifting begin"<<std::endl;
			//std::cout<<"\n before shifting _block ===="<<std::endl;
			//for(int k=0; k<_size_row; k++) std::cout<<_block[j*_size_row + k]<<'\t';

			int i_first_prev = _first_non_zero[i] ;

			//======
			//new first non zero index
			//======
			_first_non_zero[i] = j ;

			int i_gap = i_first_prev - j ;

			//======
			//prep shifting
			//with _size_col = 2*Nband
			//======
			double* d_buffer = new double[_size_col] ;
			for(int k=0; k<_size_col; k++) {
				d_buffer[k] = _block[i*_size_col + k] ;
				_block[i*_size_col + k] = 0.0 ; //cleaning
			}
			
			//=====
			//shifting with _size_col = 2*Nband
			//by pushing back
			//=====
			for(int k=0; k<_size_col; k++)
			{
				if(k<  i_gap) _block[i*_size_col + k] = 0.0 ; //cleaning
				if(k>= i_gap) _block[i*_size_col + k] = d_buffer[k - i_gap] ;  //copy the shifted values
			}
			delete[] d_buffer ; 


			//======
			//put new first one
			//======
			_block[i*_size_col + 0] = value ;

			//std::cout<<"\n after shifting _block ===="<<std::endl;
			//for(int k=0; k<_size_row; k++) std::cout<<_block[j*_size_row + k]<<'\t';
				
		
		}
       
		//=======
		//NOT a new first non zero position
		//so just put by addition into proper location
		//=======
		else
		{
			int i_loc = j - _first_non_zero[i] ;

			if(i_loc >= _size_col) { //neglect one outside 2*Nband to avoid non zero value outside bandwidth, which is not realistic case
			}
			else{
				_block[i*_size_col + i_loc] = _block[i*_size_col + i_loc] + value ;
			}
		}


	}

	//std::cout<<"\n at the end of putting ===="<<std::endl;
	//for(int k=0; k<_size_row; k++) std::cout<<_block[j*_size_row + k]<<'\t';
	//std::cout<<'\n';

    return;
}


double Compact_Matrix_forBanded_and_NonUniform_Distribution::max_value()
//Description==================================
//   a member fn  to find maximum value of present matrix
//
//OUT   : double  maximum value as a return
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	double d_temp=0.0;

	for(int i=0; i< _size_row*_size_col; i++) 
		if(_block[i] > d_temp) d_temp = _block[i] ;

    return d_temp;
}


void Compact_Matrix_forBanded_and_NonUniform_Distribution::get_block(double * d_value)
//Description==================================
//	get stored block
//  
//  Feb 1, 2013
//
//OUT  : double d_value[_size_row*_size_col] 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	for(int i=0; i<_size_row*_size_col; i++)
		d_value[i] = _block[i]; 

    return;
}


void Compact_Matrix_forBanded_and_NonUniform_Distribution::put_block(double * d_value)
//Description==================================
//	put the new into block
//  
//  Feb 6, 2013
//
//IN  : double d_value[_size_row*_size_col] 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	for(int i=0; i<_size_row*_size_col; i++)
		_block[i] = d_value[i]; 

    return;
}


void Compact_Matrix_forBanded_and_NonUniform_Distribution::get_first_non_zero(int * i_value)
//Description==================================
//	get stored _first_non_zero
//  
//  Feb 1, 2013
//
//OUT  : int i_value[_size_row] 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	for(int i=0; i<_size_row; i++)
		i_value[i] = _first_non_zero[i]; 

    return;
}


void Compact_Matrix_forBanded_and_NonUniform_Distribution::put_first_non_zero(int * i_value)
//Description==================================
//	put the new into _first_non_zero
//  
//  Feb 6, 2013
//
//IN  : int i_value[_size_row] 
//
//Note: all variable preceded by '_' are private of Grid
//=============================================
{
	for(int i=0; i<_size_row; i++)
		_first_non_zero[i] = i_value[i]; 

    return;
}

//=====================
//Constructor
//=====================
//IN   : int size_col = 2*Nband = no of row of the compact matrix
//IN   : int size_row = r       = no of row of the compact matrix. the last proc's r may be different
//
//Note: Since _size_.. are const, direct copy, e.g., _size=size, causes compile error.  
//=====================
Compact_Matrix_forBanded_and_NonUniform_Distribution::Compact_Matrix_forBanded_and_NonUniform_Distribution(int n_equation, int size_row, int size_col)
: _n_equation(n_equation), _size_row(size_row), _size_col(size_col), _block(new double[size_row*size_col]), _first_non_zero(new int[size_row])
{
    //=======
	//initialize _block with 0
	//=======
	clear(0.0) ; 

	//=======
	//initialize with n_equation as starting no for local index, which will be used in local putting
	//=======
	for(int i=0; i<size_row; i++) _first_non_zero[i] = n_equation-1 ; //first non zero element of each col

}

//=====================
//Destructor
//=====================
Compact_Matrix_forBanded_and_NonUniform_Distribution::~Compact_Matrix_forBanded_and_NonUniform_Distribution()
{
	delete[] _block ; 
	delete[] _first_non_zero ; 
}
