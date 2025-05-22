
//======================
//Compact matrix class declaration 
//for a system with banded nature which were
// NonUniformly Distributed along row direction
//
// original has [RU]local(r,n)  // where r= n/totalnodes and the last Proc may have the different one
//
// compact  has comp[RU]local(r, 2*Nband) where Nband = bandwidth
//
//
//Note: 1. only the elements within the bandwidth are stored in the compact storage
//      2. but the access can happen through the normal index (i,j) for RU_i,j 
//      3. all operation on an element outside the bandwidth is null
//
// Last update: July 11, 2010
//======================
class Compact_Matrix_forBanded_and_NonUniform_Distribution{

public:
	//==========================
	// a member fn to 
	//initialize all items in the matrix
	//==========================
	void clear(double value=0.0); 

	//======================
	//get the matrix dimension assuming rectangular shape of (_size_row, _size_col)
	//Note: _size_row and _size_col are private. 
	//      'const' after () means that this member fn is itself constant which is only meaningful for member fn's
	//      remark: only const object can call const member fn. 
	//      below member fn is necessary to access the private _size otherwise it's not possible to access to it from outside
	//======================
	int n_equation() const {return _n_equation;}
	int size_row() const {return _size_row;}
	int size_col() const {return _size_col;}


	//==========================
	// a member fn at the construction step to 
	//put items into the matrix
	//with global indices(i,j) or [RU]local 
	//to proper local indices of comp[RU]local
	//==========================
	void put_one_element(double value=0.0, const int i=0, const int j=0); 

	//==========================
	// a member fn at the construction step to 
	//put by addition of items into the matrix
	//with global indices(i,j) or [RU]local 
	//to proper local indices of comp[RU]local
	//==========================
	void putBYadd_one_element(double value=0.0, const int i=0, const int j=0); 

    //=========================
	//a member fn to 
	// find maxvalue 
	//=========================
    double max_value() ; 



    //======================
	//(i,j) operator overloading 
	//to access to the items of the array storage
	//Note: _block is one-dimensional vector for the compact matrix
	//      below operation overload of () is necessary to access the private _block otherwise it's not possible to access to it from outside
	//======================
	//probably, can be updated
	//Note: return is a reference to double 
	//      this () overload seems to be used in x=class(i,j) or cout<<class(i,j)
	//======================
	double & operator() (int i, int j) 
	{
		//std::cout<<"\n\n in &()========"; 

        double d_temp=0.0;
		double& d_element = d_temp ;  //default 

		//================
		//(i,j) is for [RU]local(i,j)
		//so proper indexing is required 
		//================
		int i_first = _first_non_zero[i] ; //get the stored index of first non zero element on ith row
        
		if(j< i_first || j >= i_first + _size_col) //outof bandwidth. _size_col= 2*Nband
		{
		    //NONE
		}
		else //within bandwidth
		{
			//d_element = _block[i* _size_col + (j - i_first)] ;  //this does not working for class(i,j) = 111 or so. 
			return _block[i*_size_col + (j - i_first)] ;   //this is working for class(i,j) = 111 or so.
		}
		
		return d_element;
	}       

    //========================
	//probably, only get the stored data
	//Note: return is a const double value 
	//========================
	double   operator() (int i, int j) const 
	{
		//std::cout<<"\n\n in () const========"; 
		
		double d_element=0.0; 

		//================
		//(i,j) is for [RU]local(i,j)
		//so proper indexing is required 
		//================
		int i_first = _first_non_zero[i] ; //get the stored index of first non zero element on ith row
        
		if(j< i_first || j >= i_first + _size_col) //outof bandwidth. _size_col= 2*Nband
		{
			d_element = 0.0;
		}
		else //within bandwidth
		{
			d_element = _block[i* _size_col + (j - i_first)] ;
		}	
		

		return d_element;
	} 


    //==========================
	//get the stored _block[_size_row*_size_col]
	//Feb 1, 2013
	//==========================
	void get_block(double* d_value); 

    //==========================
	//put the new into _block[_size_row*_size_col]
	//Feb 6, 2013
	//==========================
	void put_block(double* d_value); 

	//==========================
	//get the stored _first_non_zero[_size_row]
	//Feb 1, 2013
	//==========================
	void get_first_non_zero(int* i_value); 

	//==========================
	//put the new into _first_non_zero[_size_row]
	//Feb 6, 2013
	//==========================
	void put_first_non_zero(int* i_value); 

public:
	Compact_Matrix_forBanded_and_NonUniform_Distribution(int n_equation, int size_row, int size_col); //constructor
	~Compact_Matrix_forBanded_and_NonUniform_Distribution();        //destructor


//================
//data members
//================
private:
	const int _size_row; 
	const int _size_col; 
	const int _n_equation; //total size in one direction of the system
	
	double* _block;  
	int* _first_non_zero ; //array for INDEX of first non-zero term of each column of the original [RU]local(n,r) 
	
//================
//below is for avoiding possible error by an automatic task done by compiler
//In fact doing nothing as below makes it stable 
//================
private:
	Compact_Matrix_forBanded_and_NonUniform_Distribution(const Compact_Matrix_forBanded_and_NonUniform_Distribution &) ;
	const Compact_Matrix_forBanded_and_NonUniform_Distribution & operator = (const Compact_Matrix_forBanded_and_NonUniform_Distribution &) ;
};
