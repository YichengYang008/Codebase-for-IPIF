
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
class rbind_FHDI{

public:

	int size_col()   const {return _size_col;}  //size of total columns
	int size_row()   const {return (int)_v_block.size()/_size_col;} //total number of rows
	int size_block() const {return (int)_v_block.size();} //size of total meaningful data stored in the block

	//----------------
	//initialize all private memory
	//----------------
	void initialize(int new_size_col);
	
    //======================
	//(i,j) operator overloading 
	//to access to the items of the array storage
	//Note: _v_block is one-dimensional vector for the compact matrix of columns ncol
	//      below operation overload of () is necessary to access the private _v_block 
	//======================
	double & operator() (int i, int j) //ith row, jth col (like c++,i.e. from 0)
	{
	
        double d_temp=0.0;
		double& d_element = d_temp ;  //default 

		//================
		//(i,j) is for an entity: ith row and jth column term 
		//so proper indexing is required 
		//================
		const int i_size_block = size_block(); //total stored values
        if(_size_col*i + 1 > i_size_block ) {return d_element;} //out of total range
		
		if(j< 0 || j >= _size_col ) //out of width. 
		{
		    //NONE
		}
		else //within width
		{
			return _v_block[i*_size_col + j] ;
		}
		
		return d_element;
	}       

    //========================
	//probably, only get the stored data
	//Note: return is a const double value 
	//========================
	double   operator() (int i, int j) const //like c++ rule, i.e., from 0 
	{
	
		double d_element=0.0; 

		//================
		//(i,j) is for an entity: ith row and jth col term 
		//so proper indexing is required 
		//================
		const int i_size_block = size_block(); //total stored values
		if(_size_col*i + 1 > i_size_block ) {return d_element;} //out of total range
		
		if(j< 0 || j >= _size_col ) //out of width. 
		{
		    //NONE
		}
		else //within width
		{
			return _v_block[i*_size_col + j] ;
		}		
		
		return d_element;

	} 

	//==========================
	//	get all the stored non-null values from _v_block 
	//==========================
	void unlist(std::vector<double> & d_value);

	//==========================
	//put entire block into the storage _v_block
	//==========================
	void put_entire_block(std::vector<double> d_value);
	
	//==========================
	//append a row into storage. Column size is fixed 
	//==========================
	void append_block(double* d_value); 
	
    //==========================
	//get the stored _v_block at the i_row row of the matrix 
	//==========================
	void get_block(const int i_row, double* d_value); 

    //==========================
	//append the new matrix onto _v_block's end 
	//MUST have the same column size as _size_col
	//==========================
	void bind_blocks(const int n_row, const int n_col, double ** d_value);


    //==========================
	//return a matrix of the stored entire matrix from _v_block 
	//MUST have the same column size and row size as the stored
	//==========================
	void matrix_rbind(const int n_row, const int n_col, double ** d_value);
	
	//=====================
	//print out rbind_FHDI
	//=====================
	void print_rbind_FHDI();
	void print_rbind_FHDI_Yicheng(ofstream& TestOut);

	//void TestOut_rbind_FHDI(ofstream& TestOut)//Yicheng
public:
	rbind_FHDI(int size_col); //constructor
	~rbind_FHDI();        //destructor


//================
//data members
//================
private:
	int _size_col; //fixed number of columns of the matrix
	std::vector<double> _v_block; //store many rows * _size_col data 
	
//================
//below is for avoiding possible error by an automatic task done by compiler
//In fact doing nothing as below makes it stable 
//================
private:
	rbind_FHDI(const rbind_FHDI &) ;
	const rbind_FHDI & operator = (const rbind_FHDI &) ;
};

