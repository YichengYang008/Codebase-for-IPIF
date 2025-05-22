//======================
//Compact LIST class declaration 
//to replace "R" list class
//
//
//Note: 1. only the given number of elements are stored in the compact storage
//      2. but the access can happen through the normal index (i,j)  like c++, starting from 0
//      3. all operation on an element outside the size of each list entity is null
//
// Last update: Oct 12, 2016
//
// by Dr. Cho, I. 
// All rights reservd
//======================
class List_FHDI{

public:

	int size_row()   const {return _size_row;}  //size of total rows
	int size_block() const {return _v_block.size();} //size of total meaningful data stored in the block

	//----------------
	//initialize all private memory
	//----------------
	void initialize(int new_size_row);
	
    //======================
	//(i,j) operator overloading 
	//to access to the items of the array storage
	//Note: _v_block is one-dimensional vector for the compact LIST
	//      below operation overload of () is necessary to access the private _v_block 
	//======================
	double & operator() (int i, int j) //ith list, jth entity (like c++,i.e. from 0)
	{
	
        double d_temp=0.0;
		double& d_element = d_temp ;  //default 

		//================
		//(i,j) is for LIST entity: ith list and jth term 
		//so proper indexing is required 
		//================
		int i_size_of_list = _n_each_row_size[i] ; //get the size of the ith list row
        
		if(j< 0 || j >= i_size_of_list ) //out of width. 
		{
		    //NONE
		}
		else //within width
		{
			//get accumulated location of (i-1)th list row
			int i_sum = 0; for(int k=0; k<i; k++) {i_sum += _n_each_row_size[k];}
			return _v_block[i_sum + j] ;
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
		//(i,j) is for LIST entity: ith list and jth term 
		//so proper indexing is required 
		//================
		int i_size_of_list = _n_each_row_size[i] ; //get the size of the ith list row
        
		if(j< 0 || j >= i_size_of_list ) //out of width. 
		{
			d_element = 0.0;
		}
		else //within width
		{
			int i_sum = 0; for(int k=0; k<i; k++) {i_sum += _n_each_row_size[k];}
			//return _block[i*_size_col + (j - i_size_of_list)] ;   //this is working for class(i,j) = 111 or so.
			return _v_block[i_sum + j] ;
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
	//get the stored _v_block at the i_row row of the list 
	//==========================
	void get_block(const int i_row, double* d_value); 
	void get_block_yicheng(const int i_row, std::vector<int> &v_value);
	void get_block_yicheng(const int i_row, std::vector<double> &v_value);
	void get_block(const int i_row, const int n_size_row, const int n_size_col, 
                   double** d_value);	

    //==========================
	//put the new into _v_block's row i_row with n_size entities 
	//==========================
	void put_block(const int i_row, const int n_size, double* d_value); 

	void put_block_yicheng(const int i_row, const int n_size, std::vector<double> v_value);

	void put_block_yicheng_replace(const int i_row, const int n_size, std::vector<double> v_value);

	void remove_block_yicheng(const int i_row); // remove existing row at i_row with n_size entities 

    void put_block(const int i_row, const int n_size_row, const int n_size_col, 
                   double ** d_value);	

    void put_block(const int i_row, std::vector<double> v_value);
	
	void put_block(const int i_row, std::vector<int> v_value);
	//==========================
	//get the stored _n_each_row_size at row i_row
	//==========================
	void get_a_row_size(const int i_row, int &i_value); 

	//==========================
	//put the new row size into storage
	//==========================
	void put_a_row_size(const int i_row, int i_value); 
	
	//=====================
	//print out List_FHDI
	//=====================
	void print_List_FHDI();	
	void print_List_FHDI_yicheng(ofstream& TestOut);
	void print_List_FHDI_yicheng();
	//=====================
	//print out ONE row of List_FHDI
	//=====================
	void print_one_List_FHDI(const int i_row);	
	void print_one_List_FHDI_yicheng(const int i_row, ofstream& TestOut);

public:
	List_FHDI(int size_row); //constructor
	~List_FHDI();        //destructor


//================
//data members
//================
private:
	int _size_row; //total row of the current LIST
	std::vector<double> _v_block;
	int* _n_each_row_size ; //array for the size of each row of LIST
	
//================
//below is for avoiding possible error by an automatic task done by compiler
//In fact doing nothing as below makes it stable 
//================
private:
	List_FHDI(const List_FHDI &) ;
	const List_FHDI & operator = (const List_FHDI &) ;
};
