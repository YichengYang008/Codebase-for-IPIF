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
class List_string_FHDI{

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
	std::string & operator() (int i, int j) //ith list, jth entity (like c++,i.e. from 0)
	{
	
        std::string s_temp="";
		std::string& s_element = s_temp ;  //default 

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
		
		return s_element;
	}       

    //========================
	//probably, only get the stored data
	//Note: return is a const double value 
	//========================
	std::string   operator() (int i, int j) const //like c++ rule, i.e., from 0 
	{
	
		std::string s_element=""; 

		//================
		//(i,j) is for LIST entity: ith list and jth term 
		//so proper indexing is required 
		//================
		int i_size_of_list = _n_each_row_size[i] ; //get the size of the ith list row
        
		if(j< 0 || j >= i_size_of_list ) //out of width. 
		{
			s_element = "";
		}
		else //within width
		{
			int i_sum = 0; for(int k=0; k<i; k++) {i_sum += _n_each_row_size[k];}
			//return _block[i*_size_col + (j - i_size_of_list)] ;   //this is working for class(i,j) = 111 or so.
			return _v_block[i_sum + j] ;
		}	
		

		return s_element;
	} 

	//==========================
	//	get all the stored non-null values from _v_block 
	//==========================
	void unlist(std::vector<std::string> & s_value);

	//==========================
	//put entire block into the storage _v_block
	//==========================
	void put_entire_block(std::vector<std::string> s_value);
	
    //==========================
	//get the stored _v_block at the i_row row of the list 
	//==========================
	void get_block(const int i_row, std::string s_value[]); 

    //==========================
	//put the new into _v_block's row i_row with n_size entities 
	//==========================
	void put_block(const int i_row, const int n_size, std::string s_value[]); 
    void put_block(const int i_row, std::vector<std::string> s_value);
	
	//==========================
	//get the stored _n_each_row_size at row i_row
	//==========================
	void get_a_row_size(const int i_row, int &i_value); 

	//==========================
	//put the new row size into storage
	//==========================
	void put_a_row_size(const int i_row, int i_value); 
	
	//=====================
	//print out List_string_FHDI
	//=====================
	void print_List_string_FHDI();	
	
	//=====================
	//print out ONE row of List_string_FHDI
	//=====================
	void print_one_List_string_FHDI(const int i_row);	
	
public:
	List_string_FHDI(int size_row); //constructor
	~List_string_FHDI();        //destructor


//================
//data members
//================
private:
	int _size_row; //total row of the current LIST
	std::vector<std::string> _v_block;
	int* _n_each_row_size ; //array for the size of each row of LIST
	
//================
//below is for avoiding possible error by an automatic task done by compiler
//In fact doing nothing as below makes it stable 
//================
private:
	List_string_FHDI(const List_string_FHDI &) ;
	const List_string_FHDI & operator = (const List_string_FHDI &) ;
};
