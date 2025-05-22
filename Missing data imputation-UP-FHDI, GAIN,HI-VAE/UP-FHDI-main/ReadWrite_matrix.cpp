//
//  ReadWrite_matrix.cpp
//
#include <stdio.h>

//double type// read & write
void ReadWrite_matrix(double** matrix,
	const int nrow, const int ncol,
	ostream&  Out_matrix_binary,
	ifstream&  In_matrix_binary,

	const int i_row_target,
	double* d_read_out,
	int i_option_ReadWrite)
	//Description==============================================
	//  Read and Write "matrix[][]" (double values)
	//
	//IN   : i_option_ReadWrite = 1
	//  when writing, sequentially append the matrix[nrow][ncol] on the file
	//       i_option_ReadWrite = 2
	//  when reading, entire column of the "i_row_target" row are read
	//       return the restored row in the output vector
	//OUT  : double* d_read_out[ncol]
	//
	//before call this function in main, the following must be declared
	//for writing
	//ofstream Out_matrix_binary ; //for Complete restart of MPI_Group, for backup
	//  Out_matrix_binary.open(Text_Test_Restart_binary, ios::out|ios::binary|ios::trunc); //Write general data   //by 'trunc' discard any current contents
	//Out_matrix_binary.close();
	//
	//for reading
	//ifstream In_matrix_binary ;
	//In_matrix_binary.open(Text_Test_Restart_binary, ios::in|ios::binary); //read general data
	//In_matrix_binary.seekg (0, ios::beg); //set the get pointer to the beginning of the file
	//In_matrix_binary.close();
	//
	//==========================================================
{
	double* d_temp = NULL;
	int double_total = 0;
	int d_loc = 0;
	int block = 0;

	//----------
	//WRITE (appending the matrix onto the file)
	//----------
	if (i_option_ReadWrite == 1)
	{
		double_total = nrow*ncol; //number of total double variables
		d_temp = new double[double_total];
		d_loc = 0;//reset
		for (int i = 0; i<nrow; i++) { for (int j = 0; j<ncol; j++) d_temp[d_loc++] = matrix[i][j]; }  //On All proc's
																									   //Write---
		block = sizeof(double)*double_total;
		Out_matrix_binary.write((char*)d_temp, block);
		delete[] d_temp;
	}

	//----------
	//READ (appending the matrix onto the file)
	//----------
	if (i_option_ReadWrite == 2)
	{
		//error check
		if (i_row_target >= nrow || i_row_target < 0)
		{
			cout << "Error! Attempt to read Binary file out of range!" << endl; return;
		}

		double_total = ncol; //number of total double variables (one row is read)
		d_temp = new double[double_total];

		//move the pointer to the desired row location
		int i_new_position = sizeof(double)*i_row_target*ncol; //i_row_target is the id of the desired row [0,nrow-1]
		In_matrix_binary.seekg(i_new_position); //set the get pointer to the beginning of the file

												//Read---
		block = sizeof(double)*double_total;
		In_matrix_binary.read((char*)d_temp, block);
		for (int j = 0; j<ncol; j++) d_read_out[j] = d_temp[j];  //only one row is passed
		delete[] d_temp;
	}

	return;
}



//double// write only
void ReadWrite_matrix(double** matrix,
	const int nrow, const int ncol,
	ostream&  Out_matrix_binary)
	//Description==============================================
	//  Write "matrix[nrow][ncol]" (double values)
	//  onto Binary file Out_matrix_binary
	//
	//  Append matrix, if multiple time call
	//==========================================================
{
	double* d_temp = NULL;
	int double_total = 0;
	int d_loc = 0;
	int block = 0;

	//----------
	//WRITE (appending the matrix onto the file)
	//----------
	double_total = nrow*ncol; //number of total double variables
	d_temp = new double[double_total];
	d_loc = 0;//reset
	for (int i = 0; i<nrow; i++) { for (int j = 0; j<ncol; j++) d_temp[d_loc++] = matrix[i][j]; }  //On All proc's
																								   //Write---
	block = sizeof(double)*double_total;
	Out_matrix_binary.write((char*)d_temp, block);
	delete[] d_temp;

	return;
}


//double type// read only
void ReadWrite_matrix(const int nrow, const int ncol,

	ifstream&  In_matrix_binary,

	const int i_row_target,
	double* d_read_out)
	//Description==============================================
	//  Read one row of "matrix[nrow][ncol]" (double values)
	//  from the Binary file "In_matrix_binary"
	//  return in the vector "d_read_out[ncol]"
	//
	//OUT  : double* d_read_out[ncol]
	//
	//
	//==========================================================
{
	double* d_temp = NULL;
	int double_total = 0;
	int d_loc = 0;
	int block = 0;

	//----------
	//READ one row only
	//----------
	double_total = ncol; //number of total double variables (one row is read)
	d_temp = new double[double_total];

	//move the pointer to the desired row location
	int i_new_position = sizeof(double)*i_row_target*ncol; //i_row_target is the id of the desired row [0,nrow-1]
	In_matrix_binary.seekg(i_new_position); //set the get pointer to the beginning of the file

											//Read---
											//error check
	if (i_row_target >= nrow || i_row_target < 0)
	{
		cout << "Error! Attempt to read Binary file out of range!" << endl; return;
	}

	block = sizeof(double)*double_total;
	In_matrix_binary.read((char*)d_temp, block);
	for (int j = 0; j<ncol; j++) d_read_out[j] = d_temp[j];  //only one row is passed
	delete[] d_temp;

	return;
}

//double type// read only with specified column
void ReadWrite_matrix_column(const int nrow, const int ncol,
	ifstream&  In_matrix_binary,
	const int i_col_target,
	double* d_read_out)
	//Description==============================================
	//  Read one row of "matrix[nrow][ncol]" (double values)
	//  from the Binary file "In_matrix_binary"
	//  return in the vector "d_read_out[ncol]"
	//
	//OUT  : double* d_read_out[ncol]
	//
	//Note that the source matrix should be written column-wise in the binary file
	//==========================================================
{
	double* d_temp = NULL;
	int block = 0;
	d_temp = new double[1];

	for (int j = 0; j < nrow; j++) {

		int i_new_position = sizeof(double)*ncol*j + sizeof(double)*i_col_target; //i_row_target is the id of the desired row [0,nrow-1]

		In_matrix_binary.seekg(i_new_position); //set the get pointer to the beginning of the file
		block = sizeof(double);
		In_matrix_binary.read((char*)d_temp, block);
		d_read_out[j] = d_temp[0];
	}

	delete[] d_temp;

	return;
}



//double type// read only: this function is made to read complete colum-wise distributed z matrix (after quantile process) only

void ReadWrite_matrix_datz(const int nrow, const int ncol,
	ifstream&  In_matrix_binary, int totalnodes,
	double** matrix)
	//Description==============================================
	//  Read one row of "matrix[nrow][ncol]" (double values)
	//  from the Binary file "In_matrix_binary"
	//  return in the vector "d_read_out[ncol]"
	//
	//OUT  : double* d_read_out[ncol]
	//
	//Note that the source matrix should be written column-wise in the binary file
	//==========================================================
{

	const int L = ncol; //size of d_rw 
	int numWorkPerProc = (int)floor(1.0*L / (1.0*totalnodes - 1));
	int numWorkLocalLast = L - numWorkPerProc * (totalnodes - 2);


	double* d_temp = NULL;
	double* d_temp_last = NULL;

	int block = 0;
	d_temp = new double[numWorkPerProc];
	d_temp_last = new double[numWorkLocalLast];

	for (int i = 0; i < nrow; i++) {
		int counter = 0;

		for (int j = 0; j < (totalnodes - 1); j++) {

			int i_new_position = 0;

			if (j != (totalnodes - 2)) i_new_position = sizeof(double)*numWorkPerProc*i + sizeof(double)*nrow*numWorkPerProc*j;
			if (j == (totalnodes - 2)) i_new_position = sizeof(double)*nrow*numWorkPerProc*(totalnodes - 2) + i*numWorkLocalLast * sizeof(double);
			//if (i != 0 && j == (totalnodes - 2)) i_new_position = sizeof(double)*nrow*numWorkPerProc*(totalnodes - 2) + (i+1)* numWorkLocalLast*sizeof(double);

			//cout << "i_new_position_complete is " << i_new_position << " at i=" << i << ", j=" << j << endl;

			In_matrix_binary.seekg(i_new_position); //set the get pointer to the beginning of the file

			if (j != (totalnodes - 2)) {
				block = sizeof(double)*numWorkPerProc;
				In_matrix_binary.read((char*)d_temp, block);

				for (int k = 0; k < numWorkPerProc; k++) {
					matrix[i][counter] = d_temp[k];
					counter++;
					//cout << "counter at i=" << i << ", j=" << j << " is " << counter << endl;
				}


			}

			if (j == (totalnodes - 2)) {
				block = sizeof(double)*numWorkLocalLast;
				In_matrix_binary.read((char*)d_temp_last, block);

				for (int k = 0; k < numWorkLocalLast; k++) {
					matrix[i][counter] = d_temp_last[k];
					counter++;
					//cout << "counter at i=" << i << ", j=" << j << " is " << counter << endl;
				}


			}

		}

		//for (int m = 0; m < ncol;m++) cout << setw(20) << matrix[i][m];
		//cout << endl;

	}//end of main loop

	 //cout << "Enf of read data" << endl;

	delete[] d_temp;
	delete[] d_temp_last;

	return;
}

//double type// read only: this function is made to read ith row of colum-wise distributed z matrix (after quantile process) only
void ReadWrite_matrix_datz(const int nrow, const int ncol, const int i_row_target,
	ifstream&  In_matrix_binary, int totalnodes,
	double* d_read_out)
	//Description==============================================
	//  Read one row of "matrix[nrow][ncol]" (double values)
	//  from the Binary file "In_matrix_binary"
	//  return in the vector "d_read_out[ncol]"
	//
	//OUT  : double* d_read_out[ncol]
	//
	//Note that the source matrix should be written column-wise in the binary file
	//==========================================================
{
	cout << "totalnodes inside read is " << totalnodes << endl;
	const int L = ncol; //size of d_rw 
	int numWorkPerProc = (int)floor(1.0*L / (1.0*totalnodes - 1));
	int numWorkLocalLast = L - numWorkPerProc * (totalnodes - 2);


	double* d_temp = NULL;
	double* d_temp_last = NULL;

	int block = 0;
	d_temp = new double[numWorkPerProc];
	d_temp_last = new double[numWorkLocalLast];

	int counter = 0;

	for (int j = 0; j < (totalnodes - 1); j++) {

		int i_new_position = 0;

		if (j != (totalnodes - 2)) i_new_position = sizeof(double)*numWorkPerProc*i_row_target + sizeof(double)*nrow*numWorkPerProc*j;
		if (j == (totalnodes - 2)) i_new_position = sizeof(double)*nrow*numWorkPerProc*(totalnodes - 2) + i_row_target*numWorkLocalLast * sizeof(double);
		//if (i != 0 && j == (totalnodes - 2)) i_new_position = sizeof(double)*nrow*numWorkPerProc*(totalnodes - 2) + (i+1)* numWorkLocalLast*sizeof(double);

		cout << "i_new_position_row is " << i_new_position << " at i_row_target=" << i_row_target << ", j=" << j << endl;

		In_matrix_binary.seekg(i_new_position); //set the get pointer to the beginning of the file

		if (j != (totalnodes - 2)) {
			block = sizeof(double)*numWorkPerProc;
			In_matrix_binary.read((char*)d_temp, block);

			for (int k = 0; k < numWorkPerProc; k++) {
				d_read_out[counter] = d_temp[k];
				counter++;
				//cout << "counter at i_row_target=" << i_row_target << ", j=" << j << " is " << counter << endl;
			}


		}

		if (j == (totalnodes - 2)) {
			block = sizeof(double)*numWorkLocalLast;
			In_matrix_binary.read((char*)d_temp_last, block);

			for (int k = 0; k < numWorkLocalLast; k++) {
				d_read_out[counter] = d_temp_last[k];
				counter++;
				//cout << "counter at i_row_target=" << i_row_target << ", j=" << j << " is " << counter << endl;
			}


		}

	}


	//cout << "Enf of read data" << endl;

	delete[] d_temp;
	delete[] d_temp_last;

	return;
}


//double type// read only
void Read_matrix_binary(double** m_out, const int n_row, const int n_col,
	ifstream& ReadIn_binary)
	//Description-----------------------------
	//    print out real double matrix in binary format
	//    row-by-row
	//    non-block writing
	//----------------------------------------
	//IN   : double** m_out[n_row][n_col]  = to be read in
	//----------------------------------------
{
	int double_total = 0;

	double_total = n_col * n_row; // number of total elements

	double*   d_temp = new double[double_total];

	int block_byte = sizeof(double)*(double_total); //no of bytes of double type

	ReadIn_binary.read((char*)d_temp, block_byte);

	int d_loc = 0;//reset

	for (int i = 0; i < n_row; i++)
	{
		for (int k = 0; k < n_col; k++)
		{
			m_out[i][k] = d_temp[d_loc]; //each row
			d_loc++;
		}
	}

	//deallocate --------
	delete[] d_temp;

	return;

}

//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
//integer// read&write
void ReadWrite_matrix(int** matrix,
	const int nrow, const int ncol,
	ostream&  Out_matrix_binary,
	ifstream&  In_matrix_binary,

	const int i_row_target,
	int* i_read_out,
	int i_option_ReadWrite)
	//Description==============================================
	//  Read and Write "matrix[][]" (integer values)
	//
	//IN   : i_option_ReadWrite = 1
	//  when writing, sequentially append the matrix[nrow][ncol] on the file
	//       i_option_ReadWrite = 2
	//  when reading, entire column of the "i_row_target" row are read
	//       return the restored row in the output vector
	//OUT  : int* i_read_out[ncol]
	//
	//before call this function in main, the following must be declared
	//ofstream Out_matrix_binary ; //for Complete restart of MPI_Group, for backup
	//  Out_matrix_binary.open(Text_Test_Restart_binary, ios::out|ios::binary|ios::trunc); //Write general data   //by 'trunc' discard any current contents
	//Out_matrix_binary.close();
	//==========================================================
{
	int* i_temp = NULL;
	int integer_total = 0;
	int i_loc = 0;
	int block = 0;

	//----------
	//WRITE (appending the matrix onto the file)
	//----------
	if (i_option_ReadWrite == 1)
	{
		integer_total = nrow*ncol; //number of total int variables
		i_temp = new int[integer_total];
		i_loc = 0;//reset
		for (int i = 0; i<nrow; i++) { for (int j = 0; j<ncol; j++) i_temp[i_loc++] = matrix[i][j]; }  //On All proc's
																									   //Write---
		block = sizeof(int)*integer_total;
		Out_matrix_binary.write((char*)i_temp, block);
		delete[] i_temp;
	}

	//----------
	//READ (appending the matrix onto the file)
	//----------
	if (i_option_ReadWrite == 2)
	{
		//error check
		if (i_row_target >= nrow || i_row_target < 0)
		{
			cout << "Error! Attempt to read Binary file out of range!" << endl; return;
		}

		integer_total = ncol; //number of total int variables (one row is read)
		i_temp = new int[integer_total];

		//move the pointer to the desired row location
		int i_new_position = sizeof(int)*i_row_target*ncol; //i_row_target is the id of the desired row [0,nrow-1]
		In_matrix_binary.seekg(i_new_position); //set the get pointer to the beginning of the file

												//Read---
		block = sizeof(int)*integer_total;
		In_matrix_binary.read((char*)i_temp, block);
		for (int j = 0; j<ncol; j++) i_read_out[j] = i_temp[j];  //only one row is passed
		delete[] i_temp;
	}


	return;
}

//integer//write only
void ReadWrite_matrix(int** matrix,
	const int nrow, const int ncol,
	ostream&  Out_matrix_binary)
	//Description==============================================
	//  Write "matrix[nrow][ncol]" (integer values)
	//  onto the Binary file "Out_matrix_binary"
	//
	// Append the matrix, if multiple time call
	//==========================================================
{
	int* i_temp = NULL;
	int integer_total = 0;
	int i_loc = 0;
	int block = 0;

	//----------
	//WRITE (appending the matrix onto the file)
	//----------
	integer_total = nrow*ncol; //number of total int variables
	i_temp = new int[integer_total];
	i_loc = 0;//reset
	for (int i = 0; i<nrow; i++) { for (int j = 0; j<ncol; j++) i_temp[i_loc++] = matrix[i][j]; }  //On All proc's
																								   //Write---
	block = sizeof(int)*integer_total;
	Out_matrix_binary.write((char*)i_temp, block);
	delete[] i_temp;

	return;
}

//integer// read
void ReadWrite_matrix(const int nrow, const int ncol,
	ifstream&  In_matrix_binary,

	const int i_row_target,
	int* i_read_out)
	//Description==============================================
	//  Read from "matrix[nrow][ncol]" (integer values)
	//  one row at "i_row_target"
	//
	//OUT  : int* i_read_out[ncol]
	//
	//==========================================================
{
	int* i_temp = NULL;
	int integer_total = 0;
	int i_loc = 0;
	int block = 0;

	//----------
	//READ
	//----------
	integer_total = ncol; //number of total int variables (one row is read)
	i_temp = new int[integer_total];

	//move the pointer to the desired row location
	int i_new_position = sizeof(int)*i_row_target*ncol; //i_row_target is the id of the desired row [0,nrow-1]
	In_matrix_binary.seekg(i_new_position); //set the get pointer to the beginning of the file

											//Read---
											//error check
	if (i_row_target >= nrow || i_row_target < 0)
	{
		cout << "Error! Attempt to read Binary file out of range!" << endl; return;
	}

	block = sizeof(int)*integer_total;
	In_matrix_binary.read((char*)i_temp, block);
	for (int j = 0; j<ncol; j++) i_read_out[j] = i_temp[j];  //only one row is passed
	delete[] i_temp;

	return;
}


//integer type// read only with specified column
void ReadWrite_matrix_column(const int nrow, const int ncol,
	ifstream&  In_matrix_binary,
	const int i_col_target,
	int* d_read_out)
	//Description==============================================
	//  Read one row of "matrix[nrow][ncol]" (double values)
	//  from the Binary file "In_matrix_binary"
	//  return in the vector "d_read_out[ncol]"
	//
	//OUT  : double* d_read_out[ncol]
	//
	//
	//==========================================================
{
	cout << "ReadWrite_matrix_column for integer" << endl;

	int* d_temp = NULL;
	int block = 0;
	d_temp = new int[1];

	for (int j = 0; j < nrow; j++) {

		int i_new_position = sizeof(int)*ncol*j + sizeof(int)*i_col_target; //i_row_target is the id of the desired row [0,nrow-1]

		In_matrix_binary.seekg(i_new_position); //set the get pointer to the beginning of the file
		block = sizeof(int);
		In_matrix_binary.read((char*)d_temp, block);
		d_read_out[j] = d_temp[0];
	}

	delete[] d_temp;

	return;
}