//
//  MPI_IO.cpp
//
#include <stdio.h>


void MPI_In_raw(const int nrow, const int i_col_target,
	MPI_File fh, double* d_read_out)
	//Description==============================================
	//  Read one column of raw data "matrix[nrow][ncol]" (double values)
	//  from the Binary file
	//  return in the vector "d_read_out[nrow]"
	//
	//OUT  : double* d_read_out[ncol]
	//
	//Note that the source matrix should be written column-wise in the binary file
	//==========================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//MPI_Offset disp;
	unsigned long long int disp;

	disp = (unsigned long long)i_col_target*nrow * sizeof(double);
	//cout<<"disp is "<< disp <<" at column "<< i_col_target <<" at node "<<mynode<<endl;

	int success = 0;

	success = MPI_File_read_at(fh, disp, d_read_out, nrow, MPI_DOUBLE, &status);
	if (success != MPI_SUCCESS) cout << "MPI I/O fail to read the file!" << endl;

	//for (int i = 0; i < nrow; i++) cout<<"d_read_out["<<i<<"]: "<< d_read_out[i]<<" at column " << i_col_target << " at node " << mynode <<endl;

	return;
}

void MPI_In_raw(const int nrow, const int i_col_target,
	MPI_File fh, int* d_read_out)
	//Description==============================================
	//  Read one column of raw data "matrix[nrow][ncol]" (integer values)
	//  from the Binary file
	//  return in the vector "d_read_out[nrow]"
	//
	//OUT  : int* d_read_out[ncol]
	//
	//Note that the source matrix should be written column-wise in the binary file
	//==========================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//MPI_Offset disp;
	unsigned long long int disp;

	disp = (unsigned long long)i_col_target*nrow * sizeof(int);
	//cout<<"disp is "<< disp <<" at column "<< i_col_target <<" at node "<<mynode<<endl;

	int success = 0;

	success = MPI_File_read_at(fh, disp, d_read_out, nrow, MPI_INT, &status);
	if (success != MPI_SUCCESS) cout << "MPI I/O fail to read the file!" << endl;

	//for (int i = 0; i < nrow; i++) cout<<"d_read_out["<<i<<"]: "<< d_read_out[i]<<" at column " << i_col_target << " at node " << mynode <<endl;

	return;
}


void MPI_In_daty_row(const int nrow, const int ncol, const int i_row_target,
	MPI_File fh, double* d_read_out)
	//Description==============================================
	//  Read one row of raw daty "matrix[nrow][ncol]" (double values)
	//  from the Binary file
	//  return in the vector "d_read_out[nrow]"
	//
	//OUT  : double* d_read_out[ncol]
	//
	//Note that the source matrix should be written column-wise in the binary file
	//==========================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	unsigned long long int disp;
	int success = 0;
	double* d_temp = new double[1];// one can only read a single value eveytime

	for (int i = 0; i < ncol; i++) {
		disp = (unsigned long long)i_row_target * sizeof(double) + (unsigned long long)i * nrow * sizeof(double);
		//cout<<"disp is "<< disp <<" at i= "<<i<<endl;
		success = MPI_File_read_at(fh, disp, d_temp, 1, MPI_DOUBLE, &status);
		d_read_out[i] = d_temp[0];
	}

	if (success != MPI_SUCCESS) cout << "MPI I/O fail to read the file!" << endl;

	//for (int i = 0; i < nrow; i++) cout<<"d_read_out["<<i<<"]: "<< d_read_out[i]<<" at column " << i_col_target << " at node " << mynode <<endl;

	delete[] d_temp;
	return;
}


void MPI_In_datz_column(const int nrow, const int ncol, const int i_col_target,
	MPI_File fh, double* d_read_out)
	//Description==============================================
	//  Read one column of datz "matrix[nrow][ncol]" (double values)
	//  from the Binary file
	//  return in the vector "d_read_out[nrow]"
	//
	//OUT  : double* d_read_out[ncol]
	//
	//Note that the source matrix should be written row-wise (entire row) in the binary file
	//==========================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	unsigned long long int disp;
	int success = 0;
	double* d_temp = new double[1];// one can only read a single value eveytime

	for (int i = 0; i < nrow; i++) {
		disp = (unsigned long long)i_col_target * sizeof(double) + i*ncol * sizeof(double);
		//cout<<"disp is "<< disp <<" at i= "<<i<<endl;
		success = MPI_File_read_at(fh, disp, d_temp, 1, MPI_DOUBLE, &status);
		d_read_out[i] = d_temp[0];
	}

	if (success != MPI_SUCCESS) cout << "MPI I/O fail to read the file!" << endl;

	//for (int i = 0; i < nrow; i++) cout<<"d_read_out["<<i<<"]: "<< d_read_out[i]<<" at column " << i_col_target << " at node " << mynode <<endl;

	delete[] d_temp;
	return;
}

void MPI_In_datz(const int nrow, const int ncol, const int i_row_target,
	MPI_File fh, double* d_read_out)
	//Description==============================================
	//  Read a row of of matrix (double values)
	//  from the target Binary file 
	//  return in the array "d_read_out[ncol]"
	//
	//OUT  : double* d_read_out[ncol]
	//
	//Note that this function is especially for reading a specified row 
	//of categorized matrix written row-wisely from multiple processes
	//==========================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	// Need to know how categorized matrix is distributed colum-wisely in categorization process in advance
	const int L = ncol; //size of d_rw 
	int numWorkPerProc = (int)floor(1.0*L / (1.0*totalnodes - 1));
	int numWorkLocalLast = L - numWorkPerProc * (totalnodes - 2);
	double* d_temp = NULL;// buffer for slave processors
	double* d_temp_last = NULL;// buffer for the last slave processors

	d_temp = new double[numWorkPerProc];
	d_temp_last = new double[numWorkLocalLast];

	//MPI_Offset disp;// number of bytes from the beginning
	unsigned long long int disp;

	int success = 0;

	//-----------------------------------------------------
	//Read a row of source matrix across multiple processors
	//-----------------------------------------------------
	int counter = 0;
	for (int j = 0; j < (totalnodes - 1); j++) { // exclude processor 0

		if (j != (totalnodes - 2)) disp = (unsigned long long)sizeof(double)*numWorkPerProc*i_row_target + (unsigned long long)sizeof(double)*nrow*numWorkPerProc*j; // offset of starting reading
		if (j == (totalnodes - 2)) disp = (unsigned long long)sizeof(double)*nrow*numWorkPerProc*(totalnodes - 2) + (unsigned long long)i_row_target*numWorkLocalLast * sizeof(double); // offset of starting reading

																																				//cout << "disp is " << disp << " at i_row_target=" << i_row_target << ", j=" << j << endl;

																																				// slave processors excluding the last one
		if (j != (totalnodes - 2)) {
			success = MPI_File_read_at(fh, disp, d_temp, numWorkPerProc, MPI_DOUBLE, &status);
			if (success != MPI_SUCCESS) cout << "MPI I/O fail to read the file!" << endl;

			for (int k = 0; k < numWorkPerProc; k++) {
				d_read_out[counter] = d_temp[k];
				counter++;
				//cout << "counter at i_row_target=" << i_row_target << ", j=" << j << " is " << counter << endl;
			}


		}

		// the last slave processor
		if (j == (totalnodes - 2)) {
			success = MPI_File_read_at(fh, disp, d_temp_last, numWorkLocalLast, MPI_DOUBLE, &status);
			if (success != MPI_SUCCESS) cout << "MPI I/O fail to read the file!" << endl;

			for (int k = 0; k < numWorkLocalLast; k++) {
				d_read_out[counter] = d_temp_last[k];
				counter++;
				//cout << "counter at i_row_target=" << i_row_target << ", j=" << j << " is " << counter << endl;
			}


		}

	}


	//for (int k = 0; k < ncol; k++) cout << "d_read_out[" << k << "]: " << d_read_out[k] << " at node " << mynode << endl;


	delete[] d_temp;
	delete[] d_temp_last;

	return;
}

void MPI_In_datz(const int nrow, const int ncol, int numWorkPerProc,
	MPI_File fh, double** d_matrix) {
	//Description==============================================
	//  Read distributed z matrix
	//  from the Binary file
	//  return in the vector "d_matrix[nrow][ncol]"
	//
	//OUT  : double** d_matrix[nrow][ncol]
	//
	//Note that the source matrix should be written row-wise in the binary file
	//==========================================================

	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	if (mynode != 0) {
		int size = 0;
		size = nrow*ncol;
		double*   d_temp = new double[size];// buffer
		for (int k = 0; k < size; k++) d_temp[k] = 0.0;

		unsigned long long int disp;

		disp = (unsigned long long)(mynode - 1) * numWorkPerProc*ncol * sizeof(double);
		//cout<<"ZMAT disp is "<< disp <<" nrow is "<< nrow<<" size is "<< size <<" at node "<<mynode<<endl;

		int success = 0;
		success = MPI_File_read_at(fh, disp, d_temp, size, MPI_DOUBLE, &status);
		if (success != MPI_SUCCESS) cout << "MPI I/O fail to read the file!" << endl;

		int counter = 0;
		for (int i = 0; i < nrow; i++) {
			for (int j = 0; j < ncol; j++) {
				d_matrix[i][j] = d_temp[counter];
				counter++;
			}
		}

		delete[] d_temp;
	}

}

void MPI_Out_datz(const int nrow, const int ncol, int numWorkPerProc,
	MPI_File fh, double** d_matrix) {
	//Description==============================================
	//  Write a matrix (double values)
	//  to the target Binary file from multiple processors concurrently
	//  IN:   double** d_matrix[nrow][ncol]
	//  IN:   numWorkPerProc is particularly read for the last processor
	//
	//Note that this function is especially for writing categorized matrix 
	//row-wisely to local storage from multple processors concurrently
	//==========================================================

	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	if (mynode != 0) {
		double*   d_temp = new double[nrow*ncol];// buffer

		int counter = 0;
		for (int i = 0; i < nrow; i++) {
			for (int j = 0; j < ncol; j++) {
				d_temp[counter] = d_matrix[i][j]; // copy elements in matrix to buffer row-wisely
				counter++;
			}
		}

		//MPI_Offset disp;// number of bytes from the beginning
		unsigned long long int disp;
		
		disp = 0;
		int success = 0;

		disp = (unsigned long long)(mynode - 1)*nrow*numWorkPerProc * sizeof(double);

		success = MPI_File_write_at(fh, disp, d_temp, nrow*ncol, MPI_DOUBLE, &status);
		if (success != MPI_SUCCESS) cout << "MPI I/O fail to write the file!" << endl;

		delete[] d_temp;
	}

	return;
}

void MPI_Out_datz(const int nrow, const int ncol, MPI_File fh, double** d_matrix) {
	//Description==============================================
	//  Write a matrix (double values)
	//  to the target Binary file from multiple processors concurrently
	//  IN:   double** d_matrix[nrow][ncol]
	//  IN:   numWorkPerProc is particularly read for the last processor
	//
	//Note that this function is especially for writing categorized matrix 
	//row-wisely to local storage from multple processors concurrently
	//==========================================================

	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	if (mynode == 0) {
		double*   d_temp = new double[nrow*ncol];// buffer

		int counter = 0;
		for (int i = 0; i < nrow; i++) {
			for (int j = 0; j < ncol; j++) {
				d_temp[counter] = d_matrix[i][j]; // copy elements in matrix to buffer row-wisely
				counter++;
			}
		}

		//MPI_Offset disp;// number of bytes from the beginning
		unsigned long long int disp;
		
		disp = 0;
		int success = 0;

		success = MPI_File_write_at(fh, disp, d_temp, nrow*ncol, MPI_DOUBLE, &status);
		if (success != MPI_SUCCESS) cout << "MPI I/O fail to write the file!" << endl;

		delete[] d_temp;
	}

	return;
}

void MPI_Out_datz_rowise(const int nrow, const int ncol, const int ncol_total, int numWorkPerProc,
	MPI_File fh, double** d_matrix) {
	//Description==============================================
	//  Write a matrix (double values): non-contiguous writing 
	//  to the target Binary file from multiple processors concurrently
	//  IN:   double** d_matrix[nrow][ncol]
	//  IN:   numWorkPerProc is particularly read for the last processor
	//
	//d_matrix is generated column-wise
	//Note that this function is especially for writing categorized matrix 
	//row-wisely (entire row) to local storage from multple processors concurrently
	//==========================================================

	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	if (mynode != 0) {
		double*   d_temp = new double[ncol];// buffer

		unsigned long long int disp;
		int success = 0;

		for (int i = 0; i < nrow; i++) {

			for (int j = 0; j < ncol; j++) {
				d_temp[j] = d_matrix[i][j];
			}

			disp = 0;
			disp = (unsigned long long)(mynode - 1)*numWorkPerProc * sizeof(double) + (unsigned long long)i * ncol_total * sizeof(double);
			//cout<<"disp is "<< disp <<" at i = "<<i<<" at node "<<mynode<<endl;


			success = MPI_File_write_at(fh, disp, d_temp, ncol, MPI_DOUBLE, &status);
			if (success != MPI_SUCCESS) cout << "MPI I/O fail to write the file!" << endl;

		}


		delete[] d_temp;
	}

	return;
}

void MPI_Out_uox_mox(const int nrow, const int ncol, std::vector<int> i_size,
	MPI_File fh, double** d_matrix) {
	//Description==============================================
	//  Write a matrix (double values)
	//  to the target Binary file from multiple processors concurrently
	//  IN:   double** d_matrix[nrow][ncol]
	//  IN:   i_size: number of rows of uox or mox on all slave processors, in size of (totalnodes-1)
	//
	//Note that this function is especially for writing uox or mox matrix 
	//row-wisely to local storage from multple processors concurrently
	//==========================================================

	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;


	if (mynode != 0) {

		double* d_temp = new double[nrow*ncol];

		int counter = 0;
		for (int i = 0; i < nrow; i++) {
			for (int j = 0; j < ncol; j++) {
				d_temp[counter] = d_matrix[i][j]; // copy elements in matrix to buffer row-wisely
				counter++;
			}
		}

		//MPI_Offset disp;// number of bytes from the beginning
		unsigned long long int disp;
		
		disp = 0;

		int success = 0;
		for (int k = 0; k < mynode; k++) {
			disp = disp + (unsigned long long)i_size[k] * ncol * sizeof(double);
		}

		//cout << "disp is " << disp << " and nrow is " << nrow << " at node " << mynode << endl;

		success = MPI_File_write_at(fh, disp, d_temp, nrow*ncol, MPI_DOUBLE, &status);
		if (success != MPI_SUCCESS) cout << "MPI I/O fail to write the file!" << endl;

		delete[] d_temp;
	} //end of slave 

	return;
}


void MPI_In_uox_mox(const int ncol, const int i_row_target,
	MPI_File fh, double* d_read_out)
	//Description==============================================
	//  Read one row of uox or mox "matrix[nrow][ncol]" (double values)
	//  from the Binary file
	//  return in the vector "d_read_out[ncol]"
	//
	//OUT  : double* d_read_out[ncol]
	//
	//Note that the source matrix should be written row-wise in the binary file
	//==========================================================
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//MPI_Offset disp;
	unsigned long long int disp;

	disp = (unsigned long long)i_row_target*ncol * sizeof(double);


	int success = 0;

	success = MPI_File_read_at(fh, disp, d_read_out, ncol, MPI_DOUBLE, &status);
	if (success != MPI_SUCCESS) {
		cout << "disp is " << disp << " at row " << i_row_target << " at node " << mynode << endl;
		//cout << "MPI I/O fail to read the file in MPI_In_uox_mox at i_row_target at " << i_row_target << endl; 
	}

	//for (int i = 0; i < nrow; i++) cout<<"d_read_out["<<i<<"]: "<< d_read_out[i]<<" at column " << i_col_target << " at node " << mynode <<endl;

	return;
}

void MPI_Out_sorted_daty_column(const int nrow, const int ncol, MPI_File fh, double** d_matrix) {
	//Description==============================================
	//  Write a matrix (double values)
	//  to the target Binary file from multiple processors concurrently
	//  IN:   double** d_matrix[nrow][ncol]
	//  IN:   numWorkPerProc is particularly read for the last processor
	//
	//Note that this function is especially for writing categorized matrix 
	//row-wisely to local storage from multple processors concurrently
	//==========================================================

	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	if (mynode == 0) {
		double*   d_temp = new double[nrow*ncol];// buffer

		int counter = 0;
		for (int j = 0; j < ncol; j++) {
			for (int i = 0; i < nrow; i++) {
				d_temp[counter] = d_matrix[i][j]; // copy elements in matrix to buffer colum-wisely
				counter++;
			}
		}

		//MPI_Offset disp;// number of bytes from the beginning
		unsigned long long int disp;

		disp = 0;
		int success = 0;

		success = MPI_File_write_at(fh, disp, d_temp, nrow*ncol, MPI_DOUBLE, &status);
		if (success != MPI_SUCCESS) cout << "MPI I/O fail to write the file!" << endl;

		delete[] d_temp;
	}

	return;
}

void MPI_Out_sorted_daty_row(const int nrow, const int ncol, MPI_File fh, double** d_matrix) {
	//Description==============================================
	//  Write a matrix (double values)
	//  to the target Binary file from multiple processors concurrently
	//  IN:   double** d_matrix[nrow][ncol]
	//  IN:   numWorkPerProc is particularly read for the last processor
	//
	//Note that this function is especially for writing categorized matrix 
	//row-wisely to local storage from multple processors concurrently
	//==========================================================

	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	if (mynode == 0) {
		double*   d_temp = new double[nrow*ncol];// buffer

		int counter = 0;
		for (int i = 0; i < nrow; i++) {
			for (int j = 0; j < ncol; j++) {
				d_temp[counter] = d_matrix[i][j]; // copy elements in matrix to buffer colum-wisely
				counter++;
			}
		}


		//MPI_Offset disp;// number of bytes from the beginning
		unsigned long long int disp;

		disp = 0;
		int success = 0;

		success = MPI_File_write_at(fh, disp, d_temp, nrow*ncol, MPI_DOUBLE, &status);
		if (success != MPI_SUCCESS) cout << "MPI I/O fail to write the file!" << endl;

		delete[] d_temp;
	}

	return;
}




void MPI_Out_sorted_datr(const int nrow, const int ncol, MPI_File fh, int** d_matrix) {
	//Description==============================================
	//  Write a matrix (double values)
	//  to the target Binary file from multiple processors concurrently
	//  IN:   double** d_matrix[nrow][ncol]
	//  IN:   numWorkPerProc is particularly read for the last processor
	//
	//Note that this function is especially for writing categorized matrix 
	//row-wisely to local storage from multple processors concurrently
	//==========================================================

	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	if (mynode == 0) {
		int*   i_temp = new int[nrow*ncol];// buffer

		int counter = 0;
		for (int j = 0; j < ncol; j++) {
			for (int i = 0; i < nrow; i++) {
				i_temp[counter] = d_matrix[i][j]; // copy elements in matrix to buffer colum-wisely
				counter++;
			}
		}

		//MPI_Offset disp;// number of bytes from the beginning
		unsigned long long int disp;

		disp = 0;
		int success = 0;

		success = MPI_File_write_at(fh, disp, i_temp, nrow*ncol, MPI_INT, &status);
		if (success != MPI_SUCCESS) cout << "MPI I/O fail to write the file!" << endl;

		delete[] i_temp;
	}

	return;
}