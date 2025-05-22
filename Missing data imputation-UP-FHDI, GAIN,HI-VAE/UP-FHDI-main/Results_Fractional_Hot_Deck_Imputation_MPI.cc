//#include "order_FHDI_Yicheng.cc"
void Results_Fractional_Hot_Deck_Imputation(const int i_size_ol, 
		const int ncol,  const int nrow,
		int* id_ol, double* w_ol, double** d_oy, double** d_ox, 
		rbind_FHDI &rbind_imat_FHDI, int** r,
		
		rbind_FHDI &rbind_ipmat_FHDI, 
		rbind_FHDI &rbind_Resp_FHDI, 
		rbind_FHDI &rbind_irmat_FHDI)
//Description----------------------
// prepare output results of  FHDI 
//
//ipmat  = final imputation results
//     	col1: ID 	= unit index
//		col2: FID 	= ID of fractionally imputed value
// 		col3: WGT 	= weight 
//		col4: FWGT	= Frational weight
//		col5: Variables 
//		col6: Responses
//irmat  = imputation results related to the categorized matrix 
//     	col1: ID 	= unit index
//		col2: FID 	= ID of fractionally imputed value
//		col3: OID	= original rank of the imputed value
//		col4: ORDER = SN(selected donor)
//		col5: FEFIW	= Fefi weights 
//		col6: CELL	= cells //
//
// original R code: Dr. Im, J. and Dr. Kim, J. 
// c++ code: 		Dr. Cho, I. 
// All rights reserved
// 
// updated: Nov 21, 2016
//
//IN   : int i_size_ol = rows with observed cells 
//IN   : int ncol = number of column 
//IN   : int* id_ol [i_size_ol]; //same as oid
//IN   : double* w_ol [i_size_ol]; //same as ow
//IN   : double** d_oy (i_size_ol, ncol);  //observed oritianl matrix 
//IN   : double** d_ox (i_size_ol, ncol);  //observed categorized matrix 

//IN   : rbind_imat_FHDI (column =7+2*ncol) = accumulated result matrix of FHDI 
//IN   : int    r(nrow, ncol) = matrix of missing (0)/observed (1) index  
//OUT  : rbind_FHDI  rbind_ipmat_FHDI(4+ncol); //column size is 4+ncol
//OUT  : rbind_FHDI  rbind_Resp_FHDI(ncol+1);  //separate response matrix. Note: in R version it is 
//                                               attached to ipmat  
//OUT  : rbind_FHDI  rbind_irmat_FHDI(5+ncol); //column size is 5+ncol
//----------------------
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;
	double Result_1 = MPI_Wtime();
	double** d_ipmat0 = New_dMatrix(i_size_ol, 4+ncol);  
	for(int i=0; i<i_size_ol; i++)
	{
		//col1: oid
		d_ipmat0[i][0] = id_ol[i]; 
		//col2: 1
		d_ipmat0[i][1] = 1;
		//col3: ow
		d_ipmat0[i][2] = w_ol[i]; 		
		//col4: 1
		d_ipmat0[i][3] = 1;
		//col5 set: oy
		for(int j=0; j<ncol; j++) d_ipmat0[i][4+j] = d_oy[i][j]; 				
	}
	//testout
	//RPrint("d_ipmat0 :");
	//RPrint(d_ipmat0, i_size_ol, 4+ncol); 

	
	//-----
	//FHDI results output first NOTE: FHDI results should be considered separately!!!
	//-----
	const int i_row_imat_FHDI = rbind_imat_FHDI.size_row(); //total rows of accumulated matrix 
	double** d_ipmat1_FHDI = New_dMatrix(i_row_imat_FHDI, 4+ncol); 
	for(int i=0; i<4+ncol; i++)
	{
		for(int j=0; j< i_row_imat_FHDI; j++)
			d_ipmat1_FHDI[j][i] = rbind_imat_FHDI(j, i);
	}
	//testout
	//RPrint("d_ipmat1_FHDI :");
	//RPrint(d_ipmat1_FHDI, i_row_imat_FHDI, 4+ncol); 

	
	//---
	//final return matrix of impat of FHDI
	//---
	//rbind_FHDI  rbind_ipmat_FHDI(4+ncol); //column size is 4+ncol //defined outside 
	rbind_ipmat_FHDI.bind_blocks(i_size_ol, 4+ncol, d_ipmat0); //new addition with ipmat0
	rbind_ipmat_FHDI.bind_blocks(i_row_imat_FHDI, 4+ncol, d_ipmat1_FHDI); //append ipmat1 of FHDI
	//if (mynode == 0) cout << "Result1 FHDI not FEFI Running time is " << MPI_Wtime() - Result_1 << endl;
	//testout
	//RPrint("rbind_ipmat_FHDI after binding :");
	//rbind_ipmat_FHDI.print_rbind_FHDI(); 
	
	//---
	//re-order ipmat of FHDI with respect to oid, the first column
	//---
	double Result2 = MPI_Wtime();

	const int i_row_ipmat_FHDI = (int)rbind_ipmat_FHDI.size_row(); 
	int* i_order_ipmat_FHDI = new int[i_row_ipmat_FHDI];
	for(int i=0; i<i_row_ipmat_FHDI; i++) i_order_ipmat_FHDI[i] = rbind_ipmat_FHDI(i,0);

	if (i_row_ipmat_FHDI < totalnodes) {
		cout<<"Error!!! The size of i_order_ipmat_FHDI is too small for parallel"<<endl;
		return;
	}
	order_FHDI_Yicheng(i_order_ipmat_FHDI, i_row_ipmat_FHDI);
	//if (mynode == 0) cout << "Result2_1 time is " << MPI_Wtime() - Result2 << endl;

	//backup before re-order
	double Result2_2 = MPI_Wtime();
	double** ipmat_FHDI_backup = New_dMatrix(i_row_ipmat_FHDI, 4+ncol); //ordered ipmat FHDI
	for(int i=0; i<i_row_ipmat_FHDI; i++)
	{
		int i_row_current = i_order_ipmat_FHDI[i]-1; //-1 is for ACTUAL location
		for(int j=0; j<4+ncol; j++) ipmat_FHDI_backup[i][j] = rbind_ipmat_FHDI(i_row_current, j);
	}
	//if (mynode == 0) cout << "Result2_2 time is " << MPI_Wtime() - Result2_2 << endl;
	//----
	//re-initialize and store the ordered matrix 
	//----
	double Result2_3 = MPI_Wtime();
	rbind_ipmat_FHDI.initialize(4+ncol); //Note: not imat but ipmat!!!
	rbind_ipmat_FHDI.bind_blocks(i_row_ipmat_FHDI, 4+ncol, ipmat_FHDI_backup);	
	//if (mynode == 0) cout << "Result2_3 time is " << MPI_Wtime() - Result2_3 << endl;
	//testout
	//RPrint("======== after FHDI ============");
	//RPrint("ipmat_FHDI: "); rbind_ipmat_FHDI.print_rbind_FHDI();
	
	//----
	//make table of id (1st column) of ipmat
	//----
	double Result2_4 = MPI_Wtime();
	double* d_first_column_ipmat = new double[i_row_ipmat_FHDI];
	std::vector<double> v_table_name_1stcol_ipmat; 
	std::vector<int>    v_table_count_1stcol_ipmat; 
	
	for(int i=0; i<i_row_ipmat_FHDI; i++) d_first_column_ipmat[i] = rbind_ipmat_FHDI(i,0);
	
	table_cpp_Yicheng(d_first_column_ipmat, i_row_ipmat_FHDI, 
		      v_table_name_1stcol_ipmat, v_table_count_1stcol_ipmat);
	//if (mynode == 0) cout << "Result2_4 time is " << MPI_Wtime() - Result2_4 << endl;
	//const int i_size_table_ipmat = (int)v_table_count_1stcol_ipmat.size(); 
	//testout
	//RPrint("table name of ipmat: "); RPrint(v_table_name_1stcol_ipmat);
	//RPrint("table name of ipmat: "); RPrint(v_table_count_1stcol_ipmat);
	//if (mynode == 0) cout << "Result2 Running time is " << MPI_Wtime() - Result2 << endl;
	//----
	//make Resp matrix 
	//----
	double  Result3 = MPI_Wtime();
	double** d_Resp_FHDI = New_dMatrix(i_row_ipmat_FHDI, ncol+1); 
	for(int i=0; i<ncol; i++)
	{
		int i_temp = 0; 
		for(int k=0; k<nrow; k++) 
		{
			int i_repeat = v_table_count_1stcol_ipmat[k]; 
			for(int j=0; j<i_repeat; j++)	
			{
				d_Resp_FHDI[i_temp][i] = r[k][i]; 
				i_temp++; 
			}
		}
	}
	//-----
	//last column of Resp is prod of response: 0=at least one missing 
	//-----
	for(int i=0; i<i_row_ipmat_FHDI; i++)
	{
		int i_temp =1; 
		for(int j=0; j<ncol; j++) i_temp = i_temp * d_Resp_FHDI[i][j]; 
		//0= at least one missing; 1=all observed 
		d_Resp_FHDI[i][ncol] = i_temp; //last column at ncol+1 
	}
	rbind_Resp_FHDI.bind_blocks(i_row_ipmat_FHDI, ncol+1, d_Resp_FHDI);	//prep return 
	//if (mynode == 0) cout << "Result3 Running time is " << MPI_Wtime() - Result3 << endl;
	//---------
	//Note: in R serial version, ipmat and Resp are attached column-wise, but 
	//in c++ version, it is kept separate for better memory usage
	//Nov 21, 2016
	//---------
	//testout
	//RPrint("Resp_FHDI: "); RPrint(d_Resp_FHDI, i_row_ipmat_FHDI, ncol+1);
	
	

	//-------------------------
	//-------------------------
	//make irmat of FHDI
	//-------------------------
	//-------------------------
	double  Result4 = MPI_Wtime();
	const int nci_FHDI = 7+2*ncol; 

	double** d_irmat0 = New_dMatrix(i_size_ol, 5+ncol);  
	for(int i=0; i<i_size_ol; i++)
	{
		//col1: oid
		d_irmat0[i][0] = id_ol[i]; 
		//col2: 1
		d_irmat0[i][1] = 1;
		//col3: oid
		d_irmat0[i][2] = id_ol[i]; 		
		//col4: 1
		d_irmat0[i][3] = 1;
		//col5: 1
		d_irmat0[i][4] = 1;
		//col6 set: ox
		for(int j=0; j<ncol; j++) d_irmat0[i][5+j] = d_ox[i][j]; 				
	}
	//testout
	//RPrint("d_irmat0 :");
	//RPrint(d_irmat0, i_size_ol, 5+ncol); 
			
	//----
	//extract columns 1, 2, (nci-2):nci, -(1:ncol of ipmat0, (nci-2):nci)
	//----
	int* i_col_irmat1 = new int[5+ncol]; //columns to be extracted from imat_FHDI
	i_col_irmat1[0] = 1; //ACTUAL col id
	i_col_irmat1[1] = 2; //ACTUAL col id
	i_col_irmat1[2] = nci_FHDI-2; //ACTUAL col id
	i_col_irmat1[3] = nci_FHDI-1; //ACTUAL col id
	i_col_irmat1[4] = nci_FHDI  ; //ACTUAL col id
    //exclude 1:4+ncol, i.e., columns of ipmat0
	for(int i=0; i<ncol; i++) i_col_irmat1[5+i] = i+(4+ncol+1); //ACTUAL id 
    
	
	double** d_irmat1_FHDI = New_dMatrix(i_row_imat_FHDI, 5+ncol); 
	for(int i=0; i<(5+ncol); i++)
	{
		int i_temp_irmat1 = i_col_irmat1[i] - 1 ; //-1 actual column id
		for(int j=0; j< i_row_imat_FHDI; j++)
			d_irmat1_FHDI[j][i] = rbind_imat_FHDI(j, i_temp_irmat1);
	}
	//testout
	//RPrint("d_irmat1_FHDI :");
	//RPrint(d_irmat1_FHDI, i_row_imat_FHDI, 5+ncol); 
	
	
	//---
	//final return matrix of irmat of FHDI
	//---
	//rbind_FHDI  rbind_irmat_FHDI(5+ncol); //column size is 5+ncol //defined outside 
	rbind_irmat_FHDI.bind_blocks(i_size_ol, 5+ncol, d_irmat0); //new addition with irmat0
	rbind_irmat_FHDI.bind_blocks(i_row_imat_FHDI, 5+ncol, d_irmat1_FHDI); //append irmat1 of FHDI
	//if (mynode == 0) cout << "Result4 Running time is " << MPI_Wtime() - Result4 << endl;
	//testout
	//RPrint("rbind_irmat_FHDI after binding :");
	//rbind_irmat_FHDI.print_rbind_FHDI(); 
	
	
	//---
	//re-order irmat of FHDI with respect to oid, the first column
	//---
	double  Result5 = MPI_Wtime();
	const int i_row_irmat_FHDI = (int)rbind_irmat_FHDI.size_row(); 
	int* i_order_irmat_FHDI = new int[i_row_irmat_FHDI];
	for(int i=0; i<i_row_irmat_FHDI; i++) i_order_irmat_FHDI[i] = rbind_irmat_FHDI(i,0);

	if (i_row_irmat_FHDI < totalnodes) {
		cout<<"Error!!!! The size of i_order_irmat_FHDI is too small for parallel"<<endl;
		return;
	}
	order_FHDI_Yicheng(i_order_irmat_FHDI, i_row_irmat_FHDI);
	//if (mynode == 0) cout << "Result5 Running time is " << MPI_Wtime() - Result5 << endl;
	double  Result6 = MPI_Wtime();
	//backup before re-order
	double** irmat_FHDI_backup = New_dMatrix(i_row_irmat_FHDI, 5+ncol); //ordered irmat FHDI
	for(int i=0; i<i_row_irmat_FHDI; i++)
	{
		int i_row_current = i_order_irmat_FHDI[i]-1; //-1 is for ACTUAL location
		for(int j=0; j<5+ncol; j++) irmat_FHDI_backup[i][j] = rbind_irmat_FHDI(i_row_current, j);
	}
	//----
	//re-initialize and store the ordered matrix 
	//----
	rbind_irmat_FHDI.initialize(5+ncol); //Note: not imat but irmat!!!
	rbind_irmat_FHDI.bind_blocks(i_row_irmat_FHDI, 5+ncol, irmat_FHDI_backup);	
	//if (mynode == 0) cout << "Result6 Running time is " << MPI_Wtime() - Result6 << endl;
	//testout
	//RPrint("======== after FHDI ============");
	//RPrint("irmat_FHDI: "); rbind_irmat_FHDI.print_rbind_FHDI();
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//-------------------
	//Deallocation 
	//-------------------
	Del_dMatrix(d_ipmat0, i_size_ol, 4+ncol); //
	Del_dMatrix(d_ipmat1_FHDI, rbind_imat_FHDI.size_row(), 4+ncol); //
	delete[] i_order_ipmat_FHDI;//
	Del_dMatrix(ipmat_FHDI_backup, i_row_ipmat_FHDI, 4+ncol);//

	delete[] d_first_column_ipmat; //
	Del_dMatrix(d_Resp_FHDI, i_row_ipmat_FHDI, ncol+1);//
	
	Del_dMatrix(d_irmat0, i_size_ol, 5+ncol); 	//
	delete[] i_col_irmat1;//

	Del_dMatrix(d_irmat1_FHDI, rbind_imat_FHDI.size_row(), 5+ncol);//
	delete[] i_order_irmat_FHDI; //
	Del_dMatrix(irmat_FHDI_backup, i_row_irmat_FHDI, 5+ncol); //
	
	return; 	
}