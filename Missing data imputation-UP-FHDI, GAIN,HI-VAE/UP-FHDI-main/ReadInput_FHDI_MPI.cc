void ReadInput_FHDI_MPI(ifstream* pFileInput, int &i_option_read_data, int &i_option_ultra,
               int &i_option_perform, int &nrow_x,
               int &ncol_x, int &i_option_imputation,
               int &i_option_variance,
               int &i_option_merge,
			   int &M,
			   int &i_user_defined_datz, int &i_option_collapsing, int &i_option_SIS_type, int &top, int &i_option_cellmake, int &i_option_var_type,
	           int &memory,
               ofstream& TestOut)
//Description==========================
//  read input file
//  and save informations into proper
//  arrays
//IN   : ifstream& pFileInput = pointer to the object of ifstream class
//                              which indicate input file name
//OUT  : int    i_option_perform = 1: do all tasks; 2: Cell Make only; 3: Cell Prob only; 4: do all using user-defined datz 
//OUT  : int    nrow_x = row size of the data matrix
//OUT  : int	ncol_x = column size of the data matrix 
//OUT  : int    i_option_imputation = 1: FEFI; 2: FHDI;
//OUT  : int    i_option_variance = 0: no variance est.; 1: Jackknife var. est.;
//OUT  : int    i_option_merge = 1: fixed merge; 2: random merge;
//OUT  : int    M = number of donors for FHDI
//OUT  : int    i_user_defined_datz = 0: automatically make datz; 1: read and use the user-defined datz
//OUT  : int    i_option_SIS_type = 1: SIS with intersection; 2: SIS with union; 3: SIS with global ranking
//
//=====================================
{
  char buffer[256]; //storage for one line
  char c_temp[80];

  string s_temp ;
  int n_position=0;

  for(;;) // main loop
  {
     // read one line into buffer using member fn. getline
     pFileInput->getline(buffer, 256);
     s_temp = buffer;


     //================================================
     //INPUT INFORMATION
     //================================================
     n_position=0;
     n_position = s_temp.find("INPUT INFORMATION",0);
     if(n_position !=-1) //note: -1 means no mateched input
     {
        //TestOut <<'\n' <<'\n'<< "INPUT INFORMATION" <<'\n' << endl;

        for(;;)
        {
           // next line
           pFileInput->getline(buffer, 256);
           s_temp = buffer;
           //TestOut <<"i, s_temp     =     "<< i<<s_temp << endl;

		   //========
		   //i_option_read_data
		   //========
		   n_position = 0;
		   n_position = s_temp.find("i_option_read_data", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   istrstream one_line(buffer); // construct an istringstream class object with buffer
									  // note: new version= istringstream
			   one_line >> c_temp;
			   one_line >> i_option_read_data;

			   //TestOut << setw(30) << c_temp <<setw(10) << " = " <<setw(10) << i_option_read_data <<  "    (0:input from input.txt; 1:seperate input from daty.txt, datr.txt, and datz.txt) "  << endl;
		   }

		   //========
		   //i_option_ultra
		   //========
		   n_position = 0;
		   n_position = s_temp.find("i_option_ultra", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   istrstream one_line(buffer); // construct an istringstream class object with buffer
											// note: new version= istringstream
			   one_line >> c_temp;
			   one_line >> i_option_ultra;

			   //TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << i_option_ultra << "    (0:P-FHDI for big-n or big-p data based on memory ; 1:P-FHDI for ultra data based on memory and local storage) " << endl;
		   }


           //========
           //i_option_perform
           //========
           n_position=0;
           n_position = s_temp.find("i_option_perform",0);
           if(n_position != -1) //note: -1 means no mateched input
           {
              istrstream one_line(buffer); // construct an istringstream class object with buffer
                                     // note: new version= istringstream
              one_line >> c_temp  ;
              one_line >> i_option_perform;

			 // TestOut << setw(30) << c_temp <<setw(10) << " = " <<setw(10) << i_option_perform <<  "    (1:do all; 2:CellMake; 3:CellProb; 4: do all with user's datz) "  << endl;
           }

           //========
           //nrow
           //========
           n_position=0;
           n_position = s_temp.find("nrow",0);
           if(n_position != -1) //note: -1 means no mateched input
           {
              istrstream one_line(buffer); // construct an istringstream class object with buffer
                                     // note: new version= istringstream
              one_line >> c_temp  ;
              one_line >> nrow_x;

             // TestOut <<setw(30)<< c_temp <<setw(10)<<" = "<<setw(10)<< nrow_x << endl;
           }

           //========
           //ncol
           //========
           n_position=0;
           n_position = s_temp.find("ncol",0);
           if(n_position != -1) //note: -1 means no mateched input
           {
              istrstream one_line(buffer); // construct an istringstream class object with buffer
                                     // note: new version= istringstream
              one_line >> c_temp  ;
              one_line >> ncol_x;

             // TestOut <<setw(30)<< c_temp <<setw(10)<<" = "<<setw(10)<< ncol_x << endl;
           }

           //========
           //i_option_imputation
           //========
           n_position=0;
           n_position = s_temp.find("i_option_imputation",0);
           if(n_position != -1) //note: -1 means no mateched input
           {
              istrstream one_line(buffer); // construct an istringstream class object with buffer
                                     // note: new version= istringstream
              one_line >> c_temp  ;
              one_line >> i_option_imputation;

			 // TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << i_option_imputation <<  "    (1: FEFI; 2: FHDI)" << endl;
           }

           //========
           //i_option_variance
           //========
           n_position=0;
           n_position = s_temp.find("i_option_variance",0);
           if(n_position != -1) //note: -1 means no mateched input
           {
              istrstream one_line(buffer); // construct an istringstream class object with buffer
                                     // note: new version= istringstream
              one_line >> c_temp  ;
              one_line >> i_option_variance;

             // TestOut <<setw(30)<< c_temp << setw(10) << " = " << setw(10) << i_option_variance  <<"    (0: no variance est.; 1: Jackknife var. est.) "<< endl;
           }


           //========
           //i_option_merge
           //========
           n_position=0;
           n_position = s_temp.find("i_option_merge",0);
           if(n_position != -1) //note: -1 means no mateched input
           {
              istrstream one_line(buffer); // construct an istringstream class object with buffer
                                     // note: new version= istringstream
              one_line >> c_temp  ;
              one_line >> i_option_merge;

			  //TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << i_option_merge <<  "    (1: fixed; 2: random merge in CellMake)" << endl;
           }

		   //================================================
		   //M number of donors for FHDI 
		   //================================================
		   n_position = 0;
		   n_position = s_temp.find("i_donor", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   istrstream one_line(buffer); // construct an istringstream class object with buffer
											// note: new version= istringstream
			   one_line >> c_temp;
			   one_line >> M;

			  //TestOut << setw(30) << "M number of donors for FHDI " << setw(10) << " = " << setw(10) << M << endl;

		   }

			   

		   
		   //========
		   //i_user_defined_datz
		   //========
		   n_position = 0;
		   n_position = s_temp.find("i_user_defined_datz", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   istrstream one_line(buffer); // construct an istringstream class object with buffer
											// note: new version= istringstream
			   one_line >> c_temp;
			   one_line >> i_user_defined_datz;

			  // TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << i_user_defined_datz <<  "    (0: automatic; 1: user-defined) " << endl;
		   }

		   //================================================
		   //number of selected variables
		   //================================================
		   n_position = 0;
		   n_position = s_temp.find("i_option_collapsing", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   istrstream one_line(buffer); // construct an istringstream class object with buffer
											// note: new version= istringstream
			   one_line >> c_temp;
			   one_line >> i_option_collapsing;

			   //TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << i_option_collapsing << "    (0: old cell collapsing; !0: user-defined for new cell collapsing) " << endl;

		   }

		   //================================================
           //number of selected variables
           //================================================
		   n_position = 0;
		   n_position = s_temp.find("i_option_SIS_type", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   istrstream one_line(buffer); // construct an istringstream class object with buffer
											// note: new version= istringstream
			   one_line >> c_temp;
			   one_line >> i_option_SIS_type;

			  // TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << i_option_SIS_type << "    (1: SIS with intersection; 2: SIS with union; 3: SIS with global ranking) " << endl;

		   }

		   //================================================
		   //top ranking matrix for cell make
		   //================================================
		   n_position = 0;
		   n_position = s_temp.find("top_correlation", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   istrstream one_line(buffer); // construct an istringstream class object with buffer
											// note: new version= istringstream
			   one_line >> c_temp;
			   one_line >> top;

			 //  TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << top << "    (top correlation: 100 <= top <= ncol-1)" << endl;

		   }

		   //================================================
		   //different methods for cell make
		   //================================================
		   n_position = 0;
		   n_position = s_temp.find("i_option_cellmake", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   istrstream one_line(buffer); // construct an istringstream class object with buffer
											// note: new version= istringstream
			   one_line >> c_temp;
			   one_line >> i_option_cellmake;

			  // TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << i_option_cellmake << "    (1: Cell make with merging ; 2: Cell make with k-nearest neighbor) " << endl;

		   }

		   //================================================
		   //different methods for variance estimation
		   //================================================
		   n_position = 0;
		   n_position = s_temp.find("i_option_var_type", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   istrstream one_line(buffer); // construct an istringstream class object with buffer
											// note: new version= istringstream
			   one_line >> c_temp;
			   one_line >> i_option_var_type;

			 //  TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << i_option_var_type << "    (1: Parallel Jackkinfe variance estimation ; 2: Parallel variance estimation using linerization " << endl;

		   }

		   //================================================
		   //memory capacity of a processor
		   //================================================
		   n_position = 0;
		   n_position = s_temp.find("memory", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   istrstream one_line(buffer); // construct an istringstream class object with buffer
											// note: new version= istringstream
			   one_line >> c_temp;
			   one_line >> memory;

			 //  TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << memory << " GB for a processor" << endl;

		   }

		   //================================================
		   //max occurance of an unique pattern in uox or mox
		   //================================================
		   //n_position = 0;
		   //n_position = s_temp.find("max_pattern", 0);
		   //if (n_position != -1) //note: -1 means no mateched input
		   //{
			  // istrstream one_line(buffer); // construct an istringstream class object with buffer
					//						// note: new version= istringstream
			  // one_line >> c_temp;
			  // one_line >> max_overlap_size;

			  // TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << max_overlap_size << " for a unique pattern" << endl;

		   //}

		   //================================================
		   //max number of donors of all mox
		   //================================================
		   //n_position = 0;
		   //n_position = s_temp.find("max_donor", 0);
		   //if (n_position != -1) //note: -1 means no mateched input
		   //{
			  // istrstream one_line(buffer); // construct an istringstream class object with buffer
					//						// note: new version= istringstream
			  // one_line >> c_temp;
			  // one_line >> max_donor_size;

			  // TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << max_donor_size << " donors for a recipient" << endl;

		   //}

           //========
           //end of INPUT INFORMATION
           //========
           n_position=0;
           n_position = s_temp.find("END INPUT INFORMATION",0);
           if(n_position != -1) //note: -1 means no mateched input
           {
              return; // end of this input information
           }

        }
     }


   } // end of main loop


  return;
}

void InputData(ifstream* pFileInput,
			const int i_option_perform, const int nrow_x,
			const int ncol_x, const int i_option_imputation,
			const int i_option_variance,
			const int i_option_merge,
			const int i_user_defined_datz,

			double **daty,
			int    **datr,
			double **datz, 
			double *k,
			double *d,
			
			int    *NonCollapsible_categorical,

            ofstream& TestOut)
//Description==========================
//  read input file
//  and save informations into proper
//  arrays
//IN   : ifstream& pFileInput = pointer to the object of ifstream class
//                              which indicate input file name
//OUT  : double **daty = New_dMatrix(nrow_x, ncol_x);  original data with missing units
//OUT  : int    **datr = New_iMatrix(nrow_x, ncol_x); 	indices of 0 for missing cell and 1 for the observed cell
//OUT  : double **datz = New_dMatrix(nrow_x, ncol_x);  user-defined datz. Vaild only when i_user_defined_datz = 1
//OUT  : double *k = new double[ncol_x];		category for each column 
//OUT  : double *d = new double[nrow_x];  		weight for each row 
//OUT  : int    *id = new int[nrow_x]; 			id of each row 
//=====================================
{
  char buffer[25600000]; //storage for one line

  string s_temp ;
  int n_position=0;

  for(;;) // main loop
  {
     // read one line into buffer using member fn. getline
     pFileInput->getline(buffer, 25600000);
     s_temp = buffer;

     //================================================
     //daty original data matrix 
     //================================================
     n_position=0;
     n_position = s_temp.find("daty",0);
     if(n_position !=-1) //note: -1 means no matched input
     {

        int n_data=ncol_x;  //size of col
        double* d_out = new double[n_data]; 

		for(int j=0; j<nrow_x; j++)
        {

           // read one line into buffer using member fn. getline
           pFileInput->getline(buffer, 25600000);

           parse_dString(buffer, d_out, n_data) ;

           //save coordinates of current node
           for(int i=0; i<n_data; i++ )  daty[j][i]= d_out[i];
        }
		delete[] d_out ;

        /*TestOut <<'\n' <<'\n'<< "============= [daty] ==============" <<'\n' << endl;
		TestOut << setw(10) << "row  col:";  for (int j = 0; j < ncol_x; j++) TestOut<< setw(10) << j+1;
		TestOut << endl; 
		for(int i=0; i<nrow_x ; i++)
        {
			TestOut <<setw(10) <<i+1;
			for (int j = 0; j < ncol_x; j++) TestOut <<setw(10)<< daty[i][j] ;
			TestOut << endl;
        }*/

     }


     //================================================
     //datr indices for observed and missing cells 
     //================================================
     n_position=0;
     n_position = s_temp.find("datr",0);
     if(n_position !=-1) //note: -1 means no matched input
     {

        int n_data=ncol_x;  //size of col
        int* i_out = new int[n_data];

        for(int j=0; j<nrow_x; j++)
        {
           // read one line into buffer using member fn. getline
           pFileInput->getline(buffer, 25600000);
           
		   //int n_data=4;  //size of col
           //int i_out[n_data];
           
		   parse_iString(buffer, i_out, n_data) ;

           //save coordinates of current node
           for(int i=0; i<n_data; i++ )  datr[j][i]= i_out[i];
        }
        delete[] i_out; 

		/*TestOut << '\n' << '\n' << "============= [datr] ==============" << '\n' << endl;
		TestOut << "row  indices 0=missing; 1=observed " << endl;

		TestOut << setw(10) << "row  col:";  for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j+1;
		TestOut << endl;
        for(int i=0; i<nrow_x ; i++)
        {
			TestOut <<setw(10)<< i + 1 ;
			for (int j = 0; j < ncol_x; j++) TestOut << setw(10)<< datr[i][j] ;
			TestOut << endl;
        }*/


     }

	 //================================================
	 //k category vector 
	 //================================================
	 n_position = 0;
	 n_position = s_temp.find("category", 0);
	 if (n_position != -1) //note: -1 means no matched input
	 {

		 int n_data = ncol_x;  //size of col
		 double* d_out = new double[n_data];

		 for (int j = 0; j<1; j++)
		 {

			 // read one line into buffer using member fn. getline
			 pFileInput->getline(buffer, 25600000);

			 parse_dString(buffer, d_out, n_data);

			 //save coordinates of current node
			 for (int i = 0; i<n_data; i++)  k[i] = d_out[i];
		 }
		 delete[] d_out;


		/* TestOut << '\n' << '\n' << "k category per column" << '\n' << endl;
		 TestOut << "column : ";
		 for (int j = 0; j < ncol_x; j++) TestOut << setw(10)<<j+1 ;
		 TestOut << endl;
		 TestOut << "k[]    : ";
		 for (int j = 0; j < ncol_x; j++) TestOut << setw(10)<<k[j] ;
		 TestOut << endl;*/
		 
	 }

	 //================================================
	 //NonCollapsible_categorical
	 //================================================
	 n_position = 0;
	 n_position = s_temp.find("NonCollapsible_categorical", 0);
	 if (n_position != -1) //note: -1 means no matched input
	 {

		 int n_data = ncol_x;  //size of col
		 double* d_out = new double[n_data];

		 for (int j = 0; j<1; j++)
		 {

			 // read one line into buffer using member fn. getline
			 pFileInput->getline(buffer, 25600000);

			 parse_dString(buffer, d_out, n_data);

			 //save coordinates of current node
			 for (int i = 0; i<n_data; i++)  NonCollapsible_categorical[i] = d_out[i];
		 }
		 delete[] d_out;


		 /*TestOut << '\n' << '\n' << "NonCollapsible_categorical" << '\n' << endl;
		 TestOut << "column : ";
		 for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j + 1;
		 TestOut << endl;
		 TestOut << "NonCollapsible_categorical[]    : ";
		 for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << NonCollapsible_categorical[j];
		 TestOut << endl;*/

	 }

	 //================================================
	 //d weight vector 
	 //================================================
	 n_position = 0;
	 n_position = s_temp.find("weight", 0);
	 if (n_position != -1) //note: -1 means no matched input
	 {

		 int n_data = 1;  //size of col
		 double* d_out = new double[n_data];

		 for (int j = 0; j<nrow_x; j++)
		 {

			 // read one line into buffer using member fn. getline
			 pFileInput->getline(buffer, 25600000);

			 parse_dString(buffer, d_out, n_data);

			 //save coordinates of current node
			 for (int i = 0; i<n_data; i++)  d[j] = d_out[i];
		 }
		 delete[] d_out;


		/* TestOut << '\n' << '\n' << "weight per rows" << '\n' << endl;
		 TestOut << "row         : ";
		 for (int j = 0; j < nrow_x; j++) TestOut << setw(10)<<j+1 ;
		 TestOut << endl;
		 TestOut << "weight[]    : ";
		 for (int j = 0; j < nrow_x; j++) TestOut << setw(10)<<d[j] ;
		 TestOut << endl;*/

	 }



	 //================================================
	 //id of each row 
	 //================================================
	 //n_position = 0;
	 //n_position = s_temp.find("index", 0);
	 //if (n_position != -1) //note: -1 means no matched input
	 //{

		// int n_data = 1;  
		// int* i_out = new int[n_data];

		// for (int j = 0; j<nrow_x; j++)
		// {
		//	 // read one line into buffer using member fn. getline
		//	 pFileInput->getline(buffer, 256000);

		//	 parse_iString(buffer, i_out, n_data);

		//	 //save coordinates of current node
		//	 for (int i = 0; i<n_data; i++)  id[j] = i_out[i];
		// }
		// delete[] i_out;

		// TestOut << '\n' << '\n' << "id number of each row" << '\n' << endl;
		// TestOut << "row     : ";
		// for (int j = 0; j < nrow_x; j++) TestOut << setw(10)<<j ;
		// TestOut << endl;
		// TestOut << "id[]    : ";
		// for (int j = 0; j < nrow_x; j++) TestOut << setw(10)<<id[j] ;
		// TestOut << endl;

	 //}

	 //--------------
	 //read user defined datz, if needed
	 //--------------
	 if (i_user_defined_datz == 1)
	 {
		 n_position = 0;
		 n_position = s_temp.find("datz", 0);
		 if (n_position != -1) //note: -1 means no matched input
		 {

			 int n_data = ncol_x;  //size of col
			 double* d_out = new double[n_data];

			 for (int j = 0; j<nrow_x; j++)
			 {

				 // read one line into buffer using member fn. getline
				 pFileInput->getline(buffer, 25600000);

				 parse_dString(buffer, d_out, n_data);

				 //save coordinates of current node
				 for (int i = 0; i<n_data; i++)  datz[j][i] = d_out[i];
			 }
			 delete[] d_out;

			/* TestOut << '\n' << '\n' << "============= user-defined [datz] ==============" << '\n' << endl;
			 TestOut << setw(10) << "row  col:";  for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j;
			 TestOut << endl;
			 for (int i = 0; i<nrow_x; i++)
			 {
				 TestOut << setw(10) << i + 1;
				 for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << datz[i][j];
				 TestOut << endl;
			 }*/
		 }
	 }

	 //========
	 //end of INPUT INFORMATION
	 //========
	 n_position = 0;
	 n_position = s_temp.find("END DATA INPUT", 0);
	 if (n_position != -1) //note: -1 means no mateched input
	 {
		 return; // end of this input information
	 }

   } // end of main loop

  return;


}

void InputData_Seperate(ifstream* pFileInput, ifstream* pFileInput_daty, ifstream* pFileInput_datr,
	const int i_option_perform, const int nrow_x,
	const int ncol_x, const int i_option_imputation,
	const int i_option_variance,
	const int i_option_merge,
	const int i_user_defined_datz,

	double **daty,
	int    **datr,
	double **datz,
	double *k,
	double *d,

	int    *NonCollapsible_categorical,

	ofstream& TestOut)
	//Description==========================
	//  read input file
	//  and save informations into proper
	//  arrays
	//IN   : ifstream& pFileInput = pointer to the object of ifstream class
	//                              which indicate input file name
	//OUT  : double **daty = New_dMatrix(nrow_x, ncol_x);  original data with missing units
	//OUT  : int    **datr = New_iMatrix(nrow_x, ncol_x); 	indices of 0 for missing cell and 1 for the observed cell
	//OUT  : double **datz = New_dMatrix(nrow_x, ncol_x);  user-defined datz. Vaild only when i_user_defined_datz = 1
	//OUT  : double *k = new double[ncol_x];		category for each column 
	//OUT  : double *d = new double[nrow_x];  		weight for each row 
	//OUT  : int    *id = new int[nrow_x]; 			id of each row 
	//=====================================
{
	char buffer_daty[25600000]; //storage for one line

	//string s_temp_daty;

	// read one line into buffer using member fn. getline
	pFileInput_daty->getline(buffer_daty, 25600000);
	//s_temp_daty = buffer_daty;

	//================================================
	//daty original data matrix 
	//================================================

	int n_data = ncol_x;  //size of col
	double* d_out = new double[n_data];

	for (int j = 0; j < nrow_x; j++)
	{

		// read one line into buffer using member fn. getline
		pFileInput_daty->getline(buffer_daty, 25600000);

		parse_dString(buffer_daty, d_out, n_data);

		//save coordinates of current node
		for (int i = 0; i < n_data; i++)  daty[j][i] = d_out[i];
	}
	delete[] d_out;

	//TestOut << '\n' << '\n' << "============= [daty] ==============" << '\n' << endl;
	//TestOut << setw(10) << "row  col:";  for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j + 1;
	//TestOut << endl;
	//for (int i = 0; i < nrow_x; i++)
	//{
	//	TestOut << setw(10) << i + 1;
	//	for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << daty[i][j];
	//	TestOut << endl;
	//}
		



	char buffer_datr[25600000]; //storage for one line

	// read one line into buffer using member fn. getline
	pFileInput_datr->getline(buffer_datr, 25600000);

	//================================================
	//datr original data matrix 
	//================================================

	double* d_out2 = new double[n_data];

	for (int j = 0; j < nrow_x; j++)
	{

		// read one line into buffer using member fn. getline
		pFileInput_datr->getline(buffer_datr, 25600000);

		parse_dString(buffer_datr, d_out2, n_data);

		//save coordinates of current node
		for (int i = 0; i < n_data; i++)  datr[j][i] = d_out2[i];
	}
	delete[] d_out2;

	//TestOut << '\n' << '\n' << "============= [datr] ==============" << '\n' << endl;
	//TestOut << setw(10) << "row  col:";  for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j + 1;
	//TestOut << endl;
	//for (int i = 0; i < nrow_x; i++)
	//{
	//	TestOut << setw(10) << i + 1;
	//	for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << datr[i][j];
	//	TestOut << endl;
	//}



	int n_position = 0;

	string s_temp;

	char buffer[25600000]; //storage for one line

	for (;;) // main loop
	{
		// read one line into buffer using member fn. getline
		pFileInput->getline(buffer, 25600000);
		s_temp = buffer;

		//================================================
		//k category vector 
		//================================================
		n_position = 0;
		n_position = s_temp.find("category", 0);
		if (n_position != -1) //note: -1 means no matched input
		{

			int n_data = ncol_x;  //size of col
			double* d_out = new double[n_data];

			for (int j = 0; j < 1; j++)
			{

				// read one line into buffer using member fn. getline
				pFileInput->getline(buffer, 25600000);

				parse_dString(buffer, d_out, n_data);

				//save coordinates of current node
				for (int i = 0; i < n_data; i++)  k[i] = d_out[i];
			}
			delete[] d_out;


			 //TestOut << '\n' << '\n' << "k category per column" << '\n' << endl;
			 //TestOut << "column : ";
			 //for (int j = 0; j < ncol_x; j++) TestOut << setw(10)<<j+1 ;
			 //TestOut << endl;
			 //TestOut << "k[]    : ";
			 //for (int j = 0; j < ncol_x; j++) TestOut << setw(10)<<k[j] ;
			 //TestOut << endl;

		}

		//================================================
		//NonCollapsible_categorical
		//================================================
		n_position = 0;
		n_position = s_temp.find("NonCollapsible_categorical", 0);
		if (n_position != -1) //note: -1 means no matched input
		{

			int n_data = ncol_x;  //size of col
			double* d_out = new double[n_data];

			for (int j = 0; j < 1; j++)
			{

				// read one line into buffer using member fn. getline
				pFileInput->getline(buffer, 25600000);

				parse_dString(buffer, d_out, n_data);

				//save coordinates of current node
				for (int i = 0; i < n_data; i++)  NonCollapsible_categorical[i] = d_out[i];
			}
			delete[] d_out;


			//TestOut << '\n' << '\n' << "NonCollapsible_categorical" << '\n' << endl;
			//TestOut << "column : ";
			//for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j + 1;
			//TestOut << endl;
			//TestOut << "NonCollapsible_categorical[]    : ";
			//for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << NonCollapsible_categorical[j];
			//TestOut << endl;

		}

		//================================================
		//d weight vector 
		//================================================
		n_position = 0;
		n_position = s_temp.find("weight", 0);
		if (n_position != -1) //note: -1 means no matched input
		{

			int n_data = 1;  //size of col
			double* d_out = new double[n_data];

			for (int j = 0; j < nrow_x; j++)
			{

				// read one line into buffer using member fn. getline
				pFileInput->getline(buffer, 25600000);

				parse_dString(buffer, d_out, n_data);

				//save coordinates of current node
				for (int i = 0; i < n_data; i++)  d[j] = d_out[i];
			}
			delete[] d_out;


			/* TestOut << '\n' << '\n' << "weight per rows" << '\n' << endl;
			 TestOut << "row         : ";
			 for (int j = 0; j < nrow_x; j++) TestOut << setw(10)<<j+1 ;
			 TestOut << endl;
			 TestOut << "weight[]    : ";
			 for (int j = 0; j < nrow_x; j++) TestOut << setw(10)<<d[j] ;
			 TestOut << endl;*/

		}



		//================================================
		//id of each row 
		//================================================
		//n_position = 0;
		//n_position = s_temp.find("index", 0);
		//if (n_position != -1) //note: -1 means no matched input
		//{

		   // int n_data = 1;  
		   // int* i_out = new int[n_data];

		   // for (int j = 0; j<nrow_x; j++)
		   // {
		   //	 // read one line into buffer using member fn. getline
		   //	 pFileInput->getline(buffer, 256000);

		   //	 parse_iString(buffer, i_out, n_data);

		   //	 //save coordinates of current node
		   //	 for (int i = 0; i<n_data; i++)  id[j] = i_out[i];
		   // }
		   // delete[] i_out;

		   // TestOut << '\n' << '\n' << "id number of each row" << '\n' << endl;
		   // TestOut << "row     : ";
		   // for (int j = 0; j < nrow_x; j++) TestOut << setw(10)<<j ;
		   // TestOut << endl;
		   // TestOut << "id[]    : ";
		   // for (int j = 0; j < nrow_x; j++) TestOut << setw(10)<<id[j] ;
		   // TestOut << endl;

		//}

		//--------------
		//read user defined datz, if needed
		//--------------
		if (i_user_defined_datz == 1)
		{
			n_position = 0;
			n_position = s_temp.find("datz", 0);
			if (n_position != -1) //note: -1 means no matched input
			{

				int n_data = ncol_x;  //size of col
				double* d_out = new double[n_data];

				for (int j = 0; j < nrow_x; j++)
				{

					// read one line into buffer using member fn. getline
					pFileInput->getline(buffer, 25600000);

					parse_dString(buffer, d_out, n_data);

					//save coordinates of current node
					for (int i = 0; i < n_data; i++)  datz[j][i] = d_out[i];
				}
				delete[] d_out;

				/* TestOut << '\n' << '\n' << "============= user-defined [datz] ==============" << '\n' << endl;
				 TestOut << setw(10) << "row  col:";  for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j;
				 TestOut << endl;
				 for (int i = 0; i<nrow_x; i++)
				 {
					 TestOut << setw(10) << i + 1;
					 for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << datz[i][j];
					 TestOut << endl;
				 }*/
			}
		}

		//========
		//end of INPUT INFORMATION
		//========
		n_position = 0;
		n_position = s_temp.find("END DATA INPUT", 0);
		if (n_position != -1) //note: -1 means no mateched input
		{
			return; // end of this input information
		}

	} // end of main loop

	return;


}



void InputData_Seperate_ultra(ifstream* pFileInput,
    const int nrow_x,
	const int ncol_x, 

	double *k,
	double *d,

	int    *NonCollapsible_categorical,

	ofstream& TestOut)
	//Description==========================
	//  read input file
	//  and save informations into proper
	//  arrays
	//IN   : ifstream& pFileInput = pointer to the object of ifstream class
	//                              which indicate input file name
	//OUT  : double **daty = New_dMatrix(nrow_x, ncol_x);  original data with missing units
	//OUT  : int    **datr = New_iMatrix(nrow_x, ncol_x); 	indices of 0 for missing cell and 1 for the observed cell
	//OUT  : double **datz = New_dMatrix(nrow_x, ncol_x);  user-defined datz. Vaild only when i_user_defined_datz = 1
	//OUT  : double *k = new double[ncol_x];		category for each column 
	//OUT  : double *d = new double[nrow_x];  		weight for each row 
	//OUT  : int    *id = new int[nrow_x]; 			id of each row 
	//=====================================
{
	
	int n_position = 0;

	string s_temp;

	char buffer[90000000]; //storage for one line

	for (;;) // main loop
	{
		// read one line into buffer using member fn. getline
		pFileInput->getline(buffer, 90000000);
		s_temp = buffer;

		//================================================
		//k category vector 
		//================================================
		n_position = 0;
		n_position = s_temp.find("category", 0);
		if (n_position != -1) //note: -1 means no matched input
		{

			int n_data = ncol_x;  //size of col
			double* d_out = new double[n_data];

			for (int j = 0; j < 1; j++)
			{

				// read one line into buffer using member fn. getline
				pFileInput->getline(buffer, 90000000);

				parse_dString(buffer, d_out, n_data);

				//save coordinates of current node
				for (int i = 0; i < n_data; i++)  k[i] = d_out[i];
			}
			delete[] d_out;


			//TestOut << '\n' << '\n' << "k category per column" << '\n' << endl;
			//TestOut << "column : ";
			//for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j + 1;
			//TestOut << endl;
			//TestOut << "k[]    : ";
			//for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << k[j];
			//TestOut << endl;

		}

		//================================================
		//NonCollapsible_categorical
		//================================================
		n_position = 0;
		n_position = s_temp.find("NonCollapsible_categorical", 0);
		if (n_position != -1) //note: -1 means no matched input
		{

			int n_data = ncol_x;  //size of col
			double* d_out = new double[n_data];

			for (int j = 0; j < 1; j++)
			{

				// read one line into buffer using member fn. getline
				pFileInput->getline(buffer, 90000000);

				parse_dString(buffer, d_out, n_data);

				//save coordinates of current node
				for (int i = 0; i < n_data; i++)  NonCollapsible_categorical[i] = d_out[i];
			}
			delete[] d_out;


			//TestOut << '\n' << '\n' << "NonCollapsible_categorical" << '\n' << endl;
			//TestOut << "column : ";
			//for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j + 1;
			//TestOut << endl;
			//TestOut << "NonCollapsible_categorical[]    : ";
			//for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << NonCollapsible_categorical[j];
			//TestOut << endl;

		}

		//================================================
		//d weight vector 
		//================================================
		n_position = 0;
		n_position = s_temp.find("weight", 0);
		if (n_position != -1) //note: -1 means no matched input
		{

			int n_data = 1;  //size of col
			double* d_out = new double[n_data];

			for (int j = 0; j < nrow_x; j++)
			{

				// read one line into buffer using member fn. getline
				pFileInput->getline(buffer, 90000000);

				parse_dString(buffer, d_out, n_data);

				//save coordinates of current node
				for (int i = 0; i < n_data; i++)  d[j] = d_out[i];
			}
			delete[] d_out;


			/* TestOut << '\n' << '\n' << "weight per rows" << '\n' << endl;
			TestOut << "row         : ";
			for (int j = 0; j < nrow_x; j++) TestOut << setw(10)<<j+1 ;
			TestOut << endl;
			TestOut << "weight[]    : ";
			for (int j = 0; j < nrow_x; j++) TestOut << setw(10)<<d[j] ;
			TestOut << endl;*/

		}


		//========
		//end of INPUT INFORMATION
		//========
		n_position = 0;
		n_position = s_temp.find("END DATA INPUT", 0);
		if (n_position != -1) //note: -1 means no mateched input
		{
			return; // end of this input information
		}

	} // end of main loop

	return;


}



void InputData_Seperate_ultra_txt(ifstream* pFileInput, ifstream* pFileInput_daty, ifstream* pFileInput_datr,
	const int nrow_x,
	const int ncol_x, 
	double *k,
	double *d,

	int    *NonCollapsible_categorical,

	ofstream& TestOut)
	//Description==========================
	//  read input file
	//  and save informations into proper
	//  arrays
	//IN   : ifstream& pFileInput = pointer to the object of ifstream class
	//                              which indicate input file name
	//OUT  : double **daty = New_dMatrix(nrow_x, ncol_x);  original data with missing units
	//OUT  : int    **datr = New_iMatrix(nrow_x, ncol_x); 	indices of 0 for missing cell and 1 for the observed cell
	//OUT  : double **datz = New_dMatrix(nrow_x, ncol_x);  user-defined datz. Vaild only when i_user_defined_datz = 1
	//OUT  : double *k = new double[ncol_x];		category for each column 
	//OUT  : double *d = new double[nrow_x];  		weight for each row 
	//OUT  : int    *id = new int[nrow_x]; 			id of each row 
	//=====================================
{
	//-- MPI variables
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	MPI_File fh_binary_daty_column;// File to write daty column-wisely
	MPI_File fh_binary_daty_row;// File to write daty row-wisely
	MPI_File fh_binary_datr;// File to write datr

	int success = 0;// inidicator of accessing file system

	success = MPI_File_open(MPI_COMM_WORLD, "./daty_column_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_binary_daty_column);
	if (success != MPI_SUCCESS) cout << "MPI I/O fail to open the daty_column_binary file in ReadInput!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./daty_row_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_binary_daty_row);
	if (success != MPI_SUCCESS) cout << "MPI I/O fail to open the daty_row_binary file in ReadInput!" << endl;

	success = MPI_File_open(MPI_COMM_WORLD, "./datr_column_binary.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_binary_datr);
	if (success != MPI_SUCCESS) cout << "MPI I/O fail to open the datr_column_binary file in ReadInput!" << endl;


	string s_temp_daty;
	string s_temp_datr;
	int n_position = 0;

	if (mynode == 0) {

		double ** daty = New_dMatrix(nrow_x, ncol_x);
		int ** datr = New_iMatrix(nrow_x, ncol_x);

		//================================================
		//daty original data matrix 
		//================================================
		char buffer_daty[25600000]; //storage for one line

		// read one line into buffer using member fn. getline
		pFileInput_daty->getline(buffer_daty, 25600000);
		s_temp_daty = buffer_daty;

		n_position = 0;
		n_position = s_temp_daty.find("daty", 0);

		if (n_position != -1) {

			int n_data = ncol_x;  //size of col
			double* d_out = new double[n_data];

			for (int j = 0; j < nrow_x; j++)
			{

				// read one line into buffer using member fn. getline
				pFileInput_daty->getline(buffer_daty, 25600000);

				parse_dString(buffer_daty, d_out, n_data);

				//save coordinates of current node
				for (int i = 0; i < n_data; i++)  daty[j][i] = d_out[i];
			}
			delete[] d_out;

			//TestOut << '\n' << '\n' << "============= [daty] ==============" << '\n' << endl;
			//TestOut << setw(10) << "row  col:";  for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j + 1;
			//TestOut << endl;
			//for (int i = 0; i < nrow_x; i++)
			//{
			//	TestOut << setw(10) << i + 1;
			//	for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << daty[i][j];
			//	TestOut << endl;
			//}
		}
	
		//================================================
		//datr original data matrix 
		//================================================
		char buffer_datr[25600000]; //storage for one line

		// read one line into buffer using member fn. getline
		pFileInput_datr->getline(buffer_datr, 25600000);
		s_temp_datr = buffer_datr;

		n_position = 0;
		n_position = s_temp_datr.find("datr", 0);

		if (n_position != -1) {

			int n_data = ncol_x;  //size of col
			int* i_out = new int[n_data];

			for (int j = 0; j<nrow_x; j++)
			{
				// read one line into buffer using member fn. getline
				pFileInput_datr->getline(buffer_datr, 25600000);

				//int n_data=4;  //size of col
				//int i_out[n_data];

				parse_iString(buffer_datr, i_out, n_data);

				//save coordinates of current node
				for (int i = 0; i<n_data; i++)  datr[j][i] = i_out[i];
			}
			delete[] i_out;

			//TestOut << '\n' << '\n' << "============= [datr] ==============" << '\n' << endl;
			//TestOut << setw(10) << "row  col:";  for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j + 1;
			//TestOut << endl;
			//for (int i = 0; i < nrow_x; i++)
			//{
			//	TestOut << setw(10) << i + 1;
			//	for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << datr[i][j];
			//	TestOut << endl;
			//}
		
		}

		//-------------------------------------
		//Write daty and datr to HD
		//--------------------------------------

		MPI_Out_sorted_daty_column(nrow_x, ncol_x, fh_binary_daty_column, daty);
		MPI_Out_sorted_daty_row(nrow_x, ncol_x, fh_binary_daty_row, daty);
		MPI_Out_sorted_datr(nrow_x, ncol_x, fh_binary_datr, datr);


		//Deallocation
		Del_dMatrix(daty, nrow_x, ncol_x);
		Del_iMatrix(datr, nrow_x, ncol_x);

	}//end of mynode==0!!!!!!!!!!!!!!!!
	
	//---------------
	//Close all files
	//---------------
	success = MPI_File_close(&fh_binary_daty_column);
	if (success != MPI_SUCCESS) cout << "MPI I/O of daty fail to close the file!" << endl;

	success = MPI_File_close(&fh_binary_daty_row);
	if (success != MPI_SUCCESS) cout << "MPI I/O of daty fail to close the file!" << endl;

	success = MPI_File_close(&fh_binary_datr);
	if (success != MPI_SUCCESS) cout << "MPI I/O of datr fail to close the file!" << endl;



	n_position = 0;

	string s_temp;

	char buffer[90000000]; //storage for one line

	for (;;) // main loop
	{
		// read one line into buffer using member fn. getline
		pFileInput->getline(buffer, 90000000);
		s_temp = buffer;

		//================================================
		//k category vector 
		//================================================
		n_position = 0;
		n_position = s_temp.find("category", 0);
		if (n_position != -1) //note: -1 means no matched input
		{

			int n_data = ncol_x;  //size of col
			double* d_out = new double[n_data];

			for (int j = 0; j < 1; j++)
			{

				// read one line into buffer using member fn. getline
				pFileInput->getline(buffer, 90000000);

				parse_dString(buffer, d_out, n_data);

				//save coordinates of current node
				for (int i = 0; i < n_data; i++)  k[i] = d_out[i];
			}
			delete[] d_out;


			//TestOut << '\n' << '\n' << "k category per column" << '\n' << endl;
			//TestOut << "column : ";
			//for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j + 1;
			//TestOut << endl;
			//TestOut << "k[]    : ";
			//for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << k[j];
			//TestOut << endl;

		}

		//================================================
		//NonCollapsible_categorical
		//================================================
		n_position = 0;
		n_position = s_temp.find("NonCollapsible_categorical", 0);
		if (n_position != -1) //note: -1 means no matched input
		{

			int n_data = ncol_x;  //size of col
			double* d_out = new double[n_data];

			for (int j = 0; j < 1; j++)
			{

				// read one line into buffer using member fn. getline
				pFileInput->getline(buffer, 90000000);

				parse_dString(buffer, d_out, n_data);

				//save coordinates of current node
				for (int i = 0; i < n_data; i++)  NonCollapsible_categorical[i] = d_out[i];
			}
			delete[] d_out;


			//TestOut << '\n' << '\n' << "NonCollapsible_categorical" << '\n' << endl;
			//TestOut << "column : ";
			//for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << j + 1;
			//TestOut << endl;
			//TestOut << "NonCollapsible_categorical[]    : ";
			//for (int j = 0; j < ncol_x; j++) TestOut << setw(10) << NonCollapsible_categorical[j];
			//TestOut << endl;

		}

		//================================================
		//d weight vector 
		//================================================
		n_position = 0;
		n_position = s_temp.find("weight", 0);
		if (n_position != -1) //note: -1 means no matched input
		{

			int n_data = 1;  //size of col
			double* d_out = new double[n_data];

			for (int j = 0; j < nrow_x; j++)
			{

				// read one line into buffer using member fn. getline
				pFileInput->getline(buffer, 90000000);

				parse_dString(buffer, d_out, n_data);

				//save coordinates of current node
				for (int i = 0; i < n_data; i++)  d[j] = d_out[i];
			}
			delete[] d_out;


			/* TestOut << '\n' << '\n' << "weight per rows" << '\n' << endl;
			TestOut << "row         : ";
			for (int j = 0; j < nrow_x; j++) TestOut << setw(10)<<j+1 ;
			TestOut << endl;
			TestOut << "weight[]    : ";
			for (int j = 0; j < nrow_x; j++) TestOut << setw(10)<<d[j] ;
			TestOut << endl;*/

		}


		//========
		//end of INPUT INFORMATION
		//========
		n_position = 0;
		n_position = s_temp.find("END DATA INPUT", 0);
		if (n_position != -1) //note: -1 means no mateched input
		{
			return; // end of this input information
		}

	} // end of main loop

	return;
}