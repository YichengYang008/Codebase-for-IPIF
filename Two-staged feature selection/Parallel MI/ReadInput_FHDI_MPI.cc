void ReadInput_FHDI_MPI(ifstream* pFileInput, int &nrow, int &ncol,
               int &K, int &i_col_target, int &i_SIS)
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
		   n_position = s_temp.find("nrow", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   istrstream one_line(buffer); // construct an istringstream class object with buffer
									  // note: new version= istringstream
			   one_line >> c_temp;
			   one_line >> nrow;

			  // TestOut << setw(30) << c_temp <<setw(10) << " = " <<setw(10) << nrow <<  "    (0:input from input.txt; 1:seperate input from daty.txt, datr.txt, and datz.txt) "  << endl;
		   }

		   //========
		   //i_option_ultra
		   //========
		   n_position = 0;
		   n_position = s_temp.find("ncol", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   istrstream one_line(buffer); // construct an istringstream class object with buffer
											// note: new version= istringstream
			   one_line >> c_temp;
			   one_line >> ncol;

			   //TestOut << setw(30) << c_temp << setw(10) << " = " << setw(10) << ncol << "    (0:P-FHDI for big-n or big-p data based on memory ; 1:P-FHDI for ultra data based on memory and local storage) " << endl;
		   }


           //========
           //i_option_perform
           //========
           n_position=0;
           n_position = s_temp.find("category",0);
           if(n_position != -1) //note: -1 means no mateched input
           {
              istrstream one_line(buffer); // construct an istringstream class object with buffer
                                     // note: new version= istringstream
              one_line >> c_temp  ;
              one_line >> K;

			 //TestOut << setw(30) << c_temp <<setw(10) << " = " <<setw(10) << i_option_perform <<  "    (1:do all; 2:CellMake; 3:CellProb; 4: do all with user's datz) "  << endl;
           }

           //========
           //nrow
           //========
           n_position=0;
           n_position = s_temp.find("i_col_target",0);
           if(n_position != -1) //note: -1 means no mateched input
           {
              istrstream one_line(buffer); // construct an istringstream class object with buffer
                                     // note: new version= istringstream
              one_line >> c_temp  ;
              one_line >> i_col_target;

             // TestOut <<setw(30)<< c_temp <<setw(10)<<" = "<<setw(10)<< nrow_x << endl;
           }

           //========
           //ncol
           //========
           n_position=0;
           n_position = s_temp.find("i_selection",0);
           if(n_position != -1) //note: -1 means no mateched input
           {
              istrstream one_line(buffer); // construct an istringstream class object with buffer
                                     // note: new version= istringstream
              one_line >> c_temp  ;
              one_line >> i_SIS;

             // TestOut <<setw(30)<< c_temp <<setw(10)<<" = "<<setw(10)<< ncol_x << endl;
           }

		   //========
		   //end of INPUT INFORMATION
		   //========
		   n_position = 0;
		   n_position = s_temp.find("END INPUT INFORMATION", 0);
		   if (n_position != -1) //note: -1 means no mateched input
		   {
			   return; // end of this input information
		   }

        }//end for
     }


   } // end of main loop


  return;
}

