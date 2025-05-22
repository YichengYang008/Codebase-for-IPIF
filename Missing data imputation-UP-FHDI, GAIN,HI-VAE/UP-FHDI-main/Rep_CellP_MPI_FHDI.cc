#include "matrix_utility_FHDI.cc"
#include <vector>
#include <string>
#include <mpi.h>
#include <iostream>
//#include "Cell_Prob_Extension_cpp.cc"

void Rep_CellP_FHDI(double** d_cx, const int nrow, const int ncol,
	RepWeight_FHDI &d_rw, int*  id,
	//List_FHDI        &List_rst_prob,
	//List_string_FHDI &List_rst_name,
	std::vector<int> &List_rst_prob_size,
	std::vector<std::string> &s_ncx, ofstream& TestOut)
//Description============================================
// compute cell probability using replicate weight rw
// 
// R code: Dr. Im, J., and Dr. Kim, J. 
// C++   : Dr. Cho, I.
// All rights reserved
// Last update: March 28, 2017
//

//IN   : double d_cx[nrow, ncol] = categoraized matrix
//IN   : double d_rw[nrow, nrow] = replicate weights
//IN   : int    id[nrow] = index of rows
//
//below two lists have meaningful values up to i_nc rows  
//OUT  : List_FHDI List_rst_prob(nrow->i_nc); //list of joint probabilities for all missing patterns 
//OUT  : List_string_FHDI List_rst_name(nrow->i_nc); //names of joint probabilities for all missing patterns 
//OUT  : std::vector<std::string> s_ncx; //uniqe cn0
//======================================================== 
{
	//--------------
	//make a condensed expression "cn0" of cx, i.e. z
	//--------------
	//std::string cn0[nrow];
	//std::string cn0_backup[nrow];
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;
	//cout << "Mynode: " << mynode<<", Rep_Debug1"<<endl;
	double Rep_CellP_begin = MPI_Wtime();
	std::string *cn0 = new std::string[nrow];
	//std::string *cn0_backup = new std::string[nrow];

	Trans(d_cx, nrow, ncol, cn0);//make a condensed string expression with a given array

								 //for (int i = 0; i<nrow; i++) cn0_backup[i] = cn0[i];

								 //---------------------
								 //SORT & UNIQUE patterns of cn0
								 //---------------------
								 //std::string s_cn0_temp[nrow]; 
	std::string *s_cn0_temp = new std::string[nrow];
	for (int i = 0; i<nrow; i++) s_cn0_temp[i] = cn0[i];
	std::sort(s_cn0_temp, s_cn0_temp + nrow);
	//cout << "Mynode: " << mynode << ", Rep_Debug2" << endl;
	//------------
	//memorize observed patterns 
	//------------
	//std::vector<std::string> s_ncx; //uniqe cn0

	int i_count_cn0 = 0; //total number of unique cn0 
	std::string s_temp;
	//for (int i = 0; i<nrow; i++)
	//{
	//	s_temp = s_cn0_temp[i]; //get a string from the sorted strings 
	//	for (int j = 0; j<nrow; j++) //search all rows 
	//	{
	//		//----
	//		//below condition is needed for finding UNIQUE pattern
	//		//----
	//		if (s_temp.compare(cn0_backup[j]) == 0) //0: equal string
	//		{
	//			s_ncx.push_back(cn0_backup[j]);  //store the found observed pattern

	//											 //----------
	//											 //remove all identical string after the current string
	//											 //----------
	//			for (int k = j; k<nrow; k++)
	//			{
	//				if (s_temp.compare(cn0_backup[k]) == 0) //0: equal string
	//				{
	//					cn0_backup[k] = ""; //nullify for the next search
	//				}
	//			}

	//			i_count_cn0++;
	//			break;
	//		}

	//	}
	//}

	for (int i = 0; i < nrow; i++) {
		s_temp = s_cn0_temp[i]; //get a string from the sorted strings 
		if (i == 0 || (i > 0 && s_temp.compare(s_cn0_temp[i - 1]) != 0)) {
			s_ncx.push_back(s_temp);
			i_count_cn0++;
		}
	}

	//Now, i_count_cn0 means the total number of unique sorted strings
	const int i_nc = i_count_cn0;
	//cout << "Yang Running time of Rep_cell1 at node " << mynode << " is " << MPI_Wtime() - Rep_CellP_begin << endl;
	//testout 

	//testout
	//RPrint("=====in Rep_CellP ========");
	//RPrint("s_ncx"); RPrint(s_ncx);
	//RPrint("i_nc"); RPrint(i_nc);
	/*cout<<"=====in Rep_CellP ========"<<endl;
	cout<<"s_ncx"<<endl;
	for(int i=0; i<(int)s_ncx.size(); i++) cout<<s_ncx[i]<<" ,  ";
	cout<<endl;
	cout<<"i_nc: "<<i_nc<<endl;
	*/

	//-----------------------------
	//calculate joint probability and names of all missing patterns
	//using the Jackknife replicate weights
	//------------------------------
	//List_FHDI        List_rst_prob(i_nc); //list of joint probabilities for all missing patterns 
	//List_string_FHDI List_rst_name(i_nc); //names of joint probabilities for all missing patterns 

	std::vector<double> jp_prob_return; 
	std::vector<std::string> jp_name_return;

	std::vector<int> List_rst_prob_size_send;
	//std::vector<double> w_UserDefined; 
	double* w_UserDefined = new double[nrow]; 

	double Rep_CellP_begin2 = MPI_Wtime();
	//if (mynode == 0) cout << "The value of i_nc in RepCellP is " << i_nc << endl;
	if (i_nc < totalnodes) {
		TestOut << "i_nc is smaller than totalnodes, ERROR!!!" << endl;
		return;
	}
	//------------------------------------------------------ Job assignment
	int number_before_split = 0;
	int number_after_split = 0;
	int startpoint = 0;
	int endpoint = 0;
	int splitnode = 0;

	if (i_nc % (totalnodes - 1) != 0) { splitnode = i_nc % (totalnodes - 1); }
	if (i_nc % (totalnodes - 1) == 0) { splitnode = 1; }
	//if (mynode == 0) cout << "Split: " << splitnode << endl;
	number_after_split = floor(1.0*i_nc / (1.0*totalnodes - 1));
	number_before_split = 1.0*(i_nc - floor(1.0*i_nc / (1.0*totalnodes - 1)) *(totalnodes - splitnode - 1)) / splitnode;

	if (mynode >= 1 && mynode <= splitnode) {
		startpoint = (mynode - 1) * number_before_split;
		endpoint = mynode *number_before_split;
	}
	if (mynode > splitnode) {
		startpoint = splitnode * number_before_split + (mynode - splitnode - 1)*number_after_split;
		endpoint = splitnode * number_before_split + (mynode - splitnode)*number_after_split;
	}
	if ((number_before_split*splitnode + number_after_split* (totalnodes - splitnode - 1)) != i_nc) {
		TestOut << "Work Assignment Error!!!!" << endl;
		return;
	}

	for (int i = startpoint; i<endpoint; i++)
	{
		//---
		//search current missing pattern from all strings
		//---
		std::string s_temp = s_ncx[i]; 
		int i_loc = 0; 
		for(int j=0; j<nrow; j++) 
		{
			if(s_temp.compare(cn0[j])==0)
			{
				i_loc = j; 
				break; 
			}
		}
		//testout
		/*
		cout<<"loop i (1:i_nc) :"<<i<<"  found i_loc:"<<i_loc<<endl;
		if(i==7)
		{
			for(int j_temp=0; j_temp<nrow; j_temp++) cout<<d_rw[j_temp][i_loc]<<",  "; 
		}
		cout<<endl;
		*/
		
		//----
		//joint probability and names
		//----
		jp_prob_return.clear(); 
		jp_name_return.clear();
		//w_UserDefined.clear();
		//for(int j=0; j<nrow; j++) w_UserDefined.push_back(d_rw[j][i_loc]) ; 
		for (int j = 0; j<nrow; j++) w_UserDefined[j] = d_rw(j, i_loc);
		
		Cell_Prob_Extension_cpp(d_cx, nrow, ncol,
							    jp_prob_return,
							    jp_name_return, 
							    w_UserDefined, id, TestOut);		
								
		//---
		//prep return
		//---
		//List_rst_prob.put_block(i, jp_prob_return); //jth row has joint prob
		//List_rst_name.put_block(i, jp_name_return); //jth row has name of the joint prob
		List_rst_prob_size_send.push_back(jp_prob_return.size());
	}
	
	if (mynode != 0) {
		MPI_Send(&List_rst_prob_size_send[0], List_rst_prob_size_send.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);
	}

	std::vector<int> List_rst_prob_size_recv_before(number_before_split);
	std::vector<int> List_rst_prob_size_recv_after(number_after_split);

	if (mynode == 0) {
		for (int j = 1; j < totalnodes; j = j + 1) {
			//cout << "J is at "<< j <<endl;
			if (j >= 1 && j <= splitnode) {
				MPI_Recv(&List_rst_prob_size_recv_before[0], number_before_split, MPI_INT, j, 1, MPI_COMM_WORLD, &status);

				for (int k = 0; k < number_before_split;k++) {
					List_rst_prob_size.push_back(List_rst_prob_size_recv_before[k]);
				}
			}//end of j>1

			if (j > splitnode) {
				MPI_Recv(&List_rst_prob_size_recv_after[0], number_after_split, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				
				for (int k = 0; k < number_after_split;k++) {
					List_rst_prob_size.push_back(List_rst_prob_size_recv_after[k]);
				}

			}//end of j > splitnode
		}//end of j
	}//end of mynode=0

	if (mynode != 0) { 
		List_rst_prob_size.resize(i_nc); 
	}

	MPI_Bcast(&List_rst_prob_size[0], i_nc, MPI_INT, 0, MPI_COMM_WORLD);
	//testout
	//RPrint("List_rst_name"); List_rst_name.print_List_string_FHDI();
	//RPrint("List_rst_prob"); List_rst_prob.print_List_FHDI();

	//cout<<"List_rst_prob_size is"<< List_rst_prob_size.size() <<" at node "<<mynode<<endl;
	//for (int k1 = 0; k1 < List_rst_prob_size.size();k1++) {
	//	cout << "List_rst_prob_size[" <<k1<<"]: "<< List_rst_prob_size[k1]<< " at node "<<mynode<<endl;
	//}

	delete[] cn0; 
	//delete[] cn0_backup; 
	delete[] s_cn0_temp;
	
	delete[] w_UserDefined;
	
	return;				   
}
