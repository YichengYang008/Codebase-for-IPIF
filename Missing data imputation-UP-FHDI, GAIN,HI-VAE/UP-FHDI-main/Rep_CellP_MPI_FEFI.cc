#include "matrix_utility_FHDI.cc"
#include <vector>
#include <string>
#include <mpi.h>
#include <iostream>
#include "Cell_Prob_Extension_cpp.cc"
//#include "Cell_Prob_Extension_cpp_MPI.cc"

void RepWeight(const int n, double** d_rw)
//Description -------------------------------------
// Jackknife replicate weight for simpler random sampling
// 
// original R code: Dr. Im, J. and Dr. Kim, J. 
// c++ code: 		Dr. Cho, I. 
// All rights reserved
// 
// updated: Nov 17, 2016
//IN   : int n = matrix dimension
//OUT  : double d_rw[n,n] = replicate weights 
//--------------------------------------------------
{


	const double d_rw0 = (1.0*n) / (n - 1);

	//--------
	//initialize
	//--------
	Fill_dMatrix(d_rw, n, n, d_rw0);

	//--------
	//put 0 into diagonal terms
	//--------
	for (int i = 0; i<n; i++)
	{
		d_rw[i][i] = 0.0;
	}

	return;
}


void Rep_CellP_FEFI(double** d_cx, const int nrow, const int ncol, 
	RepWeight_FHDI &d_rw, int*  id,
	List_FHDI        &List_rst_prob,
	List_string_FHDI &List_rst_name,
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

	std::vector<double> jp_prob_return;// vector to hold the output of cell_prob_extension [one row]
	std::vector<std::string> jp_name_return; // vector to hold the output of cell_prob_extension [one row]

											 //send for prob
	std::vector<double> jp_prob_return_send; // send buffer of jp_prob [all rows]
											 //std::vector<double> jp_prob_return_recv(jp_prob_return_send.size());

											 //send for string
	std::vector<std::string> jp_name_return_send_temp;
	std::vector<char> jp_name_return_send; // send buffer of jp_name [all rows]
										   //std::vector<char> jp_name_return_recv(jp_name_return_send.size());

										   //std::vector<double> w_UserDefined; 
	double* w_UserDefined = new double[nrow];
	double Rep_CellP_begin2 = MPI_Wtime();
	//if (mynode == 0) cout << "The value of i_nc in RepCellP is " << i_nc << endl;
	if (i_nc < totalnodes) {
		TestOut<<"i_nc is smaller than totalnodes, ERROR!!!"<<endl;
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
	//if (mynode == 0) {
	//	cout << "Mynode: " << mynode << ", number_before_split: " << number_before_split << ", number_after_split:" << number_after_split << endl;
	//}
	//cout << "Mynode " << mynode << ", " << startpoint << "<= x < " << endpoint << endl;
	//---------------------------------------------------------

	for (int i = startpoint; i<endpoint; i++)
	{
		//---
		//search current missing pattern from all strings
		//---
		std::string s_temp = s_ncx[i];
		int i_loc = 0;
		for (int j = 0; j<nrow; j++)
		{
			if (s_temp.compare(cn0[j]) == 0)
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
		for (int j = 0; j<nrow; j++) w_UserDefined[j] = d_rw(j,i_loc);

		Cell_Prob_Extension_cpp(d_cx, nrow, ncol,
			jp_prob_return,
			jp_name_return,
			w_UserDefined, id, TestOut);

		//for (int i = 0;i < jp_name_return.size();i++) {
		//	jp_name_return_temp.push_back(atoi(jp_name_return[i].c_str()));
		//}
		//---
		//prep return
		//---

		//List_rst_prob.put_block(i, jp_prob_return); //jth row has joint prob
		//List_rst_name.put_block(i, jp_name_return); //jth row has name of the joint prob

		for (int j = 0;j < jp_prob_return.size();j++) {//double
			jp_prob_return_send.push_back(jp_prob_return[j]);
		}

		for (int j = 0;j < jp_name_return.size();j++) {//string
			jp_name_return_send_temp.push_back(jp_name_return[j]);
		}
		//if (mynode == 1) cout << "size of jp_name_return_send_temp is " << jp_name_return_send_temp.size() << endl;

		//for (int i = 0; i < jp_name_return_send_temp.size();i++) {//char
		//	for (int j = 0; j < jp_name_return_send_temp[i].size();j++) {
		//		jp_name_return_send.push_back(jp_name_return_send_temp[i][j]);
		//		cout <<"i_nc: "<<i<<" j: "<<j<<"wow: " << jp_name_return_send_temp[i][j] <<endl;
		//	}
		//}

	} // end of i_nc

	for (int i = 0; i < jp_name_return_send_temp.size();i++) {//char
		for (int j = 0; j < jp_name_return_send_temp[i].size();j++) {
			jp_name_return_send.push_back(jp_name_return_send_temp[i][j]);
			//cout <<"i_nc: "<<i<<" j: "<<j<<"wow: " << jp_name_return_send_temp[i][j] <<endl;
		}
	}

	//cout <<"Mynode: "<< mynode<<", jp_name_return_send(why 0): " << jp_name_return_send.size() << endl;
	//cout <<"Mynode: "<< mynode<< ", jp_prob_return_send(why 0): " << jp_prob_return_send.size() << endl;
	//cout << "Yang Running time of Rep_cell_parallel at node " << mynode << " is " << MPI_Wtime() - Rep_CellP_begin2 << endl;

	//----------------------------------------------------------------------------------------------------
	int jp_prob_size_send = jp_prob_return_send.size();
	int jp_name_size_send = jp_name_return_send.size();
	int jp_prob_size_each = jp_prob_return.size();
	int jp_name_size_each = jp_name_return.size();

	int jp_prob_size_before =0;
	int jp_name_size_before =0;
	int jp_prob_each_standard=0;
	int jp_name_each_standard=0;

	int jp_prob_size_after=0;
	int jp_name_size_after=0;

	//if (mynode != 0) {
	//	MPI_Send(&jp_prob_size_send, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	//	MPI_Send(&jp_name_size_send, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
	//}

	//if (mynode == 0) {
	//	for (int j = 1;j < totalnodes;j++) {
	//		MPI_Recv(&jp_prob_size, 1, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
	//		MPI_Recv(&jp_name_size, 1, MPI_INT, j, 2, MPI_COMM_WORLD, &status);
	//	}
	//}
	double Rep_CellP_begin_communication = MPI_Wtime();
	if (mynode == splitnode) {
		MPI_Send(&jp_prob_size_send, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&jp_name_size_send, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
		MPI_Send(&jp_prob_size_each, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
		MPI_Send(&jp_name_size_each, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);
	}

	if (mynode == (splitnode + 1)) {
		MPI_Send(&jp_prob_size_send, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&jp_name_size_send, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
		//MPI_Send(&jp_prob_size_return, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
		//MPI_Send(&jp_name_size_return, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);
	}

	if (mynode == 0) {
		MPI_Recv(&jp_prob_size_before, 1, MPI_INT, splitnode, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&jp_name_size_before, 1, MPI_INT, splitnode, 2, MPI_COMM_WORLD, &status);
		MPI_Recv(&jp_prob_each_standard, 1, MPI_INT, splitnode, 3, MPI_COMM_WORLD, &status);
		MPI_Recv(&jp_name_each_standard, 1, MPI_INT, splitnode, 4, MPI_COMM_WORLD, &status);


		MPI_Recv(&jp_prob_size_after, 1, MPI_INT, (splitnode + 1), 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&jp_name_size_after, 1, MPI_INT, (splitnode + 1), 2, MPI_COMM_WORLD, &status);
		//MPI_Recv(&jp_prob_each_last, 1, MPI_INT, (totalnodes - 1), 3, MPI_COMM_WORLD, &status);
		//MPI_Recv(&jp_name_each_last, 1, MPI_INT, (totalnodes - 1), 4, MPI_COMM_WORLD, &status);
	}
	//if (mynode == 0) cout << "jp_prob_size_before is " << jp_prob_size_before << endl;
	//if (mynode == 0) cout << "jp_name_size_before is " << jp_name_size_before << endl;
	//if (mynode == 0) cout << "jp_prob_each_standard is " << jp_prob_each_standard << endl;
	//if (mynode == 0) cout << "jp_name_each_standard is " << jp_name_each_standard << endl;

	//if (mynode == 0) cout << "jp_prob_size_after is " << jp_prob_size_after << endl;
	//if (mynode == 0) cout << "jp_name_size_after is " << jp_name_size_after << endl;
	//cout <<"Mynode: "<<mynode <<", jp_prob_size_before is " << jp_prob_size_before << endl;
	//cout << "Mynode: " << mynode << ", jp_name_size_before is " << jp_name_size_before << endl;
	//cout << "Mynode: " << mynode << ", jp_prob_each_standard is " << jp_prob_each_standard << endl;
	//cout << "Mynode: " << mynode << ", jp_name_each_standard is " << jp_name_each_standard << endl;
	//cout << "Mynode: " << mynode << ", jp_prob_size_after is " << jp_prob_size_after << endl;
	//cout << "Mynode: " << mynode << ", jp_name_size_after is " << jp_name_size_after << endl;
	//-----------------------------------------------------------------------------------------------------

	MPI_Barrier(MPI_COMM_WORLD);
	std::vector<double> jp_prob_before_return_recv(jp_prob_size_before);
	std::vector<char> jp_name_before_return_recv(jp_name_size_before);

	std::vector<double> jp_prob_after_return_recv(jp_prob_size_after);
	std::vector<char> jp_name_after_return_recv(jp_name_size_after);

	//cout << "Yang Running time of Rep_cell_communication1 at node " << mynode << " is " << MPI_Wtime() - Rep_CellP_begin_communication << endl;
	if (mynode != 0) {
		MPI_Send(&jp_prob_return_send[0], jp_prob_return_send.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&jp_name_return_send[0], jp_name_return_send.size(), MPI_CHAR, 0, 1, MPI_COMM_WORLD);
	}


	if (mynode == 0) {
		for (int j = 1; j < totalnodes; j = j + 1) {
			//cout << "J is at "<< j <<endl;
			if (j >= 1 && j <= splitnode) {
				MPI_Recv(&jp_prob_before_return_recv[0], jp_prob_size_before, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
				int counter = 0;
				std::vector<double> jp_prob_buffer;
				for (int i = (j - 1)*number_before_split; i < j*number_before_split; i++) {
					jp_prob_buffer.clear();
					for (int k = 0; k < jp_prob_each_standard;k++) {
						jp_prob_buffer.push_back(jp_prob_before_return_recv[counter]);
						//cout << "counter: " << counter << "jp_prob_before_return_recv: " << jp_prob_before_return_recv[counter] << endl;
						counter++;
					}

					//for (int k = 0;k < jp_prob_buffer.size();k++) {
					//	cout << "jp: " << jp_prob_buffer[k] << endl;
					//}
					//cout << "jp_prob_buffer size: " << jp_prob_buffer.size() << endl;
					//cout << "I'm here" << endl;
					List_rst_prob.put_block(i, jp_prob_buffer);
				}

				MPI_Recv(&jp_name_before_return_recv[0], jp_name_size_before, MPI_CHAR, j, 1, MPI_COMM_WORLD, &status);
				std::vector<std::string> jp_name_return_recv_temp;// the actual string vector for jp_name [all rows]
				int name_counter = 0;
				string jp_name_temp;
				int jp_name_return_send_temp_size = jp_name_size_before / ncol; // !!Attension, should be same as jp_name_reurn _send_temp
																		 //cout<<"jp_name_return_send_temp_size: "<< jp_name_return_send_temp_size <<endl;

				for (int k = 0;k < jp_name_return_send_temp_size;k++) {
					jp_name_temp = "";
					for (int i = 0;i < ncol;i++) {
						jp_name_temp = jp_name_temp + jp_name_before_return_recv[name_counter];
						name_counter++;
					}
					jp_name_return_recv_temp.push_back(jp_name_temp);
				}

				//cout<<"jp_name_return_recv_temp_before: "<< jp_name_return_recv_temp.size()<<endl;
				//---------------------
				int counter_yicheng = 0;
				std::vector<std::string> jp_name_buffer;
				for (int i = (j - 1)*number_before_split; i < j*number_before_split; i++) {
					jp_name_buffer.clear();
					for (int k = 0; k < jp_name_each_standard;k++) {
						jp_name_buffer.push_back(jp_name_return_recv_temp[counter_yicheng]);
						//cout << "Counter_name: " << counter << ", jp_name_return_recv_temp_before: " << jp_name_return_recv_temp[counter_yicheng]<<endl;
						counter_yicheng++;
					}
					List_rst_name.put_block(i, jp_name_buffer);
				}

			}//end of j>1

			if (j > splitnode) {
				MPI_Recv(&jp_prob_after_return_recv[0], jp_prob_size_after, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
				int counter = 0;
				std::vector<double> jp_prob_buffer;
				for (int i = splitnode * number_before_split + (j - splitnode - 1)*number_after_split; i < splitnode * number_before_split + (j - splitnode)*number_after_split; i++) {
					jp_prob_buffer.clear();
					for (int k = 0; k < jp_prob_each_standard;k++) {
						jp_prob_buffer.push_back(jp_prob_after_return_recv[counter]);
						//cout << "counter: " << counter << "jp_prob_after_return_recv: " << jp_prob_after_return_recv[counter] << endl;
						counter++;
					}

					//for (int k = 0;k < jp_prob_buffer.size();k++) {
					//	cout << "jp: " << jp_prob_buffer[k] << endl;
					//}
					//cout << "jp_prob_buffer size: " << jp_prob_buffer.size() << endl;
					//cout << "I'm here" << endl;
					List_rst_prob.put_block(i, jp_prob_buffer);
				}

				MPI_Recv(&jp_name_after_return_recv[0], jp_name_size_after, MPI_CHAR, j, 1, MPI_COMM_WORLD, &status);
				std::vector<std::string> jp_name_return_recv_temp;// the actual string vector for jp_name [all rows]
				int name_counter = 0;
				string jp_name_temp;
				int jp_name_return_send_temp_size = jp_name_size_after / ncol; // !!Attension, should be same as jp_name_reurn _send_temp
																				//cout<<"jp_name_return_send_temp_size: "<< jp_name_return_send_temp_size <<endl;

				for (int k = 0;k < jp_name_return_send_temp_size;k++) {
					jp_name_temp = "";
					for (int i = 0;i < ncol;i++) {
						jp_name_temp = jp_name_temp + jp_name_after_return_recv[name_counter];
						name_counter++;
					}
					jp_name_return_recv_temp.push_back(jp_name_temp);
				}

				//cout << "jp_name_return_recv_temp_after: " << jp_name_return_recv_temp.size() << endl;
				//---------------------
				int counter_yicheng = 0;
				std::vector<std::string> jp_name_buffer;
				for (int i = splitnode * number_before_split + (j - splitnode - 1)*number_after_split; i < splitnode * number_before_split + (j - splitnode)*number_after_split; i++) {
					jp_name_buffer.clear();
					for (int k = 0; k < jp_name_each_standard;k++) {
						jp_name_buffer.push_back(jp_name_return_recv_temp[counter_yicheng]);
						//cout << "Counter_name: " << counter << ", jp_name_return_recv_temp_after: " << jp_name_return_recv_temp[counter_yicheng] << endl;
						counter_yicheng++;
					}
					List_rst_name.put_block(i, jp_name_buffer);
				}

			}//end of j > splitnode
		}//end of j
	}//end of mynode=0

	//cout << "Yang Running time of Rep_cell_communication2 at node " << mynode << " is " << MPI_Wtime() - Rep_CellP_begin_communication << endl;
	//MPI_Barrier(MPI_COMM_WORLD);
	//testout
	//if (mynode == 0) { 
	//	cout << "List_rst_name" << endl;;
	//	List_rst_name.print_List_string_FHDI();
	//	cout << "List_rst_prob" << endl;
	//	List_rst_prob.print_List_FHDI();
	//}

	//------------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------Broadcast of jp_prob
	int List_rst_prob_size = 0;
	if (mynode == 0) {
		int temp_size = 0;
		for (int i = 0;i < i_nc;i++) {
			List_rst_prob.get_a_row_size(i, temp_size);
			List_rst_prob_size = List_rst_prob_size + temp_size;
		}
	} //end of mynode
	MPI_Bcast(&List_rst_prob_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//cout<<"Mynode: "<<mynode<<",The size of jp_prob is "<< List_rst_prob_size<<endl;
	double *List_rst_prob_broadcast = new double[List_rst_prob_size];
	//double *prob_temp = new double[temp_size];
	if (mynode == 0) {
		int temp_size_prob = 0;
		int counter_prob = 0;
		for (int i = 0;i < i_nc;i++) {
			List_rst_prob.get_a_row_size(i, temp_size_prob);
			double *prob_temp = new double[temp_size_prob];
			//vector <double> *prob_temp = new vector<double>(temp_size_prob);
			List_rst_prob.get_block(i, prob_temp);

			for (int k = 0;k < temp_size_prob;k++) {
				List_rst_prob_broadcast[counter_prob] = prob_temp[k];
				counter_prob++;
			}
			delete[] prob_temp;
		}//end of i_nc
	}
	MPI_Bcast(List_rst_prob_broadcast, List_rst_prob_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	////for (int z = 0;z < List_rst_prob_size;z++) {
	////	cout << "z: " << z<<", prob: "<< List_rst_prob_broadcast [z]<< endl;
	////}

	////cout<<"Mynode "<<mynode<<" Jp_debug1-----"<<endl;
	////----------------------------------------Broadcast of jp_name
	int List_rst_name_size = 0;
	if (mynode == 0) {
		int temp_size2 = 0;
		for (int i = 0;i < i_nc;i++) {
			List_rst_name.get_a_row_size(i, temp_size2);
			List_rst_name_size = List_rst_name_size + temp_size2;
		}
	}
	MPI_Bcast(&List_rst_name_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	////cout << "Mynode " << mynode << " Jp_debug2-----" << endl;
	char *List_rst_name_broadcast_char = new char[List_rst_name_size*ncol];

	if (mynode == 0) {
		string *List_rst_name_broadcast_string = new string[List_rst_name_size];
		int temp_size_name = 0;
		int counter_name = 0;
		for (int i = 0;i < i_nc;i++) {
			List_rst_name.get_a_row_size(i, temp_size_name);
			string *name_temp = new string[temp_size_name];
			List_rst_name.get_block(i, name_temp);

			for (int k = 0;k < temp_size_name;k++) {
				List_rst_name_broadcast_string[counter_name] = name_temp[k];
				counter_name++;
			}
			delete[] name_temp;
		}//end of i_nc
		//cout << "I'm here--" << endl;
		//--------convert to char
		int counter_trash =0;
		//char *List_rst_name_broadcast_char = new char[List_rst_name_size*ncol];
		for (int i = 0; i < List_rst_name_size;i++) {
			for (int m = 0;m < ncol;m++) {
				List_rst_name_broadcast_char[counter_trash] = List_rst_name_broadcast_string[i][m];
				//cout << "counter_trash: " << counter_trash << ", i:"<<i<<", m:"<<m<<endl;
				counter_trash++;
			}
		}
		//cout << "I'm here as well" << endl;
		delete[] List_rst_name_broadcast_string;
	}// end of if mynode
	MPI_Bcast(List_rst_name_broadcast_char, (List_rst_name_size*ncol), MPI_CHAR, 0, MPI_COMM_WORLD);
	//cout << "Mynode " << mynode << " Jp_debug3-----" << endl;
	////--------------------------Re-assemble of jp_prob and jp_name in each slave
	////----jp_prob
	if (mynode != 0) {
		int counter0 = 0;
		std::vector<double> jp_prob_buffer_again;
		for (int i = 0;i < i_nc;i++) {
			jp_prob_buffer_again.clear();
			for (int k = 0; k < jp_prob_size_each;k++) {
				jp_prob_buffer_again.push_back(List_rst_prob_broadcast[counter0]);
				//cout << "mynode " << mynode << "counter0: " << counter0 << endl;
				counter0++;
			}
			List_rst_prob.put_block(i, jp_prob_buffer_again);
		}
		//cout << "Mynode " << mynode << " Jp_debug4-----" << endl;
	//	//----------------------------jp_name
		std::vector<std::string> List_rst_name_broadcast_string_temp;
		string jp_name_temp0;
		int counter1 = 0;
		//convert char array to string vector
		for (int k = 0;k < List_rst_name_size;k++) {
			jp_name_temp0 = "";
			for (int i = 0;i < ncol;i++) {
				jp_name_temp0 = jp_name_temp0 + List_rst_name_broadcast_char[counter1];
				//cout << "mynode " << mynode << "counter1: " << counter1 << endl;
				counter1++;
			}
			List_rst_name_broadcast_string_temp.push_back(jp_name_temp0);
		}
		//cout << "Mynode " << mynode << " Jp_debug5-----" << endl;
		//---------
		std::vector<std::string> jp_name_buffer_again;
		int counter2 = 0;
		for (int i = 0;i < i_nc;i++) {
			jp_name_buffer_again.clear();
			//cout << "jp_name_size_return: " << jp_name_size_return << endl;
			for (int k = 0;k < jp_name_size_each;k++) {
				jp_name_buffer_again.push_back(List_rst_name_broadcast_string_temp[counter2]);
				//cout<<"counter2: "<<counter2<<", i: "<<i<<", k: "<<k<<endl;
				counter2++;
			}
			List_rst_name.put_block(i, jp_name_buffer_again);
		}
		//cout << "Mynode " << mynode << " Jp_debug6-----" << endl;
	}//end of mynode

	//testout------
	//cout <<"Mynode: "<<mynode <<", List_rst_name:" << endl;;
	//List_rst_name.print_List_string_FHDI();
	//cout << "Mynode: "<<mynode<<", List_rst_prob:" << endl;
	//List_rst_prob.print_List_FHDI();
	//cout << "Yang Running time of Rep_cell_communication3 at node " << mynode << " is " << MPI_Wtime() - Rep_CellP_begin_communication << endl;
	//if (mynode == 0) { cout << "Rep_CellP successfully finished--------------------------"; }
	//-------------------------------------------------------------------------------------------------------------------------------
	delete[] cn0;
	//delete[] cn0_backup;
	delete[] s_cn0_temp;

	delete[] w_UserDefined;

	delete[] List_rst_prob_broadcast;
	delete[] List_rst_name_broadcast_char;

	//delete[]

	return;
}