void order_FHDI(int* i_original, const int n)
//Description ================================
// Order the array in ascending order
//
//INOUT   : int i_original(n) returned with the ordered cell numbers  [0, ...]   
//
//=============================================
{
	int* i_source[n]; 
	int* i_order[n]; 
	
	for(int i=0; i<n; i++) 
	{
		i_source[i] = i_original[i]; //backup
	}

	//-----------
	//leverage sorting library
	//-----------
	std::sort(i_source, i_source+n);
	int i_now = 0;
	
	i_order[0] = 1; //first cell location as default
	for(int i=0; i<n; i++)
	{
		i_now = i_source[i];
		//----------------
		//comparisons from the first entiry to now 
		//----------------
		for(int j=0; j<i; j++)
		{
			if(i_now == i_original[j])
			{
				i_order[i] = j+1; //Actual location
				i_original[j] = -12345678; //dummy value
				break; 
			}				
		}
	}
/*	
	//-----------
	//basic comparison with (n-1) entities
	//-----------
	int i_now = 0; 
	int i_left = 0; 
	int i_temp1 = 0; 
	int i_temp2 = 0; 
	for(int i_col=1; i_col<n;i_col++)
	{
		//----------------
		//comparisons from now all the way to the first entity
		//----------------
		for(int i=i_col; i>0; i--)
		{
			i_now = i_source[i]; 
			i_left = i_source[i-1];
			//---
			//swap two cells 
			//---
			if(i_now<i_left) 
			{
				i_temp1 = i_source[i-1]; //contents swap
				i_temp2 = i_source[i];
				i_source[i-1] = i_temp2; 
				i_source[i] = i_temp1; 
				
				i_temp1 = i_order[i-1]; //cell no. swap 
				i_temp2 = i_order[i];
				i_order[i-1] = i_temp2; 
				i_order[i] = i_temp1;
			}
			//---
			//no need to swap two cells 
			//---
			if(i_now>=i_left) 
			{
				break; //exist inner loop
			}
		}
	}
*/
	
	//---prep return
	Copy_iVector(i_order, n, i_original);
	
	delete[] i_source; 
	delete[] i_order; 
	
	return;
}
