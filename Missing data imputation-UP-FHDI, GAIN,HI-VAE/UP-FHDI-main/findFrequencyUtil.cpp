void findFrequencyUtil(double* arr, int low, int high,
	vector<int>& freq)
{
	// If element at index low is equal to element  
	// at index high in the array 
	if (arr[low] == arr[high])
	{
		// increment the frequency of the element 
		// by count of elements between high and low 
		freq[arr[low]] += high - low + 1;
	}
	else
	{
		// Find mid and recurse for left and right  
		// subarray 
		int mid = (low + high) / 2;
		findFrequencyUtil(arr, low, mid, freq);
		findFrequencyUtil(arr, mid + 1, high, freq);
	}
}