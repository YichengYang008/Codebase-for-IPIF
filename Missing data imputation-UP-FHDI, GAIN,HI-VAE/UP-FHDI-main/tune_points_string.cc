void tune_points_string(int& startpoint, int& endpoint, string* i_source, int mynode, int totalnodes,int n) {
	if (mynode != 0 && mynode != 1 && mynode != (totalnodes - 1)) {
		if ((i_source[startpoint]).compare(i_source[startpoint - 1])==0) {
			for (int i = startpoint + 1; i < n;i++) {
				if ((i_source[startpoint]).compare(i_source[i])!=0) {
					startpoint = i;
					break;
				}
				if (i == (n - 1)) {
					startpoint = i + 1;
				}
			}
		}

		if ((i_source[endpoint - 1]).compare(i_source[endpoint])==0) {
			for (int i = endpoint;i < n;i++) {
				//cout << "i1:" << i << endl;
				if ((i_source[endpoint - 1]).compare(i_source[i]) !=0) {
					//cout << "i2:" << i << endl;
					endpoint = i;
					break;
				}
				if (i == (n - 1)) {
					endpoint = i + 1;
				}
			}
		}
	}

	if (mynode == 1) {
		//cout << "I'm here" << endl;
		//cout << "i_source[endpoint-1]: " << i_source[endpoint - 1] << ", i_source[endpoint]: " << i_source[endpoint] << endl;
		if ((i_source[endpoint - 1]).compare(i_source[endpoint])==0) {
			for (int i = endpoint;i < n;i++) {
				//cout << "i1:" << i << endl;
				if ((i_source[endpoint - 1]).compare(i_source[i])!=0) {
					//cout << "i2:" << i << endl;
					endpoint = i;
					break;
				}
				if (i == (n - 1)) {
					endpoint = i + 1;
				}
			}
		}
	}

	if (mynode == (totalnodes - 1)) {
		if ((i_source[startpoint]).compare(i_source[startpoint - 1])==0) {
			for (int i = startpoint + 1; i < n;i++) {
				if ((i_source[startpoint]).compare(i_source[i])!=0) {
					startpoint = i;
					break;
				}
				if (i == (n - 1)) {
					startpoint = i + 1;
				}
			}
		}
	}
}