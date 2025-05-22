#include <vector>
#include <map>

vector<double> vote(map< vector<double>, int >& freqMap) {
	vector<double> ret;
	int maxCount = 0; // temp for only iteration
	for (map<vector<double>, int>::iterator it = freqMap.begin(); it != freqMap.end(); ++it) {
		if (it->second > maxCount) {
			ret = it->first;
			maxCount = it->second;
		}
	}
	return ret;
}