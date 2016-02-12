#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace std;


string revComp(const string& seq);
string getCanonical(const string& seq);
string getKmer(const string& sequence, int posi, int k);
