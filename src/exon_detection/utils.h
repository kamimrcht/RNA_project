#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "correctionGraph.h"

#ifndef UTILS
#define UTILS

using namespace std;


string revComp(const string& seq);
string getCanonical(const string& seq);
string getKmer(const string& sequence, int posi, int k);
string getSequenceInWindow(const string & readSequence, uint w, uint k, uint indexWin);
string getSequenceInConsecutiveWindows(const string & readSequence, uint w, uint k, uint firstIndexWin, uint lastIndexWin);
uint absolute(int a);

#endif
