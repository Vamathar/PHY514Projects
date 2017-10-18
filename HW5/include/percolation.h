#ifndef PERCOLATION_H
#define PERCOLATION_H
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <root/TRandom3.h>

using namespace std;

class percolate{
    public:
	// default contructor
	percolate(int, double);
	// destructor
	~percolate();
	double P();
	int N();
	void SetP(double); // sets p to a new value and remakes tfgrid
	void Cluster(); // clusters tfgrid and makes clusterlabels
	vector<int> Percolating(); // returns proper labels of percolating clusters
	vector<vector<int> > Grid(); // returns clustered and relabeled tfgrid

    private:
	bool occupied(int,int); // returns whether or not a site in tfgrid is occupied
	int relabel(int); // returns proper label of input label

	int n; // size of grid
	double p; // probability
	vector< vector<double> > pgrid; // grid filled with probability values
	vector< vector<int> > tfgrid; // grid filled with 1s and 0s
	vector<int> clusterlabels; // keep track of proper cluster labels
	bool clustered;
	TRandom3 rand;

};


#endif
