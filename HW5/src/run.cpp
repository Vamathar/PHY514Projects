#include "percolation.h"
#include <string>
#include <fstream>

void savevector(vector<int>&, string);
void savegrid(vector< vector< int > >&, string);
int main(){
    int n = 12;
    double p = 0.58;
    percolate a(n,p);
    a.Cluster();

    vector<int> perc;
    string percname = "PercNames.tsv";
    perc = a.Percolating();
    savevector(perc,percname);

    vector< vector< int > > grid;
    string gridname = "Grid.tsv";
    grid = a.Grid();
    savegrid(grid,gridname);

    return 0;

}


void savevector(vector<int> &vec, string filename){
    ofstream file;
    file.open(filename.c_str(),ios::out);
    for(int i=0; i<vec.size(); i++){
	file << vec[i] << "\t";
    }
    file << endl;
    file.close();
}

void savegrid(vector< vector< int > > &grid, string filename){
    ofstream file;
    file.open(filename.c_str(),ios::out);
    for(int i=0; i<grid.size(); i++){
	for(int j=0; j<grid[i].size(); j++){
	    file << grid[i][j] << "\t";
	}
	file << endl;
    }
}

