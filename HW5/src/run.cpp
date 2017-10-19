#include "percolation.h"
#include <string>
#include <fstream>

//==============================================================================
void savevector(vector<int>&, string);
void savevector(vector<double>&, string);
void savegrid(vector< vector< int > >&, string);
//==============================================================================
int main(){

    // Part 1: make a cluster grid for visual pleasure
    int n = 50;
    double p = 0.58;
    percolate a(n,p);
    vector<int> clust;
    clust = a.Cluster();
    string clustname = "ClustNames.tsv";
    savevector(clust,clustname);

    vector< vector< int > > grid;
    string gridname = "Grid.tsv";
    grid = a.Grid();
    savegrid(grid,gridname);

    // Part 2: make P vs p coordinates
    vector<double> pvec;
    n = 200;
    double pstart = 0.53;
    double pend = 0.63;
    int repeat = 10; // number of times to average
    string xname = "Px200.tsv";
    string yname = "Py200.tsv";
    for(p = pstart; p<pend; p+=0.01){
        pvec.push_back(p);
    }
    vector<double> bigPvec(pvec.size(),0.);

    for(int i = 0; i < repeat; i++){
        cout << "i: " << i << endl;
        percolate a(n,p);
        for(int j = 0; j<pvec.size(); j++){
            a.SetP(pvec[j]);
            a.Cluster();
            if(a.Percolating().size() != 0){
                bigPvec[j] += 1.;
            }
        }
    }
    for(int i = 0; i<bigPvec.size(); i++){
        bigPvec[i] = bigPvec[i]/double(repeat);
    }
    savevector(pvec,xname);
    savevector(bigPvec,yname);

    bigPvec.clear();
    pvec.clear();
    n = 100;
    pstart = 0.45;
    pend = 0.65;
    repeat = 20; // number of times to average
    xname = "Px100.tsv";
    yname = "Py100.tsv";
    for(p = pstart; p<pend; p+=0.01){
        pvec.push_back(p);
    }
    bigPvec.assign(pvec.size(),0.);

    for(int i = 0; i < repeat; i++){
        cout << "i: " << i << endl;
        percolate a(n,p);
        for(int j = 0; j<pvec.size(); j++){
            a.SetP(pvec[j]);
            a.Cluster();
            if(a.Percolating().size() != 0){
                bigPvec[j] += 1.;
            }
        }
    }
    for(int i = 0; i<bigPvec.size(); i++){
        bigPvec[i] = bigPvec[i]/double(repeat);
    }
    savevector(pvec,xname);
    savevector(bigPvec,yname);

    bigPvec.clear();
    pvec.clear();
    n = 50;
    pstart = 0.40;
    pend = 0.70;
    repeat = 40; // number of times to average
    xname = "Px50.tsv";
    yname = "Py50.tsv";
    for(p = pstart; p<pend; p+=0.01){
        pvec.push_back(p);
    }
    bigPvec.assign(pvec.size(),0.);

    for(int i = 0; i < repeat; i++){
        cout << "i: " << i << endl;
        percolate a(n,p);
        for(int j = 0; j<pvec.size(); j++){
            a.SetP(pvec[j]);
            a.Cluster();
            if(a.Percolating().size() != 0){
                bigPvec[j] += 1.;
            }
        }
    }
    for(int i = 0; i<bigPvec.size(); i++){
        bigPvec[i] = bigPvec[i]/double(repeat);
    }
    savevector(pvec,xname);
    savevector(bigPvec,yname);
    return 0;
}
//==============================================================================
void savevector(vector<int> &vec, string filename){
    ofstream file;
    file.open(filename.c_str(),ios::out);
    for(int i=0; i<vec.size(); i++){
        file << vec[i] << "\t";
    }
    file << endl;
    file.close();
}
//==============================================================================
void savevector(vector<double> &vec, string filename){
    ofstream file;
    file.open(filename.c_str(),ios::out);
    for(int i=0; i<vec.size(); i++){
        file << vec[i] << "\t";
    }
    file << endl;
    file.close();
}
//==============================================================================
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

