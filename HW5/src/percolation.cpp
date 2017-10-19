#include "percolation.h"
#include <ctime>
//==================================================================== 
percolate::percolate(int nin, double pin){
    n = nin;
    p = pin;
    clustered = false;
    rand.SetSeed(time(NULL));
    for(int i=0; i<n; i++){
        vector<double> pvec;
        vector<int> tfvec;
        for(int j=0; j<n; j++){
            double val = double(rand.Rndm());
            pvec.push_back(val);
            if(val < p){
                tfvec.push_back(1);
            }
            else{
                tfvec.push_back(0);
            }
        }
        pgrid.push_back(pvec);
        tfgrid.push_back(tfvec);
    }	
}
//==================================================================== 
percolate::~percolate(){}
//==================================================================== 
double percolate::P(){ return p; }
//==================================================================== 
int percolate::N(){ return n; }
//==================================================================== 
void percolate::SetP(double pin){
    p=pin;
    for(int i=0; i<pgrid.size(); i++){
        for(int j=0; j< pgrid[i].size(); j++){
            if(pgrid[i][j] < p){
                tfgrid[i][j] = 1;
            }
            else{
                tfgrid[i][j] = 0;
            }
        }
    }
    clustered = false;
}
//==================================================================== 
vector<int> percolate::Cluster(){
    if(clustered){ return clusterlabels; }
    int label = 0;
    for(int i=0; i < tfgrid.size(); i++){
        for(int j=0; j < tfgrid[i].size(); j++){
            if(occupied(i,j)){ // skip empty sites
                if(!occupied(i-1,j) and !occupied(i,j-1)){ // new cluster condition
                    label++;
                    clusterlabels.push_back(label);
                    tfgrid[i][j] = label;
                }
                else if(occupied(i-1,j) and !occupied(i,j-1)){
                    tfgrid[i][j] = tfgrid[i-1][j];
                }
                else if(!occupied(i-1,j) and occupied(i,j-1)){
                    tfgrid[i][j] = tfgrid[i][j-1];
                }
                else{
                    if(tfgrid[i-1][j] == tfgrid[i][j-1]){
                        tfgrid[i][j] = tfgrid[i-1][j];
                    }
                    else if(tfgrid[i-1][j] < tfgrid[i][j-1]){
                        tfgrid[i][j] = tfgrid[i-1][j];
                        for(int k=0; k<clusterlabels.size(); k++){
                            if(clusterlabels[k] == tfgrid[i][j-1]){
                                clusterlabels[k] = tfgrid[i][j];
                            }
                        }
                    }
                    else if(tfgrid[i][j-1]<tfgrid[i-1][j]){
                        tfgrid[i][j] = tfgrid[i][j-1];
                        for(int k=0; k<clusterlabels.size(); k++){
                            if(clusterlabels[k] == tfgrid[i-1][j]){
                                clusterlabels[k] = tfgrid[i][j];
                            }
                        }
                    }
                    else{ cout << "Error" << endl;}
                }
            }
        }
    }

    for(int l=0; l < clusterlabels.size(); l++){
        if(clusterlabels[l] != l+1){
            for(int m = l+1; m < clusterlabels.size(); m++){
                if(clusterlabels[m] == l+1){clusterlabels[m]=clusterlabels[l]; }
            }
        }
    }
    clustered = true;
    return clusterlabels;
}
//==================================================================== 
bool percolate::occupied(int i, int j){
    if( i < 0 ){ return false; }
    if( j < 0 ){ return false; }
    if( i >= tfgrid.size() ){ return false; }
    if( j >= tfgrid[0].size() ){ return false; } 
    if( tfgrid[i][j] ){ return true; }
    return false;
}
//==================================================================== 
int percolate::relabel(int x){
    if(x <= 0 or x > clusterlabels.size()){ return 0; }
    return clusterlabels[x-1];
}
//==================================================================== 
vector<int> percolate::Percolating(){

    vector<int> top(tfgrid[0]);
    vector<int> bottom(tfgrid[tfgrid.size()-1]);
    vector<int> intersection;
    
    for(int i=0; i<top.size(); i++){
        top[i] = relabel(top[i]);
        bottom[i] = relabel(bottom[i]);
    }

    sort(top.begin(),top.end());
    sort(bottom.begin(),bottom.end());
    
    int ele = top[0];
    while(ele == 0){
        top.erase(top.begin());
        ele = top[0];
    }

    ele = bottom[0];
    while(ele == 0){
        bottom.erase(bottom.begin());
        ele = bottom[0];
    }

    vector<int> toplabels;
    vector<int> bottomlabels;
    
    toplabels.push_back(top[0]);
    bottomlabels.push_back(bottom[0]);
    
    for(int i=0; i<top.size(); i++){
        if(toplabels.back() != top[i]){
            toplabels.push_back(top[i]);
        }
    }

    for(int i=0; i<bottom.size(); i++){
        if(bottomlabels.back() != bottom[i]){
            bottomlabels.push_back(bottom[i]);
        }
    }

    set_intersection(toplabels.begin(),toplabels.end(),\
            bottomlabels.begin(),bottomlabels.end(),\
            back_inserter(intersection));
    return intersection;
}
//==================================================================== 
vector< vector< int > > percolate::Grid(){
   if(!clustered){ Cluster(); }
   for(int i=0; i<tfgrid.size(); i++){
        for(int j=0; j<tfgrid[i].size(); j++){
            tfgrid[i][j] = relabel(tfgrid[i][j]);
        }
    }
    return tfgrid;
}

