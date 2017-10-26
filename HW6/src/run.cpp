#include "ising.h"
#include <string>
#include <fstream>

//==============================================================================
void savegrid(int**, int,  string);
//==============================================================================

int main(){
  int N = 50;
  double Nsq = N*N;
  int iter = N*N*N*N;
  int start = iter/2; // number of iterations to allow system to equilibriate
  double t = 0.1;
  double e,m,c,xi;
  ising vanilla(N,iter,start,t);

  string LowT = "LowT.csv";
  string HiT = "HiT.csv";
  string TEM = "TEM.csv";

  ofstream file;
  file.open(TEM.c_str(),ios::out);

  savegrid(vanilla.Grid(),N,LowT);

  for(t=0.1 ; t < 10; t+=0.1){
    cout << "T = " << t << endl;
    vanilla.SetTemp(t);
    e = vanilla.Eavg/Nsq;
    m = vanilla.Mabs/Nsq;
    c = (vanilla.Esq/Nsq - e*vanilla.Eavg)/(t*t);
    xi = (vanilla.Msq/Nsq - m*vanilla.Mabs)/t;
    file << t << "," << e << "," << m << "," << c << "," << xi << endl;
  }
  savegrid(vanilla.Grid(),N,HiT);
  

  return 0;
}

//==============================================================================
void savegrid(int **grid, int N, string filename){
  ofstream file;
  file.open(filename.c_str(),ios::out);
  for(int i=0; i < N; i++){
    for(int j=0; j< N; j++){
      if(j < N-1){
        file << grid[i][j] << ",";
      }
      else{
        file << grid[i][j] << endl;
      }
    }
  }
  file.close();
}

