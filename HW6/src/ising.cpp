#include "ising.h"
#include <ctime>

//==============================================================================
// Constructor
ising::ising(int NN, int EEQQ, int strt, double TT){
  
  // Assign private variables
  n = NN;
  eq = EEQQ;
  T = TT;
  start = strt;
  ntimes = eq-start;
  equil = false;

  // Establish random number generator
  TRandom3 rand;
  rand.SetSeed(time(NULL));

  // Initialize grid
  grid = new int*[n];
  for(int i = 0; i < n; i++){
    grid[i] = new int[n];
    for(int j = 0; j < n; j++){
      grid[i][j] = (rand.Rndm() < 0.5)? 1 : -1;
    }
  }


  // Equilibriate system
  prob(T);
  EqSSF(start);
}

//==============================================================================
// Destructor
ising::~ising(){
  for(int i = 0; i < n; i++){
    delete[] grid[i];
  }
  delete[] grid;
}

//==============================================================================
// calculates possible probabilities with given temperature
void ising::prob(double t){
  pe[0] = 1.;
  pe[1] = 1.;
  pe[2] = 1.;
  pe[3] = exp(-de[3]/t);
  pe[4] = exp(-de[4]/t);
}

//==============================================================================
// returns the number of nearest neighbors with the same spin
int ising::NSame(int i, int j){
  
  int value = 0;

  // nearest neighbor indices
  int im1 = i-1;
  int ip1 = i+1;
  int jm1 = j-1;
  int jp1 = j+1;

  // periodic boundary conditions
  if(im1 < 0){ im1 = n-1; }
  if(ip1 > n-1){ ip1 = 0; }
  if(jm1 < 0){ jm1 = n-1; }
  if(jp1 > n-1){ jp1 = 0; }

  if(grid[i][j] == grid[im1][j]){ value++; }
  if(grid[i][j] == grid[ip1][j]){ value++; }
  if(grid[i][j] == grid[i][jm1]){ value++; }
  if(grid[i][j] == grid[i][jp1]){ value++; }
  
  return value;
}

//==============================================================================
// Returns the temperature
double ising::Temp(){ return T; }

//==============================================================================
// Set a new temperature
void ising::SetTemp(double TT){ 
  TRandom3 rand;
  rand.SetSeed(time(NULL));
  T = TT;
  prob(T);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      grid[i][j] = (rand.Rndm() < 0.5)? 1 : -1;
    }
  }
  EqSSF(start);
}

//==============================================================================
// Return pointer to grid
int **ising::Grid(){ return grid; }

//==============================================================================
// Returns nearest neighbor sum for Hamiltonian
int ising::SumNN(int i, int j){
  // nearest neighbor indices
  int im1 = i-1;
  int ip1 = i+1;
  int jm1 = j-1;
  int jp1 = j+1;

  // periodic boundary conditions
  if(im1 < 0){ im1 = n-1; }
  if(ip1 > n-1){ ip1 = 0; }
  if(jm1 < 0){ jm1 = n-1; }
  if(jp1 > n-1){ jp1 = 0; }
  
  int sum = 0;
  sum += grid[i][j]*grid[ip1][j];
  sum += grid[i][j]*grid[im1][j];
  sum += grid[i][j]*grid[i][jp1];
  sum += grid[i][j]*grid[i][jm1];

  return sum;
}

//==============================================================================
// Returns energy of system 
int ising::En(){
  int energy = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      energy -= SumNN(i,j);
    }
  }
  return energy;
}

//==============================================================================
// Returns magnetization of system
int ising::Mag(){
  int mag = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      mag+=grid[i][j];
    }
  }
  return mag;
}


//==============================================================================
// Average values
void ising::Avg(int Mnew, int Enew){
  Mabs += abs(Mnew);
  Msq += Mnew*Mnew;
  Eavg += Enew;
  Esq += Enew*Enew;

}
//==============================================================================
// Brings system to equilibrium using single spin flip Metropolis MC.
void ising::EqSSF(int start){


  // Establish random number generator
  TRandom3 rand;
  rand.SetSeed(time(NULL));

  // Initialize averaging quantities
  Mabs = 0.;
  Msq = 0.;
  E = 0.;
  Esq = 0.;
  bool doavg = false;
  int Mnew, Enew, count;

  int i,j,old,diff,same;
  double pflip;
  for(int ieq = 0; ieq < eq; ieq++){
    // Start averaging after enough iterations have passed
    if(ieq == start){
      count = 1;
      M = Mag();
      Mabs = abs(M);
      Msq = M*M;
      E = En();
      Eavg = E;
      Esq = E*E;
      doavg = true;
    }

    // get random point and save old value for later
    i = int(rand.Rndm()*n);
    j = int(rand.Rndm()*n);
    old = grid[i][j];
    same = NSame(i,j);
    pflip = pe[same];

    if(pflip == 1.0){ // flip spin if energy decreases
      grid[i][j]*=-1;

    }
    else{ // flip spin with probability exp(-de/T)
      if(rand.Rndm() < pflip){
        grid[i][j]*=-1;
      }
    }
    if(doavg){
      count++;
      diff = grid[i][j]-old;
      if(diff != 0){
        Mnew = M+diff;
        Enew = E+de[same];
        Avg(Mnew,Enew);
        M = Mnew;
        E = Enew;
      }
      else{
        Avg(M,E);
      }
    }
  }
  Eavg = Eavg/count;
  Esq = Esq/count;
  Mabs = Mabs/count;
  Msq = Msq/count;
  equil = true;
  return;
}
