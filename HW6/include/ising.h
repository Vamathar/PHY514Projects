#ifndef ISING_H
#define ISING_H

#include <iostream>
#include <root/TRandom3.h>
#include <functional>
#include <cmath>

using namespace std;

class ising{
  public:
    // initialized with grid size and temperature
    ising(int,int,int,double);
    ~ising();
    double Temp();
    void SetTemp(double);
    int **Grid();
    int En();
    int Mag();

    // Equilibrium Algorithms (plural for future projects)
    void EqSSF(int);

    double Mabs;
    double Msq;
    double M;
    double E;
    double Eavg;
    double Esq;

  private:

    int SumNN(int,int);
    void prob(double);
    int NSame(int,int);
    void Avg(int,int);
    
    int n; // size of grid
    int start; // number of iterations to allow system to equilibriate
    int ntimes; // iterations of averaging
    int eq; // number of iterations to come to equilibrium
    bool equil; // keeps track of whether or not system has been allowed to equilibriate
    double T; // temperature

    int **grid; // grid of 1s and -1s representing spin state of system
    double de[5] = {-8.,-4.,0.,4.,8.};
    double pe[5];

};

#endif
