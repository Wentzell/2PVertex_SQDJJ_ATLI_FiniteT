#ifndef GAMMA_H
#define GAMMA_H

// INCLUDES
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <complex>
#include <iostream>
#include <omp.h>
#include "boost/multi_array.hpp"
#include <cassert>
#include "def.h"

//const int VERT_FREQ_COUNT = 2*100;

class vertx{
   public:

      double u;
      double U[2][2][2][2];
      double Dd;
      double beta;
      double Z;
      double phi_dot;

      vertx();
      //vertx(double u_val, double Dd_val, double beta);
      vertx(double u_val, double Dd_val, double beta, double freq[], int len);
      dcomplex G1p(double w, int i, int j);
      MatQN G1p(double w);
      MatQN invG0_ATLI(double w);
      MatQN Sig_ATLI(double w);
      dcomplex G2pc(double w1, double w2, double w2p, int m, int n, int o, int p, dcomplex PVal[3][3][3][4][4][4][4]); //PVal contains Phi Values for speedup
      dcomplex G2pc_w2w3(double w1, double w2, double w2p, int m, int n, int o, int p);
      dcomplex vrtx(double w1, double w2, double w2p, int m, int n, int o, int p);
      dcomplex vrtx_w2w3(double w1, double w2, double w2p, int m, int n, int o, int p);

      void initGrid(double f[], int len);

      typedef boost::multi_array<dcomplex, 6> vertex_type;
      vertex_type Grid;
      dcomplex Gvrtx(int w1, int w2, int m, int n, int o, int p);   // vertex calculated from grid

      dcomplex phi(double w1, double w2, double w3, int i, int j, int k, int l);

   private:
      M4cd evec;
      V4d eval;
      V4d rho;
      double eta;
      M4cd nambu[2];
      dcomplex MEl[4][4][2];      // WATCH OUT, Amound of EV PUT IN EXPLICITELY here
      dcomplex MElC[4][4][2];      // WATCH OUT, Amound of EV PUT IN EXPLICITELY here
      void initU();

      double kdel(double a, double b);
      dcomplex matEl(int i, int j, M4cd op);
};

#endif
