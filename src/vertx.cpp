#include "vertx.h"
#include <stdio.h>

using namespace Eigen;
using namespace std;

vertx::vertx()
{}

vertx::vertx(double u_val, double phi_val, double beta_val, double freq[], int len)
    : Grid(boost::extents[len][len][2][2][2][2]), u(u_val), Dd(cos(phi_val/2.0)), phi_dot(phi_val), beta(beta_val)
{
//    cout << " Dd " << Dd << endl;
//    cout << " u " << u << endl;
    M4cd H;
    H <<    0.0,    0.0,        0.0,    0.0,      // SQDJJ
            0.0,    u/2.0,      -Dd,    0.0,
            0.0,    -Dd,        u/2.0,  0.0,
            0.0,    0.0,        0.0,    0.0;

    SelfAdjointEigenSolver<Matrix4cd> es;
    es.compute(H);
    eval=es.eigenvalues();
	evec=es.eigenvectors();
	rho.setZero();

	Z = 0.0;

    for(int i = 0; i < eval.size(); i++)
    {
        Z += exp(-beta * eval(i));
        rho(i) = exp(-beta * eval(i));
    }
    rho*=1.0/Z;

//    cout << " Z " << Z << endl;
//    cout << " eval " << endl << eval << endl;
//    cout << " rho " << endl << rho << endl;

    nambu[0].setZero();
    nambu[1].setZero();
	nambu[0](0,1) = 1.0; nambu[0](2,3) = 1.0; //  Nambu basis 00 10 01 11
	nambu[1](0,2) = 1.0; nambu[1](1,3) = -1.0;  // !!!!! Define in the right way!

	eta=1e-8;

	initU();


	// ----- PREPROCESS OPERATOR MATRIX ELEMENTS

    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
            for(int k = 0; k < 2; ++k){
                MEl[i][j][k] = matEl(i,j,nambu[k]);
                MElC[i][j][k] = matEl(i,j,nambu[k].adjoint());
            }
        }
    }

	//p = atli(0.5, u, 0.0, phi/PI, 0.0);

//	nambu[0](0,1) = 1.0; nambu[0](2,3) = 1.0; // Dot Basis
//	nambu[1](2,0) = 1.0; nambu[1](3,1) = -1.0;



    initGrid(freq, len);
}

dcomplex vertx::matEl(int i, int j, M4cd op)
{
    return evec.adjoint().row(i)*op*evec.col(j);
}

dcomplex vertx::G1p(double w, int m, int n)
{
    //M4cd a[2] = {nambu[m], nambu[n].adjoint()};
    dcomplex val = 0.0;
    for(int i = 0; i < eval.size(); ++i){
        for(int j = 0; j < eval.size(); ++j){
            val +=  //rho(i)*(  (  matEl(i,j,a[0]) * matEl(j,i,a[1]) ) / (eval(i)-eval(j)+I*w) + (  matEl(i,j,a[1]) * matEl(j,i,a[0]) ) / (eval(j)-eval(i)+I*w) );
                    (rho(i) + rho(j)) *  MEl[i][j][m] * MElC[j][i][n] / (eval(i)-eval(j)+I*w);
        }
    }
//    if(std::isnan(real(val)) == 1 || std::isnan(imag(val)) == 1){ cout << " G1PC NAN w1 m n"
//        << w << " "<< m << " " << n << endl; getchar(); } // Use continuous extension on problematic points!!
    return val;
}

M2cd vertx::G1p(double w)
{
    M2cd Gm;
    Gm <<    G1p(w,0,0),   G1p(w,0,1),
            G1p(w,1,0),   G1p(w,1,1);
    return Gm;
}

M2cd vertx::invG0_ATLI(double w){
    M2cd invG;
    invG <<		I*w,	Dd,
                Dd,			I*w;
    return invG;
}

M2cd vertx::Sig_ATLI(double w){
    return invG0_ATLI(w) - G1p(w).inverse();
}

dcomplex vertx::G2pc(double w1, double w2, double w2p, int m, int n, int o, int p, dcomplex PVal[3][3][3][4][4][4][4])
{
//    M4cd a[4] = {nambu[m], nambu[n], nambu[p].adjoint(), nambu[o].adjoint()};
    //double w[3] = {w1, w2, -w2p};
    dcomplex val = 0.0;

    for(int i = 0; i < eval.size(); ++i){
        for(int j = 0; j < eval.size(); ++j){
            for(int k = 0; k < eval.size(); ++k){
                for(int l = 0; l < eval.size(); ++l){
                    val +=  PVal[0][1][2][i][j][k][l] * MEl[i][j][m]  * MEl[j][k][n]  * MElC[k][l][p] * MElC[l][i][o]
                        -   PVal[0][2][1][i][j][k][l] * MEl[i][j][m]  * MElC[j][k][p] * MEl[k][l][n]  * MElC[l][i][o]
                        -   PVal[1][0][2][i][j][k][l] * MEl[i][j][n]  * MEl[j][k][m]  * MElC[k][l][p] * MElC[l][i][o]
                        +   PVal[1][2][0][i][j][k][l] * MEl[i][j][n]  * MElC[j][k][p] * MEl[k][l][m]  * MElC[l][i][o]
                        -   PVal[2][1][0][i][j][k][l] * MElC[i][j][p] * MEl[j][k][n]  * MEl[k][l][m]  * MElC[l][i][o]
                        +   PVal[2][0][1][i][j][k][l] * MElC[i][j][p] * MEl[j][k][m]  * MEl[k][l][n]  * MElC[l][i][o];

//                    val +=  phi(w[0],w[1],w[2],i,j,k,l) * MEl[i][j][m] * MEl[j][k][n] * MElC[k][l][p] * MElC[l][i][o]
//                        -   phi(w[0],w[2],w[1],i,j,k,l) * MEl[i][j][m] * MElC[j][k][p] * MEl[k][l][n] * MElC[l][i][o]
//                        -   phi(w[1],w[0],w[2],i,j,k,l) * MEl[i][j][n] * MEl[j][k][m] * MElC[k][l][p] * MElC[l][i][o]
//                        +   phi(w[1],w[2],w[0],i,j,k,l) * MEl[i][j][n] * MElC[j][k][p] * MEl[k][l][m] * MElC[l][i][o]
//                        -   phi(w[2],w[1],w[0],i,j,k,l) * MElC[i][j][p] * MEl[j][k][n] * MEl[k][l][m] * MElC[l][i][o]
//                        +   phi(w[2],w[0],w[1],i,j,k,l) * MElC[i][j][p] * MEl[j][k][m] * MEl[k][l][n] * MElC[l][i][o];

//                    val +=  phi(w[0],w[1],w[2],i,j,k,l) * matEl(i,j,a[0]) * matEl(j,k,a[1]) * matEl(k,l,a[2]) * matEl(l,i,a[3])
//                        -   phi(w[0],w[2],w[1],i,j,k,l) * matEl(i,j,a[0]) * matEl(j,k,a[2]) * matEl(k,l,a[1]) * matEl(l,i,a[3])
//                        -   phi(w[1],w[0],w[2],i,j,k,l) * matEl(i,j,a[1]) * matEl(j,k,a[0]) * matEl(k,l,a[2]) * matEl(l,i,a[3])
//                        +   phi(w[1],w[2],w[0],i,j,k,l) * matEl(i,j,a[1]) * matEl(j,k,a[2]) * matEl(k,l,a[0]) * matEl(l,i,a[3])
//                        -   phi(w[2],w[1],w[0],i,j,k,l) * matEl(i,j,a[2]) * matEl(j,k,a[1]) * matEl(k,l,a[0]) * matEl(l,i,a[3])
//                        +   phi(w[2],w[0],w[1],i,j,k,l) * matEl(i,j,a[2]) * matEl(j,k,a[0]) * matEl(k,l,a[1]) * matEl(l,i,a[3]);
                }
            }
        }
    }

    if(w2 == w2p) { val -= beta*G1p(w1,m,o)*G1p(w2,n,p);}
    if(w1 == w2p) { val += beta*G1p(w1,m,p)*G1p(w2,n,o);}

//    if(std::isnan(real(val)) == 1 || std::isnan(imag(val)) == 1){ cout << " G2PC NAN w1 w2 w2p m n o p "
//            << w1 << " " << w2 << " " << w2p << " " << m << " " << n << " " << o << " " << p << " "
//            << endl; getchar(); } // Use continuous extension on problematic points!!
    return val;
}

dcomplex vertx::G2pc_w2w3(double w1, double w2, double w2p, int m, int n, int o, int p)
{
    M4cd a[4] = {nambu[m], nambu[n], nambu[p].adjoint(), nambu[o].adjoint()};
    double w[3] = {w1, w2, -w2};
    dcomplex val = 0.0;
    if ( w2 != w2p ) {return val;}

    for(int i = 0; i < eval.size(); i++){
        for(int j = 0; j < eval.size(); j++){
            for(int k = 0; k < eval.size(); k++){
                for(int l = 0; l < eval.size(); l++){
                    val +=  -1.0/(I*w[2]+eval(k)-eval(l))*rho(j)/(I*w[0]+eval(i)-eval(j)) * kdel(eval(j),eval(l))
                        * matEl(i,j,a[0]) * matEl(j,k,a[1]) * matEl(k,l,a[2]) * matEl(l,i,a[3])
                        + 1.0/(I*w[1]+eval(k)-eval(l))*rho(j)/(I*w[0]+eval(i)-eval(j)) * kdel(eval(j),eval(l))
                        * matEl(i,j,a[0]) * matEl(j,k,a[2]) * matEl(k,l,a[1]) * matEl(l,i,a[3])
                        - 1.0/(I*w[0]+eval(k)-eval(l))*rho(i)/(I*w[2]+eval(j)-eval(k)) * kdel(eval(i),eval(k))
                        * matEl(i,j,a[1]) * matEl(j,k,a[2]) * matEl(k,l,a[0]) * matEl(l,i,a[3])
                        + 1.0/(I*w[0]+eval(k)-eval(l))*rho(i)/(I*w[1]+eval(j)-eval(k)) * kdel(eval(i),eval(k))
                        * matEl(i,j,a[2]) * matEl(j,k,a[1]) * matEl(k,l,a[0]) * matEl(l,i,a[3]);
                }
            }
        }
    }
    val -= G1p(w1,m,o)*G1p(w2,n,p);// Substract disconnected diagram without exchange
    return beta*val;
}

dcomplex vertx::vrtx(double w1, double w2, double w2p, int m, int n, int o, int p)
{
    M2cd G1Inv = G1p(w1).inverse();
    M2cd G2Inv = G1p(w2).inverse();
    M2cd G3Inv = G1p(w1+w2-w2p).inverse();
    M2cd G4Inv = G1p(w2p).inverse();

//    if( std::isnan( (G1Inv*G2Inv*G3Inv*G4Inv).norm() ) == 1 ){ cout << " G1Inv NAN " << G1p(w1) << " rho " << rho << endl; getchar();}

    // ----- PREPROCESS PHI (since independent of i,j,k,l in the following -----------
    dcomplex PVal[3][3][3][4][4][4][4];// WATCH OUT, Amound of EV PUT IN EXPLICITELY here
    double w[3] = {w1, w2, -w2p};
    for(int i = 0; i < eval.size(); ++i){
        for(int j = 0; j < eval.size(); ++j){
            for(int k = 0; k < eval.size(); ++k){
                for(int l = 0; l < eval.size(); ++l){
                    PVal[0][1][2][i][j][k][l] = phi(w[0],w[1],w[2],i,j,k,l);
                    PVal[0][2][1][i][j][k][l] = phi(w[0],w[2],w[1],i,j,k,l);
                    PVal[1][0][2][i][j][k][l] = phi(w[1],w[0],w[2],i,j,k,l);
                    PVal[1][2][0][i][j][k][l] = phi(w[1],w[2],w[0],i,j,k,l);
                    PVal[2][1][0][i][j][k][l] = phi(w[2],w[1],w[0],i,j,k,l);
                    PVal[2][0][1][i][j][k][l] = phi(w[2],w[0],w[1],i,j,k,l);
                }
            }
        }
    }
    // --------------------------------
    dcomplex val = 0.0;
    for(int i = 0; i < 2; ++i){
        for(int j = 0; j < 2; ++j){
            for(int k = 0; k < 2; ++k){
                for(int l = 0; l < 2; ++l){
                    val += G1Inv(m,i) * G2Inv(n,j) * G2pc(w1,w2,w2p,i,j,k,l,PVal) * G3Inv(k,o) * G4Inv(l,p);
                }
            }
        }
    }
    if(std::isnan(real(val)) == 1 || std::isnan(imag(val)) == 1){ cout << " VERTEX NAN w1 w2 w2p m n o p "
            << w1 << " " << w2 << " " << w2p << " " << m << " " << n << " " << o << " " << p << " "
            << endl; getchar(); } // Use continuous extension on problematic points!!

    return val;
}

dcomplex vertx::vrtx_w2w3(double w1, double w2, double w2p, int m, int n, int o, int p)
{
    M2cd G1Inv = G1p(w1).inverse();
    M2cd G2Inv = G1p(w2).inverse();
    M2cd G3Inv = G1p(w1+w2-w2p).inverse();
    M2cd G4Inv = G1p(w2p).inverse();

    dcomplex val = 0.0;
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            for(int k = 0; k < 2; k++){
                for(int l = 0; l < 2; l++){
                    val +=  G1Inv(m,i) * G2Inv(n,j) * G2pc_w2w3(w1,w2,w2p,i,j,k,l) * G3Inv(k,o) * G4Inv(l,p);
                }
            }
        }
    }
    return val;
}

void vertx::initGrid(double f[], int len)
{
    cout << " start init " << endl;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < len; ++i){
        cout << " outer loop i " << i << endl;
        for(int j = i; j < len; ++j){
            Grid[i][j][0][0][0][0] = vrtx(f[i],f[j],f[j],0,0,0,0); //-----1
            Grid[j][i][0][0][0][0] = Grid[i][j][0][0][0][0];

            Grid[i][j][0][0][0][1] = vrtx(f[i],f[j],f[j],0,0,0,1); //-----2
            Grid[j][i][0][0][1][0] = Grid[i][j][0][0][0][1];

            Grid[i][j][0][0][1][0] = vrtx(f[i],f[j],f[j],0,0,1,0); //-----3
            Grid[j][i][0][0][0][1] = Grid[i][j][0][0][1][0];

            Grid[i][j][0][0][1][1] = vrtx(f[i],f[j],f[j],0,0,1,1); //------4
            Grid[j][i][0][0][1][1] = Grid[i][j][0][0][1][1];

            Grid[i][j][0][1][0][0] = Grid[i][j][0][0][0][1];
            Grid[j][i][1][0][0][0] = Grid[i][j][0][1][0][0];

            Grid[i][j][0][1][0][1] = vrtx(f[i],f[j],f[j],0,1,0,1); //------5
            Grid[j][i][1][0][1][0] = Grid[i][j][0][1][0][1];

            Grid[i][j][0][1][1][0] = vrtx(f[i],f[j],f[j],0,1,1,0); //------6
            Grid[j][i][1][0][0][1] = Grid[i][j][0][1][1][0];

            Grid[i][j][1][0][0][0] = Grid[i][j][0][0][1][0];
            Grid[j][i][0][1][0][0] = Grid[i][j][1][0][0][0];

            Grid[i][j][0][1][1][1] = -conj(Grid[i][j][1][0][0][0]);
            Grid[j][i][1][0][1][1] = Grid[i][j][0][1][1][1];

            Grid[i][j][1][0][0][1] = conj(Grid[i][j][0][1][1][0]);
            Grid[j][i][0][1][1][0] = Grid[i][j][1][0][0][1];

            Grid[i][j][1][0][1][0] = conj(Grid[i][j][0][1][0][1]);
            Grid[j][i][0][1][0][1] = Grid[i][j][1][0][1][0];

            Grid[i][j][1][0][1][1] = -conj(Grid[i][j][0][1][0][0]);
            Grid[j][i][0][1][1][1] = Grid[i][j][1][0][1][1];

            Grid[i][j][1][1][0][0] = conj(Grid[i][j][0][0][1][1]);
            Grid[j][i][1][1][0][0] = Grid[i][j][1][1][0][0];

            Grid[i][j][1][1][0][1] = -conj(Grid[i][j][0][0][1][0]);
            Grid[j][i][1][1][1][0] = Grid[i][j][1][1][0][1];

            Grid[i][j][1][1][1][0] = -conj(Grid[i][j][0][0][0][1]);
            Grid[j][i][1][1][0][1] = Grid[i][j][1][1][1][0];

            Grid[i][j][1][1][1][1] = conj(Grid[i][j][0][0][0][0]);
            Grid[j][i][1][1][1][1] = Grid[i][j][1][1][1][1];

//            for(int m = 0; m < 2; m++){
//                for(int n = 0; n < 2; n++){
//                    for(int o = 0; o < 2; o++){
//                        for(int p = 0; p < 2; p++){
//                            cout << " ijmnop " << i << " " << j << " " << m << " " << n << " " << o << " " << p << " " << Grid[i][j][m][n][o][p] << endl;
//                        }
//                    }
//                }
//            }
//            getchar();
        }
    }
    cout << " end init " << endl;
}

dcomplex vertx::Gvrtx(int w1, int w2, int m, int n, int o, int p)
{
//    cout << w1 << " " << w2 << " " << m << " " << n << " " << o << " " << p << endl;
    return Grid[w1][w2][m][n][o][p];
}   // vertex calculated from grid

dcomplex vertx::phi(double w1, double w2, double w3, int i, int j, int k, int l) // contains only parts without delta function, also T = 0 !!
{
    dcomplex phi_val;

    if( ( w2 != -w3 || kdel(eval(j),eval(l)) != 1.0 ) && ( w1 != -w2 || kdel(eval(i),eval(k)) != 1.0 ) )
    {
        phi_val =   1.0/(I*w3+eval(k)-eval(l)) * (
                        1.0/(I*(w2+w3)+eval(j)-eval(l)) * (
                            (rho(i)+rho(j)) / (I*w1+eval(i)-eval(j)) - (rho(i)+rho(l)) / (I*(w1+w2+w3)+eval(i)-eval(l))
                        )
                        - 1.0/(I*w2+eval(j)-eval(k)) * (
                            (rho(i)+rho(j)) / (I*w1+eval(i)-eval(j)) - (rho(i)-rho(k)) / (I*(w1+w2)+eval(i)-eval(k))
                        )
                    );
    }
    else if( ( w2 == -w3 && kdel(eval(j),eval(l)) == 1.0 ) && ( w1 != -w2 || kdel(eval(i),eval(k)) != 1.0 ) )
    {
        phi_val =   1.0/(I*w3+eval(k)-eval(l)) * (
                        (rho(i) + rho(j))/(I*w1+eval(i)-eval(j))/(I*w1+eval(i)-eval(j)) - beta*rho(j)/(I*w1+eval(i)-eval(j))
                        - 1.0/(I*w2+eval(j)-eval(k)) * (
                            (rho(i)+rho(j)) / (I*w1+eval(i)-eval(j)) - (rho(i)-rho(k)) / (I*(w1+w2)+eval(i)-eval(k))
                        )
                    );
    }
    else if( ( w2 != -w3 || kdel(eval(j),eval(l)) != 1.0 ) && ( w1 == -w2 && kdel(eval(i),eval(k)) == 1.0 ) )
    {
        phi_val =   1.0/(I*w3+eval(k)-eval(l)) * (
                        1.0/(I*(w2+w3)+eval(j)-eval(l)) * (
                            (rho(i)+rho(j)) / (I*w1+eval(i)-eval(j)) - (rho(i)+rho(l)) / (I*(w1+w2+w3)+eval(i)-eval(l))
                        )
                        - 1.0/(I*w2+eval(j)-eval(k)) * (
                            (rho(i)+rho(j)) / (I*w1+eval(i)-eval(j)) + beta*rho(i)
                        )
                    );
    }
    else if( ( w2 == -w3 && kdel(eval(j),eval(l)) == 1.0 ) && ( w1 == -w2 && kdel(eval(i),eval(k)) == 1.0 ) )
    {
        phi_val =   1.0/(I*w3+eval(k)-eval(l)) * (
                        (rho(i) + rho(j))/(I*w1+eval(i)-eval(j))/(I*w1+eval(i)-eval(j)) - beta*rho(j)/(I*w1+eval(i)-eval(j))
                        - 1.0/(I*w2+eval(j)-eval(k)) * (
                            (rho(i)+rho(j)) / (I*w1+eval(i)-eval(j)) + beta*rho(i)
                        )
                    );
    }

    return phi_val;
}

double vertx::kdel(double a, double b)
{
    if (a == b){ return 1.0;}
    return 0.0;
}

void	// Assign Initial Values to U
vertx::initU(){

	for (int m = 0; m < 2; m++){ // Initialisiere zu 0
		for (int n = 0; n < 2; n++){
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					U[m][n][k][l] = 0.0;
				}
			}
		}
	}
//
	U[0][1][0][1] = -u; // Setze Anfangsbedingungen
	U[0][1][1][0] = u;
	U[1][0][0][1] = u;
	U[1][0][1][0] = -u;
}
