#include <iostream>
#include "vertx.h"
#include "def.h"
#include "../../Tools/Plot2D/plot2d.h"

using namespace std;

int main()
{
    double u = PI;//1.0;
    double beta = 100.0;
    double e = 0.0; double phi = 0.0*PI; double B=0.0; double phi_l = phi / 2.0;
    //dcomplex Dd = 0.5*exp(I*phi_l) + 0.5*exp(-I*phi_l); // 0 corresponds to SIAM
    double Dd = cos(phi/2.0);

    double w1=1 , w2=3, w2p=3;
    double w1p = w1+w2-w2p;


// ---- TEST INIT --------
    const int pos_freq_count = 100;
    const int freq_count = 2*pos_freq_count;
    double f[freq_count];

    for (int i = -pos_freq_count; i < pos_freq_count; ++i)
    {
        f[i + pos_freq_count] =    PI/beta*(2*i+1);//pow(2.0,i) ; //
    }

    vertx g(u, phi, beta, f, freq_count);

//    g.initGrid(f, freq_count);


//// ------------ PLOT VERTEX --------------
//    dcomplex *val0000 = new dcomplex[freq_count * freq_count];
//    dcomplex *val0001 = new dcomplex[freq_count * freq_count];
//    dcomplex *val0010 = new dcomplex[freq_count * freq_count];
//    dcomplex *val0011 = new dcomplex[freq_count * freq_count];
//    dcomplex *val0101 = new dcomplex[freq_count * freq_count];
//    dcomplex *val0110 = new dcomplex[freq_count * freq_count];
//////////
//      cout << "assign" << endl;
//    for (int i = 0; i < freq_count; ++i)
//    {
//        for (int j = 0; j < freq_count; ++j)
//        {
//            val0000[i + freq_count * j] = g.Gvrtx(i,j, 0, 0, 0, 0) - g.vrtx_w2w3(f[i],f[j],f[j],0, 0, 0, 0); //+ g.U[0][0][0][0]
//            val0001[i + freq_count * j] = g.Gvrtx(i,j, 0, 0, 0, 1) - g.vrtx_w2w3(f[i],f[j],f[j],0, 0, 0, 1); //+ g.U[0][0][0][1]
//            val0010[i + freq_count * j] = g.Gvrtx(i,j, 0, 0, 1, 0) - g.vrtx_w2w3(f[i],f[j],f[j],0, 0, 1, 0); //+ g.U[0][0][1][0]
//            val0011[i + freq_count * j] = g.Gvrtx(i,j, 0, 0, 1, 1) - g.vrtx_w2w3(f[i],f[j],f[j],0, 0, 1, 1); //+ g.U[0][0][1][1]
//            val0101[i + freq_count * j] = g.Gvrtx(i,j, 0, 1, 0, 1) - g.vrtx_w2w3(f[i],f[j],f[j],0, 1, 0, 1); //+ g.U[0][1][0][1]
//            val0110[i + freq_count * j] = g.Gvrtx(i,j, 0, 1, 1, 0) - g.vrtx_w2w3(f[i],f[j],f[j],0, 1, 1, 0); //+ g.U[0][1][1][0]
////            cout << " vertx 0 1 0 1 at w1 w2 " << i << " " << j << " " << val[i + freq_count * j] << endl;
//        }
//    }
////
//    plot(f, f, val0000, freq_count, freq_count, "val0000");
//    plot(f, f, val0001, freq_count, freq_count, "val0001");
//    plot(f, f, val0010, freq_count, freq_count, "val0010");
//    plot(f, f, val0011, freq_count, freq_count, "val0011");
//    plot(f, f, val0101, freq_count, freq_count, "val0101");
//    plot(f, f, val0110, freq_count, freq_count, "val0110");

    plot_file_imre("val0000");
    plot_file_imre("val0001");
    plot_file_imre("val0010");
    plot_file_imre("val0011");
    plot_file_imre("val0101");
    plot_file_imre("val0110");


//// ------------ PLOT VERTEX W2W3-------------
//    dcomplex *val0000 = new dcomplex[freq_count * freq_count];
//    dcomplex *val0001 = new dcomplex[freq_count * freq_count];
//    dcomplex *val0010 = new dcomplex[freq_count * freq_count];
//    dcomplex *val0011 = new dcomplex[freq_count * freq_count];
//    dcomplex *val0101 = new dcomplex[freq_count * freq_count];
//    dcomplex *val0110 = new dcomplex[freq_count * freq_count];
//////
////      cout << "assign" << endl;
//    for (int i = 0; i < freq_count; ++i)
//    {
//        cout << i << endl;
//        for (int j = 0; j < freq_count; ++j)
//        {
//            val0000[i + freq_count * j] = 2*PI/beta*g.vrtx_w2w3(f[i],f[j],f[j],0,0,0,0);
//            val0001[i + freq_count * j] = 2*PI/beta*g.vrtx_w2w3(f[i],f[j],f[j],0,0,0,1);
//            val0010[i + freq_count * j] = 2*PI/beta*g.vrtx_w2w3(f[i],f[j],f[j],0,0,1,0);
//            val0011[i + freq_count * j] = 2*PI/beta*g.vrtx_w2w3(f[i],f[j],f[j],0,0,1,1);
//            val0101[i + freq_count * j] = 2*PI/beta*g.vrtx_w2w3(f[i],f[j],f[j],0,1,0,1);
//            val0110[i + freq_count * j] = 2*PI/beta*g.vrtx_w2w3(f[i],f[j],f[j],0,1,1,0);
////            cout << " vertx 0 1 0 1 at w1 w2 " << i << " " << j << " " << val[i + freq_count * j] << endl;
//        }
//    }
////
//    plot(f, f, val0000, freq_count, freq_count, "val0000");
//    plot(f, f, val0001, freq_count, freq_count, "val0001");
//    plot(f, f, val0010, freq_count, freq_count, "val0010");
//    plot(f, f, val0011, freq_count, freq_count, "val0011");
//    plot(f, f, val0101, freq_count, freq_count, "val0101");
//    plot(f, f, val0110, freq_count, freq_count, "val0110");


//// ------------ 1D PLOT VERTEX --------------
//    dcomplex val0000_1D[freq_count];
//    dcomplex val0001_1D[freq_count];
//    dcomplex val0010_1D[freq_count];
//    dcomplex val0011_1D[freq_count];
//    dcomplex val0101_1D[freq_count];
//    dcomplex val0110_1D[freq_count];
//
//    int x = 60+pos_freq_count;
//
//    for (int i = 0; i < freq_count; ++i)
//    {
//            val0000_1D[i] = g.Gvrtx(x,i, 0, 0, 0, 0);
//            val0001_1D[i] = g.Gvrtx(x,i, 0, 0, 0, 1);
//            val0010_1D[i] = g.Gvrtx(x,i, 0, 0, 1, 0);
//            val0011_1D[i] = g.Gvrtx(x,i, 0, 0, 1, 1);
//            val0101_1D[i] = g.Gvrtx(x,i, 0, 1, 0, 1);
//            val0110_1D[i] = g.Gvrtx(x,i, 0, 1, 1, 0);
//    }
////
//    plot(f, val0000_1D, freq_count, "1D_val0000");
//    plot(f, val0001_1D, freq_count, "1D_val0001");
//    plot(f, val0010_1D, freq_count, "1D_val0010");
//    plot(f, val0011_1D, freq_count, "1D_val0011");
//    plot(f, val0101_1D, freq_count, "1D_val0101");
//    plot(f, val0110_1D, freq_count, "1D_val0110");

//    cout << " Sig(0) (0,1) " << g.Sig_ATLI(PI/beta)(0,1) << endl;

//    for(double w = 1e+4; w > 1e-2; w*=1e-1)
//    {
//        //cout << g.p.getSigATLI(w) - g.sig_EOM(w) << endl;
//        cout << g.p.getSigATLI(w) << endl << endl;
//        cout << g.sig_EOM(w) << endl << endl << endl;
//    }

    //ofstream data; string ffname;

//    ffname = "data/vertx.dat"; //sprintf(ffname, "data/%s/dat%s%s", par[i], par[i], fname);
//    data.open(ffname.c_str(), ios::out | ios::trunc);
//
//    for(int i = 1; i < 100; i++)
//    {
//        double u_val = u / 100.0 * i;
//        g = vertx(u_val, phi);
//        char temp[1000];
//        sprintf(temp, "%.5e %.5e\n", u_val, real(g.vrtx(0.0, 0.0, 0.0, 0, 1, 0, 1)/u_val));  // Writes Frequency
//        data << temp;
//    }
//    data.close();

//
//    cout << "g.vrtx(w1,w2,w2p, 0,1,0,1) " << g.vrtx(w1,w2,w2p, 0,1,0,1)<< endl;
//    cout << "g.vrtx(w1,w2,w2p, 0,1,1,1) " << g.vrtx(w1,w2,w2p, 0,1,1,1) << endl;
//    cout << "g.vrtx(w1,w2,w2p, 1,1,1,1) " << g.vrtx(w1,w2,w2p, 1,1,1,1) << endl;
//    cout << "g.vrtx(w1,w2,w2p, 0,0,1,1) " << g.vrtx(w1,w2,w2p, 0,0,1,1) << endl;
//    cout << "g.vrtx(w1,w2,w2p, 0,0,0,0) " << g.vrtx(w1,w2,w2p, 0,0,0,0) << endl;
//    cout << "g.vrtx(w1,w2,w2p, 0,0,0,1) " << g.vrtx(w1,w2,w2p, 0,0,0,1) << endl;

//
//    double gamUUSiamD22P = PI*u*u/2.0/(w1*w1*w2p*w2p)*(w1*w1+u*u/4.0)*(w2p*w2p+u*u/4.0);
//
//    cout << " gamUUSiamD22P " << gamUUSiamD22P << endl;

    // EQUATION OF MOTION CHECK


}
