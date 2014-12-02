
/************************************************************************************************//**
 *  		
 * 	file: 		test.cpp
 * 	contents:   	See file.h
 * 
 ****************************************************************************************************/

//#include "file.h"
#include <vector>
#include <Eigen/Core>
#include <memory>
#include <iostream>
#include <boost/multi_array.hpp>
#include <chrono>

using namespace std;
using namespace Eigen;

/**
 *	Description of main method	
 */
int main (int argc, char * argv[])
{
   int len = 10000;
  
   boost::multi_array<double,2>*  arr = new boost::multi_array<double,2>(boost::extents[len][len]);
   
   cout << " size of multiarray " << sizeof(arr) << endl;
  


   auto start_time = chrono::steady_clock::now();
   /********************* Start Timing********************/

   double temp = 0.0;
   for (int i = 0; i < len; ++i)
   {
      temp += (*arr)[0][i];
   }
   
   /********************* End Timing ********************/
   auto end_time = chrono::steady_clock::now();
   auto diff =  end_time - start_time;




   cout << " Execution time " << chrono::duration <double, milli> (diff).count() << " ms "  << endl;

   return 0;
}


