#include "vmc.h"
#include "constants.h"
#include <cv.h>
#include <highgui.h>
#include <cxcore.h>
#include<cvaux.h>


int main()
{
  SOM s(NX,NY,NZ, NO, NI);

  float*  x = allocate<float>(NI);
  float*  xn = allocate<float>(NI);
  uint nc = NX * NY * NZ;
  float** w = allocate<float>(nc, NI);
  //==================================//
 
  float* y= allocate<float>(NO);               //  output at each node 
  float* y_d=allocate<float>(NO);              //desired output at

  float* eta =allocate<float>(3);

  float* Er = allocate<float>(NO);

  float* res_count = allocate<float>(NO); //to count the total number of
  float* percentage = allocate<float>(NO);				      //output expressions that give 1
  float ErNorm = 0.0f;
  
  float teta=0.0;
  //each node


  //================================//
  int* index=allocate<int>(3);
  int count=0;
  float sigma=0;
  zero(y, NO);//initialized to zero
  zero (res_count, NO);
  zero (percentage, NO);
#ifndef DEBUG
  //s.readweight();
  s.readSOM();// read all the three files 
  //  display("should not display");
#endif

  std::ifstream fr;

  std::ofstream fy;


  std::ifstream fyd;

  std::ofstream ferr;
  //======================================================//
 
 

  uint   fcount = 0;
  //  do{

  open2read(fr, "./tr_data_normalized/au_test_data.txt");//read a data for test

  open2write(fy, "results.txt");//to write the results data
  open2write(ferr, "test_err_norm.txt");

  open2read(fyd, "./tr_data_normalized/yd_test_data.txt");
  while (!fr.eof())// number of lines in the file
    {

      read(fr, x, NI);
      read(fyd, y_d, NO);

      copy(x, xn, NI);
      s.selectWinner(xn);
      s.get_winner_index(index);//not needed 
      //====================================//
      // get the count using the obtained winner index.
      // get the y[count] and pass it to the function
     

      //neighbourhood width
      //  sigma = SIGMA_INIT * pow( (SIGMA_FINAL / SIGMA_INIT), (float(count)/NOD) );
      //       sigma = SIGMA_INIT * exp( (-float(count)/NOD) );

      

      s.Spread() = SIGMA;//works if we keep a fixed sigma ( say sigma=5)

	   

      s.Network_output(x, y);//check-- y is calculated using this function


      // write the results data into a file

      print(fy, y, NO);

     
      //===========this part is added to get the result count==============================//
      for ( int i=0; i<NO; i++)
	{

	  if ( y[i]==1)

	    res_count[i]++;
	}
      //==================end of result count=======================//

      //read(fyd, y_d, NO);
      sub(y_d, y, Er, NO);
      ErNorm = norm2(Er, NO);

      ferr<<ErNorm<<std::endl;


      count++;

    }


  display (res_count, NO);
  float total_count =0;

  for (int i=0; i<NO; i++)
    total_count= total_count + res_count[i];

  std::cout << "total count =\t"<< total_count<<std::endl;

  //===============calculate the percentage here=======================//
  display ("percentage are = \t");
  for (int i=0; i<NO; i++)
    {
      percentage[i] = (res_count[i] / total_count )*100.0;

      std::cout <<percentage[i] <<"\t";
    }



  
  // display (percentage, NO);
  //=========percentage calculation is done here==================//

  
  fr.close();
  

  fcount++;
  //  }while(fcount < 3);
  //=========================================//


 
  //do proper termination of the code
  fy.close();
  fyd.close();
  ferr.close();
  del(Er);
  del (x);
  del (xn);
  del (w,nc);
  del (index);
  del (y);
  del (eta);
  del (percentage);
  return(0);
}
