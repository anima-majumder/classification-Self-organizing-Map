#include "vmc.h"
#include "constants.h"
#include <cv.h>
#include <highgui.h>
#include <cxcore.h>
#include<cvaux.h>
#include <time.h>

int main()
{
  SOM s(NX,NY,NZ, NO, NI);


  uint nc = NX * NY * NZ;
  
  float*  x = allocate<float>(NI);
  float*  xn = allocate<float>(NI);
 
  float** w = allocate<float>(nc, NI);
 
 
  float* y= allocate<float>(NO);               //  output at each node 
  float* y_d=allocate<float>(NO);              //desired output at

  float* eta =allocate<float>(3);

  float* Er = allocate<float>(NO);
  float ErNorm = 0.0f;
  
  float teta=0.0;
 
  //================================//
  int* index=allocate<int>(3);
  int count=0;
  float sigma=0;
  zero(y, NO);//initialized to zero

#ifndef DEBUG
  //  s.readweight();
  //  display("should not display");
#endif

  
  //  s.getWeight(w);

  std::ifstream fr;

  //=====read desired output from file ====================//
  std::ifstream f_yd;


  //======================================================//
  std::ofstream fo;
  open2write(fo,"trData.txt");

  uint   fcount = 0;
  // do{

    open2read(fr, "./tr_data_normalized/au_training_data.txt");//read a data for test
    open2read(f_yd, "tr_data_normalized/yd_tr_data.txt");
   
    
    //===start execution time ====================//

    clock_t begin, end;
    double time_spent;

    begin = clock();

   // for (int itr=0; itr<NOD; itr++)// number of lines in the file
   while ( !fr.eof())// number of lines in the file
      {

	zero(y, NO);//initialized to zero
	

	read(fr, x, NI);

	read(f_yd, y_d, NO);
// #ifdef DEBUG
// 	//this is for case of 2D data
// 	x[0]=1.0;
// 	x[1]=1.0;

// 	y_d[0]=5.0;
// 	y_d[0]=5.0;
    
// #endif

	copy(x, xn, NI);
	s.selectWinner(xn);
	s.get_winner_index(index);
	//====================================//
 

	//neighbourhood width
	  sigma = SIGMA_INIT * pow( (SIGMA_FINAL / SIGMA_INIT), (float(count)/NOD) );
	  s.Spread() = sigma;

      
	//	s.Spread() = 1.50;//works if we keep a fixed sigma ( say sigma=5)
  

		s.Spread() = SIGMA;//works if we keep a fixed sigma ( say sigma=5)
	   
  
	teta = ETA1_INIT * pow( (ETA1_FINAL / ETA1_INIT), (float(count)/(10*NOD)));
	// teta = ETA1_INIT * exp( -(float(count)/NOD));
 
	set(eta, teta, 3);

	s.setRate(eta);
	s.Network_output(x, y);//check-- y is calculated using this function
	s.Parameters_training( y,  x,  y_d);

	/*	sub(y_d, y, Er, NO);
	ErNorm = norm2(Er, NO);

	if( count > 1)
	  fo<<ErNorm<<std::endl;
	*/

	count++;

      }
    //=========check the execution time ================//

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    std::cout<<"execution time =" <<time_spent <<std::endl; 

    
     
    fr.close();
    f_yd.close();

    fcount++;
    // }while(fcount < 1);
  //=========================================//

#ifndef DEBUG
  s.saveSOM();

#endif

 
  //do proper termination of the code

  fo.close();
  // fr.close();
  // f_yd.close();
  del(Er);
  del (x);
  del (xn);
  del (w,nc);
  del (index);
 
  del (y);
  del (eta);
  del (y_d);
  return(0);
}
