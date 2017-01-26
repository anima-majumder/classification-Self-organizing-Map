#include<iostream>
#include "constants.h"

int main()
{
  
  double** clsCnt;
  double*  origClasCnt;
  const  size_t N=332; // total number of data in the file 
  const size_t dimData=1; 
  const size_t numClass=6;
  int temp =0;
  std::ifstream ft;
  std::ifstream fr;
   open2read(ft,"./tr_data/yd_test1d.txt");   
   open2read(fr,"./results.txt");

  // for(size_t i=0; i<N; i++){
  //   dataset[i] = new int[1+1]; 
  
   // } 


   clsCnt = new double* [numClass];
   for(size_t i=0; i<numClass; i++){
     
     clsCnt[i] = new double[numClass]; 
  
  }
    origClasCnt = new double[numClass];
    int c;
    int op;
    int count;
  for(size_t i=0; i<N; i++){
    //read response data to bla
    count =0;
    if (dimData ==6)
      {

	for (size_t i =0; i<dimData ; i++)
	  {
	    ft>>temp;

	    if(temp ==1)
	      count++;
	    
	  }
	ft>>c; 

	fr>>op;


    
      }

    else
      {

   ft>>c; 

    fr>>op;

      }
    origClasCnt[c-1]=origClasCnt[c-1]+1;
    clsCnt[c-1][op-1]++;
  }

  for(size_t t=0; t<numClass; t++){
      for(size_t v=0; v<numClass; v++){
	std::cout<<clsCnt[t][v]<<"\t";
      }
      std::cout<<std::endl;
    }

  for(size_t z=0; z<numClass; z++){
    for(size_t r=0; r<numClass; r++){
      clsCnt[z][r]=100*clsCnt[z][r]/origClasCnt[z];
      }
  }
  std::cout<<"====================================================="<<std::endl;
  
    for(size_t t=0; t<numClass; t++){
      for(size_t v=0; v<numClass; v++){
	std::cout<<clsCnt[t][v]<<"\t";
      }
      std::cout<<std::endl;
    }
    
      
  

  ft.close();
  fr.close();


  delete  origClasCnt;

  
 return 0;
 

}
