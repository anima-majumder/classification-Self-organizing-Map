
#include "constants.h"

#include "vmc.h"



int main()
{


  std::ifstream ft;
  std::ifstream fr;
  std::ofstream fw;
 
  open2read(ft,"./tr_data_normalized/yd_test_data.txt");   
  open2read(fr,"./results.txt");

  float *yd = NULL;
  yd= allocate <float>(NO);
  float *y = NULL;
  y= allocate <float>(NO);
  float temp=0.0;
  // results for happiness
  //1=happiness [ 1 0 0 0 0 0]
  // 2= sadness, [ 0 1 0 0 0 0]
  //3= disgust,  [ 0 0 1 0 0 0]
  //4= anger  [ 0 0 0 1 0 0]
  //, 5= surprise  [ 0 0 0 0 1 0]
  //and 6= fear  [ 0 0 0 0 0 1]


  //total count
  float h_count =0, sa_count =0, d_count =0, a_count=0, sur_count= 0, fe_count =0;

  
  

  // count after test
  float *h_count_test =allocate <float>(NO);
  float *sa_count_test =allocate <float>(NO);
  float *d_count_test =allocate <float>(NO);
  float *a_count_test =allocate <float>(NO);
  float *sur_count_test =allocate <float>(NO);
  float *fe_count_test =allocate <float>(NO);

  //each pointer loads the count of detected emotion
  // if the emotion is happy then h_count_test[0] will have correctly
  // recognition emotion count and other will have confused detection
  // count like wise for sadness sa_count_test[1] will have correctly
  // recognized emotion. 


  float* h_percent = NULL;
  float* sa_percent = NULL;
  float* d_percent = NULL;
  float* a_percent = NULL;
  float* sur_percent = NULL;
  float* fe_percent = NULL;

  h_percent = allocate <float> (NO);
  sa_percent = allocate <float> (NO);
  d_percent = allocate <float> (NO);
  a_percent = allocate <float> (NO);
  sur_percent = allocate <float> (NO);
  fe_percent = allocate <float> (NO);
  

  
  // assign zero to all the counts 
  for ( size_t i=0; i < NO; i++)

    {

      h_count_test[ i] =0;
      sa_count_test [i] =0;

      d_count_test[ i] =0;
      a_count_test [i] =0;

      sur_count_test[ i] =0;
      fe_count_test [i] =0;



    }

  
  int input=0; 
 
  int count =0; 
  open2write(fw, "err_percentage.txt");
  //******************************************************//

 
  //read all the data from different files, combine all and write back all the data into a single file
  while (!ft.eof() )
   
    {
   
      count =0;
      for (size_t i= 0; i<NO; i++)
	{
	  ft >>yd[i];
	  fr >>y[i];

	

	  if ( yd [i] == 1)
	    
	      input = i+1;
	    
	      else
		input = input;
	  }
      //===============//

      
      switch (input)
	{

	case 1:
	  //input data is for happiness
	  h_count++;

	      
	  for ( size_t i =0; i < NO; i++)

	    {

	      if( y[i] == 1 )
		{
		count++;
	     
		//	if(count ==1)
		  h_count_test[i]++;
		}

	    }// end of for loop in case 1
		
	    

	  break;
	  //========================================//
	case 2:

	  //input data is for sadness

	  sa_count++;

	      
	  for ( size_t i =0; i < NO; i++)

	    {

	      if( y[i] == 1 )
		{
		  count++;
	     
		  //	if(count ==1)

		  sa_count_test[i]++;
		}

	    }// end of for loop in case 1
		

	  
	  break;
	  //===============================//
	case 3:

	  //input data is for disgust
	  d_count++;

	      
	  for ( size_t i =0; i < NO; i++)

	    {

	      if(  y[i] == 1 )
		{
		  count++;
	     
		  //	  if(count ==1)
		  d_count_test[i]++;
		}
	    }// end of for loop in case 1
		


	  
	  break;
	  //=============================================//

	case 4:

	  // input data is for anger

	  a_count++;

	      
	  for ( size_t i =0; i < NO; i++)

	    {

	      if( y[i] == 1 )
		{
		  count++;
	     
		  //	if(count ==1)
		a_count_test[i]++;

	     
		}
		  
	    }// end of for loop in case 1
		
	  
	  break;

	  //==========================================//


	  
	case 5:

	  //input data is for surprise

	  sur_count++;

	      
	  for ( size_t i =0; i < NO; i++)

	    {

	      if(  y[i] == 1 )
		{
		  count++;
	     
		  //	if(count ==1)
		sur_count_test[i]++;
		}

		  
	    }// end of for loop in case 1
		

	  
	  break;
	  //==========================================//
	case 6:

	  //input data is for fear

	  fe_count++;

	      
	  for ( size_t i =0; i < NO; i++)

	    {

	      if(  y[i] == 1 )
		{
		  count++;
	     
		  //	if(count ==1)
		fe_count_test[i]++;

		}
		  
	    }// end of for loop in case 1
		
	  
	  break;

	  //=======================================// 
	default :

	  //no input detected
	  std::cout <<" no emotion level detected in this input vector " <<std::endl;
	  break;
	      

	} // end of switch case


      std::cout<<"count = " <<count <<std::endl;

  //===========================================================//
      
}// end of while loop 

// calculate the percentage here.  Basically it calculates the confusion
// matrix here 
for (size_t i =0; i< NO; i++)

  {
    h_percent[i] = 100 * (h_count_test[i] / h_count);
    sa_percent[i] = 100 * (sa_count_test[i] / sa_count);
    d_percent[i] = 100 * (d_count_test[i] / d_count);
    a_percent[i] = 100 * (a_count_test[i] / a_count);
    sur_percent[i] = 100 * (sur_count_test[i] / sur_count);
    fe_percent[i] = 100 * (fe_count_test[i] / fe_count);
    // every column of the file writes recognition percentages for that
    // emotion. say for sa_percentage .. row 0 : recognition result of
    // happiness data when sadness is input ..row tow is the recognition
    // result of sadness when sadness is input  
    // fw<< h_percent[i] <<"\t " <<sa_percent[i] <<"\t" <<d_percent[i] <<"\t" <<a_percent[i] <<"\t"
    //   <<sur_percent[i]<<"\t" <<fe_percent[i] <<std::endl;

  }




  for (size_t k =0; k<NO; k++)
    fw<<h_percent[k] <<"\t";
  fw<<std::endl;

  for (size_t k =0; k<NO; k++)
    fw<<sa_percent[k] <<"\t";
  fw<<std::endl;

  for (size_t k =0; k<NO; k++)
    fw<<d_percent[k] <<"\t";
  fw<<std::endl;

  for (size_t k =0; k<NO; k++)
    fw<<a_percent[k] <<"\t";
  fw<<std::endl;


  for (size_t k =0; k<NO; k++)
    fw<<sur_percent[k] <<"\t";
  fw<<std::endl;


  for (size_t k =0; k<NO; k++)
    fw<<fe_percent[k] <<"\t";
  fw<<std::endl;


 

// std::cout<<"percentage= "<< h_count_conf/ sur_count <<"\t"<< sa_count_conf/ sur_count <<"\t"<< d_count_conf/ sur_count<<"\t"<< a_count_conf/ sur_count <<"\t"<< sur_count_cor/ sur_count <<"\t"<< fe_count_conf/ sur_count <<std::endl;
  


//fw<< h_per<<"\t"<< sa_per<<"\t";
	    
//fw<<"\n";



 
//deallocate all the memories
//===============================//
del (yd);
del (y);

del (h_count_test);
del (sa_count_test);

del (d_count_test);
del (a_count_test);

del (sur_count_test);
del (fe_count_test);

del (h_percent);
del (sa_percent);
del (d_percent);
del (a_percent);
del (sur_percent);
del (fe_percent);


ft.close();
fr.close();
fw.close();



return(0);
}
