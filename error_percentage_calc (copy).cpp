
#include "constants.h"
//#include "utils.h"
#include "vmc.h"



int main()
{


  std::ifstream ft;
  std::ifstream fr;
  std::ofstream fw;
 
  open2read(ft,"./tr_data/yd_test_data.txt");   
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

  
  float h_count_cor =0, sa_count_cor =0, d_count_cor =0, a_count_cor=0, sur_count_cor= 0, fe_count_cor =0;

  // count after test
  float *h_count_test =allocate <float>(NO);
  float *sa_count_test =allocate <float>(NO);
  float *d_count_test =allocate <float>(NO);
  float *a_count_test=allocate <float>(NO);
  float *sur_count_test= allocate <float>(NO);
  float *fe_count_test =allocate <float>(NO);

  // assign zero to all the counts 
  for ( size_t i=0; i < NO; i++)

    {

      h_count_test[ i] =0;
      sa_count_test [i] =0;

      d_count_test[ i] =0;
      a_count_test [i] =0;

      su_count_test[ i] =0;
      fe_count_test [i] =0;



    }

  
  int input=0; 
  
  open2write(fw, " err_percentage.txt");
  //******************************************************//

 
  //read all the data from different files, combine all and write back all the data into a single file
  while (!ft.eof() )
   
    {
   

      for (size_t i= 0; i<NO; i++)
	{
	  ft >>yd[i];
	  fr >>y[i] ;

	

	  if ( yd [i] == 1)
	    {
	      input = i+1;

	      else
		input = input;
	  
	    }

	}
      //===============//

      
      switch (input)
	{

	case 1:
	  //input data is for happiness
	  h_count++;

	      
	  for ( size_t i =0; i < NO; i++)

	    {

	      if( i == 0 && y[i] == 1 )

		h_count_test[i]++;

	     }// end of for loop in case 1
		
	    

	  break;

	case 2:

	  //input data is for sadness

	  sa_count++;

	      
	  for ( size_t i =0; i < NO; i++)

	    {

	      if( i == 0 && y[i] == 1 )

		sa_count_test++;

	     }// end of for loop in case 1
		

	  
	  break;
	  //===============================//
	case 3:

	  //input data is for disgust
	  d_count++;

	      
	  for ( size_t i =0; i < NO; i++)

	    {

	      if( i == 0 && y[i] == 1 )

		h_count_conf++;

	      else if ( i == 1 && y[i] == 1 )

		sa_count_conf++;

	      else if ( i == 2 && y[i] == 1 )

		d_count_cor++;

	      else if ( i == 3 && y[i] == 1 )

		a_count_conf++;

	      else if ( i == 4 && y[i] == 1 )

		su_count_conf++;
		    
	      else if ( i == 5 && y[i] == 1 )

		fe_count_conf++;

		  
	    }// end of for loop in case 1
		


	  
	  break;
	  //=============================================//

	case 4:

	  // input data is for anger

	  a_count++;

	      
	  for ( size_t i =0; i < NO; i++)

	    {

	      if( i == 0 && y[i] == 1 )

		h_count_conf++;

	      else if ( i == 1 && y[i] == 1 )

		sa_count_conf++;

	      else if ( i == 2 && y[i] == 1 )

		d_count_conf++;

	      else if ( i == 3 && y[i] == 1 )

		a_count_cor++;

	      else if ( i == 4 && y[i] == 1 )

		su_count_conf++;
		    
	      else if ( i == 5 && y[i] == 1 )

		fe_count_conf++;

		  
	    }// end of for loop in case 1
		
	  
	  break;

	  //==========================================//


	  
	case 5:

	  //input data is for surprise

	  su_count++;

	      
	  for ( size_t i =0; i < NO; i++)

	    {

	      if( i == 0 && y[i] == 1 )

		h_count_conf++;

	      else if ( i == 1 && y[i] == 1 )

		sa_count_conf++;

	      else if ( i == 2 && y[i] == 1 )

		d_count_conf++;

	      else if ( i == 3 && y[i] == 1 )

		a_count_conf++;

	      else if ( i == 4 && y[i] == 1 )

		su_count_cor++;
		    
	      else if ( i == 5 && y[i] == 1 )

		fe_count_conf++;

		  
	    }// end of for loop in case 1
		

	  
	  break;
	  //==========================================//
	case 6:

	  //input data is for fear

	  fe_count++;

	      
	  for ( size_t i =0; i < NO; i++)

	    {

	      if( i == 0 && y[i] == 1 )

		h_count_conf++;

	      else if ( i == 1 && y[i] == 1 )

		sa_count_conf++;

	      else if ( i == 2 && y[i] == 1 )

		d_count_conf++;

	      else if ( i == 3 && y[i] == 1 )

		a_count_conf++;

	      else if ( i == 4 && y[i] == 1 )

		su_count_conf++;
		    
	      else if ( i == 5 && y[i] == 1 )

		fe_count_cor++;

		  
	    }// end of for loop in case 1
		
	  
	  break;

	  //=======================================// 
	default :

	  //no input detected

	  break;
	      

	} // end of switch case


    } // end of for loop

 //===========================================================//
      
}


// std::cout<< fe_count<< "\t" <<std::endl;

// std::cout<< h_count_cor<< "\t" << sa_count_conf<< "\t" <<d_count_conf <<"\t" <<a_count_conf<<"\t"<<sur_count_conf<<"\t"<<fe_count_conf<<std::endl;

// std::cout<< h_count_conf<< "\t" << sa_count_conf<< "\t" <<d_count_conf <<"\t" <<a_count_conf<<"\t"<<sur_count_conf<<"\t"<<fe_count_cor<<std::endl;



std::cout<<"percentage= "<< h_count_conf/ sur_count <<"\t"<< sa_count_conf/ sur_count <<"\t"<< d_count_conf/ sur_count<<"\t"<< a_count_conf/ sur_count <<"\t"<< sur_count_cor/ sur_count <<"\t"<< fe_count_conf/ sur_count <<std::endl;
  


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

del (su_count_test);
del (fe_count_test);



ft.close();
fr.close();
fw.close();



return(0);
}
