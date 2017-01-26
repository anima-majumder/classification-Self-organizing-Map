#include"SOM.h"

//===================================================================//
SOM::SOM()
{

  s_neuron_Nx = 0;  
  s_neuron_Ny = 0;
  s_neuron_Nz = 0;
  s_op_dimn_NO = 0;
  s_ip_dimn_NI = 0;

  s_spread=0.0;
  
  s_w = NULL;
  s_b = NULL;
  s_A = NULL;

  s_delta_b = NULL; 
  s_delta_w = NULL; 
  s_delta_A = NULL;
  a_del_x = NULL;
  yj = NULL;
 
  s_y = NULL;
  s_bmu_w_index = NULL;
  s_eta = NULL;
  
  s_del_xw = NULL;

 
  s_y_err = NULL;
  s_temp = NULL;
  sigma_dvt = NULL;
  
  s_test_sigma = NULL;
  srand(time(NULL));


 
  return ;
}

//===================================================================//
void SOM::createBuffers()
{

  uint NO = s_op_dimn_NO;
  uint NI = s_ip_dimn_NI;

  uint nc = s_neuron_Nx * s_neuron_Ny * s_neuron_Nz;


  s_w = allocate<float>(nc , NI);
  s_b = allocate<float>(nc, NO);
  s_A = allocate<float>(nc , (NO*NI) );

  s_delta_A = allocate<float>(nc , (NO*NI) );
  s_delta_b = allocate<float>(nc,NO);
  s_delta_w = allocate<float>(nc,NI);
  a_del_x = allocate<float>(nc, NO);
  yj = allocate<float>(nc, NO);


  s_y = allocate<float>(NO);
  s_bmu_w_index = allocate<int>(3);
  s_eta = allocate<float>(3);

 
  s_del_xw = allocate<float>(NI);
  


  s_y_err = allocate<float>(NO);
  s_temp = allocate<float>(NO);
  sigma_dvt = allocate <float>(1);
  s_test_sigma = allocate <float>(NO);
  

  return;

}
//===================================================================
SOM::SOM(const uint& Nx, const uint& Ny, const uint& Nz, const uint& NO,const uint& NI)
{

  srand(time(NULL));

  s_neuron_Nx = Nx;  
  s_neuron_Ny = Ny;
  s_neuron_Nz = Nz;

  s_op_dimn_NO = NO;
  s_ip_dimn_NI = NI;

  uint nc = Nx*Ny*Nz;

  createBuffers();

  
  for(size_t i = 0;i<nc;i++)
    {
      randomize(s_A[i], (NO*NI), float(0.05) );
      
      randomize(s_b[i], NO, float(0.05) );
      randomize (s_w[i], NI, float (0.001)); // weight is initialized
      // here 17th june 2013
    }

   
  
  set(s_bmu_w_index,-1, 3); // set the index of bmu to initial value -1
  // (which is actually not possible)

}
//===================================================================
SOM::~SOM()
{

  uint nc = s_neuron_Nx * s_neuron_Ny * s_neuron_Nz;


  del(s_w, nc);
  del(s_b, nc);
  del(s_A, nc);

  
  del(s_delta_b, nc);
  del(s_delta_w, nc);
  del(s_delta_A, nc);
  del (yj, nc);
  del (a_del_x, nc);

  del(s_bmu_w_index);

  del (s_y);
  del(s_eta);

  del(s_del_xw);
 

  del(s_y_err);
  del(s_temp);
  del(s_test_sigma);

#ifdef DEBUG
  display("SOM destructed ");
#endif


}
//===================================================================
void SOM::getIndex(const uint& count, uint& x, uint& y, uint& z)
{
  //   routine to get the lattice index from Parameter index count.
  //     count = (x*Ny + y)*Nz + z
  // z = count%Nz
  //     let count1 = (count - z)/Nz
  // then y = count1 %Ny
  //     and x  = (count1 - y)/Ny

  uint count1;

  z = count% s_neuron_Nz;
  count1 = (count - z)/s_neuron_Nz;

  y = count1 % s_neuron_Ny;

  x = (count1-y)/s_neuron_Ny;

  return;

}
//=========================================//
void SOM::getWeight(float** const& w)const{

  uint nc = s_neuron_Nx * s_neuron_Ny * s_neuron_Nz;

  for(uint n=0; n < nc; n++){
    copy(s_w[n], w[n], s_ip_dimn_NI);
  }

}
//=========================//


void SOM::get_winner_index(int *winner_index)
{



  for (uint n=0; n < 3; n++)
    {
      winner_index[n]=  s_bmu_w_index[n];
    }

}



//===================================================================
void SOM::selectWinner(const float* const& inp_vect)
{

  uint NI = s_ip_dimn_NI;
  uint NC = s_neuron_Nx * s_neuron_Ny * s_neuron_Nz;

  uint i=0;

  uint x = 0;
  uint y = 0;
  uint z = 0;

  float minD=0;
  float dist = 0;
  float sum=0;

  uint count = 0;

  set(s_bmu_w_index, -1, 3);

  for(count = 0;count < NC; count++)
    {

      sum = 0;

      for(i = 0; i < NI; i++)
	sum += pow((s_w[count][i] - inp_vect[i]), float(2.0) );

      dist= sqrt(sum);


      if(count > 0)
	{
          if(dist < minD)
	    {
	      minD = dist;
	      getIndex(count, x, y, z);
	      s_bmu_w_index[0] = x;
	      s_bmu_w_index[1] = y;
	      s_bmu_w_index[2] = z;
	    }
	}

      else
	{
	  minD = dist;
	  s_bmu_w_index[0] = 0;
	  s_bmu_w_index[1] = 0;
	  s_bmu_w_index[2] = 0;

	}

    }

#ifdef DEBUG
  display("Winner Index::", s_bmu_w_index,3);
#endif

  return ;

}
//===================================================================
void SOM::setRate(const float* const &eta)
{

  copy(eta,s_eta,3);

  return;
}
//===================================================================
void SOM::getRate(float* const &eta)const
{

  copy(s_eta,eta,3);

  return ;
}
//===================================================================
void SOM::Network_output(const float* const &inp_vect, float* const &y0)
{
  uint nc = s_neuron_Nx * s_neuron_Ny * s_neuron_Nz;
  float sigma=Spread();
 
  int i=0;
  int j=0;
  int k=0;
  int l=0;
  int m=0;

  int count = 0;

  float dnorm=0.0;                       //Lattice distance Norm
  float sumh=0.0;                        //Summation of spread
  float h=0.0;                           //Membership
  float p=0.0; //to set the threshold value
  
  //=============================================
  //getting the winner index
  int x = s_bmu_w_index[0];
  int y = s_bmu_w_index[1];
  int z = s_bmu_w_index[2];
  
  // //=================================//
  
  sumh = 0.0;  
  zero(y0,s_op_dimn_NO);
  zero(yj, nc, s_op_dimn_NO);  
  zero(a_del_x, nc, s_op_dimn_NO);
  zero (s_y, s_op_dimn_NO);
  //===============================//
  for(i = 0; i < s_neuron_Nx; i++)  //moving x lattice
    {
      for(j = 0; j < s_neuron_Ny; j++)//moving y lattice
	{
	  for(k = 0; k < s_neuron_Nz; k++) // moving along z lattice
	    {

	      //computing the distance in lattice space
	      dnorm = pow((x -i), 2.0) + pow((y-j), 2.0) + pow((z-k), 2.0);
	      
	       
	      h = exp(-dnorm / (2 * sigma * sigma));

	     
	      zero(yj[count], s_op_dimn_NO);  

	      for(l = 0; l < s_op_dimn_NO; l++)
		{

		  zero(a_del_x[count],  s_op_dimn_NO);
		  
		  //A*(X-W)
		  for(m = 0; m < s_ip_dimn_NI; m++)
		    {

		      a_del_x[count][l] +=  s_A[count][l*s_ip_dimn_NI + m] * (inp_vect[m] - s_w[count][m]);
		    }

		  //h*(bj + A*(X-W))
		  yj[count][l] = s_b[count][l] + a_del_x[count][l] ;
		  
		  s_y[l] += yj[count][l] *h;
		}
	     
	      count++;
	      sumh += h;

	    }
	}
    }


  
  
  divide(s_y, sumh, y0, s_op_dimn_NO);

 // s_test_sigma = y0;

  //=====================================================//

  //this part is only for action_unit data
  //sigmoid function to threshold the network output value to 1 when
  //positive data is coming from the network else threshold to 0
  
  /*
    for ( int i=0; i<s_op_dimn_NO; i++)
    {

    if (y0[i]>=0)
	

    y0[i]=1;
    else
    y0[i]=0;

    }

  */
  //===========sigmoid function ====================//
 //find the largest value of p and the respective index. (date 26th june 2013)

  float largest=0.0;
  int index=0;
  // std::cout << " network output = " <<std::endl; 
  for ( int i=0; i<s_op_dimn_NO; i++)
    {

      p= 1/(1+ exp (-y0[i]));
      s_test_sigma[i] = p;
      if(p>largest)
	{
	largest =p;
	index =i;
	}

    }

  for(int i=0; i<s_op_dimn_NO; i++)
    {
      // if (p>=0.5)
      if(i==index)
	y0[i]=1;

	 
      else
  	y0[i]=0;
       
      // std::cout << y0[i] <<"\t"; 
       
      //  y0[i] =6.0*p; // remove this 

    }

  //==============old version ======================//
  /*
   for ( int i=0; i<s_op_dimn_NO; i++)
    {

      p= 1/(1+ exp (-y0[i]));

    
     if (p>=0.5)
     
	y0[i]=1;

	 
      else
  	y0[i]=0;
      }
  */
  //============sigmoid function end ==============//


  /*
  
  display("network output is : \t");
  display(y0, s_op_dimn_NO);
  */
  //========================================================//

 
  return;

}



 
//===================================================================//
void SOM:: saveSOM(const char* const pname, const char* const wname, const char* const tname, const char* const jname)
{

  uint nc = s_neuron_Nx * s_neuron_Ny * s_neuron_Nz;

  uint count = 0;

#ifdef DEBUG
  display("\nSaving network parameters ...\n");
#endif


  std::ofstream fp;
  std::ofstream fw;
  std::ofstream ft;
  std::ofstream fj;

  //===============================

  if(pname == NULL)
    {

      open2write(fp,"parameter.txt");
    }
  else
    {
      open2write(fp,pname);
    }

  //===============================


  if(wname == NULL)
    {

      open2write(fw,"weight.txt");
    }
  else
    {
      open2write(fw,wname);
    }

  //===============================


  if(tname == NULL)
    {
      open2write(ft,"theta.txt");
    }
  else
    {
      open2write(ft,tname);

    }

  //===============================

  if(jname == NULL)
    {
      open2write(fj,"Jacobian.txt");
    }
  else
    {
      open2write(fj,jname);
    }


  //===============================

  fp<<s_neuron_Nx<<"\t"<<s_neuron_Ny<<"\t"<<s_neuron_Nz<<std::endl;
  fp<<std::endl;
  fp<<s_op_dimn_NO<<"\t"<<s_ip_dimn_NI<<std::endl;

  //===============================

  for(count = 0; count< nc;count++)
    {
      print(fj, s_A[count], s_op_dimn_NO, s_ip_dimn_NI);
      fj<<std::endl;
      print(ft, s_b[count], s_op_dimn_NO);
      print(fw, s_w[count], s_ip_dimn_NI);

    }

  //===============================
  ft.close();
  fj.close();
  fw.close();
  fp.close();
  //===============================
#ifdef DEBUG
  display("final value");
  for(int i=0; i < nc;i++){
    display (s_A[i], s_op_dimn_NO * s_ip_dimn_NI);
    
  }
  for(int i=0;i<nc;i++){
    display(s_b[i], s_op_dimn_NO);}
  
  printline();

  #endif
  return;

}
//===================================================================
void SOM:: readSOM(const char* const pname, const char* const wname, const char* const tname, const char* const jname)
{


  uint count = 0;

#ifdef DEBUG
  display("\nSaving network parameters ...\n");
#endif


  std::ifstream fp;
  std::ifstream fw;
  std::ifstream ft;
  std::ifstream fj;

  //===============================

  if(pname == NULL)
    {

      open2read(fp,"parameter.txt");
    }
  else
    {
      open2read(fp,pname);
    }

  //===============================


  if(wname == NULL)
    {

      open2read(fw,"weight.txt");
    }
  else
    {
      open2read(fw,wname);
    }

  //===============================


  if(tname == NULL)
    {
      open2read(ft,"theta.txt");
    }
  else
    {
      open2read(ft,tname);

    }

  //===============================

  if(jname == NULL)
    {
      open2read(fj,"Jacobian.txt");
    }
  else
    {
      open2read(fj,jname);
    }


  //===============================
  //reading the network dimensions


  fp>>s_neuron_Nx>>s_neuron_Ny>>s_neuron_Nz;
  fp>>s_op_dimn_NO>>s_ip_dimn_NI;

#ifdef DEBUG
  display(" Read SOM parameters");
  display("Number of Neurons in x-direction", s_neuron_Nx);
  display("Number of Neurons in y-direction", s_neuron_Ny);
  display("Number of Neurons in z-direction", s_neuron_Nz);
  display("Output dimension", s_op_dimn_NO);
  display("Input dimension", s_ip_dimn_NI);
#endif

  uint nc = s_neuron_Nx * s_neuron_Ny * s_neuron_Nz;

  // createBuffers();

  //===============================
  for(count = 0; count< nc;count++)
    {
      read(fj, s_A[count], s_op_dimn_NO, s_ip_dimn_NI);
      read(ft, s_b[count], s_op_dimn_NO);
      read(fw, s_w[count], s_ip_dimn_NI);

    }

  //===============================
  ft.close();
  fj.close();
  fw.close();
  fp.close();
  //===============================

  return;

}


//================anima===20th july 2012===========================//


//===================================================================
void SOM::Parameters_training( const float* const& y0, const float* const& inp_vect, const float* const& y_d)
{
  //y0 is the output obtained from the network  i.e. y . 
  float sigma=Spread();
  float dvnorm;                                    //Change in position norm
  float eta[3];

  int i=0;
  int j=0;
  int k=0;
  int l=0;
  int m=0;


  
  float dnorm=0.0;                       //Lattice distance Norm
  float sumh=0.0;                        //Summation of spread
  float h=0.0;                           //Membership
  float p=-1.0;
  uint count = 0;
  bool if_sigmoid= false; // when sigmoid funct is used at o/p node
  //===============================
  //getting the winner index (\mu is the winner)
  int x = s_bmu_w_index[0];
  int y = s_bmu_w_index[1];
  int z = s_bmu_w_index[2];
  uint nc = s_neuron_Nx * s_neuron_Ny * s_neuron_Nz;
  getRate(eta);
 
  float eta_w = eta[0];
  zero(s_delta_w, nc, s_ip_dimn_NI);// date 17th june 2013
  zero(s_delta_b, nc, s_op_dimn_NO);
  zero(s_delta_A, nc, s_op_dimn_NO * s_ip_dimn_NI); 
  // C= A-B  
  //sub(y_d,y0,s_y_err, s_op_dimn_NO);  //output err :  s_y_err=Y_d- Y_out

  //edited for test
  sub(y_d,s_test_sigma,s_y_err, s_op_dimn_NO);  //output err :  s_y_err=Y_d- y_sig

  //========================//

  //this part does : y0(1-y0)
  if(if_sigmoid == true)
    {

      for (int i=0; i<s_op_dimn_NO; i++){
	s_temp[i] =1.0;
      }
  
      //1-y0
     // sub (s_temp, y0, s_temp, s_op_dimn_NO);

      //edited for test
//(1-y_sig)
      sub (s_temp, s_test_sigma, s_temp, s_op_dimn_NO);

      //a= a*b
 
      // y0^T (1-y0)


      //edited for test
      // y_sig^T(1-y_sig)
      multiply(s_test_sigma, s_temp, sigma_dvt,  s_op_dimn_NO, 1, s_op_dimn_NO, TRANS, NO_TRANS  );

    }

  
  //===========================// 



  // display("error=");
  // display(s_y_err,s_op_dimn_NO);

  for(i = 0; i < s_neuron_Nx; i++)  //moving x lattice
    {
      
      for(j = 0; j < s_neuron_Ny; j++)//moving y lattice
	{
	  for(k = 0; k < s_neuron_Nz; k++) // moving along z lattice
	    {
	      dnorm = pow(float(x -i), 2.0) + pow(float(y-j), 2.0) + pow(float(z-k), 2.0);
	      h = exp(-dnorm / (2 * sigma * sigma));

	     
	      sumh = sumh+h;
	      //weight updation is done here =============================//
	      //c- a-b
	      sub(inp_vect, s_w[count], s_delta_w[count], s_ip_dimn_NI); //date
									 //17th
									 //june
									 //2013


	      multiply(s_delta_w[count], (eta_w*h), s_ip_dimn_NI); // data
								   // 17th
								   // june 2013


	      



	   
	      //============================================//
	      for (l=0; l<s_op_dimn_NO; l++){


		s_delta_b[count][l] =  s_y_err[l] * h ; 
		

		if(if_sigmoid == true){


		  s_delta_b[count][l] =  s_delta_b[count][l] *  sigma_dvt[0]; 
		}
		
		for(m=0; m<s_ip_dimn_NI; m++){

		  s_delta_A[count] [ l * s_ip_dimn_NI +m ] = s_y_err[l] * h * (inp_vect[m] - s_w[count][m]);


		  if(if_sigmoid == true){
		    
		    s_delta_A[count] [ l * s_ip_dimn_NI +m ] = s_delta_A[count] [ l * s_ip_dimn_NI +m ] * sigma_dvt[0];
		  }

		}
		
	      }


	      count++;
	    
	    }// Neuron index ... for each neuron
	  
	}
    }

  // //normalize by dividing with /sum hi
  // parameter update is done here
  
  for(int c=0; c<count; c++)
    {

      // update the parameters only if there is any loss . or if there
      // is any mismatch of the network output and desired output component.
      for (l=0; l<s_op_dimn_NO; l++){

	
      
	s_delta_b[c][l] = s_delta_b[c][l]/sumh;


	
      	for(m=0; m<s_ip_dimn_NI; m++){
	  s_delta_A[c][l * s_ip_dimn_NI +m] = s_delta_A[c][l * s_ip_dimn_NI +m]/ sumh;
	}


	

	
	if (y0[l] == y_d[l]) {
	  //no update

	  s_b[c][l] = s_b[c][l];


	  for(m=0; m<s_ip_dimn_NI; m++){
	    s_A[c][l * s_ip_dimn_NI +m] =  s_A[c][l * s_ip_dimn_NI +m]; 

	  }

	     
	}

	else {
	  // a= a+ eta * delta b
	  s_b[c][l] +=  s_delta_b[c][l] * eta_w;

	  for(m=0; m<s_ip_dimn_NI; m++){
	    s_A[c][l * s_ip_dimn_NI +m] +=  s_delta_A[c][l * s_ip_dimn_NI +m] * eta_w;

	  }


	}




      }//update of each output node 


      //==========update the weights here ===============//

      divide (s_delta_w[c], sumh , s_ip_dimn_NI);

      add(s_w[c], s_delta_w[c], s_ip_dimn_NI);

      //=========weight is updated====data 17th june 2013===//

    }





  
  return;
}
//===================================================================




void SOM:: readweight(const char* const wname)
{


  uint count = 0;



  std::ifstream fw;
 
  //===============================


  if(wname == NULL)
    {

      open2read(fw,"weight.txt");
    }
  else
    {
      open2read(fw,wname);
    }

  //===============================



  uint nc = s_neuron_Nx * s_neuron_Ny * s_neuron_Nz;

  //===============================
  for(count = 0; count< nc;count++)
    {
     
      read(fw, s_w[count], s_ip_dimn_NI);

    }

  //===============================

  fw.close();
 
  //===============================

  return;

}



//===================================================================


// //EOF
