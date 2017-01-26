/*
//================================================================

*/
#ifndef __SOM_
#define __SOM_
#include "constants.h"


class SOM
{


  private:
  uint s_neuron_Nx;                                   //Number of x Neurons
  uint s_neuron_Ny;                                   //No. of y Neurons
  uint s_neuron_Nz;                                   //No. of z Neurons

  uint s_op_dimn_NO;                                   //No. of Links 
  uint s_ip_dimn_NI;                                   //No. of Inputs

  float s_spread;


  float** s_w; //weight
  float** s_b; //bj
  float** s_A; //Aj

 

  float** s_delta_b;
  float** s_delta_w;
  float** s_delta_A;
  float** a_del_x;

  float** yj;   // output at each neuron

  float* s_y;
  int* s_bmu_w_index;                            //Pointer to Winner Index

   
  float* s_del_xw;
  float* s_eta;                            //Learning Rate

  float* s_y_err;
  float* s_temp;
  float*  sigma_dvt;
  float* s_test_sigma;

  void getIndex(const uint& count, uint& x, uint& y, uint& z);

  void createBuffers();


 public:
  SOM();
  SOM(const uint& Nx, const uint& Ny, const uint& Nz, const uint& NO,const uint& NI);

 

  ~SOM();


  float& Spread(){return(s_spread);}


  void getRate(float* const &eta)const;
  void setRate(const float* const& eta);

  void getWeight(float** const& w)const;
  void get_winner_index(int *winner_index);
  
  void selectWinner(const float* const& vIn);
 
 

  void Network_output(const float* const &inp_vect, float* const &y0);
  


  void Parameters_training( const float* const& y0, const float* const& inp_vect, const float* const& y_d);

 


  void saveSOM(const char* const pname=NULL, const char* const wname=NULL, const char* const tname=NULL, const char* const jname=NULL);


  void readSOM(const char* const pname=NULL, const char* const wname=NULL, const char* const tname=NULL, const char* const jname=NULL);

  void  readweight(const char* const wname=NULL);



};

#endif

//===================================================================


//EOF 
