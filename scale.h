/* ----------------------------------------------------------
created on 05.09.2009

class to change the domain of  data from one range to another range


 * ---------------------------------------------------------- */
#ifndef _SCALE_
#define _SCALE_

#include "constants.h"


//==============================================================


class scale
{

 private:

  float *m_inMax;
  float *m_inMin;

  float *m_outMax;
  float *m_outMin;

  uint n;

  void getInputLimit(float* const& inMax, float* const& inMin)const;
  void getOutputLimit(float* const& outMax, float* const& outMin)const;
  const uint& size()const;


  void createBuffers();

  void clearBuffers();

 public:

  scale();

  scale(const float* const& inMax, const float* const& inMin, const float* const& outMax, const float* const& outMin, const uint& n);

  scale(const scale& b);

  scale& operator=(const scale& b);

  ~scale();

  void transform(const float* const& in, float* const& out);

};



#endif
//==============================================================
//==============================================================
//EOF
