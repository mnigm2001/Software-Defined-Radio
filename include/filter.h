/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>
#include <cmath>

// declaration of a function prototypes
void allPass(const std::vector<float> &,std::vector<float> &,std::vector<float> &);

void bandPass(float, float, float, unsigned short int, std::vector<float> &);

void fmPLL(const std::vector<float> &,std::vector<float> &,std::vector<float> & ,float , float , float , float , float);

void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);

void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);

void convolveBlockFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &, \
  std::vector<float> &);

void convolveBlockFastFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &,\
  std::vector<float> &, const unsigned int, const bool);

void convolveBlockResampleFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, \
  std::vector<float> &state, const unsigned int audio_decim, const unsigned int audio_upsamp, bool printData);

void upsample (const std::vector<float> &, std::vector<float> &, const int);

void downsample(std::vector<float> &, const std::vector<float> &, const unsigned short int);

void fmDemod(std::vector<float> &, const std::vector<float> &, const std::vector<float> &, float &, float &);

void setVec(const std::vector<float> &, std::vector<float> &, int , int , int mode=1);



#endif // DY4_FILTER_H
