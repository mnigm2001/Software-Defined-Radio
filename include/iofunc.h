/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_IOFUNC_H
#define DY4_IOFUNC_H

// add headers as needed
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>

// declaration of a function prototypes
void printRealVector(const std::vector<float> &);

void printComplexVector(const std::vector<std::complex<float>> &);

void readBinData(const std::string, std::vector<float> &);

void writeBinData(const std::string, const std::vector<float> &);

void write_audio_data(const std::string out_fname, const std::vector<float> &audio_left, const std::vector<float> &audio_right);

void printFirstLast(const std::vector<float> &, const char*);

void printLastN(const std::vector<float> &vector, const unsigned int N, const char* vector_name);

void printFirstLastN(const std::vector<float> &vector, const unsigned int firstN, const unsigned int lastN, const char* vector_name);

void readStdinBlockData(unsigned int, unsigned int, std::vector<float> &);
#endif // DY4_IOFUNC_H
