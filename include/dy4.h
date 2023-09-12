/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_DY4_H
#define DY4_DY4_H

//Libraries for project.cpp
#include <cstring>
#include <algorithm>
#include <chrono>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>

// some general and reusable stuff
// our beloved PI constant
#define PI 3.14159265358979323846

// although we use DFT (no FFT ... yet), the number of points for a
// Fourier transform is defined as NFFT (same as matplotlib)
#define NFFT 512

// Queue size for our program
#define QUEUE_ELEMS 6

#endif // DY4_DY4_H
