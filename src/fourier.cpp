/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"
#include "iofunc.h"

// just DFT function (no FFT yet)
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
	Xf.resize(x.size(), static_cast<std::complex<float>>(0));
	for (unsigned int m = 0; m < Xf.size(); m++) {
		for (unsigned int k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<float>> &Xf, std::vector<float> &Xmag)
{
	// only the positive frequencies
	Xmag.resize(Xf.size(), static_cast<float>(0));
 	for (unsigned int i = 0; i < Xf.size(); i++) {
		Xmag[i] = std::abs(Xf[i])/Xf.size();
	}
}

// add your own code to estimate the PSD

void LinearSpacedArray(std::vector<float> &range, float max, float min=0.0, float step=1.0) {
    float N = (max - min)/step;
    for(int i=0; i<N; i++) {
        range.push_back(min + i*step);
    }
}

void estimatePSD(std::vector<float> &freq, std::vector<float> &psd_est, const std::vector<float> &samples, const float Fs){

	// rename the NFFT argument (notation consistent with matplotlib.psd)
 	// to freq_bins (i.e., frequency bins for which we compute the spectrum)
	int freq_bins = NFFT;
	// frequency increment (or resolution of the frequency bins)
	float df = Fs/freq_bins;

	// create the frequency vector to be used on the X axis
 	// for plotting the PSD on the Y axis (only positive freq)
	 LinearSpacedArray(freq, Fs/2, 0.0 , df);

	// design the Hann window used to smoothen the discrete data in order
 	// to reduce the spectral leakage after the Fourier transform
	std::vector<float> hann(freq_bins);
	for(unsigned int i=0;i<hann.size();i++){
		hann[i] = std::pow(std::sin(i*PI/freq_bins),2.0);
	}

	// create an empty list where the PSD for each segment is computed
	std::vector<float> psd_list;

	// samples should be a multiple of frequency bins, so
 	// the number of segments used for estimation is an integer
	// note: for this to work you must provide an argument for the
  // number of frequency bins not greater than the number of samples!
	int no_segments = int(floor(samples.size()/float(freq_bins)));

	// iterate through all the segments
	for (int k=0;k<no_segments;k++){
		// apply the hann window (using pointwise multiplication)
		// before computing the Fourier transform on a segment
		std::vector<float> windowed_samples;
		int temp = 0;
		for (int i=k*freq_bins; i<(k+1)*freq_bins;i++){
			windowed_samples.push_back(samples[i] * hann[temp]);
			temp++;
		}

		// compute the Fourier transform using the built-in FFT from numpy
		std::vector<std::complex<float>> Xf;
		//std::vector<std::complex<float>> twiddles;
		//compute_twiddles(twiddles);
		//FFT_optimized(windowed_samples, Xf, twiddles);

		DFT(windowed_samples, Xf);

		// since input is real, we keep only the positive half of the spectrum
		// however, we will also add the signal energy of negative frequencies
		// to have a better a more accurate PSD estimate when plotting
		Xf.resize(int(freq_bins/2));
		std::vector<float> psd_seg;
		for (unsigned int i=0; i<Xf.size();i++){
			psd_seg.push_back((1/(Fs*freq_bins/2)) * std::pow(std::abs(Xf[i]), 2.0));
			psd_seg[i] = 2*psd_seg[i];
		}

		// translate to the decibel (dB) scale
		for(unsigned int i=0; i<psd_seg.size(); i++){
			psd_seg [i] = 10*std::log10(psd_seg[i]);
		}

		// append to the list where PSD for each segment is stored
		// in sequential order (first segment, followed by the second one, ...)
		for (unsigned int i=0; i<psd_seg.size(); i++){
			psd_list.push_back(psd_seg[i]);
		}
	}

	// compute the estimate to be returned by the function through averaging
	psd_est.resize(int(freq_bins/2));

	// iterate through all the frequency bins (positive freq only)
	// from all segments and average them (one bin at a time ...)
	for(int k=0; k<int(freq_bins/2); k++){
		for (int l=0; l<no_segments; l++){									// iterate through all the segments
			psd_est[k] += psd_list[k + l*int(freq_bins/2)];
		}
		// compute the estimate for each bin
		psd_est[k] = psd_est[k]/no_segments;
	}
	//std::cerr << "psd_size = " << psd_est.size() << std::endl;
	//printRealVector(psd_est);

}

// added IDFT

void IDFT(const std::vector<std::complex<float>> &Xf, std::vector<std::complex<float>> &x) {
	x.resize(Xf.size(), static_cast<std::complex<float>>(0));
	for (unsigned int k = 0; k < x.size(); k++) {
		for (unsigned int m = 0; m < x.size(); m++) {
			std::complex<float> expval(0, 2*PI*(k*m) / Xf.size());
			x[k] += Xf[m] * std::exp(expval);
		}
		x[k] /= Xf.size();
	}
}

// added FFT

unsigned int swap_bits(unsigned int x, unsigned char i, unsigned char j) {

  unsigned char bit_i = (x >> i) & 0x1L;
  unsigned char bit_j = (x >> j) & 0x1L;

  unsigned int val = x;
  val = bit_i ? (val | (0x1L << j)) : (val & ~(0x1L << j));
  val = bit_j ? (val | (0x1L << i)) : (val & ~(0x1L << i));

  return val;
}

unsigned int bit_reversal(unsigned int x, unsigned char bit_size) {

  unsigned int val = x;

  for (int i=0; i < int(bit_size/2); i++)
    val = swap_bits(val, i, bit_size-1-i);

  return val;
}

void compute_twiddles(std::vector<std::complex<float>> &twiddles) {
	std::cout << "(int)twiddles.size() = " << (int)twiddles.size() << std::endl;
  for (int k=0; k<(int)twiddles.size(); k++) {
      std::complex<float> expval(0.0, -2*PI*float(k)/ NFFT);
      twiddles[k] = std::exp(expval);
  }
}

void FFT_recursive(const std::vector<std::complex<float>> &x, \
  std::vector<std::complex<float>> &Xf) {

  if (x.size() > 1) {
    // declare vectors and allocate space for the even and odd halves
    std::vector<std::complex<float>> xe(int(x.size()/2)), xo(int(x.size()/2));
    std::vector<std::complex<float>> Xfe(int(x.size()/2)), Xfo(int(x.size()/2));

    // split into even and odd halves
    for (int k=0; k<(int)x.size(); k++)
      if ((k%2) == 0) xe[k/2] = x[k];
      else xo[k/2] = x[k];

    // call recursively FFT of half size for even and odd halves respectively
    FFT_recursive(xe, Xfe);
    FFT_recursive(xo, Xfo);

    // merge the results from the odd/even FFTs (each of half the size)
    for (int k=0; k<(int)xe.size(); k++) {
        std::complex<float> expval(0.0, -2*PI*float(k)/ x.size());
        std::complex<float> twiddle = std::exp(expval);
        Xf[k]           = Xfe[k] + twiddle * Xfo[k];
        Xf[k+xe.size()] = Xfe[k] - twiddle * Xfo[k];
    }
  } else {
    // end of recursion - copy time domain samples to frequency bins (default values)
    Xf[0] = x[0];
  }
}

void FFT_improved(const std::vector<std::complex<float>> &x, \
  std::vector<std::complex<float>> &Xf, \
  const std::vector<std::complex<float>> &twiddles, \
  const unsigned char recursion_level) {

  if (x.size() > 1) {
    int half_size = int(x.size()/2);
    std::vector<std::complex<float>> xe(half_size), xo(half_size);
    std::vector<std::complex<float>> Xfe(half_size), Xfo(half_size);

    for (int k=0; k<half_size; k++) {
      xe[k] = x[k*2];
      xo[k] = x[k*2+1];
    }

    FFT_improved(xe, Xfe, twiddles, recursion_level+1);
    FFT_improved(xo, Xfo, twiddles, recursion_level+1);

    for (int k=0; k<half_size; k++) {
        Xf[k]           = Xfe[k] + twiddles[k*(1<<(recursion_level-1))] * Xfo[k];
        Xf[k+half_size] = Xfe[k] - twiddles[k*(1<<(recursion_level-1))] * Xfo[k];
    }
  } else {
    Xf[0] = x[0];
  }
}

void FFT_optimized(const std::vector<std::complex<float>> &x, \
  std::vector<std::complex<float>> &Xf, \
  const std::vector<std::complex<float>> &twiddles) {

	std::cout << "9" << std::endl;
	//printComplexVector(x);
	//printComplexVector(Xf);
	//printComplexVector(twiddles);

  unsigned char no_levels = (unsigned char)std::log2((float)x.size());
  for (unsigned int i=0; i<x.size(); i++) {
    Xf[i] = x[bit_reversal(i, no_levels)];
  }
	std::cout << "10" << std::endl;
  unsigned int step_size = 1;

  std::complex<float> tmp;
  for (unsigned char l=0; l<no_levels; l++) {
    for (unsigned int p=0; p<x.size(); p+=2*step_size) {
      for (unsigned int k=p; k<p+step_size; k++) {
        tmp             = Xf[k] + twiddles[(k-p)*(1<<(no_levels-1-l))] * Xf[k+step_size];
        Xf[k+step_size] = Xf[k] - twiddles[(k-p)*(1<<(no_levels-1-l))] * Xf[k+step_size];
        Xf[k]           = tmp;
      }
    }
    step_size *= 2;
  }
	std::cout << "11" << std::endl;
}
