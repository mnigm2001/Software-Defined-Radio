#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <random>
using namespace std;

#define PI 3.14159265358979323846

//Print and set functions
void printRealVector(const std::vector<int> &x);
void printRealVectorFloat(const std::vector<float> &x);
void printLastN(const std::vector<float> &vector, const unsigned int N, const char* vector_name);
void printFirstLast(const std::vector<int> &vector, const unsigned int firstN, const unsigned int lastN, \
	 const char* vector_name);
void setVec(const std::vector<float> &vec1, std::vector<float> &vec2, int begin, int end, int mode=1);

// Filtering functions
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h);
void downsample(std::vector<float> &output, const std::vector<float> &input, int ds_coeff);
void upsample (const std::vector<float> &x, std::vector<float> &xu, const int up_rate);
void convolveBlockFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, \
	std::vector<float> &state, const unsigned int audio_decim);
void convolveBlockUpsampleFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, \
	std::vector<float> &state, const unsigned int audio_upsamp);
void convolveBlockUpWithDownFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, \
	std::vector<float> &state, const unsigned int audio_decim, const unsigned int audio_upsamp);


//This function generates random int values in the range min < number < max
static std::vector<int> generate_data(size_t size, int min, int max)
{
    using value_type = int;
    // We use static in order to instantiate the random engine
    // and the distribution once only.
    // It may provoke some thread-safety issues.
    static std::uniform_int_distribution<value_type> distribution(
        min,
        max);
    static std::default_random_engine generator;

    std::vector<value_type> data(size);
    std::generate(data.begin(), data.end(), []() { return distribution(generator); });
    return data;
}



int main()
{

	// Resampling Testing
	int audio_upsamp = 147;
	int audio_decim = 800;
	int audio_taps = audio_upsamp * 101;

	/*
	//SLIDES TEST
	int audio_upsamp = 3;
	int audio_decim = 4;
	int num_taps = 17;
	int audio_taps = audio_upsamp * num_taps;

	std::vector<int> x_gen = generate_data(36, -100, 100);
	printRealVector(x_gen);
	std::vector<float> x(x_gen.begin(), x_gen.end());
	std::vector<float> h;
	impulseResponseLPF(240000, 16000, audio_taps, h);*/

	/*
	//MODE2 TEST
	std::vector<float> h;
	impulseResponseLPF(240000, 16000, audio_taps, h);
	std::vector<float> y;
	std::vector<float> state;
	state.resize(h.size()-1, 0.0);

	//FIRST BLOCK
	std::vector<int> x_gen = generate_data(5600, -100, 100);
	printFirstLast(x_gen, 3, 100, "x");
	std::vector<float> x(x_gen.begin(), x_gen.end());
	convolveBlockUpWithDownFIR(y, x, h, state, audio_decim, audio_upsamp);
	//printRealVectorFloat(state);

	//SECOND BLOCK
	x_gen = generate_data(5600, -500, 500);
	printFirstLast(x_gen, 3, 100, "x");
	std::vector<float> x2(x_gen.begin(), x_gen.end());
	convolveBlockUpWithDownFIR(y, x2, h, state, audio_decim, audio_upsamp);
	cout << " h size = " << h.size() << ", state size = " << state.size() << endl;*/




}

// FIltering Functions
// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
  // allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);

	float NormFc = Fc/(Fs/2);
	for (int i=0; i<num_taps; i++){
		if (i==(num_taps-1)/2){	h[i] = NormFc;}
		else{ h[i] = NormFc * ((std::sin(PI*NormFc*(i-(num_taps-1)/2)))/(PI*NormFc*(i-(num_taps-1)/2)));}
		h[i] *= std::pow(std::sin(i*PI/num_taps),2);
	}
}

// function to downsample a vector
void downsample(std::vector<float> &output, const std::vector<float> &input, int ds_coeff)
{

  /*
  // allocate memory for the output (filtered) data
  output.clear(); output.resize(ceil(input.size()/static_cast<float> (ds_coeff)), 0.0);
  cout << "input size = " << input.size() << "  output size = " << output.size() << endl;

  int k=0;
  for (unsigned int i=0;i<input.size();i=i+ds_coeff){
    cout << "i= " << i << "  k= " << k << endl;
    cout << "input =" << input[i] << endl;
    output[k] = input[i];
    k++;
  }*/

  output.clear(); output.resize(ceil(input.size()/static_cast<float> (ds_coeff)), 0.0);
  cout << "input size = " << input.size() << "  output size = " << output.size() << endl;
  for(unsigned int i = 0; i < output.size(); i++)
  {
    cout << "i= " << i << endl;
    cout << "input =" << input[i] << endl;
    output[i] = input[i*ds_coeff];
  }
}

// function to upsample a vector
void upsample (const std::vector<float> &x, std::vector<float> &xu, const int up_rate)
{
	xu.clear(); xu.resize(x.size()*up_rate, static_cast<float>(0));
	for(int i=0; i<(x.size()*up_rate); i++){
		if((i%up_rate) == 0){ xu[i] = x[i/up_rate]; }
		else{ xu[i] = 0; }
	}
}

// Convolution of y using x and h with state saving
void convolveBlockFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, const unsigned int audio_decim)
{
	// allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()/audio_decim, static_cast<float>(0));
	for (unsigned int m=0; m<(int)x.size(); m+=audio_decim){
			for (unsigned int n=0; n<(int)h.size(); n++){
				if (m-n >= 0){ y[m/audio_decim]	+= h[n] * x[m-n]; }
				else{	y[m/audio_decim] += h[n] * state[m-n+state.size()]; }
			}
	}
	setVec(x, state, x.size()-state.size(), x.size(), 1);
}

// Convolution with only upsampling, NO downsampling (no state saving either)
void convolveBlockUpsampleFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, \
	const unsigned int audio_upsamp)
{
	// allocate memory for the output (filtered) data
	y.clear(); y.resize((x.size()*audio_upsamp), static_cast<float>(0));

	for (int m=0; m<(int)(x.size()*audio_upsamp); m++){
		int phase = m%audio_upsamp;
		int count = 0;
		bool dotFlag = true;
		if((phase == 146) && (m<2000)){ printf(" &&&\n");}
			for (int n=phase; n<(int)h.size(); n+=audio_upsamp){

				if(((count < 10) || (count > 95)) && (m<2000)) { printf( "m = %d\t n = %d\t m-n = %d  \t phase = %d \t\t", m, n, m-n, phase); }
				else{ if(dotFlag && (m<2000)) {printf("\t...\n");} dotFlag = false;}

				if ((m-n) >= 0){
					//y[m]	+= h[n] * x[(m-n)/audio_upsamp];
					if(((count < 10) || (count > 95)) && (m<2000)) { printf("X\ty[%d] = h[%d] * x[%d]\n",m, n, m-n); }
					else{ if(dotFlag && (m<2000)) {printf("\t...\n");} dotFlag = false;}
				}
				else{
					//y[m] += h[n] * state[m-n+state.size()];
					if(((count < 10) || (count > 95)) && (m<2000)) { printf("STATE\ty[%d] = h[%d] * state[%d]\n", m, n, m-n+state.size()); }
					else { if(dotFlag && (m<2000)) {printf("\t...\n");} dotFlag = false;}
				}
				count++;
			}
			if(m<2000) { printf("____________________________________\n"); }
			//Amplify
	}
	//Prepare next state
}

// Working Resampling convolution function
void convolveBlockUpWithDownFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, \
	const unsigned int audio_decim, const unsigned int audio_upsamp)
{
	// allocate memory for the output (filtered) data
	y.clear(); y.resize((x.size()*audio_upsamp)/audio_decim, static_cast<float>(0));

	for (int m=0; m<(int)(x.size()*audio_upsamp); m+=audio_decim){
		int phase = (m)%audio_upsamp;

		int count = 0;
		bool dotFlag = true;
		int mless = 1600;
		int numbegin = 10;
		int numend = 94;

		//if((phase == 146) && (m<2000)){ printf(" &&&\n");}
			for (int n=phase; n<(int)h.size(); n+=audio_upsamp){

				if(((count < numbegin) || (count > numend)) && (m<mless)) { printf( "m = %d\t n = %d\t  m-n = %d\t phase = %d \t", m, n, m-n, phase); }
				else{ if(dotFlag && (m<mless)) {printf("\t...\n");} dotFlag = false;}

				if ((m-n) >= 0){
					y[m/audio_decim]	+= h[n] * x[(m-n)/audio_upsamp];
					if(((count < numbegin) || (count > numend)) && (m<mless)) { printf("X\ty[%d] = h[%d] * x[%d]\t y[%d] = h[%d] * x[%d]\t%f = %f * %f\n",m, n, (m-n), m/audio_decim, n, (m-n)/audio_upsamp, y[m/audio_decim], h[n], x[(m-n)/audio_upsamp]); }
					else{ if(dotFlag && (m<mless)) {printf("\t...\n");} dotFlag = false;}
				}
				else{
					y[m/audio_decim] += h[n] * state[m-n+state.size()];
					if(((count < numbegin) || (count > numend)) && (m<mless)) { printf("STATE\ty[%d] = h[%d] * state[%d]\t%f = %f * %f\n", m/audio_decim, n, m-n+state.size(), y[m/audio_decim], h[n], state[m-n+state.size()]); }
					else { if(dotFlag && (m<mless)) {printf("\t...\n");} dotFlag = false;}
				}
				count++;
			}
			if(((count < numbegin) || (count > numend)) && (m<mless)) { printf("y = %f\n",y[m/audio_decim]); }
			if((m<mless)) { printf("____________________________________\n"); }

			//Amplify
			y[m/audio_decim] += y[m/audio_decim] * audio_upsamp;
	}

	//Prepare next state
	printf("saved state values = [ ");
	int k=audio_upsamp-1;
	for (int i=audio_upsamp*x.size()-state.size(); i<(int)(audio_upsamp*x.size() - audio_upsamp);i+=audio_upsamp){
		//fprintf(stdout, "k = %d\t i = %d\ti/upsamp+1 = %d\tstate = %f\tx = %f\n", k, i, i/audio_upsamp + 1, state[k], x[i/audio_upsamp + 1]);
		state[k] = x[(i/audio_upsamp )+ 1];
		printf("%f, ",state[k]);
		k+=audio_upsamp;
	}
	printf("]\n");
}


// Print and set functions
void setVec(const std::vector<float> &vec1, std::vector<float> &vec2, int begin, int end, int mode)
{
	if(mode == 1){
		int k=0;
		for (int i=begin; i<end;i++){
			vec2[k] = vec1[i];
			k++;
		}
	}
	else if(mode == 2){
		int k=0;
		for (int i=begin; i<end;i++){
			vec2[i] = vec1[k];
			k++;
		}
	}
}

void printRealVector(const std::vector<int> &x)
{
	std::cout << "Printing int vector of size " << x.size() << "\n";
	for (unsigned int i = 0; i < x.size(); i++)
		std::cout << x[i] << " ";
	std::cout << "\n";
}

void printRealVectorFloat(const std::vector<float> &x)
{
	std::cout << "Printing float vector of size " << x.size() << "\n";
	for (unsigned int i = 0; i < x.size(); i++)
		std::cout << x[i] << " ";
	std::cout << "\n";
}

// This function prints last n values in a vector
void printLastN(const std::vector<float> &vector, const unsigned int N, const char* vector_name)
{

	if(vector.size() > N){
		std::cout << vector_name << " = [ ... , ";
		for (unsigned int i=vector.size()-N; i<vector.size()-1; i++){
			std::cout << vector[i] << ", ";
		}
		std::cout << vector[vector.size()-1] << " ]" << std::endl;
	}
	else { std::cout << "Vector size < N" << std::endl; }
}

void printFirstLast(const std::vector<int> &vector, const unsigned int firstN, const unsigned int lastN, \
	const char* vector_name)
{

		std::cout << vector_name << " = [ ";
		for(unsigned int i=0; i<firstN; i++){
			std::cout << vector[i] << ", ";
		}
		std::cout << " ... , ";
		for (unsigned int i=vector.size()-lastN; i<vector.size()-1; i++){
			std::cout << vector[i] << ", ";
		}
		std::cout << vector[vector.size()-1] << " ]" << std::endl;
}
