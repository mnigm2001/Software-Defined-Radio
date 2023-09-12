/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "iofunc.h"

// All pass filter function
void allPass(const std::vector<float> &input_block,std::vector<float> &state_block,std::vector<float> &output_block)
{
	std::vector<float> tmp;
	output_block.clear();
	tmp.clear();

	int size= state_block.size()+1;
	tmp=state_block;

	for (int unsigned k=0; k<state_block.size(); k++){
		state_block[state_block.size()-1-k]=input_block[input_block.size()-1-k];
	}

	tmp.insert(tmp.end(), input_block.begin(), (input_block.end()-tmp.size()));
	output_block=tmp;
}

// PLL function
void fmPLL(const std::vector<float> &PLLIn, std::vector<float> &ncoOut,std::vector<float> &state ,float freq, \
	float Fs, float ncoScale, float phaseAdjust, float normBandwidth)
{
	float Cp=2.666;
	float Ci=3.555;

	float Kp= normBandwidth*Cp;
	float Ki= (normBandwidth*normBandwidth)*Ci;

	ncoOut.clear(); ncoOut.resize(PLLIn.size()+1, 0.0);

	float integrator = state[0];
	float phaseEst = state[1];
	float feedbackI = state[2];
	float feedbackQ = state[3];
	ncoOut[0]=state[4];
	float trigOffset = state[5];

	float errorI,errorQ,errorD,trigArg;

	for (int unsigned k=0; k<PLLIn.size(); k++){
		// phase detector
		errorI = PLLIn[k] * (feedbackI);
		errorQ = PLLIn[k] * (-1*feedbackQ);

		//four-quadrant arctangent discriminator for phase error detection
		errorD = std::atan2(errorQ, errorI);

		//loop filter
		integrator = integrator + Ki*errorD;

		//update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator;

		//internal oscillator
		trigOffset += 1;
		trigArg = 2*PI*(freq/Fs)*(trigOffset) + phaseEst;
		feedbackI = std::cos(trigArg);
		feedbackQ = std::sin(trigArg);
		ncoOut[k+1] = std::cos(trigArg*ncoScale + phaseAdjust);
	}

	state[0]=integrator;
	state[1]=phaseEst;
	state[2]=feedbackI;
	state[3]=feedbackQ;
	state[4]=ncoOut.back();
	state[5]=trigOffset;
}

// THis function computes the bandpass filter coefficients
void bandPass(float Fs, float Fb, float Fe, unsigned short int N_taps, std::vector<float> &coeff)
{
	coeff.clear(); coeff.resize(N_taps, 0.0);
	float Normcenter = ((Fe+Fb)/2)/(Fs/2);
	float Normpass = (Fe-Fb)/(Fs/2);

	for (int i=0; i<N_taps; i++){
		if (i==(N_taps-1)/2){
			coeff[i] = Normpass;
		}
		else{
			coeff[i] = Normpass*((std::sin(PI*Normpass/2*(i-(N_taps-1)/2)))/(PI*Normpass/2*(i-(N_taps-1)/2)));
		}
		coeff[i]= coeff[i]*std::cos(i*(PI)*Normcenter);
		coeff[i] = coeff[i]*(std::sin(i*PI/N_taps))*(std::sin(i*PI/N_taps));
	}
}


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

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
  // allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);

	for (unsigned int m=0; m<y.size(); m++){
		for (unsigned int n=0; n<h.size(); n++){
			if ((m-n >= 0)&&(m-n < x.size())){
				y[m]	+= h[n]*x[m-n];
			}
		}
	}
}

// Convolution of y using x and h with state saving
void convolveBlockFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, \
	std::vector<float> &state)
{
	// allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size(), static_cast<float>(0));

	for (int n=0; n<(int)y.size(); n++){
		for (int k=0; k<(int)h.size(); k++){
			if (n-k >= 0){
				y[n]	+= h[k] * x[n-k]; }									//Using current block audio data i.e. xb for convolution
			else{
				y[n] += h[k] * state[n-k+state.size()]; }					//Using pervious block audio data for convolution
		}
	}

	//Prepare next state
	int k=0;
	for (int i=x.size()-state.size(); i<(int)x.size();i++){
		state[k] = x[i];
		k++;
	}
}


// Convolution with downsampling of y using x and h with state saving
void convolveBlockFastFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, \
	const unsigned int audio_decim, const bool printData)
{

	// allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()/audio_decim, static_cast<float>(0));
	int multCount = 0;

	for (int m=0; m<=(int)x.size(); m+=audio_decim){
			for (int n=0; n<(int)h.size(); n++){
				if (m-n >= 0){
					y[m/audio_decim]	+= h[n] * x[m-n];
					//if(m<1) { multCount++; }
					multCount++;
				}
				else{
					y[m/audio_decim] += h[n] * state[m-n+state.size()];
					//if(m<1) { multCount++; }
					multCount++;
				}
			}
	}
	//std::cout << " multCount = " << multCount << std::endl;

	//Prepare next state
	int k=0;
	for (int i=x.size()-state.size(); i<(int)x.size();i++){
		state[k] = x[i];
		k++;
	}
}

// Convolution with resampling of y using x and h with state saving
void convolveBlockResampleFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, \
	const unsigned int audio_decim, const unsigned int audio_upsamp, bool printData)
{
	// allocate memory for the output (filtered) data
	y.clear(); y.resize((x.size()*audio_upsamp)/audio_decim, static_cast<float>(0));
	int multCount = 0;

	for (int m=0; m<(int)(x.size()*audio_upsamp); m+=audio_decim){
		int phase = (m)%audio_upsamp;
			for (int n=phase; n<(int)h.size(); n+=audio_upsamp){
				if ((m-n) >= 0){
					y[m/audio_decim]	+= h[n] * x[(m-n)/audio_upsamp];
					//if(m<1) { multCount++; }
					multCount++;
				}
				else{
					y[m/audio_decim] += h[n] * state[m-n+state.size()];
					//if(m<1) { multCount++; }
					multCount++;
				}
			}
			//Amplify output signal
			y[m/audio_decim] += y[m/audio_decim] * audio_upsamp;
	}
	//std::cout << " multCount = " << multCount << std::endl;

	//Prepare next state
	int k=audio_upsamp-1;
	for (int i=audio_upsamp*x.size()-state.size(); i<(int)(audio_upsamp*x.size() - audio_upsamp);i+=audio_upsamp){
		state[k] = x[(i/audio_upsamp )+ 1];
		k+=audio_upsamp;
	}
}


//UP sampling
void upsample (const std::vector<float> &x, std::vector<float> &xu, const int up_rate)
{
	xu.clear(); xu.resize(x.size()*up_rate, static_cast<float>(0));
	for(int i=0; i<(int)(x.size()*up_rate); i++){
		if((i%up_rate) == 0){ xu[i] = x[i/up_rate]; }
		else{ xu[i] = 0; }
	}
}

// Function to downsample a vector by a factor ds_coeff
void downsample(std::vector<float> &output, const std::vector<float> &input, const unsigned short int ds_coeff)
{
  // allocate memory for the output
	output.clear(); output.resize(ceil(input.size()/static_cast<float> (ds_coeff)), 0.0);
  for(unsigned int i = 0; i < output.size(); i++)
  {
    output[i] = input[i*ds_coeff];
  }
}

// fm demodulation function
void fmDemod(std::vector<float> &fm_demod, const std::vector<float> &I, const std::vector<float> &Q, float &prev_i, float &prev_q)
{
	// allocate memory for the output
	fm_demod.clear(); fm_demod.resize(I.size(), static_cast<float>(0));

	for (unsigned int k=0; k<I.size(); k++){
		if (I[k]*I[k]+Q[k]*Q[k]==0){
			fm_demod[k]=0;
		}else if(k!=0){
			fm_demod[k]=(I[k]*(Q[k]-Q[k-1])-Q[k]*(I[k]-I[k-1]))/(I[k]*I[k]+Q[k]*Q[k]);
		}else if(k==0){
			fm_demod[k]=(I[k]*(Q[k]-prev_q)-Q[k]*(I[k]-prev_i))/(I[k]*I[k]+Q[k]*Q[k]);
		}
	}

	//Prepare next state
	prev_i=I.back();
	prev_q=Q.back();
}

//Function to append sections of data from vector to another
//The second vector start index is set to zero for mode 1, and the start index of first is zero for mode 2
void setVec(const std::vector<float> &vec1, std::vector<float> &vec2, int begin, int end, unsigned short int mode)
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
