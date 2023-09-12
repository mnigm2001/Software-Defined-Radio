#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <unistd.h>
#define PI 3.14159265358979323846


// Queue size for our program
#define QUEUE_ELEMS 6

struct PARAMS{
	int mode;
	int rf_Fs;
	int if_fs;
	float audio_Fs;
	int rf_decim;
	int audio_decim;
	int audio_upsamp;
	int audio_taps;
};

//Print and set functions
void printRealVector(const std::vector<int> &x);
void printRealVectorFloat(const std::vector<float> &x);
void printFirstLastN(const std::vector<float> &vector, const unsigned int firstN, const unsigned int lastN, \
	 const char* vector_name);
void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data);
void setVec(const std::vector<float> &vec1, std::vector<float> &vec2, int begin, int end, int mode);

// Filtering functions
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h);
void convolveBlockFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, \
	std::vector<float> &state, const unsigned int audio_decim);
void convolveBlockFastFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h,\
   std::vector<float> &state, const unsigned int audio_decim, const bool printData);
void convolveBlockUpsampleFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, \
	std::vector<float> &state, const unsigned int audio_upsamp);
void convolveBlockResampleFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, \
  std::vector<float> &state, const unsigned int audio_decim, const unsigned int audio_upsamp, bool printData);
void fmDemod(std::vector<float> &fm_demod, const std::vector<float> &I, const std::vector<float> &Q, float &prev_i, float &prev_q);

// Example Functions
void my_generate(std::vector<int>& elem);
void my_compute(std::vector<int>& elem);
void my_produce(std::queue<std::vector<int>> &my_queue, std::mutex& my_mutex, \
  std::condition_variable& my_cvar);
void my_consume(std::queue<std::vector<int>> &my_queue, std::mutex& my_mutex, \
  std::condition_variable& my_cvar);


// Consumer function
void RF_FrontEnd(std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, \
  std::condition_variable& my_cvar, struct PARAMS &audio_params)
{
	auto cpu_id = sched_getcpu();
	std::cerr << " Starting FE thread " << std::this_thread::get_id() << " processor " << cpu_id << std::endl;

	int mode = audio_params.mode;
	std::cerr << "mode = " << mode << std::endl;
	int rf_Fc = 100000;
	int rf_taps = 151;

  // coefficients for the front-end low-pass filter
	std::vector<float> rf_coeff;
	impulseResponseLPF(audio_params.rf_Fs, rf_Fc, rf_taps, rf_coeff);

  // select a block_size that is a multiple of KB
	// and a multiple of decimation factors
	unsigned int block_size;
	if ( (mode == 0) || (mode == 1) ) { block_size = (1024 * audio_params.rf_decim * audio_params.audio_decim * 2); }
	else if (mode == 2) { block_size = 7 * audio_params.audio_decim * audio_params.rf_decim * 2; }
	else if (mode == 3) { block_size = 7 * audio_params.audio_decim * audio_params.rf_decim * 2; }
  std::cerr << "block_size = " << block_size << std::endl;

  // states for rf front-end processing
  float prev_i=0.0, prev_q=0.0;
  std::vector<float> I_state;
  std::vector<float> Q_state;
  I_state.resize(rf_taps-1, 0.0);
  Q_state.resize(rf_taps-1, 0.0);

  // Parameters storing each I and Q block
  std::vector<float> I_in_block;
  std::vector<float> Q_in_block;
  bool flag = false;

  // if the number of samples in the last block is less than the block size
	// it is fine to ignore the last few samples from the raw IQ file
  for (unsigned int block_id=0; ; block_id++){
		std::vector<float> iq_data(block_size);
		readStdinBlockData(block_size, block_id, iq_data);
		if ((std::cin.rdstate()) != 0){
			std::cerr << "End of input stream reached" << std::endl;
			exit(1);
		}
	 //if( block_id < 300){ std::cerr << "\n___________________Read block " << block_id << "_________________" << std::endl; }
   //printFirstLastN(iq_data, 5, 5, "iq_data");
	 auto old_id = cpu_id;

   I_in_block.clear();Q_in_block.clear();
   flag = false;
   //Seperate the IQ data into 2 vectors
   std::partition_copy(iq_data.begin(),
                        iq_data.end(),
                        std::back_inserter(I_in_block),
                        std::back_inserter(Q_in_block),
                        [&flag](int) { return flag = !flag; });

    // Front-end I data filtering and downsampling
    std::vector<float> I_filt;
    std::vector<float> I_decim;
    convolveBlockFastFIR(I_filt, I_in_block, rf_coeff, I_state, audio_params.rf_decim, false);
    //printFirstLastN(I_in_block, 5, 5, "I_in_block");

    // Front-end Q data filtering and downsampling
    std::vector<float> Q_filt;
    std::vector<float> Q_decim;
    convolveBlockFastFIR(Q_filt, Q_in_block, rf_coeff, Q_state, audio_params.rf_decim, false);

    // FM demod
    std::vector<float> fm_demod;
    fmDemod(fm_demod, I_filt, Q_filt, prev_i, prev_q);
    //printFirstLastN(fm_demod, 5, 5, "fm_demod");

		cpu_id = sched_getcpu();
		if(cpu_id != old_id){
			std::cerr << " changed FE thread " << std::this_thread::get_id() <<  " on processor " << cpu_id << std::endl;
		}

    //pass data to the queue using mutexes and condition variables as a sync mechanism
    std::unique_lock<std::mutex> my_lock(my_mutex);
    //check if the queue is FULL - needed to be done in a custom way
    while(my_queue.size() >= QUEUE_ELEMS){
      my_cvar.wait(my_lock);
    }
		std::cerr << "\t PRODUCED FmDemod\n";
    my_queue.push(fm_demod);
    my_cvar.notify_one();
    my_lock.unlock();

  }
}

void RF_MONO(std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, \
  std::condition_variable& my_cvar, struct PARAMS &audio_params, std::vector<float> &state_mono)
{
	auto cpu_id = sched_getcpu();
	std::cerr << " Starting MONO thread " << std::this_thread::get_id() << " processor " << cpu_id << std::endl;
	int audio_Fc = 16000;
	int mode = audio_params.mode;

	// coefficients for the filter to extract mono audio
	std::vector<float> audio_coeff;
	if ( (mode == 0) || (mode == 1) ) { impulseResponseLPF(audio_params.if_fs, audio_Fc, audio_params.audio_taps, audio_coeff); }
	else if ( (mode == 2) || (mode == 3) ) { impulseResponseLPF(audio_params.if_fs*audio_params.audio_upsamp, \
		audio_Fc, audio_params.audio_taps, audio_coeff); }

  while(1)
  {
		auto old_id = cpu_id;
    // extract data from the queue using mutexes and cond variables as a sync mechanism
    std::unique_lock<std::mutex> my_lock(my_mutex);
    while(my_queue.empty()){
      my_cvar.wait(my_lock);
    }
		std::cerr << " CONSUMED FmDemod\n";
    std::vector<float> fm_demod = my_queue.front();
    my_queue.pop();
    my_cvar.notify_one();
    my_lock.unlock();

    //Extracting Mono Audio
		// extract the mono audio data through filtering
		std::vector<float> audio_filt;
		if((mode == 0) || (mode == 1)) { convolveBlockFastFIR(audio_filt, fm_demod, audio_coeff, state_mono, audio_params.audio_decim, false); }
		else if ((mode == 2) || (mode == 3)) { convolveBlockResampleFIR(audio_filt, fm_demod, audio_coeff, state_mono, \
			audio_params.audio_decim, audio_params.audio_upsamp, false); }
		//printFirstLastN(audio_filt, 5, 5, "audio_filt");

		//Convert audio data to 16-bits and write to standard output
		std::vector<short int> wav_data(audio_filt.size());
		for ( unsigned int k=0; k<audio_filt.size(); k++){
			if (std::isnan(audio_filt[k])) wav_data[k] = 0;
			else wav_data[k] = static_cast<short int> (audio_filt[k] * 16384);
		}
		fwrite(&wav_data[0], sizeof(short int), wav_data.size(), stdout);

		cpu_id = sched_getcpu();
		if(cpu_id != old_id){
			std::cerr << " changed mono thread " << std::this_thread::get_id() <<  " on processor " << cpu_id << std::endl;
		}

    // terminate if queue is empty and no input on
    if (((std::cin.rdstate()) != 0) && (my_queue.empty())){
      break;
    }
  }
}

//MAIN
int main(int argc, char* argv[])
{
  int mode = 0;
	//Set operating mode based on cmd line arguments
	if(argc < 2){
		std::cerr << "Operating in default mode 0:" << std::endl;
	}else if (argc == 2){
		mode = atoi(argv[1]);
		if(mode > 3){
			std::cerr << "Wrong mode " << mode << std::endl;
			exit(1);}
		else { std::cerr << "Operating in mode:" << mode << std::endl; }
	}else {
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or " << std::endl;
		std::cerr << "Usage: " << argv[0] << " <mode>" << std::endl;
		std::cerr << "\t\t <mode> is a value from 0 to 3" << std::endl;
		exit(1);
	}

  int rf_Fs, if_fs, rf_decim, audio_decim, audio_upsamp, audio_taps;	//SHOULD WE MAKE THEM ALL FLOAT?
	float audio_Fs;
	//Set Fs for each processing stage and mode
	if (mode == 1){ rf_Fs = 1440000; if_fs = 288000; audio_Fs = 48000.0; rf_decim = 5; audio_decim = 6; audio_upsamp = 0; audio_taps = 101;}
	else if (mode == 2){ rf_Fs = 2400000; if_fs = 240000; audio_Fs = 44100.0; rf_decim = 10; audio_decim = 800; audio_upsamp = 147; audio_taps = 101*audio_upsamp;}
	else if (mode == 3){ rf_Fs = 960000; if_fs = 320000; audio_Fs = 44100.0; rf_decim = 3; audio_decim = 3200; audio_upsamp = 441; audio_taps = 101*audio_upsamp;}
	else { rf_Fs = 2400000; if_fs = 240000; audio_Fs = 48000.0; rf_decim = 10; audio_decim = 5; audio_upsamp = 0; audio_taps = 101;}

	// Structure to be passed to threads
	struct PARAMS my_params;
	my_params = PARAMS{
		mode,
		rf_Fs,
		if_fs,
		audio_Fs,
		rf_decim,
		audio_decim,
		audio_upsamp,
		audio_taps
	};

  // states for Mono processing
	std::vector<float> state_mono;
	state_mono.resize(audio_taps-1, static_cast<float>(0));

  // context used for thread synch through queue
  std::queue<std::vector<float>> my_queue;
  std::mutex my_mutex;
  std::condition_variable my_cvar;

  // producer thread -- RF Frond-end
  std::thread t_FRE = std::thread(RF_FrontEnd, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar), \
		std::ref(my_params));

  // consumer thread1 -- Mono and Stereo
  std::thread t_consPrint = std::thread(RF_MONO, std::ref(my_queue), std::ref(my_mutex), \
    std::ref(my_cvar), std::ref(my_params), std::ref(state_mono));

  t_FRE.join();
  t_consPrint.join();

  return 0;

  /*
  // EXAMPLE
  // initialize the random generator (illustrative)
  srand(RANDOM_SEED);

  // context used for thread synch through queue
  std::queue<std::vector<int>> my_queue;
  std::mutex my_mutex;
  std::condition_variable my_cvar;

  // producer thread
  std::thread tp = std::thread(my_produce, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar));

  // consumer thread
  std::thread tc = std::thread(my_consume, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar));

  tp.join();
  tc.join();*/

}


///////////////////////////////////////// FILTERING FUNCTIONS //////////////////////////////////////////////////////////////
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

// Convolution with downsampling of y using x and h with state saving
void convolveBlockFastFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, \
	const unsigned int audio_decim, const bool printData)
{

	// allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()/audio_decim, static_cast<float>(0));

	for (int m=0; m<=(int)x.size(); m+=audio_decim){
			for (int n=0; n<(int)h.size(); n++){
				if (m-n >= 0){
					y[m/audio_decim]	+= h[n] * x[m-n]; }
				else{
					y[m/audio_decim] += h[n] * state[m-n+state.size()]; }
			}
	}

	//Prepare next state
	int k=0;
	for (int i=x.size()-state.size(); i<(int)x.size();i++){
		state[k] = x[i];
		k++;
	}
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
		//if((phase == 146) && (m<2000)){ printf(" &&&\n");}
			for (int n=phase; n<(int)h.size(); n+=audio_upsamp){

				//if(((count < 10) || (count > 95)) && (m<2000)) { printf( "m = %d\t n = %d\t m-n = %d  \t phase = %d \t\t", m, n, m-n, phase); }
				//else{ if(dotFlag && (m<2000)) {printf("\t...\n");} dotFlag = false;}

				if ((m-n) >= 0){
					//y[m]	+= h[n] * x[(m-n)/audio_upsamp];
					//if(((count < 10) || (count > 95)) && (m<2000)) { printf("X\ty[%d] = h[%d] * x[%d]\n",m, n, m-n); }
					//else{ if(dotFlag && (m<2000)) {printf("\t...\n");} dotFlag = false;}
				}
				else{
					//y[m] += h[n] * state[m-n+state.size()];
					//if(((count < 10) || (count > 95)) && (m<2000)) { printf("STATE\ty[%d] = h[%d] * state[%d]\n", m, n, m-n+state.size()); }
					//else { if(dotFlag && (m<2000)) {printf("\t...\n");} dotFlag = false;}
				}
				count++;
			}
			//if(m<2000) { printf("____________________________________\n"); }
			//Amplify
	}
	//Prepare next state
}

// Convolution with resampling of y using x and h with state saving
void convolveBlockResampleFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, \
	const unsigned int audio_decim, const unsigned int audio_upsamp, bool printData)
{
	// allocate memory for the output (filtered) data
	y.clear(); y.resize((x.size()*audio_upsamp)/audio_decim, static_cast<float>(0));

	for (int m=0; m<(int)(x.size()*audio_upsamp); m+=audio_decim){
		int phase = (m)%audio_upsamp;
			for (int n=phase; n<(int)h.size(); n+=audio_upsamp){
				if ((m-n) >= 0){
					y[m/audio_decim]	+= h[n] * x[(m-n)/audio_upsamp]; }
				else{
					y[m/audio_decim] += h[n] * state[m-n+state.size()]; }
			}
			//Amplify output signal
			y[m/audio_decim] += y[m/audio_decim] * audio_upsamp;
	}

	//Prepare next state
	int k=audio_upsamp-1;
	for (int i=audio_upsamp*x.size()-state.size(); i<(int)(audio_upsamp*x.size() - audio_upsamp);i+=audio_upsamp){
		state[k] = x[(i/audio_upsamp )+ 1];
		k+=audio_upsamp;
	}
}

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

////////////////////////////////////////////// PRINTING FUNCTIONs /////////////////////////////////////////////////
void printRealVector(const std::vector<int> &x)
{
	std::cerr << "Printing int vector of size " << x.size() << "\n";
	for (unsigned int i = 0; i < x.size(); i++)
		std::cerr << x[i] << " ";
	std::cerr << "\n";
}

void printRealVectorFloat(const std::vector<float> &x)
{
	std::cerr << "Printing float vector of size " << x.size() << "\n";
	for (unsigned int i = 0; i < x.size(); i++)
		std::cerr << x[i] << " ";
	std::cerr << "\n";
}

void printFirstLastN(const std::vector<float> &vector, const unsigned int firstN, const unsigned int lastN, \
	const char* vector_name)
{

		std::cerr << vector_name << " = [ ";
		for(unsigned int i=0; i<firstN; i++){
			std::cerr << vector[i] << ", ";
		}
		std::cerr << " ... , ";
		for (unsigned int i=vector.size()-lastN; i<vector.size()-1; i++){
			std::cerr << vector[i] << ", ";
		}
		std::cerr << vector[vector.size()-1] << " ]" << std::endl;
}

//data in the redirected input is interpreted as bytes
void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data){
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for(int k=0; k<(int)num_samples; k++){
		//automaticallt normalizes the data to the range - 1 to +1
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
	}
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

//////////////////////////////// TEMPORARY STORAGE //////////////////////////////////////
/*
void RF_FE(std::queue<std::vector<int>> &my_queue, std::mutex& my_mutex, \
  std::condition_variable& my_cvar, int &mode, std::vector<float> &iq_data)
{
  // thread runs until a custom end cond occurs
  while(generated_elems < TOTAL_ELEMS)
  {
    //generation and custom processing of elements before
    //passing them to another thread through a queue
    //any processing idependent of other threads should not use locks
    // Parameters storing each I and Q block
    std::vector<float> I_in_block;
    std::vector<float> Q_in_block;
    bool flag = false;

   I_in_block.clear();Q_in_block.clear();
   flag = false;
   //Seperate the IQ data into 2 vectors
   std::partition_copy(iq_data.begin(),
                        iq_data.end(),
                        std::back_inserter(I_in_block),
                        std::back_inserter(Q_in_block),
                        [&flag](int) { return flag = !flag; });

    // Front-end I data filtering and downsampling
    std::vector<float> I_filt;
    std::vector<float> I_decim;
    convolveBlockFastFIR(I_filt, I_in_block, rf_coeff, I_state, rf_decim, false);
    printFirstLastN(I_in_block, 5, 5, "I_in_block");

    // Front-end Q data filtering and downsampling
    std::vector<float> Q_filt;
    std::vector<float> Q_decim;
    convolveBlockFastFIR(Q_filt, Q_in_block, rf_coeff, Q_state, rf_decim, false);

    // FM demod
    std::vector<float> fm_demod;
    fmDemod(fm_demod, I_filt, Q_filt, prev_i, prev_q);
    printFirstLastN(fm_demod, 5, 5, "fm_demod");
    //bookkeeping for the end cond (very custom)
    generated_elems++;

    //pass data to the queue using mutexes and
    // condition variables as a sync mechanism
    std::unique_lock<std::mutex> my_lock(my_mutex);
    //check if the queue is FULL - needed to be done in a custom way
    while(my_queue.size() >= QUEUE_ELEMS){
      my_cvar.wait(my_lock);
    }

    my_queue.push(fm_demod);
    my_cvar.notify_one();
    my_lock.unlock();
  }
}*/
