/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"

struct PARAMS{
	int mode;
	int rf_Fs;
	int if_fs;
	float audio_Fs;
	int rf_decim;
	int audio_decim;
	int audio_upsamp;
	int audio_taps;
	int stereo_taps;
};

struct STATES {
  std::vector<float> state_mono;
	std::vector<float> state_stereo;
	std::vector<float> state_carrier;
	std::vector<float> state_stereofilt;
	std::vector<float> state_allpass;
	std::vector<float> state_PLL;
};


// Consumer function
void RF_FrontEnd(std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, \
  std::condition_variable& my_cvar, struct PARAMS &audio_params)
{
	int mode = audio_params.mode;
	std::cout << "mode = " << mode << std::endl;
	int rf_Fc = 100000;
	int rf_taps = 13;	//XX

  // coefficients for the front-end low-pass filter
	std::vector<float> rf_coeff;
	impulseResponseLPF(audio_params.rf_Fs, rf_Fc, rf_taps, rf_coeff);

  // select a block_size that is a multiple of KB
	// and a multiple of decimation factors
	unsigned int block_size;
	if ( (mode == 0) || (mode == 1) ) { block_size = (1024 * audio_params.rf_decim * audio_params.audio_decim * 2); }
	else if (mode == 2) { block_size = 7 * audio_params.audio_decim * audio_params.rf_decim * 2; }
	else if (mode == 3) { block_size = 7 * audio_params.audio_decim * audio_params.rf_decim * 2; }
  std::cout << "block_size = " << block_size << std::endl;

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

	std::chrono::duration<double, std::milli> I_runtime, Q_runtime, fmDemod_runtime, program_runtime;
	std::cout << "I Filter BEFORE " << I_runtime.count() << " milliseconds" << "\n";
	std::cout << "Q Filter BEFORE  " << Q_runtime.count() << " milliseconds" << "\n";
	std::cout << "fmDemod BEFORE " << fmDemod_runtime.count() << " milliseconds" << "\n";

	auto program_start_time = std::chrono::high_resolution_clock::now();
  // if the number of samples in the last block is less than the block size
	// it is fine to ignore the last few samples from the raw IQ file
  for (unsigned int block_id=0; ; block_id++){
		std::vector<float> iq_data(block_size);
		readStdinBlockData(block_size, block_id, iq_data);
		if ((std::cin.rdstate()) != 0){
			std::cout << "End of input stream reached" << std::endl;
			std::cout << "I Filter FINAL " << I_runtime.count() << " milliseconds" << "\n";
			std::cout << "Q Filter FINAL  " << Q_runtime.count() << " milliseconds" << "\n";
			std::cout << "fmDemod FINAL " << fmDemod_runtime.count() << " milliseconds" << "\n";

			auto program_stop_time = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double, std::milli> program_runtime = program_stop_time-program_start_time;
			std::cout << "Program ran for " << program_runtime.count() << " milliseconds" << "\n";
			exit(1);
		}

	 if( block_id < 300){ std::cout << "\n___________________Read block " << block_id << "_________________" << std::endl; }
	 //if(block_id > 0){ break; }

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
		auto start_time_I = std::chrono::high_resolution_clock::now();
    convolveBlockFastFIR(I_filt, I_in_block, rf_coeff, I_state, audio_params.rf_decim, false);
		auto stop_time_I = std::chrono::high_resolution_clock::now();
		I_runtime += stop_time_I-start_time_I;

    //printFirstLastN(I_in_block, 5, 5, "I_in_block");

    // Front-end Q data filtering and downsampling
    std::vector<float> Q_filt;
    std::vector<float> Q_decim;
		auto start_time_Q = std::chrono::high_resolution_clock::now();
    convolveBlockFastFIR(Q_filt, Q_in_block, rf_coeff, Q_state, audio_params.rf_decim, false);
		auto stop_time_Q = std::chrono::high_resolution_clock::now();
		Q_runtime += stop_time_Q-start_time_Q;

    // FM demod
    std::vector<float> fm_demod;
		auto start_time_fmDemod = std::chrono::high_resolution_clock::now();
    fmDemod(fm_demod, I_filt, Q_filt, prev_i, prev_q);
		auto stop_time_fmDemod = std::chrono::high_resolution_clock::now();
		fmDemod_runtime += stop_time_fmDemod-start_time_fmDemod;

		//std::cout << "I Filter ran for " << I_runtime.count() << " milliseconds" << "\n";
		//std::cout << "Q Filter ran for " << Q_runtime.count() << " milliseconds" << "\n";
		//std::cout << "fmDemod ran for " << fmDemod_runtime.count() << " milliseconds" << "\n";

		//std::cout << "fm_demod size = " << fm_demod.size() << std::endl;
		//std::cout << "\t PRODUCED FmDemod\n";

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
}

void RF_STEREO(std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, \
  std::condition_variable& my_cvar, struct PARAMS &audio_params, struct STATES &my_states)
{
	int audio_Fc = 16000;
	int mode = audio_params.mode;
	int audio_decim = audio_params.audio_decim;
	int audio_upsamp = audio_params.audio_upsamp;
	int if_fs = audio_params.if_fs;

	// coefficients for the filter to extract mono audio
	std::vector<float> audio_coeff;
	if ( (mode == 0) || (mode == 1) ) { impulseResponseLPF(if_fs, audio_Fc, audio_params.audio_taps, audio_coeff); }
	else if ( (mode == 2) || (mode == 3) ) { impulseResponseLPF(if_fs*audio_upsamp, audio_Fc, audio_params.audio_taps, \
		 audio_coeff); }

	//stereo coefficients
	std::vector<float> carrier_coeff;
	std::vector<float> stereo_coeff;
	bandPass(if_fs, 18.5e3, 19.5e3, audio_params.stereo_taps, carrier_coeff);
	bandPass(if_fs, 22e3, 54e3, audio_params.stereo_taps, stereo_coeff);

	std::chrono::duration<double, std::milli> allpass_runtime, channel_recov_runtime, carrier_recov_runtime, \
		PLL_runtime, mixer_runtime, mono_stereo_runtime, stereo_runtime;

  while(1)
  {
    // extract data from the queue using mutexes and cond variables as a sync mechanism
    std::unique_lock<std::mutex> my_lock(my_mutex);
    while(my_queue.empty()){
      my_cvar.wait(my_lock);
    }

    std::vector<float> fm_demod = my_queue.front();
    my_queue.pop();
    my_cvar.notify_one();
    my_lock.unlock();

		//apply all pass filter
		std::vector<float> audio_allpass;
		auto start_time_allpass = std::chrono::high_resolution_clock::now();
		allPass(fm_demod, my_states.state_allpass, audio_allpass);
		auto stop_time_allpass = std::chrono::high_resolution_clock::now();
		allpass_runtime += stop_time_allpass-start_time_allpass;

		//filter stereo_taps
		std::vector<float> stereo_filt;
		std::vector<float> carrier_filt;
		auto start_time_channel = std::chrono::high_resolution_clock::now();
		convolveBlockFIR(stereo_filt, fm_demod, stereo_coeff, my_states.state_stereo);
		auto stop_time_channel = std::chrono::high_resolution_clock::now();
		channel_recov_runtime += stop_time_channel-start_time_channel;

		auto start_time_carrier = std::chrono::high_resolution_clock::now();
		convolveBlockFIR(carrier_filt, fm_demod, carrier_coeff, my_states.state_carrier);
		auto stop_time_carrier = std::chrono::high_resolution_clock::now();
		carrier_recov_runtime += stop_time_carrier-start_time_carrier;
		//std::cout << " stereo_filt = " << stereo_filt.size() << std::endl;
		//std::cout << " carrier_filt = " << carrier_filt.size() << std::endl;

    //Extracting Mono Audio
		// extract the mono audio data through filtering
		//Mono Filtering
		std::vector<float> audio_filt;
		if((mode == 0) || (mode == 1)) {
			auto start_time_mono = std::chrono::high_resolution_clock::now();
			convolveBlockFastFIR(audio_filt, audio_allpass, audio_coeff, my_states.state_mono, \
			audio_decim, false);
			auto stop_time_mono = std::chrono::high_resolution_clock::now();
			mono_stereo_runtime += stop_time_mono-start_time_mono;

		}
		else if ((mode == 2) || (mode == 3)) {
			auto start_time_mono = std::chrono::high_resolution_clock::now();
			convolveBlockResampleFIR(audio_filt, audio_allpass, audio_coeff, my_states.state_mono, \
				audio_decim, audio_upsamp, false);
			auto stop_time_mono = std::chrono::high_resolution_clock::now();
			mono_stereo_runtime += stop_time_mono-start_time_mono;
		}
		//std::cout << " audio_filt = " << audio_filt.size() << std::endl;

		//PLL
		std::vector<float> PLL;
		auto start_time_PLL = std::chrono::high_resolution_clock::now();
		fmPLL(carrier_filt, PLL, my_states.state_PLL, 19e3, if_fs, 2.0, 0.0, 0.01);
		auto stop_time_PLL = std::chrono::high_resolution_clock::now();
		PLL_runtime += stop_time_PLL-start_time_PLL;
		//std::cout << " PLL size = " << PLL.size() << std::endl;

		//mixer
		std::vector<float> mixer;
		mixer.resize(stereo_filt.size(), 0.0);
		auto start_time_mixer = std::chrono::high_resolution_clock::now();
		for (unsigned int z=0; z<mixer.size(); z++){
			mixer[z] = stereo_filt[z]*PLL[z]*2;
		}
		auto stop_time_mixer = std::chrono::high_resolution_clock::now();
		mixer_runtime += stop_time_mixer-start_time_mixer;
		//std::cout << " mixer size = " << mixer.size() << std::endl;

		//filter after mixing
		std::vector<float> stereo_final;
		if((mode == 0) || (mode == 1)) {
			auto start_time_stereo = std::chrono::high_resolution_clock::now();
			convolveBlockFastFIR(stereo_final, mixer, audio_coeff, my_states.state_stereofilt, \
			audio_decim, false);
			auto stop_time_stereo = std::chrono::high_resolution_clock::now();
			stereo_runtime += stop_time_stereo-start_time_stereo;
		}
		else if ((mode == 2) || (mode == 3)) {
			auto start_time_stereo = std::chrono::high_resolution_clock::now();
			convolveBlockResampleFIR(stereo_final, mixer, audio_coeff, \
			my_states.state_stereofilt, audio_decim, audio_upsamp, false);
			auto stop_time_stereo = std::chrono::high_resolution_clock::now();
			stereo_runtime += stop_time_stereo-start_time_stereo;
		}

		//std::cout << " stereo_final size = " << stereo_final.size() << std::endl;
		//std::cout << " CONSUMED FmDemod\n";
		//combination
		std::vector<float> audiodata_L;
		std::vector<float> audiodata_R;
		audiodata_L.resize(stereo_final.size(), 0.0);
		audiodata_R.resize(stereo_final.size(), 0.0);
		for (unsigned int s=0; s<audiodata_L.size(); s++){
			audiodata_L[s] = stereo_final[s]+audio_filt[s];
			audiodata_R[s] = audio_filt[s]-stereo_final[s];
		}
		//allpass_runtime, channel_recov_runtime, carrier_recov_runtime, \
			//PLL_runtime, mixer_runtime, mono_stereo_runtime, stereo_runtime;

		std::cout << " allpass_runtime = " << allpass_runtime.count() << " milliseconds" << "\n";
		std::cout << " channel_recov_runtime = " << channel_recov_runtime.count() << " milliseconds" << "\n";
		std::cout << " carrier_recov_runtime = " << carrier_recov_runtime.count() << " milliseconds" << "\n";
		std::cout << " PLL_runtime = " << PLL_runtime.count() << " milliseconds" << "\n";
		std::cout << " mixer_runtime = " << mixer_runtime.count() << " milliseconds" << "\n";
		std::cout << " mono_stereo_runtime = " << mono_stereo_runtime.count() << " milliseconds" << "\n";
		std::cout << " stereo_runtime = " << stereo_runtime.count() << " milliseconds" << "\n";

		/*
		//Convert audio data to 16-bits and write to standard output
		std::vector<short int> wav_data(audiodata_L.size()*2);
		for ( unsigned int k=0; k<audiodata_L.size(); k++){
			if (std::isnan(audiodata_L[k])) wav_data[k*2] = 0;
			else wav_data[k*2] = static_cast<short int> (audiodata_L[k] * 16384);

			if (std::isnan(audiodata_R[k])) wav_data[k*2+1] = 0;
			else wav_data[k*2+1] = static_cast<short int> (audiodata_R[k] * 16384);
		}
		fwrite(&wav_data[0], sizeof(short int), wav_data.size(), stdout);*/

    // check the end cond (custom) to terminate the thread
    if (((std::cin.rdstate()) != 0) && my_queue.empty()){
      break;
    }
  }
}

void RF_MONO(std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, \
  std::condition_variable& my_cvar, struct PARAMS &audio_params, struct STATES &my_states)
{
	//auto cpu_id = sched_getcpu();
	//std::cout << " Starting MONO thread " << std::this_thread::get_id() << " processor " << cpu_id << std::endl;
	int audio_Fc = 16000;
	int mode = audio_params.mode;

	// coefficients for the filter to extract mono audio
	std::vector<float> audio_coeff;
	if ( (mode == 0) || (mode == 1) ) { impulseResponseLPF(audio_params.if_fs, audio_Fc, audio_params.audio_taps, audio_coeff); }
	else if ( (mode == 2) || (mode == 3) ) { impulseResponseLPF(audio_params.if_fs*audio_params.audio_upsamp, \
		audio_Fc, audio_params.audio_taps, audio_coeff); }

	std::chrono::duration<double, std::milli> mono_runtime;

  while(1)
  {
		//auto old_id = cpu_id;
    // extract data from the queue using mutexes and cond variables as a sync mechanism
    std::unique_lock<std::mutex> my_lock(my_mutex);
    while(my_queue.empty()){
      my_cvar.wait(my_lock);
    }

    std::vector<float> fm_demod = my_queue.front();
    my_queue.pop();
    my_cvar.notify_one();
    my_lock.unlock();

    //Extracting Mono Audio
		// extract the mono audio data through filtering
		std::vector<float> audio_filt;
		if((mode == 0) || (mode == 1)) {
			auto start_time_mono = std::chrono::high_resolution_clock::now();
			convolveBlockFastFIR(audio_filt, fm_demod, audio_coeff, my_states.state_mono, \
			audio_params.audio_decim, false);
			auto stop_time_mono = std::chrono::high_resolution_clock::now();
			mono_runtime += stop_time_mono-start_time_mono;
		}
		else if ((mode == 2) || (mode == 3)) {
			auto start_time_mono = std::chrono::high_resolution_clock::now();
			convolveBlockResampleFIR(audio_filt, fm_demod, audio_coeff, \
			my_states.state_mono, audio_params.audio_decim, audio_params.audio_upsamp, false);
			auto stop_time_mono = std::chrono::high_resolution_clock::now();
			mono_runtime += stop_time_mono-start_time_mono;
		}
		//std::cout << " CONSUMED FmDemod\n";
		//std::cout << " audio_filt size = " << audio_filt.size() << std::endl;
		std::cout << "Mono conv ran for " << mono_runtime.count() << " milliseconds" << "\n";

		/*
		//Convert audio data to 16-bits and write to standard output
		std::vector<short int> wav_data(audio_filt.size());
		for ( unsigned int k=0; k<audio_filt.size(); k++){
			if (std::isnan(audio_filt[k])) wav_data[k] = 0;
			else wav_data[k] = static_cast<short int> (audio_filt[k] * 16384);
		}
		fwrite(&wav_data[0], sizeof(short int), wav_data.size(), stdout);*/

		/*
		cpu_id = sched_getcpu();
		if(cpu_id != old_id){
			std::cout << " changed mono thread " << std::this_thread::get_id() <<  " on processor " << cpu_id << std::endl;
		}*/

    // terminate if queue is empty and no input on
    if (((std::cin.rdstate()) != 0) && (my_queue.empty())){
      break;
    }
  }
}

//////////////////////////////////////////////// MAIN ///////////////////////////////////////////////
int main(int argc, char* argv[])
{
	int mode = 0;
	int audio_channels = 1;
	//Set operating mode based on cmd line arguments
	if(argc < 3){
		std::cout << "Operating in default mode 0: 1 Channel(Mono)" << std::endl;
	}
	else if (argc == 3){
		// Parameter for mode
		mode = atoi(argv[1]);
		if(mode > 3){
			std::cout << "Wrong mode " << mode << std::endl;
			exit(1);
		}
		else {
			std::cout << "Operating in mode:" << mode << std::endl;
		}
		// Parameter for number of channels; Mono(1) or stereo(2)
		audio_channels = atoi(argv[2]);
		if(audio_channels > 2){
			std::cout << "Wrong mode " << mode << std::endl;
			exit(1);
		}
		else {
			std::cout << " Number of channels set to: " << audio_channels << std::endl;
		}
	}
	else {
		std::cout << "Usage: " << argv[0] << std::endl;
		std::cout << "or " << std::endl;
		std::cout << "Usage: " << argv[0] << " <mode>" << std::endl;
		std::cout << "\t\t <mode> is a value from 0 to 3" << std::endl;
		exit(1);
	}

  int rf_Fs, if_fs, rf_decim, audio_decim, audio_upsamp, audio_taps;	//SHOULD WE MAKE THEM ALL FLOAT?
	float audio_Fs;
	//Set Fs for each processing stage and mode
	if (mode == 1){ rf_Fs = 1440000; if_fs = 288000; audio_Fs = 48000.0; rf_decim = 5; audio_decim = 6; audio_upsamp = 0; audio_taps = 13;}
	else if (mode == 2){ rf_Fs = 2400000; if_fs = 240000; audio_Fs = 44100.0; rf_decim = 10; audio_decim = 800; audio_upsamp = 147; audio_taps = 13*audio_upsamp;}
	else if (mode == 3){ rf_Fs = 960000; if_fs = 320000; audio_Fs = 44100.0; rf_decim = 3; audio_decim = 3200; audio_upsamp = 441; audio_taps = 13*audio_upsamp;}
	else { rf_Fs = 2400000; if_fs = 240000; audio_Fs = 48000.0; rf_decim = 10; audio_decim = 5; audio_upsamp = 0; audio_taps = 13;}

	int stereo_taps = 13;

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
		audio_taps,
		stereo_taps
	};

  // states for Mono processing
	std::vector<float> state_mono;
	state_mono.resize(audio_taps-1, static_cast<float>(0));

	//states for stereo processing
	std::vector<float> state_stereo;
	std::vector<float> state_carrier;
	std::vector<float> state_stereofilt;
	std::vector<float> state_allpass;
	state_stereo.resize(stereo_taps-1, static_cast<float>(0));
	state_carrier.resize(stereo_taps-1, static_cast<float>(0));
	state_stereofilt.resize(audio_taps-1, static_cast<float>(0));
	state_allpass.resize(int((stereo_taps-1)/2), static_cast<float>(0));
	std::vector<float> state_PLL{0.0,0.0,1.0,0.0,1.0,0};

	struct STATES my_states;
	my_states = STATES{
	  state_mono,
		state_stereo,
		state_carrier,
		state_stereofilt,
		state_allpass,
		state_PLL
	};

  // context used for thread synch through queue
  std::queue<std::vector<float>> my_queue;
  std::mutex my_mutex;
  std::condition_variable my_cvar;

  // producer thread -- RF Frond-end
  std::thread t_FRE = std::thread(RF_FrontEnd, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar), \
		std::ref(my_params));

	/*
	std::thread t_audio = std::thread(RF_STEREO, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar), \
			std::ref(my_params), std::ref(my_states));*/

	std::thread t_audio;
	if(audio_channels == 1){
		// consumer thread1 -- Mono
		t_audio = std::thread(RF_MONO, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar), \
			std::ref(my_params), std::ref(my_states));
	}
	else if(audio_channels == 2){
		// consumer thread1 -- Stereo
		t_audio = std::thread(RF_STEREO, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar), \
			std::ref(my_params), std::ref(my_states));
	}

  t_FRE.join();
  t_audio.join();

  return 0;

}
