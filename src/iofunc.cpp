/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "iofunc.h"

// some basic functions for printing information from vectors
// or to read from/write to binary files in 32-bit float format
void printRealVector(const std::vector<float> &x)
{
	std::cerr << "Printing float vector of size " << x.size() << "\n";
	for (unsigned int i = 0; i < x.size(); i++)
		std::cerr << x[i] << " ";
	std::cerr << "\n";
}

void printComplexVector(const std::vector<std::complex<float>> &X)
{
	std::cerr << "Printing complex vector of size " << X.size() << "\n";
	for (unsigned int i = 0; i < X.size(); i++)
		std::cerr << X[i] << " ";
	std::cerr << "\n";
}

// assumes data in the raw binary file is in 32-bit float format
void readBinData(const std::string in_fname, std::vector<float> &bin_data)
{
	std::ifstream fdin(in_fname, std::ios::binary);
	if(!fdin) {
		std::cerr << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	} else {
		std::cerr << "Reading raw binary from \"" << in_fname << "\"\n";
	}
	fdin.seekg(0, std::ios::end);
	const unsigned int num_samples = fdin.tellg() / sizeof(float);

	bin_data.resize(num_samples);
	fdin.seekg(0, std::ios::beg);
	fdin.read(reinterpret_cast<char*>(&bin_data[0]), num_samples*sizeof(float));
	fdin.close();
}

// assumes data in the raw binary file is 32-bit float format
void writeBinData(const std::string out_fname, const std::vector<float> &bin_data)
{
	std::cerr << "Writing raw binary to \"" << out_fname << "\"\n";
	std::ofstream fdout(out_fname, std::ios::binary);
	for (unsigned int i=0; i<bin_data.size(); i++) {
		fdout.write(reinterpret_cast<const char*>(&bin_data[i]),\
								sizeof(bin_data[i]));
	}
	fdout.close();
}

// function to write audio data to a binary file that contains raw samples
// represented as 32-bit floats; we also assume two audio channels
// note: check the python script that can read this type of files
// and then reformat them to .wav files to be run on third-party players
void write_audio_data(const std::string out_fname, const std::vector<float> &audio_left, const std::vector<float> &audio_right)
{
	// file descriptor for the output to be written
	if (audio_left.size() != audio_right.size()) {
		std::cerr << "Something got messed up with audio channels\n";
		std::cerr << "They must have the same size ... exiting\n";
		exit(1);
	} else {
		std::cerr << "Writing raw audio to \"" << out_fname << "\"\n";
	}
	std::ofstream fdout(out_fname, std::ios::binary);
	for (unsigned int i=0; i<audio_left.size(); i++) {
		// we assume we have handled a stereo audio file
		// hence, we must interleave the two channels
		// (change as needed if testing with mono files)
		fdout.write(reinterpret_cast<const char*>(&audio_left[i]),\
								sizeof(audio_left[i]));
		fdout.write(reinterpret_cast<const char*>(&audio_right[i]),\
								sizeof(audio_right[i]));
	}
	fdout.close();
}

//This function prints the first and last 3 values within a given vector
void printFirstLast(const std::vector<float> &vector, const char* vector_name)
{
	std::cerr << vector_name << " = ["<< vector[0];
	std::cerr << ", "<< vector[1];
	std::cerr << ", "<< vector[2];
	std::cerr << ", ... , "<< vector[vector.size()-3];
	std::cerr << ", "<< vector[vector.size()-2];
	std::cerr << ", "<< vector[vector.size()-1] << " ]" << std::endl;
}

// This function prints last N values in a vector
void printLastN(const std::vector<float> &vector, const unsigned int N, const char* vector_name)
{
	if(vector.size() > N){
		std::cerr << vector_name << " = [ ... , ";
		for (unsigned int i=vector.size()-N; i<vector.size()-1; i++){
			std::cerr << vector[i] << ", ";
		}
		std::cerr << vector[vector.size()-1] << " ]" << std::endl;
	}
	else { std::cerr << "Vector size < N" << std::endl; }
}

//Prints firstN values from the beginning of vector and lastN from the end
void printFirstLastN(const std::vector<float> &vector, const unsigned int firstN, const unsigned int lastN, const char* vector_name)
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

//
//data in the redirected input is interpreted as bytes
void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data){
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for(int k=0; k<(int)num_samples; k++){
		//automaticallt normalizes the data to the range - 1 to +1
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
	}
}
