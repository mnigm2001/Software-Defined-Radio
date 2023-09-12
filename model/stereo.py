#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
import sys
import time

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD, compEffDemod, lpBlockFilter, convolveBlockResampleFIR, convolveBlockFastFIR, fmPll,bandPass, allPass
# for take-home add your functions
from fmMonoBasic import impResponse

#rf_Fs = 2.4e6
#rf_Fc = 100e3
#rf_taps = 151
#rf_decim = 10

#audio_Fs = 48e3
#audio_decim = 5
# add other settings for audio, like filter taps, ...
#audio_taps = 101
#audio_Fc = 16e3

# flag that keeps track if your code is running for
# in-lab (il_vs_th = 0) vs takehome (il_vs_th = 1)
il_vs_th = 1


if __name__ == "__main__":

	#Mode assignment
	mode = 0
	if (int(sys.argv[1]) > 0):
		mode = int(sys.argv[1])
		if (mode == 1):
			print("Operating in mode 1")
		elif (mode == 2):
			print("Operating in mode 2")
		elif (mode == 3):
			print("Operating in mode 3")
		else:
			print("mode is value from 0 to 3")
			exit(0)
	else:
		print("Operating in default mode 0")

	print ("mode = " + str(mode))

	# Set parameters based on mode
	rf_Fs, if_fs, rf_decim, audio_decim, audio_upsamp, audio_taps = 0,0,0,0,0,0
	audio_Fs = 0.0
	if (mode == 0):
		rf_Fs,if_fs,audio_Fs,rf_decim,audio_decim,audio_upsamp,audio_taps = 2400000,240000,48000.0,10,5,0,101
	elif (mode == 1):
		rf_Fs,if_fs,audio_Fs,rf_decim,audio_decim,audio_upsamp,audio_taps = 1440000,288000,48000.0,5,6,0,101
	elif (mode == 2):
		rf_Fs,if_fs,audio_Fs,rf_decim,audio_decim,audio_upsamp = 2400000,240000,44100.0,10,800,147
		audio_taps=101*audio_upsamp
	elif (mode == 3):
		rf_Fs,if_fs,audio_Fs,rf_decim,audio_decim,audio_upsamp= 960000,320000,44100.0,3,3200,441
		audio_taps = 101*audio_upsamp

	#Assigning constant global variables
	rf_Fc = 100000;
	rf_taps = 151;
	audio_Fc = 16000;
	#Stereo
	stereo_taps = 151;

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/lab3_iq_samples/stereo_l0_r9.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

	# coefficients for the filter to extract mono audio
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on firwin
		# same principle  as for rf_coeff (but different arguments, of course)
		audio_coeff = signal.firwin(audio_taps, audio_Fc/(240e3/2), window=('hann'))

	else:
		carrier_coeff = bandPass(stereo_taps,if_fs,18.5e3,19.5e3)
		stereo_coeff = bandPass(stereo_taps,if_fs,22e3,54e3)
		# Generate impulse response
		if (mode == 0 or mode == 1):
			audio_coeff = impResponse(audio_taps,if_fs, audio_Fc)
		elif (mode == 2 or mode == 3):
			audio_coeff = impResponse(audio_taps,if_fs*audio_upsamp, audio_Fc)

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB and a multiple of decimation factors
	if (mode == 0 or mode == 1):
		block_size =  1024 * rf_decim * audio_decim * 2
	elif (mode == 2):
		block_size = 7 * rf_decim * audio_decim * 2
	elif (mode == 3):
		block_size = 3 * rf_decim * audio_decim * 2
	print("Block size =",block_size)

	# states needed for continuity in block processing
	# 	Front-end states
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_phase = 0				#state for atan fmdemod
	state_IQ = np.zeros(2)		#state for compEff fmdemod
	#	Mono path states
	state_mono = np.zeros(audio_taps-1)
	#	Stereo Path states
	state_stereo = np.zeros(stereo_taps-1)
	state_carrier = np.zeros(stereo_taps-1)
	state_stereofilt = np.zeros(audio_taps-1)
	state_allpass = np.zeros(int((stereo_taps-1)/2))
	state_Pll=[0.0,0.0,1.0,0.0,1.0,0]

	# audio buffers that stores all the audio blocks
	#	Mono data arrays
	mono_data = np.array([])
	#	stereo data arrays
	stereo_data = np.array([])
	audio_data_left = np.array([]) # used to concatenate filtered blocks (audio data)
	audio_data_right = np.array([])

	printData = False
	test_block = 10
	block_count = 0

	start = time.time()
	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count+1)*block_size < len(iq_data):
		print('Processing block ' + str(block_count))
		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
		'''
		if(block_count == test_block):
			print('Processing block ' + str(block_count))
			printData = False
		else:
			printData = False
		'''

		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
				zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
				zi=state_q_lpf_100k)
		'''
		i_filt, state_i_lpf_100k = convolveBlockFastFIR(iq_data[(block_count)*block_size:(block_count+1)*block_size:2], rf_coeff, state_i_lpf_100k, rf_decim)
		q_filt, state_q_lpf_100k = convolveBlockFastFIR(iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2], rf_coeff, state_q_lpf_100k, rf_decim)
		if(block_count == test_block):
			print("Filtered and Downsampled FR-end ", i_filt, end='\n')
		'''

		# downsample the I/Q data from the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		if(block_count == test_block):
			print("Down Sampled FR-end", i_ds, end='\n')

		# FM demodulator
		if il_vs_th == 0:
			#fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)
			fm_demod, state_IQ = compEffDemod(i_ds, q_ds, state_IQ)
		else:
			fm_demod, state_IQ = compEffDemod(i_ds, q_ds, state_IQ)

		if(block_count == test_block):
			print("fmDemod", fm_demod, end='\n')

		# extract the stereo audio data through filtering
		if il_vs_th == 0:
			audio_filt, state_mono = signal.lfilter(audio_coeff, 1.0, fm_demod,zi=state_mono)
		else:
            #all pass filter
			audio_allpass, state_allpass = allPass(fm_demod,state_allpass)

			#filter stereo
			stereo_filt, state_stereo = signal.lfilter(stereo_coeff, 1.0, fm_demod, zi=state_stereo)
			carrier_filt, state_carrier = signal.lfilter(carrier_coeff, 1.0, fm_demod, zi=state_carrier)
			#stereo_filt, state_stereo = convolveBlockFastFIR(fm_demod, stereo_coeff, state_stereo, 1)
			#carrier_filt, state_carrier = convolveBlockFastFIR(fm_demod, carrier_coeff, state_carrier, 1)

			# with your own code for BLOCK convolution
			if (mode == 0 or mode == 1):
				audio_filt, state_mono = convolveBlockFastFIR(audio_allpass, audio_coeff, state_mono, audio_decim)
			elif (mode == 2 or mode == 3):
				audio_filt, state_mono = convolveBlockResampleFIR(audio_allpass, audio_coeff, state_mono, audio_decim, audio_upsamp)

		if(block_count == test_block):
			print("Filtered Mono Audio", audio_filt, end='\n')
			print("Audio_filt data size =" + str(len(audio_filt)))

		#Stereo PLL
		PLL, state_Pll = fmPll(carrier_filt, 19e3, if_fs, state_Pll)

		#print("PLL START:", PLL[0])
		#print("PLL END:",PLL[-1])
		#print("state:", state_Pll)

		#Stereo Processing Mixer
		mixer = PLL[:-1] * stereo_filt * 2
		if (mode == 0 or mode == 1):
			stereo_final, state_stereofilt = convolveBlockFastFIR(mixer, audio_coeff, state_stereofilt, audio_decim)
		elif (mode == 2 or mode == 3):
			stereo_final, state_stereofilt = convolveBlockResampleFIR(mixer, audio_coeff, state_stereofilt, audio_decim,\
			 audio_upsamp)

		#audio_dataL=stereo_final+audio_filt
		#print("m=",audio_filt[500])
		#print("s=",stereo_final[500])
		#print("l=",audio_dataL[500])
		#audio_dataR=audio_filt-stereo_final
		#print("r=",audio_dataR[500])

		# Extract stereo data
		audio_dataL = stereo_final+audio_filt
		audio_dataR = audio_filt-stereo_final

		# concatenate processed audio blocks
		audio_data_left = np.concatenate((audio_data_left, audio_dataL))
		audio_data_right = np.concatenate((audio_data_right, audio_dataR))

		# to save runtime select the range of blocks to log data
		# this includes both saving binary files as well plotting PSD

		if block_count >= 229 and block_count < 230:

			# plot PSD of selected block after FM demodulation
			ax0.clear()
			fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
					'Demodulated FM (block ' + str(block_count) + ')')
			# output binary file name (where samples are written from Python)
			fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			fm_demod.astype('float32').tofile(fm_demod_fname)

			# plot PSD of selected block after extracting mono audio
			# ... change as needed
			fmPlotPSD(ax1, stereo_filt, (rf_Fs/rf_decim)/1e3, subfig_height[1], 'Stereo Data')
			fmPlotPSD(ax1, carrier_filt, (rf_Fs/rf_decim)/1e3, subfig_height[1], 'Stereo Data')
			# plot PSD of selected block after downsampling mono audio
			# ... change as needed
			#fmPlotPSD(ax2, audio_filt, audio_Fs/1e3, subfig_height[2], 'Downsampled Mono Audio')
			#fmPlotPSD(ax2, stereo_final, audio_Fs/1e3, subfig_height[2], 'Downsampled Mono Audio')
			#fmPlotPSD(ax2, audio_dataL, audio_Fs/1e3, subfig_height[2], 'Downsampled Mono Audio')
			# save figure to file
			fmPlotPSD(ax2, mixer, (rf_Fs/rf_decim)/1e3, subfig_height[1], 'Stereo Data')
			fmPlotPSD(ax2, stereo_final,(rf_Fs/rf_decim)/1e3, subfig_height[1], 'Stereo Data')
			fmPlotPSD(ax2, PLL, (rf_Fs/rf_decim)/1e3, subfig_height[1], 'Stereo Data')
			#plt.plot(stereo_final)
			#plt.plot(audio_filt)
			#plt.plot(PLL)

			fig.savefig("../data/stereo_testing" + str(block_count) + ".png")

		block_count += 1

	end = time.time()
	print("runtime measured: ", end - start, end='\n')

	print('Block count =' + str(block_count))
	print("audio_data_left size = ", len(audio_data_left), end='\n')
	print("audio_data_left = ", audio_data_left, end='\n')
	print("audio_data_right size = ", len(audio_data_right), end='\n')
	print("audio_data_right = ", audio_data_right, end='\n')
	print('Finished processing all the blocks from the recorded I/Q samples')

	#concatenate both left and right channels
	audio_data = np.zeros(len(audio_data_left)+len(audio_data_right))
	print("audio_data size = ", len(audio_data), end='\n')
	for i in range(0,len(audio_data_left)-1):
		audio_data[2*i] = audio_data_left[i]
		audio_data[2*i+1] = audio_data_right[i]


	# write audio data to file
	out_fname = "../data/fmLeftMode2.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data_left/2)*32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	out_fname = "../data/fmRightMode2_2.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data_right/2)*32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	out_fname = "../data/fmRecombined.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	plt.show()
