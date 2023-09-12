#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
import matplotlib.pyplot as plt2
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
import sys


from fmSupportLib import *
from fmMonoBasic import impResponse

# flag that keeps track if your code is running for
il_vs_th = 1

if __name__ == "__main__":

    rf_Fs, if_fs, rf_decim, audio_decim, audio_upsamp, audio_taps = 0,0,0,0,0,0
    audio_Fs = 0.0
    rds_resampler_taps = 0
    rds_upsamp = 0
    rds_decim = 0
    rds_taps = 0
    rds_SPS = 0
    rds_RRC_taps = 0
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

    if (mode == 0):
      rf_Fs,if_fs,audio_Fs,rf_decim,audio_decim,audio_upsamp,audio_taps = 2400000,240000,48000.0,10,5,0,101
      rds_upsamp = 247
      rds_decim = 960
      rds_resampler_taps = 101*rds_upsamp
      rds_SPS = 26
      rds_RRC_taps = 101
    elif (mode == 1):
      rf_Fs,if_fs,audio_Fs,rf_decim,audio_decim,audio_upsamp,audio_taps = 1440000,288000,48000.0,5,6,0,101
    elif (mode == 2):
      rf_Fs,if_fs,audio_Fs,rf_decim,audio_decim,audio_upsamp = 2400000,240000,44100.0,10,800,147
      audio_taps = 101*audio_upsamp
      rds_upsamp = 817
      rds_decim = 1920
      rds_resampler_taps = 101*rds_upsamp
      rds_SPS = 43
      rds_RRC_taps = 101
    elif (mode == 3):
      rf_Fs,if_fs,audio_Fs,rf_decim,audio_decim,audio_upsamp= 960000,320000,44100.0,3,3200,441
      audio_taps = 101*audio_upsamp

    print("________________INPUT PARAMS_______________")
    # print ("rf_fs = " + str(rf_Fs))
    # print ("if_fs = " + str(if_fs))
    # print ("audio_fs = " + str(audio_Fs))
    # print ("rf_decim = " + str(rf_decim))
    # print ("audio decim = " + str(audio_decim))
    # print ("audio_upsamp = " + str(audio_upsamp))
    # print ("audio taps = " + str(audio_taps))
    print("rds_upsamp = " + str(rds_upsamp))
    print("rds_decim = " + str(rds_decim))
    print("rds_taps = " + str(rds_taps))
    print("rds_SPS = " + str(rds_SPS))
    print("rds_RRC_taps = " + str(rds_RRC_taps))
    print("__________________________________________\n")


    #Assigning constant global variables
    rf_Fc = 100000
    rf_taps = 151
    audio_Fc = 16000

    #Stereo
    stereo_taps = 151

    rds_taps = 151
    # read the raw IQ data from the recorded file
    # IQ data is assumed to be in 8-bits unsigned (and interleaved)
    in_fname = "../data/lab3_iq_samples/samples8.raw"
    raw_data = np.fromfile(in_fname, dtype='uint8')
    print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
    # IQ data is normalized between -1 and +1 in 32-bit float format
    iq_data = (np.float32(raw_data) - 128.0)/128.0
    print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

    # coefficients for the front-end low-pass filter
    #rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))
    rf_coeff = impResponse(rf_taps, rf_Fs, rf_Fc)

    # coefficients for the filters
    if il_vs_th == 0:
      audio_coeff = signal.firwin(audio_taps, audio_Fc/(240e3/2), window=('hann'))
    else:
      # Stereo coefficients
      carrier_coeff = bandPass(stereo_taps, if_fs, 18.5e3, 19.5e3)
      stereo_coeff = bandPass(stereo_taps, if_fs, 22e3, 54e3)
      # RDS coefficients
      rds_channel_coeff = bandPass(rds_taps, if_fs, 54e3, 60e3)
      rds_carrier_coeff = bandPass(rds_taps, if_fs, 113.5e3, 114.5e3)
      rds_resampler_coeff = impResponse(rds_resampler_taps, if_fs*rds_upsamp, 3e3)
      rds_RRC_coeff = impulseResponseRootRaisedCosine(2375*rds_SPS, rds_RRC_taps)	# XX Fs?
      # Mono coefficients
      if (mode == 0 or mode == 1):
        audio_coeff = impResponse(audio_taps, if_fs, audio_Fc)
      elif (mode == 2 or mode == 3):
        audio_coeff = impResponse(audio_taps, if_fs*audio_upsamp, audio_Fc)


    # set up the subfigures for plotting
    subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
    plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
    fig.subplots_adjust(hspace = .6)

    # IQ constellation figure
    plt2.rc('figure', figsize=(6, 5))	# the size of the entire figure
    fig2, ax_fig2 = plt2.subplots()
    ax_fig2.set_title("IQ Constellation Graph")


    # select a block_size that is a multiple of KB
    # and a multiple of decimation factors
    if (mode == 0 or mode == 1):
      #block_size =  1024 * rf_decim * audio_decim *  2
      block_size =  2 * rf_decim * audio_decim * rds_decim * 2
    elif (mode == 2):
      #block_size = 7 * rf_decim * audio_decim * 2
      block_size = rf_decim * audio_decim * rds_decim * 2
    elif (mode == 3):
      block_size = 3 * rf_decim * audio_decim * 2

    print("Block size =",block_size)

    # States needed for continuity in block processing
    #	Front-end states
    state_i_lpf_100k = np.zeros(rf_taps-1)
    state_q_lpf_100k = np.zeros(rf_taps-1)
    state_phase = 0
    state_IQ = np.zeros(2)
    #	Mono path states
    state_mono = np.zeros(audio_taps-1)
    #	Stereo Path states
    state_stereo = np.zeros(stereo_taps-1)
    state_carrier = np.zeros(stereo_taps-1)
    state_stereofilt = np.zeros(audio_taps-1)
    state_allpass = np.zeros(int((stereo_taps-1)/2))
    state_Pll = [0.0, 0.0, 1.0, 0.0, 1.0, 0]
    #	RDS states
    state_rds_pll = [0.0, 0.0, 1.0, 0.0, 1.0, 0, 1.0]
    state_rds_channel = np.zeros(rds_taps-1)
    state_rds_carrier = np.zeros(rds_taps-1)
    state_rds_allpass = np.zeros(int((rds_taps-1)/2))
    state_rds_resampler = np.zeros(rds_resampler_taps-1)
    state_rds_resampler2 = np.zeros(rds_resampler_taps-1)
    state_rds_RRC = np.zeros(rds_RRC_taps-1)
    state_rds_RRC2 = np.zeros(rds_RRC_taps-1)

    # audio buffer that stores all the audio blocks
    #	Mono data array
    mono_data = np.array([])
    #	stereo data arrays
    stereo_data = np.array([])
    audio_data_left = np.array([]) # used to concatenate filtered blocks (audio data)
    audio_data_right = np.array([])

    printData = False
    test_block = 10
    block_count = 0

    decoded_data = np.array([])

    # if the number of samples in the last block is less than the block size
    # it is fine to ignore the last few samples from the raw IQ file
    while (block_count+1)*block_size < len(iq_data):
      print('________________ Processing block ' + str(block_count) + "______________")

      # filter to extract the FM channel (I samples are even, Q samples are odd)
      i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
          iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
          zi=state_i_lpf_100k)
      q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
          iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
          zi=state_q_lpf_100k)

      # downsample the I/Q data from the FM channel
      i_ds = i_filt[::rf_decim]
      q_ds = q_filt[::rf_decim]
      #print("  i_ds = ", i_ds, end='\n')

      # FM demodulator
      if il_vs_th == 0:
        fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)
      else:
        fm_demod, state_IQ = compEffDemod(i_ds, q_ds, state_IQ)

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Front-End FINISHED ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      #RDS Channel Extraction
      rds_channel_filt, state_rds_channel = signal.lfilter(rds_channel_coeff, 1.0, fm_demod, zi=state_rds_channel)

      #RDS Carrier Recovery Block
          #	RDS all pass filter
      rds_allpass, state_rds_allpass = allPass(rds_channel_filt, state_rds_allpass)

      #	Squaring Non-linearity
      rds_sqr_nonlin = rds_channel_filt * rds_channel_filt

      #	Bandpass filter
      rds_carrier_filt, state_rds_carrier = signal.lfilter(rds_carrier_coeff, 1.0, rds_sqr_nonlin, zi=state_rds_carrier)

      #	RDS PLL
      rds_PLL, rds_PLL_Q, state_rds_pll = fmPll(rds_carrier_filt, 114e3, if_fs, state_rds_pll, ncoScale=0.5, \
        phaseAdjust=(3*math.pi/8), normBandwidth=0.002)

      #RDS Demodulation
      #	Mixer
      rds_mixer = rds_PLL[:-1] * rds_allpass * 2

      #	Rational Resampler
      rds_resampler_out, state_rds_resampler = convolveBlockResampleFIR(rds_mixer, rds_resampler_coeff, \
        state_rds_resampler, rds_decim, rds_upsamp)

      # Root-raised Cosine Filter
      rds_RRC_out, state_rds_RRC = signal.lfilter(rds_RRC_coeff, 1.0, rds_resampler_out, zi=state_rds_RRC)

      #Quad-component - FOR DEBUGGING
      rds_mixer2 = rds_PLL_Q[:-1] * rds_allpass * 2
      rds_resampler_out2, state_rds_resampler2 = convolveBlockResampleFIR(rds_mixer2, rds_resampler_coeff, \
        state_rds_resampler2, rds_decim, rds_upsamp)
      rds_RRC_out2, state_rds_RRC2 = signal.lfilter(rds_RRC_coeff, 1.0, rds_resampler_out2, zi=state_rds_RRC2)

      # Clock and Data Recovery
      pair = np.zeros(2)
      prev_size = 0
      start_init = 158
      to_pass_on_state = [pair, start_init, prev_size]

      #	Sample the data points
      sampling_point_array = np.zeros(len(rds_RRC_out))
      sampling_point_array2 = np.zeros(len(rds_RRC_out2))
      for i in range(start_init, len(sampling_point_array), rds_SPS):
        sampling_point_array[i] = rds_RRC_out[i]
        sampling_point_array2[i] = rds_RRC_out2[i]

      # 	Compute the CDR data and the manchester decoded data
      samples, to_pass_on_state = CDR(rds_RRC_out, rds_SPS, to_pass_on_state, block_count)

      #	Differential decoding of manchesterd data
      diff_decod = diff_decoding(samples)
      decoded_data = np.concatenate((decoded_data, diff_decod))

      #	Frame Synchronize of decoded data
      offset_type, state_rds_frame =  framesync(decoded_data)
      decoded_data = decoded_data[state_rds_frame:]	#cutdown all used data

      print("offset_type = ", offset_type)

      # Graphing
      if block_count >= 20 and block_count < 21:

        # Plot PSD of Front-end Data
        ax0.clear()
        fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
            'Demodulated FM (block ' + str(block_count) + ')')
        # output binary file name (where samples are written from Python)
        fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
        # create binary file where each sample is a 32-bit float
        fm_demod.astype('float32').tofile(fm_demod_fname)

        # Plot Resampler output, RRC output, and the sampling points
        ax1.plot(rds_RRC_out)
        ax1.plot(rds_resampler_out)
        ax1.plot(sampling_point_array, linestyle='-')
        ax1.set_xlim([-10, 1000])

        # PLot the I and Q data
        ax2.plot(rds_RRC_out)
        ax2.plot(rds_RRC_out2)
        ax2.set_xlim([-10, 1000])

        # Plot the IQ Constellation Graph
        ax_fig2.scatter(sampling_point_array, sampling_point_array2, s=10)
        ax_fig2.set_ylim([-1.25,1.25])
        #fig.savefig("../data/rds" + str(block_count) + ".png")

      block_count += 1

    #uncomment assuming you wish to show some plots
    plt.show()
