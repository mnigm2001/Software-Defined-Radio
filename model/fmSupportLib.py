#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math, cmath

# Galois Field Matrix Multiplication
def matrixMult(decoded_data, parity_matrix):
	output = np.empty(10)
	temp_mat = np.empty(len(parity_matrix))
	for k in range(0, len(output)):
		ones_count = 0
		for i in range(0, len(temp_mat)):
			temp_mat[i] = decoded_data[i]*parity_matrix[i][k]
			if(temp_mat[i] == 1):
				ones_count += 1
		if(ones_count%2 == 1):
			output[k] = 1
		else:
			output[k] = 0
	return output

# Frame Synchronization of decoded data
def framesync(diff_data):
	# Parity Matrix
	h=[[1,0,0,0,0,0,0,0,0,0],
	[0,1,0,0,0,0,0,0,0,0],
	[0,0,1,0,0,0,0,0,0,0],
	[0,0,0,1,0,0,0,0,0,0],
	[0,0,0,0,1,0,0,0,0,0],
	[0,0,0,0,0,1,0,0,0,0],
	[0,0,0,0,0,0,1,0,0,0],
	[0,0,0,0,0,0,0,1,0,0],
	[0,0,0,0,0,0,0,0,1,0],
	[0,0,0,0,0,0,0,0,0,1],
	[1,0,1,1,0,1,1,1,0,0],
	[0,1,0,1,1,0,1,1,1,0],
	[0,0,1,0,1,1,0,1,1,1],
	[1,0,1,0,0,0,0,1,1,1],
	[1,1,1,0,0,1,1,1,1,1],
	[1,1,0,0,0,1,0,0,1,1],
	[1,1,0,1,0,1,0,1,0,1],
	[1,1,0,1,1,1,0,1,1,0],
	[0,1,1,0,1,1,1,0,1,1],
	[1,0,0,0,0,0,0,0,0,1],
	[1,1,1,1,0,1,1,1,0,0],
	[0,1,1,1,1,0,1,1,1,0],
	[0,0,1,1,1,1,0,1,1,1],
	[1,0,1,0,1,0,0,1,1,1],
	[1,1,1,0,0,0,1,1,1,1],
	[1,1,0,0,0,1,1,0,1,1]]

	n=0
	offset_type = " "
	while (n<len(diff_data)-26):
		# Compute the Galois Field Matrix Multiplication of the decoded data and the parity matrix
		s = matrixMult(diff_data[n:26+n], h)
		# Check for Syndrome words
		if (s==[1,1,1,1,0,1,1,0,0,0]).all():
			offset_type = 'A'
			# If the remaining data cannot accomidate another syndrome word, break
			if(len(diff_data)-(n+26) < 26):
				break
			n+=26
		elif ( s == [1,1,1,1,0,1,0,1,0,0]).all():
			offset_type = 'B'
			if(len(diff_data)-(n+26) < 26):
				break
			n+=26
		elif (s==[1,0,0,1,0,1,1,1,0,0]).all():
			offset_type = 'C'
			if(len(diff_data)-(n+26) < 26):
				break
			n+=26
		elif(s == [1,1,1,1,0,0,1,1,0,0]).all():
			offset_type = 'C_apos'
			if(len(diff_data)-(n+26) < 26):
				break
			n+=26
		elif (s==[1,0,0,1,0,1,1,0,0,0]).all():
			offset_type = 'D'
			if(len(diff_data)-(n+26) < 26):
				break
			n+=26
		else:
			n+=1

	#Save next index for concatention
	if(offset_type == " "):
		state_index = n
	else:
		state_index = n+26

	return offset_type , state_index


def CDR(input1, rds_SPS, to_pass_on_state, block_count):
	pair = to_pass_on_state[0]#0: value from last block, 1:sampling point from current block
	start = to_pass_on_state[1]
	prev_size = to_pass_on_state[2]
	output = np.array([])
	unpaired = True
	limit = 0.3
	while(unpaired):
		#store sampling points from input
		sampling_point_array = np.zeros(len(input1))
		size=0
		for i in range(start, len(sampling_point_array), rds_SPS):
			#Pairs last sampling point of prev block with first sampling point
			if((i == start) and (start == to_pass_on_state[1]) and (prev_size%2 == 1)):
				pair[1] = input1[i]
				bit = symbolToBit(pair)
				output = np.append(output,bit)
				#shift to the left to prepare for next input
				pair[0] = pair[1]
				#Adjusts the start sampling point
				start = start + rds_SPS
				continue

			#Takes care of 3 consecutive H's or L's
			if((i>= start + 2*rds_SPS) and (sampling_point_array[i-2*rds_SPS]>0) and (sampling_point_array[i-1*rds_SPS]>0) \
				and (input1[i]>0)):
				#print("INV High ",end='\t' )
				sampling_point_array[i] = -1*input1[i]
			elif((i>= start + 2*rds_SPS) and (sampling_point_array[i-2*rds_SPS]<0) and (sampling_point_array[i-1*rds_SPS]<0) \
				and (input1[i]<0)):
				#print("INV Low ",end='\t' )
				sampling_point_array[i] = -1*input1[i]
			else:
				sampling_point_array[i] = input1[i]

			size += 1

		#print("size =", size)
		#print("___________________START VALUE =",start, "__________________\n")
		#Removing zeros
		samples = copyFrom(sampling_point_array, rds_SPS, start, size)

		for i in range(0,len(samples), 2):
			if((i+1)<len(samples)):
				if((samples[i]<0 and samples[i+1]<0)or(samples[i]>0 and samples[i+1]>0)):
					#print("irregular pair detected")
					#checks the limit condition
					if(abs(samples[i])<limit or abs(samples[i+1])<limit):
						if(abs(samples[i])<limit):
							samples[i] = -1*samples[i]
							#print("dealt with by inverting input[",i,"]")
						elif(abs(samples[i+1])<limit):
							samples[i+1] = -1*samples[i+1]
							#print("dealt with by inverting input[",i+1,"]")
					#changes starting position
					else:
						start = start + rds_SPS
						#print("could not invert either samples, restarting from sample ", start, "..")
						if(block_count != 0):
							pair[1] = samples[0]
							bit = symbolToBit(pair)
							#print("pairing ", pair[0], "with ", pair[1], "...")
							#print("output bit = ", bit, "\n")
							#output = np.concatenate((output, bit))
							output = np.append(output,bit)
							pair[0] = pair[1]
						#break_count += 1
						unpaired = True
						break
				else:
					#print("paired=false")
					unpaired = False

	#print("No irregular pairs found")
	#print("size =", size)
	#print(samples)

	#Stores the last sampling point of current block
	pair[0] = samples[-1]
	to_pass_on_state[0] = pair

	#Computes the starting point of next block
	last_index = ((size-1)*rds_SPS) + start
	next_start = rds_SPS - (len(input1) - last_index)
	to_pass_on_state[1] = next_start
	#print("last_index = ", last_index)
	#print("next_start = ", next_start)

	#store current size for next block pairing
	to_pass_on_state[2] = size

	#Updates final manchestered_data
	manchestered_data = manchestering(samples)
	#print("manchestered_data = ", manchestered_data)

	output = np.concatenate((output, manchestered_data))
	#print("Output = ", output)

	return output, to_pass_on_state

def manchestering(input1):
	manchestered_data = np.zeros(int(len(input1)/2))
	#manchestered_data = np.array([])
	last_bit = np.zeros(2)
	last_bit[0] = input1[-1]

	for i in range(0, len(input1), 2):
		if((i+1)<len(input1)):
			if((input1[i]<0) and (input1[i+1]>0)):
				#print("input[",i,"] = ", input1[i], "and input[",i+1,"] = ", input1[i+1]," ADDING ZERO")
				#output = np.concatenate((output, 0))
				manchestered_data[int(i/2)] = 0
			elif((input1[i]>0) and (input1[i+1]<0)):
				#output = np.concatenate((output, 1))
				#print("input[",i,"] = ", input1[i], "and input[",i+1,"] = ", input1[i+1]," ADDING ONE")
				manchestered_data[int(i/2)] = 1
	#print("manchestered_data size = ", len(manchestered_data))
	return manchestered_data

def symbolToBit(pair):
	bit = 0
	#In the case of LH or LL
	if(pair[0] < 0):
		bit = 0
	#In the case of HL or HH
	elif(pair[0] > 0):
		bit = 1
	return bit

# This functions removes interleaved zeros within an array
def copyFrom(sampling_point_array, rds_SPS, start, size):
	output = np.zeros(size)
	for i in range(start,len(sampling_point_array), rds_SPS):
		output[int((i-start)/rds_SPS)] = sampling_point_array[i]
		#print("for i = ",i,"output[",int((i-start)/rds_SPS),"] = ", sampling_point_array[i])
	return output

# This function computes the differentially decoded data
def diff_decoding(manch_data):
	diff_data = np.empty(len(manch_data))
	diff_data[0] = manch_data[0]
	for i in range(1, len(manch_data)):
		if(manch_data[i] != manch_data[i-1]):
			diff_data[i] = 1
		elif(manch_data[i] == manch_data[i-1]):
			diff_data[i] = 0
	return diff_data

def impulseResponseRootRaisedCosine(Fs, N_taps):

	"""
	Root raised cosine (RRC) filter

	Fs  		sampling rate at the output of the resampler in the RDS path
				sampling rate must be an integer multipler of 2375
				this integer multiple is the number of samples per symbol

	N_taps  	number of filter taps

	"""

	# duration for each symbol - should NOT be changed for RDS!
	T_symbol = 1/2375.0

	# roll-off factor (must be greater than 0 and smaller than 1)
	# for RDS a value in the range of 0.9 is a good trade-off between
	# the excess bandwidth and the size/duration of ripples in the time-domain
	beta = 0.90

	# the RRC inpulse response that will be computed in this function
	impulseResponseRRC = np.empty(N_taps)

	for k in range(N_taps):
		t = float((k-N_taps/2))/Fs
		# we ignore the 1/T_symbol scale factor
		if t == 0.0: impulseResponseRRC[k] = 1.0 + beta*((4/math.pi)-1)
		elif t == -T_symbol/(4*beta) or t == T_symbol/(4*beta):
			impulseResponseRRC[k] = (beta/np.sqrt(2))*(((1+2/math.pi)* \
					(math.sin(math.pi/(4*beta)))) + ((1-2/math.pi)*(math.cos(math.pi/(4*beta)))))
		else: impulseResponseRRC[k] = (math.sin(math.pi*t*(1-beta)/T_symbol) +  \
					4*beta*(t/T_symbol)*math.cos(math.pi*t*(1+beta)/T_symbol))/ \
					(math.pi*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol)

	# returns the RRC impulse response to be used by convolution
	return impulseResponseRRC


#All pass filter to match phase delay for split streams
def allPass(input_block, state_block):
    output_block= np.concatenate((state_block,input_block[:-len(state_block)]))
    state_block=input_block[-len(state_block):]

    return output_block,state_block

def fmPll(pllIn, freq, Fs, state, ncoScale = 2.0, phaseAdjust = 0.0, normBandwidth = 0.01,):

	# scale factors for proportional/integrator terms
	# these scale factors were derived assuming the following:
	# damping factor of 0.707 (1 over square root of 2)
	# there is no oscillator gain and no phase detector gain
	Cp = 2.666
	Ci = 3.555

	# gain for the proportional term
	Kp = (normBandwidth)*Cp
	# gain for the integrator term
	Ki = (normBandwidth*normBandwidth)*Ci

	# output array for the NCO
	ncoOut = np.empty(len(pllIn)+1)
	ncoOutQ = np.empty(len(pllIn)+1)

	# initialize internal state
	integrator = state[0]
	phaseEst = state[1]
	feedbackI = state[2]
	feedbackQ = state[3]
	ncoOut[0] = state[4]
	trigOffset = state[5]
	ncoOutQ[0] = state[6]
	# note: state saving will be needed for block processing

	for k in range(len(pllIn)):

		# phase detector
		errorI = pllIn[k] * (+feedbackI)  # complex conjugate of the
		errorQ = pllIn[k] * (-feedbackQ)  # feedback complex exponential

		# four-quadrant arctangent discriminator for phase error detection
		errorD = math.atan2(errorQ, errorI)

		# loop filter
		integrator = integrator + Ki*errorD

		# update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator

		# internal oscillator
		trigOffset += 1
		trigArg = 2*math.pi*(freq/Fs)*(trigOffset) + phaseEst
        #I and Q for the real/imaginary parts of the oscillator
		feedbackI = math.cos(trigArg)
		feedbackQ = math.sin(trigArg)
		ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)
		ncoOutQ[k+1] = math.sin(trigArg*ncoScale + phaseAdjust)

    #Store values for next state, the last ncoOut is stored for next block
	stateout = [integrator, phaseEst, feedbackI, feedbackQ, ncoOut[-1], trigOffset, ncoOutQ[-1]]
	# for stereo only the in-phase NCO component should be returned
	# for block processing you should also return the state
	return ncoOut, ncoOutQ, stateout
	# for RDS add also the quadrature NCO component to the output


#Function to generate coeffcients for a bandpass filter of N_taps with Fb and Fe cutoffs (beginning and ending)
def bandPass(N_taps,Fs, Fb,Fe):
	coeff = np.zeros(N_taps)
	Normcenter = ((Fe+Fb)/2)/(Fs/2)
	Normpass = (Fe-Fb)/(Fs/2)

	for i in range(0,N_taps):
		if (i==(N_taps-1)/2):
			coeff[i] = Normpass
		else:
			coeff[i] = Normpass*((np.sin(np.pi*Normpass/2*(i-(N_taps-1)/2)))/(np.pi*Normpass/2*(i-(N_taps-1)/2)))

		coeff[i]= coeff[i]*np.cos(i*(np.pi)*Normcenter)
		coeff[i] = coeff[i]*(np.sin(i*np.pi/N_taps))**2
	return coeff



#Function to generate impulse response coefficients
def impResponse(N_taps,Fs, Fc):
	coeff = np.zeros(N_taps)
	NormFc = Fc/(Fs/2)
	for i in range(0,N_taps):
		if (i==(N_taps-1)/2):
			coeff[i] = NormFc
		else:
			coeff[i] = NormFc*((np.sin(np.pi*NormFc*(i-(N_taps-1)/2)))/(np.pi*NormFc*(i-(N_taps-1)/2)))
		coeff[i] = coeff[i]*(np.sin(i*np.pi/N_taps))**2
	return coeff

#Function to resample data through upsampling and downsampling
def convolveBlockResampleFIR(x, h, state, audio_decim, audio_upsamp):
	output_filt = np.zeros(int(len(x)*audio_upsamp/audio_decim))
	#print("x size = ", len(x), " h size = ", len(h), " output_filt size = ", len(output_filt), end='\n')
	for m in range(0, len(x)*audio_upsamp, audio_decim):
		phase = m%audio_upsamp
		for n in range(phase, len(h), audio_upsamp):
			if (m-n >= 0):
				#print("m = ", m, " n = ", n, " m-n ", m-n," m/decim = ", m/audio_decim, " m-n/decim = ", (m-n)/audio_upsamp, end='\n')
				output_filt[int(m/audio_decim)] += h[n] * x[int((m-n)/audio_upsamp)]
			else:
				output_filt[int(m/audio_decim)] += h[n] * state[m-n+len(state)]
		output_filt[int(m/audio_decim)] = output_filt[int(m/audio_decim)] * audio_upsamp

	k = audio_upsamp-1
	for i in range(audio_upsamp*len(x)-len(state), audio_upsamp*len(x)-audio_upsamp, audio_upsamp):
		state[k] = x[int((i/audio_upsamp)) + 1]
		k += audio_upsamp

	return output_filt, state

#Function that performs convolution and downsamples by audio_decim
def convolveBlockFastFIR(x, h, state, audio_decim):
	#print("audio_decim =",audio_decim)
	#print("x.size = ",len(x))
	output_filt = np.zeros((int)(np.size(x)/audio_decim))
	#output_filt = np.zeros(1024)
	for m in range(0, len(x), audio_decim):
		for n in range(0, len(h)):
			if (m-n >= 0):
				output_filt[(int)(m/audio_decim)]	+= h[n] * x[m-n]
			else:
				output_filt[(int)(m/audio_decim)] += h[n] * state[m-n+len(state)]

	#prepare the next state
	k = 0
	for i in range(len(x)-len(state), len(x)):
		state[k] = x[i]
		k+=1

	return output_filt, state


#Function to filter data in a single pass
def lpFilter(h, x):
	y = np.zeros(len(x)+len(h)-1)
	for n in range(len(y)):
		for k in range(len(h)):
			if (n-k >= 0) and (n-k < len(x)):
				y[n] += h[k] * x[n-k]
	return y

#This function computes the convolution for a given block of data
def lpBlockFilter(h, x, state, printData):
	y = np.zeros(len(x))

	#if(printData):
	#	print("x = ", x, end='\n')
	#	print("h = ", h, end='\n')
	for n in range(len(y)):
		for k in range(len(h)):

			#if (printData and (n==0) and (k<10)):
			#	print("n = ", n, "\tk = ", k, "\tn-k = ", n-k, end='\n')

			if n-k >= 0:									#Using current block audio data i.e. xb for convolution
				y[n] += h[k] * x[n-k]
			#	if (printData and (n==0) and (k<10)):
			#		print("X\ty = ", y[n], "\th = ", h[k], "\tx = ", x[n-k], end='\n\n')
			else:
				y[n] += h[k] * state[n-k+len(state)]		#Using pervious block audio data for convolution
			#	if (printData and (n==0) and (k<10)):
			#		print("STATE\ty = ", y[n], "\th = ", h[k], "\tstate = ", state[n-k+len(state)], end='\n\n')

	#Prepare state for next block
	state = x[len(x)-len(state):len(x)]
	return y, state

#prev_IQ is the last IQ data from previous block in form [I,Q]
def compEffDemod(I, Q, prev_IQ):

	# empty vector to store the demodulated samples
	fm_demod = np.empty(len(I))
	q_sqr = np.square(Q)
	i_sqr = np.square(I)

	# iterate through each of the I and Q pairs
	for k in range(len(I)):

		#Differentiate both I and Q samples
		i_diff = I[k] - prev_IQ[0]
		q_diff = Q[k] - prev_IQ[1]
		fm_numer = ((I[k]*q_diff) - (Q[k]*i_diff))
		fm_denom = (q_sqr[k] + i_sqr[k])

		#Scale the differences
		if(math.isnan(fm_numer/fm_denom)):
			fm_demod[k] = 0.0
		else:
			fm_demod[k] = fm_numer/fm_denom

		# if(k<=2):
		# 	print("k = ", k, end='\n')
		# 	print("prev_IQ =", prev_IQ, end='\n')
		# 	print("I= ",I[k], "Q= ", Q[k], end='\n')
		# 	print("i_sqr= ",i_sqr[k], "q_sqr= ", q_sqr[k], end='\n')
		# 	print("i_diff= ",i_diff, "q_diff= ", q_diff, end='\n')
		# 	print("fm_demod[k] =", fm_demod[k], end='\n\n')

		#Save the state of the current IQ samples to compute the next
		prev_IQ = [I[k], Q[k]]

	#WE NEED TO RETURN THE LAST I AND Q SINCE THOSE ARE WHAT WE USE IN I[K-1]
	return fm_demod, prev_IQ

def fmDemodArctan(I, Q, prev_phase = 0.0):
#
# the default prev_phase phase is assumed to be zero, however
# take note in block processing it must be explicitly controlled

	# empty vector to store the demodulated samples
	fm_demod = np.empty(len(I))

	# iterate through each of the I and Q pairs
	for k in range(len(I)):

		# use the atan2 function (four quadrant version) to detect angle between
		# the imaginary part (quadrature Q) and the real part (in-phase I)
		current_phase = math.atan2(Q[k], I[k])

		# we need to unwrap the angle obtained in radians through arctan2
		# to deal with the case when the change between consecutive angles
		# is greater than Pi radians (unwrap brings it back between -Pi to Pi)
		[prev_phase, current_phase] = np.unwrap([prev_phase, current_phase])

		# take the derivative of the phase
		fm_demod[k] = current_phase - prev_phase

		# save the state of the current phase
		# to compute the next derivative
		prev_phase = current_phase

	# return both the demodulated samples as well as the last phase
	# (the last phase is needed to enable continuity for block processing)
	return fm_demod, prev_phase

# custom function for DFT that can be used by the PSD estimate
def DFT(x):

	# number of samples
	N = len(x)

	# frequency bins
	Xf = np.zeros(N, dtype='complex')

	# iterate through all frequency bins/samples
	for m in range(N):
		for k in range(N):
			Xf[m] += x[k] * cmath.exp(1j * 2 * math.pi * ((-k) * m) / N)

	# return the vector that holds the frequency bins
	return Xf

# custom function to estimate PSD based on the Bartlett method
# this is less accurate than the Welch method from matplotlib
# however, as the visual inspections confirm, the estimate gives
# the user a "reasonably good" view of the power spectrum
def estimatePSD(samples, NFFT, Fs):

	# rename the NFFT argument (notation consistent with matplotlib.psd)
	# to freq_bins (i.e., frequency bins for which we compute the spectrum)
	freq_bins = NFFT
	# frequency increment (or resolution of the frequency bins)
	df = Fs/freq_bins

	# create the frequency vector to be used on the X axis
	# for plotting the PSD on the Y axis (only positive freq)
	freq = np.arange(0, Fs/2, df)

	# design the Hann window used to smoothen the discrete data in order
	# to reduce the spectral leakage after the Fourier transform
	hann = np.empty(freq_bins)
	for i in range(len(hann)):
		hann[i] = pow(math.sin(i*math.pi/freq_bins),2)

	# create an empty list where the PSD for each segment is computed
	psd_list = []

	# samples should be a multiple of frequency bins, so
	# the number of segments used for estimation is an integer
	# note: for this to work you must provide an argument for the
	# number of frequency bins not greater than the number of samples!
	no_segments = int(math.floor(len(samples)/float(freq_bins)))

	# iterate through all the segments
	for k in range(no_segments):

		# apply the hann window (using pointwise multiplication)
		# before computing the Fourier transform on a segment
		windowed_samples = samples[k*freq_bins:(k+1)*freq_bins] * hann

		# compute the Fourier transform using the built-in FFT from numpy
		Xf = np.fft.fft(windowed_samples, freq_bins)

		# note, you can check how MUCH slower is DFT vs FFT by replacing the
		# above function call with the one that is commented below
		#
		# Xf = DFT(windowed_samples)
		#
		# note: the slow impelementation of the Fourier transform is not as
		# critical when computing a static power spectra when troubleshooting
		#
		# note also: time permitting a custom FFT can be implemented

		# since input is real, we keep only the positive half of the spectrum
		# however, we will also add the signal energy of negative frequencies
		# to have a better a more accurate PSD estimate when plotting
		Xf = Xf[0:int(freq_bins/2)] # keep only positive freq bins
		psd_seg = (1/(Fs*freq_bins/2)) * (abs(Xf)**2) # compute signal power
		psd_seg = 2*psd_seg # add the energy from the negative freq bins

		# translate to the decibel (dB) scale
		for i in range(len(psd_seg)):
			psd_seg[i] = 10*math.log10(psd_seg[i])

		# append to the list where PSD for each segment is stored
		# in sequential order (first segment, followed by the second one, ...)
		psd_list.extend(psd_seg)

	# compute the estimate to be returned by the function through averaging
	psd_est = np.zeros(int(freq_bins/2))

	# iterate through all the frequency bins (positive freq only)
	# from all segments and average them (one bin at a time ...)
	for k in range(int(freq_bins/2)):
		# iterate through all the segments
		for l in range(no_segments):
			psd_est[k] += psd_list[k + l*int(freq_bins/2)]
		# compute the estimate for each bin
		psd_est[k] = psd_est[k] / no_segments
	print("psd_est size = ",len(psd_est), end='\n')
	#print("psd_est = ",psd_est, end='\n')

	# the frequency vector and PSD estimate
	return freq, psd_est

# custom function to format the plotting of the PSD
def fmPlotPSD(ax, samples, Fs, height, title):

	x_major_interval = (Fs/12)		# adjust grid lines as needed
	x_minor_interval = (Fs/12)/4
	y_major_interval = 20
	x_epsilon = 1e-3
	x_max = x_epsilon + Fs/2		# adjust x/y range as needed
	x_min = 0
	y_max = 10
	y_min = y_max-100*height
	#ax.psd(samples, NFFT=512, Fs=Fs)
	#
	# below is the custom PSD estimate, which is based on the Bartlett method
	# it less accurate than the PSD from matplotlib, however it is sufficient
	# to help us visualize the power spectra on the acquired/filtered data
	#
	freq, my_psd = estimatePSD(samples, NFFT=512, Fs=Fs)
	ax.plot(freq, my_psd)
	#
	ax.set_xlim([x_min, x_max])
	ax.set_ylim([y_min, y_max])
	ax.set_xticks(np.arange(x_min, x_max, x_major_interval))
	ax.set_xticks(np.arange(x_min, x_max, x_minor_interval), minor=True)
	ax.set_yticks(np.arange(y_min, y_max, y_major_interval))
	ax.grid(which='major', alpha=0.75)
	ax.grid(which='minor', alpha=0.25)
	ax.set_xlabel('Frequency (kHz)')
	ax.set_ylabel('PSD (db/Hz)')
	ax.set_title(title)

if __name__ == "__main__":

	# do nothing when this module is launched on its own
	pass
