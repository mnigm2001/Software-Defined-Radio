
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD, impResponse, lpFilter

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10
audio_Fs = 48e3

audio_Fc = 	16e3			# change as needed (see spec in lab document)
audio_decim = 5			# change as needed (see spec in lab document)
audio_taps = 151			# change as you see fit

il_vs_th = 0

def diff_decod(manch_data):
	diff_data = np.empty(len(manch_data))
	diff_data[0] = manch_data[0]
	for i in range(1, len(manch_data)):
		if(manch_data[i] != manch_data[i-1]):
			diff_data[i] = 1
		elif(manch_data[i] == manch_data[i-1]):
			diff_data[i] = 0
	return diff_data

#mult is AND, addition is XOR
def matrixMult(decoded_data, parity_matrix):
	output = np.empty(10)
	temp_mat = np.empty(len(parity_matrix))
	for k in range(0, len(output)):
		ones_count = 0
		for i in range(0, len(temp_mat)):
			print("i = ", i, " k = ", k, end='\t')
			temp_mat[i] = decoded_data[i]*parity_matrix[i][k]
			print(" temp_mat[i] = ", temp_mat[i], " decoded_data[i] = ", decoded_data[i], " parity_matrix[i][k] =", \
				parity_matrix[i][k])
			if(temp_mat[i] == 1):
				ones_count += 1
		print("ones_count = ", ones_count)
		if(ones_count%2 == 1):
			output[k] = 1
		else:
			output[k] = 0

	return output


if __name__ == "__main__":

	#manch_data = np.array([0, 1, 0, 1, 0 ,1, 0, 0, 1, 1])
	#diff_data = diff_decod(manch_data)
	#print("manch_data = ", manch_data)
	#print("diff_data = ", diff_data)

	h = np.array([[1,0,0,0,0,0,0,0,0,0],
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
	[1,1,0,0,0,1,1,0,1,1]])

	data = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
	data2 = np.zeros(26)
	data3 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0])
	dataC = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,0,0])
	datatest = np.array([1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1])
	 [0. 1. 1. 0. 1. 1. 1. 1. 0. 1.]

	print("h = ", h)
	print(len(data3))
	output = matrixMult(datatest, h)
	print("output = ", output)
