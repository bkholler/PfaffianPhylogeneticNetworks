#!/usr/bin/env python

import sys, getopt, errno, time
import numpy as np
import itertools
import mpmath
import skbio

from multiprocessing import Pool
from skbio import DNA, TabularMSA
from itertools import permutations

start = time.time()

#Metric Functions:
def get1Norm(values):
	total = mpmath.mpf(0)
	for val in values:
		total = mpmath.fadd(total, mpmath.fabs(val))
	return total

def isCanonical(perm, n):
	return perm[1] < perm[n-1]

def getFrequenciesFromPermutation(frequencyArray, permutation, n):
	returnArray = np.zeros(shape=(2,)*n, dtype=np.longdouble)
	for index, freq in np.ndenumerate(frequencyArray):
		permutedIndex = [0] * n
		for i in range(n):
			permutedIndex[i] = index[permutation[i]]
		returnArray[tuple(permutedIndex)] = freq

	return returnArray

def evaluatePermutedAlignment(frequencies, perm, n):
	#find the permutation under which this leaf-labelled sunlet is invariant
	#print("Started calculation for perm: " + str(perm))
	invariantPermList = list()
	invariantPermList.append(perm[0])
	for i in range(1, n):
		invariantPermList.append(perm[n-i])
	invariantPerm = tuple(invariantPermList)
	assert(not isCanonical(invariantPerm, n))

	score = 0
	for p in [perm, invariantPerm]:

		#get permuted frequencies array
		permutedFrequencies = getFrequenciesFromPermutation(frequencies, p, n)
		# Construct M matrix
		M = np.zeros((n, n))
		for i in range(n):
			for j in range(n):
				if i > j:
					qIndex = [0]*n
					qIndex[i] = 1
					qIndex[j] = 1
					# Perform Fourier transformation to get q_{e_i + e_j}
					transformedValue = mpmath.mpf(0)
					for pIndex, freq in np.ndenumerate(permutedFrequencies):
						if freq > 0:
							summand = freq
							for l in range(n):
								summand *= Chi(qIndex[l], pIndex[l])
							transformedValue = mpmath.fadd(transformedValue, summand)
					val = mpmath.fdiv(transformedValue, mpmath.power(2,n))
					M[i,j] = val
					M[j,i] = -val

		results = list()
		# evaluate the 2x2 minors
		for i in range(1,n):
			for j in range(i+1, n):
				for k in range(j+1, n):
					for l in range(k+1, n):
						#print("{" + str(i) + "," + str(j) + "},{" + str(k) + "," + str(l) + "}")
						results.append(np.linalg.det(M[np.ix_([i,j],[k,l])]))

		# evaluate the 3x3 minors
		for i in range(1,n):
			for j in range(i+1, n):
				for k in range(j+1, n):
					for l in range(k+1, n):
						for m in range(l+1, n):
							#print("{ 1," + str(i) + "," + str(m) + "},{" + str(j) + "," + str(k) + "," + str(l) + "}")
							results.append(np.linalg.det(M[np.ix_([0,i,m],[j,k,l])]))
		score += get1Norm(results)
	#print("complete perm: " + str(perm))
	return score

PP = {
  "A": 0, 
  "G": 0, 
  "C": 1, 
  "T": 1 
}

# Function to return character values when calculating Fourier transform
def Chi(i, j):
	if i == 1 and j ==1: 
		return -1
	else:
		return 1

if __name__ == '__main__':
	MSAFilename = ""
	scoringFunction = get1Norm
	numProcesses = 16

	try:
		opts, args = getopt.getopt(sys.argv[1:],"ha:t:")
	except getopt.GetoptError:
		print("Option not recognised.")
		print("python evaluate_pfaffians.py -a <MSA file> -t <num threads>")
		print("python evaluate_pfaffians.py -h for further usage instructions.")
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-h":
			print("python evaluate.py -a <MSA file> -t <num threads>")
			print("-a <MSA file>\t\t Multiple sequence alignment file.")
			sys.exit()
		elif opt in ("-a"):
			MSAFilename = arg
		elif opt in ("-t"):
			numProcesses = int(arg)

	if len(MSAFilename) == 0:
		print("Error: You must provide an MSA file with -a.")
		sys.exit(2)

	try:
		alignment = TabularMSA.read(MSAFilename, constructor=DNA)
	except (ValueError, TypeError, skbio.io.UnrecognizedFormatError) as e:
		print(e)
		sys.exit(2)
	if not alignment or len(alignment) < 4:
		print("Error: MSA file must be a multiple sequence alignment of at least 4 sequences.")
		sys.exit(2)

	n = len(alignment)
	scores = dict()

	count = 0
	frequencies = np.zeros(shape=(2,)*n, dtype=np.longdouble)
	for col in alignment.iter_positions(ignore_metadata=True):
		if "-" not in col:
			index = [-1]*n
			for i in range(n):
				index[i] = PP[str(col[i])]
			frequencies[tuple(index)] += 1
			count += 1

	for index, freq in np.ndenumerate(frequencies):
		if freq != 0:
			frequencies[index] = mpmath.fdiv(freq, count)

	# get the permutations under which this sunlet is not invariant.
	results = dict()
	perms = permutations(range(n))
	pool = Pool(processes=numProcesses)
	for perm in perms:
		if not isCanonical(perm, n):
			continue
		results[perm] = pool.apply_async(evaluatePermutedAlignment, args=(frequencies, perm, n))
	pool.close()
	pool.join()

	for key in results.keys():
		scores[key] = results[key].get()

	rank = 1
	sortedScores = dict(sorted(scores.items(), key=lambda item: item[1]))
	for key in sortedScores.keys():
		score = scores[key]
		print(str(rank) + ":\t" + str(key) + ":\t" + str(score))
		rank += 1

	end = time.time()
	print("evaluate_pfaffians_v2.py time: " + str(end - start))


