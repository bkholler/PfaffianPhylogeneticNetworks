#!/usr/bin/env python

import sys, getopt, errno, time
import numpy as np
import itertools
import random
import skbio
from skbio import TabularMSA, DNA

start = time.time()

msa_length = 1000
outputFilename = "msa.phylip"

nucleotides = ["A", "C", "G", "T"]
substitutionRates = list()
generageEdgeRates = True
generateGamma = True
numberLeaves = 6

try:
	opts, args = getopt.getopt(sys.argv[1:],"hs:g:l:o:n:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python simulateSunletSequenceK3P.py -l <MSA length> -s <seed> -g <gamma parameter> -n <number of leaves> -o <output filename>")
	print("python simulateSunletSequenceK3P.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("python simulateSunletSequenceK3P.py -l <MSA length> -s <seed> -g <gamma parameter> -o <output filename>")
		print("-l <MSA length>\t\t Length of MSA to output")
		print("-s <seed>\t\t Integer value for seeding random numbers")
		print("-g <gamma parameter>\t\t Floating point value giving the probability of a site evolving along edge e in the network.")
		print("-n <number of leaves>\t\t Positive integer giving the number of leaves in the sunlet.")
		print("-o <output filename>\t Filename for output MSA in phylip format")
		sys.exit()
	elif opt in ("-o"):
		outputFilename = arg 
	elif opt in ("-l"):
		msa_length = int(arg)
	elif opt in ("-s"):
		random.seed(int(arg))
	elif opt in ("-g"):
		val = float(arg)
		if val < 0 or val > 1:
			print ("Error: User supplied gamma parameter is not in interval [0,1]: " + param)
			sys.exit(2)
		gamma = val
		generateGamma = False
	elif opt in ("-n"):
		numberLeaves = int(arg)
		if numberLeaves < 3:
			print("Error: Number of leaves must be at least 3.")
			sys.exit(-1)


if generateGamma:
	gamma = random.random()


for i in range(2*numberLeaves):
	random1 = random.uniform(0.95,1)
	random2 = random.uniform(0,0.02)
	random3 = random.uniform(0,0.02)
	random4 = random.uniform(0,0.02)
	randomSum = random1 + random2 + random3 + random4
	substitutionRates.append( [random1/randomSum, random2/randomSum, random3/randomSum, random4/randomSum])

#print("Simulating alignments with parameters: ")
#edge_params = ""
#for edgeRates in substitutionRates:
#	edge_params += "(" + str(edgeRates[0]) + ","+ str(edgeRates[1]) + ","+ str(edgeRates[2]) + ","+ str(edgeRates[3]) + "), "
#print("Edge Parameters: " + edge_params[0:-1])
#print("Gamma: " + str(gamma))

def mutate_K3P(start_nucl, edge_param):
	rand = random.random()
	if start_nucl == "A":
		if rand < edge_param[0]:
			return "A"
		elif rand < edge_param[0] + edge_param[1]:
			return "G"
		elif rand < edge_param[0] + edge_param[1] + edge_param[2]:
			return "C"
		else:
			return "T"
	elif start_nucl == "G":
		if rand < edge_param[1]:
			return "A"
		elif rand < edge_param[1] + edge_param[0]:
			return "G"
		elif rand < edge_param[1] + edge_param[0] + edge_param[3]:
			return "C"
		else:
			return "T"
	elif start_nucl == "C":
		if rand < edge_param[2]:
			return "A"
		elif rand < edge_param[2] + edge_param[3]:
			return "G"
		elif rand < edge_param[2] + edge_param[3] + edge_param[0]:
			return "C"
		else:
			return "T"
	else: # start_nucl == "T":
		if rand < edge_param[3]:
			return "A"
		elif rand < edge_param[3] + edge_param[2]:
			return "G"
		elif rand < edge_param[3] + edge_param[2] + edge_param[1]:
			return "C"
		else:
			return "T"

sequences = [""] * numberLeaves
for i in range(msa_length):

	cycleVerts = [None] * numberLeaves
	# start with a uniform distribution at the root = cycleVerts[1]
	cycleVerts[1] = nucleotides[random.randrange(4)]
	for i in range(2,numberLeaves):
		cycleVerts[i] = mutate_K3P(cycleVerts[i-1], substitutionRates[numberLeaves + i-1])

	if random.random() < gamma:
		cycleVerts[0] = mutate_K3P(cycleVerts[1], substitutionRates[numberLeaves])
	else:
		cycleVerts[0] = mutate_K3P(cycleVerts[-1], substitutionRates[-1])

	for i in range(numberLeaves):
		sequences[i] += mutate_K3P(cycleVerts[i], substitutionRates[i])

skbioSeqs = list()
for i in range(numberLeaves):
	skbioSeqs.append(DNA(sequences[i], metadata={"id":"Taxon" + str(i)}))

aln = TabularMSA(skbioSeqs, minter="id")

f = open(outputFilename, "w")
aln.write(f, format='phylip')
f.close()

end = time.time()
print("simulateSunletSequenceK3P.py time: " + str(end - start))


