#!/usr/bin/python

import ctypes

cgplib = ctypes.CDLL('./libcgp.so')

inputs = 1
nodes = 20
outputs = 1
arity = 2

params = cgplib.initialiseParameters(inputs,nodes,outputs,arity)

cgplib.addNodeFunction(params, "add,sub,mul,div,sin")

cgplib.printParameters(params)

trainingData = cgplib.initialiseDataSetFromFile("symbolic.data")

chromo = cgplib.runCGP(params, trainingData, 100)

cgplib.printChromosome(chromo)


cgplib.getChromosomeFitness.restype = ctypes.c_float

f = cgplib.getChromosomeFitness(chromo)





print f
