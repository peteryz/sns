#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 20:52:09 2018

@author: pyzhang
"""
import numpy

# Beware:
# When reading node source, sink, the indexing in python needs to be i = source-1


instance = 1

# some groups might have empty set of arcs, filter them out.

numNodes = range(instance)

numNodes[0] = 104

cTree = {}
cCostTree = {}
cArcsTree = {}
capacityTreeFile = open("IO/instance_"+str(instance)+"_capacityTree.txt", "r")
index = 0
for line_in in capacityTreeFile:

    line = line_in.strip().split(',')    
    arcString = line[2].strip().split(';')    
    if len(arcString) == 1: # no arcs
        continue

    cTree[index] = float(line[0])    
    cCostTree[index] = float(line[1])
    cArcsTree[index] = () # tuple is hashable
    
    for i in range(len(arcString)):
        pairs = arcString[i].strip().split('-')
        if len(pairs) < 2:
            continue
        source = int(pairs[0]) - 1
        sink = int(pairs[1]) - 1
        cArcsTree[index] = cArcsTree[index] + ([source, sink],)
    index += 1

numGroupsTree = index
pathGroupTree = numpy.zeros((numNodes[0],numNodes[0],numGroupsTree))
for g in range(len(cArcsTree)):
    for e in range(len(cArcsTree[g])):
        print cArcsTree[g], cArcsTree[g][e], cArcsTree[g][e][0], cArcsTree[g][e][1], g
        i = cArcsTree[g][e][0]
        j = cArcsTree[g][e][1]
        pathGroupTree[i,j,g] = 1


cFull = {}
cCostFull = {}
cArcsFull = {}
capacityFullFile = open("IO/instance_"+str(instance)+"_capacityFull.txt", "r")
index = 0
for line_in in capacityFullFile:

    line = line_in.strip().split(',')    
    arcString = line[2].strip().split(';')        
    if len(arcString) == 1: # no arcs
        continue

    cFull[index] = float(line[0])    
    cCostFull[index] = float(line[1])
    cArcsFull[index] = () # tuple is hashable
    
    for i in range(len(arcString)):
        pairs = arcString[i].strip().split('-')
        if len(pairs) < 2:
            continue
        source = int(pairs[0]) - 1
        sink = int(pairs[1]) - 1
        cArcsFull[index] = cArcsFull[index] + ([source, sink],)
    index += 1

numGroupsFull = index
pathGroupFull = numpy.zeros((numNodes[0],numNodes[0],numGroupsFull))
for g in range(len(cArcsFull)):
    for e in range(len(cArcsFull[g])):
        i = cArcsFull[g][e][0]
        j = cArcsFull[g][e][1]
        pathGroupFull[i,j,g] = 1


print numpy.nonzero(pathGroupTree)

print numpy.nonzero(pathGroupFull)
