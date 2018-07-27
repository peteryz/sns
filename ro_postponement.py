##############################################################################

## Code for postponement in application section of He and Yehua's IJOO paper

#  received from Yehua March 2018.

#  Row and column generation.

## Modification by Peter from March to April 2018:

# Take input from txt files (to sync input with Julia code)

# Run a slightly different formulation: (1) uncertainty set, (2) capacity (inv)

# as decision variable, (3) demand loss penalty \in integers from the start.

##############################################################################

## TODO ##
# Capacity var
# Capacity costs
# Flow costs (already there?)
# Groups of capacity and flows, indexed by k
# pp var replicate p, indexed by k


from gurobipy import *

from itertools import combinations

from math import sqrt

from math import floor

import time

import random

import numpy

import copy

from scipy.stats import norm

#import pylab as pl

import csv

import os

import sys


def ro_opt_rowcol(plants, products, allNodes, dmean, dvar, underage_margin, bits, invCost, flowCost, arcs, paths, instance):

    for i,j in paths:
      print i,j,flowCost[(i,j)], "\n"

    #trajectory_file = open('IO/log_instance_' + str(instance) + '.txt','w')
    #trajectory_file.write("m.objVal, ")
    #trajectory_file.write("second stage cost, ")
    #trajectory_file.write("bi.objVal, ")
    #trajectory_file.write("bi.MIPGap, ")
    #trajectory_file.write("Time limit, ")
    #trajectory_file.write("Cut iteration, ")
    #trajectory_file.write('\n')

    start = time.time()

    postponed_cap = 0.2

    # Create master optimization model

    m = Model('constraint_generate')

    m.setParam('OutputFlag', 0)

    m.setAttr("ModelSense", GRB.MINIMIZE)



    # Inventory variables, for plant and product inv

    s_plant = {}

    s_prod = {}

    # delta represents sum of first stage inv cost, 2nd stage underage cost, and flow costs

    delta = m.addVar(name='miss_demand', obj=1)

    for i in plants:

        s_plant[i] = m.addVar(lb=0.0, name='sn_%s' % i, obj=invCost[i])

        m.update()

    for j in products:

        s_prod[j] = m.addVar(lb=0.0, name='s_%s' % j, obj=invCost[j+len(plants)])

        m.update()

#    m.addConstr(quicksum(s[j] for j in products) == preallocated_cap*total_capacity)

    m.optimize()

    flow = {}

    lost = {}



    # second stage optimization (constraint generation)

    bi = Model('bilinear_sub')

    bi.setParam('OutputFlag', 0)

    p = {}; q = {}; y = {}; d = {}; xi = {};

    for j in products:

        for n in bits:

            y[(j,n)] = bi.addVar(obj=2**n)

            q[(j,n)] = bi.addVar(vtype=GRB.BINARY)

        #variables for demand and absolute demand deviation

#        w[j] = bi.addVar(name='w_%s' % j)

        d[j] = bi.addVar(name='d_%s' % j, ub=dmean[j]+dvar[j])

    for i in plants:

        p[i] = bi.addVar()

    for i in allNodes:

        xi[i] = bi.addVar(lb=0.0, ub=1.0)

    bi.update()



    for i,j in paths:

        bi.addConstr(quicksum(q[(j,n)]*2**n for n in bits) - p[i]

                    <= flowCost[(i,j)]) # flexibility constraint

    for j in products:

        for n in bits:

            bi.addConstr(y[(j,n)] <= (dmean[j]+dvar[j])*q[(j,n)])

            bi.addConstr(y[(j,n)] <= d[j])

            bi.addConstr(y[(j,n)] >= d[j] + (dmean[j]+dvar[j])*q[(j,n)] - (dmean[j]+dvar[j])) #redudant

        bi.addConstr(quicksum(q[(j,n)]*2**n for n in bits) <= underage_margin[j])

        #bi.addConstr(demand[j] - d[j] <= sigma[j]*w[j])

        #bi.addConstr(d[j] - demand[j] <= sigma[j]*w[j])



    #bi.addConstr(quicksum(w[j] for j in products) <= beta)



    # Demand uncertainty description

    # 1. xi[i] \in [0,1]. Done. Constructed as this.

    # 2. xi[0] == 1

    # 3. xi[i] >= xi[j] for all i,j in adj

    # 4. for i in allNodes, xi[i]*gamma[i] >= sum(for j in allNodes * adj[i,j],

    # 5. d[j] == dmean[j] + dvar[j] * xi[j+num_plants]

    bi.addConstr(xi[0] == 1)

    for i,j in arcs:

        bi.addConstr(xi[i] >= xi[j])

    for i in allNodes:

        bi.addConstr(xi[i] * gamma[i] >= quicksum(xi[j] for i,j in arcs.select(i,'*')))

    for j in products:

        bi.addConstr(d[j] == dmean[j] + dvar[j] * xi[j + len(plants)])

    time_limit = 60
    master_tolerance = 0.001

    is_opt = False

    count_constr = 0

    bi.setAttr("ModelSense", GRB.MAXIMIZE)

    bi.setParam("MIPGap", 0.05)

    bi.setParam("TimeLimit", time_limit)

    while is_opt == False:

        is_opt = True

        # Solve bilinear programming subproblem

        for j in products:

            for n in bits:

                q[(j,n)].setAttr("obj", -s_prod[j].x*2**n)

        for i in plants:

            p[i].setAttr("obj",-s_plant[i].x)

        #bi.setAttr("objCon", -delta.x)

        olddelta = delta.x

        bi.optimize()

        #if (bi.objVal > delta.x) or (bi.objVal <= delta.x and bi.objVal + bi.MIPGap * abs(bi.objVal) - delta.x > master_tolerance * m.objVal) or m.objVal == 0:
        if (bi.objVal > delta.x):

            #add columns (variables) to the master problem

            for j in products:

                lost[j, count_constr] = m.addVar()

            for i,j in paths:

                flow[i,j, count_constr] = m.addVar()

            m.update()

            #add flow constraints to the master problem

            for i in plants:

                m.addConstr(quicksum(flow[i,j,count_constr] for i,j in paths.select(i,'*')) <= s_plant[i])

            for j in products:

                m.addConstr(quicksum(flow[i,j, count_constr] for i,j in paths.select('*',j)) + lost[j, count_constr] + s_prod[j] >= d[j].x)

                # Lost sales constaint

            m.addConstr(quicksum(flow[i,j, count_constr] * flowCost[(i,j)] for i,j in paths) + quicksum(lost[j, count_constr]*underage_margin[j] for j in products)<= delta)

            count_constr += 1

            is_opt = False

            if bi.status == GRB.Status.TIME_LIMIT:
              time_limit = time_limit + 60
              bi.setParam("TimeLimit", time_limit)
            if bi.status == GRB.Status.OPTIMAL:
              time_limit = 60
              bi.setParam("TimeLimit", time_limit)
            #trajectory_file.write(str(m.objVal)+', ')
            #trajectory_file.write(str(delta.x)+', ')
            #trajectory_file.write(str(bi.objVal)+', ')
            #trajectory_file.write(str(bi.MIPGap)+', ')
            #trajectory_file.write(str(time_limit)+', ')
            #trajectory_file.write(str(count_constr))
            #trajectory_file.write('\n')

        print ("Constr: %3d" % count_constr)

        print ("MILP")

        print ("violate: %8.2f" % (bi.objVal))

        print ("Delta:")

        print (olddelta)

            # Compute optimal solution

        m.optimize()

        if m.status != GRB.status.OPTIMAL:

            print ("Iteration Error")

            break



    end = time.time()

#    print ("Time:", end - start)

#    print ("Flow and demand loss cost: ", delta.x)

#    print ("Inventory cost: ", sum(s_plant[i].x*invCost[i] for i in plants) + sum(s_prod[j].x*invCost[j+len(plants)] for j in products))

#    print ("Total cost: ", m.objVal)


    totalInvCost = sum(s_plant[i].x*invCost[i] for i in plants) + sum(s_prod[j].x*invCost[j+len(plants)] for j in products)

    totalDemandLossCost = 0.0
    for k in range(count_constr):
      totalDemandLossCost = max(sum(lost[j,k].x * underage_margin[j] for j in products), totalDemandLossCost)

    totalFlowCost = m.objVal - totalInvCost - totalDemandLossCost


    sol = [s_plant[i].x for i in plants] + [s_prod[j].x for j in products]

    #print "Total inventory", sum(sol)

    #trajectory_file.close()

    return end-start, count_constr, sol, m.objVal, totalInvCost, totalFlowCost, totalDemandLossCost, bi.objVal, bi.MIPGap



# Start Script



num_instances = 1 ## Specify the total number of inputs here

instances = xrange(num_instances)

numNodes = range(num_instances)

write_data = open('IO/colrow_output_instance_' + sys.argv[1] + '.txt','w')

constrs_diffinstance = range(num_instances)

time_diffinstance = range(num_instances)

objVal = range(num_instances)



constrs_diffinstance_tree = range(num_instances)

time_diffinstance_tree = range(num_instances)

objVal_tree = range(num_instances)


instance = int(sys.argv[1])

print ("instance: ", instance)

### Read in the following parameters

# 0 - numNodes, ",",

# 1 - numTreeArcs, ",",

# 2 - numFullArcs, ",",

# 3 - numLevels, ",",

# 4 - length(find(supplyNodes)), ",",

# 5 - length(find(demandNodes)) ),

sizeFile  = open("./IO/instance_"+str(instance)+"_size.txt", "r")

sizeInput = sizeFile.read().strip().split(',')




## new input by me

num_plants = int(sizeInput[4])

plants = xrange(num_plants)

num_products = int(sizeInput[5])

products = xrange(num_products)

allNodes = xrange(num_plants+num_products)

numNodes[0] = num_plants+num_products

#samples = xrange(1000) # not used anywhere



## new input by me

# In my formulation, inventories (He's capacities) are decision var.

# I will remove the use of this parameter in my modification, but set dummy

# for now

#capacity = {}

#for i in range(num_plants):

#    capacity[i] = 1000000

#total_capacity = sum(capacity.values())



## new input by me

# demand.txt format:

# 0 - node index i (for all nodes)

# 1 - mean demand \bar{d}_i

# 2 - surge demand \hat{d}_i

# 3 - demand loss penalty at i

# 4 - budget for this node \Gamma_i



dmean = {} # final length = num_products

dvar = {} # final length = num_products

gamma = {} # final length = num_plants + num_products

underage_margin = {} # final length = num_products



demandFile = open("IO/instance_"+str(instance)+"_demand.txt", "r")

for line_in in demandFile:

    line = line_in.strip().split(',')

    node_id = int(line[0]) # node_id starts from 1

    gamma[node_id-1] = int(line[4])

    if node_id > num_plants:

        dmean[node_id - 1 - num_plants] = float(line[1])

        dvar[node_id - 1 - num_plants] = float(line[2])

        ## Assumption: input is integer, so He's code runs to full accuracy

        underage_margin[node_id - 1 - num_plants] = float(line[3])



## new input by me

# underage_margin already populated when reading demand input file.

bits = xrange(16) # how many digits do we need? for dual var q in subproblem,

## Additional input by me

# Need to add: inventory holding costs for all nodes, flow costs

invCost = {}

invCostFile = open("IO/instance_"+str(instance)+"_invCost.txt", "r")

for line_in in invCostFile:

    line = line_in.strip().split(',')

    node_id = int(line[0]) # node_id starts from 1

    invCost[node_id - 1] = float(line[1])





paths = [] ## Note: path (i,j), i = plant index, j = product index, following He's legacy code

flowCost = {} ##

flowCostFile = open("IO/instance_"+str(instance)+"_pathFull.txt", "r")

for line_in in flowCostFile:

    line = line_in.strip().split(',')

    i = int(line[0]) - 1

    j = int(line[1]) - 1

    if i < num_plants and j >= num_plants: ## Note: only include paths between supply and demand nodes

        paths += [(i,j-num_plants)]

        flowCost[(i,j-num_plants)] = float(line[2])

paths = tuplelist(paths)



paths_tree = [] ## Note: path (i,j), i = plant index, j = product index, following He's legacy code

flowCost_tree = {} ##

flowCostFile_tree = open("IO/instance_"+str(instance)+"_pathTree.txt", "r")

for line_in in flowCostFile_tree:

    line = line_in.strip().split(',')

    i = int(line[0]) - 1

    j = int(line[1]) - 1

    if i < num_plants and j >= num_plants: ## Note: only include paths between supply and demand nodes

        paths_tree += [(i,j-num_plants)]

        flowCost_tree[(i,j-num_plants)] = float(line[2])

paths_tree = tuplelist(paths_tree)



## arcs just tree arcs, base tree network

arcs = [] ## Note: arc (i,j): both i and j are unique *node id*. Need to postprocess in formulation. It is only used to construct demand uncertainty set.

arcFile = open("IO/instance_"+str(instance)+"_adjTree.txt", "r")

for line_in in arcFile:

    line = line_in.strip().split(',')

    i = int(line[0]) - 1

    j = int(line[1]) - 1

    arcs += [(i,j)]

arcs = tuplelist(arcs)


## new optimization call by me

#    [time_diffinstance[k], constrs_diffinstance[k], sol, objVal[k]] = ro_opt_rowcol(plants, products, allNodes, dmean, dvar, underage_margin, bits, invCost, flowCost, arcs, paths)


[time_diffinstance_tree[0], constrs_diffinstance_tree[0], sol_tree, objVal_tree[0], totalInvCost_tree, totalFlowCost_tree, totalDemandLossCost_tree, subObj_tree, subGap_tree] = ro_opt_rowcol(plants, products, allNodes, dmean, dvar, underage_margin, bits, invCost, flowCost_tree, arcs, paths_tree, instance)

[time_diffinstance[0], constrs_diffinstance[0], sol, objVal[0], totalInvCost, totalFlowCost, totalDemandLossCost, subObj, subGap] =                                              ro_opt_rowcol(plants, products, allNodes, dmean, dvar, underage_margin, bits, invCost, flowCost,      arcs, paths, instance)

write_data.write('instance#, ')

write_data.write('Num nodes, ')

write_data.write('Obj tree, ')

write_data.write('time tree, ')

write_data.write('#constraints tree, ')

write_data.write("inv cost tree, ")

write_data.write("flow cost tree, ")

write_data.write("demand loss cost tree, ")

write_data.write("Subproblem objective tree, ")

write_data.write("Subproblem gap tree, ")

write_data.write('Obj general, ')

write_data.write('time general, ')

write_data.write('#constraints general, ')

write_data.write("inv cost general, ")

write_data.write("flow cost general, ")

write_data.write("demand loss cost general, ")

write_data.write("Subproblem objective general, ")

write_data.write("Subproblem gap general")

write_data.write('\n')

write_data.write(str(instance)+',')

write_data.write(str(numNodes[0])+',')

write_data.write(str(objVal_tree[0])+',')

write_data.write(str(time_diffinstance_tree[0])+',')

write_data.write(str(constrs_diffinstance_tree[0])+',')

write_data.write(str(totalInvCost_tree)+',')

write_data.write(str(totalFlowCost_tree)+',')

write_data.write(str(totalDemandLossCost_tree)+',')

write_data.write(str(subObj_tree)+',')

write_data.write(str(subGap_tree)+',')

write_data.write(str(objVal[0])+',')

write_data.write(str(time_diffinstance[0])+',')

write_data.write(str(constrs_diffinstance[0])+',')

write_data.write(str(totalInvCost)+',')

write_data.write(str(totalFlowCost)+',')

write_data.write(str(totalDemandLossCost)+',')

write_data.write(str(subObj)+',')

write_data.write(str(subGap))

write_data.write('\n')

write_data.close()
