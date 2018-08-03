using JuMP
using JuMPeR
using Combinatorics
using Gurobi

include("affineFunctionsEngaging.jl")

## Engaging Directory
engagingDir = "IO/"

  ############################################################################
  ######################### User specify 1: meta info ########################
  ############################################################################
srand(314);
tempDir = engagingDir;
cleanDir = engagingDir;
numInstances = 3600;
numNodesList = zeros(numInstances);
affineObj = zeros(numInstances,2); # first column for tree network, second for general network based on this tree
affineSolTime = zeros(numInstances);
enumObj = zeros(numInstances,2); # objective values from enumeration
enumSolTime = zeros(numInstances);
affineLocalObj = zeros(numInstances,2); # first column for tree network, second for general network based on this tree
affineLocalSolTime = zeros(numInstances);
affinePathObj = zeros(numInstances,2);
affinePathSolTime = zeros(numInstances);
affineXiTreeObj = zeros(numInstances,2);
affineXiTreeSolTime = zeros(numInstances);
affineXiObj = zeros(numInstances,2);
affineXiSolTime = zeros(numInstances);

for instance = 1:numInstances

  println(string("instance: ", instance))
  dir = string(tempDir,"instance_",instance,"_");

  # Read size file
  sizes = zeros(6)
  sizeFile = open(string(dir,"size.txt"))
  lines = readlines(sizeFile)
  sizes = [parse(Int, s) for s in split(strip(lines[1]),',')] # for local machine

  numNodes = convert(Int16, sizes[1])
  numNodesList[instance] = numNodes
  numTreeArcs = convert(Int16, sizes[2])
  numFullArcs = convert(Int16, sizes[3])
  numLevels = convert(Int16, sizes[4])
  numSupplyNodes = convert(Int16, sizes[5])
  numDemandNodes = convert(Int16, sizes[6])

  # Based on tree: parent, kids, offspr, offsprNr,
  # parent: Dict(), kid -> parent
  # kids: Dict(), parent -> [kids]
  # offspr: Dict(), root -> [offsprings]
  # offsprNr: contains root too
  parent = Dict()
  kids = Dict()
  offspr = Dict()
  offsprNr = Dict()

  #Base on full network
  offsprFull = Dict()

  # Read adjacency file tree: i, j, arc cost
  Atree = spzeros(Int8, numNodes, numNodes)
  Ctree = spzeros(numNodes, numNodes)
  open(string(dir,"adjTree.txt")) do adjTreeFile
    lines = readlines(adjTreeFile)

    for l in lines
      #input = parse.(split(strip(l),',')) # can do this in 0.5 and beyond
      input = [parse(Float64, s) for s in split(strip(l),',')] # for local machine
      i = convert(Int16, input[1])
      j = convert(Int16, input[2])
      cost = input[3]
      Atree[i,j] = 1
      Ctree[i,j] = cost
      parent[j] = i
      kids[i] = push!(get(kids,i,[]),j)
    end
  end

  # Read adjacency file: i, j, arc cost
  Afull = spzeros(Int8, numNodes, numNodes)
  Cfull = spzeros(numNodes, numNodes)
  open(string(dir,"adjFull.txt")) do adjFullFile
    lines = readlines(adjFullFile)

    for l in lines
      #input = parse.(split(strip(l),',')) # can do this in 0.5 and beyond
      input = [parse(Float64, s) for s in split(strip(l),',')] # for local machine
      i = convert(Int16, input[1])
      j = convert(Int16, input[2])
      cost = input[3]
      Afull[i,j] = 1
      Cfull[i,j] = cost
    end
  end

  # pathTree and cost
  Ptree = spzeros(Int8, numNodes, numNodes)
  PCtree = spzeros(numNodes, numNodes)
  open(string(dir,"pathTree.txt")) do file
    lines = readlines(file)
    for l in lines
      #input = parse.(split(strip(l),',')) # can do this in 0.5 and beyond
      input = [parse(Float64, s) for s in split(strip(l),',')] # for local machine
      i = convert(Int16, input[1])
      j = convert(Int16, input[2])
      cost = input[3]
      Ptree[i,j] = 1
      PCtree[i,j] = cost
      offspr[i] = push!(get(offspr,i,[]),j)
      offsprNr[i] = push!(get(offsprNr,i,[]),j)
    end
  end
  for i=1:numNodes
    offsprNr[i] = push!(get(offsprNr,i,[]),i)
  end

  # pathFull files, i,j, path cost
  Pfull = spzeros(Int8, numNodes, numNodes)
  PCfull = spzeros(numNodes, numNodes)
  open(string(dir,"pathFull.txt")) do file
    lines = readlines(file)
    for l in lines
      #input = parse.(split(strip(l),',')) # can do this in 0.5 and beyond
      input = [parse(Float64, s) for s in split(strip(l),',')] # for local machine
      i = convert(Int16, input[1])
      j = convert(Int16, input[2])
      cost = input[3]
      Pfull[i,j] = 1
      PCfull[i,j] = cost
      offsprFull[i] = push!(get(offsprFull,i,[]),j)
    end
  end

  # inv file: node index, inv
  invCost = zeros(numNodes)
  open(string(dir,"invCost.txt")) do file
    lines = readlines(file)
    for l in lines
      #input = parse.(split(strip(l),',')) # can do this in 0.5 and beyond
      input = [parse(Float64, s) for s in split(strip(l),',')] # for local machine
      i = convert(Int16, input[1])
      cost = input[2]
      invCost[i] = cost
    end
  end


  # demand file: index, dmean, dvar, losspenalty, Gamma
  dmean = zeros(numNodes)
  dvar = zeros(numNodes)
  dLossPenalty = zeros(numNodes)
  Γ = zeros(numNodes)
  open(string(dir,"demand.txt")) do file
    lines = readlines(file)
    for l in lines
      #input = parse.(split(strip(l),',')) # can do this in 0.5 and beyond
      input = [parse(Float64, s) for s in split(strip(l),',')] # for local machine
      i = convert(Int16, input[1])
      dm = input[2]
      dv = input[3]
      b = input[4]
      γ = convert(Int16, input[5]) # imposing integer attack scales
      dmean[i] = dm
      dvar[i] = dv
      dLossPenalty[i] = b
      Γ[i] = γ
    end
  end

  # read in capacity parameters
    # Two capacity files, one for tree and one for nontree: capacityTree.txt and capacityFull.txt
    # Each line is a group, in the form of capacity, cost per capacity, i-j;a-b;
  # Store cost per capacity in an array, index is the group ID
  # Store arcs of capacity groups in a Dict(), index by group ID and

  println("")
  println("****")
  dir = "IO/instance_1_" # For local testing

  capCostTree = Float64[]
  capPathsTree = Dict()

  open(string(dir,"capacityTree.txt")) do file
    index = 1
    lines = readlines(file)
    for l in lines
      input = split(strip(l),',') # can do this in 0.5 and beyond
      if input[3] == ""
        continue
      end
      cap = parse(Float64, input[1])
      cost = parse(Float64, input[2])
      pathString = split(input[3],';')
      numPaths = length(pathString)-1
      #capCostTree = vcat(capCostTree, cost)
      push!(capCostTree, cost)
      paths = zeros(Int64, numPaths,2)
      for i=1:numPaths
        nodeString = split(pathString[i],'-')
        paths[i,1] = parse(Int64, nodeString[1])
        paths[i,2] = parse(Int64, nodeString[2])
      end
      capPathsTree[index] = paths
      index += 1
    end
  end

  numCapacityGroups = length(capCostTree)
  println("numCapacityGroups")
  println(numCapacityGroups)
  println("capCostTree:")
  println(capCostTree)
  println("capPathsTree:")
  println(capPathsTree)

  #### Generate helper data structures

  ###################### Solve affine formulation on this instance #######################

  ## Study 1: tree structure, with delay cost
  #F0return = spzeros(numNodes,numNodes);
  #Fdreturn = spzeros(numNodes,numNodes,numNodes);
  F0return = zeros();
  Fdreturn = zeros();

  affine_xreturn = zeros(numNodes);
  xreturn = zeros(numNodes);
  println("**********")
  #affineObj[instance,1] = affine(numNodes, Γ, Atree, Ptree, invCost, dLossPenalty, PCtree, dmean, dvar, F0return, Fdreturn, affine_xreturn, affineSolTime, instance);
  #affineObj[instance,2] = affine(numNodes, Γ, Atree, Pfull, invCost, dLossPenalty, PCfull, dmean, dvar, F0return, Fdreturn, affine_xreturn, affineSolTime, instance);
  #affineLocalObj[instance,1] = affineLocal(numNodes, Γ, Atree, Ptree, invCost, dLossPenalty, PCtree, dmean, dvar, F0return, Fdreturn, affine_xreturn, affineSolTime, instance);
  #affineLocalObj[instance,2] = affineLocal(numNodes, Γ, Atree, Pfull, invCost, dLossPenalty, PCfull, dmean, dvar, F0return, Fdreturn, affine_xreturn, affineSolTime, instance);
  #affinePathObj[instance,1] = affinePathLocalFS(numNodes, numSupplyNodes, Γ, Atree, Ptree, invCost, dLossPenalty, PCtree, dmean, dvar, F0return, Fdreturn, affine_xreturn, affinePathSolTime, instance);
  #affinePathObj[instance,2] = affinePathLocalFS(numNodes, numSupplyNodes, Γ, Atree, Pfull, invCost, dLossPenalty, PCfull, dmean, dvar, F0return, Fdreturn, affine_xreturn, affinePathSolTime, instance);

  affineXiTreeObj[instance,1] = affineXiTree(numNodes, numSupplyNodes, Γ, Atree, Ptree, Ptree, parent, kids,             offspr, offsprNr, invCost, dLossPenalty, PCtree, dmean, dvar, F0return, Fdreturn, affine_xreturn, affineXiTreeSolTime, instance);
  affineXiObj[instance,1] =         affineXi(numNodes, numSupplyNodes, Γ, Atree, Pfull, Ptree, parent, kids, offsprFull, offspr, offsprNr, invCost, dLossPenalty, PCfull, dmean, dvar, F0return, Fdreturn, affine_xreturn, affineXiSolTime,     instance);

  # Testing: passing tree into affineXi
  #affineXiObj[instance,1] =         affineXi(numNodes, numSupplyNodes, Γ, Atree, Ptree, Ptree, parent, kids, offspr, offspr, offsprNr, invCost, dLossPenalty, PCtree, dmean, dvar, F0return, Fdreturn, affine_xreturn, affineXiSolTime,     instance);

  #enumObj[instance,1] = fullyAdapt(numNodes, Γ,        Ptree, invCost, dLossPenalty, PCtree, dmean, dvar, xreturn, enumSolTime);
  #enumObj[instance,2] = fullyAdapt(numNodes, Γ,        Pfull, invCost, dLossPenalty, PCfull, dmean, dvar, xreturn, enumSolTime);

end


 ############# Output affine file ##################
affineOutputFile = open(string(cleanDir,"affineOutput.txt"), "w");
write(affineOutputFile,"instance#, #nodes, obj tree, time tree, obj, time", "\n");
for i=1:numInstances
  write(affineOutputFile, string(i, ",", numNodesList[i], ",", affineXiTreeObj[i,1], ",", affineXiTreeSolTime[i], ",", affineXiObj[i,1], ",", affineXiSolTime[i,1]));
  write(affineOutputFile, "\n");
end
close(affineOutputFile);

## Close all files
