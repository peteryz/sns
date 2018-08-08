# April 7. Generate data, run affine, output .txt file for He's python code.
# Combining generateData20180406.jl and study1-20180402.jl

using JuMP
using JuMPeR
using Combinatorics
using Gurobi

include("affineFunctionsEngaging.jl")

engagingDir = "IO/"

  ############################################################################
  ######################### User specify 1: meta info ########################
  ############################################################################
srand(314);
tempDir = engagingDir;
cleanDir = engagingDir;
numInstances = 200;
numNodesList = zeros(numInstances);
affineObj = zeros(numInstances,2); # first column for tree network, second for general network based on this tree
affineSolTime = zeros(numInstances);
enumObj = zeros(numInstances,2); # objective values from enumeration
enumSolTime = zeros(numInstances);
affineLocalObj = zeros(numInstances,2); # first column for tree network, second for general network based on this tree
affineLocalSolTime = zeros(numInstances);

for instance=1:numInstances

  println(string("###### Instance: ", instance, " ######"))
  ##################################################################################
  ########################### User specify 2: instance data ########################
  ##################################################################################
  println(string("instance: ", instance, ". Line 34."))

  ### April 17 new network generaiton method: generate by n ###
  ndiff = 7
  targets = zeros(ndiff)
  targets[1] = 100 # 100
  targets[2] = 200 # 200
  targets[3] = 500 # 500
  targets[4] = 1000 # 1000
  targets[5] = 2000 # 2000
  targets[6] = 5000 # 5000
  targets[7] = 10000 # 10000

  case = mod(instance, ndiff) == 0 ? ndiff : mod(instance, ndiff)
  numNodesTarget = targets[case]

  numLevels = rand(2:5); # randomize 2-5
  intensity = rand()*0.2; # 0 to 1, fraction of kids being impacted
  numChildrenPerLevel = zeros(numLevels);
  numChildrenPerLevel[1] = 1; # numChildrenPerLevel[i] is the number of children per node in level i-1
  numChildrenPerLevel[2] = rand(5:15);
  #numChildrenPerLevel[3] = 3 + rand(1:3);
  if numLevels >= 3
    numChildrenPerLevel[3] = rand(5:15);
    if numLevels >= 4
      numChildrenPerLevel[4:numLevels] = rand(5:15);
    end
  end
  # implement this function
  (numNodes, Atree, nodeLevel, pc, supplyNodes, demandNodes) = randomTreeByN(numNodesTarget, numChildrenPerLevel);
  numTreeArcs = sum(Atree);

  numNodesList[instance] = numNodes

  b = 100 # demand loss penalty
  h_low = (b-1)*rand()*0.1
  #h_high = h_low + (b-h_low)*rand()
  h_high = (b-1)*rand()*0.1 + (b-1)*0.1
  #f_root = (b-h_low)*rand()
  f_root = (b-1)*rand()*1.0
  survivability = rand()

  pArc = rand()*0.01

  #pcFull = pc;
  #Afull = copy(Atree);
  #numFullArcs = sum(Afull);


  # For generalized constraints: old, before Aug 7
  #numCapacityGroups = rand(1:100) # Change later for test
  #capProb = rand(numCapacityGroups) * 100.0 / numNodes # determines how many arcs per capacity group
  #capPerGroup = rand(numCapacityGroups) * 1000 * 10 / 2 # demand is 1000 at each node # determines the capacity for this group
  #capCostPerGroup = rand(numCapacityGroups) * b

  # Capacity group generation: new, after Aug 7
  # Parameter for instance is α ∈ [0,1], fraction of paths in each group
  # For each demand node d, get the number of paths ending at d: g_d
  # Randomly partition the g_d paths into min(g_d, Int16(1/α)) groups
  # Affine policy tractability:
  # In the special strucutre,all paths in one group go to the same demand node,
  # therefore we can actually have very efficient dualization with local policy
  numCapacityGroups = 0
  capCostPerGroup = rand() * b
  α = rand()


  println(string("instance: ", instance, ". Line 82."))

  pcFull = Dict();
  for i in keys(pc)
    pcFull[i] = collect(copy(pc[i]));
  end

  Afull = extraArcs(pcFull, pArc, numNodes);
  numFullArcs = sum(Afull);

  println(string("instance: ", instance, ". Line 92."))

  invCost = zeros(numNodes);
  for i=1:numNodes
    invCost[i] = h_low + (h_high - h_low) * (0.5 + 0.5*rand() + nodeLevel[i] - 1) / numLevels
  end

  # arc cost
  Ctree = spzeros(numNodes, numNodes);
  C = spzeros(numNodes, numNodes);
  for i in keys(pcFull)
    for j in pcFull[i]
      cost = (0.5 + 0.5 * rand()) * f_root / numLevels; # root to demand path cost is 0.5 to 1.0 f_root
      C[i,j] = cost;
      if get(pc,i,-1) != -1 && find(pc[i] .== j) != []
        Ctree[i,j] = cost;
      end
    end
  end
  dmean = zeros(numNodes);
  dvar = zeros(numNodes);
  for i in find(demandNodes)
    dvar[i] = 1000*rand();
  end
  dLossPenalty = b * ones(numNodes)
  Γ = zeros(Int16, numNodes);
  for i in keys(pc)
    Γ[i] = round(Int16, size(pc[i],1) * min(1, intensity * (0.8+0.4*rand())))
  end

  println(string("instance: ", instance, ". Line 122."))

#=
  ### General network generation ###
  # Format: supply nodes first, demand nodes follow. Root node first.
  numLevels = 3; # 3 to 6
  numChildren = 3; # max number of children per node.
  pArc = 0.1; # percentage addition on arcs. extra 0.1*numNodes^2 arcs
  intensity = 0.5; # to calibrate Γ[i]. higher demand intensity, more demand

  # Atree: adjacency matrix for tree base
  # nodeLeve: dict, node id => level of node. leaf is 1, root is numLevels
  (numNodes, Atree, nodeLevel, pc, supplyNodes, demandNodes) = randomTree(numLevels, numChildren);
  numTreeArcs = sum(Atree);

  println(string("Number of nodes: ",numNodes))

  numNodesList[instance] = numNodes;

  pcFull = Dict();
  for i in keys(pc)
    pcFull[i] = copy(pc[i]);
  end

  Afull = extraArcs(pcFull, pArc, numNodes);
  numFullArcs = sum(Afull);

  invCost = zeros(numNodes);
  for i=1:numNodes
    invCost[i] = 5 * (1 + numLevels - nodeLevel[i]) - 2*rand() # from root: 3 to 5, 8 to 10, 13 to 15, etc.
  end

  # arc cost
  Ctree = spzeros(Int16, numNodes, numNodes);
  C = spzeros(Int16, numNodes, numNodes);
  for i in keys(pcFull)
    for j in pcFull[i]
      cost = round(Int16, (10 + 4 * rand())/numLevels); # root to demand path has roughly cost of 10
      C[i,j] = cost;
      if get(pc,i,-1) != -1 && find(pc[i] .== j) != []
        Ctree[i,j] = cost;
      end
    end
  end

  dmean = zeros(numNodes);
  dvar = zeros(numNodes);
  for i in find(demandNodes)
    dvar[i] = 1000 * rand();
  end

  dLossPenalty = 100 * ones(numNodes)

  Γ = zeros(Int16, numNodes);
  for i in keys(pc)
    Γ[i] = max(1, round(Int8, (0.5 + rand() * intensity) * size(pc[i],1)));
  end

  ### End: General network generation ###
=#

  ##################################################################################
  ############################### Auto generate area ###############################
  ##################################################################################
  ### Generate path matrix from the adjanency matrix
  Ptree = sparse(zeros(Int8, numNodes, numNodes))
  Ptree = copy(Atree); # Each entry = number of numLevels-step paths
  PCtree = copy(Ctree); # initialize the path costs on tree as arc costs
  for i=1:(numLevels-1)
    makeOne(Ptree);
    updateCost(Ptree, Atree, PCtree, Ctree)
    Ptree = Ptree + Ptree * Atree;
  end
  makeOne(Ptree)
  Pfull = sparse(zeros(Int8, numNodes, numNodes))
  Pfull = copy(Afull); # Each entry = number of numLevels-step paths
  PC = copy(C) # initialize path costs as arc costs
  for i=1:(numLevels-1)
    makeOne(Pfull)
    updateCost(Pfull, Afull, PC, C)
    Pfull = Pfull + Pfull * Afull;
  end
  makeOne(Pfull);
  ### Construct path cost from arc costs
  # !!!!! Assumption 1: There is no flow capacity,
  # therefore the path from i to j is virtually unique, i.e., the cheapest one
  # (Otherwise my formulation as-is would work since the current implementation is path-based)

  println(string("instance: ", instance, ". Line 210."))

  ###################### Output instance data #######################

  dir = string(tempDir,"instance_",instance,"_");

  # Create files for input
  sizeFile = open(string(dir,"size.txt"), "w");
  numSupplyNodes = length(find(supplyNodes))
  numDemandNodes = length(find(demandNodes))
  write(sizeFile,
        string(numNodes, ",",
               numTreeArcs, ",",
               numFullArcs, ",",
               numLevels, ",",
               numSupplyNodes, ",",
               numDemandNodes, ",",
               numCapacityGroups),
        "\n");
  close(sizeFile);
  # i,j,flow cost for arc (i,j)
  adjTreeFile = open(string(dir,"adjTree.txt"), "w");
  nz = findn(Atree)
  for i=1:length(nz[1])
    write(adjTreeFile, string(nz[1][i], ",", nz[2][i], ",", C[nz[1][i],nz[2][i]]),"\n")
  end
  close(adjTreeFile);

  adjFullFile = open(string(dir,"adjFull.txt"), "w");
  nz = findn(Afull)
  for i=1:length(nz[1])
    write(adjFullFile, string(nz[1][i], ",", nz[2][i], ",", C[nz[1][i],nz[2][i]]),"\n")
  end
  close(adjFullFile);

  pathTreeFile = open(string(dir,"pathTree.txt"), "w");
  nz = findn(Ptree)
  for i=1:length(nz[1])
    write(pathTreeFile, string(nz[1][i], ",", nz[2][i], ",", PCtree[nz[1][i],nz[2][i]]),"\n")
  end
  close(pathTreeFile);

  pathFullFile = open(string(dir,"pathFull.txt"), "w");
  nz = findn(Pfull)
  for i=1:length(nz[1])
    write(pathFullFile, string(nz[1][i], ",", nz[2][i], ",", PC[nz[1][i],nz[2][i]]),"\n")
  end
  close(pathFullFile);

  #
  # capacity groups
  capacityTreeFile = open(string(dir,"capacityTree.txt"), "w");
  capacityTree = Dict()
  capacityFullFile = open(string(dir,"capacityFull.txt"), "w");
  capacityFull = Dict()
  ind = 0
  for j=numSupplyNodes+1:numNodes
    numJGroups = floor(α*sum(Ptree[:,j]))+1
    for k=1:numJGroups
      # capacity per group, capacity cost per group
      capThisGroup = 1000
      capCostThisGroup = rand() * b
      capacityTree[ind+k] = string(capThisGroup, ",", capCostThisGroup, ",")
      capacityFull[ind+k] = string(capThisGroup, ",", capCostThisGroup, ",")
    end
    for i=1:sum(Pfull[:,j])
      ancestor = findn(Pfull[:,j])[i]
      toGroup = rand(1:numJGroups)
      if Ptree[ancestor,j] == 1
        capacityTree[ind+toGroup] = string(capacityTree[ind+toGroup], ancestor, "-", j, ";")
      end
      capacityFull[ind+toGroup] = string(capacityFull[ind+toGroup], ancestor, "-", j, ";")
    end
    ind += numJGroups
    numCapacityGroups += numJGroups
  end
  for j=1:numCapacityGroups
    write(capacityTreeFile, string(capacityTree[j],"\n"))
    write(capacityFullFile, string(capacityFull[j],"\n"))
  end
  close(capacityTreeFile)
  close(capacityFullFile)


  # inv cost: node index (for all supply nodes), cost
  invCostFile = open(string(dir,"invCost.txt"), "w");
  for i=1:length(invCost)
    write(invCostFile, string(i,",",invCost[i]),"\n")
  end
  close(invCostFile)

  # Demand: node index, dmean, dvar, lossPenalty, Γ
  demandFile = open(string(dir,"demand.txt"), "w");
  for i=1:length(dmean)
    write(demandFile, string(i,",",dmean[i],",",dvar[i],",",dLossPenalty[i],",",Γ[i],",",survivability),"\n")
  end
  close(demandFile)

  ###################### Solve affine formulation on this instance #######################

  ## Study 1: tree structure, with delay cost
  #F0return = spzeros(numNodes,numNodes);
  #Fdreturn = spzeros(numNodes,numNodes,numNodes);
  println(string("instance: ", instance, ". Line 275."))

  affine_xreturn = zeros(numNodes);
  xreturn = zeros(numNodes);
  println("**********")
  #affineObj[instance,1] = affine(numNodes, Γ, Atree, Ptree, invCost, dLossPenalty, PCtree, dmean, dvar, F0return, Fdreturn, affine_xreturn, affineSolTime, instance);
  #affineObj[instance,2] = affine(numNodes, Γ, Atree, Pfull, invCost, dLossPenalty, PC,     dmean, dvar, F0return, Fdreturn, affine_xreturn, affineSolTime, instance);
  #affineLocalObj[instance,1] = affineLocal(numNodes, Γ, Atree, Ptree, invCost, dLossPenalty, PCtree, dmean, dvar, F0return, Fdreturn, affine_xreturn, affineSolTime, instance);
  #affineLocalObj[instance,2] = affineLocal(numNodes, Γ, Atree, Pfull, invCost, dLossPenalty, PC,     dmean, dvar, F0return, Fdreturn, affine_xreturn, affineSolTime, instance);

  #enumObj[instance,1] = fullyAdapt(numNodes, Γ,        Ptree, invCost, dLossPenalty, PCtree, dmean, dvar, xreturn, enumSolTime);
  #enumObj[instance,2] = fullyAdapt(numNodes, Γ,        Pfull, invCost, dLossPenalty, PC, dmean, dvar, xreturn, enumSolTime);

  println(string("instance, numNodes, numLevels, intensity, h_low, h_high, f_root, survivability, pArc, numFullArcs, numTreeArcs, numCapacityGroups, alpha, ",
                 instance, ", ", numNodes, ", ", numLevels, ", ", intensity, ", ", h_low, ", ",
                 h_high, ", ", f_root, ", ", survivability, ", ", pArc, ", ", numTreeArcs, ", ",
                 numFullArcs, ", ", numCapacityGroups, ", ", α))

end

 ############# Output affine file ##################
affineOutputFile = open(string(cleanDir,"affineOutput.txt"), "w");
write(affineOutputFile,"instance, affineObj Tree, affineObj nontree, affine local tree, affine local nontree, affine local nontree time, numNodes", "\n");
for i=1:numInstances
  write(affineOutputFile, string(i,",",affineObj[i,1],",",affineObj[i,2],",",affineLocalObj[i,1],",",affineLocalObj[i,2],",",affineSolTime[i],","));
  write(affineOutputFile, string(numNodesList[i]));
  write(affineOutputFile, "\n");
end
close(affineOutputFile);

## Close all files
