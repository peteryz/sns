############### Helper functions ###############

function randPlusMinus1()
  number = rand()
  if number > 0.5
    return 1
  else
    return -1
  end
end


function makeOne(M) # return the 0/1 version of this matrix. Note: passing matrix by reference!
  nz = findn(M)
  for i=1:length(nz[1])
    M[nz[1][i],nz[2][i]] = 1;
  end
end

# Helper function for constructing path costs from arc costs.
# Pfull is a matrix (for k steps), Pij = 1 if there is a path from i to j within k steps, 0 otherwise.
# Afull is an adjacency matrix
# PC is the current cost matrix associated with Pfull
# After the function call, PC is the path matrix cost
function updateCost(P, A, PC, C)
  nzP = findn(P)
  nzA = findn(A)
  for p=1:length(nzP[1])
    for q=1:length(nzA[1])
      i = nzP[1][p];
      j = nzP[2][p];
      k = nzA[1][q];
      l = nzA[2][q];
      if j==k
        if PC[i,l] > 0 ## It's already initialized, safe to operate as is
          PC[i,l] = min(PC[i,l], PC[i,j]+C[k,l]);
        else ## path cost from i to l not initialized yet, so effectively infinite
          PC[i,l] = PC[i,j]+C[k,l];
        end
      end
    end
  end
end

# Tree creation
# Recursive function to add children
# resulting indices for nodes are depth first, need to sort later to put leaves at the end
# dict: parentIndex => children index array
# p: current parent index
# c: max number of children this parent has, rand(1:c)
# e: total number of nodes so far
# l: number of levels remaining to create
# leaf: current set of leaf indices
# return: leaf node set, end index, i.e., number of elements
function addC(dict, nl, p, c, e, l, leaf) # dict, parent index, mean number of child (rand 1:2c), last index so far, number of levels left
  nl[p] = l+1;
  if l > 0
    nc = rand(1:c);
    dict[p] = [e+1:e+nc];
    np = e;
    e = e+nc;
    for pp = np+1:np+nc
      (leaf, e) = addC(dict,nl,pp,c,e,l-1,leaf);
    end
  elseif l == 0
    push!(leaf, p)
  end
  return leaf, e;
end


# (numNodes, Atree, nodeLevel, pc, supplyNodes, demandNodes) = randomTreeByN(numNodesTarget, numChildrenPerLevel);
# New network generation function after April 17 meeting with Nikos (refer to scanned notes for details).
# The key is, given total number of nodes (target), and average number of children per node (level-dependent),
# generate a tree
# Atree - adjacency matrix.
# nodeLevel - hash table for looking up the level given a node index
# Input: target total number of nodes, f[i] - average children per node on level i-1
function randomTreeByN(N, f)
  L = size(f,1); # number of levels
  factors = zeros(L);
  factors[1] = 1;
  for i=2:L
    factors[i] = prod(f[1:i]);
  end

  if sum(factors) < N
    ub = 1;
    lb = (sum(factors)-1)/N;
  else
    lb = 1;
    ub = (sum(factors)-1)/N;
  end

  F = (ub+lb)/2;

  count = 1000;

  while count > 0 && abs(sum([factors[i]/F^(i-1) for i=1:L]) - N) > 1
    if sum([factors[i]/F^(i-1) for i=1:L]) > N
      lb = F;
    else
      ub = F;
    end
    F = (ub+lb)/2;
    count = count - 1;
  end

  for i=2:L
    f[i] = round(f[i]/F);
  end
  #println(string("Target number of nodes: ", N, ". Achievable number of nodes: ", sum([prod(f[1:i]) for i=1:L])));

  # now f is calibrated, construct pc and nodeLevel
  pc = Dict();
  nodeLevel = Dict();
  nodeLevel[1] = 1;

  q = [1]; # initialize queue with root node
  numNodes = 1;

  while !isempty(q)
    id = shift!(q);
    l = nodeLevel[id];
    if l+1 > L
      numKids = 0;
    else
      numKids = round(Int64, f[l+1] * (0.9 + 0.2 * rand()));
      kids = numNodes+1:numNodes+numKids;
      pc[id] = kids;
      append!(q,kids);
      for k in kids
        nodeLevel[k] = l+1;
      end
      numNodes = numNodes + numKids;
    end
  end

  Atree = spzeros(Int8,numNodes,numNodes); # tree adjacency matrix

  for i in keys(pc)
    c = pc[i];
    for k=1:size(c,1)
      Atree[i,c[k]] = 1;
    end
  end

  supplyNodes = zeros(Int8, numNodes);
  demandNodes = zeros(Int8, numNodes);

  for i=1:numNodes
    if get(pc,i,-1) != -1 # i is a parent
      supplyNodes[i] = 1;
    else
      demandNodes[i] = 1;
    end
  end

  return numNodes, Atree, nodeLevel, pc, supplyNodes, demandNodes;
end

# numChildren: max number of children per node. Random from 1 to numChildren
# numLevels: number of levels in this tree
function randomTree(numLevels, numChildren)
  # Create parent=>[children] dictionary
  parentChildren = Dict();
  nl = Dict(); # node index => level of node. leaf = 1, root = numLevels. Re-index later.
  leaf = Int16[]; # dynamic array to keep indices of leaves
  (leaf, numNodes) = addC(parentChildren, nl, 1, numChildren, 1, numLevels-1, leaf); # first 1: root index, seond 1: current total count of nodes
  invCost = zeros(numNodes);

  # Sort index in parentChildren so that increase with level.
  oldToNew = Dict(); # oldIndex => newIndex
  q = Int16[];
  sortedId = Int16[];
  push!(q,1); # root index is 1
  while !isempty(q)
    h = shift!(q);
    push!(sortedId,h);
    if get(parentChildren, h, -1) != -1
      append!(q,parentChildren[h]);
    end
  end
  for i=1:size(sortedId,1)
    oldToNew[sortedId[i]] = i;
  end

  nodeLevel = Dict();
  for i in keys(nl)
    newI = oldToNew[i];
    nodeLevel[newI] = nl[i];
  end

  # New parent to children dict based on sorted index
  pc = Dict();
  for i in keys(parentChildren)
    p = oldToNew[i];
    c = parentChildren[i];
    nc = Int16[];
    for j=1:size(c,1)
      push!(nc,oldToNew[c[j]])
    end
    pc[p] = nc;
  end

  Atree = spzeros(Int8,numNodes,numNodes); # tree adjacency matrix

  for i in keys(pc)
    c = pc[i];
    for k=1:size(c,1)
      Atree[i,c[k]] = 1;
    end
  end

  supplyNodes = zeros(Int8, numNodes);
  demandNodes = zeros(Int8, numNodes);

  for i=1:numNodes
    if get(pc,i,-1) != -1 # i is a parent
      supplyNodes[i] = 1;
    else
      demandNodes[i] = 1;
    end
  end

  return numNodes, Atree, nodeLevel, pc, supplyNodes, demandNodes;
end

# Add extra arcs to the pc dictionary according to probability p
function extraArcs(pc, p, n)
  for i in keys(pc)
    c = pc[i];
    for j in c
      for k in c
        if j != k && get(pc,k,-1) != -1 && get(pc,j,-1) != -1
          kc = pc[k];
          for l in kc
            if rand() < p
              jc = pc[j];
              push!(jc,l);
              pc[j] = jc;
            end
          end
        end
      end
    end
  end

  Afull = spzeros(Int8,n,n); # tree adjacency matrix

  for i in keys(pc)
    c = pc[i];
    for k=1:size(c,1)
      Afull[i,c[k]] = 1;
    end
  end

  return Afull;
end


# solving affine formulation assuming A (and P? maybe not  necessary) describe a tree
function affine(n, Γ, Atree, P, h, b, c, dmean, dvar, F0return, Fdreturn, xreturn, solveTime, instance)

  # All 1's mimics the total budget=1 scenario
  # Γ1=k, downstream all 1, mimics the total budget=k scenario
  #Γ1 = Γs[1]
  #Γ2 = Γs[2]
  #Γ3 = Γs[3]
  #Γ4 = Γs[4]
  #Γ5 = Γs[5]

  # We need to use RobustModel, not Model
  model = RobustModel(solver=GurobiSolver(Method=2, Crossover=0, BarConvTol=0.000001));
  #model = RobustModel();

  # Set up variables
  @variable(model, x[1:n] >= 0); # inventory decisions
  FdIndices = String[]
  F0Indices = String[]
  for i=1:n
    for j=1:n
      if P[i,j] == 1
        push!(F0Indices,"$(i)_$(j)")
        for k=1:n
          push!(FdIndices,"$(i)_$(j)_$(k)")
        end
      end
    end
  end
  @variable(model, F0[F0Indices]) # offset terms for affine policy flows
  @variable(model, Fd[FdIndices]) # affine coefficients for flows; first dimension is for the uncertainty
  @variable(model, S0[1:n]); # offset terms for affine policy flows
  @variable(model, Sd[1:n,1:n]); # affine coefficients for flows
  @variable(model, z >= 0); # second stage cost for epigraph format

  # Set up uncertainty polytope
  @uncertain(model, ξ[1:n])
  @constraint(model, ξ[1] >= 1);
  @constraint(model, ξ[1] <= 1);
  for i=1:n
    @constraint(model, ξ[i] >= 0);
    @constraint(model, ξ[i] <= 1);
    for j=1:n
      if (Atree[i,j] > 0) # A is a tree
        @constraint(model, ξ[i] >= ξ[j]);
      end
    end
    if (sum(Atree[i,:]) > 0) # this node has children. A is a tree
      @constraint(model, sum([ξ[j]*Atree[i,j] for j=1:n]) <= Γ[i] * ξ[i]); # A is a tree
    end
  end

  # flow constraints
  for j=1:n
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=1:n]) + F0["$(i)_$(j)"] : 0 ) for i=1:n]) + sum([ξ[s]*Sd[s,j] for s=1:n]) + S0[j] >= dmean[j] + dvar[j] * ξ[j] + sum([(P[j,k] == 1 ? sum([ξ[s]*Fd["$(j)_$(k)_$(s)"] for s=1:n]) + F0["$(j)_$(k)"] : 0) for k=1:n]))
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=1:n]) + F0["$(i)_$(j)"] : 0 ) for i=1:n]) >= sum([ (P[j,k] == 1 ? sum([ξ[s]*Fd["$(j)_$(k)_$(s)"] for s=1:n]) + F0["$(j)_$(k)"] : 0) for k=1:n]))
  end

  for j=1:n
    for i = 1:n
      if P[i,j] == 1
        @constraint(model, sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=1:n]) + F0["$(i)_$(j)"] >= 0)
      end
    end
    @constraint(model, sum([ξ[s]*Sd[s,j] for s=1:n]) + S0[j] >= 0)
  end

  #### YES delay penalty ####
  @constraint(model, z >= sum([(sum([ξ[s]*Sd[s,j] for s=1:n]) + S0[j]) * b[j] for j=1:n]) +  sum([sum([  (P[i,j]==1 ? (sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=1:n]) + F0["$(i)_$(j)"])*c[i,j] : 0) for j=1:n]) for i=1:n])) # with delay penalty

  #### NO delay penalty ####
  #@constraint(model, z >= sum([(sum([γ[s]*Sd[s,j] for s=1:n]) + S0[j]) * b[j] for j=1:n])) # without delay penalty


  # We will use the non-macro version of setObjective, although
  # we could use the macro version if we wanted to.
  # Remember: objectives must be certain.
  @objective(model, :Min, sum([h[i]*x[i] for i=1:n]) + z)

  tic()
  status = solve(model)
  solveTime[instance] = toc()

#  for i=1:n
#    for j=1:n
#      if P[i,j] == 1
#        F0return[i,j] = getvalue(F0["$(i)_$(j)"])
#        for k=1:n
#          Fdreturn[i,j,k] = getvalue(Fd["$(i)_$(j)_$(k)"])
#        end
#      end
#    end
#  end

  for i=1:n
    xreturn[i] = getvalue(x[i])
  end

  aObj = getobjectivevalue(model)
  return sum([h[i]*getvalue(x[i]) for i=1:n]) + getvalue(z)
end


# fully adaptive policy by scenario enumeration
function fullyAdapt(n, Γs, P, h, b, c, dmean, dvar, xreturn, solveTime)

  # We need to use RobustModel, not Model
  model = Model();

  Γ1 = Γs[1]
  Γ2 = Γs[2]
  Γ3 = Γs[3]
  Γ4 = Γs[4]
  Γ5 = Γs[5]
  Γsum = Γ1+Γ2+Γ3+Γ4+Γ5
  ######### If change the budgets, also have to change line 117, combinations(1:20,sum(gammas))

  # K is the number of scenarios; we are iterating through all the scenarios for the full adaptibility solution
  K = multinomial(Γ1,4-Γ1)*multinomial(Γ2,4-Γ2)*multinomial(Γ3,4-Γ3)*multinomial(Γ4,4-Γ4)*multinomial(Γ5,4-Γ5)

  # Set up variables
  @variable(model, x[1:n] >= 0); # inventory decisions
  #@variable(model, f[1:n,1:n,1:K] >= 0); # offset terms for affine policy flows
  fIndices = String[]
  for i=1:n
    for j=1:n
      if P[i,j] == 1
        for k=1:K
          push!(fIndices,"$(i)_$(j)_$(k)")
        end
      end
    end
  end
  @variable(model, f[fIndices] >= 0) # offset terms for affine policy flows
  @variable(model, s[1:n,1:K] >= 0); # offset terms for affine policy flows
  @variable(model, z >= 0); # second stage cost for epigraph format

  # flow constraints
  indx = 0
  for scenario in combinations(1:20,Γsum) # master list of scenarios; may not be a valid scenario

    γ = zeros(n)
    γ[1+scenario] = 1

    ## we do not immediately enforce the hierarchical consistence (e.g., γ[6] <= γ[2]) because we are slightly abusing the notation for γ[i]; a scenario on leaf i
    ## is active only when γ[i] * γ[p_i] = 1
    if (sum([γ[i] for i=2:5]) == Γ1) && (sum([γ[i] for i=6:9]) == Γ2) && (sum([γ[i] for i=10:13]) == Γ3) && (sum([γ[i] for i=14:17]) == Γ4) && (sum([γ[i] for i=18:21]) == Γ5)
      indx = indx + 1
      for j=1:n
        #if sum(P[j,:]) > 0
        #  @constraint(model, x[j] + sum([f[i,j,indx] * P[i,j] for i=1:n]) >= sum([P[j,k] * f[j,k,indx] for k=1:n]))
        #else
        #  @constraint(model, x[j] + sum([f[i,j,indx] * P[i,j] for i=1:n]) + s[j,indx] >= dmean[j] + dvar[j] * γ[j] * γ[ceil((j-5)/4)+1])
        #end
        @constraint(model, x[j] + sum([(P[i,j] == 1 ? f["$(i)_$(j)_$(indx)"] : 0) for i=1:n]) >= sum([(P[j,k] == 1 ? f["$(j)_$(k)_$(indx)"] : 0) for k=1:n]))
        if j >= 6
          @constraint(model, x[j] + sum([(P[i,j] == 1 ? f["$(i)_$(j)_$(indx)"] : 0) for i=1:n]) + s[j,indx] >= dmean[j] + dvar[j] * γ[j] * γ[ceil((j-5)/4)+1] + sum([(P[j,k] == 1 ? f["$(j)_$(k)_$(indx)"] : 0) for k=1:n]))
        end
      end
      #### YES delay penalty ####
      @constraint(model, z >= sum([s[j,indx] * b[j] for j=1:n]) + sum([ sum([(P[i,j] == 1 ? f["$(i)_$(j)_$(indx)"]*c[i,j] : 0) for j=1:n]) for i=1:n])) # with delay penalty

      #### NO delay penalty ####
      #@constraint(model, z >= sum([s[j,indx] * b[j] for j=1:n])) # without delay penalty
    end
  end

  # We will use the non-macro version of setObjective, although
  # we could use the macro version if we wanted to.
  # Remember: objectives must be certain.
  @objective(model, :Min, sum([h[i]*x[i] for i=1:n]) + z)

  tic()
  status = solve(model)
  solveTime[1] = toc()

  for i=1:n
    xreturn[i] = getvalue(x[i])
  end

  fObj = getobjectivevalue(model)
  return sum([h[i]*getvalue(x[i]) for i=1:n]) + getvalue(z)
end

# Local: flow ij is a function of \xi_j only
function affineLocal(n, Γ, Atree, P, h, b, c, dmean, dvar, F0return, Fdreturn, xreturn, solveTime, instance)

  # All 1's mimics the total budget=1 scenario
  # Γ1=k, downstream all 1, mimics the total budget=k scenario
  #Γ1 = Γs[1]
  #Γ2 = Γs[2]
  #Γ3 = Γs[3]
  #Γ4 = Γs[4]
  #Γ5 = Γs[5]

  # We need to use RobustModel, not Model
  #model = RobustModel();
  model = RobustModel(solver=GurobiSolver(Method=2, Crossover=0, BarConvTol=0.000001));

  # Set up variables
  @variable(model, x[1:n] >= 0); # inventory decisions
  FdIndices = String[]
  F0Indices = String[]
  for i=1:n
    for j=1:n
      if P[i,j] == 1
        push!(F0Indices,"$(i)_$(j)")
        for k=j # affine local policy
          push!(FdIndices,"$(i)_$(j)_$(k)")
        end
      end
    end
  end
  @variable(model, F0[F0Indices]) # offset terms for affine policy flows
  @variable(model, Fd[FdIndices]) # affine coefficients for flows; first dimension is for the uncertainty
  @variable(model, S0[1:n]); # offset terms for affine policy flows
  @variable(model, Sd[1:n,1:n]); # affine coefficients for flows
  @variable(model, z >= 0); # second stage cost for epigraph format

  # Set up uncertainty polytope
  @uncertain(model, ξ[1:n])
  @constraint(model, ξ[1] >= 1);
  @constraint(model, ξ[1] <= 1);
  for i=1:n
    @constraint(model, ξ[i] >= 0);
    @constraint(model, ξ[i] <= 1);
    for j=1:n
      if (Atree[i,j] > 0) # A is a tree
        @constraint(model, ξ[i] >= ξ[j]);
      end
    end
    if (sum(Atree[i,:]) > 0) # this node has children. A is a tree
      @constraint(model, sum([ξ[j]*Atree[i,j] for j=1:n]) <= Γ[i] * ξ[i]); # A is a tree
    end
  end

  # flow constraints
  for j=1:n
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=j]) + F0["$(i)_$(j)"] : 0 ) for i=1:n]) + sum([ξ[s]*Sd[s,j] for s=1:n]) + S0[j] >= dmean[j] + dvar[j] * ξ[j] + sum([(P[j,k] == 1 ? sum([ξ[s]*Fd["$(j)_$(k)_$(s)"] for s=k]) + F0["$(j)_$(k)"] : 0) for k=1:n]))
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=j]) + F0["$(i)_$(j)"] : 0 ) for i=1:n]) >= sum([ (P[j,k] == 1 ? sum([ξ[s]*Fd["$(j)_$(k)_$(s)"] for s=k]) + F0["$(j)_$(k)"] : 0) for k=1:n]))
  end

  for j=1:n
    for i = 1:n
      if P[i,j] == 1
        @constraint(model, sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=j]) + F0["$(i)_$(j)"] >= 0)
      end
    end
    @constraint(model, sum([ξ[s]*Sd[s,j] for s=1:n]) + S0[j] >= 0)
  end

  #### YES delay penalty ####
  @constraint(model, z >= sum([(sum([ξ[s]*Sd[s,j] for s=1:n]) + S0[j]) * b[j] for j=1:n]) +  sum([sum([  (P[i,j]==1 ? (sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=j]) + F0["$(i)_$(j)"])*c[i,j] : 0) for j=1:n]) for i=1:n])) # with delay penalty

  #### NO delay penalty ####
  #addConstraint(model, z >= sum([(sum([γ[s]*Sd[s,j] for s=1:n]) + S0[j]) * b[j] for j=1:n])) # without delay penalty


  # We will use the non-macro version of setObjective, although
  # we could use the macro version if we wanted to.
  # Remember: objectives must be certain.
  @objective(model, :Min, sum([h[i]*x[i] for i=1:n]) + z)

  tic()
  status = solve(model)
  solveTime[instance] = toc()

#  for i=1:n
#    for j=1:n
#      if P[i,j] == 1
#        F0return[i,j] = getvalue(F0["$(i)_$(j)"])
#        for k=j
#          Fdreturn[i,j,k] = getvalue(Fd["$(i)_$(j)_$(k)"])
#        end
#      end
#    end
#  end

  for i=1:n
    xreturn[i] = getvalue(x[i])
  end

  aObj = getobjectivevalue(model);
  return sum([h[i]*getvalue(x[i]) for i=1:n]) + getvalue(z)
end

# Path-aware implementation, reducing the number of constraints and constraint lengths
# Refer to April 19 scanned notes in the SNS revision folder.
# nS is the number of supply nodes, assumed to have ID smaller than all demand nodes
function affinePath(n, nS, Γ, Atree, P, h, b, c, dmean, dvar, F0return, Fdreturn, xreturn, solveTime, instance)
  model = RobustModel(solver=GurobiSolver(Method=2, Crossover=0, BarConvTol=0.000001));
  #model = RobustModel();

  # Set up variables
  @variable(model, x[1:n] >= 0); # inventory decisions
  FdIndices = String[]
  F0Indices = String[]
  for i=1:n
    for j=1:n
      if P[i,j] == 1
        push!(F0Indices,"$(i)_$(j)")
        for k=1:n
          push!(FdIndices,"$(i)_$(j)_$(k)")
        end
      end
    end
  end
  @variable(model, F0[F0Indices]) # offset terms for affine policy flows
  @variable(model, Fd[FdIndices]) # affine coefficients for flows; first dimension is for the uncertainty
  @variable(model, S0[1:n]); # offset terms for affine policy flows
  @variable(model, Sd[1:n,1:n]); # affine coefficients for flows
  @variable(model, z >= 0); # second stage cost for epigraph format

  # Set up uncertainty polytope
  @uncertain(model, ξ[1:n])
  @constraint(model, ξ[1] >= 1);
  @constraint(model, ξ[1] <= 1);
  for i=1:n
    @constraint(model, ξ[i] >= 0);
    @constraint(model, ξ[i] <= 1);
    for j=1:n
      if (Atree[i,j] > 0) # A is a tree
        @constraint(model, ξ[i] >= ξ[j]);
      end
    end
    if (sum(Atree[i,:]) > 0) # this node has children. A is a tree
      @constraint(model, sum([ξ[j]*Atree[i,j] for j=1:n]) <= Γ[i] * ξ[i]); # A is a tree
    end
  end

  # demand satisfaction constraints:
  for j=(nS+1):n
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=1:n]) + F0["$(i)_$(j)"] : 0 ) for i=1:n]) + sum([ξ[s]*Sd[s,j] for s=1:n]) + S0[j] >= dmean[j] + dvar[j] * ξ[j])
  end

  # supply node capacity/conservation constraints:
  for i=1:n
    @constraint(model, x[i] >= sum([(P[i,j] == 1 ? sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=1:n]) + F0["$(i)_$(j)"] : 0) for j=(nS+1):n]))
  end

  for j=1:n
    for i = 1:n
      if P[i,j] == 1
        @constraint(model, sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=1:n]) + F0["$(i)_$(j)"] >= 0)
      end
    end
    @constraint(model, sum([ξ[s]*Sd[s,j] for s=1:n]) + S0[j] >= 0)
  end

  #### YES delay penalty ####
  @constraint(model, z >= sum([(sum([ξ[s]*Sd[s,j] for s=1:n]) + S0[j]) * b[j] for j=(nS+1):n]) +  sum([sum([(P[i,j]==1 ? (sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=1:n]) + F0["$(i)_$(j)"])*c[i,j] : 0) for j=(nS+1):n]) for i=1:n])) # with delay penalty

  # We will use the non-macro version of setObjective, although
  # we could use the macro version if we wanted to.
  # Remember: objectives must be certain.
  @objective(model, :Min, sum([h[i]*x[i] for i=1:n]) + z)

  tic()
  status = solve(model)
  solveTime[instance] = toc()

#  for i=1:n
#    for j=1:n
#      if P[i,j] == 1
#        F0return[i,j] = getvalue(F0["$(i)_$(j)"])
#        for k=1:n
#          Fdreturn[i,j,k] = getvalue(Fd["$(i)_$(j)_$(k)"])
#        end
#      end
#    end
#  end

  for i=1:n
    xreturn[i] = getvalue(x[i])
  end

  aObj = getobjectivevalue(model)
  return sum([h[i]*getvalue(x[i]) for i=1:n]) + getvalue(z)
end


# Path-aware implementation, reducing the number of constraints and constraint lengths
# Refer to April 19 scanned notes in the SNS revision folder.
# nS is the number of supply nodes, assumed to have ID smaller than all demand nodes
function affinePathLocal(n, nS, Γ, Atree, P, h, b, c, dmean, dvar, F0return, Fdreturn, xreturn, solveTime, instance)
  model = RobustModel(solver=GurobiSolver(Method=2, Crossover=0, BarConvTol=0.000001));
  #model = RobustModel();

  # Set up variables
  @variable(model, x[1:n] >= 0); # inventory decisions
  FdIndices = String[]
  F0Indices = String[]
  for i=1:n
    for j=1:n
      if P[i,j] == 1
        push!(F0Indices,"$(i)_$(j)")
        for k=j
          push!(FdIndices,"$(i)_$(j)_$(k)")
        end
      end
    end
  end
  @variable(model, F0[F0Indices]) # offset terms for affine policy flows
  @variable(model, Fd[FdIndices]) # affine coefficients for flows; first dimension is for the uncertainty
  @variable(model, S0[1:n]); # offset terms for affine policy flows
  @variable(model, Sd[1:n,1:n]); # affine coefficients for flows
  @variable(model, z >= 0); # second stage cost for epigraph format

  # Set up uncertainty polytope
  @uncertain(model, ξ[1:n])
  @constraint(model, ξ[1] >= 1);
  @constraint(model, ξ[1] <= 1);
  for i=1:n
    @constraint(model, ξ[i] >= 0);
    @constraint(model, ξ[i] <= 1);
    for j=1:n
      if (Atree[i,j] > 0) # A is a tree
        @constraint(model, ξ[i] >= ξ[j]);
      end
    end
    if (sum(Atree[i,:]) > 0) # this node has children. A is a tree
      @constraint(model, sum([ξ[j]*Atree[i,j] for j=1:n]) <= Γ[i] * ξ[i]); # A is a tree
    end
  end

  # demand satisfaction constraints:
  for j=(nS+1):n
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=j]) + F0["$(i)_$(j)"] : 0 ) for i=1:n]) + sum([ξ[s]*Sd[s,j] for s=1:n]) + S0[j] >= dmean[j] + dvar[j] * ξ[j])
  end

  # supply node capacity/conservation constraints:
  for i=1:n
    @constraint(model, x[i] >= sum([(P[i,j] == 1 ? sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=j]) + F0["$(i)_$(j)"] : 0) for j=(nS+1):n]))
  end

  for j=1:n
    for i = 1:n
      if P[i,j] == 1
        @constraint(model, sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=j]) + F0["$(i)_$(j)"] >= 0)
      end
    end
    @constraint(model, sum([ξ[s]*Sd[s,j] for s=1:n]) + S0[j] >= 0)
  end

  #### YES delay penalty ####
  @constraint(model, z >= sum([(sum([ξ[s]*Sd[s,j] for s=1:n]) + S0[j]) * b[j] for j=(nS+1):n]) +  sum([sum([(P[i,j]==1 ? (sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=j]) + F0["$(i)_$(j)"])*c[i,j] : 0) for j=(nS+1):n]) for i=1:n])) # with delay penalty

  # We will use the non-macro version of setObjective, although
  # we could use the macro version if we wanted to.
  # Remember: objectives must be certain.
  @objective(model, :Min, sum([h[i]*x[i] for i=1:n]) + z)

  tic()
  status = solve(model)
  solveTime[instance] = toc()

#  for i=1:n
#    for j=1:n
#      if P[i,j] == 1
#        F0return[i,j] = getvalue(F0["$(i)_$(j)"])
#        for k=1:n
#          Fdreturn[i,j,k] = getvalue(Fd["$(i)_$(j)_$(k)"])
#        end
#      end
#    end
#  end

  for i=1:n
    xreturn[i] = getvalue(x[i])
  end

  aObj = getobjectivevalue(model)
  return sum([h[i]*getvalue(x[i]) for i=1:n]) + getvalue(z)
end


# Path-aware implementation, reducing the number of constraints and constraint lengths
# Refer to April 19 scanned notes in the SNS revision folder.
# nS is the number of supply nodes, assumed to have ID smaller than all demand nodes
# April 21 improvement: localize everything including s, and use deterministic version of demand when possible (April 20 scanned notes)
# (a) Only index flow i,j if P[i,j] == 1 AND j is in demand nodes
# (b) index loss s by string too, not 2D array with n^2 vars
function affinePathLocalFS(n, nS, Γ, Atree, P, h, b, c, dmean, dvar, F0return, Fdreturn, xreturn, solveTime, instance)
  model = RobustModel(solver=GurobiSolver(Method=2, Crossover=0, BarConvTol=0.000001));
  #model = RobustModel();

  # Set up variables
  @variable(model, x[1:n] >= 0); # inventory decisions
  FdIndices = String[]
  F0Indices = String[]
  SdIndices = String[]
  S0Indices = String[]
  for j=(nS+1):n # for every demand node
    push!(S0Indices,"$(j)")
    push!(SdIndices,"$(j)")
    for i=1:n
      if P[i,j] == 1
        push!(F0Indices,"$(i)_$(j)")
        for k=j
          push!(FdIndices,"$(i)_$(j)_$(k)")
        end
      end
    end
  end
  @variable(model, F0[F0Indices]) # offset terms for affine policy flows
  @variable(model, Fd[FdIndices]) # affine coefficients for flows; first dimension is for the uncertainty
  @variable(model, S0[S0Indices]); # offset terms for affine policy flows
  @variable(model, Sd[SdIndices]); # affine coefficients for flows
  @variable(model, z >= 0); # second stage cost for epigraph format

  # Set up uncertainty polytope
  @uncertain(model, ξ[1:n])
  @constraint(model, ξ[1] >= 1);
  @constraint(model, ξ[1] <= 1);
  for i=1:n
    @constraint(model, ξ[i] >= 0);
    @constraint(model, ξ[i] <= 1);
    for j=1:n
      if (Atree[i,j] > 0) # A is a tree
        @constraint(model, ξ[i] >= ξ[j]);
      end
    end
    if (sum(Atree[i,:]) > 0) # this node has children. A is a tree
      @constraint(model, sum([ξ[j]*Atree[i,j] for j=1:n]) <= Γ[i] * ξ[i]); # A is a tree
    end
  end

  # demand satisfaction constraints, deterministic version since only involve \xi_j
  for j=(nS+1):n
    # ξ_j = 1
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? Fd["$(i)_$(j)_$(j)"] + F0["$(i)_$(j)"] : 0 ) for i=1:n]) + Sd["$(j)"] + S0["$(j)"] >= dmean[j] + dvar[j])
    # ξ_j = 0
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? F0["$(i)_$(j)"] : 0 ) for i=1:n]) + S0["$(j)"] >= dmean[j])
  end

  # supply node capacity/conservation constraints:
  for i=1:nS
    @constraint(model, x[i] >= sum([(P[i,j] == 1 ? sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=j]) + F0["$(i)_$(j)"] : 0) for j=(nS+1):n]))
  end

  for j=(nS+1):n
    for i = 1:n
      if P[i,j] == 1
        @constraint(model, Fd["$(i)_$(j)_$(j)"] + F0["$(i)_$(j)"] >= 0) # deterministic version
        @constraint(model, F0["$(i)_$(j)"] >= 0)
      end
    end
    @constraint(model, Sd["$(j)"] + S0["$(j)"] >= 0)
    @constraint(model, S0["$(j)"] >= 0)
  end

  #### YES delay penalty ####
  @constraint(model, z >= sum([(ξ[j]*Sd["$(j)"] + S0["$(j)"]) * b[j] for j=(nS+1):n]) +  sum([sum([(P[i,j]==1 ? (ξ[j]*Fd["$(i)_$(j)_$(j)"] + F0["$(i)_$(j)"])*c[i,j] : 0) for j=(nS+1):n]) for i=1:n])) # with delay penalty

  # We will use the non-macro version of setObjective, although
  # we could use the macro version if we wanted to.
  # Remember: objectives must be certain.
  @objective(model, :Min, sum([h[i]*x[i] for i=1:n]) + z)

  tic()
  status = solve(model)
  solveTime[instance] = toc()

#  for i=1:n
#    for j=1:n
#      if P[i,j] == 1
#        F0return[i,j] = getvalue(F0["$(i)_$(j)"])
#        for k=1:n
#          Fdreturn[i,j,k] = getvalue(Fd["$(i)_$(j)_$(k)"])
#        end
#      end
#    end
#  end

  for i=1:n
    xreturn[i] = getvalue(x[i])
  end

  aObj = getobjectivevalue(model)
  return sum([h[i]*getvalue(x[i]) for i=1:n]) + getvalue(z)
end


# Based on affinePathLocalFS, with the additional feature of dualizing the node inventory capacity x >= sum(f) constraints
# Refer to scanned notes April 21, 2018
# Additional input argument: all based on tree
# Ptree,
# parent: Dict(), kid -> parent
# kids: Dict(), parent -> [kids]
# offspr: Dict(), root -> [offsprings]
# offsprNr: contains root too
#
function affineXiTree(n, nS, Γ, Atree, P, Ptree, parent, kids, offspr, offsprNr, h, b, c, dmean, dvar, F0return, Fdreturn, xreturn, solveTime, instance)
  model = RobustModel(solver=GurobiSolver(Method=2, Crossover=0, BarConvTol=0.000001));
  #model = RobustModel();

  # Set up variables
  @variable(model, x[1:n] >= 0); # inventory decisions
  FdIndices = String[]
  F0Indices = String[]
  SdIndices = String[]
  S0Indices = String[]
  for j=(nS+1):n # for every demand node
    push!(S0Indices,"$(j)")
    push!(SdIndices,"$(j)")
    for i=1:n
      if P[i,j] == 1
        push!(F0Indices,"$(i)_$(j)")
        for k=j
          push!(FdIndices,"$(i)_$(j)_$(k)")
        end
      end
    end
  end
  @variable(model, F0[F0Indices]) # offset terms for affine policy flows
  @variable(model, Fd[FdIndices]) # affine coefficients for flows; first dimension is for the uncertainty
  @variable(model, S0[S0Indices]); # offset terms for affine policy flows
  @variable(model, Sd[SdIndices]); # affine coefficients for flows
  @variable(model, z >= 0); # second stage cost for epigraph format

  # dual variables for flow constraints
  @variable(model, p[1:n]) # 1:nS if demand nodes cannot serve demand nodes.
  qIndices = String[]
  rIndices = String[]
  for t=1:n
    for i in offsprNr[t]
      push!(rIndices,"$(t)_$(i)")
      for j in get(kids,i,[])
        push!(qIndices, "$(t)_$(i)_$(j)")
      end
    end
  end
  @variable(model, q[qIndices] >= 0)
  @variable(model, r[rIndices] >= 0)

  # Set up uncertainty polytope
  @uncertain(model, ξ[1:n])
  @constraint(model, ξ[1] >= 1);
  @constraint(model, ξ[1] <= 1);
  for i=1:n
    @constraint(model, ξ[i] >= 0);
    @constraint(model, ξ[i] <= 1);
    for j=1:n
      if (Atree[i,j] > 0) # A is a tree
        @constraint(model, ξ[i] >= ξ[j]);
      end
    end
    if (sum(Atree[i,:]) > 0) # this node has children. A is a tree
      @constraint(model, sum([ξ[j]*Atree[i,j] for j=1:n]) <= Γ[i] * ξ[i]); # A is a tree
    end
  end

  # demand satisfaction constraints, deterministic version since only involve \xi_j
  for j=(nS+1):n
    # ξ_j = 1
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? Fd["$(i)_$(j)_$(j)"] + F0["$(i)_$(j)"] : 0 ) for i=1:n]) + Sd["$(j)"] + S0["$(j)"] >= dmean[j] + dvar[j])
    # ξ_j = 0
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? F0["$(i)_$(j)"] : 0 ) for i=1:n]) + S0["$(j)"] >= dmean[j])
  end

  # supply node capacity/conservation constraints:
  #for i=1:n
  #  @constraint(model, x[i] >= sum([(P[i,j] == 1 ? sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=j]) + F0["$(i)_$(j)"] : 0) for j=(nS+1):n]))
  #end

  # robustified supply node capacity constraints. Refer to scanned notes April 21.
  for t=1:nS
    @constraint(model, p[t] <= x[t] - sum([(j > nS ? F0["$(t)_$(j)"] : 0) for j in get(offspr,t,[])]))
    @constraint(model, p[t] - sum([q["$(t)_$(t)_$(j)"] for j in get(kids,t,[])]) - Γ[t]*r["$(t)_$(t)"] >= 0)
    for i in get(offspr,t,[])
      if i <= nS
        @constraint(model, q["$(t)_$(parent[i])_$(i)"] - sum([q["$(t)_$(i)_$(j)"] for j in get(kids,i,[])]) + r["$(t)_$(parent[i])"] - Γ[i]*r["$(t)_$(i)"] >= 0)
      else
        @constraint(model, q["$(t)_$(parent[i])_$(i)"] + r["$(t)_$(parent[i])"] >= Fd["$(t)_$(i)_$(i)"])
      end
    end
  end

  # non-negativity of flow and loss variables
  # deterministic version given 0 <= ξ[i] <= 1, and worst case happens at extreme point
  for j=(nS+1):n
    for i = 1:n
      if P[i,j] == 1
        @constraint(model, Fd["$(i)_$(j)_$(j)"] + F0["$(i)_$(j)"] >= 0)
        @constraint(model, F0["$(i)_$(j)"] >= 0)
      end
    end
    @constraint(model, Sd["$(j)"] + S0["$(j)"] >= 0)
    @constraint(model, S0["$(j)"] >= 0)
  end

  #### YES delay penalty ####
  @constraint(model, z >= sum([(ξ[j]*Sd["$(j)"] + S0["$(j)"]) * b[j] for j=(nS+1):n]) +  sum([sum([(P[i,j]==1 ? (ξ[j]*Fd["$(i)_$(j)_$(j)"] + F0["$(i)_$(j)"])*c[i,j] : 0) for j=(nS+1):n]) for i=1:n])) # with delay penalty

  # We will use the non-macro version of setObjective, although
  # we could use the macro version if we wanted to.
  # Remember: objectives must be certain.
  @objective(model, :Min, sum([h[i]*x[i] for i=1:n]) + z)

  tic()
  status = solve(model)
  solveTime[instance] = toc()

#  for i=1:n
#    for j=1:n
#      if P[i,j] == 1
#        F0return[i,j] = getvalue(F0["$(i)_$(j)"])
#        for k=1:n
#          Fdreturn[i,j,k] = getvalue(Fd["$(i)_$(j)_$(k)"])
#        end
#      end
#    end
#  end

  for i=1:n
    xreturn[i] = getvalue(x[i])
  end

  aObj = getobjectivevalue(model)
  return sum([h[i]*getvalue(x[i]) for i=1:n]) + getvalue(z)
end



# Based on affineXiTree, now for LLG, life loss guarantee formulation.
# PDF Notes: 20180512 affineXiTreeLLG
function affineXiTreeLLG(n, nS, Γ, Atree, P, Ptree, parent, kids, offspr, offsprNr, h, b, c, rho, surv, dmean, dvar, F0return, Fdreturn, xreturn, solveTime, instance)
  model = RobustModel(solver=GurobiSolver(Method=2, Crossover=0, BarConvTol=0.000001));
  #model = RobustModel();

  # Set up variables
  @variable(model, x[1:n] >= 0); # inventory decisions
  FdIndices = String[]
  F0Indices = String[]
  SdIndices = String[]
  S0Indices = String[]
  for j=(nS+1):n # for every demand node
    push!(S0Indices,"$(j)")
    push!(SdIndices,"$(j)")
    for i=1:n
      if P[i,j] == 1
        push!(F0Indices,"$(i)_$(j)")
        for k=j
          push!(FdIndices,"$(i)_$(j)_$(k)")
        end
      end
    end
  end
  @variable(model, F0[F0Indices]) # offset terms for affine policy flows
  @variable(model, Fd[FdIndices]) # affine coefficients for flows; first dimension is for the uncertainty
  @variable(model, S0[S0Indices]); # offset terms for affine policy flows
  @variable(model, Sd[SdIndices]); # affine coefficients for flows
  @variable(model, z >= 0); # second stage cost for epigraph format

  # dual variables for flow constraints
  @variable(model, p[1:n]) # 1:nS if demand nodes cannot serve demand nodes.
  qIndices = String[]
  rIndices = String[]
  for t=1:n
    for i in offsprNr[t]
      push!(rIndices,"$(t)_$(i)")
      for j in get(kids,i,[])
        push!(qIndices, "$(t)_$(i)_$(j)")
      end
    end
  end
  @variable(model, q[qIndices] >= 0)
  @variable(model, r[rIndices] >= 0)

  # Set up uncertainty polytope
  @uncertain(model, ξ[1:n])
  @constraint(model, ξ[1] >= 1);
  @constraint(model, ξ[1] <= 1);
  for i=1:n
    @constraint(model, ξ[i] >= 0);
    @constraint(model, ξ[i] <= 1);
    for j=1:n
      if (Atree[i,j] > 0) # A is a tree
        @constraint(model, ξ[i] >= ξ[j]);
      end
    end
    if (sum(Atree[i,:]) > 0) # this node has children. A is a tree
      @constraint(model, sum([ξ[j]*Atree[i,j] for j=1:n]) <= Γ[i] * ξ[i]); # A is a tree
    end
  end

  # demand satisfaction constraints, deterministic version since only involve \xi_j (constraint 2a and 2b in scanned notes)
  for j=(nS+1):n
    # ξ_j = 1
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? (Fd["$(i)_$(j)_$(j)"] + F0["$(i)_$(j)"])*rho[i,j] : 0 ) for i=1:n]) + 0*(Sd["$(j)"] + S0["$(j)"]) >= surv[j] * (dmean[j] + dvar[j]))
    # ξ_j = 0
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? F0["$(i)_$(j)"]*rho[i,j] : 0 ) for i=1:n]) + 0*S0["$(j)"] >= surv[j] * dmean[j])
  end

  # flow cap constraint (constraints 2.1 in scanned notes)
  for j=(nS+1):n
    # ξ_j = 1
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? Fd["$(i)_$(j)_$(j)"] + F0["$(i)_$(j)"] : 0 ) for i=1:n]) <= dmean[j] + dvar[j])
    # ξ_j = 0
    @constraint(model,        sum([(P[i,j] == 1 ? F0["$(i)_$(j)"] : 0 ) for i=1:n]) <= dmean[j])
  end

  # supply node capacity/conservation constraints:
  #for i=1:n
  #  @constraint(model, x[i] >= sum([(P[i,j] == 1 ? sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=j]) + F0["$(i)_$(j)"] : 0) for j=(nS+1):n]))
  #end

  # robustified supply node capacity constraints. Refer to scanned notes April 21.
  for t=1:nS
    @constraint(model, p[t] <= x[t] - sum([(j > nS ? F0["$(t)_$(j)"] : 0) for j in get(offspr,t,[])]))
    @constraint(model, p[t] - sum([q["$(t)_$(t)_$(j)"] for j in get(kids,t,[])]) - Γ[t]*r["$(t)_$(t)"] >= 0)
    for i in get(offspr,t,[])
      if i <= nS
        @constraint(model, q["$(t)_$(parent[i])_$(i)"] - sum([q["$(t)_$(i)_$(j)"] for j in get(kids,i,[])]) + r["$(t)_$(parent[i])"] - Γ[i]*r["$(t)_$(i)"] >= 0)
      else
        @constraint(model, q["$(t)_$(parent[i])_$(i)"] + r["$(t)_$(parent[i])"] >= Fd["$(t)_$(i)_$(i)"])
      end
    end
  end

  # non-negativity of flow and loss variables
  # deterministic version given 0 <= ξ[i] <= 1, and worst case happens at extreme point
  for j=(nS+1):n
    for i = 1:n
      if P[i,j] == 1
        @constraint(model, Fd["$(i)_$(j)_$(j)"] + F0["$(i)_$(j)"] >= 0)
        @constraint(model, F0["$(i)_$(j)"] >= 0)
      end
    end
    @constraint(model, Sd["$(j)"] + S0["$(j)"] >= 0)
    @constraint(model, S0["$(j)"] >= 0)
  end

  # We will use the non-macro version of setObjective, although
  # we could use the macro version if we wanted to.
  # Remember: objectives must be certain.
  @objective(model, :Min, sum([h[i]*x[i] for i=1:n]))

  tic()
  status = solve(model)
  solveTime[instance] = toc()

#  for i=1:n
#    for j=1:n
#      if P[i,j] == 1
#        F0return[i,j] = getvalue(F0["$(i)_$(j)"])
#        for k=1:n
#          Fdreturn[i,j,k] = getvalue(Fd["$(i)_$(j)_$(k)"])
#        end
#      end
#    end
#  end

  for i=1:n
    xreturn[i] = getvalue(x[i])
  end

  aObj = getobjectivevalue(model)
  return sum([h[i]*getvalue(x[i]) for i=1:n]) + getvalue(z)
end

function subsetOf(s, S)
  isSubset = true
  for i in s
    if findfirst(S,i) == 0
      isSubset = false
    end
  end
  return isSubset
end

# Return t, the root of smallest subtree that contains all nodes connected to k via Pfull
# offsprFull: Dict(): node -> [offsprings] based on Pfull
# offspr: Dict(), root -> [offsprings] (based on tree)
function generalRoot(k,offsprFull,offspr)
  for t=k:-1:1
    if subsetOf(get(offsprFull,k,[]), get(offspr,t,[]))
      return t
    end
  end
  println("General root did not find root for ")
  println("offsprFull")
  println(offsprFull)
  println("offspr")
  println(offspr)
  println("k")
  println(k)
  println("get(offsprFull,k,[])")
  println(get(offsprFull,k,[]))
  println("get(offspr,k,[])")
  println(get(offspr,k,[]))
  return -1 # this condition shouldn't be activated
end

# TODO add offsprFull,
# TODO implement generalRoot(k,Pfull,Ptree)
# affineXi, based on affineXitree, now general network
# Main new edit: for each node i, extract the set of demand nodes and it's
# Refer to scanned notes May 16, 2018
# P is Pfull
# offsprFull: Dict(): node -> [offsprings] based on Pfull
# All following based on tree structure:
# Ptree,
# parent: Dict(), kid -> parent
# kids: Dict(), parent -> [kids]
# offspr: Dict(), root -> [offsprings]
# offsprNr: contains root too
function affineXi(n, nS, Γ, Atree, P, Ptree, parent, kids, offsprFull, offspr, offsprNr, h, b, c, dmean, dvar, F0return, Fdreturn, xreturn, solveTime, instance)
  model = RobustModel(solver=GurobiSolver(Method=2, Crossover=0, BarConvTol=0.000001));
  #model = RobustModel();

  # Set up variables
  @variable(model, x[1:n] >= 0); # inventory decisions
  FdIndices = String[]
  F0Indices = String[]
  SdIndices = String[]
  S0Indices = String[]
  for j=(nS+1):n # for every demand node
    push!(S0Indices,"$(j)")
    push!(SdIndices,"$(j)")
    for i=1:n
      if P[i,j] == 1
        push!(F0Indices,"$(i)_$(j)")
        for k=j
          push!(FdIndices,"$(i)_$(j)_$(k)")
        end
      end
    end
  end
  @variable(model, F0[F0Indices]) # offset terms for affine policy flows
  @variable(model, Fd[FdIndices]) # affine coefficients for flows; first dimension is for the uncertainty
  @variable(model, S0[S0Indices]); # offset terms for affine policy flows
  @variable(model, Sd[SdIndices]); # affine coefficients for flows
  @variable(model, z >= 0); # second stage cost for epigraph format

  # Set up uncertainty polytope
  @uncertain(model, ξ[1:n])
  @constraint(model, ξ[1] >= 1);
  @constraint(model, ξ[1] <= 1);
  for i=1:n
    @constraint(model, ξ[i] >= 0);
    @constraint(model, ξ[i] <= 1);
    for j=1:n
      if (Atree[i,j] > 0) # A is a tree
        @constraint(model, ξ[i] >= ξ[j]);
      end
    end
    if (sum(Atree[i,:]) > 0) # this node has children. A is a tree
      @constraint(model, sum([ξ[j]*Atree[i,j] for j=1:n]) <= Γ[i] * ξ[i]); # A is a tree
    end
  end

  # demand satisfaction constraints, deterministic version since only involve \xi_j
  for j=(nS+1):n
    # ξ_j = 1
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? Fd["$(i)_$(j)_$(j)"] + F0["$(i)_$(j)"] : 0 ) for i=1:n]) + Sd["$(j)"] + S0["$(j)"] >= dmean[j] + dvar[j])
    # ξ_j = 0
    @constraint(model, x[j] + sum([(P[i,j] == 1 ? F0["$(i)_$(j)"] : 0 ) for i=1:n]) + S0["$(j)"] >= dmean[j])
  end

  # supply node capacity/conservation constraints:
  #for i=1:n
  #  @constraint(model, x[i] >= sum([(P[i,j] == 1 ? sum([ξ[s]*Fd["$(i)_$(j)_$(s)"] for s=j]) + F0["$(i)_$(j)"] : 0) for j=(nS+1):n]))
  #end

  pIndices = String[]
  qIndices = String[]
  rIndices = String[]
  # robustified supply node capacity constraints. Refer to scanned notes May 16.
  for k=1:nS
    t = generalRoot(k,offsprFull,offspr)
    push!(pIndices,"$(k)_$(t)")
    for i in offsprNr[t]
      push!(rIndices,"$(k)_$(t)_$(i)")
      for j in get(kids,i,[])
        push!(qIndices, "$(k)_$(t)_$(i)_$(j)")
      end
    end
  end
  # dual variables for supply conservation constraints
  @variable(model, p[pIndices]) # 1:nS if demand nodes cannot serve demand nodes.
  @variable(model, q[qIndices] >= 0)
  @variable(model, r[rIndices] >= 0)

  for k=1:nS
    t = generalRoot(k,offsprFull,offspr)
    @constraint(model, p["$(k)_$(t)"] <= x[k] - sum([(j > nS ? F0["$(k)_$(j)"] : 0) for j in get(offsprFull,k,[])]))
    @constraint(model, p["$(k)_$(t)"] - sum([q["$(k)_$(t)_$(t)_$(j)"] for j in get(kids,t,[])]) - Γ[t]*r["$(k)_$(t)_$(t)"] >= 0)
    for i in get(offspr,t,[])
      if i <= nS
        @constraint(model, q["$(k)_$(t)_$(parent[i])_$(i)"] - sum([q["$(k)_$(t)_$(i)_$(j)"] for j in get(kids,i,[])]) + r["$(k)_$(t)_$(parent[i])"] - Γ[i]*r["$(k)_$(t)_$(i)"] >= 0)
      else
        if findfirst(get(offsprFull,k,[]),i) == 0
          @constraint(model, q["$(k)_$(t)_$(parent[i])_$(i)"] + r["$(k)_$(t)_$(parent[i])"] >= 0)
        else
          @constraint(model, q["$(k)_$(t)_$(parent[i])_$(i)"] + r["$(k)_$(t)_$(parent[i])"] >= Fd["$(k)_$(i)_$(i)"])
        end
      end
    end
  end

  # non-negativity of flow and loss variables
  # deterministic version given 0 <= ξ[i] <= 1, and worst case happens at extreme point
  for j=(nS+1):n
    for i = 1:n
      if P[i,j] == 1
        @constraint(model, Fd["$(i)_$(j)_$(j)"] + F0["$(i)_$(j)"] >= 0)
        @constraint(model, F0["$(i)_$(j)"] >= 0)
      end
    end
    @constraint(model, Sd["$(j)"] + S0["$(j)"] >= 0)
    @constraint(model, S0["$(j)"] >= 0)
  end

  #### YES delay penalty ####
  @constraint(model, z >= sum([(ξ[j]*Sd["$(j)"] + S0["$(j)"]) * b[j] for j=(nS+1):n]) +  sum([sum([(P[i,j]==1 ? (ξ[j]*Fd["$(i)_$(j)_$(j)"] + F0["$(i)_$(j)"])*c[i,j] : 0) for j=(nS+1):n]) for i=1:n])) # with delay penalty

  # We will use the non-macro version of setObjective, although
  # we could use the macro version if we wanted to.
  # Remember: objectives must be certain.
  @objective(model, :Min, sum([h[i]*x[i] for i=1:n]) + z)

  tic()
  status = solve(model)
  solveTime[instance] = toc()

#  for i=1:n
#    for j=1:n
#      if P[i,j] == 1
#        F0return[i,j] = getvalue(F0["$(i)_$(j)"])
#        for k=1:n
#          Fdreturn[i,j,k] = getvalue(Fd["$(i)_$(j)_$(k)"])
#        end
#      end
#    end
#  end

  for i=1:n
    xreturn[i] = getvalue(x[i])
  end

  aObj = getobjectivevalue(model)
  return sum([h[i]*getvalue(x[i]) for i=1:n]) + getvalue(z)
end
