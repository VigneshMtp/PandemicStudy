# Pandemic Disease Integer Programming
workspace()
include("GraphConstruction.jl")
using JuMP
using Gurobi
using .GraphConstruction

# Read Basic Data for pandemic study
InputData =open("Inputfile.txt")
lines =readdlm(InputData)
NI =lines[1]; # No of Indiviuals
TH =lines[2]; # Time Horizon
tw =lines[3]; # Time period in transistion state
ta =lines[4]; # Time period in sick state (alive)
td =lines[5]; # Time period in sick state (dead)
alpha = lines[11]; # Percentage of sick nodes selected for isolation

# Output file declaration
RiskDegreeOutFile = open("RiskDegree.csv","w")
PandemicOutFile = open("PandemicResults.csv", "w")
[print(PandemicOutFile, "S.No\t",",", "EV", ",", "Inf 1 ; 2", ",", "Time", ",", "Explored Nodes",
       ",", "Gap in %", ",", "Objective value", ",", "Infected", ",", "Death", ",", "Isolated\n")]
close(PandemicOutFile)

# Indiviual Interaction and Risk graph construction
EdgeMat = GraphConstruction.undirectedGraph(NI)
WeightMatrix = rand(1:10, NI, NI)
for j = 1:NI
  [print(RiskDegreeOutFile, EdgeMat[i,j]*WeightMatrix[i,j],",") for i in 1:NI-1]
  [print(RiskDegreeOutFile, EdgeMat[NI,j]*WeightMatrix[NI,j],"\n")]
end
close(RiskDegreeOutFile)
w = readcsv("RiskDegree.csv")
d = readcsv("ContDegree.csv");d1=d[1] # contagious degree

# EV = 0 No Enhancement
# EV = 1 All variable fixing
# EV = 2 All variable fixing and moderated big-M coefficients
# EV = 3 All variable fixing, moderated big-M coefficients, and symmetry breaking

for EV in 0:3
 for IIF in 1:30   # No of instances
   NIF1 = reshape([14	44	21	24	12	50	31	56	28	58	31	13	16	48	8
                    3	57	55	44	39	51	52	39	41	12	17	58	58	38	40], 1,30) # Infected indiviual node
   NIF2 = reshape([25	59	30	28	29	39	48	24	22	25	18	23	1	26	21
                    7	18	45	41	32	20	45	58	56	1	25	15	54	15	13], 1,30) # Infected indiviual node
PILP=Model(solver = GurobiSolver(MIPGap = 1e-4, Threads =1, TimeLimit = 1800))
tic()

InitialData =zeros(Int64, NI)
InitialData[NIF1[IIF]] = Int(maximum(w)*maximum(d)); # maximum degree of risk
InitialData[NIF2[IIF]] = Int(maximum(w)*maximum(d));

N1= Int64[];append!(N1,find(x ->(x!=0), InitialData))
println(N1)

# Maximum possible amount of risk in the entire network for person i at time j
u = maximum(d)
Zmax =zeros(NI, TH)
[Zmax[i,j] = Int(sum(w[:,i]*u) - w[i,i]*u) for i=1:NI, j in 1:TH]
Zmax1 = maximum(Zmax)

# S --> set of potential indiviuals to be infected at time period j
S = Array(Vector{Int64}, TH)
S[1] = N1 ; [S[j]=[] for j in 2:TH]
for j in tw+2 : TH
  for k in 1:max(ta,td)
    if j-tw-k>=1
      [issubset(ii,S[j-tw-k]) && !issubset(i,S[j]) && w[ii,i]>0 ? push!(S[j],i) :"" for i in 1:NI, ii in 1:NI]
    end
  end
end

# Big M Coefficients ; set of potential indiviuals to be sick at time period j
SCap = Array(Vector{Int64}, TH)
[SCap[j]=[] for j in 1:TH]
for j in 1:TH
  for k in 1:max(ta,td)
    if j-tw-k>=1
      [issubset(i,S[j-tw-k]) && !issubset(i,SCap[j]) ? push!(SCap[j],i) :""  for i in 1:NI]
    end
  end
end

ZmaxCap = zeros(NI,TH) # Big-M coefficients
SCap1 = Array(Vector{Int64}, TH)
for j in 1:TH
  for i in 1:NI
    SCap1[j] = setdiff(SCap[j], i)
    length(SCap1[j])>0 ? ZmaxCap[i,j]= sum(w[SCap1[j][ii],i]*u for ii in 1:length(SCap1[j])) :""
  end
end
ZmaxCap1 = maximum(ZmaxCap) # Big-M coefficients

beta2 = 0; temp=0; Ztemp =0
SCap1 = Array(Vector{Int64}, TH)
for j in 1:TH
  for i in 1:NI
    SCap1[j] = setdiff(SCap[j], i)
    length(SCap1[j])>0 ? temp = max(sum(ceil(w[SCap1[j][ii],i]/10) for ii in 1:length(SCap1[j])),1) : temp=1
    Ztemp = Ztemp +ZmaxCap[i,j]/temp
  end
end
beta2 =floor(Ztemp/(NI*TH))
beta3 = floor(sum(sum(ZmaxCap,2),1)/(NI*TH))
beta = Int64[0  beta2  beta3]  # Represents amount of virus that has infected an indiviual

A=zeros(NI,3)
[InitialData[i]>=beta[1] && InitialData[i]<=beta[2]-1 ? A[i,1] =1 :A[i,1]=0 for i in 1:NI]
[InitialData[i]>=beta[2] && InitialData[i]<=beta[3]-1 ? A[i,2] =1 :A[i,2]=0 for i in 1:NI]
[InitialData[i]>=beta[3] && InitialData[i]<=Zmax1 ? A[i,3] =1 :A[i,3]=0 for i in 1:NI]

if EV==2 || EV==3
  tempZmax = Zmax
  Zmax = ZmaxCap
  Zmax1 = ZmaxCap1
  [InitialData[i,1]>0 ? Zmax[i,1]=tempZmax[i,1] : "" for i in 1:NI]
end

@defVar(PILP, x[1:NI, 1:TH], Bin)
@defVar(PILP, y[1:NI, 1:TH, 1:3], Bin)
@defVar(PILP, yb[1:NI, 1:TH, 1:max(ta,td), 1:3], Bin)
@defVar(PILP, z[1:NI, 1:TH], Int)
@defVar(PILP, zb[1:NI, 1:TH], Int)

@setObjective(PILP, Min, sum{y[i,j,2]+25*y[i,j,3], i=1:NI, j=1:TH})

@addConstraint(PILP, MaxSelect[j=1:TH], sum{x[i,j], i=1:NI} <= 0 + alpha*sum{sum{y[i,jj,2], jj= max(j-tw-ta,1):j-tw-1}
+ sum{y[i,jj,3], jj= max(j-tw-td,1):j-tw-1}, i = 1:NI})

@addConstraint(PILP, SelectIsolate[i=1:NI], sum{x[i,j], j=1:TH} <= 1)
@addConstraint(PILP, Category[i=1:NI, j=1:TH], sum{y[i,j,l], l=1:3} == 1)
@addConstraint(PILP, Category2[i=1:NI], sum{y[i,j,l], l=2:3, j=1:TH} <=1)

@addConstraint(PILP, SelectIsolate2[i=1:NI, j=1:TH], x[i,j] <= 0 + sum{y[i,jj,2], jj= max(j-tw-ta,1):j-tw-1}
+ sum{y[i,jj,3], jj= max(j-tw-td,1):j-tw-1})

@addConstraint(PILP, Category21[i=1:NI, k=1:ta,j=tw+k+1:TH,  l=2], yb[i,j,k,l]<=y[i,j-tw-k,l])
@addConstraint(PILP, Category22[i=1:NI,  k=1:ta, j=tw+k+1:TH, l=2], yb[i,j,k,l]<=1-sum{x[i,jj], jj=1:j})
@addConstraint(PILP, Category23[i=1:NI, k=1:ta, j=tw+k+1:TH,  l=2], yb[i,j,k,l]>=y[i,j-tw-k,l] - sum{x[i,jj], jj=1:j})

@addConstraint(PILP, Category31[i=1:NI, k=1:td, j=tw+k+1:TH, l=3], yb[i,j,k,l]<=y[i,j-tw-k,l])
@addConstraint(PILP, Category32[i=1:NI, k=1:td, j=tw+k+1:TH,  l=3], yb[i,j,k,l]<=1-sum{x[i,jj], jj=1:j})
@addConstraint(PILP, Category33[i=1:NI, k=1:td, j=tw+k+1:TH,  l=3], yb[i,j,k,l]>=y[i,j-tw-k,l] - sum{x[i,jj], jj=1:j})

@addConstraint(PILP, AuxConstraint1[i=1:NI, j=2:TH], z[i,j] == sum{(sum{w[ii,i]*d[k]*yb[ii,j,k,2], k=1:min(ta,j-tw-1)}+sum{w[ii,i]*d[k]*yb[ii,j,k,3], k=1:min(td,j-tw-1)}), ii=1:i-1}
+ sum{(sum{w[ii,i]*d[k]*yb[ii,j,k,2], k=1:min(ta,j-tw-1)}+sum{w[ii,i]*d[k]*yb[ii,j,k,3], k=1:min(td,j-tw-1)}), ii=i+1:NI})

@addConstraint(PILP, AuxConstraint2[i=1:NI, j=1:TH, l=1:3], y[i,j,l]*beta[l]<=zb[i,j])
@addConstraint(PILP, AuxConstraint3[i=1:NI, j=1:TH, l=1:2], zb[i,j]<=(beta[l+1]-1-Zmax1)*y[i,j,l]+Zmax1)
@addConstraint(PILP, AuxConstraint4[i=1:NI, j=1:TH], zb[i,j]<=z[i,j])
@addConstraint(PILP, AuxConstraint5[i=1:NI, j=1:TH], zb[i,j]<=Zmax[i,j]*(1-sum{y[i,jj,2], jj=1:j-1}-sum{y[i,jj,3], jj=1:j-1}))
@addConstraint(PILP, AuxConstraint6[i=1:NI, j=1:TH], zb[i,j]>=z[i,j]-Zmax1*(sum{y[i,jj,2], jj=1:j-1}+sum{y[i,jj,3], jj=1:j-1}))
@addConstraint(PILP, InitialDConstraint[i=1:NI], z[i,1]==InitialData[i])

iT = Array(Vector{Int64}, TH)
[iT[j]=[] for j in 1:TH]
[iT[j] = setdiff(1:NI,S[j]) for j in 1:TH]

if EV == 1 || EV==2 || EV==3
  @addConstraint(PILP, InitialVConstraint[i=1:NI, j=1:tw+1], x[i,j]==0)
  @addConstraint(PILP, InitialV1Constraint[i=1:NI, j=2:tw+1], y[i,j,1]==1)
  @addConstraint(PILP, InitialV2Constraint[i=1:NI, k=1:3], y[i,1,k]==A[i,k])
  @addConstraint(PILP, InitialV3Constraint[j=tw+2:TH, i = 1:length(iT[j])], z[iT[j][i],j]==0)
  @addConstraint(PILP, InitialV4Constraint[j=tw+2:TH, i = 1:length(iT[j])], y[iT[j][i],j,1]==1)
end

TT = []
[length(S[j-1])>0 ? push!(TT,j) :"" for j in tw+3:TH]

if EV == 3
  @addConstraint(PILP, ValidIneq[i=1:NI, j=1:length(TT)], x[i,TT[j]]-y[i,TT[j]-tw-1,2]-y[i,TT[j]-tw-1,3] - (1-(alpha*sum{sum{y[i,jj,2], jj= max(TT[j]-tw-ta-1,1):TT[j]-tw-1-1}
  + sum{y[i,jj,3], jj= max(TT[j]-tw-td-1,1):TT[j]-tw-1-1}, i = 1:NI} - sum{x[i,TT[j]-1], i=1:NI})/length(S[TT[j]-1]))<=x[i,TT[j]-1])
end
solve(PILP)
println("\n")
println("Total Infected\t", sum(getValue(y[:,:,2])))
println("Total Death\t", sum(getValue(y[:,:,3])))
println("Total Isolated\t", sum(getValue(x[:,:])))
println()
println(getObjectiveValue(PILP))
println(getnodecount(PILP))
println(getsolvetime(PILP))
println(getobjbound(PILP))
println(abs((getobjbound(PILP)-getObjectiveValue(PILP))/getObjectiveValue(PILP))*100)
for j in 1:TH
  print(j,"\t")
  [getValue(x[i,j])>0 ? print(i," "):"" for i in 1:NI];print("\t\t")
  [getValue(y[i,j,2])>0 ? print(i," "):"" for i in 1:NI];print("\t\t")
  [getValue(y[i,j,3])>0 ? print(i," "):"" for i in 1:NI];print("\t\t")
  println()
end


PandemicOutFile= open("PandemicResults.csv", "a")
print(PandemicOutFile, IIF ,",", EV, ",", NIF1[IIF]," ; ",NIF2[IIF], ",", getsolvetime(PILP), ",",getnodecount(PILP), ",",abs((getobjbound(PILP)-getObjectiveValue(PILP))/getObjectiveValue(PILP))*100, ",",getObjectiveValue(PILP), ",")
print(PandemicOutFile,sum(getValue(y[:,:,2])), ",",sum(getValue(y[:,:,3])),",", sum(getValue(x[:,:])))
print(PandemicOutFile,"\n")
close(PandemicOutFile)
end
PandemicOutFile= open("PandemicResults.csv", "a")
print(PandemicOutFile,"\n\n")
close(PandemicOutFile)
end
