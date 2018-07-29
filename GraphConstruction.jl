module GraphConstruction
using StatsBase
export undirectedGraph
function undirectedGraph(Nvertex)
  VertexArray = Int64[]
  [push!(VertexArray,i) for i in 1: Nvertex]

  MDVertex = rand(1:0.05*Nvertex, Nvertex)
  EdgeMatrix = zeros(Int64, Nvertex, Nvertex)

  i = 1
  while length(VertexArray)>1
    NewEdge = Int64[]
    j = VertexArray[rand(1:length(VertexArray))]
    println(j)
    if j!=i
      EdgeMatrix[i,j]=1; EdgeMatrix[j,i]=1
      NewEdge = sample(1:Nvertex, rand(Int(round(MDVertex[i])):Int(round(0.05*Nvertex))), replace=false)
      [EdgeMatrix[i,NewEdge[k]]=1 for k=1:length(NewEdge)]
      [EdgeMatrix[NewEdge[k],i]=1 for k=1:length(NewEdge)]
      deleteat!(VertexArray,findin(VertexArray,i))
      println(VertexArray)
    end
    i=j;
  end
  return(EdgeMatrix)
end
end
