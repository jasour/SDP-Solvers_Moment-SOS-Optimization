#####################################################
# https://github.com/wangjie212/TSSOS
#####################################################
using TSSOS
using DynamicPolynomials
using SparseArrays
using MultivariatePolynomials


# dimension of NLP
n=20
# order of NLP
POW=20
# relaxation order (order>=POW)
order=POW

@polyvar x[1:n]
f=sum((x-0.5*ones(n)).^2)
pop=[f,0.5^2-sum((x-0.5*ones(n)).^2)-(x[1]-0.5)^(2*POW)]

opt,sol,data=cs_tssos_first(pop,x,order,numeq=0,TS="MD",solution=true)
#opt,sol,data=cs_tssos_higher!(data,TS="MD",solution=true)

println("optimum $(opt)")
println("solution $(sol)")

##################################################
#Options:
#nb: specify the first nb variables to be binary variables (satisfying xi^2=1)
#CS (correlative sparsity): "MD" or "MF" (an approximately minimum chordal extension), "NC" (no chordal extension)
#TS: "block" (using the maximal chordal extension), "MD" or "MF" (using an approximately minimum chordal extension), false (without term sparsity)
#order: d (the relaxation order), "multi" (applying the lowest relaxation orders for each variable clique)
#extra_sos: true (adding a first-order moment matrix for each variable clique), false
#solution: true (extract a solution), false (don't extract a solution)