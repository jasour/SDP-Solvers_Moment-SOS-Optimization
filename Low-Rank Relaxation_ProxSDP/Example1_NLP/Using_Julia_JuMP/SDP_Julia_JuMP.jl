###########################################################################
#This code solves the standard SDP given the matrices A,b,C
###########################################################################
# Ashkan Jasour, Research Scientist, MIT 2020
# jasour.mit.edu  rarnop.mit.edu
###########################################################################

println("Calling the packages")

using ProxSDP
using JuMP
using LinearAlgebra
using DelimitedFiles
using DynamicPolynomials
using MosekTools

filepath=@__DIR__;
cd(filepath)

println("Reading the SDP Data")
# read the matrices
C=readdlm("C.txt"); b=readdlm("b.txt");
nc=size(b,1)-1;nx=size(C,1)
A=Array{Float64}(undef, nx,nx,nc);
t_p=readdlm("tol_primal.txt");
m_i =readdlm("max_iter.txt");

#SDP solver ProxSDP
model = Model(with_optimizer(ProxSDP.Optimizer, tol_primal=t_p[1], max_iter=Int(m_i[1]), log_verbose=true))
#model = Model(with_optimizer(Mosek.Optimizer))

@variable(model, X[1:nx, 1:nx], PSD)
@objective(model, Max,  -1*dot(C, X))

println("Constraints of the SDP")
for i in 1:nc
print("--$(i)--")
A[:,:,i]=readdlm("A$i.txt");
@constraint(model, dot(A[:,:,i],X) .== b[i+1])
end

println("Construct the SDP")
# Solve the SDP
JuMP.optimize!(model)

# SDP Solution
Xsol = JuMP.value.(X)

# write the obtained solution to text file
open("Sol.txt", "w") do io; writedlm(io, Xsol, ','); end;

Optimum=[b[1]-tr(C*Xsol)]
