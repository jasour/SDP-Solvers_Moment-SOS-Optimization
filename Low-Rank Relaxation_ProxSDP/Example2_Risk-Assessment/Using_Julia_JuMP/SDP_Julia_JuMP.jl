###########################################################################
#This code solves the standard SDP given the matrices A,b,C
###########################################################################
# Risk = probability (g(x)>=0) where x is a random vector and g(x)>=0 is a nonlinear safety constraint.
#############################################################################
# measure-moment SDP formulation:
# max_{mu} vol(mu)                     : obj
#   s.t mu supported on g(x)>=0        : support constraint
#       mu <= given Lebesgue Measure   : measure constraint
# where mu is unknown measure defined on safety set g(x)>=0.
# Risk = Expected value of {Polynomial Indicator function of the set {x:g(x)>=0} } with respect to the probability distribution of random vector x.
# Solution of the dual SDP is the coefficients of the polynomial indicator function.
# Lecture 10: Probabilistic Nonlinear Safety Verification, rarnop.mit.edu 
# Ashkan Jasour, Research Scientist, MIT 2020
# jasour.mit.edu  rarnop.mit.edu
###########################################################################

println("Calling the packages")

using ProxSDP
using JuMP
using LinearAlgebra
using DelimitedFiles
using DynamicPolynomials

filepath=@__DIR__;
cd(filepath)

println("Reading the SDP Data")
# read the matrices
C=readdlm("C.txt"); b=readdlm("b.txt");
nc=size(b,1);nx=size(C,1)
A=Array{Float64}(undef, nx,nx,nc);
Mind=readdlm("Mind.txt",Int);
yx1x2=readdlm("yx1x2.txt");
t_p=readdlm("tol_primal.txt");
m_i =readdlm("max_iter.txt");

#SDP solver ProxSDP
model = Model(with_optimizer(ProxSDP.Optimizer, tol_primal=t_p[1], max_iter=Int(m_i[1]), log_verbose=true))

@variable(model, X[1:nx, 1:nx], PSD)
@objective(model, Max,  -1*dot(C, X))

println("Constraints of the SDP")
for i in 1:nc
print("--$(i)--")
A[:,:,i]=readdlm("A$i.txt");
@constraint(model, dot(A[:,:,i],X) .== b[i])
end

println("Construct the SDP")
# Solve the SDP
JuMP.optimize!(model)

# SDP Solution
Xsol = JuMP.value.(X)

# write the obtained solution to text file
open("Sol.txt", "w") do io; writedlm(io, Xsol, ','); end;

# Calculated Risk
Risk=tr(yx1x2[Mind]*Xsol[1:size(Mind,1),1:size(Mind,1)]);
println("Risk: $(Risk)")

