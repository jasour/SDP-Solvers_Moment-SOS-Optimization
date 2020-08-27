using MomentOpt
using DynamicPolynomials
using MosekTools
using ProxSDP


# dimension of NLP
n=2
# degree of the NLP
POW=2
# relaxation order: it uses (2*d- maxdegree(Obj,Con)) number of moments
d=3


# NLP, optimal sol x*=0.5*In
@polyvar x[1:n]
Obj=sum((x-0.5*ones(n)).^2)
Con=0.5^2-sum((x-0.5*ones(n)).^2)-(x[1]-0.5)^(POW);

# SDP solver
#model = GMPModel(with_optimizer(Mosek.Optimizer))
model = GMPModel(with_optimizer(ProxSDP.Optimizer,tol_primal=0.001, max_iter=10^5 ,log_verbose=true))

# define probability measure mu and its support
@variable model mu Meas(x, support = @set(Con >= 0))
@constraint model Mom(1, mu) == 1

# objective function of  moment SDP
@objective model Min Mom(Obj, mu)

# relaxation order of Opt
set_approximation_degree(model, d)

# solve moment SDP
optimize!(model)

# Obtained results in moments
mu_sol = measure(mu)
v = objective_value(model)
y = moment_value.(moments(mu_sol))
L = monomials(mu_sol)

# Extracted solution of NLP
x_sol = atomic(mu);
println("Objective value: $(v)")
println("Extracted solution of NLP: $(x_sol)")

#if x_sol isa Nothing
#    println("Using MultivariateSeries for extraction")
#    using MultivariateSeries
#    _, x_sol = ms_decompose(sum(series(yi,l) for (l,yi) in zip(L, y)))
#end
#x_sol