%% NLP: min p(x) s.t. g(x)>=0
%% Ashkan Jasour, Research Scientist, MIT 2020
% jasour.mit.edu  rarnop.mit.edu
%%
clc;clear all;close all
%% Parameters
nx=10; % number of design parameters
POW=2; % degree of polynomila
d=2; % relaxation order of SDP: SDP uses 2*d number of the moments.
% Parameters of the SDP solver : ProxSDP
max_iter=10^5;
tol_primal=0.008;
%% NLP  min p(x) s.t. g(x)>=0
x=mpvar('x',[1 nx]);               
p=sum((x-0.5).^2);
g=0.5^2-sum((x-0.5).^2)-(x(1)-0.5)^(2*POW) ;

%% 1: Generate standard SDP min <C,X> s.t. <A,X>=b
display('Generate standard SDP')
[A,C,b]=func_Standard_SDP_Gen(nx,d,p,g);

%% 2: Save the Data for Julia
for i=1:size(A,2); 
    fname = sprintf('A%d.txt', i);
    dlmwrite(fname,full(A{i}),'delimiter',' '); 
end
dlmwrite('b.txt',b,'delimiter',' ')
dlmwrite('C.txt',full(C),'delimiter',' ')
dlmwrite('tol_primal.txt',tol_primal,'delimiter',' ')
dlmwrite('max_iter.txt',max_iter,'delimiter',' ')

%% Solve with CVX
m=size(A,2);  n=size(A{1},1);
cvx_begin
    variable XX( n, n ) symmetric
    dual variables y{m}
sdpsettings('solver','mosek')

    minimize( trace(C*XX) );
       for k = 1 : m;
        trace(A{k}*XX)==b(k+1) : y{k};  % coeff of x^k , b=[x^0,x,...]
    end
    XX == semidefinite(n);
cvx_end

p_0=b(1); 
Optimum_cvx=[p_0-trace(C*XX)];
Sol_cvx_dual=-cell2mat(y(1:nx));
 

%% Solve with GloptiPoly
 
mset('yalmip',true);mset(sdpsettings('solver','mosek'));
mpol('x',1,nx); p = sum((x-0.5).^2); g=0.5^2-sum((x-0.5).^2)-(x(1)-0.5)^(2*POW);

P = msdp(min(p), g>=0, d)
[status,Optimum_Glopti,~,dua1] = msol(P);
y_mom=double(mvec(meas)); Sol_Glopti=y_mom(2:2+nx-1);

%% Compare
[Sol_cvx_dual,Sol_Glopti]
[Optimum_cvx, Optimum_Glopti]

