function [A,C,b]=func_Standard_SDP_Gen(nx,d,p,g)

%% NLP
% min_{x} p(x)    :obj
% s.t. g(x)>=0    :con

% Formulation in terms of nonnegative polynomials:
% max_{gamma} gamma                  :obj
% s.t. p(x)-gamma >=0  on g(x)>=0    :con

% SOS Optimization:
% max_{gamma,sigma_0(x),sigma_1(x)}  gamma        :obj
% s.t. P(x)-gamma = sigma_0(x)+sigma_1(x)g(x)     :cons(1)  ----> Coeffs of right and left habd sides match
%      sigma_0(x)>=0 , sigma_1(x)>=0              :cons(2)  ----> Gram matrix of sigma_0 , sigma_1 is PSD

% SDP Optimization:
% max_{Q,Q0,Q1} gamma         :obj
% s.t. <As0_i,Q0>+<As1g_i,Q1> =p_i       for i>0  : cons(1)
%      <As0_0,Q0>+<As1g_0,Q1> =p_0-gamma for i=0       
%       Q0 >=0, Q1>=0                       :cons(2)
% where 
% Q0:Gram matrix of sigma_0(x), and Q1: Gram matrix of sigma_1(x)

% SDP Optimization: From con(1):  gamma =p_0 - <As0_0,Q0> - <As1g_0,Q1> 
% min_{Q,Q0,Q1} <As0_0,Q0>+<As1g_0,Q1>          :obj
% s.t. <As0_i,Q0>+<As1g_i,Q1> =p_i       for i>0  : cons(1)    
%       Q0 >=0, Q1>=0                       :cons(2)
% in obj, p_0 is constant;Hence, is removed. Also, max is replaced with min.
%Lecture 5: Duality of SOS and Moment based Semidefinite Programs(SDPs), Appendix II
% https://rarnop.mit.edu/Lectures-Codes
%%
clc;
nvar=nx; 
Ny= round(factorial(nvar+2*d)/(factorial(2*d)*factorial(nvar))); % Number of coefficients of P(x): Poly of order 2*d in x  

%% 1: sigma_0(x)=Sum c_i*x^i = x'Q0x,  Q0: Gram matrix
% coeffs of sigma_0(x): <As0_i,Q0> where Q0 is the Gram matrix of sigma_0(x) polynomial
[As0]=func_standard_moment(nvar,d);

%% 2: sigma_1(x)*g(x) 
% coeffs of sigma_1(x)*g(x): <As1g_i,Q1> where Q1 is the Gram matrix of sigma_1(x) polynomial
[As1g]=func_standard_localizing(g,nvar,d);

%% Generate A, C, b in Standard SDP min <C,X> s.t. <A,X>=b
% where X=blkdiag{Q0,Q1} 

for i=1:Ny-1; A{i}=   blkdiag(As0{i+1},As1g{i+1}); end

C=   blkdiag(As0{1},As1g{1});

m=size(As0,2); 
Cpc=p.coef; Cc=0+Cpc(:,:); % Cc1 : coefficent, 
Dpc=p.deg;  Dc=0+Dpc(:,:); % Bc1 : degrees
b=zeros(m,1);
for i=1:size(Dc,1); b(glex2num(Dc(i,:)))=Cc(i) ; end


end

 
