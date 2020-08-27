%% Risk = probability (g(x)>=0)  where x: random vector
% Eaxample:  Probability ( -x(1)^4+0.5*(x(1)^2-x(2)^2)+0.1 >=0 ) 
% where x1 hase Uniform probability distribution over [0,1]
%       x2 has Beta(1.3,3) probability distribution over [0,1]
% Lecture 10: Probabilistic Nonlinear Safety Verification, rarnop.mit.edu 
%% Ashkan Jasour, Research Scientist, MIT 2020
% jasour.mit.edu  rarnop.mit.edu
%%
clc;clear all;close all
%% Parameters
nx=2; % number of uncertain parameters
d=5; % relaxation order: SDP uses 2*d number of the moments of uncertainties to calculate the Risk.
% Parameters of the SDP solver : ProxSDP
max_iter=10^5;
tol_primal=0.008;
%% Safety Constraint : Risk= Prob(g(x)>=0)  x:random vector
x=mpvar('x',[1 nx]);               
g=-x(1)^4+0.5*(x(1)^2-x(2)^2)+0.1;
%% Moments Information

%moments of Lebesgue Measure over [-1,1]^2 to calculate the integral
u=1;l=-1; yL=[2];for i=1:2*d ;yL(i+1,1)= ( u^(i+1) - l^(i+1) )/(i+1);end 
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(nx,k)]; end; 
yL=prod(yL(vpow+1),2);

%moments of Uniform probability distribution over [0,1] for uncertain variable x1
u=1;l=0;yx1=[1];for i=1:2*d ;yx1(i+1,1)=(1/(u-l))*((u^(i+1) - l^(i+1))/(i+1));end 

%moments of Beta(aB,bB) probability distribution over [0,1] for uncertain variable x2
aB=1.3;bB=3; yx2=[1];for k=1:2*d; yx2=[yx2;(aB+k-1)/(aB+bB+k-1)*yx2(end) ]; end;

% moments of joint distribution of x1 and x2
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(nx,k)]; end; 
yx1x2=yx1(vpow(:,1)+1).*yx2(vpow(:,2)+1);

%% 1: Generate standard SDP min <C,X> s.t. <A,X>=b
display('Generate standard SDP')
[A,C,b,Mind]=func_Standard_SDP_Gen(nx,d,g,yL);

%% 2: Save the Data for Julia
for i=1:size(A,2); 
    fname = sprintf('A%d.txt', i);
    dlmwrite(fname,full(A{i}),'delimiter',' '); 
end
dlmwrite('b.txt',b,'delimiter',' ')
dlmwrite('C.txt',C,'delimiter',' ')
dlmwrite('Mind.txt',Mind,'delimiter',' ')
dlmwrite('yx1x2.txt',yx1x2,'delimiter',' ')
dlmwrite('tol_primal.txt',tol_primal,'delimiter',' ')
dlmwrite('max_iter.txt',max_iter,'delimiter',' ')
