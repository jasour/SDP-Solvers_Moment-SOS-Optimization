clc;
X = dlmread('Sol.txt');
% Rank Test
S=svd(X); ind=find(S>5);size(ind);

% Polynomial Indicator function
z=sym('x',[1 nx]); xv=prod(z.^vpow,2); xc=xv(Mind);
P=trace(xc*X(1:size(Mind,1),1:size(Mind,2))); % polynomial indicator fuction

% plot polynomial indicator fuction
[x1,x2] = meshgrid(-0.95:0.01:0.95,-0.95:0.01:0.95); P=eval(P);
close all;surfc(x1,x2,P,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7);
camlight; lighting gouraud; hold on;grid on;set(gca,'fontsize',25)
title('polynomial indicator function');

