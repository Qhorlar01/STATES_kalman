clear all
close all
clc
%% JOINT DISTRIBUTION
N=200;
sig1=2;
sig2=5;

var1=(sig1)^2;
var2=(sig2)^2;

rho=0.9;
m=[0;0];
po=0.91; %degree of confidence 

Cxy=rho*sig1*sig2;
C=[var1 Cxy; Cxy var2];
% To Generate N = 200 realizations of this random variable.
x= chol(C,"lower")*randn(length(C),N)+m*ones(1,N);
mean_x1=[mean(x(1,:));mean(x(2,:))]
cov_x=cov(x(1,:),x(2,:));
% To plot confidence ellipse
t=linspace(0,2*pi,200);
X=sqrt(-2*log(1-po))*chol(C,'lower')*[cos(t);sin(t)]+m*ones(1,length(t));



plot(X(1,:),X(2,:))
hold on
plot (x(1,:),x(2,:),'x')


