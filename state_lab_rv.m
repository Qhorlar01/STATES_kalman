clear all
close all
clc

m=0;
sig=1;


%%random variables
N1 = 100; % number of realisation

N2 = 4000; % number of realisation

nb1= sqrt (N1); %number of bins
nb2 = sqrt (N2); %number of bins
%% normal distribution
x_norm1 = randn (N1,1);
x_norm2 = randn (N2,1);
pdfGauss = @(x_norm1, m, sig) 1/sqrt(2*pi)/sig*exp(-0.5*((x_norm1-m)/sig).^2);
x = linspace (-4,4,100)
x_norm1mean= mean (x_norm1)
fprintf('1.1. mean of the normal distribution: \n');
disp(x_norm1)
x_norm1std= std (x_norm1)
fprintf('1.1. standard deviation of the normal distribution: \n');
disp(x_norm1)

x_norm2mean= mean (x_norm2)
fprintf('1.2. mean of the normal distribution: \n');
disp(x_norm2)
x_norm1std= std (x_norm1)
fprintf('1.2. standard deviation of the normal distribution: \n');
disp(x_norm2)

figure(1)
subplot(1,2,1)
[n,b] = hist(x_norm1,nb1);
bar (b,n/(b(2)-b(1))/sum(n), 1)
hold on
plot (x,pdfGauss(x,m,sig),"r")

subplot(1,2,2)
xn=[0:1:length(x_norm1)-1]
plot(x_norm1,"x")

figure(2)      
subplot(1,2,1)
[n,b] = hist(x_norm2,nb2);
bar (b,n/(b(2)-b(1))/sum(n), 1)
hold on
plot (x,pdfGauss(x,m,sig),"r")
subplot(1,2,2)
xn=[0:1:length(x_norm2)-1]
plot(x_norm2,"x")
%%Uniform Distribution
x_uni1 = 2*sqrt(3)*(rand(N1,1)-0.5);

x_uni2 = 2*sqrt(3)*(rand(N2,1)-0.5);

xs= [-sqrt(3) -sqrt(3) sqrt(3) sqrt(3)];
ys = [0 1/(2*sqrt(3)) 1/(2*sqrt(3)) 0];
x_uni1mean= mean (x_uni1)
x_uni2mean= mean (x_uni2)

x_uni1std= std (x_uni1)
x_uni2std= std (x_uni2)

figure(3)
subplot(1,2,1)
[n,b] = hist(x_uni1,nb1);
bar (b,n/(b(2)-b(1))/sum(n), 1)
hold on 
plot (xs,ys,"r")
subplot(1,2,2)
xn=[0:1:length(x_uni1)-1]
plot(x_uni1,"x")

figure(4)
subplot(1,2,1)
[n,b] = hist(x_uni2,nb2);
bar (b,n/(b(2)-b(1))/sum(n), 1)
hold on
plot (xs,ys,"r")
subplot(1,2,2)
xn=[0:1:length(x_uni2)-1]
plot(x_uni1,"x")
%gaussian mixture
m = 0.95;
sig= sqrt (1-m^2);
x_gauss1 = randn (N1,1)*sqrt(1-m*m)+m;
k1 = find(rand(N1,1)>0.5);
x_gauss1(k1) = x_gauss1(k1)-2*m;

pdfGauss = @(x_gauss1, m, sig) 1/sqrt(2*pi)/sig*exp(-0.5*((x_gauss1-m)/sig).^2);

x_gauss1mean= mean (x_gauss1)

x_gauss1std= std (x_gauss1)


%plot
figure(5)
subplot(1,2,1)
[n,b] = hist(x_gauss1,nb1);
bar (b,n/(b(2)-b(1))/sum(n), 1)
hold on
plot (x, 0.5*pdfGauss(x, -m, sig)+0.5*pdfGauss(x, m, sig))
subplot(1,2,2)
xn=[0:1:length(x_gauss1)-1]
plot(x_gauss1,"x")

x_gauss2 = randn (N2,1)*sqrt(1-m*m)+m;
k2 = find(rand(N2,1)>0.5);
x_gauss2(k2) = x_gauss2(k2)-2*m;

pdfGauss = @(x_gauss2, m, sig) 1/sqrt(2*pi)/sig*exp(-0.5*((x_gauss2-m)/sig).^2);

x_gauss2mean= mean (x_gauss2)
x_norm2std= std (x_norm2)


%plot
figure(6)
subplot(1,2,1)
[n,b] = hist(x_gauss2,nb2);
bar (b,n/(b(2)-b(1))/sum(n), 1);
hold on
plot (x, 0.5*pdfGauss(x, -m, sig)+0.5*pdfGauss(x, m, sig))

subplot(1,2,2)
xn=[0:1:length(x_gauss2)-1]
plot(x_gauss2,"x")

