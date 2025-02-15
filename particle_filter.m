%Particle filter

P1      =  0.1;
Q       =  10;
R       =  1;
g       =  @(v, nu) 1/sqrt(2*pi*nu)*exp(-v.*v/2/nu) ;   % zero mean variance nu gaussian PDF


N       = 100;
x_0     = 0;
m       = 0;
w       = sqrt(R)*randn(N,1) + m*ones(N,1);
v       = sqrt(Q)*randn(N,1) + m*ones(N,1);
t       = 1:N;
%%

[x,y]   =  simulate(N,x_0,w,v);

 f10      = figure(10);set(f10,'name','state trajectory');
 %subplot(211)
% plot(t,x,'-r')
 %leg1    = legend('$x_n$');
 %set(leg1,'Interpreter','latex');
 %grid on
 %title('State Trajectory of $x$','Interpreter','latex')

  plot(t,y,'-b')
 leg1    = legend('$y_n$');
 set(leg1,'Interpreter','latex');
 grid on
 title('State Trajectory of output $y$','Interpreter','latex')

xekf = ekf(y,P1,Q,R,N,x_0);

xCorrectedUKF = UKF(R, Q,y,N);

xEst    = Particle(particleFilter( @f ,@LikelihoodFcn), 1000000,0,Q,y,N);

f11      = figure(11);set(f11,'name','Extended Kalman Filter');
plot(t,x,'-k',t,xEst(1:N),'-r', t,xekf(1:N),'-b', t,xCorrectedUKF(1:N),'-g','Linewidth',1.1)
leg1    = legend('$x_n$','$x_{n\{Particle\}}$','$x_{n\{EKF\}}$','$x_{n\{UKF\}}$');
set(leg1,'Interpreter','latex');
grid on
title('Comaprison of State  Trajectory $x$, with the $x_{n\{Particle\}}$, $x_{n\{UKF\}}$, and $x_{n\{EKF\}}$  ','Interpreter','latex')

f12      = figure(12);set(f12,'name','Extended Kalman Filter');
plot(t,y,'-k',t,xEst(1:N),'-r', t,xekf(1:N),'-b', t,xCorrectedUKF(1:N),'-g','Linewidth',1.1)
leg1    = legend('$y_n$','$y_{n\{Particle\}}$','$y_{n\{EKF\}}$','$y_{n\{UKF\}}$');
set(leg1,'Interpreter','latex');
grid on
title('Comaprison of State  Trajectory $y$, with the $y_{n\{Particle\}}$, $y_{n\{UKF\}}$, and $y_{n\{EKF\}}$  ','Interpreter','latex')

ll = 0;









%plot(x,'k'), hold on, plot(xEst,'r'), plot(xCorrectedUKF,'b'),  plot(xekf,'g'), hold off


function xe = ekf(y,P1,Q,R,N,x_0)

     xe  = zeros(N,1);
     xe(1) = x_0;
     P = P1;
     for n = 1:N
         H = dh(xe(n),n);

         y_pred = h(xe(n),n);
         Cxy = P*H;
         Cyy = H*P*H + R;
         
         y_observed = y(n);
         
         xe(n) = xe(n) + (y_observed - y_pred)*Cxy/Cyy;
         P = P - Cxy*Cxy/Cyy;
         
         F = df(xe(n),n);
         
         xe(n+1) = f(xe(n),n);
         
         P = F*P*F + Q;
     end

end


function xCorrectedUKF = UKF(R, Q,y,N)

    initialStateGuess = 0; 

    ukf = unscentedKalmanFilter(...
         @(x, n)(0.5*x + 25*x./(1+x.*x)+8*cos(1.2*n)),... % State transition function
         @(x)(x.*x/20),... % Measurement function
        initialStateGuess,...
        'HasAdditiveMeasurementNoise',true,...
        'HasAdditiveProcessNoise',true, 'alpha', 0.55);
    ukf.MeasurementNoise = R;
    ukf.ProcessNoise = Q;

    for k=1:N
        % Let k denote the current time.
        %
        % Residuals (or innovations): Measured output - Predicted output
        e(k) = y(k) - f(ukf.State,k); % ukf.State is x[k|k-1] at this point
        % Incorporate the measurements at time k into the state estimates by
        % using the "correct" command. This updates the State and StateCovariance
        % properties of the filter to contain x[k|k] and P[k|k]. These values
        % are also produced as the output of the "correct" command.
        [xCorrectedUKF(k), PCorrected(k,:,:)] = correct(ukf,y(k));
        % Predict the states at next time step, k+1. This updates the State and
        % StateCovariance properties of the filter to contain x[k+1|k] and
        % P[k+1|k]. These will be utilized by the filter at the next time step.
        predict(ukf,k);
    end

end
%%
function xEst = Particle(myPF,particles, initial, Q,y,N)
    initialize(myPF,particles,initial,Q);
    myPF.StateEstimationMethod = 'mean';
    myPF.ResamplingMethod = 'systematic';

    xEst = zeros(N,1);
    for k=1:N
        xEst(k) = correct(myPF,y(k));
        predict(myPF,k);
    end
end


%%

function likelihood = LikelihoodFcn(predictedParticles,measurement)
% vdpMeasurementLikelihoodFcn Example measurement likelihood function for 
%                             particle filter
%
% The measurement is the first state.
%
% likelihood = vdpMeasurementLikelihoodFcn(predictedParticles, measurement)
%
% Inputs:
%    predictedParticles - NumberOfStates-by-NumberOfParticles matrix that
%                         holds the predicted particles
%
% Outputs:
%    likelihood - A vector with NumberOfParticles elements whose n-th
%                 element is the likelihood of the n-th particle
%
% See also extendedKalmanFilter, unscentedKalmanFilter

%   Copyright 2017 The MathWorks, Inc.

%#codegen

% The tag %#codegen must be included if you wish to generate code with 
% MATLAB Coder.

numberOfMeasurements = 1; % Expected number of measurements

% Validate the measurement
validateattributes(measurement, {'double'}, {'vector', 'numel', numberOfMeasurements}, ...
    'vdpMeasurementLikelihoodFcn', 'measurement');

% Assume that measurements are subject to Gaussian distributed noise with
% variance 0.016
% Specify noise as covariance matrix
measurementNoise = 1 * eye(numberOfMeasurements);
  
% The measurement contains the first state variable. Get the first state of
% all particles
predictedMeasurement = h(predictedParticles(1,:),0);

% Calculate error between predicted and actual measurement
measurementError = bsxfun(@minus, predictedMeasurement, measurement(:)');

% Use measurement noise and take inner product
measurementErrorProd = dot(measurementError, measurementNoise \ measurementError, 1);

% Convert error norms into likelihood measure. 
% Evaluate the PDF of the multivariate normal distribution. A measurement
% error of 0 results in the highest possible likelihood.
likelihood = 1/sqrt((2*pi).^numberOfMeasurements * det(measurementNoise)) * exp(-0.5 * measurementErrorProd);
end

function [x,y] =  simulate(N,x_0,w,v)

    x = zeros(N,1);
    y = zeros(N,1);
    
    x(1) = x_0;
    y(1) = h(x(1),1) + w(1);

    for n = 2:N
        x(n) = f(x(n-1),n-1) + v(n);
        y(n) = h(x(n),n) + w(n);  
    end

end

function y    = f(x, n) 

y = 0.5*x + 25*x./(1+x.*x)+8*cos(1.2*n);

end

function y    = h(x, n) 

y = x.*x/20;
end

function y    = df(x, n) 

y = 0.5 + 25*(1-x.*x)./(1+x.*x)./(1+x.*x);
end

function y    = dh(x, n) 

y = x/10;
end