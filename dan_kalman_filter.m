%Parameters for simulation
%====================================
Duration    = 0.5;
Ts          = 0.001;
delta       = 0.1;
amplitude   = 0.1;
T           = 20e-3;
G           = 50;
L           = 512;
x1          = [0.1,0];
q           = 1e2;

T_filter          = 25e-3;

%=====================================

%Test Suite
%====================================
%[sysd,A,B,C]     = getdigitizedsystem(G,T,Ts);
%[yy,tOut ,xx]    = lsim(sysd,u,t,x1) ;
 %y = round(y*L/2/pi)*2*pi/L;


[u,t]             = inputvoltage(Duration,amplitude,delta,Ts);

 f10      = figure(10);set(f10,'name','Input Voltage');
 plot(t,u,'-k')
 leg1    = legend('$u(t)$');
 set(leg1,'Interpreter','latex');
 grid on
 title('Input Voltage u(t) with $\Delta$ = 100ms and $V_{peak-to-peak}$ = 0.1V','Interpreter','latex')


[yt, y,x]             = simulate(u,G,T,Ts,L,[0,0]);
f11      = figure(11);set(f11,'name','Simulation Result, ');
 subplot(111), plot(t,y,'-b',t,yt,'-r','LineWidth',1.5)
 leg1    = legend('$y<\theta(t)>$','$y_n<\theta_n>$');
 set(leg1,'Interpreter','latex');
 grid on
 title('Output $y <\theta(t)> $, and quantized $y_n <\theta_n>$','Interpreter','latex')


 f12      = figure(12);set(f12,'name','Simulation Result, ');
 
 subplot(211), plot(t,x(1,:),'-b','LineWidth',1.2)
 leg1    = legend('$x_1 = \theta(t)$');
 set(leg1,'Interpreter','latex');
 grid on
 title('State Variable $x_1 = \theta(t)$','Interpreter','latex')
 
 subplot(212), plot(t,x(2,:),'-g','LineWidth',1.2)
 leg1    = legend('$x_2 = \Omega(t)$');
 set(leg1,'Interpreter','latex');
 grid on
 title('State Variable $x_2 = \Omega(t)$','Interpreter','latex')


%for z = 1:5
    q = 25;
    
xe                = kal(y,u,G,T_filter,Ts,L,x1,[pi*pi/3 0; 0 0],q);
xf                = stationary_kal(y,u,G,T_filter,Ts,L,x1,q);


%f13      = figure(13);set(f13,'name','Filter Result, ');

%subplot(211), plot(t,y,'-b',t,xe(1,:),'-r',t,xf(1,:),'--','LineWidth',1), hold on
%leg1    = legend('$\theta_n$','$\theta_n(kalman_filter)$', '$\theta_n(Stationary_kalman)$');
%set(leg1,'Interpreter','latex');
%grid on
%title('ouput of Kalman Filter and stationary kalman filter given $q \rightarrow [0.001]$','Interpreter','latex')

%subplot(212), plot(t,x(2,:),'-b',t,xe(2,:),'-r',t,xf(2,:),'--','LineWidth',1)
%leg1    = legend('$\Omega_n$', '$\Omega_n(kalman_filter)$', '$\Omega_n(Stationary_kalman)$');
%set(leg1,'Interpreter','latex');
%grid on
%title('Speed ($\Omega_n$) of Kalman Filter and stationary kalman filter given $q \rightarrow [0.001] $','Interpreter','latex')
%end


%The model of the system is perfect
%=======================================================================


x1                =  [0.1 0];
x1_stat           =  [0.1 0];
xe                = kal(y,u,G,T,Ts,L,x1,[pi*pi/3 0; 0 0],q);
xf                = stationary_kal(y,u,G,T,Ts,L,x1_stat,q);

figure(3), subplot(211), hold on,plot(t, xe(1,:)),plot(t, xf(1,:)), plot(t, y, t, x(1,:)),grid on  
leg1    = legend('$\theta_n(kalman_filter)$', '$\theta_n(Stationary_kalman)$','$\theta_n$');
set(leg1,'Interpreter','latex');
grid on
title('Perfect model-ouput of Kalman Filter and stationary kalman filter given $q \rightarrow [25]$','Interpreter','latex')

figure(3), subplot(212), hold on,plot(t, xe(2,:)),plot(t, xf(2,:)), plot(t,x(2,:)),grid on  
leg1    = legend( '$\Omega_n(kalman_filter)$', '$\Omega_n(Stationary_kalman)$','$\Omega_n$');
set(leg1,'Interpreter','latex');
grid on
title('Perfect-model-Speed ($\Omega_n$) of Kalman Filter and stationary kalman filter given $q \rightarrow [25] $','Interpreter','latex')


%The model of the system is rough
%========================================================================
%q                 = 6*1e3;
G_filter          = 50;


xe                = kal(y,u,G_filter,T_filter,Ts,L,x1,[pi*pi/3 0; 0 0],q);
xf                = stationary_kal(y,u,G_filter,T_filter,Ts,L,x1,q);

figure(4), subplot(211), hold on,plot(t, xe(1,:)),plot(t, xf(1,:)), plot(t, y, t, x(1,:)),grid on 
leg2    = legend('$\theta_n(kalman_filter)$', '$\theta_n(Stationary_kalman)$','$\theta_n$');
set(leg2,'Interpreter','latex');
grid on
title('model rough ouput of Kalman Filter and stationary kalman filter given $q \rightarrow [25]$','Interpreter','latex')

figure(4), subplot(212), hold on,plot(t, xe(2,:)),plot(t, xf(2,:)), plot(t,x(2,:)),grid on  
leg2    = legend( '$\Omega_n(kalman_filter)$', '$\Omega_n(Stationary_kalman)$','$\Omega_n$');
set(leg2,'Interpreter','latex');
grid on
title('model rough Speed ($\Omega_n$) of Kalman Filter and stationary kalman filter given $q \rightarrow [25] $','Interpreter','latex')




%%

function xe = stationary_kal(y,u,G,T,Ts,L,x1_0,q)  
    y_actual        = y;
    [sysd,A,B,C]    = getdigitizedsystem(G,T,Ts);
    N               = length(u);
    x               = zeros(2,N);
    x(:,1)          = x1_0';
    r               = 4*pi*pi/(12*L*L);
    
    F               = A;
    H               = C;
    R               = r;
    Q               = B*q*B';
    syms P [2 2];
    eqn             = P == F*P*F' - F*P*(H')*H*P*F'/(H*P*H' + R) + Q;
    
    %being a covariance matrix, P must be positive definite hence ->
    assume(P,'real');
    assume(P,'positive');
    
    P_sol           = vpasolve(eqn,[P1_1 P1_2 P2_1 P2_2]);
    P_inf           = [P_sol.P1_1(1) P_sol.P1_2(1); P_sol.P2_1(1) P_sol.P2_2(1)];
    
    
    K = P_inf*H'/(H*P_inf*H' + R);
    
    %now the Stationary Kalman Loop
    
    for n = 1:N-1
        
        y_pred      = C* x(:,n);
        
        y_observed  = y_actual(n) ;
        x(:,n)      = x(:,n) + K*(y_observed - y_pred);
        
        
        x(:,n+1)    = A*x(:,n) + B*(u(n)); 
    end
    
    xe = x;


end

function xe = kal(y,u,G,T,Ts,L,x1_0,p1_0,q)   
    y_actual        = y;
    [sysd,A,B,C]    = getdigitizedsystem(G,T,Ts);
    N               = length(u);
    x               = zeros(2,N);
    x(:,1)          = x1_0';
    P               = p1_0';
    r               = 4*pi*pi/(12*L*L);

 
    for n = 1:N-1
        
        y_pred      = C* x(:,n);
        Cxy         = P*C';
        Cyy         = C*P*C' + r;
        
        y_observed  = y_actual(n) ;
      
        x(:,n)      = x(:,n) + Cxy/Cyy*(y_observed - y_pred);
        P           = P   - Cxy*Cxy'/Cyy;
        
        
        x(:,n+1)    = A*x(:,n) + B*(u(n));
        P           = A*P*A'   + B*q*B';
        
        
    end
    xe = x;

end

function [yt, y,x] =  simulate(u,G,T,Ts,L,x1)
    [sysd,A,B,C]    = getdigitizedsystem(G,T,Ts);
    N               = length(u);
    x               = zeros(2,N);
    y               = zeros(1,N);
    x(:,1)          = x1';
    
  
    for n = 1:N-1
        x(:,n+1)    = A* x(:,n) + B*u(n);
        y(n)        = C* x(:,n);
        
    end
    y(N) = y(N-1);
    yt             = y;
     y             = round(y*L/2/pi)*2*pi/L;


end

function [sysd,A,B,C] = getdigitizedsystem(G,T,Ts)

    A_       = [0 1; 0 -1/T];
    B_       = [0; G/T];
    C_       = [1 0];
    sysc    = ss(A_,B_,C_,0);
    sysd    = c2d(sysc,Ts);
    A       = sysd.A;
    B       = sysd.B;
    C       = sysd.C;
    
end

function [u,t] = inputvoltage(D,A,Delta,Ts)
    t =(0:Ts:D);
    u = A/2*square(2*pi*t/Delta);
end


function r = get_quant_error()

r = abs(randn); 

end
