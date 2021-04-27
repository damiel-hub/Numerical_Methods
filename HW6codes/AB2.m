close all;
clc;
h = 0.5;     % Step size
Tmax = 6;    % Maximum time
N = Tmax / h;  % Maximum number of steps
alpha=0.5;
t = linspace(0,Tmax,N+1);  % Time range
g = 9.81;
l = 0.6;

% Analytical solution of the differential equation
treal = 0:0.01:Tmax;
theta_real = 10*cos(sqrt(g/l)*treal);
plot(treal,theta_real);
hold on

%Numerical solution
f=@(t,theta) [theta(2); (-g/l)*theta(1)]; % Governing system of equations

% Initial Conditions
Theta = [10; 0];
% Initialization with second order Runge-Kutta method
k1 = h.*f(t(1),Theta(:,1));
k2 = h.*f(t(1)+alpha.*h, Theta(:,1)+alpha.*k1);
Theta(:,2) = Theta(:,1) + (1-1/2/alpha).*k1 + k2/2/alpha;

% AB2 method steps
for i=2:N
    Theta(:,i+1) = Theta(:,i) + (3/2).*h.*f(t(i),Theta(:,i)) - (1/2).*h.*f(t(i-1),Theta(:,i-1));
end
plot(t,Theta(1,:),'.:')
legend('Exact Solution','Leapfrog Solution','Location','NorthEast')
title(['AB2 method steps, When h = ' num2str(h)])
xlabel('t')
ylabel('\theta')