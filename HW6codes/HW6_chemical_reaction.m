close all;
clc;
h = 0.001;     % Step size
Tmax = 100000;    % Maximum time
N = Tmax / h;  % Maximum number of steps
alpha=0.5;
t = linspace(0,Tmax,N+1);  % Time range
k1 = 2000;
k2 = 0.001;
k3 = 10;


%Numerical solution
f=@(t,c) [-k1*c(1)*c(2)+k2*c(3); -k1*c(1)*c(2)+(k2+k3)*c(3);k1*c(1)*c(2)-(k2+k3)*c(3);k3*c(3)]; % Governing system of equations

% Initial Concentrations
c = [1; 0.00005;0;0];

if 0
    % Fourth order Runge-Kutta method
    for i=1:N
        RKk1 = f(t(i),c(:,i));
        RKk2 = f(t(i)+0.5.*h, c(:,i)+0.5.*RKk1*h);
        RKk3 = f(t(i)+0.5.*h, c(:,i)+0.5.*RKk2*h);
        RKk4 = f(t(i)+h,c(:,i)+RKk3*h);
        c(:,i+1) = c(:,i) + (1/6)*(RKk1+2*RKk2+2*RKk3+RKk4)*h;
    end
    figure
    set(gcf,'position',[400,0.5,600,750])
    subplot(4,1,1)
    loglog(t,c(1,:),'.:')
    title(['C_{S} concentration, RK4 method steps, When h = ' num2str(h)])
    xlabel('t')
    ylabel('C_{S}')

    subplot(4,1,2)
    loglog(t,c(2,:),'.:')
    title(['C_{E} concentration, RK4 method steps, When h = ' num2str(h)])
    xlabel('t')
    ylabel('C_{E}')

    subplot(4,1,3)
    loglog(t,c(3,:),'.:')
    title(['C_{ES} concentration, RK4 method steps, When h = ' num2str(h)])
    xlabel('t')
    ylabel('C_{ES}')

    subplot(4,1,4)
    loglog(t,c(4,:),'.:')
    title(['C_{P} concentration, RK4 method steps, When h = ' num2str(h)])
    xlabel('t')
    ylabel('C_{P}')
end

if 1
    %ode23s
    [t,y] = ode23s(f,[0 Tmax],c);
    figure
    set(gcf,'position',[400,0.5,600,750])
    subplot(4,1,1)
    loglog(t,y(:,1),'.:')
    title(['C_{S} concentration, ode23s method'])
    xlabel('t')
    ylabel('C_{S}')

    subplot(4,1,2)
    loglog(t,y(:,2),'.:')
    title(['C_{E} concentration, ode23s method'])
    xlabel('t')
    ylabel('C_{E}')

    subplot(4,1,3)
    loglog(t,y(:,3),'.:')
    title(['C_{ES} concentration, ode23s method'])
    xlabel('t')
    ylabel('C_{ES}')

    subplot(4,1,4)
    loglog(t,y(:,4),'.:')
    title(['C_{P} concentration, ode23s method'])
    xlabel('t')
    ylabel('C_{P}')
end
