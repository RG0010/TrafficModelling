%% File information
%{
MIT License

Copyright (c) 2023 RG0010

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Author:  RG0010
Title:  GadunovSolver1_3_1_FINAL.m
Version:  3.0
Date:  13 May 2023
Description:

    This program solves the following partial differential equation:
            rho_t + f(rho)_x = s(x,t)

    This program uses Gadunov's method.  The algorithm is:
        1) Define delta_x, delta_y, x_max, & t_max
        2) Define initial & boundary conditions
            i) Initial conditions (rho_x_t_min)
           ii) Boundary conditions (rho_x_min_t & rho_x_max_t)
        3) Compute the apprx. flow at cell boundaries
        4) Calculate the next density in time using values from step 3
        5) Repeat for all values of time and x

    INPUTS (some in initial conditions file):
        rho_x_t_min - distance boundary conditions (veh/m)
        rho_x_min_t - time boundary conditions when x=x_min (veh/m)
        rho_x_max_t - time boundare conditions when x=x_max (veh/m)
        dx - width of distance step (m)
        dt - width of time step (s)
        delta_p - step size for rho to compute Riemann problem solution
        rho_max - traffic jame density
        v_f - maximum traffic speed
        x_min - minimum x value (m)
        x_max - maximum x value (m)
        t_min - minimum time value (s)
        t_max - maximum time value (s)
    OUTPUTS:
        rho - nxj matrix of solutions to LWR model (veh/m)
    FUNCTIONS:
        F(rho_l,rho_r) - function which solves the Riemann Problem at each
                        cell interface
        f(rho) - function for the flow (a function of density only) (veh/s)
        CFL(delta_t, delta_x) - function to check CFL condition
        V(rho) - fundamental diagram (in this case, Greenshields)

    HOW TO USE:
        To use this code, call a separate file of initial conditions
        defining rho_x_t_min, rho_x_min_t, and rho_x_max_t. If using one of
        mine, uncomment one case and ensure all other cases are commented.
%}
%% Initial tidying
clear;clc;clf;
%% Define consts. for bigger practice problem
%{%}
global v_f;
global rho_max;
    v_f = 30;       % m/s
    rho_max = 0.15;    % density of traffic at full jam (veh/m)

global delta_p;
    
    delta_p = 0.0001;
% Integration constants
dx = 32;     % define step value along x axis (km)
dt = 1;     % define stap value along t axis (hrs)
x_min = 0;
x_max = 12000;       % define max x value (mi)
t_min = 0;
t_max = 360;      % define max t value (hrs)

x = x_min:dx:x_max;    % initialize x values array
t = t_min:dt:t_max;            % initialize t values array


% Initialize matrix P
rho = zeros(length(x),length(t));
q = zeros(length(x),length(t));

Init_Conds_2;
%Init_Conds_LWR_Selivanov

%% Check CFL condition
C = CFL(dt,dx);
if C > 1
    error("Error: Solution will not converge.  CFL condition is " + C + ".");
else    
    disp("CFL condition is " + C + ".");
end

%% Computation of P
% Solve equation P(n+1,j) = P(n,j) - (delta_t/delta_x)*(F_pos - F_neg)
for n = 1:1:length(t)-1               % iterate thru all t values
    for j = 2:1:length(x)-1 % iterate thru all x values
        
        F_pos = F(rho(j,n),rho(j+1,n));
        F_neg = F(rho(j-1,n),rho(j,n));
        rho(j,n+1) = rho(j,n) - (dt/dx)*(F_pos - F_neg);
    end
    rho(length(x),n) = rho(length(x)-1,n);
        rho(1,n) = rho(2,n);
        q(length(x),n) = q(length(x)-1,n);
        q(1,n) = q(2,n);
end

%% Plotting code - uncomment as needed
%{
surf(x,t,rho','LineStyle','none'); 
limits1 = [0,rho_max];
c1 = colorbar;
c1.Label.String = "Traffic Density (m/s)";
c1.Label.FontSize = 12;
%set(gca,'clim',limits1);
title("Traffic density by time and distance - LWR", 'FontSize', 12);
axis([x_min x_max t_min t_max 0 rho_max]);
xlabel("Distance (m)", 'FontSize', 12);
ylabel("Time (s)", 'FontSize', 12);
zlabel("Density (veh/m)", 'FontSize', 12);
%}
%{
% Plots the traffic density for all t and x
figure(1)
subplot(1,2,1)
surf(x,t,rho','LineStyle','none'); 
limits1 = [0,rho_max];
c1 = colorbar;
c1.Label.String = "Traffic Density (m/s)";
c1.Label.FontSize = 12;
set(gca,'clim',limits1);
title("Traffic density by time and distance - LWR", 'FontSize', 12);
axis([x_min x_max t_min t_max 0 rho_max]);
xlabel("Distance (m)", 'FontSize', 12);
ylabel("Time (s)", 'FontSize', 12);
zlabel("Density (veh/m)", 'FontSize', 12);
% Traffic speed by time and distance plot
subplot(1,2,2)
surf(x,t,(f(rho)./rho)','LineStyle','none');
limits2 = [0,v_f];
c2 = colorbar;
c2.Label.String = "Traffic Velocity (m/s)";
c2.Label.FontSize = 12;
set(gca,'clim',limits2);
title("Traffic velocity by time and distance - LWR", 'FontSize', 12);
axis([x_min x_max t_min t_max 0 v_f]);
xlabel("Distance (m)", 'FontSize', 12);
ylabel("Time (s)", 'FontSize', 12);
zlabel("Speed (mph)", 'FontSize', 12);
%}
%{
% Traffic flow by time and distance plot
figure(3);
subplot(2,1,1);
plot(x,rho(:,find(t==50)));
title("LWR Traffic Density at t=50s", 'FontSize', 12);
xlabel("Position (m)", 'FontSize', 12);
ylabel("Traffic density (veh/m)", 'FontSize', 12);
axis([x_min x_max 0 (rho_max+0.005)]);
grid on
subplot(2,1,2);
plot(x,rho(:,find(t==150)));
title("LWR Traffic Density at t=150s", 'FontSize', 12);
xlabel("Position (m)", 'FontSize', 12);
ylabel("Traffic density (veh/m)", 'FontSize', 12);
axis([x_min x_max 0 (rho_max+0.005)]);
grid on

% Plots fundamental diagram
p = 0:0.001:rho_max;
figure(10);
%{
subplot(3,1,1)
plot(f(p),f(p)./p);
title("Fundamental Diagrams");
xlabel("Flow (veh/hr)");
ylabel("Velocity (km/hr)");
subplot(3,1,2)
plot(p,f(p)./p);
xlabel("Density (veh/km)");
ylabel("Velocity (km/hr)");
subplot(3,1,3)
%}
plot(p,f(p));
grid on;
title("Greenshields Fundamental Diagram Used for All Modelling",'FontSize',12);
xlabel("Density (veh/km)",'FontSize',12);
ylabel("Flow (veh/hr)",'FontSize',12);
axis([0 rho_max 0 (max(f(p))+0.08*max(f(p)))]);
%}
%{
% Plots the traffic density for all t and x
figure(1)
subplot(1,2,1);
surf(x,t,P,'LineStyle','none'); 
title("Traffic density by time and distance");
axis([1 5 0 5 0 1]);
xlabel("Distance (km)");
ylabel("Time (hr)");
zlabel("Density (veh/km)");
subplot(1,2,2);
surf(x,t,P,'LineStyle','none'); 
title("Traffic density by time and distance");
axis([1 5 0 5 0 1]);
xlabel("Distance (km)");
ylabel("Time (hr)");
zlabel("Density (veh/km)");
% Traffic speed by time and distance plot
figure(2)
surf(x,t,f(P)./P,'LineStyle','none'); 
title("Traffic speed by time and distance");
axis([1 5 0 5 0 1]);
xlabel("Distance (km)");
ylabel("Time (hr)");
zlabel("Speed (km/hr)");

% Traffic flow by time and distance plot
figure(3)
surf(x,t,f(P),'LineStyle','none'); 
title("Traffic flow by time and distance");
axis([1 5 0 5 0 1]);
xlabel("Distance (km)");
ylabel("Time (hr)");
zlabel("Flow (veh/hr)");

% Plots fundamental diagram
p = 0:0.001:1;
figure(10);
subplot(3,1,1)
plot(f(p),f(p)./p);
title("Fundamental Diagrams");
xlabel("Flow (veh/hr)");
ylabel("Velocity (km/hr)");
subplot(3,1,2)
plot(p,f(p)./p);
xlabel("Density (veh/km)");
ylabel("Velocity (km/hr)");
subplot(3,1,3)
plot(p,f(p));
xlabel("Density (veh/km)");
ylabel("Flow (veh/hr)");
%}
%% Functions
function soln = F(P_l, P_r)         % This function solves the Riemann Problems
    global delta_p;
    if P_r < P_l
        p = P_r:delta_p:P_l;
        f_lr = f(p);
        soln = max(f_lr);
    elseif P_l < P_r
        p = P_l:delta_p:P_r;
        f_lr = f(p);
        soln = min(f_lr);
    else
        soln = f(P_l);
    end
end

function Q = f(p)                   % This is the function for the flow
                            % This fcn assumes the Q(p) = p * V(p)
    global v_f;
    global rho_max;
    %V = v_f * (1 - (p./rho_max));         % Greenshields model
    Q = p .* V(p);
end

function C = CFL(delta_t, delta_x)  % Checks CFL condition
    global rho_max;
    global delta_p;
    p = 0:delta_p:rho_max;
    Q = f(p);
    for index = 1:1:length(Q)-1
        Q_p = (Q(index+1) - Q(index)) / delta_p;
    end
    a = max(abs(Q_p));
    C = a * (delta_t /delta_x);
end
% Greenshields model
function vel = V(rho)
    n = 1;
    global v_f;
    global rho_max;
    vel = v_f .* (1 - (rho ./ rho_max).^n);
end