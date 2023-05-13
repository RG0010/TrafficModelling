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
Title:  GadunovSolver1_2_FINAL.m
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
        rho_x_t_min - distance boundary traffic density conditions (veh/m)
        rho_x_min_t - time boundary traffic density conditions when x=x_min (veh/m)
        rho_x_max_t - time boundare traffic density conditions when x=x_max (veh/m)
        q_x_t_min - distance boundary traffic flow deviation conditions (veh/s)
        q_x_min_t - time boundary traffic flow deviation conditions when x=x_min (veh/s)
        q_x_max_t - time boundare traffic flow deviation conditions when x=x_max (veh/s)
        S - source term definitions
        tau - relaxation constant
        dx - width of distance step (m)
        dt - width of time step (s)
        x_min - minimum x value (m)
        x_max - maximum x value (m)
        t_min - minimum time value (s)
        t_max - maximum time value (s)
    OUTPUTS:
        rho - nxj matrix of traffic density solutions to ARZ model (veh/m)
        v - nxj matrix of traffic velocity solutions to ARZ model (m/s)
        q - nxj matrix of traffic flow deviation solutions to ARZ model (veh/s)
    FUNCTIONS:
        F_1(U_L, U_R) - function which solves the Riemann Problem at each
                        cell interface for solving for rho
        F_2(U_L, U_R) - function which solves the Riemann Problem at each
                        cell interface for solving for q
        f_1(U) - function for the traffic density flux (a function of density & flow deviation) (veh/s)
        f_2(U) - function for the traffic flow deviation flux (a function of density & flow deviation) (veh*m/s/s)
        V(rho) - fundamental diagram (in this case, Greenshields)
        V_prime(rho) - derivative of fundamental diagram w.r.t. rho 
    HOW TO USE:
        To use this code, call a separate file of initial conditions
        defining rho_x_t_min, rho_x_min_t, rho_x_max_t, q_x_t_min, q_x_min_t,
        and q_x_max_t. If using one of mine, uncomment one case and ensure 
        all other cases are commented.

%}
%% setup
clc;clear;clf;
format long
%% Define constants and axes
global v_f;
global rho_max;
    v_f = 30;       % m/s
    rho_max = 0.15;    % density of traffic at full jam (veh/m)

% step values (units being meters and seconds)
dx = 32;     % meters    
dt = 1;   
tau = 1e100;
% Max and min x values
x_min = 0.0;
x_max = 12000.0;

% some remarks on macroscopic traffic flow modelling
% Max and min t values
t_min = 0.0;
t_max = 360.0;

% definitions for x and t axes
t = t_min:dt:t_max;
x = x_min:dx:x_max;

%% Set initial conditions
% Page 1 is rho values, page 2 is flow deviation
rho = zeros(length(x),length(t));
v = zeros(length(x),length(t));
q = zeros(length(x),length(t));

Init_Conds_2;
%Init_Conds_LWR_Selivanov

for n = 1:1:(length(t) - 1)
    for i = 2:1:(length(x) - 1)
        
        U_L = [rho(i-1,n), q(i-1,n)];
        U_M = [rho(i,n)  , q(i,n)];
        U_R = [rho(i+1,n), q(i+1,n)];
        C_r = max(eig_1(U_M),eig_2(U_M))*(dt/dx);
        if (C_r > 1) || (isnan(C_r)) 
            disp(C_r);
            error("CFL Violation.");
        end
        F_1_L = F_1(U_L, U_M);
        F_1_R = F_1(U_M, U_R);
        F_2_L = F_2(U_L, U_M);
        F_2_R = F_2(U_M, U_R);
        rho(i,n+1) = U_M(1) - (F_1_R - F_1_L)*(dt/dx) + (dt/2)*(S(i,n+1) + S(i,n));
        q(i,n+1)   = (2*tau/(dt+2*tau))*((U_M(2) - (F_2_R - F_2_L)*(dt/dx)) - U_M(2)*(dt/(2*tau))); 
        if (isnan(rho(i,n+1)))
            error("Answer is NaN!!");
        end
    end
        rho(length(x),n) = rho(length(x)-1,n);
        rho(1,n) = rho(2,n);
        q(length(x),n) = q(length(x)-1,n);
        q(1,n) = q(2,n);
end

v = (q)./(rho) + V(rho);
figure(1)
contourf(t,x,rho,'LineStyle','none');
title("Traffic density by time and distance - ARZ", 'FontSize', 18);
ylabel("Distance (m)", 'FontSize', 18);
xlabel("Time (s)", 'FontSize', 18);
zlabel("Traffic Density (veh/m)", 'FontSize', 12);
c3 = colorbar;
%ylabel(c3, "Traffic Density (veh/m)", "Rotation", 270, 'FontSize', 12);
c3.Label.String = "Traffic Densite (veh/m)";
c3.Label.FontSize = 12;
limits1=[0.0,0.15];
set(gca,'clim',limits1);
ylim([x_min x_max]);
xlim([t_min t_max]);

figure(50)
subplot(1,2,1);
plot(x,rho(:,find(t==50)));
title("ARZ Traffic Density at t=50s", 'FontSize', 12);
xlabel("Position (m)", 'FontSize', 12);
ylabel("Traffic density (veh/m)", 'FontSize', 12);
axis([3000 9600 0 0.16]);
grid on
subplot(1,2,2);
plot(x,rho(:,find(t==150)));
title("ARZ Traffic Density at t=150s", 'FontSize', 12);
xlabel("Position (m)", 'FontSize', 12);
ylabel("Traffic density (veh/m)", 'FontSize', 12);
axis([1500 10300 0 0.16]);
grid on

figure(2);
subplot(1,2,1);
surf(x,t,rho','LineStyle','none');
title("Traffic density by time and distance - ARZ", 'FontSize', 12);
xlabel("Distance (m)", 'FontSize', 12);
ylabel("Time (s)", 'FontSize', 12);
zlabel("Traffic Density (veh/m)", 'FontSize', 12);
c3 = colorbar;
%ylabel(c3, "Traffic Density (veh/m)", "Rotation", 270, 'FontSize', 12);
c3.Label.String = "Traffic Densite (veh/m)";
c3.Label.FontSize = 12;
limits1=[0,rho_max];
set(gca,'clim',limits1);
xlim([x_min x_max]);
ylim([t_min t_max]);
subplot(1,2,2)
surf(x,t,v','LineStyle','none');
title("Traffic velocity by time and distance - ARZ", 'FontSize', 12);
xlabel("Distance (m)", 'FontSize', 12);
ylabel("Time (s)", 'FontSize', 12);
zlabel("Speed (m/s)", 'FontSize', 12);
colorbar
limits2 = [0,v_f];
c4 = colorbar; %('Ticks',[0:(v_f/5):v_f]);
%ylabel(c4, "Traffic Velocity (m/s)", "Rotation", 270, 'FontSize', 12);
c4.Label.FontSize = 12;
c4.Label.String = "Traffic Velocity (m/s)";
%c4.Label.Rotation = 0;
set(gca,'clim',limits2);
xlim([x_min x_max]);
ylim([t_min t_max]);

figure(3);
subplot(2,1,1);
plot(x,rho(:,find(t==50)));
title("ARZ Traffic Density at t=50s", 'FontSize', 12);
xlabel("Position (m)", 'FontSize', 12);
ylabel("Traffic density (veh/m)", 'FontSize', 12);
axis([x_min x_max 0 (rho_max+0.005)]);
grid on
subplot(2,1,2);
plot(x,rho(:,find(t==150)));
title("ARZ Traffic Density at t=150s", 'FontSize', 12);
xlabel("Position (m)", 'FontSize', 12);
ylabel("Traffic density (veh/m)", 'FontSize', 12);
axis([x_min x_max 0 (rho_max+0.005)]);
grid on

figure(20);
plot(x,rho(:,find(t==t_min)),'LineWidth',2);
hold on;
plot(x,rho(:,find(t==100)),'LineWidth',2);
plot(x,rho(:,find(t==200)),'LineWidth',2);
plot(x,rho(:,find(t==300)),'LineWidth',2);
hold off;
legend('T=0','T=100','T=200','T=300');
title("ARZ Highway Traffic Density at Various Times", 'FontSize', 12);
xlabel("Position (m)", 'FontSize', 12);
ylabel("Traffic density (veh/m)", 'FontSize', 12);
axis([x_min x_max 0 (rho_max+0.005)]);
grid on

function f_LR = F_1(U_L, U_R)
    lambda_1_L = eig_1(U_L);
    lambda_1_R = eig_1(U_R);
    lambda_2_L = eig_2(U_L);
    lambda_2_R = eig_2(U_R);
    
    %{
    % A different method of approximating wave speeds
    a_L = lambda_1_L; 
    a_R = lambda_2_R; 
    %}
    %{
    % A different method of approximating wave speeds
    a_L = min(lambda_1_L, lambda_1_R);
    a_R = max(lambda_2_L, lambda_2_R);
    %}
    a_R = max([lambda_1_L,lambda_1_R,lambda_2_L,lambda_2_R]);
    a_L = -a_R;
    
    c_1 = 0.5*((a_R - abs(a_R)) - (a_L - abs(a_L))) / (a_R - a_L);
    c_2 = 0.5*((a_R + abs(a_R)) - (a_L + abs(a_L))) / (a_R - a_L);
    c_3 = -0.5*((a_R*abs(a_L)) - (a_L*abs(a_R))) / (a_R - a_L);
    
    f_L = f_1(U_L);
    f_R = f_1(U_R);
    
    f_LR = c_1*f_R + c_2*f_L + c_3*(U_R(1) - U_L(1));
end

function f_LR = F_2(U_L, U_R) % p. 451 of Davis 1988 **WORKS**
    lambda_1_L = eig_1(U_L);
    lambda_1_R = eig_1(U_R);
    lambda_2_L = eig_2(U_L);
    lambda_2_R = eig_2(U_R);
    
    %{
    % A different method of approximating wave speeds
    a_L = lambda_1_L; 
    a_R = lambda_2_R; 
    %}
    %{
    % A different method of approximating wave speeds
    a_L = min(lambda_1_L, lambda_1_R);
    a_R = max(lambda_2_L, lambda_2_R);
    %}
    a_R = max([lambda_1_L, lambda_1_R,lambda_2_L, lambda_2_R]);
    a_L = -a_R;
    
    c_1 = 0.5*((a_R - abs(a_R)) - (a_L - abs(a_L))) / (a_R - a_L);
    c_2 = 0.5*((a_R + abs(a_R)) - (a_L + abs(a_L))) / (a_R - a_L);
    c_3 = -0.5*((a_R*abs(a_L)) - (a_L*abs(a_R))) / (a_R - a_L);
    
    f_L = f_2(U_L);
    f_R = f_2(U_R);
    
    f_LR = c_1*f_R + c_2*f_L + c_3*(U_R(2) - U_L(2));
end

function flux = f_1(U)
    flux = U(2) + U(1)*V(U(1));
end

function flux = f_2(U)
    flux = ((U(2)^2)/U(1)) + U(2)*V(U(1));
end

function lambda = eig_1(U)
    lambda = V(U(1)) + U(1)*V_prime(U(1)) + (U(2)/U(1));
end

function lambda = eig_2(U)
    lambda = V(U(1)) + (U(2)/U(1));
end

% Greenshields model
function vel = V(rho)
    n = 1;
    global v_f;
    global rho_max;
    vel = v_f .* (1 - (rho ./ rho_max).^n);
end

% Partial derivative of Greenshields model w.r.t. rho
function vel_prime = V_prime(rho)
    n = 1;
    global v_f;
    global rho_max;
    vel_prime = -n*v_f*(rho/rho_max)^n;
end