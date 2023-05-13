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
Title:  Init_Conds_2.m
Version:  3.0
Date:  13 May 2023

    TO USE:
        Uncomment the desired test case and then run either
        ARZsolver1_2_1_FINAL.m or GadunovSolver1_3_1_FINAL.m
%}
%% Case 1
% Equilibrium velocity, and situation where faster, sparser traffic catches
% up to slower traffic.  I.e., a transition from low to high traffic
% density with free flow traffic
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.1*rho_max*ones(1,length(x_min:dx:((x_max/2)))), ... % 0.1 if x < 1
               0.45*rho_max*ones(1,length((x_max/2):dx:x_max))];           % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

q_x_t_min = 0*ones(1,length(x_min:dx:x_max)); % 0 v(x,0) = V(rho)
q_x_min_t = q_x_t_min(1)*ones(length(t),1);
q_x_max_t = q_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;

q(:,1) = q_x_t_min;
q(1,:) = q_x_min_t;
q(length(x),:) = q_x_max_t;
%}
%% Case 2
% Equilibrium velocity, and situation where faster, sparser traffic catches
% up to slower traffic.  I.e., a transition from low to high traffic
% density with congested traffic (dt = 0.1)
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.55*rho_max*ones(1,length(x_min:dx:((x_max/2)))), ... % 0.1 if x < 1
               0.9*rho_max*ones(1,length((x_max/2):dx:x_max))];           % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

q_x_t_min = 0*ones(1,length(x_min:dx:x_max)); % 0 v(x,0) = V(rho)
q_x_min_t = q_x_t_min(1)*ones(length(t),1);
q_x_max_t = q_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;

q(:,1) = q_x_t_min;
q(1,:) = q_x_min_t;
q(length(x),:) = q_x_max_t;
%}
%% Case 3
% Equilibrium velocity, and situation where slower traffic expands into 
% faster, sparser traffic.  I.e., a transition from high to low traffic
% density with free flow traffic
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.45*rho_max*ones(1,length(x_min:dx:((x_max/2)))), ... % 0.1 if x < 1
               0.1*rho_max*ones(1,length((x_max/2):dx:x_max))];           % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

q_x_t_min = 0*ones(1,length(x_min:dx:x_max)); % 0 v(x,0) = V(rho)
q_x_min_t = q_x_t_min(1)*ones(length(t),1);
q_x_max_t = q_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;

q(:,1) = q_x_t_min;
q(1,:) = q_x_min_t;
q(length(x),:) = q_x_max_t;
%}
%% Case 4
% Equilibrium velocity, and situation where slower traffic expands into 
% faster, sparser traffic.  I.e., a transition from high to low traffic
% density with congested traffic traffic
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.9*rho_max*ones(1,length(x_min:dx:((x_max/2)))), ... % 0.1 if x < 1
               0.55*rho_max*ones(1,length((x_max/2):dx:x_max))];           % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

q_x_t_min = 0*ones(1,length(x_min:dx:x_max)); % 0 v(x,0) = V(rho)
q_x_min_t = q_x_t_min(1)*ones(length(t),1);
q_x_max_t = q_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;

q(:,1) = q_x_t_min;
q(1,:) = q_x_min_t;
q(length(x),:) = q_x_max_t;
%}
%% Case 5
% Equilibrium velocity, and situation where there is a section of jammed
% traffic that expands into regions of lowere traffic density.  I.e.,
% traffic jam dissolution. (dt=0.1)
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.35*rho_max*ones(1,length(x_min:dx:(x_max/3)-dx)), ... % 0.1 if x < 1
               0.99*rho_max*ones(1,length(((x_max/3)-dx):dx:(2*x_max/3)-2*dx)), ...
               0.35*rho_max*ones(1,length((2*x_max/3):dx:x_max)),];           % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

q_x_t_min = 0*ones(1,length(x_min:dx:x_max)); % 0 v(x,0) = V(rho)
q_x_min_t = q_x_t_min(1)*ones(length(t),1);
q_x_max_t = q_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;

q(:,1) = q_x_t_min;
q(1,:) = q_x_min_t;
q(length(x),:) = q_x_max_t;
%}
%% Case 6
% Equilibrium velocity, and situation where there is a section of non-jammed
% traffic in which regions of high traffic density expand.  I.e.,
% traffic jam dissolution into a small region of non-jammed traffic.
% (dt=0.1)
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.99*rho_max*ones(1,length(x_min:dx:(x_max/3)-dx)), ... % 0.1 if x < 1
               0.35*rho_max*ones(1,length(((x_max/3)-dx):dx:(2*x_max/3)-2*dx)), ...
               0.99*rho_max*ones(1,length((2*x_max/3):dx:x_max)),];           % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

q_x_t_min = 0*ones(1,length(x_min:dx:x_max)); % 0 v(x,0) = V(rho)
q_x_min_t = q_x_t_min(1)*ones(length(t),1);
q_x_max_t = q_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;

q(:,1) = q_x_t_min;
q(1,:) = q_x_min_t;
q(length(x),:) = q_x_max_t;
%}
%% Case 7
% A middle hump with initial velocity as defined in Zhang (2002) in
% free-flow area of fundamental diagram
% I.e. velocity before the hump is equal to the velocity during the hump
% I.e. v1 = v2 = v
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.1*rho_max*ones(1,length(x_min:dx:((1*x_max)/6))), ... % 0.1 if x < 1
           0.4*rho_max*ones(1,length(((1*x_max)/6):dx:((2*x_max)/6)-dx)), ...     % 0.9 if 1 <= x < 2
           0.1*rho_max*ones(1,length(((2*x_max)/6):dx:x_max)),];           % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;

v = V(min(min(rho_x_t_min)))*ones(1,length(x_min:dx:x_max));

q_x_t_min = rho_x_t_min.*(v-V(rho_x_t_min));
q_x_min_t = q_x_t_min(1)*ones(length(t),1);
q_x_max_t = q_x_t_min(length(x))*ones(length(t),1);

q(:,1) = q_x_t_min;
q(1,:) = q_x_min_t;
q(length(x),:) = q_x_max_t;
%}
%% Case 8
% A middle hump with initial velocity as defined in Zhang (2002) in
% congested area of fundamental diagram
% I.e. velocity before the hump is equal to the velocity during the hump
% I.e. v1 = v2 = v
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.55*rho_max*ones(1,length(x_min:dx:((1*x_max)/6))), ... % 0.1 if x < 1
           0.9*rho_max*ones(1,length(((1*x_max)/6):dx:((2*x_max)/6)-dx)), ...     % 0.9 if 1 <= x < 2
           0.55*rho_max*ones(1,length(((2*x_max)/6):dx:x_max)),];           % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;

v = V(min(min(rho_x_t_min)))*ones(1,length(x_min:dx:x_max));

q_x_t_min = rho_x_t_min.*(v-V(rho_x_t_min));
q_x_min_t = q_x_t_min(1)*ones(length(t),1);
q_x_max_t = q_x_t_min(length(x))*ones(length(t),1);

q(:,1) = q_x_t_min;
q(1,:) = q_x_min_t;
q(length(x),:) = q_x_max_t;
%}
%% Case 9
% Phantom traffic jam that interacts with equilibrium traffic upstream
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.30*rho_max*ones(1,length(x_min:dx:((1*x_max)/6))), ... % 0.1 if x < 1
           0.6*rho_max*ones(1,length(((1*x_max)/6):dx:((2*x_max)/6)-dx)), ...     % 0.9 if 1 <= x < 2
           0.30*rho_max*ones(1,length(((2*x_max)/6):dx:((4*x_max)/6)-dx)),...
           0.60*rho_max*ones(1,length(((4*x_max)/6):dx:x_max)),];           % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;

v = [V(min(min(rho_x_t_min)))*ones(1,length(x_min:dx:((4*x_max)/6)-dx)), ...
     V(rho_x_t_min(find(x==(4*x_max/6)):end))];


q_x_t_min = rho_x_t_min.*(v-V(rho_x_t_min));
q_x_min_t = q_x_t_min(1)*ones(length(t),1);
q_x_max_t = q_x_t_min(length(x))*ones(length(t),1);

q(:,1) = q_x_t_min;
q(1,:) = q_x_min_t;
q(length(x),:) = q_x_max_t;
%}
%% Case 10
% Phantom traffic jam in congested flow that interacts with congested equilibrium
% traffic upstream
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.55*rho_max*ones(1,length(x_min:dx:((1*x_max)/6))), ... % 0.1 if x < 1
           0.7*rho_max*ones(1,length(((1*x_max)/6):dx:((2*x_max)/6)-dx)), ...     % 0.9 if 1 <= x < 2
           0.55*rho_max*ones(1,length(((2*x_max)/6):dx:((4*x_max)/6)-dx)),...
           0.75*rho_max*ones(1,length(((4*x_max)/6):dx:x_max)),];           % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;

v = [V(min(min(rho_x_t_min)))*ones(1,length(x_min:dx:((4*x_max)/6)-dx)), ...
     V(rho_x_t_min(find(x==(4*x_max/6)):end))];


q_x_t_min = rho_x_t_min.*(v-V(rho_x_t_min));
q_x_min_t = q_x_t_min(1)*ones(length(t),1);
q_x_max_t = q_x_t_min(length(x))*ones(length(t),1);

q(:,1) = q_x_t_min;
q(1,:) = q_x_min_t;
q(length(x),:) = q_x_max_t;
%}
%% Case 11
% 
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.4*rho_max*ones(1,length(x))];           % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;

v = [V(0.3*rho_max)*ones(1,length(x_min:dx:((x_max)/3)-dx)), ...
     V(0.5*rho_max)*ones(1,length(((1*x_max)/3):dx:((2*x_max)/3)-dx)), ...
     V(0.7*rho_max)*ones(1,length(((2*x_max)/3):dx:x_max))];

q_x_t_min = rho_x_t_min.*(v-V(rho_x_t_min));
q_x_min_t = q_x_t_min(1)*ones(length(t),1);
q_x_max_t = q_x_t_min(length(x))*ones(length(t),1);

q(:,1) = q_x_t_min;
q(1,:) = q_x_min_t;
q(length(x),:) = q_x_max_t;
%}
%% Case 12
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.6*rho_max*ones(1,length(x))];           % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;

v = [V(0.3*rho_max)*ones(1,length(x_min:dx:((x_max)/3)-dx)), ...
     V(0.5*rho_max)*ones(1,length(((1*x_max)/3):dx:((2*x_max)/3)-dx)), ...
     V(0.7*rho_max)*ones(1,length(((2*x_max)/3):dx:x_max))];

q_x_t_min = rho_x_t_min.*(v-V(rho_x_t_min));
q_x_min_t = q_x_t_min(1)*ones(length(t),1);
q_x_max_t = q_x_t_min(length(x))*ones(length(t),1);

q(:,1) = q_x_t_min;
q(1,:) = q_x_min_t;
q(length(x),:) = q_x_max_t;
%}
% Greenshields model
function vel = V(rho)
    n = 1;
    global v_f;
    global rho_max;
    vel = v_f .* (1 - (rho ./ rho_max).^n);
end