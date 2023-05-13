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
Title:  Init_Conds_LWR_Selivanov.m
Version:  3.0
Date:  13 May 2023

    TO USE:
        Uncomment the desired test case and then run either
        ARZsolver1_2_1_FINAL.m or GadunovSolver1_3_1_FINAL.m
%}
%% Case 1
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0*rho_max*ones(1,length(x_min:dx:0)), ... % 0.1 if x < 1
           0.5*rho_max*ones(1,length((0+dx:dx:x_max)))];     % 0.9 if 1 <= x < 2          % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;
%}
%% Case 2
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.8*rho_max*ones(1,length(x_min:dx:0)), ... % 0.1 if x < 1
           0*rho_max*ones(1,length((0+dx:dx:x_max)))];     % 0.9 if 1 <= x < 2          % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;
%}
%% Case 3
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.4*rho_max*ones(1,length(x_min:dx:1)), ... % 0.1 if x < 1
               0.7*rho_max*ones(1,length(1+dx:dx:2)), ... % 0.1 if x < 1
           0.9*rho_max*ones(1,length((2+dx:dx:x_max)))];     % 0.9 if 1 <= x < 2          % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;
%}
%% Case 4
%{
S = zeros(length(x),length(t)); % Source always equals 0

rho_x_t_min = [0.1*rho_max*ones(1,length(x_min:dx:1)), ... % 0.1 if x < 1
               0.9*rho_max*ones(1,length(1+dx:dx:2)), ... % 0.1 if x < 1
               0.2*rho_max*ones(1,length(2+dx:dx:3)), ... % 0.1 if x < 1
               0.4*rho_max*ones(1,length(3+dx:dx:4)), ... % 0.1 if x < 1
               0.7*rho_max*ones(1,length(4+dx:dx:5)), ... % 0.1 if x < 1
           0.9*rho_max*ones(1,length((5+dx:dx:x_max)))];     % 0.9 if 1 <= x < 2          % 0.9 if x >= 5
rho_x_min_t = rho_x_t_min(1)*ones(length(t),1);
rho_x_max_t = rho_x_t_min(length(x))*ones(length(t),1);

rho(:,1) = rho_x_t_min;
rho(1,:) = rho_x_min_t;
rho(length(x),:) = rho_x_max_t;