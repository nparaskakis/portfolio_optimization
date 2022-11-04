%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Advanced Topics in Convex Optimization %%
%% Nikolaos Paraskakis - A.M.: 2018030027 %%
%%       Portfolio Optimization           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%
%% Initialization %%
%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;
rng(6799);





%%%%%%%%%%%%%%%%%%%%%%
%% Add needed paths %%
%%%%%%%%%%%%%%%%%%%%%%

addpath('./Financial_Data');
addpath('./Projected_Subgradient_Method');
addpath('./Interior_Point_Method_Phase1');
addpath('./Interior_Point_Method_Phase2');
addpath('./FDPG_Method_For_Projection');





%%%%%%%%%%%%%%%%%%%%%%%
%% Load problem data %%
%%%%%%%%%%%%%%%%%%%%%%%

mu = load('meanReturns.mat');
mu = mu.x';

Sigma = load('covMat.mat');
Sigma = Sigma.x;

n = length(mu);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1: Optimize portfolio for minimum volatility %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('./Minimum_Volatility');

v1 = zeros(n,1);
v2 = ones(n,1);

s1 = struct('n',n, 'v1',v1, 'v2',v2, 'cost_func', @(w) p2_cost_func(w,Sigma));

fprintf('\nCVX for solving "1: Minimum volatility":\n');

tic;
[w_opt_1_cvx,f_opt_1_cvx] = p2_cvx(s1);
toc;

fprintf('f_opt_1_cvx: %g\n',f_opt_1_cvx);

print_val = 0;
c = 1.0001;
w_initial = rand(n,1);
Lf = max(abs(eig(Sigma)));

s2 = struct('print_val', print_val, 'c',c, 'Lf',Lf, 'f_opt_cvx',f_opt_1_cvx, 'x_initial', w_initial, 'cost_func',@(w) p2_cost_func(w,Sigma) , 'subgradient',@(w) p2_subgradient(w,Sigma), 'projection',@(d) p1_projection(d,n,v1,v2), 'stepsize',@(f_x_k,df_x_k,f_opt_cvx,k,Lf) dynamic_step(f_x_k,df_x_k,f_opt_cvx,k,Lf));

fprintf('\nProjected subgradient for solving "1: Minimum volatility":\n');

tic;
[w_opt_1_proj_subgrad,f_opt_1_proj_subgrad,~,~,iter_1_proj_subgrad] = projected_subgradient_method(s2);
toc;

fprintf('Iterations are %g.\n',iter_1_proj_subgrad);
fprintf('f_opt_1_proj_subgrad: %g\n',f_opt_1_proj_subgrad);

rmpath('./Minimum_Volatility');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2: Optimize portfolio by minimizing     %%
%%    volatility for a given target return %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('./Minimum_Volatility_On_Target_Return');

v1 = zeros(n,1);
v2 = ones(n,1);

ret = 0.15;

s3 = struct('n',n, 'mu',mu, 'ret',ret, 'v1',v1, 'v2',v2, 'cost_func', @(w) p2_cost_func(w,Sigma));

fprintf('\nCVX for solving "2: Minimum volatility for a given target return":\n');

tic;
[w_opt_2_cvx,f_opt_2_cvx] = p2_cvx(s3);
toc;

fprintf('f_opt_2_cvx: %g\n',f_opt_2_cvx);

print_val = 0;
c = 1.0001;
w_initial = rand(n,1);
Lf = max(abs(eig(Sigma)));

s4 = struct('print_val', print_val, 'c',c, 'Lf',Lf, 'f_opt_cvx',f_opt_2_cvx, 'x_initial', w_initial, 'cost_func',@(w) p2_cost_func(w,Sigma) , 'subgradient',@(w) p2_subgradient(w,Sigma), 'projection',@(d) p1_projection(d,n,mu,ret,v1,v2), 'stepsize',@(f_x_k,df_x_k,f_opt_cvx,k,Lf) dynamic_step(f_x_k,df_x_k,f_opt_cvx,k,Lf));

fprintf('\nProjected subgradient for solving "2: Minimum volatility for a given target return":\n');

tic;
[w_opt_2_proj_subgrad,f_opt_2_proj_subgrad,~,~,iter_2_proj_subgrad] = projected_subgradient_method(s4);
toc;

fprintf('Iterations are %g.\n',iter_2_proj_subgrad);
fprintf('f_opt_2_proj_subgrad: %g\n',f_opt_2_proj_subgrad);

rmpath('./Minimum_Volatility_On_Target_Return');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3: Optimize portfolio by maximizing utility %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('./Maximum_Utility');

v1 = zeros(n,1);
v2 = ones(n,1);

risk_aversion = 1;

s5 = struct('n',n, 'v1',v1, 'v2',v2, 'cost_func', @(w) p2_cost_func(w,mu,Sigma,risk_aversion));

fprintf('\nCVX for solving "3: Maximum utility":\n');

tic;
[w_opt_3_cvx,f_opt_3_cvx] = p2_cvx(s5);
toc;

fprintf('f_opt_3_cvx: %g\n',f_opt_3_cvx);

print_val = 0;
c = 1.00001;
w_initial = rand(n,1);
Lf = max(abs(eig(Sigma)));

offset = 0;

if (f_opt_3_cvx < 0)
    offset = abs(f_opt_3_cvx) + 1;
end

s6 = struct('print_val', print_val, 'c',c, 'Lf',Lf, 'f_opt_cvx',f_opt_3_cvx + offset, 'x_initial', w_initial, 'cost_func',@(w) p2_cost_func(w,mu,Sigma,risk_aversion) + offset , 'subgradient',@(w) p2_subgradient(w,mu,Sigma,risk_aversion), 'projection',@(d) p1_projection(d,n,v1,v2), 'stepsize',@(f_x_k,df_x_k,f_opt_cvx,k,Lf) dynamic_step(f_x_k,df_x_k,f_opt_cvx,k,Lf));

fprintf('\nProjected subgradient for solving "3: Maximum utility":\n');

tic;
[w_opt_3_proj_subgrad,f_opt_3_proj_subgrad,~,~,iter_3_proj_subgrad] = projected_subgradient_method(s6);
toc;

fprintf('Iterations are %g.\n',iter_3_proj_subgrad);
fprintf('f_opt_3_proj_subgrad: %g\n',f_opt_3_proj_subgrad - offset);

rmpath('./Maximum_Utility');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4: Optimize portfolio by maximizing     %%
%%    return for a given target volatility %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('./Maximum_Return_On_Target_Volatility');

v1 = zeros(n,1);
v2 = ones(n,1);

A = ones(1,n);
b = 1;

alpha = 0.05;
beta = 0.5;
mu_param = 1.5;
t_bm = rand;

p = 1;
m = 2*n + 1;

e = 10^(-4);

vol = 0.3;
var = vol^2;

s7 = struct('n',n, 'Sigma',Sigma, 'var',var, 'v1',v1, 'v2',v2, 'cost_func', @(w) f0_nbm_phase2(w,mu));

fprintf('\nCVX for solving "4: Maximum return for a given target volatility":\n');

tic;
[w_opt_4_cvx,f_opt_4_cvx] = cvx_phase2(s7);
toc;

fprintf('f_opt_4_cvx: %g\n',f_opt_4_cvx);

%w_0 = findFeasibleCVX(n,v1,v2,Sigma,var);

s8 = struct('n',n, 'p',p, 'm',m, 'A',A, 'b',b, 'alpha',alpha, 'beta',beta, 'mu',mu_param, 't_bm', t_bm, 'e',e, 'f0',@(w) f0_nbm_phase1(w), 'grad_f0',@(w) grad_f0_nbm_phase1(w), 'hessian_f0',@(w) hessian_f0_nbm_phase1(w), 'phi',@(w) phi_nbm_phase1(w,Sigma,var), 'grad_phi',@(w) grad_phi_nbm_phase1(w,Sigma,var), 'hessian_phi',@(w) hessian_phi_nbm_phase1(w,Sigma,var), 'is_feasible',@(w) is_feasible(w,A,b,Sigma,var), 'find_s_initial',@(w) find_s_initial(w,Sigma,var));
w_0 = interior_point_method_phase1(s8);

s9 = struct('n',n, 'p',p, 'm',m, 'A',A, 'b',b, 'alpha',alpha, 'beta',beta, 'mu',mu_param, 't_bm', t_bm, 'e',e, 'x_initial',w_0, 'f0',@(w) f0_nbm_phase2(w,mu), 'grad_f0',@(w) grad_f0_nbm_phase2(w,mu), 'hessian_f0',@(w) hessian_f0_nbm_phase2(w), 'phi',@(w) phi_nbm_phase2(w,Sigma,var), 'grad_phi',@(w) grad_phi_nbm_phase2(w,Sigma,var), 'hessian_phi',@(w) hessian_phi_nbm_phase2(w,Sigma,var), 'is_feasible',@(w) is_feasible(w,A,b,Sigma,var));

fprintf('\nInterior point for solving "4: Maximum return for a given target volatility":\n');

tic;
[w_opt_4_ip,f_opt_4_ip,~,~,iter_4_ip] = interior_point_method_phase2(s9);
toc;

fprintf('Iterations are %g.\n',iter_4_ip);
fprintf('f_opt_4_ip: %g\n',f_opt_4_ip);

rmpath('./Maximum_Return_On_Target_Volatility');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5: Optimize portfolio by   %%
%%    maximizing Sharpe ratio %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('./Maximum_Sharpe_Ratio');

v1 = zeros(n,1);
v2 = ones(n,1);

A = mu';
b = 1;

alpha = 0.05;
beta = 0.5;
mu_param = 1.5;
t_bm = rand;

p = 1;
m = 2*n;

e = 10^(-4);

C = (-1)*eye(n);
d = (-1)*ones(n,1);

C_hat = zeros(size(C));

for i = (1:1:n)
    for j = (1:1:n)
        C_hat(i,j) = C(i,j) - d(i);
    end
end

s10 = struct('n',n, 'mu',mu, 'C_hat',C_hat, 'cost_func', @(y) f0_nbm_phase2(y,Sigma));

fprintf('\nCVX for solving "5: Maximum Sharpe ratio":\n');

tic;
[w_opt_5_cvx,f_opt_5_cvx] = cvx_phase2(s10);
toc;

fprintf('f_opt_5_cvx: %g\n',f_opt_5_cvx);

%y_0 = findFeasibleCVX(n,A,b,C_hat);

s11 = struct('n',n, 'p',p, 'm',m, 'A',A, 'b',b, 'alpha',alpha, 'beta',beta, 'mu',mu_param, 't_bm', t_bm, 'e',e, 'f0',@(y) f0_nbm_phase1(y), 'grad_f0',@(y) grad_f0_nbm_phase1(y), 'hessian_f0',@(y) hessian_f0_nbm_phase1(y), 'phi',@(y) phi_nbm_phase1(y,C_hat), 'grad_phi',@(y) grad_phi_nbm_phase1(y,C_hat), 'hessian_phi',@(y) hessian_phi_nbm_phase1(y,C_hat), 'is_feasible',@(y) is_feasible(y,A,b,C_hat), 'find_s_initial',@(y) find_s_initial(y,C_hat));
y_0 = interior_point_method_phase1(s11);

s12 = struct('n',n, 'p',p, 'm',m, 'A',A, 'b',b, 'alpha',alpha, 'beta',beta, 'mu',mu_param, 't_bm', t_bm, 'e',e, 'x_initial',y_0, 'f0',@(y) f0_nbm_phase2(y,Sigma), 'grad_f0',@(y) grad_f0_nbm_phase2(y,Sigma), 'hessian_f0',@(y) hessian_f0_nbm_phase2(y,Sigma), 'phi',@(y) phi_nbm_phase2(y,C_hat), 'grad_phi',@(y) grad_phi_nbm_phase2(y,C_hat), 'hessian_phi',@(y) hessian_phi_nbm_phase2(y,C_hat), 'is_feasible',@(y) is_feasible(y,A,b,C_hat));

fprintf('\nInterior point for solving "5: Maximum Sharpe ratio":\n');

tic;
[w_opt_5_ip,f_opt_5_ip,~,~,iter_5_ip] = interior_point_method_phase2(s12);
toc;

w_opt_5_ip = w_opt_5_ip/sum(w_opt_5_ip);

fprintf('Iterations are %g.\n',iter_5_ip);
fprintf('f_opt_5_ip: %g\n',f_opt_5_ip);

rmpath('./Maximum_Sharpe_Ratio');





%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove added paths %%
%%%%%%%%%%%%%%%%%%%%%%%%

rmpath('./Financial_Data');
rmpath('./Projected_Subgradient_Method');
rmpath('./Interior_Point_Method_Phase1');
rmpath('./Interior_Point_Method_Phase2');
rmpath('./FDPG_Method_For_Projection');