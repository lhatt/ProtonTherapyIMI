% file init_data.m
% in this file all problem data and adaptivity parameter are initialized
% this file is called from afem.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  problem  data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we declare a "struct" for storing the problem data
global prob_data

% folder or directory where the domain mesh is described
domain = 'square_all_dirichlet';

% initial global refinements
global_refinements = 3;

prob_data.f = inline('f_poly(x)','x ');

% Dirichlet data, function g_D
prob_data.gD = inline('u_poly(x)','x ');

% Neumann data, function g_N
prob_data.gN = inline('0', 'x');

%prob_data.max_iterations = 5;

prob_data.grd_u_exact = inline('grdu_poly(x)','x ');

prob_data.epsilon = 0.01;
prob_data.advection = [-1,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  data for a posteriori estimators and adaptive strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We declare a data structure for storing all 
% the adaptive strategy parameters
global adapt

% weight in front of interior residual
adapt.C(1) = 1.0;
% weight in front of jump residual
adapt.C(2) = 1.0;

% tolerance for the adaptive strategy
adapt.tolerance = 1e-5;

% maximum number of iterations of adaptive strategy
adapt.max_iterations = 5;

% marking_strategy, possible options are
% GR: global (uniform) refinement,  
% MS: maximum strategy,  
% GERS: guaranteed error reduction strategy (D\"orfler's)
% ES: equidistribution strategy,           (not implemented yet) 
% MES: modified equidistribution strategy, (not implemented yet)
adapt.strategy = 'GR';

% n_refine, number of refinements of each marked element
adapt.n_refine = 2;

% parameters of the different marking strategies
% MS: Maximum strategy
adapt.MS_gamma = 0.75;

% GERS: guaranteed error reduction strategy (D\"orfler's)
adapt.GERS_theta_star = 0.8;
adapt.GERS_nu = 0.1;