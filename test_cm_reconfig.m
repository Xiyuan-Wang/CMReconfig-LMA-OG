% Clear workspace and close all figures
clear
close all

% Specify the coupling matrix file and algorithm to use
cm_file = 'cm/cm4'; % Path to the coupling matrix file
% Algorithm choice: 
% 'gen_iso_flow' or 'leven_marq' for LM algorithm
algo = 'gen_iso_flow'; 

%%
% Load the canonical cross coupling matrix M0 and weight matrix W
load(cm_file)
cm_reduce = str2func(algo); % Convert algorithm name to function handle

% Set algorithm options
opt.lambda = .1; % Regularization parameter
opt.alpha = 0.5; % Step size parameter
opt.max_iter = 50000; % Maximum number of iterations
opt.tolopt = -inf; % Tolerance for optimization
opt.tolfun = 1e-6; % Tolerance for function value change
% opt.timing_interval = 5 % Uncomment to set timing interval
opt.verbose = true; % Display iteration details
opt.lossless = true; % Ensure lossless reconfiguration
opt.rand_init = true; % Use random initialization; set to false for identity matrix initialization

%% Perform coupling matrix reconfiguration
% Execute the reconfiguration algorithm with specified options
[M, obj_val, Q, timing] = cm_reduce(M0, W, opt);



