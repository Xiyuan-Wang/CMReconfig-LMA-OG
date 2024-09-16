clear
close all

% canonical cross coupling matrix M_can and weight matrix W
load cm8.mat

% set LM options
opt.lambda = 0.1;
opt.max_iter = 100;
opt.tolopt = -inf;
opt.verbose = true;
opt.lossless = true;
opt.rand_init = true;

%% Perform coupling matrix reconfiguration
[M, obj_val, Q] = leven_marq(M_can, W, opt);



