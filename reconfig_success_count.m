clear % Clear all variables from the workspace
close all % Close all open figure windows

cm_file = 'cm/cm4'; % path to the coupling matrix file
algo = 'gen_iso_flow'; % Algorithm to use, can be "leven_marq" or "gen_iso_flow"

%%
load(cm_file) % Load the configuration matrix from the specified file
cm_reduce = str2func(algo); % Convert the algorithm name to a function handle

% Set optimization options
opt.alpha = 0.5; % Set the alpha parameter for the algorithm
opt.max_iter = 50000; % Set the maximum number of iterations
opt.tolopt = -inf; % Set the tolerance for optimization
opt.tolfun = 1e-6; % Set the tolerance for function value
opt.verbose = false; % Disable verbose output
opt.lossless = true; % Enable lossless option

opt.rand_init = true; % Enable random initialization

%%

N_test = 1000; % Set the number of tests to run
run_times = inf(N_test, 2); % Initialize an array to store run times

%% use parfor to accelerate if possible
for n_test = 1:N_test % Loop over the number of tests
    fprintf("No. of Test = %d\n", n_test); % Print the current test number
    [~, obj_val, Q, run_time] = cm_reduce(M0, W, opt); % Run the algorithm and capture outputs
    
    if obj_val(run_time(1)) <= opt.tolfun % Check if the objective value meets the tolerance
      run_times(n_test,:) =  run_time; % Store the run time if condition is met
    end    
end

%% Save the results to the data directory
save('-v6', ['data/reconfig_success_count_', algo, '_', cm_file, '.mat'], ...
     'N_test', 'run_times') 

