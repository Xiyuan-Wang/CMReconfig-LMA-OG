clear % Clear all variables from the workspace
close all % Close all open figure windows

cm_file = 'cm/cm4'; % Path to the coupling matrix file
algo = 'gen_iso_flow'; % Algorithm to use, can be 'leven_marq' or 'gen_iso_flow'

%%
load(cm_file) % Load the configuration matrix from the specified file
cm_reduce = str2func(algo); % Convert the algorithm name to a function handle

%% Set optimization options
opt.max_iter = 1e3; % Set the maximum number of iterations
opt.alpha = 0.5; % Set alpha parameter for the algorithm
opt.tolopt = -inf; % Set the tolerance for optimization
opt.tolfun = -inf; % Set the tolerance for function value
opt.verbose = false; % Disable verbose output
opt.lossless = true; % Enable lossless optimization

%% Set timing options
opt.timing_interval = 1; % Interval for timing measurements
opt.rand_init = true; % Enable random initialization

%% Initialize test parameters
N_test = 10; % Number of tests to run
n_timing = floor(opt.max_iter / opt.timing_interval); % Number of timing intervals
timings = zeros(n_timing, N_test); % Preallocate timing results
obj_vals = zeros(opt.max_iter+1, N_test); % Preallocate objective values
timing_iters = [1; opt.timing_interval * (1:n_timing)' + 1]; % Timing iteration indices

% Use parfor to accelerate the loop if possible
for n_test = 1:N_test
    fprintf("No. of Test = %d\n", n_test); % Display the current test number
    [~, obj_val, Q, timing] = cm_reduce(M0, W, opt); % Run the reduction algorithm
    obj_vals(:,n_test) = obj_val; % Store objective values
    timings(:,n_test) = timing(:,2); % Store timing results
end

%% Adjust timings for plotting
timings = [zeros(1,N_test); timings]; % Add zero timing for the initial point

%% Save results to data directory
% save('-v6', ['data/reconfig_converge_', algo, '_', cm_file, '.mat'], ...
%     'timings','obj_vals','timing_iters')

%% Plot results
blue_str = '#0095EF'; % Define a color for plotting
blue = sscanf(blue_str(2:end),'%2x%2x%2x',[1 3])/255; % Convert color to RGB
% Plot objective vs. time
loglog(timings, obj_vals(timing_iters,:), 'color', blue, ... 
       'linewidth', 1)
xlim([1e-3, max(timings(:))]) % Set x-axis limits
ylim([min(obj_vals(:)),10]) % Set y-axis limits
grid on % Enable grid
xlabel('Time (s)') % Label x-axis
ylabel('Objective') % Label y-axis
figure % Create a new figure
% Plot objective vs. iterations
loglog(1:size(obj_vals,1), obj_vals, 'color', blue, ... 
	 'linewidth', 1)
ylim([min(obj_vals(:)),10]) % Set y-axis limits
grid on % Enable grid
xlabel('No. of Iterations') % Label x-axis
ylabel('Objective') % Label y-axis
