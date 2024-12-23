function [Mr, varargout] = gen_iso_flow(M, W, varargin)
% GEN_ISO_FLOW Performs reconfiguration of a coupling matrix M using general isospectral flow algorithm.

%   Inputs:
%       M - Input coupling matrix.
%       W - Binary weight matrix whose nonzero elements mark the couplings to be reduced.
%       varargin - Additional options as a struct.
%
%   Outputs:
%       Mr - reduced coupling matrix.
%       varargout - Additional outputs including objective values,
%                   transformation matrix Q, and timing information.

    % Set default options for the algorithm
    option.alpha = .3; % Backtracking parameter (Armijo rule)
    option.beta = .7; % Backtracking parameter
    option.sigma = .3; % Backtracking parameter
    option.max_iter = 100; % Maximum number of iterations
    option.tolfun = 1e-3; % Tolerance for objective function value
    option.tolopt = 1e-3; % Tolerance for optimality (first order condition)
    option.tolstep = 1e-3; % Tolerance for step size
    option.verbose = false; % Verbose output flag
    option.plot = false; % Plot objective function values
    option.timing_interval = inf; % Interval for timing output
    option.retraction = 'gs'; % Default retraction method (Gram-Schmidt)
    option.rand_init = false; % Random initialization flag

    % Override default options with user-specified options
    if nargin > 2
        user_option = varargin{1};
        for fn = fieldnames(user_option)'
            option.(fn{1}) = user_option.(fn{1});
        end
    end

    % Define retraction methods
    switch option.retraction
        case 'exp'
            retr = @(U, X) U*expm(-X);
        case 'cayley'
            retr = @(U, X) U * ((eye(size(X)) + X) \ (eye(size(X)) - X));
        case 'gs'
            retr = @(U, X) corth(U * (eye(size(X)) - X)); 
        case 'euclid'
            retr = @(U, X) U - X;
        otherwise
            retr = @(U, X) U*expm(-X);
    end

    % Initialize variables
    n = size(M, 1); % Size of the matrix M
    D = diag([1 zeros(1, n-2) 1]); % Diagonal matrix for boundary conditions

    % Define transformation and projection functions
    T = @(X) blkdiag(0, X, 0);
    Tt = @(X) X(2:size(X,1)-1,2:size(X,2)-1);
    Pw = @(X) W .* X; % Element-wise multiplication with weight matrix W
    obj = @(Q) norm(Pw((T(Q.') + D) * M * (T(Q) + D)), 'fro')^2; % Objective function

    % Initialize Q
    Q = eye(n-2);
    if option.rand_init
        Q = rand_corth_mat(n-2, option.lossless);
    end

    % Initialize objective values and gradient
    obj_values = zeros(1, option.max_iter + 1);
    obj_values(1) = obj(Q);

    Mq = (T(Q.') + D) * M * (T(Q) + D);
    A = (T(Q.') + D) * M * T(Q);
    G = Tt(A' * Pw(Mq) - Pw(Mq) * conj(A));
    grad_norm = norm(G, 'fro');

    % Verbose output header
    if option.verbose
        fprintf('\n')
        fprintf('------------------------------------------\n')
        fprintf(' Iter | Obj Value | Grad Norm | Step Norm \n')
        fprintf('------------------------------------------\n')
    end

    % Initialize timing and iteration variables
    start_time = tic;
    timing = [];
    iter = 1;
    t = option.alpha;

    if ~isinf(option.timing_interval)
        tic;
    end

    % Main optimization loop
    while (iter <= option.max_iter) && (grad_norm > option.tolopt) ...
          && (obj_values(iter) > option.tolfun)

        grad_norm = norm(G, 'fro');
        obj_diff = grad_norm^2;

        % Compute new Q using retraction
        Q_new = retr(Q, t * G);    
        obj_value_new = obj(Q_new);    
        obj_diff = t * option.sigma * obj_diff;    

        % Backtracking line search
        if obj_value_new + obj_diff > obj_values(iter)
            step_norm = 0;
            t = option.beta * t;
            obj_values(iter+1) = obj_values(iter);    
        else
            step_norm = norm(Q - Q_new, 'fro');
            Q = Q_new;

            Mq = (T(Q.') + D) * M * (T(Q) + D);
            A = (T(Q.') + D) * M * T(Q);
            G = Tt(A' * Pw(Mq) - Pw(Mq) * conj(A));

            obj_values(iter+1) = obj_value_new;
            t = option.alpha;
        end

        % Verbose output for each iteration
        if option.verbose
            fprintf(' %4d | %-9.3e | %-9.3e | %-9.3e\n', ...
                    iter, obj_value_new, grad_norm, step_norm)
        end

        % Timing output
        if ~isinf(option.timing_interval) && ...
            (mod(iter, option.timing_interval)==0)
            timing = [timing;iter toc];
            tic;
        end

        iter = iter + 1;
    end

    % Calculate total runtime
    run_time = toc(start_time);

    % Verbose output footer
    if option.verbose
        fprintf('-----------------------------------------------------\n')
    end

    % Plot objective function values if requested
    if option.plot
        figure
        semilogy(1:option.max_iter+1, obj_values)
        grid on
        xlabel('No. of iterations')
        ylabel('Objective function value')
    end

    % Compute final reconfigured matrix
    Mr = T(Q.') * M * T(Q);

    % Return additional outputs if requested
    if nargout > 1
        varargout{1} = obj_values';
    end

    if nargout > 2 
        varargout{2} = T(Q);
    end

    if nargout > 3
        if ~isempty(timing)
            varargout{3} = [timing(:,1) cumsum(timing(:,2))];
        else
            varargout{3} = [iter run_time];
        end
    end
end



