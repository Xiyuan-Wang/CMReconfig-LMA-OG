function [Mr, varargout] = leven_marq(M, W, varargin)
    % LEVEN_MARQ Perform Levenberg-Marquardt optimization on orthogonal
    % matrices for coupling matrix reduction.
    %
    %   [Mr, varargout] = LEVEN_MARQ(M, W, varargin) performs the
    %   Levenberg-Marquardt optimization on orthogonal matrices for
    %   reconfiguration of the input coupling matrix M with structure
    %   specified by the binary weight matrix W.
    %   Additional options can be provided via varargin.
    %
    %   Inputs:
    %       M - Input coupling matrix.
    %       W - Binary weight matrix whose nonzero elements mark the couplings to be reduced.
    %       varargin - Additional options as a struct.
    %
    %   Outputs:
    %       Mr - reduced coupling matrix.
    %       varargout - Additional outputs including objective values,
    %                   transformation matrix Q, and timing information.
    
    % Set default options
    option.lambda = .1; % damping parameter for the Levenberg-Marquardt method
    option.max_iter = 100; % maximum number of iterations
    option.tolfun = 1e-3; % tolerance for objective function value
    option.tolopt = 1e-3; % tolerance for optimality (first order condition) 
    option.tolstep = 1e-3; % tolerance for step
    option.verbose = false; % set true to print objective and gradient values during iteration 
    option.plot = false; % set true to plot the objective function values
    option.timing_interval = inf; % set to a positive integer for recording timing information
    option.retraction = 'exp'; % retraction method for the orthogonal matrix
    option.rand_init = false; % set true to initialize the transformation matrix with a random orthogonal matrix

    % Override default options with user-provided options
    if nargin > 2
        user_option = varargin{1};
        for fn = fieldnames(user_option)'
            option.(fn{1}) = user_option.(fn{1});
        end
    end

    % Define retraction methods, the default is the matrix exponential
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
    n = size(M, 1);
    U0 = eye(n-2);
    if option.rand_init
        U0 = rand_corth_mat(n-2, option.lossless);
    end

    lambda = option.lambda;

    % Initialization
    U = U0;
    obj_values = zeros(1, option.max_iter + 1);
    obj_values(1) = obj(U, M, W);
    [~, grad_norm] = lm_step(U, lambda, M, W);    

    % Print header if verbose
    if option.verbose
        fprintf('\n')
        fprintf('-----------------------------------------------------\n')
        fprintf(' Iter | Obj Value | Grad Norm | Step Norm | Sol Norm\n')
        fprintf('-----------------------------------------------------\n')
    end

    start_time = tic;
    timing = [];
    iter = 1;
    if ~isinf(option.timing_interval)
        tic;
    end

    % Levenberg-Marquardt iteration
    Lu = 11;
    Ld = 9;
    while (iter <= option.max_iter) && (grad_norm > option.tolopt) ...
	  && (obj_values(iter) > option.tolfun)
        [G, grad_norm] = lm_step(U, lambda, M, W);    
        U_new = retr(U, G);    
        obj_value_new = obj(U_new, M, W);

        % Update lambda and step norm based on objective value
        if obj_value_new >= obj_values(iter)
            lambda = min(Lu * lambda, 1e7);
            step_norm = 0;
            obj_values(iter+1) = obj_values(iter);    
        else
            lambda = max(lambda / Ld, 1e-7);
            step_norm = norm(U - U_new, 'fro');
            U = U_new;
            obj_values(iter+1) = obj_value_new;    
        end

        % Print iteration details if verbose
        if option.verbose
            fprintf(' %4d | %-9.3e | %-9.3e | %-9.3e | %-9.3e\n', ...
                    iter, obj_values(iter+1), grad_norm, step_norm, norm(U, 'fro'))
        end

        % Record timing if required
        if ~isinf(option.timing_interval) && (mod(iter, option.timing_interval)==0)
            timing = [timing;iter toc];
            tic;
        end

        iter = iter + 1;
    end

    run_time = toc(start_time);

    % Print footer if verbose
    if option.verbose
        fprintf('-----------------------------------------------------\n')
    end

    % Plot objective function values if required
    if option.plot
        figure
        semilogy(1:option.max_iter+1, obj_values)
        grid on
        xlabel('No. of iterations')
        ylabel('Objective function value')
    end

    % Compute final result
    Q = blkdiag(1, U, 1);
    Mr = Q.' * M * Q;

    % Return additional outputs if requested
    if nargout > 1
        varargout{1} = obj_values';
    end

    if nargout > 2 
        varargout{2} = Q;
    end

    if nargout > 3
      if ~isempty(timing)
        varargout{3} = [timing(:,1) cumsum(timing(:,2))];
      else
	varargout{3} = [iter run_time];
      end      
    end    
end

% Objective function
function obj_val = obj(U, M, W)
    % OBJ Compute the objective function value.
    %
    %   obj_val = OBJ(U, M, W) computes the objective function value for
    %   the given transformation matrix U, input matrix M, and binary weight
    %   matrix W.
    %
    %   Inputs:
    %       U - Transformation matrix.
    %       M - Input coupling matrix.
    %       W - Binary weight matrix.
    %
    %   Outputs:
    %       obj_val - Objective function value.
    
    Q = blkdiag(1, U, 1);
    obj_val = norm(W .* (Q.' * M * Q), 'fro')^2;
end

% Levenberg-Marquardt step
function [G, grad_norm] = lm_step(U, lambda, M0, W)
    % LM_STEP Perform a single Levenberg-Marquardt step.
    %
    %   [G, grad_norm] = LM_STEP(U, lambda, M0, W) computes the gradient
    %   and step direction for the given transformation matrix U, damping
    %   parameter lambda, input coupling matrix M0, and binary weight matrix W.
    %
    %   Inputs:
    %       U - Transformation matrix.
    %       lambda - Damping parameter.
    %       M0 - Input coupling matrix.
    %       W - Binary weight matrix.
    %
    %   Outputs:
    %       G - Gradient matrix.
    %       grad_norm - Norm of the gradient.
    
    n = size(M0, 1);
    m = n - 2;

    Q = blkdiag(0, U, 0);
    D = diag([1, zeros(1,m), 1]);
    E = [zeros(1,m); eye(m); zeros(1,m)];
    
    M = (Q.' + D) * M0 * (Q + D);
    A = Q.' * M0 * Q + D * M0 * Q;
    
    J = W(:) .* (kron(E, A * E) - kron(A * E, E));

    indicator_lower = tril(ones(m), -1);
    [cols, rows] = find(indicator_lower);
    id_lo = sub2ind([m,m], cols, rows);
    id_up = sub2ind([m,m], rows, cols);
    
    Jc = J(:,id_lo) - J(:,id_up);

    % Compute gradient and update step
    grad = Jc' * M(:);
		% g = pinv(Jc' * Jc + lambda * eye((m^2-m)/2)) * grad;
    g = (Jc' * Jc + lambda * eye((m^2-m)/2)) \ grad;

    G = zeros(m);
    G(id_lo) = g;
    G(id_up) = -g;
    
    grad_norm = norm(grad);
end
