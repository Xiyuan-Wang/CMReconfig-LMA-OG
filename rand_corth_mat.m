function R = rand_corth_mat(n, varargin)
%% RAND_CORTH_MAT Generate a random complex or real orthogonal matrix
%   R = RAND_CORTH_MAT(n) generates an n-by-n random complex orthogonal matrix.
%   R = RAND_CORTH_MAT(n, true) generates an n-by-n random real orthogonal matrix.
%
%   Inputs:
%       n       - Size of the matrix (n-by-n)
%       varargin- Optional argument to specify if the matrix should be real
%                 (true for real, false for complex; default is false)
%
%   Outputs:
%       R       - n-by-n orthogonal matrix (real or complex)
%
%   The function uses QR decomposition to generate the orthogonal matrix.
%
%   Example:
%       R = rand_corth_mat(5);       % Generates a 5x5 complex orthogonal matrix
%       R = rand_corth_mat(5, true); % Generates a 5x5 real orthogonal matrix

REAL = false;
if nargin > 1
    REAL = varargin{1};
end

if REAL
    R = randn(n);
else
    R = randn(n) + 1i*randn(n);
end

R = corth(R);



