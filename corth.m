function Q = corth(A)
% CORTH Performs real/complex orthogonalization of the column vectors of A using the Gram-Schmidt method.
%   This function takes a matrix A and returns a matrix Q such that:
%   - Q.' * Q = I, where I is the identity matrix.
%
%   Input:
%       A - A complex matrix with column vectors to be orthogonalized.
%
%   Output:
%       Q - A complex matrix with orthonormal column vectors.

n = size(A, 2);

Q = A;
Q(:,1) = Q(:,1) / sqrt(sum(Q(:,1).^2));
for k = 2 : n
    Q(:,k) = Q(:,k) - Q(:,1:k-1) * Q(:,1:k-1).' * Q(:,k);
    Q(:,k) = Q(:,k) / sqrt(sum(Q(:,k).^2));
end
