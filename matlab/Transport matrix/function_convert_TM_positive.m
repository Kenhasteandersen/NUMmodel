%
% Makes the off-diagonal elements of a transport matrix non-negative by
% moving each negative element to its transpose position. The diagonal is
% then set such that the volume-weighted column sums are unchanged:
%
%   dv' * A3 == dv' * A
%
% This is the condition for mass conservation when the matrix acts on
% concentrations in boxes of different volume (mass = dv'*c). Correcting the
% unweighted column sums instead is only conservative for boxes of equal
% volume and makes the global model leak tracer.
%
% In:
%  A  - transport matrix (rate matrix or propagator; both work since the
%       weighted column sums are preserved, whatever they are).
%  dv - volume of each box (m3), in the same ordering as A.
%
% Out:
%  A3 - matrix with non-negative off-diagonal elements.
%
function A3 = function_convert_TM_positive(A, dv)

n = size(A,1);
w = reshape(dv, 1, n);

N = A - spdiags(diag(A), 0, n, n); % Off-diagonal part
N(N>0) = 0;                        % ... of which only the negative elements
A2 = A - N - N';                   % Move them to the transpose position
A2 = A2 - spdiags(diag(A2), 0, n, n);

DA = (w*A - w*A2) ./ w;            % Diagonal that keeps dv'*A3 == dv'*A
A3 = A2 + spdiags(DA(:), 0, n, n);
