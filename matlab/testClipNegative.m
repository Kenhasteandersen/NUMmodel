%
% Tests function_clip_negative on a small hand-made case with three water
% columns (boxes 1-4, 5-8 and 9-11) and unequal box volumes. The value
% argument is not used; it is there to fit testFunction.
%
function bSuccess = testClipNegative(~)

addpath(strcat(fileparts(mfilename('fullpath')),'/Transport matrix'));

cols.order = (1:11)';
cols.pos = (1:11)';
cols.first = [1 1 1 1 5 5 5 5 9 9 9]';
cols.last = [4 4 4 4 8 8 8 8 11 11 11]';
dv = [1 2 1.5 1 1 3 2 1 2 1 1]';
%
% (a) column 1: one negative box with positive neighbours above and below
% (b) column 2: [-1 0.1 0.1 0.1]; the column cannot cover the deficit
% (c) variable 2, column 1: the worked example [-0.5 0.01 1 1]
% (d) variable 3: all positive
%
U = [ [0.3 -0.2 0.4 0.5,  -1.0 0.1 0.1 0.1,  1 1 1]', ...
      [-0.5 0.01 1.0 1.0,  1 1 1 1,  1 1 1]', ...
      ones(11,1) ];
[u, st] = function_clip_negative(U, dv, cols);

tol = 1e-12;
bSuccess = true;
% No negative values, and the mass changes only by what is added:
bSuccess = bSuccess && all(u(:) >= 0);
bSuccess = bSuccess && all(abs(dv'*u - dv'*U - st(2,:)) <= tol*(dv'*U));
% (a) the neighbours give the same fraction; the next box down is untouched:
bSuccess = bSuccess && abs(u(1,1)/U(1,1) - u(3,1)/U(3,1)) < tol && u(4,1) == U(4,1);
% (b) the whole column ends at 0 and nothing outside it changes:
bSuccess = bSuccess && all(u(5:8,1) == 0) && all(u(9:11,1) == U(9:11,1));
bSuccess = bSuccess && abs(st(2,1) - (dv(5) - 0.1*sum(dv(6:8)))) < tol;
% (c) the worked example:
take = 0.5*dv(1) - 0.01*dv(2);
bSuccess = bSuccess && all(u(1:2,2) == 0) && abs(u(3,2) - (1 - take/dv(3))) < tol ...
    && u(4,2) == 1 && st(2,2) == 0 && abs(st(3,2) - take) < tol;
% (d) unchanged:
bSuccess = bSuccess && isequal(u(:,3), U(:,3)) && all(st(:,3) == 0);

if ~bSuccess
    fprintf(2, 'function_clip_negative does not give the expected result\n');
end
