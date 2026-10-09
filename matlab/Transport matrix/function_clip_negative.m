%
% Removes negative concentrations after a transport step with the raw
% transport matrices. The mass missing in a negative box is taken from the
% nearest boxes above and below in the same water column: first from the
% boxes one step away, then two steps away, etc. The boxes at a given
% distance give in proportion to their mass. If the whole water column
% cannot cover the deficit, the remainder is added as new mass (the box is
% left at zero). No mass is moved between water columns, and each state
% variable is treated separately.
%
% Within each water column the negative boxes are corrected from the
% surface down. Positive boxes are only ever scaled down, so no box ends
% up negative.
%
% In:
%  u    - state (nb x n concentrations)
%  dv   - volume of each box (nb x 1, m3)
%  cols - the water columns (made in simulateGlobal), with fields:
%         order: all boxes, column by column, each from the surface down
%         pos:   position of each box in order
%         first, last: position in order of the top and bottom box of
%                each box's water column
%
% Out:
%  u  - state without negative values
%  st - 3 x n statistics for each state variable (concentration * m3):
%       st(1,:): negative mass removed
%       st(2,:): mass added because the water column could not cover it
%       st(3,:): mass taken from boxes more than one step away
%
function [u, st] = function_clip_negative(u, dv, cols)

st = zeros(3, size(u,2));
[ib, jt] = find(u < 0);
if isempty(ib)
    return
end
%
% Sort the negative entries by state variable and from the surface down:
%
[~, s] = sortrows([jt, cols.pos(ib)]);
ib = ib(s);
jt = jt(s);

for m = 1:numel(ib)
    b = ib(m);
    j = jt(m);
    d = -u(b,j)*dv(b); % The deficit
    st(1,j) = st(1,j) + d;
    u(b,j) = 0;

    q = cols.pos(b);
    q0 = cols.first(b);
    q1 = cols.last(b);
    r = 1;
    while d > 0 && (q-r >= q0 || q+r <= q1)
        % The positive boxes r steps above and below:
        bAbove = 0;
        bBelow = 0;
        avail = 0;
        if q-r >= q0
            bAbove = cols.order(q-r);
            if u(bAbove,j) > 0
                avail = avail + u(bAbove,j)*dv(bAbove);
            else
                bAbove = 0;
            end
        end
        if q+r <= q1
            bBelow = cols.order(q+r);
            if u(bBelow,j) > 0
                avail = avail + u(bBelow,j)*dv(bBelow);
            else
                bBelow = 0;
            end
        end
        % Take the deficit in proportion to their mass:
        if avail > 0
            take = min(d, avail);
            f = 1 - take/avail;
            if bAbove > 0
                u(bAbove,j) = u(bAbove,j)*f;
            end
            if bBelow > 0
                u(bBelow,j) = u(bBelow,j)*f;
            end
            if r > 1
                st(3,j) = st(3,j) + take;
            end
            d = d - take;
        end
        r = r + 1;
    end
    st(2,j) = st(2,j) + d; % What could not be covered is added
end
