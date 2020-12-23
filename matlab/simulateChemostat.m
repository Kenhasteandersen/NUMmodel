function sim = simulateChemostat(p, L)
%
% Find ix for nutrients and unicellulars:
%
ix = [1,2]; % Nutrients and DOC
for i = 1:p.nGroups
    if (p.typeGroups(i)==1)
        ix = [ix, p.ixStart(i):p.ixEnd(i)];
    end
end
%
% Concentrations in the deep layer:
%
uDeep = ix*0;
uDeep(1:2) = p.u0(1:2);
%
% Simulate:
%
[t,u] = ode23(@deriv, [0 365], p.u0);
%
% Assemble result:
%
sim.t = t;
sim.N = u(:,1);
sim.DOC = u(:,2);
sim.B = u(:,3:end);
sim.p = p;
sim.rates = calcDerivatives(p,u(end,:),L);
for iGroup = 1:p.nGroups
    sim.Bgroup(:,iGroup) = sum( u(:, p.ixStart(iGroup):p.ixEnd(iGroup)),2);
end

    %
    % Function to assemble derivative for chemostat:
    %
    function dudt = deriv(t,u)
        rates = calcDerivatives(p,u',L);
        dudt = rates.dudt;
        %
        % Chemostat dynamics for nutrients and unicellulars:
        %
        dudt(ix) = dudt(ix) + p.d*(uDeep-u(ix)');
        dudt = dudt';
    end
    
    %function dudt = deriv(t,u)
    %    dudt = zeros(p.n+2,1);
    %    [u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcrates',T, L, 2+p.n, u, 1.0, 1.0, dudt);
    %    dudt(1) = dudt(1) + p.d*(150-u(1));
    %    dudt(2) = dudt(2) - p.d*u(2);
    %    dudt(3:end) = dudt(3:end) - p.d*u(3:end);
    %end
    
    

end