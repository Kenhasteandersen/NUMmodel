function sim = simulateChemostat(p, T, L)

u = [p.N0, p.DOC0, p.B0];

[t,u] = ode23(@deriv, [0 365], u);
sim.t = t;
sim.N = u(:,1);
sim.DOC = u(:,2);
sim.B = u(:,3:end);


    function dudt = deriv(t,u)
        dudt = zeros(p.n+2,1);
        [u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcrates',T, L, 2+p.n, u, 1.0, 1.0, dudt);
        dudt(1) = dudt(1) + p.d*(150-u(1));
        dudt(2) = dudt(2) - p.d*u(2);
        dudt(3:end) = dudt(3:end) - p.d*u(3:end);
    end

end