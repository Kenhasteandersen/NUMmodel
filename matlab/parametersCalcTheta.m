%
% Calculate the size-based feeding preferences theta(predator,prey).
%
function theta = parametersCalcTheta(m,mLower,mDelta,beta,sigma)
n = length(m);

theta = zeros(n,n);
for i = 1:n % Predators
    %p.theta[i,] = phi(p.m[i]/p.m, p.beta, p.sigma)   % (predator, prey)
    for j = 1:n % Prey
        % THIS NEEDS TO BE UPDATED BECAUSE M(2)/M(1) IS NOT CONSTANT
        %theta(i,j) = Phi(m(i)/m(j), 1+mDelta/mLower, beta(i), sigma(i));
        
        % Just use the basic size-selection function for now. Needs to 
        % be updated to correct for bin widts.
        theta(i,j) = exp( -(log(m(i)/m(j)/beta(i)))^2/(2*sigma(i)^2));
    end
end

    %
    % Prey size function integrated over size groups:
    %
    function res = Phi(z, Delta,beta,sigma)
        m = logspace(-12,3,1000);
        dm = diff(m);
        m = m(2:1000);
        
        phi = @(z, beta, sigma) exp( -(log(z/beta)).^2/(2*sigma^2) );
        
        fPrey = @(m, w0, Delta) ...
            integral( @(logw) phi(m./exp(logw), beta, sigma), ...
            log(w0/sqrt(Delta)), log(w0*sqrt(Delta)));
        
        function int = fTot(m0,w0, Delta)
            ix = find((m>m0/sqrt(Delta)) & (m<m0*sqrt(Delta)));
            int = 0;
            for k = 1:length(ix)
                int = int + fPrey(m(ix(k)), w0, Delta)/m(ix(k))*dm(ix(k)) / log(Delta)^2;
            end
        end
        
        res = fTot(1,1/z,Delta);
    end

end