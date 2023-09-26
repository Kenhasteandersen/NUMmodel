%
% Returns a vector of strategies determined by the uptake rates
% Second return argument is colors corresponding to each strategy
% 
function [strategy, col] = calcTrophicStrategy(rates)

arguments
    rates struct
end

rhoCN = 5.68;

col = zeros(length(rates.jL),3);

for i = 1:length(rates.jL)
    strategy{i} = 'Unknown';
    
    if rates.jN(i)*rhoCN > rates.jL(i)
        strategy{i} = 'Light limited';
        col(i,:) = [0,1,0];
    else
        strategy{i} = 'Nutrient limited';
        col(i,:) = [0,0,1];
    end
    
    if rates.jDOC(i) > rates.jLreal(i)
        strategy{i} = 'Osmoheterotroph';
        col(i,:) = [0.5,0,0.5];
    end
    
    if (rates.jFreal(i)/rates.jL(i) > 0.25) || (rates.jF(i) > rates.jN(i)*rhoCN)
        strategy{i} = 'Mixotroph';
        col(i,:) = [1,0.5,0.5];
    end
    
    if (rates.jNloss(i) > 1e-5) && (rates.jN(i) < rates.jF(i)/rhoCN)
        strategy{i} = 'Heterotroph';
        col(i,:) = [1,0,0];
    end
      
end


% strategy = rep('Unknown', p$n)
%   strategy[r$jN*p$rhoCN>r$jL] = "Light limited"
%   strategy[r$jL>=r$jN*p$rhoCN] = "Nutrient limited"
%   strategy[r$jDOC > r$jLreal] = "Osmoheterotroph"
%   strategy[(r$jNloss>1e-5) & (r$JN<r$JF/p$rhoCN)] = "Heterotroph"
%   strategy[((r$jFreal/r$jL > 0.25) | (r$jF>r$jN*p$rhoCN) )& !strategy=="Heterotroph"] = "Mixotroph"  
%   return(strategy)
% }
