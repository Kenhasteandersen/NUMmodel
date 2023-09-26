% Converts generalists and diatoms mass to radius in mum
% Converts copepod body-mass to prosome length/2 in mum, based on
% prosome length(um) to body-mass relationship for copepods
% from Chisholm and Roff (1990)
% -------treat POM as generalists------- MUST BE ADDRESSED
function r = calcRadiusGroups(p)

rho = 0.4*1e6*1e-12; % mug/cm3 (Andersen et al 2016; rho = m/V = 0.3*d^3/(4/3*pi*(d/2)^3) )

for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    m(ix) = p.m(ix+p.idxB-1);
    % ix = (p.ixStart(iGroup):p.ixEnd(iGroup));
    % m = p.m(ix);

    %
    % Generalists:
    %
    if ((p.typeGroups(iGroup)==1) || (p.typeGroups(iGroup)==5)|| (p.typeGroups(iGroup)==100))
        r(ix) = (3/(4*pi)*m(ix)/rho).^(1/3);
    end
    %
    % Diatoms:
    %
    if ((p.typeGroups(iGroup)==3) || (p.typeGroups(iGroup)==4))
        v = 0.6;
        r(ix) = (3/(4*pi)*m(ix)/rho/(1-v)).^(1/3); 
    end
    if ((p.typeGroups(iGroup)==10) || (p.typeGroups(iGroup)==11))
        % prosome length divided by 2 to be equivalent to radius
        r(ix)= ( m(ix)/(0.73*0.48*exp(-16.41))^(1/2.74) )/2; 
    end

end

