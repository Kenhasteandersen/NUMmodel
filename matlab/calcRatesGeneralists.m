function rates = calcRatesGeneralists(ixStart, ixEnd, u, rates, p, L)
ix = ixStart:ixEnd;
%B = pmax(0,u(ixStart:ixEnd));
N = max(0,u(1));
DOC = max(0, u(2));
%
% Temperature corrections:
%
ANmT = p.ANm*fTemp(1.5,p.T);
f2 = fTemp(2,p.T);
JmaxT = p.Jmax*f2;
JR = p.Jrespm*f2;
JFmaxmT = p.JFmaxm*f2;
%
% Uptakes
%
rates.JN(ix) =   ANmT*N*p.rhoCN; % Diffusive nutrient uptake in units of C/time
%
rates.JDOC(ix) = ANmT*DOC; % Diffusive DOC uptake, units of C/time

rates.JL(ix) =   p.epsilonL * p.ALm*L;  % Photoharvesting

% Light acclimation:
%JLreal = pmax( 0, JL - pmax(0,(JL+JDOC-JR - JN)) )

% Feeding as type II with surface limitation:
%F = theta %*% B
%JF = epsilonF * JFmaxmT * AFm*F / (AFm*F + JFmaxmT) %        % Feeding

% Passive losses:
%Jloss_passive = p.cLeakage * m^(2/3) % in units of C

% Total nitrogen uptake:
rates.JNtot(ix) = rates.JN(ix)+rates.JF(ix)-p.Jloss_passive_m; % In units of C
% Total carbon uptake
rates.JCtot(ix) = rates.JL(ix)+rates.JF(ix)+rates.JDOC(ix)-JR-p.Jloss_passive_m;

% Liebig + synthesis limitation:
%Jtot = pmin( JNtot, JCtot, JmaxT )
rates.Jtot(ix) = min( rates.JNtot(ix), rates.JCtot(ix) );
f = rates.Jtot(ix)./(rates.Jtot(ix) + JmaxT);

% If synthesis-limited then down-regulate feeding:
%JFreal = pmax(0, JF - pmax( 0,  pmax(0, Jtot-JmaxT) ))
rates.JFreal(ix) = max(zeros(1,length(ix)), rates.JF(ix) - ...
    (rates.Jtot(ix)-f.*JmaxT).*(rates.Jtot(ix)>0));
rates.Jtot(ix) = f.*JmaxT;
rates.JLreal(ix) = rates.JL(ix) - max(zeros(1,length(ix)), ...
    min((rates.JCtot(ix) - (rates.JF(ix)-rates.JFreal(ix))-rates.Jtot(ix)), rates.JL(ix)));

% Actual uptakes:
rates.JCtot(ix) = rates.JLreal(ix) + rates.JDOC(ix) + rates.JFreal(ix) - JR - ....
    p.Jloss_passive_m;
rates.JNtot(ix) = rates.JN(ix) + rates.JFreal(ix) - p.Jloss_passive_m;
%
% Losses:
%
rates.JCloss_feeding(ix) = (1-p.epsilonF)/p.epsilonF*rates.JFreal(ix); % Incomplete feeding (units of carbon per time)
rates.JCloss_photouptake(ix) = (1-p.epsilonL)/p.epsilonL*rates.JLreal(ix);
rates.JNlossLiebig(ix) = max(zeros(1,length(ix)), rates.JNtot(ix)-rates.Jtot(ix));  % In units of C
rates.JClossLiebig(ix) = max(zeros(1,length(ix)), rates.JCtot(ix)-rates.Jtot(ix)); % C losses from Liebig, not counting losses from photoharvesting

rates.JNloss(ix) = rates.JCloss_feeding(ix) + rates.JNlossLiebig(ix) + p.Jloss_passive_m; % In units of C
rates.JCloss(ix) = rates.JCloss_feeding(ix) + rates.JCloss_photouptake(ix) + ...
    rates.JClossLiebig(ix) + p.Jloss_passive_m;
