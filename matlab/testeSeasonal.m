clear all;
sim=testMyChemostatSeasonal;

for i=1:length(sim.t);jDOC(i,:) = cell2mat({sim.rates(i).jDOC})';end
for i=1:length(sim.t);jLreal(i,:) = cell2mat({sim.rates(i).jLreal})';end
for i=1:length(sim.t);jFreal(i,:) = cell2mat({sim.rates(i).jFreal})';end

figure;
subplot(4,1,1)
plot(sim.t,sum(sim.B,2))


subplot(4,1,2)
plot(sim.t,sum(sim.B.*(jFreal./(jFreal+jLreal+jDOC)),2));subtitle('phagotrophy')
subplot(4,1,3)
plot(sim.t,sum(sim.B.*(jLreal./(jFreal+jLreal+jDOC)),2));subtitle('phototrophy')

subplot(4,1,4)
plot(sim.t,sum(jFreal,2),'-r')
hold on
plot(sim.t,sum(jLreal,2),'-g')


function sim = testMyChemostatSeasonal

% p = setupNUMmodel();
p = setupGeneralistsSimpleOnly();

p = parametersChemostat(p, 'lat_lon', [60 -15]);
p.tEnd=365*4;
sim = simulateChemostat(p, 'bUnicellularloss', false);
% 
% p_bis = parametersChemostat(p, 'seasonalAmplitude', 1);
% sim_bis = simulateChemostat(p_bis, 'bUnicellularloss', false);

% sumB = sum(sim.B(:));
% sumB_bis = sum(sim_bis.B(:));
% if ( sumB > 2693652  && sumB < 2693653 ) && ( sumB_bis > 33999 && sumB_bis < 34000 )
%     bSuccess = true;
% else
%     bSuccess = false;
%     fprintf(2,"sum(B) = %f; sum(B_bis) = %f\n", [sumB, sumB_bis]);
end