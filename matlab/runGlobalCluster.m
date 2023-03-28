%
% This should be run on the cluster. Execute via submitGlobalCluster.m
%
%c = parcluster('dcc R2019a');
%c.AdditionalProperties.EmailAddress = 'kha@aqua.dtu.dk';
%c.AdditionalProperties.MemUsage = '64GB';
%c.AdditionalProperties.ProcsPerNode = 0;
%c.AdditionalProperties.WallTime = '4:00';
%c.saveProfile
%
clust = parcluster('dcc R2019a');
numW=8;    % Exactly the number of nodes times the number of processors per cores requested
parpool(clust, numW);

%load('tmpparameters');
%p = parametersGlobal(parameters([]),2);
p = parametersGlobal( setupNUMmodel(bParallel=true) );

sim = simulateGlobal(p);

save('tmp','sim','-v7.3');
