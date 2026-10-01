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
% The cells are threaded inside the library with OpenMP, so no parallel pool is
% needed. Set OMP_NUM_THREADS to the number of cores requested from the queue
% before matlab starts.
%
%load('tmpparameters');
%p = parametersGlobal(parameters([]),2);
p = parametersGlobal( setupNUMmodel );

sim = simulateGlobal(p);

save('tmp','sim','-v7.3');
