function sim = simulateGlobalCluster(p,sim)
%
% Set up cluster:
%
c = parcluster('dcc R2019a');
c.AdditionalProperties.EmailAddress = 'kha@aqua.dtu.dk';
c.AdditionalProperties.MemUsage = '15GB';
c.AdditionalProperties.ProcsPerNode = 0;
c.AdditionalProperties.WallTime = '4:00';
c.saveProfile

%clust = parcluster('dcc R2019a');
numW=10;    % Exactly the number of nodes times the number of processors per cores requested
parpool(c, numW);

sim = simulateGlobal(p,sim);

%
% 80 cores: 11 minutes
% 20 cores: 4:01
% 10 cores: 4:18 min
% 5 cores:  5:24
