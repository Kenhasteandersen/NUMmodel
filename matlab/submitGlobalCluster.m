%
% Submit to DTU cluster. The results can be picked up 
% via a call to loadGlobalCluster
%
function submitGlobalCluster()
%
% Remove old log files:
%system('rm Error_*txt Output_*txt');
% Get parameters:
%save('tmpparameters','p');
% Submit run:
system('bsub < clusterrun.sh');


%
% 80 cores: 11 minutes
% 20 cores: 4:01
% 10 cores: 4:18 min
% 5 cores:  5:24
