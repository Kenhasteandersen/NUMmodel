%% Decide which parameters to use in the sensitivity study
nrRandIter=6;
theParam={'bL','bN','bDOC','bF','bg'};
% randParam=rand(nrRandIter,length(theParam));
% randParam=repmat(linspace(0.1,0.8,nrRandIter),length(theParam),1)';
vector=linspace(0.02,0.5,nrRandIter);
vectors={vector,vector,vector,vector,vector};
n = numel(vectors); %// number of vectors
combs = cell(1,n); %// pre-define to generate comma-separated list
[combs{end:-1:1}] = ndgrid(vectors{end:-1:1}); %// the reverse order in these two
%// comma-separated lists is needed to produce the rows of the result matrix in
%// lexicographical order 
combs = cat(n+1, combs{:}); %// concat the n n-dim arrays along dimension n+1
combs = reshape(combs,[],n); %// reshape to obtain desired matrix
randParam=combs;

for itnr=1:length(randParam)%nrRandIter
    copyfile ../input/input_orig.nlm ../input/input.nlm

        system(['C:\cygwin64\bin\sed.exe "s/.*bL.*/      bL = ',num2str(randParam(itnr,1)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*bN.*/      bN = ',num2str(randParam(itnr,2)),'/" ../input/tmp1.nlm > ../input/tmp.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*bDOC.*/      bDOC = ',num2str(randParam(itnr,3)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*bF.*/      bF = ',num2str(randParam(itnr,4)),'/" ../input/tmp1.nlm > ../input/tmp.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*bg.*/      bg = ',num2str(randParam(itnr,5)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);

        movefile ../input/tmp1.nlm ../input/input.nlm

% end
    
    p = setupGeneralistsOnly;
    p = parametersChemostat(p);
    p.tEnd = 2000;
    p.d = 0.001;   
    p.rand =randParam(itnr,:);
    tic
    sim = simulateChemostat(p, 100);
    toc
    time = datestr(clock,'YYYY_mm_dd_HH_MM_SS');
    save(['outputsim_',time,'.mat'],'sim');
end
