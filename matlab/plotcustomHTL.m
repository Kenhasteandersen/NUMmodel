%creates a bifurcation diagramm for custom HTL mortalitty values
function plotcustomHTL(nmort,maxmort,mAdult,fcr,opt)

arguments
    nmort double = 10 %mortality values tested (linspace 0-maxmort)
    maxmort double = 0.4 %maximum mortality value tested day^-1
    mAdult double = [20 1000] %copepod max sizes μg C
    fcr double = 10 %times smaller than smallest zoo adult size, opt=2
    opt int = 2 %type of mortality, 1:constant, 2:pHTL, 3: (coming soon)
end

mort=linspace(0,maxmort,nmort);
Bzoo=zeros(nmort,length(mAdult));
p = setupGeneric(mAdult);
p = parametersChemostat(p);
mmort=min(mAdult./fcr);
for i=1:nmort
   sim= testHTLmort(mort(i),mAdult,opt,mmort);
  for ii=1:length(mAdult)
      Bin = floor(0.8*length(sim.u));
      yend = mean(sim.u(Bin:end,p.ixStart(ii+1):p.ixEnd(ii+1)));
      Bzoo(i,ii)=sum(yend);
  end

end

semilogy(mort,Bzoo,'LineWidth',2)
xlabel('HTL mortality (day^{-1})')
ylabel('Biomass (μg l^{-1})')