ixD=11:20;
mD=sim.p.m(14:23); % diatom mass
bL=0.08; bN=0.045; bg=0.02;
bJx_tot=bL*sim.rates.jLreal(ixD)+bN*(sim.rates.jN(ixD)+sim.rates.jSi(ixD)+...
    sim.rates.jDOC(ixD));
%%
L=linspace(100,500,6);
for i=1:length(L)
   for j=ixD 
    jR(i,j)=simL(i).rates.jR(j);
    jTot(i,j)=bg*simL(i).rates.jTot(j);
    jXtots(i,j)= bL*simL(i).rates.jLreal(j)+bN*(simL(i).rates.jN(j)+...
        simL(i).rates.jSi(j)+simL(i).rates.jDOC(j));
    BdL(i,j)=simL(i).B(j);
       end
end
%%
for i=1:length(d)
   for j=ixD 
    jRd(i,j)=simL(i).rates.jR(j);
    jTotd(i,j)=bg*simL(i).rates.jTot(j);
    jXtotsd(i,j)= bL*simL(i).rates.jLreal(j)+bN*(simL(i).rates.jN(j)+...
        simL(i).rates.jSi(j)+simL(i).rates.jDOC(j));
        Bdd(i,j)=simd(i).B(j);
   end
end

ix=[11,15,20];
%%
figure(7)
clf(7)
subplot(231)
% for i=1:length(L)
hold on
% set(gca,'LineWidth',3)
semilogx(L,jR(:,ix),'LineWidth',3)
semilogx(L,jTot(:,ix),'--','LineWidth',3)
semilogx(L,jXtots(:,ix),':','LineWidth',2)
xlabel('Light')
ylabel('rates (/day)')
axis tight

subplot(232)
hold on
% set(gca,'LineWidth',3)
% semilogx(L,jR,'LineWidth',3)
semilogx(L,jTot(:,11:20),'--','LineWidth',3)
% semilogx(L,jXtots,':','LineWidth',2)
% semilogx(L(i),bg*simL(i).rates.jTot(1),'bd')
% semilogx(L(i),bJx_tot(1),'gd')
xlabel('Light')
ylabel('rates (/day)')
axis tight


subplot(232)
hold on
% set(gca,'LineWidth',3)
% semilogx(L,jR,'LineWidth',3)
semilogx(L,jTot(:,11:20),'--','LineWidth',3)
% semilogx(L,jXtots,':','LineWidth',2)
% semilogx(L(i),bg*simL(i).rates.jTot(1),'bd')
% semilogx(L(i),bJx_tot(1),'gd')
xlabel('Light')
ylabel('rates (/day)')
axis tight
title('Jtot all sizes')

subplot(233)
% for i=1:length(L)
hold on
% set(gca,'LineWidth',3)
semilogx(L,BdL,'LineWidth',3)
xlabel('Light')
ylabel('rates (/day)')
title('biomass')
axis tight

subplot(234)
% for i=1:length(L)
hold on
% set(gca,'LineWidth',3)
semilogx(d,jRd(:,ix),'LineWidth',3)
semilogx(d,jTotd(:,ix),'--','LineWidth',3)
semilogx(d,jXtotsd(:,ix),':','LineWidth',2)
xlabel('mixing rate')
ylabel('rates (/day)')
axis tight

subplot(235)
hold on
% set(gca,'LineWidth',3)
% semilogx(L,jR,'LineWidth',3)
% semilogx(d,jTotd(:,11),'--','LineWidth',3)
semilogx(d,jRd(:,11),'LineWidth',3)
semilogx(d,jTotd(:,11),'--','LineWidth',3)
semilogx(d,jXtotsd(:,11),':','LineWidth',2)
% semilogx(L,jXtots,':','LineWidth',2)
% semilogx(L(i),bg*simL(i).rates.jTot(1),'bd')
% semilogx(L(i),bJx_tot(1),'gd')
xlabel('mixing rate')
ylabel('rates (/day)')
legend('jR','jTot','\betaJx')
axis tight
title('respiration for smallest')

subplot(236)
% for i=1:length(L)
hold on
% set(gca,'LineWidth',3)
semilogx(d,Bdd,'LineWidth',3)
xlabel('mixing rate')
ylabel('rates (/day)')
title('biomass')
axis tight