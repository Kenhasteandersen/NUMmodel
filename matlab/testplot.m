figure('color','w')
tiledlayout(3,2)
t1=nexttile(1);
pcolor(rand(100,100))
t3=nexttile(3);
pcolor(rand(100,100))
cbar=colorbar;
cbar.Position(4)=(t1.Position(2)-t3.Position(2)+t1.Position(4));

nexttile(2)
plot(1:10)
nexttile(4)
plot(1:20)
nexttile(6)
plot(1:30)

