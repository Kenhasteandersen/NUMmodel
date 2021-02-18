%
% Returns the biomass of pico, nano and micro plankton in mugC/liter.
% Only works with generalists.
% Example: call calcPicoNanoMicro(sim.B(end,:), sim.p.pGeneralists);
%
function Bpnm = calcPicoNanoMicro(B,p)

ESD = 10000 * 1.5 * (p.m*1e-6).^(1/3);

ixPico = find(ESD<2);
ixNano = find(ESD>=2 & ESD<20);
ixMicro = find(ESD>=20);

lnESD = log(ESD);
delta = diff(lnESD);
delta = delta(1);
lnESDmid = lnESD(1:end-1)+0.5*delta;

Bpico = sum(B(ixPico));
Bnano = sum(B(ixNano));
Bmicro = sum(B(ixMicro));

if (lnESDmid(ixPico(end))<log(2))
    tmp = B(ixPico(end)+1)*(lnESD(ixPico(end)+1)-log(2))/delta;
    Bpico = Bpico + tmp;
    Bnano = Bnano - tmp;
else
    tmp = B(ixPico(end))*(lnESD(ixPico(end))-log(2))/delta;
    Bpico = Bpico - tmp;
    Bnano = Bnano + tmp;
end

if (lnESDmid(ixNano(end))<log(20))
    tmp = B(ixNano(end)+1)*(lnESD(ixNano(end)+1)-log(20))/delta;
    Bnano = Bnano + tmp;
    Bmicro = Bmicro - tmp;
else
    tmp = B(ixNano(end))*(lnESD(ixNano(end))-log(20))/delta;
    Bnano = Bnano - tmp;
    Bmicro = Bmicro + tmp;
end

Bpnm = [Bpico, Bnano, Bmicro];

