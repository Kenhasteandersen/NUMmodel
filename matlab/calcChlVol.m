%
% Calculate the functions from any simulation.
%
% In:
%   sim - simulation structure
%
% Out:
%   sim - same as input but with fields of global function added:
%
%         For global simulations additional fields are:
%         sim.Ntotal    - total dissolved N as a function of time (mugN)
%         sim.DOCtotal  - total dissolved DOC (mugC)

%
function sim = calcChlVol(sim)

% switch sim.p.nameModel
%         
%     case 'global'
        %
        % Primary production in gC per m2/year::
        %
        if ~isfield(sim,'ProdGross')
            % Get grid volumes:
            load(sim.p.pathGrid,'dv','dz','dx','dy');
            ix = ~isnan(sim.N(:,:,1,1)); % Find all relevant grid cells
            
            aRatio = zeros(length(sim.x), length(sim.y),1,sim.p.n-sim.p.idxB+1,length(sim.t)); %formulation that can work with diatoms
            jLreal = zeros(length(sim.t),length(sim.x), length(sim.y), length(sim.z),sim.p.n-sim.p.idxB+1);
            jL = zeros(length(sim.t),length(sim.x), length(sim.y), length(sim.z),sim.p.n-sim.p.idxB+1);
            sim.BB= zeros(length(sim.x), length(sim.y),1,sim.p.n-sim.p.idxB+1,length(sim.t));

            ChlArea = 0*dx;
            ChlVolume = zeros(length(sim.t),length(sim.x), length(sim.y),length(sim.z));
            
            nTime = length(sim.t);
            nX = length(sim.x);
            nY = length(sim.y);
            nZ = length(sim.z);
            for iTime = 1:nTime
                for i = 1:nX
                    for j = 1:nY
                       % for w = 1:sim.p.n-2
                        k=1;% for k = 1:nZ
                            if ~isnan(sim.N(i,j,k,iTime))
                                if sim.p.nNutrients==3
                                u = [squeeze(sim.N(i,j,k,iTime)), ...
                                    squeeze(sim.DOC(i,j,k,iTime)), ...
                                    squeeze(sim.Si(i,j,k,iTime)), ...
                                    squeeze(sim.B(i,j,k,:,iTime))'];
                                else
                                u = [squeeze(sim.N(i,j,k,iTime)), ...
                                    squeeze(sim.DOC(i,j,k,iTime)), ...
                                    squeeze(sim.B(i,j,k,:,iTime))'];
                                end
                                % Chl:
                                rates = getRates(sim.p, u, sim.L(i,j,k,iTime), sim.T(i,j,k,iTime));
                                tmp =  calcChl( squeeze(sim.B(i,j,k,:,iTime)), rates, sim.L(i,j,k,iTime));
%                                 sim.BB(i,j,k,:,iTime) =  calcChl( sim.B(i,j,k,:,iTime), rates, sim.L(i,j,k,iTime)); 
%                                 aRatio(i,j,k,:,iTime) = rates.jLreal./rates.jL;
%                                 jL(iTime,i,j,k,:) = rates.jL;
%                                 jLreal(iTime,i,j,k,:) = rates.jLreal;
%                                 sim.jLreal(i,j,k,:,iTime) = rates.jLreal;
%                                 sim.jL(i,j,k,:,iTime) = rates.jL;
                                if ~isnan(tmp)
%                                     ChlArea(i,j,:) = ChlArea(i,j,:) + tmp .* dz(i,j,k);
                                    ChlVolume(iTime,i,j,k) = ChlVolume(iTime,i,j,k) + tmp;
%                                       sim.BB(i,j,k,:,iTime) = sim.BB(i,j,k,:,iTime) + tmp;
                                end
                            end
%                         end
%                         sim.aRatio(i,j,k,:,iTime)=aRatio;
                        %end
                    end
                end
            end
            sim.ChlArea = ChlArea/length(sim.t);
            sim.ChlVolume = ChlVolume/length(sim.t);
            sim.jLreal = jLreal;
            sim.jL = jL;
%             sim.aRatio=aRatio;
            %
            % Global totals
            %
            calcTotal = @(u) sum(u(ix(:)).*dv(ix(:))); % mug/day
            
            for i = 1:length(sim.t)
                sim.Ntotal(i) = calcTotal(sim.N(:,:,:,i));
                sim.DOCtotal(i) = calcTotal(sim.DOC(:,:,:,i)); % mugC
                sim.Btotal(i) = 0;
                for j = 1:(sim.p.n-sim.p.idxB+1)
                    sim.Btotal(i) = sim.Btotal(i) + calcTotal(sim.B(:,:,:,j,i)); % mugC
                end
            end
            
            sim.ChlTotal = ChlArea.*dx.*dy;
            sim.ChlTotal = sum(sim.ChlTotal(:))/ 1000; %gChl
        end

end
