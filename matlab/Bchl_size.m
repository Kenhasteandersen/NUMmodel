        for i= 1:(sim.p.n-sim.p.nNutrients)
         Bchl_size(:,:,:,i,:) = squeeze( sim.B(:,1:61,1,i,:).* sim.jLreal(:,:,:,i,:) )...
             ./squeeze(sim.L(:,1:61,1,:)); % in units of mu g Chl per l 
        end