function ix = findIxUnicellular(p)

ix = [];
for iGroup = 1:p.nGroups
    if p.typeGroups(iGroup)<10
        ix = [ix p.ixStart(iGroup):p.ixEnd(iGroup)];
    end
end
      
ix = ix - p.idxB+1;