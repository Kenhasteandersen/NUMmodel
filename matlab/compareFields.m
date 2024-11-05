% 
% Compare fields from 2 differents simulations (used in compareSimulations)
% 
% In:
%   field1 & field2 - fields to be compared (must be the same size)
%
% Out :
%   Fields - variation between field1 and field2 in %


function Field = compareFields(field1,field2)
arguments
    field1 double
    field2 double
end

Field = (field2-field1);
notZero = find(field2);
Field(notZero) = Field(notZero)./field2(notZero)*100;
Zero = find(field2==0 & field1~=0);
Field(Zero) = 100; %if field1 is not zero whereas field2 is zero the variation is set to 100%