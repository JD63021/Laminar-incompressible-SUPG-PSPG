function rowP = globalPressureRow(pNode, Nxy)
% GLOBALPRESSUREROW  Maps pressure node ID pNode (in [1..Npr])
%  to the global row index if the velocity DOFs go first in blocks x,y.
%
% Suppose x-block = 1..Nxy, y-block= Nxy+1..2*Nxy,
% so pressure block = (2*Nxy+1) .. (2*Nxy+Npr).
%
    rowP = 2*Nxy + pNode;
end
