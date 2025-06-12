function rowIndex = globalRow(vNode, comp, Nxy)
% GLOBALROW  Maps a velocity node ID vNode (in [1..Nxy]) and a component
%  ('x' or 'y') to the corresponding row index in the global vector.
%
% EXAMPLE:
%   rowX = globalRow(17, 'x', 81);   % => 17
%   rowY = globalRow(17, 'y', 81);   % => 17 + 81 = 98
%
    switch comp
        case 'x'
            rowIndex = vNode;
        case 'y'
            rowIndex = vNode + Nxy;
        otherwise
            error('Invalid velocity component (use ''x'' or ''y'')');
    end
end
