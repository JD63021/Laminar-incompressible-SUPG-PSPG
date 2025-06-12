function uv_new = update_bc(uv, boundaryInfo, nodeInfo, Nxy, t, corner, bcFlags,inletProfile)
% update_bc applies Dirichlet boundary conditions for the driven cavity.
%
% Inputs:
%   uv           : solution vector [2*Nxy + Npr]
%   boundaryInfo : structure with boundary info (fields like 'flag_2', etc.)
%   nodeInfo     : structure with node coordinates.
%   Nxy          : number of velocity DOFs per component.
%   t            : current time.
%   corner       : index for pressure correction.
%   bcFlags      : structure with fields:
%                    .top   - flag number for top boundary.
%                    .sides - array of flag numbers for side/bottom boundaries.
%
% Output:
%   uv_new       : solution vector with boundary conditions applied.

uv_new = uv;
H = max(nodeInfo.velocity.y) - min(nodeInfo.velocity.y);

% --- Apply top boundary condition (Dirichlet inlet)
topFlag = bcFlags.inlet;
if isfield(boundaryInfo, ['flag_' num2str(topFlag)])
    inletNodes = boundaryInfo.(['flag_' num2str(topFlag)])(:).';
else
    inletNodes = [];
end
for nodeID = inletNodes
    rowX = globalRow(nodeID, 'x', Nxy);
    rowY = globalRow(nodeID, 'y', Nxy);
    Uin = inletProfile(t, nodeInfo.velocity.y(nodeID), H);
    uv_new(rowX) = Uin;
    uv_new(rowY) = 0;
end

% --- Apply side and bottom boundary conditions (zero velocity)
sides = bcFlags.wall;
sideNodes = [];
for i = 1:length(sides)
    flagStr = ['flag_' num2str(sides(i))];
    if isfield(boundaryInfo, flagStr)
        sideNodes = [sideNodes; boundaryInfo.(flagStr)(:)];
    end
end
sideNodes = unique(sideNodes);
for nodeID = sideNodes'
    rowX = globalRow(nodeID, 'x', Nxy);
    rowY = globalRow(nodeID, 'y', Nxy);
    uv_new(rowX) = 0;
    uv_new(rowY) = 0;
end

% --- Fix corner pressure DOF
cornerGlobal = 2*Nxy + corner;
uv_new(cornerGlobal) = 0;
end

function row = globalRow(vNode, comp, Nxy)
% globalRow computes the global index for a given velocity node and component.
switch comp
    case 'x'
        row = vNode;
    case 'y'
        row = vNode + Nxy;
    otherwise
        error('Invalid component. Use "x" or "y".');
end
end

