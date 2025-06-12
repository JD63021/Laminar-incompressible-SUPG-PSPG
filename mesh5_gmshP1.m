function [nodeInfo, elemInfo, boundaryInfo] = mesh5_gmshP1( ...
    gmshVelFile, gmshPresFile, ...
    excludedVelFlags, pBoundFlags)
% mesh5_gmsh_p1p1 reads two Gmsh .m files for a P1–P1 triangular approach:
%   (1) Velocity mesh => P1 => TRIANGLES
%       boundary => LINES
%   (2) Pressure mesh => P1 => TRIANGLES
%       boundary => LINES
%
% And builds:
%   nodeInfo.velocity (P1),
%   nodeInfo.pressure (P1),
%   elemInfo.velElements (3 nodes each),
%   elemInfo.presElements (3 nodes each),
%   boundaryInfo.velLine2Elements.flag_<X> => Nx2 lines from velocity mesh,
%   boundaryInfo.presLine2Elements.flag_<X> => Nx2 lines from pressure mesh.
%
% 'excludedVelFlags' (if given) are velocity boundary flags we ignore.
% 'pBoundFlags' can be used similarly if you want to handle or store
% additional data for the pressure boundary.
%
% Usage:
%   [nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh_p1p1('velMesh.m','presMesh.m')
%
% This function assumes that each .m file from Gmsh (via "Save as MATLAB") 
% loads a structure 'msh' with fields:
%   msh.POS         -> node coordinates
%   msh.TRIANGLES   -> Nx4: the first 3 columns are node IDs of the triangle, 
%                             the 4th is a region or partition index (optional)
%   msh.LINES       -> Nx3: the first 2 columns are node IDs of the line segment,
%                             the 3rd is the boundary flag
%
% ----------------------------------------------------------------------

if ~exist('excludedVelFlags','var') || isempty(excludedVelFlags)
    excludedVelFlags = [];
end
if ~exist('pBoundFlags','var') || isempty(pBoundFlags)
    pBoundFlags = [];
end

%% (A) Read Velocity Mesh => P1 => TRIANGLES
if ischar(gmshVelFile)
    run(gmshVelFile);  % loads struct 'msh' for velocity
    msh_vel = msh;
else
    error('Velocity mesh file not provided as a string.');
end

if ~isfield(msh_vel,'TRIANGLES')
    error('Expected TRIANGLES in the velocity mesh (P1). Check your Gmsh export.');
end

triElementsVel = msh_vel.TRIANGLES(:,1:3);
numElemsVel    = size(triElementsVel,1);

allVelCoords = msh_vel.POS;  % Nx(2 or 3)
numVelNodes  = size(allVelCoords,1);

vel_nodes.id = (1:numVelNodes).';
vel_nodes.x  = allVelCoords(:,1);
vel_nodes.y  = allVelCoords(:,2);

%% (B) Read Pressure Mesh => P1 => TRIANGLES
if ischar(gmshPresFile)
    run(gmshPresFile); % loads struct 'msh' for pressure
    msh_pres = msh;
else
    error('Pressure mesh file not provided as a string.');
end

if ~isfield(msh_pres,'TRIANGLES')
    error('Expected TRIANGLES in the pressure mesh (P1). Check your Gmsh export.');
end

triElementsPres = msh_pres.TRIANGLES(:,1:3);
numElemsPres    = size(triElementsPres,1);

allPresCoords = msh_pres.POS;
numPresNodes  = size(allPresCoords,1);

pres_nodes.id = (1:numPresNodes).';
pres_nodes.x  = allPresCoords(:,1);
pres_nodes.y  = allPresCoords(:,2);

%% (C) Store element connectivity
elemInfo.velElements  = triElementsVel;   % 3 nodes each
elemInfo.presElements = triElementsPres;  % 3 nodes each

%% (D) Fill nodeInfo
nodeInfo.velocity = vel_nodes;    % P1 for velocity
nodeInfo.pressure = pres_nodes;   % P1 for pressure

%% (E) Build boundaryInfo
boundaryInfo = struct();

% 1) Velocity boundaries => LINES => Nx3 [nodeA nodeB flag]
if ~isfield(msh_vel,'LINES') || isempty(msh_vel.LINES)
    warning('No LINES found in velocity mesh => boundaryInfo.velLine2Elements = empty.');
    boundaryInfo.velLine2Elements = struct();
else
    lines2_vel = msh_vel.LINES; % Nx3 => [nA nB flag]
    uniqueVFlags = unique(lines2_vel(:,3));
    velLine2Struct = struct();
    for iF = 1:numel(uniqueVFlags)
        gFlag = uniqueVFlags(iF);
        % Optionally ignore certain flags
        if ismember(gFlag, excludedVelFlags)
            continue;
        end
        mask = (lines2_vel(:,3) == gFlag);
        theseLines = lines2_vel(mask, 1:2);  % columns 1..2 => node IDs
        fieldName  = sprintf('flag_%d', gFlag);
        velLine2Struct.(fieldName) = theseLines;
    end
    boundaryInfo.velLine2Elements = velLine2Struct;
end

% 2) Pressure boundaries => LINES => Nx3 [nodeA nodeB flag]
if ~isfield(msh_pres,'LINES') || isempty(msh_pres.LINES)
    warning('No LINES found in pressure mesh => boundaryInfo.presLine2Elements = empty.');
    boundaryInfo.presLine2Elements = struct();
else
    lines2_pres = msh_pres.LINES; % Nx3 => [nA nB flag]
    uniquePFlags = unique(lines2_pres(:,3));
    presLine2Struct = struct();
    for iF = 1:numel(uniquePFlags)
        pFlag = uniquePFlags(iF);
        maskP = (lines2_pres(:,3) == pFlag);
        theseLinesP = lines2_pres(maskP,1:2);
        fNameP = sprintf('flag_%d', pFlag);
        presLine2Struct.(fNameP) = theseLinesP;
    end
    boundaryInfo.presLine2Elements = presLine2Struct;
end

%% (F) Optionally collect all velocity boundary nodes
boundaryInfo.allVelNodes = [];
if isfield(boundaryInfo,'velLine2Elements')
    fnames = fieldnames(boundaryInfo.velLine2Elements);
    for iF = 1:numel(fnames)
        thisField = fnames{iF}; % e.g. 'flag_5'
        lineEls   = boundaryInfo.velLine2Elements.(thisField);
        if isempty(lineEls), continue; end
        boundaryNodes = unique(lineEls(:));
        boundaryInfo.(thisField) = boundaryNodes;  % store them
        boundaryInfo.allVelNodes = [boundaryInfo.allVelNodes; boundaryNodes];
    end
    boundaryInfo.allVelNodes = unique(boundaryInfo.allVelNodes);
end

% If you wish to store pressure boundary flags similarly, you can do so:
% boundaryInfo.allPresNodes = [];
% if isfield(boundaryInfo,'presLine2Elements')
%     % etc. ...
% end

end


% function [nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh( ...
%     gmshVelFile, gmshPresFile, ...
%     excludedVelFlags, pBoundFlags)
% % mesh6_gmsh reads two Gmsh files for a P2–P1 triangular approach:
% %   (1) Velocity mesh => P2 => TRIANGLES6, boundary => LINES3
% %   (2) Pressure mesh => P1 => TRIANGLES,  boundary => LINES
% %
% % And it builds:
% %   nodeInfo.velocity (P2) + nodeInfo.pressure (P1),
% %   elemInfo.velElements (6 nodes) + elemInfo.presElements (3 nodes),
% %   boundaryInfo.velLine3Elements.flag_<X> => Nx3 lines from velocity mesh,
% %   boundaryInfo.presLine2Elements.flag_<X> => Nx2 lines from pressure mesh.
% %
% % If 'excludedVelFlags' is specified, those velocity boundary segments are ignored.
% % If pBoundFlags/vals are specified, the code sets pressure boundary info as well.
% 
% if ~exist('excludedVelFlags','var') || isempty(excludedVelFlags)
%     excludedVelFlags = [];
% end
% if ~exist('pBoundFlags','var') || isempty(pBoundFlags)
%     pBoundFlags = [];
% end
% if ~exist('pBoundVals','var') || isempty(pBoundVals)
%     pBoundVals = [];
% end
% 
% %% (A) Read Velocity Mesh => P2 => TRIANGLES6
% if ischar(gmshVelFile)
%     run(gmshVelFile);  % loads struct 'msh' into workspace
%     msh_vel = msh;
% else
%     error('Velocity mesh file not provided as a string.');
% end
% 
% if ~isfield(msh_vel,'TRIANGLES6')
%     error('Expected TRIANGLES6 in the velocity mesh (P2).');
% end
% 
% triElementsVel = msh_vel.TRIANGLES6(:,1:6);
% numElemsVel    = size(triElementsVel,1);
% 
% allVelCoords = msh_vel.POS;    % Nx3 or Nx2
% numVelNodes  = size(allVelCoords,1);
% 
% vel_nodes.id = (1:numVelNodes).';
% vel_nodes.x  = allVelCoords(:,1);
% vel_nodes.y  = allVelCoords(:,2);
% 
% %% (B) Read Pressure Mesh => P1 => TRIANGLES
% if ischar(gmshPresFile)
%     run(gmshPresFile);
%     msh_pres = msh;
% else
%     error('Pressure mesh file not provided as a string.');
% end
% 
% if ~isfield(msh_pres,'TRIANGLES')
%     error('Expected TRIANGLES in the pressure mesh (P1).');
% end
% 
% triElementsPres = msh_pres.TRIANGLES(:,1:3);
% numElemsPres    = size(triElementsPres,1);
% 
% allPresCoords = msh_pres.POS;
% numPresNodes  = size(allPresCoords,1);
% 
% pres_nodes.id = (1:numPresNodes).';
% pres_nodes.x  = allPresCoords(:,1);
% pres_nodes.y  = allPresCoords(:,2);
% 
% %% (C) Store element connectivity
% elemInfo.velElements  = triElementsVel;  % 6 nodes each
% elemInfo.presElements = triElementsPres; % 3 nodes each
% 
% %% (D) Fill nodeInfo
% nodeInfo.velocity = vel_nodes;
% nodeInfo.pressure = pres_nodes;
% 
% %% (E) Build boundaryInfo
% boundaryInfo = struct();
% 
% % 1) Velocity boundaries => LINES3 => Nx4 [nodeA nodeMid nodeB flag]
% if ~isfield(msh_vel,'LINES3') || isempty(msh_vel.LINES3)
%     warning('No LINES3 found in velocity mesh => boundaryInfo.velLine3Elements = empty.');
%     boundaryInfo.velLine3Elements = struct();
% else
%     lines3_vel = msh_vel.LINES3; % Nx4 => [nA nM nB flag]
%     uniqueVFlags = unique(lines3_vel(:,4));
%     velLine3Struct = struct();
%     for iF = 1:numel(uniqueVFlags)
%         gFlag = uniqueVFlags(iF);
%         if ismember(gFlag, excludedVelFlags)
%             continue;
%         end
%         mask = (lines3_vel(:,4) == gFlag);
%         theseLines = lines3_vel(mask, 1:3);  % columns 1..3 => node IDs
%         fieldName  = sprintf('flag_%d', gFlag);
%         velLine3Struct.(fieldName) = theseLines;
%     end
%     boundaryInfo.velLine3Elements = velLine3Struct;
% end
% 
% % 2) Pressure boundaries => LINES => Nx3 [nodeA nodeB flag]
% if ~isfield(msh_pres,'LINES') || isempty(msh_pres.LINES)
%     warning('No LINES found in pressure mesh => boundaryInfo.presLine2Elements = empty.');
%     boundaryInfo.presLine2Elements = struct();
% else
%     lines2_pres = msh_pres.LINES; % Nx3 => [nA nB flag]
%     uniquePFlags = unique(lines2_pres(:,3));
%     presLine2Struct = struct();
%     for iF = 1:numel(uniquePFlags)
%         pFlag = uniquePFlags(iF);
%         maskP = (lines2_pres(:,3) == pFlag);
%         theseLinesP = lines2_pres(maskP,1:2);
%         fNameP = sprintf('flag_%d', pFlag);
%         presLine2Struct.(fNameP) = theseLinesP;
%     end
%     boundaryInfo.presLine2Elements = presLine2Struct;
% end
% 
% %% (F) Optionally store velocity boundary node sets
% boundaryInfo.allVelNodes = [];
% if isfield(boundaryInfo,'velLine3Elements')
%     fnames = fieldnames(boundaryInfo.velLine3Elements);
%     for iF = 1:numel(fnames)
%         thisField = fnames{iF}; % e.g. 'flag_5'
%         lineEls   = boundaryInfo.velLine3Elements.(thisField);
%         if isempty(lineEls), continue; end
%         boundaryNodes = unique(lineEls(:));
%         boundaryInfo.(thisField) = boundaryNodes;  % store them
%         boundaryInfo.allVelNodes = [boundaryInfo.allVelNodes; boundaryNodes];
%     end
%     boundaryInfo.allVelNodes = unique(boundaryInfo.allVelNodes);
% end
