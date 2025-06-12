function [nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh( ...
    gmshVelFile, gmshPresFile, ...
    excludedVelFlags, pBoundFlags)
% mesh6_gmsh reads two Gmsh files for a P2–P1 triangular approach:
%   (1) Velocity mesh => P2 => TRIANGLES6, boundary => LINES3
%   (2) Pressure mesh => P1 => TRIANGLES,  boundary => LINES
%
% And it builds:
%   nodeInfo.velocity (P2) + nodeInfo.pressure (P1),
%   elemInfo.velElements (6 nodes) + elemInfo.presElements (3 nodes),
%   boundaryInfo.velLine3Elements.flag_<X> => Nx3 lines from velocity mesh,
%   boundaryInfo.presLine2Elements.flag_<X> => Nx2 lines from pressure mesh.
%
% If 'excludedVelFlags' is specified, those velocity boundary segments are ignored.
% If pBoundFlags/vals are specified, the code sets pressure boundary info as well.

if ~exist('excludedVelFlags','var') || isempty(excludedVelFlags)
    excludedVelFlags = [];
end
if ~exist('pBoundFlags','var') || isempty(pBoundFlags)
    pBoundFlags = [];
end
if ~exist('pBoundVals','var') || isempty(pBoundVals)
    pBoundVals = [];
end

%% (A) Read Velocity Mesh => P2 => TRIANGLES6
if ischar(gmshVelFile)
    run(gmshVelFile);  % loads struct 'msh' into workspace
    msh_vel = msh;
else
    error('Velocity mesh file not provided as a string.');
end

if ~isfield(msh_vel,'TRIANGLES6')
    error('Expected TRIANGLES6 in the velocity mesh (P2).');
end

triElementsVel = msh_vel.TRIANGLES6(:,1:6);
numElemsVel    = size(triElementsVel,1);

allVelCoords = msh_vel.POS;    % Nx3 or Nx2
numVelNodes  = size(allVelCoords,1);

vel_nodes.id = (1:numVelNodes).';
vel_nodes.x  = allVelCoords(:,1);
vel_nodes.y  = allVelCoords(:,2);

%% (B) Read Pressure Mesh => P1 => TRIANGLES
if ischar(gmshPresFile)
    run(gmshPresFile);
    msh_pres = msh;
else
    error('Pressure mesh file not provided as a string.');
end

if ~isfield(msh_pres,'TRIANGLES')
    error('Expected TRIANGLES in the pressure mesh (P1).');
end

triElementsPres = msh_pres.TRIANGLES(:,1:3);
numElemsPres    = size(triElementsPres,1);

allPresCoords = msh_pres.POS;
numPresNodes  = size(allPresCoords,1);

pres_nodes.id = (1:numPresNodes).';
pres_nodes.x  = allPresCoords(:,1);
pres_nodes.y  = allPresCoords(:,2);

%% (C) Store element connectivity
elemInfo.velElements  = triElementsVel;  % 6 nodes each
elemInfo.presElements = triElementsPres; % 3 nodes each

%% (D) Fill nodeInfo
nodeInfo.velocity = vel_nodes;
nodeInfo.pressure = pres_nodes;

%% (E) Build boundaryInfo
boundaryInfo = struct();

% 1) Velocity boundaries => LINES3 => Nx4 [nodeA nodeMid nodeB flag]
if ~isfield(msh_vel,'LINES3') || isempty(msh_vel.LINES3)
    warning('No LINES3 found in velocity mesh => boundaryInfo.velLine3Elements = empty.');
    boundaryInfo.velLine3Elements = struct();
else
    lines3_vel = msh_vel.LINES3; % Nx4 => [nA nM nB flag]
    uniqueVFlags = unique(lines3_vel(:,4));
    velLine3Struct = struct();
    for iF = 1:numel(uniqueVFlags)
        gFlag = uniqueVFlags(iF);
        if ismember(gFlag, excludedVelFlags)
            continue;
        end
        mask = (lines3_vel(:,4) == gFlag);
        theseLines = lines3_vel(mask, 1:3);  % columns 1..3 => node IDs
        fieldName  = sprintf('flag_%d', gFlag);
        velLine3Struct.(fieldName) = theseLines;
    end
    boundaryInfo.velLine3Elements = velLine3Struct;
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

%% (F) Optionally store velocity boundary node sets
boundaryInfo.allVelNodes = [];
if isfield(boundaryInfo,'velLine3Elements')
    fnames = fieldnames(boundaryInfo.velLine3Elements);
    for iF = 1:numel(fnames)
        thisField = fnames{iF}; % e.g. 'flag_5'
        lineEls   = boundaryInfo.velLine3Elements.(thisField);
        if isempty(lineEls), continue; end
        boundaryNodes = unique(lineEls(:));
        boundaryInfo.(thisField) = boundaryNodes;  % store them
        boundaryInfo.allVelNodes = [boundaryInfo.allVelNodes; boundaryNodes];
    end
    boundaryInfo.allVelNodes = unique(boundaryInfo.allVelNodes);
end

%% (G) Pressure boundary conditions (if requested)
% boundaryInfo.pressureConditions = struct('flag',{},'nodes',{},'value',{});
% if ~isempty(pBoundFlags)
%     for iPf = 1:numel(pBoundFlags)
%         thisFlag  = pBoundFlags(iPf);
%         thisValue = pBoundVals(iPf);
%         boundaryInfo.pressureConditions(iPf).flag  = thisFlag;
%         boundaryInfo.pressureConditions(iPf).nodes = [];  % or find from LINES
%         boundaryInfo.pressureConditions(iPf).value = thisValue;
%     end
% end

end


% function [nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh( ...
%     gmshVelFile, gmshPresFile, ...
%     excludedVelFlags, pBoundFlags, pBoundVals)
% % mesh5_gmsh reads two Gmsh files for a Q2–Q1 approach:
% %   (1) Velocity mesh: Q2 => QUADS9, boundary lines => LINES3
% %   (2) Pressure mesh: Q1 => QUADS4, boundary lines => LINES
% %
% % and it builds:
% %   nodeInfo.velocity (Q2) + nodeInfo.pressure (Q1),
% %   elemInfo.quadElements (9 nodes for velocity) + elemInfo.presElements (4 nodes),
% %   boundaryInfo.velLine3Elements.flag_<X> => Nx3 lines from velocity mesh,
% %   boundaryInfo.presLine2Elements.flag_<X> => Nx2 lines from pressure mesh.
% %
% % This way, you can directly integrate each boundary line in your force integral,
% % avoiding partial edges or double counting.
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
% %% (A) Read Velocity Mesh => Q2
% if ischar(gmshVelFile)
%     run(gmshVelFile);  % populates variable 'msh'
%     msh_vel = msh;
% else
%     error('Velocity mesh file not provided as a string.');
% end
% 
% % Q2 means volume elements are QUADS9
% % boundary lines are typically LINES3 => (nodeA, midNode, nodeB, flag)
% 
% allVelCoords = msh_vel.POS;
% numVelNodes  = size(allVelCoords,1);
% vel_nodes.id = (1:numVelNodes).';
% vel_nodes.x  = allVelCoords(:,1);
% vel_nodes.y  = allVelCoords(:,2);
% 
% % For Q2 volume elements:
% if ~isfield(msh_vel,'QUADS9')
%     error('Expected QUADS9 in the velocity mesh (Q2).');
% end
% quadElementsVel = msh_vel.QUADS9(:,1:9);
% numElemsVel     = size(quadElementsVel,1);
% 
% %% (B) Read Pressure Mesh => Q1
% if ischar(gmshPresFile)
%     run(gmshPresFile); % populates variable 'msh'
%     msh_pres = msh;
% else
%     error('Pressure mesh file not provided as a string.');
% end
% 
% % Q1 means volume elements are QUADS4
% % boundary lines typically LINES => (nodeA, nodeB, flag)
% 
% allPresCoords = msh_pres.POS;
% numPresNodes  = size(allPresCoords,1);
% pres_nodes.id = (1:numPresNodes).';
% pres_nodes.x  = allPresCoords(:,1);
% pres_nodes.y  = allPresCoords(:,2);
% 
% if ~isfield(msh_pres,'QUADS')
%     error('Expected QUADS in the pressure mesh (Q1).');
% end
% quadElementsPres = msh_pres.QUADS(:,1:4);
% numElemsPres     = size(quadElementsPres,1);
% 
% %% (C) Construct final connectivity
% % For velocity => Q2
% elemInfo.quadElements = quadElementsVel;  % 9 velocity nodes
% % For pressure => Q1
% elemInfo.presElements = quadElementsPres; % 4 pressure nodes
% 
% %% (D) Fill nodeInfo
% nodeInfo.velocity = vel_nodes;
% nodeInfo.pressure = pres_nodes;
% 
% %% (E) Build boundaryInfo
% 
% boundaryInfo = struct();
% 
% % 1) Velocity boundaries => look for LINES3, ignoring excludedVelFlags
% if ~isfield(msh_vel,'LINES3') || isempty(msh_vel.LINES3)
%     warning('No LINES3 found in the velocity mesh. boundaryInfo.velLine3Elements = empty.');
%     boundaryInfo.velLine3Elements = struct();
% else
%     lines3_vel = msh_vel.LINES3; % Nx4 => [nodeA nodeMid nodeB flag]
%     uniqueVFlags = unique(lines3_vel(:,4));
% 
%     velLine3Struct = struct();
%     for iF = 1:numel(uniqueVFlags)
%         gFlag = uniqueVFlags(iF);
%         if ismember(gFlag, excludedVelFlags)
%             continue;
%         end
% 
%         mask = (lines3_vel(:,4) == gFlag);
%         theseLines = lines3_vel(mask,1:3);  % columns 1..3 are node IDs
%         fieldName  = sprintf('flag_%d', gFlag);
%         velLine3Struct.(fieldName) = theseLines;
%     end
%     boundaryInfo.velLine3Elements = velLine3Struct;
% end
% 
% % 2) Pressure boundaries => look for LINES
% if ~isfield(msh_pres,'LINES') || isempty(msh_pres.LINES)
%     warning('No LINES found in the pressure mesh. boundaryInfo.presLine2Elements = empty.');
%     boundaryInfo.presLine2Elements = struct();
% else
%     lines2_pres = msh_pres.LINES; % Nx3 => [nodeA nodeB flag]
%     uniquePFlags = unique(lines2_pres(:,3));
% 
%     presLine2Struct = struct();
%     for iF = 1:numel(uniquePFlags)
%         pFlag = uniquePFlags(iF);
%         maskP = (lines2_pres(:,3) == pFlag);
%         theseLinesP = lines2_pres(maskP, 1:2);  % nodeA, nodeB
%         fNameP = sprintf('flag_%d', pFlag);
%         presLine2Struct.(fNameP) = theseLinesP;
%     end
%     boundaryInfo.presLine2Elements = presLine2Struct;
% end
% 
% %% (F) Possibly build older 'flag_<X>' lists for velocity boundary
% %   if you want a single boundary chain for reference.
% %   We gather node IDs from the lines3 => ignoring the 'excludedVelFlags'.
% boundaryInfo.allNodes = [];
% if isfield(boundaryInfo,'velLine3Elements')
%     fnames = fieldnames(boundaryInfo.velLine3Elements);
%     for iF = 1:numel(fnames)
%         thisField = fnames{iF}; % 'flag_X'
%         lineEls = boundaryInfo.velLine3Elements.(thisField);
%         if isempty(lineEls), continue; end
%         boundaryNodes = unique(lineEls(:));
%         boundaryInfo.(thisField) = boundaryNodes;  % store them as well
%         boundaryInfo.allNodes = [boundaryInfo.allNodes; boundaryNodes];
%     end
%     boundaryInfo.allNodes = unique(boundaryInfo.allNodes);
% end
% 
% %% (G) Pressure boundary conditions (if provided)
% if ~isempty(pBoundFlags)
%     % If you want to gather pressure BC from the Q1 lines => 
%     % you'd do so similarly here. We'll skip or adapt.
%     boundaryInfo.pressureConditions = struct('flag',{},'nodes',{},'value',{});
%     for iPf = 1:numel(pBoundFlags)
%         thisFlag  = pBoundFlags(iPf);
%         thisValue = pBoundVals(iPf);
%         % Not fully implementing, but you could do lines2 search:
%         % pNodes = ...
%         boundaryInfo.pressureConditions(iPf).flag  = thisFlag;
%         boundaryInfo.pressureConditions(iPf).nodes = [];
%         boundaryInfo.pressureConditions(iPf).value = thisValue;
%     end
% end
% 
% end
