function testCase = testcase()
% fpc - Test case configuration for the NS problem.
%
% This file specifies:
%   (i) the gmsh files to run (now two files: one for velocity and one for pressure),
%   (ii) any velocity boundary flags to exclude (if desired),
%   (iii) pressure boundary flags and values (if desired), and 
%   (iv) the boundary flags to use for enforcing Dirichlet conditions.
%
% 

    testCase.gmshVelFile  = 'quadcp2.m';  % gmsh file for velocity grid (P2 or P1)
    testCase.gmshPresFile = 'quadcp1.m';  % gmsh file for pressure grid (P1)
    
    testCase.excludedVelFlags = [];             % velocity boundary flags to exclude
    testCase.pBoundFlags      = [];                % (no pressure BCs in this case)
    testCase.pBoundVals       = [];
    testCase.boundaryFlags.inlet = 3;              % top boundary flag (for inlet)
    testCase.boundaryFlags.wall  = [1,2,4]; % side and bottom boundaries

    % Specify the inlet velocity profile (function handle).
    % This function should accept: (t, y, H) and return the x-velocity.
    % testCase.inletProfile   = @(t, y, H) 4*1.5*sin(pi*t/8)*(y*(H-y))/(H^2);
    testCase.inletProfileSS = @(t, y, H) 1;%4*0.3*(y*(H-y))/(H^2);
    testCase.corner = 2;%113;
end


