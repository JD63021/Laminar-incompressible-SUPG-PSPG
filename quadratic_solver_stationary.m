% quadratic_solver.m
% =========================================================================
% This driver sets up the mesh, initial conditions, and calls the stationary
% solver (main1_SS) for the Navier–Stokes problem.
% =========================================================================

tic
clc
warning off MATLAB:NonScalarInput
warning off MATLAB:nearlySingularMatrix

%% Load test case configuration6
% testCase = fpc();

%% User parameters
D    = 1;      % Diffusivity (or similar)
Re   = 1;
mu   = 1/1;
rho  = 1;
o    = 1;
t1   = 0;     % initial time (for stationary solver, this can be arbitrary)

%% Generate Mesh & Setup

testCase = testcase();

%% Generate Mesh & Setup
[nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh(...
    testCase.gmshVelFile, testCase.gmshPresFile, ...
    testCase.excludedVelFlags, testCase.pBoundFlags);

Nxy = length(nodeInfo.velocity.x);  % Number of velocity DOFs
Npr = max(elemInfo.presElements(:));


% Set initial conditions.
if exist('U','var')
    U1 = U(1:Nxy); 
    U2 = U(Nxy+1:2*Nxy); 
    U3 = U(2*Nxy+1:end);
else
    U1 = zeros(Nxy,1);
    U2 = zeros(Nxy,1);
    U3 = zeros(Npr,1);
end

U_ms = zeros(Nxy,1); V_ms = zeros(Nxy,1); fu_ms = zeros(Nxy,1); fv_ms = zeros(Nxy,1);
P_ms = zeros(Npr,1);



%% Stationary solver parameters (for a steady-state problem)
dt    = 1;   % Not used here
nt    = 1;   % (Single "time-step" for steady state)
gamma = 0;

%% Call the stationary solver.


[U, x, y,x1,y1, z1, z2,z3, z, nodeInfo, elemInfo, boundaryInfo,K] = ...
    main_stationary_picard( ...
      U1, U2, U3, ...           % initial guesses for x-vel, y-vel, pressure
      D, Re, o, ...             % PDE parameters
      gamma, mu, rho, ...       % PDE parameters
      nodeInfo, elemInfo, boundaryInfo, ...
      testCase.corner, ...          % which pressure node to pin
      200, 1e-6,testCase.boundaryFlags,testCase.inletProfileSS );



% Update velocity and pressure parts from the solution.
U1 = U(1:Nxy);
U2 = U(Nxy+1:2*Nxy);
U3 = U(2*Nxy+1:end);


%% --- Write VTP Files for Diagnostic Snapshots ---

    filename = sprintf('x velocity.vtp');
    writeVTP(x, y, U1, filename);
    fprintf('Wrote snapshot %d to file %s (mu = %f)\n', 1, filename, mu);

% %% --- Generate a PVD File Linking the Snapshots ---
% pvdFilename = 'solutions.pvd';
% baseFilename = 'solution_%04d.vtp';
% writePVD(recordedTimes, baseFilename, pvdFilename);






%% (Optional) Visualize the final velocity magnitude.
figure;
scatter3(x, y, sqrt(U1.^2 + U2.^2), 20, sqrt(U1.^2 + U2.^2), 'filled');
title('Final Velocity Magnitude');
xlabel('x'); ylabel('y'); zlabel('Velocity Magnitude');
colorbar;

toc

figure;
scatter3(x1, y1, U3, 20, U3, 'filled');
title('Pressure Magnitude');
xlabel('x'); ylabel('y'); zlabel('Pressure Magnitude');
colorbar;

figure;
scatter3(x, y, U1 , 20, U1, 'filled');
title('x Velocity Magnitude');
xlabel('x'); ylabel('y'); zlabel('x Velocity Magnitude');
colorbar;

figure;
scatter3(x, y, U2 , 20, U2, 'filled');
title('y Velocity Magnitude');
xlabel('x'); ylabel('y'); zlabel('y Velocity Magnitude');
colorbar;


