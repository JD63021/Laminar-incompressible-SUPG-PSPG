function [U,x,y,x1,y1,z1,z2,z3,z,nodeInfo,elemInfo,boundaryInfo,K_new] = ...
    main_stationary_picard( ...
      U1,U2,U3, ...                 % initial guess
      D,Re,o, ...                   % PDE parameters
      gamma,mu,rho, ...             % PDE parameters
      nodeInfo,elemInfo,boundaryInfo,...
      cornerIndex, ...              % pinned pressure node
      maxPicardIters,tolRes,bcFlags,inletProfile)

% ================================================================
% 0) basic sizes / helper handles
% ================================================================
Nxy   = numel(nodeInfo.velocity.x);
Npr   = max(elemInfo.presElements(:));
x     = nodeInfo.velocity.x;   y = nodeInfo.velocity.y;
x1    = nodeInfo.pressure.x;   y1 = nodeInfo.pressure.y;

if nargin<11||isempty(cornerIndex);  cornerIndex=1;      end
if nargin<12||isempty(maxPicardIters); maxPicardIters=60; end
if nargin<13||isempty(tolRes);       tolRes=1e-4;         end

fprintf('\n=== Stationary Picard – ω–relaxed ===\n');
fprintf('Re = %.0f   mu = %.3e   maxIter = %d   tol = %.1e\n',...
        Re,mu,maxPicardIters,tolRes);

% ------------------------------------------------
% USER PARAM – under-relaxation factor  (0<ω≤1)
% typical 0.5–0.8 for high Re; 1.0 → pure Picard
omega = 1;
% ------------------------------------------------

% initial guess ---------------------------------------------------
U_current = [U1;U2;U3];
 U_current = update_bc(U_current,boundaryInfo,nodeInfo,...
                          Nxy,1,cornerIndex,bcFlags,inletProfile);
 iterpicard = 0;
% ================================================================
% Picard loop
% ================================================================
for iter = 1:maxPicardIters
    tNow = 1;  
    iterpicard = iterpicard +1;
    % dummy “time”
  
    % --- build matrix + extra Robin RHS --------------------------------
    [K_new] = build_stiffness_operator_stationary( ...
        nodeInfo,elemInfo,boundaryInfo, ...
        D,Re,o, ...
        U_current(1:Nxy),U_current(Nxy+1:2*Nxy),U_current(2*Nxy+1:end), ...
        gamma,mu,rho,cornerIndex);


    

    % --- right hand side (Dirichlet rows already handled inside K) -----
    rhs = update_bc( zeros(size(U_current)), ...
                     boundaryInfo,nodeInfo,Nxy,tNow, ...
                     cornerIndex,bcFlags,inletProfile );
    

    % --- linear solve --------------------------------------------------
    U_raw = K_new \ rhs;                      % Newton/Picard update
    

    % --- ω-relaxation --------------------------------------------------
    U_next = U_current + omega*(U_raw - U_current);
   
    % --- convergence check --------------------------------------------
    dU = norm(U_next - U_current);
    fprintf('  it %2d : ||ΔU|| = %.3e   (ω = %.2f)\n',iter,dU,omega);

    if dU < tolRes
        disp('  >>> converged');   U_current = U_next;   break;
    end
    U_current = U_next;                       % iterate
end

% -----------------------------------------------------------------
% 4)  pack remaining outputs (unchanged)
% -----------------------------------------------------------------
U        = U_current;
z1       = U(1:Nxy);       z2 = U(Nxy+1:2*Nxy);
z3       = U(2*Nxy+1:end); z  = hypot(z1,z2);

end

