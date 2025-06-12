function [K] = ...
          build_stiffness_operator_stationary( ...
              nodeInfo, elemInfo, boundaryInfo, ...
              D, Re, o, U1, U2, U3, gamma, mu, rho, corner)
%BUILD_STIFFNESS_OPERATOR_P2P1  (no-SUPG, no Jacobian)
%  Stationary Navier–Stokes stiffness matrix for P2–P1 triangles.
%
%  INPUT  --------
%    nodeInfo.velocity.x/y     : coords of P2 velocity nodes
%    nodeInfo.pressure.x/y     : coords of P1 pressure nodes
%    elemInfo.velElements      : (nEl×6)   velocity connectivity
%    elemInfo.presElements     : (nEl×3)   pressure connectivity
%    boundaryInfo.allNodes     : velocity DOFs with Dirichlet BC
%    D,Re,o, gamma, mu, rho    : physical/weak-form params
%    U1,U2 (size = #velNodes)  : current velocity guess
%    U3   (size = #presNodes)  : current pressure guess
%    corner                    : one pressure node to pin (p=0)
%
%  OUTPUT --------
%    K   global sparse stiffness (2*Nxy+Npr) × (2*Nxy+Npr)
%    Jx1 … Jc2  sub-blocks (optional)


% ------------------------------------------------------------------
% 0) basic sizes
% ------------------------------------------------------------------
numVel  = length(nodeInfo.velocity.x);    % P2 nodes
numPres = length(nodeInfo.pressure.x);    % P1 nodes
Nxy     = numVel;                         % velocity DOFs / comp.

nEl  = size(elemInfo.velElements,1);      % # triangles

% triplet pre-allocation (6×6=36, 6×3=18, 3×6=18 entries / element)
nP2P2 = 36*nEl;    nP2P1 = 18*nEl;    nP1P2 = 18*nEl; nP1P1 = 9*nEl;  

I_jx1 = zeros(nP2P2,1);  J_jx1 = I_jx1;  V_jx1 = I_jx1;
I_jx2 = I_jx1;           J_jx2 = I_jx1;  V_jx2 = I_jx1;
I_jx3 = zeros(nP2P1,1);  J_jx3 = I_jx3;  V_jx3 = I_jx3;

I_jy1 = I_jx1;  J_jy1 = I_jx1;  V_jy1 = I_jx1;
I_jy2 = I_jx1;  J_jy2 = I_jx1;  V_jy2 = I_jx1;
I_jy3 = I_jx3;  J_jy3 = I_jx3;  V_jy3 = I_jx3;

I_jc1 = zeros(nP1P2,1);  J_jc1 = I_jc1;  V_jc1 = I_jc1;
I_jc2 = I_jc1;           J_jc2 = I_jc1;  V_jc2 = I_jc1;

                % 3×3 entries per element
I_jc3 = zeros(nP1P1,1);  J_jc3 = zeros(nP1P1,1); V_jc3 = zeros(nP1P1,1);


ptr = struct('jx1',1,'jx2',1,'jx3',1,'jy1',1,'jy2',1,'jy3',1,'jc1',1,'jc2',1,'jc3',1);



% ------------------------------------------------------------------
% 1) Gauss & shape functions (reference triangle)
% ------------------------------------------------------------------
[NP2, dNxiP2, dNetaP2, wG, ~]   = precomputeShapeFunctionsP2_Tri();
[NP1, dNxiP1, dNetaP1, ~ ,  ~]  = precomputeShapeFunctionsP1_Tri();

nGp = length(wG);             % same rule for all sub-blocks

% ------------------------------------------------------------------
% 2) boundary velocity DOFs
% ------------------------------------------------------------------
bVel = boundaryInfo.allVelNodes(:);

% ------------------------------------------------------------------
% 3) element loop
% ------------------------------------------------------------------
% ---- choose a target wall shear once, before the element loop -----

% at top of build_stiffness_operator_stationary.m


for e = 1:nEl
    Kvel = elemInfo.velElements(e,:);   % 1×6
    Kpr  = elemInfo.presElements(e,:);  % 1×3

    xV = nodeInfo.velocity.x(Kvel);
    yV = nodeInfo.velocity.y(Kvel);
    xP = nodeInfo.pressure.x(Kpr);
    yP = nodeInfo.pressure.y(Kpr);




    % local copies of unknowns
    U1el = U1(Kvel);   U2el = U2(Kvel);   Pel = U3(Kpr);


    % ------------------------------------------------------------------
    % τ_m  (SUPG / momentum VMS)  and  τ_p  (PSPG)   --  P2-P1 version
    % ------------------------------------------------------------------
    % ❶ geometric length
    vec1    = [xV(2)-xV(1) ; yV(2)-yV(1)];
    vec2    = [xV(3)-xV(1) ; yV(3)-yV(1)];
    A_elem  = 0.5 * abs( vec1(1)*vec2(2) - vec1(2)*vec2(1) );   % scalar area
    h_lin  = 2*sqrt( A_elem / pi );         % “circle with same area”
    h      = h_lin / 3;                    % divide by (p+1)  (p = 2)

    % ❷ convective velocity magnitude  –– use element‐average you already have
    u_mag = sqrt( mean( U1el.^2 + U2el.^2 ) ) + 1e-12;  % add eps to avoid 0/0

    nu    = mu/(rho);
    % ❸ Tezduyar/Shakib τ formula for triangles
    C_I   = 3.0;                            % 3–4 for triangles, 1/√3 for quads
    inv2  = ( 2*u_mag / h )^2 + ( C_I*nu / h^2 )^2;
    tau_m = 1.0 / sqrt( inv2 );  % SUPG / VMS parameter
    



    % ❺ plug into your code
    taus1 = 1*tau_m;
    % taus1 = 0;
    taup = 0;


    % local 6×6 / 6×3 / 3×6 blocks
    Kx1e = zeros(6,6); Kx2e = zeros(6,6); Kx3e = zeros(6,3);
    Ky1e = zeros(6,6); Ky2e = zeros(6,6); Ky3e = zeros(6,3);
    Kc1e = zeros(3,6); Kc2e = zeros(3,6); Kc3e = zeros(3,3);
    
    % ------- gauss loop ------------------------------------------
    for gp = 1:nGp
         % --- geometry & first derivatives ------------------------------------
        [dNx, dNy, detJ, d2Nx2, d2Ny2] = p2ShapeDerivativesAllNodes( ...   
                                       xV, yV, ...
                                       dNxiP2(:,gp), dNetaP2(:,gp));

      
       % -------------------------------------------------------------
      

    % ---------------------------------------------------------------------
    % *** all the existing viscosity / convection code follows here ***
    % ---------------------------------------------------------------------
    %                                               
        % P1 derivatives for velocity-pressure
        [dNxp, dNyp, ~] = p1ShapeDerivativesAllNodes( ...
                                xP, yP, dNxiP1(:,gp), dNetaP1(:,gp));

        Ni   = NP2(:,gp);                      % 6×1
        Np   = NP1(:,gp);                      % 3×1
        wdet = wG(gp)*detJ;

        % element-average velocities & their derivatives
        a1 = Ni.'  * U1el;      a2 = Ni.'  * U2el;
       

           
        
        % Fill up matrix A: symmetric shear terms + non-linear velocity
        % terms + SUPG terms
        
       
                supg_test = a1*dNx + a2*dNy; lap_u_gp = (d2Nx2 + d2Ny2); lap_v_gp = (d2Nx2 + d2Ny2); shap_i = Ni;
                visc_xx = mu*(2*(dNx*dNx.') + dNy*dNy.') ...
                          + gamma*(dNx*dNx.') - taus1*mu*(supg_test*lap_u_gp.');
                visc_xy = mu*(dNy*dNx.') ...
                          + gamma*(dNx*dNy.');
                visc_yx = mu*(dNx*dNy.') ...
                          + gamma*(dNy*dNx.');
                visc_yy = mu*(dNx*dNx.' + 2*(dNy*dNy.')) ...
                          + gamma*(dNy*dNy.') - taus1*mu*(supg_test*lap_v_gp.');

                conv_x = rho*Re*(shap_i*(a1*dNx + a2*dNy).') + taus1*rho*Re*(supg_test*(a1*dNx + a2*dNy).');
                conv_y = rho*Re*(shap_i*(a1*dNx + a2*dNy).') + taus1*rho*Re*(supg_test*(a1*dNx + a2*dNy).');

                Kx1e = Kx1e + (visc_xx + conv_x)*wdet;
                Kx2e = Kx2e + (visc_xy          )*wdet;
                Ky1e = Ky1e + (visc_yx          )*wdet;
                Ky2e = Ky2e + (visc_yy + conv_y)*wdet;


        % Fill up matrix B^T: pressure terms + SUPG terms

                supg_test = a1*dNx + a2*dNy;
                Kx3e = Kx3e - (dNx*Np.' - taus1*(supg_test*(dNxp).')) * wdet;
                Ky3e = Ky3e - (dNy*Np.' - taus1*(supg_test*(dNyp).')) * wdet;

        
          % Fill up matrix B: continuity terms       
               
                Kc1e = Kc1e - (Np*dNx.')*wdet  ;
                Kc2e = Kc2e - (Np*dNy.')*wdet  ;
        
      
    end % gp

    % ------------------------------------------------------------------
    % accumulate in triplets
    % ------------------------------------------------------------------
    % helpers
    [i6,i6T] = ndgrid(Kvel,Kvel);   % 36
    range36  = 0:35;

    % Jx1
    idx = ptr.jx1 + range36;       ptr.jx1 = ptr.jx1 + 36;
    I_jx1(idx)=i6(:); J_jx1(idx)=i6T(:); V_jx1(idx)=Kx1e(:);

    % Jx2
    idx = ptr.jx2 + range36;       ptr.jx2 = ptr.jx2 + 36;
    I_jx2(idx)=i6(:); J_jx2(idx)=i6T(:); V_jx2(idx)=Kx2e(:);

    % Jy1
    idx = ptr.jy1 + range36;       ptr.jy1 = ptr.jy1 + 36;
    I_jy1(idx)=i6(:); J_jy1(idx)=i6T(:); V_jy1(idx)=Ky1e(:);

    % Jy2
    idx = ptr.jy2 + range36;       ptr.jy2 = ptr.jy2 + 36;
    I_jy2(idx)=i6(:); J_jy2(idx)=i6T(:); V_jy2(idx)=Ky2e(:);

    % velocity–pressure 6×3
    [iv,jp] = ndgrid(Kvel,Kpr);    range18 = 0:17;

    idx = ptr.jx3 + range18;       ptr.jx3 = ptr.jx3 + 18;
    I_jx3(idx)=iv(:); J_jx3(idx)=jp(:); V_jx3(idx)=Kx3e(:);

    idx = ptr.jy3 + range18;       ptr.jy3 = ptr.jy3 + 18;
    I_jy3(idx)=iv(:); J_jy3(idx)=jp(:); V_jy3(idx)=Ky3e(:);

    % continuity 3×6
    [ip,jv] = ndgrid(Kpr,Kvel);

    idx = ptr.jc1 + range18;       ptr.jc1 = ptr.jc1 + 18;
    I_jc1(idx)=ip(:); J_jc1(idx)=jv(:); V_jc1(idx)=Kc1e(:);

    idx = ptr.jc2 + range18;       ptr.jc2 = ptr.jc2 + 18;
    I_jc2(idx)=ip(:); J_jc2(idx)=jv(:); V_jc2(idx)=Kc2e(:);


    % Kpr = elemInfo.presElements(e,:)  % 1×3 global pressure dof indices
    [Ipairs,Jpairs] = ndgrid(Kpr,Kpr); % each is 3×3

    % linear indices into the triplet arrays
    range9 = 0:8;
    idx   = ptr.jc3 + range9;
    
    % scatter element‐matrix into the triplet buffers
    I_jc3(idx) = Ipairs(:);
    J_jc3(idx) = Jpairs(:);
    V_jc3(idx) = Kc3e(:);
    
    % advance the pointer
    ptr.jc3 = ptr.jc3 + 9;
end % element loop ends

% ------------------------------------------------------------------
% 4) build sparse sub-blocks once
% ------------------------------------------------------------------
Jx1 = sparse(I_jx1,J_jx1,V_jx1,Nxy,Nxy);
Jx2 = sparse(I_jx2,J_jx2,V_jx2,Nxy,Nxy);
Jx3 = sparse(I_jx3,J_jx3,V_jx3,Nxy,numPres);

Jy1 = sparse(I_jy1,J_jy1,V_jy1,Nxy,Nxy);
Jy2 = sparse(I_jy2,J_jy2,V_jy2,Nxy,Nxy);
Jy3 = sparse(I_jy3,J_jy3,V_jy3,Nxy,numPres);

Jc1 = sparse(I_jc1,J_jc1,V_jc1,numPres,Nxy);
Jc2 = sparse(I_jc2,J_jc2,V_jc2,numPres,Nxy);
% Jc3 = spalloc(numPres,numPres,0);
% assemble the global pressure‐pressure stiffness block
Jc3 = sparse( I_jc3, J_jc3, V_jc3, numPres, numPres );

% % ------------------------------------------------------------------
% 5) impose Dirichlet BCs (velocity rows only)
% ------------------------------------------------------------------
Jx1(bVel,:) = 0; 
Jx2(bVel,:) = 0;  Jx3(bVel,:) = 0;
Jy1(bVel,:) = 0;  Jy2(bVel,:) = 0;  Jy3(bVel,:) = 0;
Jx1(sub2ind(size(Jx1),bVel,bVel)) = 1;
Jy2(sub2ind(size(Jy2),bVel,bVel)) = 1;


% % ------------------------------------------------------------------
% 6) fix one pressure DOF to zero
% ------------------------------------------------------------------
Jc3(corner,corner) = 1;



% ------------------------------------------------------------------
% 7) build global K
% ------------------------------------------------------------------
K = [Jx1 Jx2 Jx3;
     Jy1 Jy2 Jy3;
     Jc1 Jc2 Jc3];




end

% ------------------------------------------------------------------
% ------------------------------------------------------------------
function [dNxAll,dNyAll,detJ,d2Nx2,d2Ny2] = p2ShapeDerivativesAllNodes( ...
           x,y,dNxi,dNeta)
% 6-node quadratic triangle – returns first *and* second derivatives
% (second rows identical to what you already used in the residual builder)

% Jacobian and its inverse
dX_dxi  = sum(x.*dNxi);   dX_deta = sum(x.*dNeta);
dY_dxi  = sum(y.*dNxi);   dY_deta = sum(y.*dNeta);
J       = [dX_dxi dY_dxi; dX_deta dY_deta];
detJ    = det(J);
invJ    = inv(J);

% first derivatives wrt x,y
dNxAll  = invJ(1,1)*dNxi + invJ(1,2)*dNeta;
dNyAll  = invJ(2,1)*dNxi + invJ(2,2)*dNeta;

% reference-element second derivatives (constant)
d2N_xi2_ref   = [ 4;  4; 0; -8; 0;  0];
d2N_eta2_ref  = [ 4;  0; 4;  0; 0; -8];
d2N_xieta_ref = [ 4;  0; 0; -4; 4; -4];

a = invJ(1,1);  b = invJ(1,2);
c = invJ(2,1);  d = invJ(2,2);

% Laplacian bits for every local node
d2Nx2 = a*a*d2N_xi2_ref  + 2*a*b*d2N_xieta_ref + b*b*d2N_eta2_ref;
d2Ny2 = c*c*d2N_xi2_ref  + 2*c*d2N_xieta_ref   + d*d*d2N_eta2_ref;
end
% ------------------------------------------------------------------

% ------------------------------------------------------------------
function [dNxAll,dNyAll,detJ] = p1ShapeDerivativesAllNodes( ...
           x,y,dNxi,dNeta)
% 3-node linear triangle
dX_dxi  = sum(x.*dNxi);   dX_deta = sum(x.*dNeta);
dY_dxi  = sum(y.*dNxi);   dY_deta = sum(y.*dNeta);
J = [dX_dxi dY_dxi; dX_deta dY_deta];     detJ = det(J);
invJ = inv(J);
dNxAll = invJ(1,1)*dNxi + invJ(1,2)*dNeta;
dNyAll = invJ(2,1)*dNxi + invJ(2,2)*dNeta;
end
