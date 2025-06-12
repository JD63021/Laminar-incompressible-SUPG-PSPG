function [K] = build_stiffness_operator_stationaryP1( ...
          nodeInfo, elemInfo, boundaryInfo, ...
          D, Re, o, U1, U2, U3, gamma, mu, rho, corner)
%BUILD_STIFFNESS_OPERATOR_P1P1  PSPG+SUPG‐stabilized stationary Navier–Stokes for P1–P1 triangles
%
%  INPUT
%    nodeInfo.velocity.x/y    : coords of P1 velocity nodes
%    nodeInfo.pressure.x/y    : coords of P1 pressure nodes
%    elemInfo.velElements     : (nEl×3) velocity connectivity
%    elemInfo.presElements    : (nEl×3) pressure connectivity
%    boundaryInfo.allVelNodes : Dirichlet velocity DOFs
%    D,Re,o,gamma,mu,rho      : physical / weak‐form params
%    U1,U2  (#vel nodes)      : velocity guess
%    U3     (#pres nodes)     : pressure guess
%    corner                   : one pressure node to pin (p=0)
%
%  OUTPUT
%    K  global sparse stiffness (2*Nvel + Npres) × (2*Nvel + Npres)

% 0) sizes & preallocate triplets
numVel  = length(nodeInfo.velocity.x);  % P1 nodes
numPres = length(nodeInfo.pressure.x);
Nxy     = numVel;
nEl     = size(elemInfo.velElements,1);

nPP = 9*nEl;  % 3×3 per element
I_jx1 = zeros(nPP,1); J_jx1 = I_jx1; V_jx1 = I_jx1;
I_jx2 = I_jx1;  J_jx2 = I_jx1; V_jx2 = I_jx1;
I_jx3 = zeros(nPP,1); J_jx3 = I_jx3; V_jx3 = I_jx3;
I_jy1 = I_jx1;  J_jy1 = I_jx1; V_jy1 = I_jx1;
I_jy2 = I_jx1;  J_jy2 = I_jx1; V_jy2 = I_jx1;
I_jy3 = I_jx3;  J_jy3 = I_jx3; V_jy3 = I_jx3;
I_jc1 = zeros(nPP,1); J_jc1 = I_jc1; V_jc1 = I_jc1;
I_jc2 = I_jc1;           J_jc2 = I_jc1; V_jc2 = I_jc1;
I_jc3 = zeros(nPP,1); J_jc3 = I_jc3; V_jc3 = I_jc3;

ptr = struct('jx1',1,'jx2',1,'jx3',1,'jy1',1,'jy2',1,'jy3',1,'jc1',1,'jc2',1,'jc3',1);

% 1) shape functions on reference triangle (linear)
[NP1, dNxiP1, dNetaP1, wG, ~] = precomputeShapeFunctionsP1_Tri();
nGp = length(wG);

% 2) Dirichlet velocity rows
bVel = boundaryInfo.allVelNodes(:);

% 3) element loop
for e = 1:nEl
    Kvel = elemInfo.velElements(e,:);   % 1×3
    Kpr  = elemInfo.presElements(e,:);  % 1×3

    xV = nodeInfo.velocity.x(Kvel);
    yV = nodeInfo.velocity.y(Kvel);
    xP = nodeInfo.pressure.x(Kpr);
    yP = nodeInfo.pressure.y(Kpr);

    U1el = U1(Kvel);  U2el = U2(Kvel);  Pel = U3(Kpr);

    % compute stabilization parameters (p=1)
    vec1   = [xV(2)-xV(1); yV(2)-yV(1)];
    vec2   = [xV(3)-xV(1); yV(3)-yV(1)];
    A_elem = 0.5*abs(vec1(1)*vec2(2) - vec1(2)*vec2(1));
    h_lin  = 2*sqrt(A_elem/pi);
    h      = h_lin/2;                % divide by (p+1), p=1
    u_mag  = sqrt(mean(U1el.^2 + U2el.^2)) + 1e-12;
    nu     = mu/rho;
    C_I    = 3;
    inv2   = (2*u_mag/h)^2 + (C_I*nu/h^2)^2;
    tau_m  = 1.0/sqrt(inv2);

    taus1  = tau_m; % tau for SUPG (turn on)
    taup   = tau_m; % tau for PSPG
    % taus1 =0;         % turn off SUPG
    % taup = 0;       % turn off PSPG

    % init local blocks
    Kx1e = zeros(3,3); Kx2e = zeros(3,3); Kx3e = zeros(3,3);
    Ky1e = zeros(3,3); Ky2e = zeros(3,3); Ky3e = zeros(3,3);
    Kc1e = zeros(3,3); Kc2e = zeros(3,3); Kc3e = zeros(3,3);

    % Gauss loop
    for gp = 1:nGp
        % geometry & derivatives (P1)
        [dNx, dNy, detJ] = p1ShapeDerivativesAllNodes( ...
                              xV, yV, dNxiP1(:,gp), dNetaP1(:,gp));
        [dNxp, dNyp, ~]  = p1ShapeDerivativesAllNodes( ...
                              xP, yP, dNxiP1(:,gp), dNetaP1(:,gp));
        Ni   = NP1(:,gp);               % 3×1
        Np   = NP1(:,gp);
        wdet = wG(gp)*detJ;

        % velocity approx & grads
        a1   = Ni.'*U1el;  a2   = Ni.'*U2el;
        % a1x, a1y not needed except for laplacian
        

        

         supg_test = a1*dNx + a2*dNy; 
                visc_xx = mu*(2*(dNx*dNx.') + dNy*dNy.') ...
                          + gamma*(dNx*dNx.') ;
                visc_xy = mu*(dNy*dNx.') ...
                          + gamma*(dNx*dNy.');
                visc_yx = mu*(dNx*dNy.') ...
                          + gamma*(dNy*dNx.');
                visc_yy = mu*(dNx*dNx.' + 2*(dNy*dNy.')) ...
                          + gamma*(dNy*dNy.') ;

                conv_x = Re*(Ni*(a1*dNx + a2*dNy).') + taus1*Re*(supg_test*(a1*dNx + a2*dNy).');
                conv_y = Re*(Ni*(a1*dNx + a2*dNy).') + taus1*Re*(supg_test*(a1*dNx + a2*dNy).');

                Kx1e = Kx1e + (visc_xx + conv_x)*wdet;
                Kx2e = Kx2e + (visc_xy          )*wdet;
                Ky1e = Ky1e + (visc_yx          )*wdet;
                Ky2e = Ky2e + (visc_yy + conv_y)*wdet;


        % Fill up matrix B^T: pressure terms + SUPG terms

                supg_test = a1*dNx + a2*dNy; 
                Kx3e = Kx3e - (dNx*Np.' - taus1*(supg_test*(dNxp).')) * wdet;
                Ky3e = Ky3e - (dNy*Np.' - taus1*(supg_test*(dNyp).')) * wdet;

        
        % Fill up matrix B: continuity terms + PSPG       
                pspg_x = dNxp;  pspg_y = dNyp;
                Kc1e = Kc1e - (Np*dNx.')*wdet  - taup*(pspg_x* (Re*(a1*dNx + a2*dNy).' ))*wdet;
                Kc2e = Kc2e - (Np*dNy.')*wdet  - taup*(pspg_y* (Re*(a1*dNx + a2*dNy).' ))*wdet;
                
        % Fill up null matrix: PSPG  
                Kc3e = Kc3e - taup*((pspg_x*dNxp.') + (pspg_y*dNyp.'))*wdet;
        % 
        % % continuity + PSPG
        % for p = 1:3
        %     pspg_x = dNxp(p);  pspg_y = dNyp(p);
        %     for k = 1:3
        %         dNxk = dNx(k);  dNyk = dNy(k);
        %         Kc1e(p,k) = Kc1e(p,k) + dNx(k)*NP1(p,gp)*wdet ...
        %                     + taup*pspg_x*(Re*(a1*dNxk + a2*dNyk) )*wdet;
        %         Kc2e(p,k) = Kc2e(p,k) + dNy(k)*NP1(p,gp)*wdet ...
        %                     + taup*pspg_y*(Re*(a1*dNxk + a2*dNyk) )*wdet;
        %     end
        %     for k = 1:3
        %         dNxp_k = dNxp(k); dNyp_k = dNyp(k);
        %         Kc3e(p,k) = Kc3e(p,k) + taup*(dNxp_k*pspg_x + dNyp_k*pspg_y)*wdet;
        %     end
        % end
    end  % gp

    % scatter element‐matrices into triplets
    [i3,i3T] = ndgrid(Kvel,Kvel);   % 9 entries
    idx = ptr.jx1 + (0:8); I_jx1(idx)=i3(:); J_jx1(idx)=i3T(:); V_jx1(idx)=Kx1e(:);
    ptr.jx1 = ptr.jx1+9;
    idx = ptr.jx2 + (0:8); I_jx2(idx)=i3(:); J_jx2(idx)=i3T(:); V_jx2(idx)=Kx2e(:);
    ptr.jx2 = ptr.jx2+9;

    idx = ptr.jy1 + (0:8); I_jy1(idx)=i3(:); J_jy1(idx)=i3T(:); V_jy1(idx)=Ky1e(:);
    ptr.jy1 = ptr.jy1+9;
    idx = ptr.jy2 + (0:8); I_jy2(idx)=i3(:); J_jy2(idx)=i3T(:); V_jy2(idx)=Ky2e(:);
    ptr.jy2 = ptr.jy2+9;

    [iv,jp] = ndgrid(Kvel,Kpr);  % 9 entries
    idx = ptr.jx3 + (0:8); I_jx3(idx)=iv(:); J_jx3(idx)=jp(:); V_jx3(idx)=Kx3e(:);
    ptr.jx3 = ptr.jx3+9;
    idx = ptr.jy3 + (0:8); I_jy3(idx)=iv(:); J_jy3(idx)=jp(:); V_jy3(idx)=Ky3e(:);
    ptr.jy3 = ptr.jy3+9;

    [ip,jv] = ndgrid(Kpr,Kvel);
    idx = ptr.jc1 + (0:8); I_jc1(idx)=ip(:); J_jc1(idx)=jv(:); V_jc1(idx)=Kc1e(:);
    ptr.jc1 = ptr.jc1+9;
    idx = ptr.jc2 + (0:8); I_jc2(idx)=ip(:); J_jc2(idx)=jv(:); V_jc2(idx)=Kc2e(:);
    ptr.jc2 = ptr.jc2+9;

    [Ip,Jp] = ndgrid(Kpr,Kpr);
    idx = ptr.jc3 + (0:8); I_jc3(idx)=Ip(:); J_jc3(idx)=Jp(:); V_jc3(idx)=Kc3e(:);
    ptr.jc3 = ptr.jc3+9;
end  % element

% 4) build sparse sub‐blocks
Jx1 = sparse(I_jx1,J_jx1,V_jx1,Nxy,Nxy);
Jx2 = sparse(I_jx2,J_jx2,V_jx2,Nxy,Nxy);
Jx3 = sparse(I_jx3,J_jx3,V_jx3,Nxy,numPres);
Jy1 = sparse(I_jy1,J_jy1,V_jy1,Nxy,Nxy);
Jy2 = sparse(I_jy2,J_jy2,V_jy2,Nxy,Nxy);
Jy3 = sparse(I_jy3,J_jy3,V_jy3,Nxy,numPres);
Jc1 = sparse(I_jc1,J_jc1,V_jc1,numPres,Nxy);
Jc2 = sparse(I_jc2,J_jc2,V_jc2,numPres,Nxy);
Jc3 = sparse(I_jc3,J_jc3,V_jc3,numPres,numPres);

% 5) impose Dirichlet BCs on velocity rows
Jx1(bVel,:)=0; Jx2(bVel,:)=0; Jx3(bVel,:)=0; Jy1(bVel,:)=0; Jy2(bVel,:)=0; Jy3(bVel,:)=0;
Jx1(sub2ind(size(Jx1),bVel,bVel))=1;
Jy2(sub2ind(size(Jy2),bVel,bVel))=1;

% 6) pin one pressure DOF
Jc3(corner,corner)=1;

% 7) assemble global K
K = [Jx1 Jx2 Jx3;
     Jy1 Jy2 Jy3;
     Jc1 Jc2 Jc3];
end


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
