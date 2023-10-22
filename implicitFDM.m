function [uN,speedX,speedY,x,y,Nx,Ny] = implicitFDM(Nmod,varK,phiVect,wavenum0,wavenum1,KMean,Nx,Ny,dx,dy,x,y)
%% 2D implicit FDM flow-solver for realizations of random velocity fields
%  hydr. conductivity K is a Kraichnan field with Nmod cosine mods defined by wawenumbers 'wavenum' and phases 'phi'
%  2D-Homogenuous Stationary BVP with Variable Coefficients and Dirichlet-Neumann BC

%%   Variable Initialization
phi = phiVect(1:Nmod) ;
C1 = wavenum0(1:Nmod) ;
C2 = wavenum1(1:Nmod) ;
%%   Matrix Initializations
f = zeros(Nx,Ny) ;
CM1 = K(x + dx/2,y,Nmod,KMean,varK,C1,C2,phi,Nx,Ny) ;
CM2 = K(x,[y(1) - dy, y] + dy/2,Nmod,KMean,varK,C1,C2,phi,Nx,Ny+1) ;
CM_Realizari_X = K(x(2:end-1),y,Nmod,KMean,varK,C1,C2,phi,Nx-2,Ny) ;

A = CM1(1:end-2,:)*(1/dx^2) ;
D = CM1(2:end-1,:)*(1/dx^2) ;

B = CM2(2:end-1,1:end-1)*(1/dy^2) ;
E = CM2(2:end-1,2:end)*(1/dy^2) ;

C = (-1) * (A+B+D+E) ;
%%   The Left Hand Side of the Numerical Solution
Diag0 = C(:) ;

auxLW = A(2:end,:) ;
auxLW = [ auxLW ; zeros(1,Ny) ] ;
DiagLW = auxLW(:) ;
DiagLW = DiagLW(1:end-1) ;

auxUP = D(1:end-1,:) ;
auxUP = [ auxUP ; zeros(1,Ny) ] ;
DiagUP = auxUP(:) ;
DiagUP = DiagUP(1:end-1) ;

DiagUP2 = E(:,1:end-1) ;
DiagUP2 = DiagUP2(:) + [ B(:,1) ; zeros(((Nx-2)*Ny-(Nx-2))-length(B(:,1)),1)] ;

DiagLW2 = B(:,2:end);
DiagLW2 = DiagLW2(:) + [zeros(((Nx-2)*Ny-(Nx-2))-length(E(:,end)),1) ; E(:,end)] ;

LHS = sparse((Nx-2)*Ny,(Nx-2)*Ny) ;
LHS = spdiags(Diag0,0,LHS) ;                   % main diagonal
LHS = spdiags(DiagLW,-1,LHS) ;                 % lower diagonal
LHS = spdiags(DiagLW2,-(Nx-2),LHS) ;           % -(Nx-2) diagonal

DiagUP = [0 DiagUP']';                         % upper diagonal
LHS = spdiags(DiagUP,1,LHS) ;

DiagUP2 = [zeros(Nx-2,1) ; DiagUP2] ;          % (Nx-2) upper diagonal
LHS = spdiags(DiagUP2,Nx-2,LHS) ;
%%   Boundary Conditions
BCXL = @(v) 1 ;
BCXR = @(v) 0 ;

BCYB = @(u) 0 ;
BCYU = @(u) 0 ;

uN = zeros(Nx,Ny) ;
uN(1,:) = BCXL(y) ;
uN(Nx,:) = BCXR(y) ;

%%   The Right Hand Side of the Matrix Solution

vRHS_F = f(2:end-1,:) ;
vRHS_F = vRHS_F(:) ;


auxRHS_A = [ A(1,:).*BCXL(y)  ; zeros(Nx-3,length(A(1,:))) ] ;
auxRHS_A = auxRHS_A(:) ;
vRHS_A = (-1) *auxRHS_A ;

auxRHS_D = [ zeros(Nx-3,length(D(end,:))) ; D(end,:).*BCXR(y) ] ;
vRHS_D = (-1) * auxRHS_D(:) ;

vRHS_Start = [2*dy * (BCYB(x(2:Nx-1)) .* B(:,1)')  zeros(1,(Nx-2)*Ny-(Nx-2))]';
vRHS_End = [zeros(1,(Nx-2)*Ny-(Nx-2))  (-2)*dy * (BCYU(x(2:Nx-1)) .* E(:,Ny)')]';
RHS = vRHS_F + vRHS_A + vRHS_D + vRHS_Start + vRHS_End ;


%%   Numerical Solution

rez = (LHS\RHS)' ;
rez = reshape(rez,Nx-2,Ny) ;
uN(2:Nx-1,:) = rez ;

speedX = zeros(Nx-2,Ny-2) ;
speedY = zeros(Nx-2,Ny-2) ;

for j = 2 : Ny - 1

    uN_SpeedX = (- 1/(2*dx)) * (CM_Realizari_X(:,j) .* (uN(3:Nx,j) - uN(1:Nx-2,j))) ;
    uN_SpeedY = (- 1/(2*dy)) * (CM_Realizari_X(:,j) .* (uN(2:Nx-1,j+1) - uN(2:Nx-1,j-1))) ;

    speedX(:,j-1) = uN_SpeedX ;
    speedY(:,j-1) = uN_SpeedY ;

end
end