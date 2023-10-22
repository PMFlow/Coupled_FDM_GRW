function [p,Vx,Vy] = explicitFDM(Nmod,varK,phi,wavenum,KMean,I,J,dx,dy,x,y,x2,y2,p0)
%% 2D explicit FDM flow-solver for realizations of random velocity fields
%  hydr. conductivity K is a Kraichnan field with Nmod cosine mods defined by wawenumbers 'wavenum' and phases 'phi'
%  2D-Homogenuous Stationary BVP with Variable Coefficients and Dirichlet-Neumann BC

close all

%%   Assigning Parameters
C1 = wavenum(:,1);
C2 = wavenum(:,2);
%%   Boundary Conditions
%homogeneous:
BCXL = @(v) 1 ;
BCXR = @(v) 0 ;
BCYB = @(u) 0 ;
BCYU = @(u) 0 ;
%%  Coefficients and Forcing Term
Dx = K(x2(1:I-1),y(2:J-1),Nmod,KMean,varK,C1,C2,phi,I-1,J-2); 
Dy = K(x(2:I-1),y2(1:J-1),Nmod,KMean,varK,C1,C2,phi,I-2,J-1); 
D  = K(x(2:I-1),y(2:J-1),Nmod,KMean,varK,C1,C2,phi,I-2,J-2);  

%% solution
%%%% Parameters
S=1e5; 
maxr=0.8;
Tolerance=1e-6; 

dt=8*(max(max(Dx))/dx^2+max(max(Dy))/dy^2); dt=maxr*1/dt
rx=dt*Dx/dx^2; ry=dt*Dy/dy^2;   % r<=1/4 such that 1-4*r>0 !!!
rloc=1-(rx(1:I-2,:)+rx(2:I-1,:)+ry(:,1:J-2)+ry(:,2:J-1));
%%%%% explicit iterative flow-solver
p = p0; pa=p;
pp=zeros(I,J); eps=zeros(1,S); 
for s=1:S
    pp(2:I-1,2:J-1)=rloc.*p(2:I-1,2:J-1) ...
        +rx(1:I-2,:).*p(1:I-2,2:J-1)+rx(2:I-1,:).*p(3:I,2:J-1) ...
        +ry(:,1:J-2).*p(2:I-1,1:J-2)+ry(:,2:J-1).*p(2:I-1,3:J);
    %%%% boundary conditions
    %%%% BCX Left/Right
    pp(1,:)=p0(1,:); pp(I,:)=p0(I,:);
    %%%% BCY Bottom/Up
    derB=BCYB(x(2:I-1));
    pp(2:I-1,1)=pp(2:I-1,2)-derB*dy;
    derU=BCYU(x(2:I-1));
    pp(2:I-1,J)=pp(2:I-1,J-1)+derU*dy;    
    p=pp;
    eps(s)=dx*norm(p-pa)+norm(p-pa)/norm(p);
    if eps(s) <= Tolerance
        fprintf('Number of iterations = %d Error = %e\n',s,eps(s)) ;
        break
    end
    pa=p;   
end
%%%% Velocities
Vx=-(p(3:I,2:J-1)-p(1:I-2,2:J-1)).*D/(2*dx);
Vy=-(p(2:I-1,3:J)-p(2:I-1,1:J-2)).*D/(2*dy);
