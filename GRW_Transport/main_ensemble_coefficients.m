%% Ensemble dispersion coefficients computed with the 2D unbiased GRW algorithm for advection-diffusion in random fields
close all; clear all;

tic

ifdm=1; % expicit FDM
% ifdm=2; % impicit FDM

N=10^24;
R= 100;
Lx=199; Ly=99; Tt=10;
D1=0.01; D2=D1; d1=2; d2=d1;
stepU=5; U_MEAN=0.7134;
dx=0.1; dy=dx;
dt=stepU*dx/U_MEAN;
T=round(Tt/dt);
rx=2*D1*dt/(d1^2*dx^2);
ry=2*D2*dt/(d2^2*dy^2);
r=rx+ry; % rx+ry<=1 !
n=zeros(Lx,Ly); nn=zeros(Lx,Ly);
Dx=zeros(1,T); Mx=zeros(1,T); Varx=zeros(1,T);u=zeros(Lx,Ly);v=zeros(Lx,Ly);
Dy=zeros(1,T); My=zeros(1,T); Vary=zeros(1,T);
restr=0; restjump=0; restjumpx=0; restjumpy=0;
x0=round(Lx*dx/10);
y0=round(Ly*dy/2);
iIC=1;
if iIC==0
    n(x0/dx,y0/dy)=N;
    n0=n; % Dirac IC
else
    ti=0.01;
    gss=Gauss_IC(ti,dx,dy,x0,y0,Lx,Ly,U_MEAN,D1);
    n0 = N*gss/sum(sum(gss)); % Gaussian IC
end
tp=round(T/3)-1; tsel=tp:tp:T;
X=1:1:Lx; Y=1:1:Ly;
xr=(X*dx-x0); yr=(Y*dy-y0);
Mx0=sum(xr*n0)/N; Varx0=sum(xr.^2*n0)/N-Mx0*Mx0;
My0=sum(yr*n0')/N; Vary0=sum(yr.^2*n0')/N-My0*My0;
Varx0=Varx0-Mx0*Mx0; Vary0=Vary0-My0*My0;
meanVx=0; meanVy=0;
for k=1:R
    if ifdm==1
        load(['..\FDM_Flow\test_explicitFDM','\dataGAUSS_',num2str(k),'.mat']);
    else
        load(['..\FDM_Flow\test_implicitFDM','\dataGAUSS_',num2str(k),'.mat']);
    end
    meanVx=meanVx+mean(mean(Vx)); meanVy=meanVy+mean(mean(Vy));
    U=Vx-U_MEAN; V=Vy;
    if mod(k,10)==0 && k<100
        disp(k)
    end
    for i=1:Lx
        for j=1:Ly
            ur=U(i,j); vr=V(i,j);
            u(i,j)=floor((ur+1)*stepU+0.5);
            v(i,j)=floor(vr*stepU+0.5);
        end
    end
    n=n0;
    for t=1:T
        ntot=0;
        for x=d1+1:(Lx-d1-1)
            for y=d2+1:(Ly-d2-1)
                if n(x,y) > 0
                    xa=x+u(x,y); ya=y+v(x,y);
                    restr=n(x,y)*(1-r)+restr; nsta=floor(restr);
                    restr=restr-nsta; njump=n(x,y)-nsta;
                    nn(xa,ya)=nn(xa,ya)+nsta;
                    restjump=njump*rx/r+restr;
                    njumpx=floor(restjump); restjump=restjump-njumpx;
                    njumpy=njump-njumpx;
                    if(njumpx)>0
                        restjumpx=njumpx/2+restjumpx;
                        nj(1)=floor(restjumpx); restjumpx=restjumpx-nj(1);
                        nj(2)=njumpx-nj(1);
                        for i=1:2
                            if nj(i)>0
                                xd=xa+(2*i-3)*d1;
                                nn(xd,ya)=nn(xd,ya)+nj(i);
                            end
                        end
                    end
                    if(njumpy)>0
                        restjumpy=njumpy/2+restjumpy;
                        nj(1)=floor(restjumpy); restjumpy=restjumpy-nj(1);
                        nj(2)=njumpy-nj(1);
                        for i=1:2
                            if nj(i)>0
                                yd=ya+(2*i-3)*d2;
                                nn(xa,yd)=nn(xa,yd)+nj(i);
                            end
                        end
                    end
                end
            end
        end
        for x=d1:(Lx-d1)
            for y=d2:(Ly-d2)
                n(x,y)=nn(x,y); nn(x,y)=0; ntot=ntot+n(x,y);
                xr=x*dx-x0; yr=y*dy-y0;
                Mx(t)=Mx(t)+xr*n(x,y); Varx(t)=Varx(t)+xr^2*n(x,y);
                My(t)=My(t)+yr*n(x,y); Vary(t)=Vary(t)+yr^2*n(x,y);
            end
        end
        if k==1 && sum(ismember(t,tsel))
            [lin,col]=find(n>0);
            [nl,nc]=size(n);
            linind=sub2ind([nl,nc],lin,col);
            figure(1); plot3(lin*dx,col*dx,n(linind));
            if t==tsel(1);hold
            end
            grid on; view(40,20)
            xlabel('x'); ylabel('y'); zlabel('n(x,y)'); xlim([0 12]);
            set(gca, 'XTick', [0:2:12]); set(gca, 'YTick', [0:2:10])
        end

    end
    if mod(k,100)==0
        disp(k)
    end
end
for t=1:T
    Mx(t)=Mx(t)/ntot/R; My(t)=My(t)/ntot/R;
    Varx(t)=Varx(t)/ntot/R-Mx(t)*Mx(t)-Varx0;
    Vary(t)=Vary(t)/ntot/R-My(t)*My(t)-Vary0;
    Dx(t)=Varx(t)/(2.0*t*dt); Dy(t)=Vary(t)/(2.0*t*dt);
    Mx(t)=Mx(t)/(t*dt); My(t)=My(t)/(t*dt);
end
meanVx=meanVx/R
meanVy=meanVy/R
Pe=U_MEAN*dx/D1
t=1:T;
if ifdm==1
    save ('ensemble_coefficients_explicitFEM','t','dt','Mx','My','Dx','Dy')
else
    save ('ensemble_coefficients_implicitFEM','t','dt','Mx','My','Dx','Dy')
end

figure;
subplot(1,2,1)
plot(t*dt,Mx,t*dt,My,'LineWidth',1.1); ylim([-0.4 1.2]);
legend('V_{x}(t)','M_{y}(t)','Location','best'); legend('boxoff');
xlabel('t');
subplot(1,2,2)
plot(t*dt,Dx,t*dt,Dy,'LineWidth',1.1);
legend('D_{x}(t)','D_{y}(t)','Location','best'); legend('boxoff');
xlabel('t');
figure;
[lin,col]=find(n>0);
[nl,nc]=size(n);
linind=sub2ind([nl,nc],lin,col);
plot3(lin,col,n(linind));
grid on; view(40,20)
xlabel('x'); ylabel('y'); zlabel('n(x,y)');

toc
%% Elapsed time is 2.575515 seconds.

