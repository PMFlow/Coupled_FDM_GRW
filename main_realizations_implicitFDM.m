%% Computes ensembles of random velocity fields with the implicit FDM 

clear ;
close all ;

tic ;
global state;
initstate; state;

Nmod = 10; %% Nmod = 10^2 ; % for comparison with the 'linear_approximation'
NRealiz = 1; % 10^2; % %% NRealiz = 10^4 ; % for comparison with the 'linear_approximation'
varK = 0.1;
ZC1 = 1.0;
ZC2 = 1.0;
KMean = 15;
Nx=201; Ny=101;
a = 0 ;
b = 20 ;
c = 0 ;
d = 10 ;
dx = (b-a)/(Nx-1) ;
dy = (d-c)/(Ny-1) ;
x = a + ( (1:Nx)-1 )*dx ;
y = c + ( (1:Ny)-1 )*dy ;

%% GAUSSIAN CORRELATION

for m = 1 : length(varK)
    for n = 1 : NRealiz
        [wavenum, phi] = Kraichnan_Gauss_param(Nmod,ZC1,ZC2);
        [p,Vx,Vy,x,y,Nx,Ny] = implicitFDM(Nmod,varK(m),phi,wavenum(:,1),wavenum(:,2),KMean,Nx,Ny,dx,dy,x,y) ;
        if mod(n,1000)==0
            fprintf('Number of realization : %d \n',n) ;
        end
        % save(['implicitFDM','\dataGAUSS_',num2str(n),'.mat'],'p','Vx','Vy') ;
    end
end

toc ;
%% Elapsed time is 0.098 seconds / 1 realization (= 11.60 / 100 realiations)
