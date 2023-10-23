function sol = K(x,y,Nmod,KMean,varK,C1,C2,phi,Nx,Ny)
%% computes a realization of the Kraichnan fild K

coeff = sqrt(varK*2/Nmod) ;       % amplitude
sol = zeros(Nx,Ny) ;
for i = 1 : Nx
    for j = 1 : Ny
        phase = 2*pi*(C1*x(i)+C2*y(j))+phi;
        ak = sum(coeff.*cos(phase));
        sol(i,j) = KMean*exp(-varK/2)*exp(ak); %hydraulic conductivity:
    end
end
end