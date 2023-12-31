function [wavenum, phi] = Kraichnan_Gauss_param(NMOD,ZC1,ZC2)
%% computes phases and wavenumbers of the Kraichnan fild K

sq2Pi=1/sqrt(2)/pi;
wavenum(:,1)=sq2Pi*randn(NMOD,1)/ZC1;
wavenum(:,2)=sq2Pi*randn(NMOD,1)/ZC2;
phi = 2.*pi*rand(NMOD,1) ;
