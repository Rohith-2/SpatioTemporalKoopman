clear all
close all
clc

x = linspace(0,10,100);
t = linspace(0,100,100);
U = zeros(length(x));

varepsilon1 = 10^-8;
varepsilon2 = 10^-7;

for i=1:length(x)
    for j=1:length(t)
        X = x(i); T = t(j);
        U(i,j) = (0.5 + sin(X)) * (2*cos(2*X - (2*pi/45)*T) + 0.5*cos(10*X - sqrt(10)*T));
    end
end

contourf(U)

[Vreconst,Modes,Amplitudes,Amplitudesx,GrowthRatex,Frequencyx,Amplitudest,GrowthRatet,Frequencyt]=...
CalculateDMDdSdT(10,1,x,t,U,varepsilon1,varepsilon2);
figure

contourf(abs(Vreconst))
shading('interp')
colormap(jet)