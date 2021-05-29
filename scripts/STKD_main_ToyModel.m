clc
clear all
close all

addpath '/Volumes/Rohith/College/SpatioTemporalKoopman/HODMD'

%load V.mat
load sent_snp.mat
%load momo.mat
V = sent;
%load Time.mat
%load Exis.mat
s = size(V);
Time = linspace(0,1,s(2));
Exis = linspace(0,1,s(1));


% APPLY STKD TO V
[J,K]=size(V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% STKD  PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SET PARAMETERS STKD
% Tolerances
% SVD
varepsilon1=1e-8;
% DMD
varepsilon2=1e-7;
% Set index d: dSpace for spatial analysis and dTime for temporal analysis
dSpace=2;
dTime=1;

[Vreconst,Modes,Amplitudes,Amplitudesx,GrowthRatex,Frequencyx,Amplitudest,GrowthRatet,Frequencyt]=...
CalculateDMDdSdT(dSpace,dTime,Time,Exis,V,varepsilon1,varepsilon2);

% CALCULATE RRMS ERROR
dif=Vreconst-V;
errorRMS=norm(dif(:))/norm(V(:));
disp('Reconstruction error: RRMSE')
disp(errorRMS)


% PLOT RECONSTRUCTION
h9=figure;
axes9 = axes('Parent',h9,'FontSize',14,'FontName','Agency FB');
box(axes9,'on');
pcolor(abs(Vreconst))
shading('interp')
colormap(jet)
xlabel('Date')
ylabel('Sentiment')

figure
pcolor(V)
shading('interp')
colormap(jet)
xlabel('Date')
ylabel('Sentiment')
title("OG Data")

% Plot: frequency/wavenumber vs. amplitude, vs. growth rate
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
box(axes1,'on');
semilogy(Frequencyt,GrowthRatet,'o','linewidth',2,'color','k','MarkerSize',8);
semilogy(Frequencyx,GrowthRatex,'o','linewidth',2,'color','b','MarkerSize',8);
set(axes1,'YMinorTick','on');
xlabel('Frequency/Wavenumber')
ylabel('Growth Rate')


figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
box(axes1,'on');
semilogy(Frequencyt,Amplitudest,'+','linewidth',2,'color','k','MarkerSize',8);
semilogy(Frequencyx,Amplitudesx,'+','linewidth',2,'color','b','MarkerSize',8);
set(axes1,'YMinorTick','on','YScale','log');
xlabel('Frequency/Wavenumber')
ylabel('Amplitude')