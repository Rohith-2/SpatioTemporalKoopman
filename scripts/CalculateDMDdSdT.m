
function  [Vreconst,Modes,Amplitudes,Amplitudesx,GrowthRatex,Frequencyx,Amplitudest,GrowthRatet,Frequencyt]...
    =CalculateDMDdSdT(d1,d2,Tiempos,Exis,V,varepsilon1,varepsilon2)

% dimension reduction: SVD
[U,Sigma,T]=svd(V);
sigmas=diag(Sigma);
n=length(sigmas);

SigmNorm=norm(sigmas,2);
kk=0;
for k=1:n
    if norm(sigmas(k:n),2)/SigmNorm>varepsilon1
        kk=kk+1;
    end
end
disp('Spatial Dimension: singular values')
disp(kk)

% Escale with singular values
hatx=Sigma(1:kk,1:kk)*U(:,1:kk)';
hatT=Sigma(1:kk,1:kk)*T(:,1:kk)';
Sigma=Sigma(1:kk,1:kk);
%[N,~]=size(hatT);

% HODMD: SPACE
if d1>1
    [GrowthRatex,Frequencyx,hatModesx] = DMDd_STKD(d1,hatx,Exis,varepsilon1,varepsilon2);
else
    [GrowthRatex,Frequencyx,hatModesx] = DMD1_STKD(hatx,Exis,varepsilon1,varepsilon2);
end
[NNx,MMx]=size(hatModesx);

Modesx0=hatModesx;
Modesx=zeros(NNx,MMx);
for pp=1:MMx
    Amplitudesx(pp)=norm(Modesx0(:,pp),2)/sqrt(NNx);
    Modesx(:,pp)=Modesx0(:,pp)/Amplitudesx(pp);
end
%GrowthRateOmegAmplx=[(1:MMx)',GrowthRatex',Frequencyx',Amplitudesx']

% HODMD: Time
if d2>1
    [GrowthRatet,Frequencyt,hatModest] =...
        DMDd_STKD(d2,hatT,Tiempos,varepsilon1,varepsilon2);
else
    [GrowthRatet,Frequencyt,hatModest] =DMD1_STKD(hatT,Tiempos,varepsilon1,varepsilon2);
end
[NNt,MMt]=size(hatModest);
Modest0=hatModest;
Modest=zeros(NNt,MMt);
for pp=1:MMt
    Amplitudest(pp)=norm(Modest0(:,pp),2)/sqrt(NNt);
    Modest(:,pp)=Modest0(:,pp)/Amplitudest(pp);
end

% DMD MODES: Q
q=(Modesx*diag(Amplitudesx))'*inv(Sigma)*(Modest*diag(Amplitudest));

Amplitudes=abs(q);
Modes=q./Amplitudes;
[MM,NN]=size(Amplitudes);
ModesNumber0=MM*NN;
ModesNumber=0;
for n=1:NN
    for m=1:MM
        if Amplitudes(m,n)/max(max(Amplitudes))<varepsilon2
            Amplitudes(m,n)=0;
        else
            ModesNumber= ModesNumber+1;
        end
    end
end
disp(ModesNumber0)
disp(ModesNumber)

qq=Amplitudes.*Modes;

% RECONSTRUCTION
Vreconst=Reconst_STKD(Tiempos,Exis,qq,GrowthRatex,Frequencyx,GrowthRatet,Frequencyt);
% MODES in X and T
disp('Number of spatial modes')
disp(MMx)
disp('Number of temporal modes')
disp(MMt)
