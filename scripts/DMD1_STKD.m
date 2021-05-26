
function  [GrowthRates,Frequency,hatmodos] =DMD1_STKD(hatT,Tiempos,varepsilon1,varepsilon2)

[Jprima,K]=size(hatT);
N=Jprima;
Deltat=Tiempos(2)-Tiempos(1);
%
[hatU1,hatSigma,hatU2]=svd(hatT(:,1:K-1),'econ');
%% Calculate Koopman operator
hatR=hatT(:,2:K)*hatU2*inv(hatSigma)*hatU1';
[Q,tildeMM]=eig(hatR);
eigenvalues=diag(tildeMM);

M=length(eigenvalues);
qq=log(eigenvalues);
GrowthRates=real(qq)/Deltat;
Frequency=imag(qq)/Deltat;
%% Normalize eigenvectors
[NN,MMM]=size(Q);
for m=1:MMM
    Qnorm=Q(:,m);
    Q(:,m)= Q(:,m)/norm(Qnorm(:),2);
end

%% Calculate amplitudes
Mm=zeros(NN*K,M);
Bb=zeros(NN*K,1);
Id=eye(MMM);
for k=1:K
    Mm(1+(k-1)*NN:k*NN,:)=Q*Id;
    Id=Id*tildeMM;
    Bb(1+(k-1)*NN:k*NN,1)=hatT(:,k);
end

[Ur,Sigmar,Vr]=svd(Mm,'econ');
a=Vr*(Sigmar\(Ur'*Bb));

u=zeros(NN,M);
for m=1:M
    u(:,m)=a(m)*Q(:,m);
end

hatamplitudes=zeros(M,1);

for m=1:M
    aca=u(:,m);
    hatamplitudes(m)=norm(aca(:),2);
end

UU=[u;GrowthRates';Frequency';hatamplitudes']';
UU1=sortrows(UU,-(NN+3));

UU=UU1';
u=UU(1:NN,:);
GrowthRates=UU(NN+1,:);
Frequency=UU(NN+2,:);
hatamplitudes=UU(NN+3,:);

kk3=0;
for m=1:M
    if hatamplitudes(m)/hatamplitudes(1)>varepsilon2
        kk3=kk3+1;
    else
    end
end
%% Spectral complexity: number of DMD modes.
('Spectral complexity')
kk3

%amplitudes=UU(NN+3,:);
hatmodos=u(:,1:kk3);
Frequency=Frequency(1:kk3);
GrowthRates=GrowthRates(1:kk3);
