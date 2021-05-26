
function  [GrowthRate,Frequency,hatmodos] =DMDd_STKD(d,hatT,Tiempos,varepsilon1,varepsilon2)
[Jprima,K]=size(hatT);
N=Jprima;

% Create the modified snapshot matrix
tildeT=zeros(d*N,K-d+1);
for ppp=1:d
 tildeT((ppp-1)*N+1:ppp*N,:)=hatT(:,ppp:ppp+K-d);   
end

% Dimension reduction
[U1,Sigma1,T1]=svd(tildeT);
sigmas1=diag(Sigma1);

Deltat=Tiempos(2)-Tiempos(1);
n=length(sigmas1);

NormS=norm(sigmas1,2);
kk1=0;


for k=1:n
    RRMSEE(k)=norm(sigmas1(k:n),2)/NormS;
        if RRMSEE(k)>varepsilon1
        kk1=kk1+1;
        end
end

disp('Spatial dimension reduction')
disp(kk1)

U1=U1(:,1:kk1);
hatT1=Sigma1(1:kk1,1:kk1)*T1(:,1:kk1)';

% Reduced modified snapshot matrix
[~,K1]=size(hatT1);
[tildeU1,tildeSigma,tildeU2]=svd(hatT1(:,1:K1-1),'econ');

% Reduced modified Koopman matrix
tildeR=hatT1(:,2:K1)*tildeU2/tildeSigma*tildeU1';
[tildeQ,tildeMM]=eig(tildeR);
eigenvalues=diag(tildeMM);

M=length(eigenvalues);
qq=log(eigenvalues);
GrowthRate=real(qq)/Deltat;
Frequency=imag(qq)/Deltat;

Q=U1*tildeQ;
Q=Q((d-1)*N+1:d*N,:);
[NN,MMM]=size(Q);
 
 for m=1:MMM
    NormQ=Q(:,m);
    Q(:,m)= Q(:,m)/norm(NormQ(:),2);
 end

Mm=zeros(NN*K,M);
Bb=zeros(NN*K,1);
aa=eye(MMM);
for k=1:K
 Mm(1+(k-1)*NN:k*NN,:)=Q*aa; 
 aa=aa*tildeMM;
 Bb(1+(k-1)*NN:k*NN,1)=hatT(:,k);
end

[Ur,Sigmar,Vr]=svd(Mm,'econ');
a=Vr*(Sigmar\(Ur'*Bb));
disp("size of a")
size(a)
disp("size of Q")
size(Q)
u=zeros(NN,M);
for m=1:M
    u(:,m)=a(m)*Q(:,m);
end
hatamplitudes=zeros(M,1);

for m=1:M
    aca=u(:,m);
    hatamplitudes(m)=norm(aca(:),2);    
end

UU=[u;GrowthRate';Frequency';hatamplitudes']';
UU1=sortrows(UU,-(NN+3));

UU=UU1';
u=UU(1:NN,:);
GrowthRate=UU(NN+1,:);
Frequency=UU(NN+2,:);
hatamplitudes=UU(NN+3,:);

kk3=0;
for m=1:M
    if hatamplitudes(m)/hatamplitudes(1)>varepsilon2
        kk3=kk3+1;
    end
end

% Spectral complexity: number of DMD modes.
disp('Spectral complexity')
disp(kk3)

%amplitudes=UU(NN+3,:);
hatmodos=u(:,1:kk3);
Frequency=Frequency(1:kk3);
GrowthRate=GrowthRate(1:kk3);
