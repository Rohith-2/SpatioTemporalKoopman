function Vreconst=Reconst_STKD(Tiempos,Exis,q,deltasx,omegasx,deltast,omegast)

[N1,N2]=size(q);
J=length(Exis);
K=length(Tiempos);
%
vvx=zeros(N1,J);
for jj=1:J
    for m=1:N1
        vvx(m,jj)=exp((deltasx(m)+1i*omegasx(m))*Exis(jj));
    end
end
%
vvt=zeros(N2,K);
for kk=1:K
    for m=1:N2
        vvt(m,kk)=exp((deltast(m)+1i*omegast(m))*Tiempos(kk));
    end
end
%
Vreconst=vvx'*(q*vvt);
