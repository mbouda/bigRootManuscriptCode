function [QR,psiX]=predBigRootQRGC(X,psiS,GC,nL,nT,sc)

    kr=X(1:nL)/sc(1);
    Kx=X(nL+(1:nL))/sc(2);
    L=X(2*nL+(1:nL))/sc(3);
    
	[b2,c1,c2,c3,c4,c5,b12,b21,c12,c21,k12,k21,CnA,CpA,CnB,CpB]=formParams5(kr,Kx,L);
    
    A1=formPentAnAKG(b2,b12,b21,c12,c21,k12,k21,CnA,CpA,nL);
    
    psiX=zeros(nL,nT);
    for i=1:nT
        B1=formPentAnBKG(b2,c1,b12,b21,c12,c21,k12,k21,CnA,CpA,nL,psiS(:,i),GC(i));
        psiX(:,i)=A1\B1;
    end
    
    QR=(repmat(kr.*L,[1 nT]).*(psiS-psiX));    
    
end