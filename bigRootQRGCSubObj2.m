function J=bigRootQRGCSubObj2(X,psiS,dS,GC,nL,nT,dT,sc)

    X=abs(X);

    kr=X(1:nL)/sc(1);
    Kx=X(nL+(1:nL))/sc(2);
    L=X(2*nL+(1:nL))/sc(3);    
    
    QR=-dS/dT;

	[b2,c1,c2,c3,c4,c5,b12,b21,c12,c21,k12,k21,CnA,CpA,CnB,CpB]=formParams5(kr,Kx,L);
 
    A1=formPentAnAKG(b2,b12,b21,c12,c21,k12,k21,CnA,CpA,nL);
    A2=formPentAn2AKG(c2,c3,c5,CnB,CpB,k12,k21,nL);
    
    psiX1=zeros(nL,nT);
    psiX2=zeros(nL,nT);
    for i=1:nT
        B1=formPentAnBKG(b2,c1,b12,b21,c12,c21,k12,k21,CnA,CpA,nL,psiS(:,i),GC(i));
        B2=formPentAn2BKG(c1,c2,c3,c4,c5,k12,k21,CnB,CpB,nL,psiS(:,i),GC(i));
        psiX1(:,i)=A1\B1;
        psiX2(:,i)=A2\B2;
    end
    
    QR1=(repmat(kr.*L,[1 nT]).*(psiS-psiX1));
    QR2=(repmat(kr.*L,[1 nT]).*(psiS-psiX2));
    
    [G0,G1]=gradJG(b2,c1,c12,c21,CnA,CpA,psiX1,psiS,GC,nT);
    QR3=-repmat(Kx,[1 nT]).*(G1-G0);
        
    [G0,G1]=gradJG(b2,c1,c12,c21,CnA,CpA,psiX2,psiS,GC,nT);
    QR4=-repmat(Kx,[1 nT]).*(G1-G0);
    
    j1=(QR1-QR)';
    j2=(QR2-QR)'; 
    j3=(QR3-QR)'; 
    j4=(QR4-QR)'; 
            
    J=5e5*sum(sum((cat(1,j1,j2,j3,j4)).^2));
   
end