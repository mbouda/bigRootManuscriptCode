function J=forwPentAnQRObj(X,psiS,psiC,QR,nL,nSol,sc)

    kr=X(1:nL)/sc(1);
    Kx=X(nL+(1:nL))/sc(2);
    L=X(2*nL+(1:nL))/sc(3);
    
	[b2,c1,c2,c3,c4,c5,b12,b21,c12,c21,k12,k21,CnA,CpA,CnB,CpB]=formParams5(kr,Kx,L);
    
    A1=formPentAnAK2(b2,c4,b12,b21,c12,c21,k12,k21,CnA,CpA,nL);
    A2=formPentAn2AKOld(c2,c3,c5,CnB,CpB,k12,k21,nL);
    
    psiX1=zeros(nL,nSol);
    psiX2=zeros(nL,nSol);
    for i=1:nSol
        B1=formPentAnBK2(b2,c1,c3,c4,b12,b21,c12,c21,k12,k21,CnA,CpA,nL,psiS(:,i),psiC(i));
        B2=formPentAn2BKOld(c1,c2,c3,c5,k12,k21,CnB,CpB,nL,psiS(:,i),psiC(i));
        psiX1(:,i)=A1\B1;
        psiX2(:,i)=A2\B2;
    end
    
    QR1=(repmat(kr.*L,[1 nSol]).*(psiS-psiX1));
    QR2=(repmat(kr.*L,[1 nSol]).*(psiS-psiX2));
    
    [G0,G1]=gradJOld(b2,c1,c2,c3,c4,c12,c21,CnA,CpA,psiX1,psiS,psiC,nSol);
    QR3=-repmat(Kx,[1 nSol]).*(G1-G0);
        
    [G0,G1]=gradJOld(b2,c1,c2,c3,c4,c12,c21,CnA,CpA,psiX2,psiS,psiC,nSol);
    QR4=-repmat(Kx,[1 nSol]).*(G1-G0);
    
    j1=(QR1-QR)';
    j2=(QR2-QR)'; 
    j3=(QR3-QR)'; 
    j4=(QR4-QR)'; 
    
    J=cat(1,j1,j2,j3,j4);

end