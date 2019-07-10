function [j1,j2,j3,j4]=fit1BGEval(Kr,Kx,L,psiS,psiX,QR,GC,nSol,nL)

    [b2,c1,c2,c3,c4,c5,b12,b21,c12,c21,k12,k21,CnA,CpA,CnB,CpB]=formParams5(Kr,Kx,L);
    
    A1=formPentAnAKG(b2,b12,b21,c12,c21,k12,k21,CnA,CpA,nL);
    A2=formPentAn2AKG(c2,c3,c5,CnB,CpB,k12,k21,nL);
    
    pred1=zeros(nL,nSol);
    pred2=zeros(nL,nSol);
    for i=1:nSol
        B1=formPentAnBKG(b2,c1,b12,b21,c12,c21,k12,k21,CnA,CpA,nL,psiS(:,i),GC(i));
        B2=formPentAn2BKG(c1,c2,c3,c4,c5,k12,k21,CnB,CpB,nL,psiS(:,i),GC(i));
        pred1(:,i)=A1\B1;
        pred2(:,i)=A2\B2;
    end
    j1=(pred1-psiX)';
    j2=(pred2-psiX)';
    
    QR1=(repmat(Kr.*L,[1 nSol]).*(psiS-pred1));
    QR2=(repmat(Kr.*L,[1 nSol]).*(psiS-pred2));
       
    j3=((QR1-QR))';
    j4=((QR2-QR))'; 


end