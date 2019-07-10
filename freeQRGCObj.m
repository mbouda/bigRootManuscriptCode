function J=freeQRGCObj(X,psiS,QR,GC,KR,nL,nT)
   
    coeffs=reshape(X,[nL 8]);
    
    A=formPentAGC(coeffs);
    psiX=zeros(nL,nT);
    for i=1:nT
        B=formPentBGC(coeffs,psiS(:,i),GC(i));
        psiX(:,i)=A\B;
    end
    QRM=(repmat(KR,[1 nT]).*(psiS-psiX));
    
    J=1e6*(QRM-QR)';
        
end