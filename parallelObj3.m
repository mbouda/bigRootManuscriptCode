function J=parallelObj3(X,psiS,QR,QC,sc,nT,nL)

    Kr=abs(X)*sc(1);
    psiC=(sum(repmat(Kr,[1 nT]).*psiS)-QC')./sum(Kr);
    A=diag(Kr);
    
    qrMod=zeros(nL,nT);
    for i=1:nT
        b=Kr/2.*(psiC(i)+psiS(:,i));
        qrMod(:,i)=-(Kr*2).*(repmat(psiC(i),[nL 1])-A\b);
    end
    J=2e6*sum(sum((QR-qrMod).^2));    
    
end