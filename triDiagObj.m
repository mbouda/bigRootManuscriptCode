function J=triDiagObj(X,psiS,QR,nT,nL,sc)

    X=abs(X);
    Kr=X(1:nL)*sc(1);
    Kx=X(nL+(1:nL-1))*sc(2);
    QR=QR*sc(3);
    
    QC=sum(QR)';
   
    psiX=QR'./repmat(-Kr',[nT 1])+psiS';
    QX=repmat(Kx',[nT 1]).*(psiX(:,2:end)-psiX(:,1:end-1));
    
    QX=[QC QX zeros(nT,1)];
    
    r=(QX(:,1:end-1)-(QR'+QX(:,2:end)));
    
    J=sum(sum(r.^2));
end 