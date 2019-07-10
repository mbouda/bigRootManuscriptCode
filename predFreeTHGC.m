function thMod=predFreeTHGC(X,psiS,ch,GC,KR,nL,nT,dZ,dT)

    coeffs=reshape(X,[nL 8]);
    A=formPentAGC(coeffs);
    
    opts=optimoptions('fsolve','Display','off');
    psiMod=zeros(nL,nT);
    psiMod(:,1)=psiS(:,1);
    for i=1:nT-1
        psiMod(:,i+1)=fsolve(@(x)crankNicholsonFreeTH(x,psiMod(:,i),ch,A,coeffs,KR,GC,dT,dZ),psiMod(:,i),opts);
    end
    
    thMod=clappHornTh(psiMod,ch);
    