function thMod=parallelPred(Kr,th,ch,QC,nL,nT,dT,dZ,sc)
      
    opts=optimoptions('fsolve','Display','off','FiniteDifferenceType','forward','UseParallel',false);
    psiMod=zeros(nL,nT);
    psiMod(:,1)=clappHornPsi(th(1,:)',ch);
    for i=1:nT-1
        psiMod(:,i+1)=fsolve(@(x)crankNicholsonParallelTH2(x,psiMod(:,i),ch,Kr,QC(i),sc,dT,dZ,nL),psiMod(:,i),opts);

    end
    
    thMod=clappHornTh(psiMod,ch);
    
end