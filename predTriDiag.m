function thMod=predTriDiag(Kr,Kx,th,ch,QC,nL,nT,dT,dZ,sc)

    ci=Kr;
    ci(1)=ci(1)+Kx(1);
    ci(2:end-1)=ci(2:end-1)+Kx(1:end-1)+Kx(2:end);
    ci(end)=ci(end)+Kx(end);
    A=diag(ci)-diag(Kx,1)-diag(Kx,-1);

    opts=optimoptions('fsolve','Display','off','FiniteDifferenceType','forward','UseParallel',false);
  
    psiMod=zeros(nL,nT);
    psiMod(:,1)=clappHornPsi(th(1,:)',ch);
    for i=1:nT-1
          psiMod(:,i+1)=fsolve(@(x)crankNicholsonSeriesTH2(x,psiMod(:,i),ch,A,Kr,QC(i),sc,dT,dZ),psiMod(:,i),opts);
    end
    thMod=clappHornTh(psiMod,ch);
    
end