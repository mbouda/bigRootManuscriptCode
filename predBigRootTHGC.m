function thMod=predBigRootTHGC(X,th,ch,GC,nL,nT,dZ,dT,sc)

    kr=X(1:nL)/sc(1);
    Kx=X(nL+(1:nL))/sc(2);
    L=X(2*nL+(1:nL))/sc(3);
    
	[b2,c1,~,~,~,~,b12,b21,c12,c21,k12,k21,CnA,CpA,~,~]=formParams5(kr,Kx,L);
    
    A1=formPentAnAKG(b2,b12,b21,c12,c21,k12,k21,CnA,CpA,nL);

    opts=optimoptions('fsolve','Display','off','FiniteDifferenceType','forward','UseParallel',false);
    psiMod=zeros(nL,nT);
    psiMod(:,1)=clappHornPsi(th(:,1)',ch);

    for i=1:nT-1
        psiMod(:,i+1)=fsolve(@(x)crankNicholsonQRTH(x,psiMod(:,i),ch,A1,L,kr,b2,c1,b12,b21,c12,c21,k12,k21,CnA,CpA,nL,GC(i:i+1),dT,dZ),psiMod(:,i),opts);
    end
    
    thMod=clappHornTh(psiMod,ch);
    
    
end