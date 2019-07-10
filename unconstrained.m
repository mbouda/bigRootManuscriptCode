function [QRMod1,thMod1,QRMod2,thMod2]=unconstrained()

    dat=load('uswrcDroughtDataSet.mat');
    QR=dat.QR;
    TH=dat.TH;
	ch=dat.ch;
    psiS=dat.psiS;
    subset=dat.subset;
    t=dat.t;
    dZ=dat.dZ;
    rootFrac=dat.rootFrac;
    
    [nL,nT]=size(QR);
    dT=1800;
	
    QC=sum(QR);
    Kx1=(0.0015/1800)/0.5e6;
    GC=-QC/Kx1;
    
    myPool=parpool(24);
    warning off
    parfevalOnAll(myPool,@warning,0,'off');

    baseRes=load('baseLineRes.mat','XBR');
    sc=[1e8 1e11 1e2]*0.1;
    KR=(baseRes.XBR(1:nL)*sc(1)).*(baseRes.XBR((2*nL+1):3*nL)*sc(3));
    
    X0=ones(nL*8,1);
    
    
    
    opt=optimoptions('lsqnonlin','Display','iter','UseParallel',true,'MaxIter',2500,...
        'MaxFunEval',2e5,'Algorithm','Levenberg-Marquardt','TolFun',1e-10);
    
    [X1,F]=lsqnonlin(@(X)freeQRGCObj(X,psiS(:,subset),QR(:,subset),GC(subset),KR,nL,sum(subset)),X0,[],[],opt);

    [X2,F]=lsqnonlin(@(X)freeTHGCObj(X,TH,psiS,GC,KR,ch,nL,nT,dT,dZ),X0,[],[],opt);


    QRMod1=predFreeQRGC(X1,psiS,GC,KR,nL,nT);
    thMod1=predFreeTHGC(X1,psiS,ch,GC,KR,7,nT,dZ,dT);

    QRMod2=predFreeQRGC(X2,psiS,GC,KR,nL,nT);
    thMod2=predFreeTHGC(X2,psiS,ch,GC,KR,7,nT,dZ,dT);

end