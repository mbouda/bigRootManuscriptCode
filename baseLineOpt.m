function [qrModPar,qrModSer,qrModBR]=baseLineOpt()


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
	
    KR=((1e-4)/1800)/5e5;
    Kr0=rootFrac*KR/sum(rootFrac);
    
    QC=sum(QR);
    Kx1=(0.0015/1800)/0.5e6;
    GC=-QC/Kx1;
    
    init=load('inits'); %load feasible starting points and scaling factors for each model
    X0BR=init.X0BR;
    scBR=init.scBR;
    X0Par=init.X0Par;
    scPar=init.scPar;
    X0Ser=init.X0Ser;
    scSer=init.scSer;

    opts=optimoptions('fmincon','Display','iter','MaxIter',250,...
        'MaxFunEval',25000,'UseParallel',true,'Algorithm','interior-point');
    
    myPool=parpool(24);
    warning off
    parfevalOnAll(myPool,@warning,0,'off');


	XBR=fmincon(@(X)bigRootQRGCSubObj2(X,psiS(:,subset),dS(:,subset),GC(subset),nL,sum(subset),dT,scBR),X0BR,...
            [],[],[],[],[],[],@(X)bigRootQRGCSubSumCon(X,scBR,Kr0,nL),opts);
    XBR=abs(XBR);

    XPar=fmincon(@(X)parallelObj3(X,psiS(:,subset),QR(:,subset),QC(subset),scPar,sum(subset),nL),X0Par,...
        [],[],[],[],[],[],@(X)parallelSumCon(X,KR,scPar),opts);

    XSer=fmincon(@(X)triDiagObj(X,psiS(:,subset),QR(:,subset),sum(subset),nL,scSer),X0Ser,...
        [],[],[],[],[],[],@(X)seriesSumCon(X,KR,scSer,nL),opts);
    XSer=abs(XSer);
                
    qrModBR=predBigRootQRGC(XBR,psiS,GC,nL,nT,scBR);
    
    KrPar=abs(XPar)*scPar;
    psiC=(sum(repmat(KrPar,[1 nT]).*psiS)-QC')./sum(KrPar);
    A=diag(KrPar);
    qrModPar=zeros(nL,nT);
    for i=1:nT
        b=KrPar/2.*(psiC(i)+psiS(:,i));
        qrModPar(:,i)=-(KrPar*2).*(repmat(psiC(i),[nL 1])-A\b);
    end
    
    KrSer=XSer(1:nL)*scSer(1);
    KxSer=XSer(nL+(1:nL-1))*scSer(2);
    ci=KrSer;
    ci(1)=ci(1)+KxSer(1);
    ci(2:end-1)=ci(2:end-1)+KxSer(1:end-1)+KxSer(2:end);
    ci(end)=ci(end)+KxSer(end);
    A=diag(ci)-diag(KxSer,1)-diag(KxSer,-1);
    psiX=zeros(nL,nT);
    for i=1:nT
        b=KrSer.*psiS(:,i);
        b(1)=b(1)-QC(i)*scSer(3);
        psiX(:,i)=A\b;
    end
    qrModSer=-repmat(KrSer,[1 nT]).*(psiX-psiS)/scSer(3);


end