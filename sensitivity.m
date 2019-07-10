function sensitivity()

    dat=load('uswrcDroughtDataSet.mat');
    psiS=dat.psiS;
    QR=dat.QR;
    TH=dat.TH; %these data missing from release due to rights
               %can be obtained from US-Wrc PI
    
    rootFrac=dat.rootFrac;
    dZ=dat.dZ;
    ch=dat.ch;
   
    [nL,nT]=size(QR);
    dT=1800;
	
    KR=((1e-4)/1800)/5e5;
    Kr0=rootFrac*KR/sum(rootFrac);
    
    QC=sum(QR);
    Kx1=(0.0015/1800)/0.5e6;
    GC=-QC/Kx1;
    
    prevOpt=load('baseLineRes','XBR','XPar','XSer');
    
    XBR=abs(prevOpt.XBR);
    cBR=0.1;
    scBR=[1e8 1e11 1e2]'*cBR; 
    
    cPar=1e6;
    scPar=1/cPar;
    XPar=abs(prevOpt.XPar);
    KrPar=abs(XPar)*scPar;
    
    cr=1e-14;
    cx=1e-10;
    cq=1e0;
    scSer=[cr cx cq];
    XSer=abs(prevOpt.XSer);
    KrSer=XSer(1:nL)*scSer(1);
    KxSer=XSer(nL+(1:nL-1))*scSer(2);

       
    myPool=parpool(24);
    warning off
    parfevalOnAll(myPool,@warning,0,'off');
    
    eBR=zeros(nL,nTrials);
    ePar=zeros(nL,nTrials);
    eSer=zeros(nL,nTrials);
    
    pertB=[-0.2 0 0.2];
    pertP=[-500 0 500];
    
    perturb=combvec(pertB,pertP,...
        pertB,pertP,...
        pertB,pertP,...
        pertB,pertP,...
        pertB,pertP,...
        pertB,pertP,...
        pertB,pertP);
    for i=1:7
        k=(2*(i-1))+(1:2);
        del=any(perturb(k,:)==0) & ~all(perturb(k,:)==0);
        perturb(:,del)=[];
    end
    
    nTrials=size(perturb,2);
    
    parfor i=1:nTrials
       
        ch2=ch;
        d=perturb(:,i);
        d=reshape(d,[2,7]);
        dB=d(1,:)';
        dP=d(2,:)';
        ch2.B=ch2.B+dB;
        ch2.pSat=ch2.pSat+dP;
        
        thModBR=predBigRootTHGC(XBR,TH,ch2,GC,nL,nT,dZ,dT,scBR);
        eBR(:,i)=sqrt(mean((thModBR-TH).^2,2));

        thModPar=parallelPred(KrPar,TH',ch2,QC,nL,nT,dT,dZ); 
        ePar(:,i)=sqrt(mean((thModPar-TH).^2,2));
        
        thModSer=predTriDiag(KrSer,KxSer,TH',ch2,QC,nL,nT,dT,dZ,scSer);
        eSer(:,i)=sqrt(mean((thModSer-TH).^2,2));
    end
    
    save('errorSens','eBR','ePar','eSer')
    
    exit 
end
