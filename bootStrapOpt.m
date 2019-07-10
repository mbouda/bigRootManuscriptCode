function [eBR,ePar,eSer,thModBR,thModPar,thModSer]=bootStrapOpt()

    dat=load('uswrcDroughtDataSet.mat');
    psiS=dat.psiS;
    QR=dat.QR;
    TH=dat.TH; %these data missing from release due to rights
               %can be obtained from US-Wrc PI
    
    rootFrac=dat.rootFrac;
    dZ=dat.dZ;
    ch=dat.ch;
    
    t=dat.t;
    dT=1800;
    [nL,nT]=size(QR);
    
    KR=((1e-4)/1800)/5e5;
    Kr0=rootFrac*KR/sum(rootFrac);

    Kx1=(0.0015/1800)/0.5e6;
    QC=sum(QR)';
    GC=QR/Kx1;
    
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
    
    
    nTrials=1000;
    nDay=3;
    
    midNight=find(t==round(t)); midNight(end)=[];
    totDays=length(midNight);
    
    eBR=zeros(nL,nTrials);
    ePar=zeros(nL,nTrials);
    eSer=zeros(nL,nTrials);
    
    thModBR=zeros(nL,nT,nTrials);
    thModPar=zeros(nL,nT,nTrials);
    thModSer=zeros(nL,nT,nTrials);
    
    for i=1:nTrials
        
        idx = datasample(1:totDays,nDay,'Replace',false);
        subset=false(nT,1);
        for j=1:nDay
            subset(midNight(idx(j)):(midNight(idx(j))+47))=true;
        end

        XBR=fmincon(@(X)bigRootQRGCSubObj2(X,psiS(:,subset),dS(:,subset),GC(subset),nL,sum(subset),dT,scBR),X0BR,...
            [],[],[],[],[],[],@(X)bigRootQRGCSubSumCon(X,scBR,Kr0,nL),opts);
        XBR=abs(XBR);

        XPar=fmincon(@(X)parallelObj3(X,psiS(:,subset),QR(:,subset),QC(subset),scPar,sum(subset),nL),X0Par,...
            [],[],[],[],[],[],@(X)parallelSumCon(X,KR,scPar),opts);

        XSer=fmincon(@(X)triDiagObj(X,psiS(:,subset),QR(:,subset),sum(subset),nL,scSer),X0Ser,...
            [],[],[],[],[],[],@(X)seriesSumCon(X,KR,scSer,nL),opts);
        XSer=abs(XSer);
                
        thModBR(:,:,i)=predBigRootTHGC(XBR,TH,ch,GC,nL,nT,dZ,dT,scBR);
        eBR(:,i)=sqrt(mean((thModBR(:,:,i)-TH).^2,2));

        KrPar=abs(XPar)*scPar(1);
        thModPar(:,:,i)=parallelPred(KrPar,TH',ch,QC,nL,nT,dT,dZ,scPar); 
        ePar(:,i)=sqrt(mean((thModPar(:,:,i)-TH).^2,2));
        
        KrSer=XSer(1:nL)*scSer(1);
        KxSer=XSer(nL+(1:nL-1))*scSer(2);
        thModSer(:,:,i)=predTriDiag(KrSer,KxSer,TH',ch,QC,nL,nT,dT,dZ,scSer); 
        eSer(:,i)=sqrt(mean((thModSer(:,:,i)-TH).^2,2));

%       Reporting:
%         if any([0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]==i/nTrials)
%             fprintf(1,'------------\n',i);
%             fprintf(1,'\n',i);
%             fprintf(1,'Finished %d \n',i);
%             fprintf(1,'\n',i);
%             fprintf(1,'------------\n',i);
%         end
        
    end
end
