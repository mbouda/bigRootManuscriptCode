function fitCase1BG(dataDir)

    addpath ~/oldCode/
    inputDir=strcat('~/work/synth/',dataDir,filesep);

    inp=load(strcat(inputDir,'case1bResSM'));
    GC=inp.GC;
    psiS=inp.psiBC(2:end,:);
    psiX=inp.psiX;
    QR=-inp.QR;

    clear inp
    
    [nLayers,nSol]=size(psiX);

    sc=[1e8 1e11 1e2 1e9]';
    
    X0=ones(3*nLayers,1);

    myPool=parpool(24);
    
    warning off
    parfevalOnAll(myPool,@warning,0,'off');
    
    opts=optimoptions('lsqnonlin','Display','iter','Algorithm','Levenberg-Marquardt',...
        'UseParallel',true,'FiniteDifferenceType','forward',...
        'MaxIter',150,'StepTolerance',1e-10,'FunctionTolerance',1e-24);
    [X,F,~,~,out]=lsqnonlin(@(X)forwPentAnCombObjConG(X,psiS,GC,psiX,QR,nLayers,nSol,sc),X0,[],[],opts);

    save(strcat(inputDir,'optOutsConG'),'X','F','out');
    
    delete(gcp('nocreate'))

end

