function fitCase1BQR(dataDir)

    addpath ~/oldCode/
    inputDir=strcat('~/work/synth/',dataDir,filesep);

    inp=load(strcat(inputDir,'case1bResSM'));
    psiC=inp.psiBC(1,:);
    psiS=inp.psiBC(2:end,:);
    QR=-inp.QR;
    clear inp
    
    [nLayers,nSol]=size(QR);
    sc=[1e8 1e11 1e2]';
    X0=ones(3*nLayers,1);

    myPool=parpool(24);
    
    warning off
    parfevalOnAll(myPool,@warning,0,'off');
    
    opts=optimoptions('lsqnonlin','Display','iter','Algorithm','Levenberg-Marquardt',...
        'UseParallel',true,'FiniteDifferenceType','forward',...
        'MaxIter',150,'StepTolerance',1e-10,'FunctionTolerance',1e-24);
    [X,F,~,~,out]=lsqnonlin(@(X)forwPentAnQRObj(X,psiS,psiC,QR,nLayers,nSol,sc),X0,[],[],opts);

    save(strcat(inputDir,'optOutsQR'),'X','F','out');
    
    delete(gcp('nocreate'))

end


