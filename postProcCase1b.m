%This script post-processes case1b fits for four boundary conditions


%want a table of errors in 
%psiX
%QR
%kr * S
%Kx / S

baseDir='~/work/synth/';
files=dir(baseDir);
files=struct2cell(files)';
caseDirs=files(contains(files(:,1),'case1BCon'),1);
nDirs=size(caseDirs,1);

sc=[1e8 1e11 1e2 1e9]';

JPX=cell(nDirs,6);
JPQ=cell(nDirs,6);
JGX=cell(nDirs,6);
JGQ=cell(nDirs,6);

RPX=cell(nDirs,6);
RPQ=cell(nDirs,6);
RGX=cell(nDirs,6);
RGQ=cell(nDirs,6);

myPool=parpool(24);
warning off
parfevalOnAll(myPool,@warning,0,'off');

parfor i=1:nDirs
    
    inp=load(strcat(baseDir,caseDirs{i},filesep,'case1bResSM'));
    
    
    psiC=inp.psiBC(1,:);
    GC=inp.GC;
    psiS=inp.psiBC(2:end,:);
    psiX=inp.psiX;
    QR=-inp.QR;
 	rc=inp.rc;
    [nL,nSol]=size(psiS);
    
    [L,Kr,Kx,~,~,~,~,~,~]=constructCase1();
	L=L.*rc(:,1); 
	Kr=Kr.*rc(:,2);
    Kx=Kx.*rc(:,3);
    
    KrS=Kr.*L;
    KxS=Kx./L;
    
    res=load(strcat(baseDir,caseDirs{i},filesep,'optOutsCon2'));
    
    KrX=res.X(1:nL)/sc(1);
    KxX=res.X(nL+(1:nL))/sc(2);
    LX=res.X(2*nL+(1:nL))/sc(3);
    
    modKrS=KrX.*LX;
    modKxS=KxX./LX;
    
%     [JPX{i,1},JPX{i,2},JPX{i,3},JPX{i,4}]=fit1BEval(KrX,KxX,LX,psiS,psiX,QR,psiC,nSol,nL);
%     JPX{i,5}=modKrS-KrS;
%     JPX{i,6}=modKxS-KxS;
    [j1,j2,j3,j4]=fit1BEval(KrX,KxX,LX,psiS,psiX,QR,psiC,nSol,nL);
    j5=modKrS-KrS;
    j6=modKxS-KxS;
    JPX(i,:)={j1,j2,j3,j4,j5,j6};
    RPX(i,:)={psiX' psiX' QR' QR' KrS KxS};
    
    res=load(strcat(baseDir,caseDirs{i},filesep,'optOutsConG'));
    
    KrX=res.X(1:nL)/sc(1);
    KxX=res.X(nL+(1:nL))/sc(2);
    LX=res.X(2*nL+(1:nL))/sc(3);
    
    modKrS=KrX.*LX;
    modKxS=KxX./LX;
    
%     [JGX{i,1},JGX{i,2},JGX{i,3},JGX{i,4}]=fit1BGEval(KrX,KxX,LX,psiS,psiX,QR,GC,nSol,nL);
%     JGX{i,5}=modKrS-KrS;
%     JGX{i,6}=modKxS-KxS;
    [j1,j2,j3,j4]=fit1BGEval(KrX,KxX,LX,psiS,psiX,QR,GC,nSol,nL);
    j5=modKrS-KrS;
    j6=modKxS-KxS;
    JGX(i,:)={j1,j2,j3,j4,j5,j6};
    RGX(i,:)={psiX' psiX' QR' QR' KrS KxS};
    
    res=load(strcat(baseDir,caseDirs{i},filesep,'optOutsQR'));
    
    KrX=res.X(1:nL)/sc(1);
    KxX=res.X(nL+(1:nL))/sc(2);
    LX=res.X(2*nL+(1:nL))/sc(3);
    
    modKrS=KrX.*LX;
    modKxS=KxX./LX;
    
%     [JPQ{i,1},JPQ{i,2},JPQ{i,3},JPQ{i,4}]=fit1BEval(KrX,KxX,LX,psiS,psiX,QR,psiC,nSol,nL);
%     JPQ{i,5}=modKrS-KrS;
%     JPQ{i,6}=modKxS-KxS;
    [j1,j2,j3,j4]=fit1BEval(KrX,KxX,LX,psiS,psiX,QR,psiC,nSol,nL);
    j5=modKrS-KrS;
    j6=modKxS-KxS;
    JPQ(i,:)={j1,j2,j3,j4,j5,j6};
    RPQ(i,:)={psiX' psiX' QR' QR' KrS KxS};

    
    res=load(strcat(baseDir,caseDirs{i},filesep,'optOutsQRUnc'));
    
    KrX=res.X(1:nL)/sc(1);
    KxX=res.X(nL+(1:nL))/sc(2);
    LX=res.X(2*nL+(1:nL))/sc(3);
    
    modKrS=KrX.*LX;
    modKxS=KxX./LX;
    
%     [JGQ{i,1},JGQ{i,2},JGQ{i,3},JGQ{i,4}]=fit1BGEval(KrX,KxX,LX,psiS,psiX,QR,GC,nSol,nL);
%     JGQ{i,5}=modKrS-KrS;
%     JGQ{i,6}=modKxS-KxS;
    
    [j1,j2,j3,j4]=fit1BGEval(KrX,KxX,LX,psiS,psiX,QR,GC,nSol,nL);
    j5=modKrS-KrS;
    j6=modKxS-KxS;
    JGQ(i,:)={j1,j2,j3,j4,j5,j6};
    RGQ(i,:)={psiX' psiX' QR' QR' KrS KxS};
end

save('case1BFits','JPX','JGX','JPQ','JGQ','RPX','RGX','RPQ','RGQ');

%%

load case1BFits

epX=zeros(20,2);
rpX=zeros(20,2);
for i=1:20
    for j=1:2
        epX(i,j)=sqrt(mean(mean(JPX{i,j}.^2)));
        rpX(i,j)=sqrt(mean(mean(RPX{i,j}.^2)));
    end
end
%12 and 17 not converged
%figure; histogram(log(epX)) %shows that all others have

eQX=zeros(20,2);
rQX=zeros(20,2);
for i=1:20
    for j=3:4
        eQX(i,j-2)=sqrt(mean(mean(JPX{i,j}.^2)));
        rQX(i,j-2)=sqrt(mean(mean(RPX{i,j}.^2)));
    end
end

erX=zeros(20,1);
rrX=zeros(20,1);
for i=1:20
    for j=5
        erX(i,j-4)=sqrt(mean(JPX{i,j}.^2));
        rrX(i,j-4)=sqrt(mean(RPX{i,j}.^2));
    end
end

exX=zeros(20,1);
rxX=zeros(20,1);
for i=1:20
    for j=6
        exX(i,j-5)=sqrt(mean(JPX{i,j}.^2));
        rxX(i,j-5)=sqrt(mean(RPX{i,j}.^2));
    end
end

epX([12 17],:)=[];
eQX([12 17],:)=[];
erX([12 17],:)=[];
exX([12 17],:)=[];

rpX([12 17],:)=[];
rQX([12 17],:)=[];
rrX([12 17],:)=[];
rxX([12 17],:)=[];

table(:,1)=[ mean(mean(epX)); mean(mean(eQX)); mean(mean(erX)); mean(mean(exX))];
relTable(:,1)=[ mean(mean(rpX)); mean(mean(rQX)); mean(mean(rrX)); mean(mean(rxX))];

epX=zeros(20,2);
for i=1:20
    for j=1:2
        epX(i,j)=sqrt(mean(mean(JGX{i,j}.^2)));
    end
end
%8 not converged
%figure; histogram(log(epX)) %shows that all others have

eQX=zeros(20,2);
for i=1:20
    for j=3:4
        eQX(i,j-2)=sqrt(mean(mean(JGX{i,j}.^2)));
    end
end

erX=zeros(20,1);
for i=1:20
    for j=5
        erX(i,j-4)=sqrt(mean(JGX{i,j}.^2));
    end
end

exX=zeros(20,1);
for i=1:20
    for j=6
        exX(i,j-5)=sqrt(mean(JGX{i,j}.^2));
    end
end

epX(8,:)=[];
eQX(8,:)=[];
erX(8,:)=[];
exX(8,:)=[];
table(:,2)=[ mean(mean(epX)); mean(mean(eQX)); mean(mean(erX)); mean(mean(exX))];

%averaging across layers is a problem because n'th layer is wrong.
epX=zeros(20,2);
for i=1:20
    for j=1:2
        epX(i,j)=sqrt(mean(mean(JPQ{i,j}(:,1:end-1).^2)));
        epN(i,j)=sqrt(mean(JPQ{i,j}(:,end).^2));
    end
end

eQX=zeros(20,2);
for i=1:20
    for j=3:4
        eQX(i,j-2)=sqrt(mean(mean(JPQ{i,j}(:,1:end-1).^2)));
        eQN(i,j-2)=sqrt(mean(JPQ{i,j}(:,end).^2));
    end
end

erX=zeros(20,1);
for i=1:20
    for j=5
        erX(i,j-4)=sqrt(mean(JPQ{i,j}(1:end-1).^2));
        erN(i,j-4)=sqrt(mean(JPQ{i,j}(end).^2));
    end
end

exX=zeros(20,1);
for i=1:20
    for j=6
        exX(i,j-5)=sqrt(mean(JPQ{i,j}(1:end-1).^2));
        exN(i,j-5)=sqrt(mean(JPQ{i,j}(end).^2));
    end
end

table(:,3)=[ mean(mean(epX)); mean(mean(eQX)); mean(mean(erX)); mean(mean(exX))];
table(:,4)=[ mean(mean(epN)); mean(mean(eQN)); mean(mean(erN)); mean(mean(exN))];



epX=zeros(20,2);
for i=1:20
    for j=1:2
        epX(i,j)=sqrt(mean(mean(JGQ{i,j}.^2)));
    end
end

eQX=zeros(20,2);
for i=1:20
    for j=3:4
        eQX(i,j-2)=sqrt(mean(mean(JGQ{i,j}.^2)));
    end
end

erX=zeros(20,1);
for i=1:20
    for j=5
        erX(i,j-4)=sqrt(mean(JGQ{i,j}.^2));
    end
end

exX=zeros(20,1);
for i=1:20
    for j=6
        exX(i,j-5)=sqrt(mean(JGQ{i,j}.^2));
    end
end

table(:,5)=[ mean(mean(epX)); mean(mean(eQX)); mean(mean(erX)); mean(mean(exX))];

save('errorTable','table')


psiXAll=cell2mat(cat(2,RPX(:,1:2),RPQ(:,1:2),RGX(:,1:2),RGQ(:,1:2)));
psiXMean=mean(mean(psiXAll));
psiXSTD=std(reshape(psiXAll,[numel(psiXAll) 1]));

QRAll=cell2mat(cat(2,RPX(:,3:4),RPQ(:,3:4),RGX(:,3:4),RGQ(:,3:4)));
QRMean=mean(mean(abs(QRAll)));
QRSTD=std(reshape(QRAll,[numel(QRAll) 1]));

KrAll=cell2mat(cat(2,RPX(:,5),RPQ(:,5),RGX(:,5),RGQ(:,5)));
KrMean=mean(mean(KrAll));
KrSTD=std(reshape(KrAll,[numel(KrAll) 1]));

KxAll=cell2mat(cat(2,RPX(:,6),RPQ(:,6),RGX(:,6),RGQ(:,6)));
KxMean=mean(mean(KxAll));
KxSTD=std(reshape(KxAll,[numel(KxAll) 1]));

save('moments1b','psiXMean','psiXSTD','QRMean','QRSTD','KrMean','KrSTD','KxMean','KxSTD')


