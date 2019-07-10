function  singleRootLibrary()

%This function builds a library of analytical solutions for different sets of
%moisture conditions, using a simple synthetic system: one root passing
%through seven layers, with uniform but different properties in each layer

jobID=getenv('PBS_JOBID');
C=regexp(jobID,'\d+','match');
jobNr=C{1};
rng(str2double(jobNr));

[L,Kr,Kx,a,~,~,arch,inLayer,nL,rc]=constructCase1b();

PSIC=-(0.4:0.2:1.2)*1e6;
PSIL=-(0.2:0.2:1.0)*1e6;
psiBC=combvec(PSIC,PSIL,PSIL,PSIL,PSIL,PSIL,PSIL,PSIL); %for 7 layers
nBC=size(psiBC,2);

psiX=zeros(nL,nBC);
QR=zeros(nL,nBC);
myPool=parpool(24);
parfor i=1:nBC
    
    BC=psiBC(:,i);
    psiC=BC(1);
    psiLayer=BC(2:end);
    
    psiS=psiLayer(inLayer);
    hydro=[Kr Kx psiS];
    C=linSysLF2(arch,hydro,psiC);

    G1=a.*C(:,1).*exp(a.*L)-a.*C(:,2).*exp(-a.*L);
    G0=a.*C(:,1)-a.*C(:,2);

    psiX(:,i)=psiS+(1./(L.*a.^2)).*(G1-G0);
    QR(:,i)=-Kx.*(G1-G0);
    
end

outDir=strcat('~/work/synth/case1BCon_',jobNr,filesep);
mkdir(outDir);
save(strcat(outDir,'case1bResSM'),'psiBC','psiX','QR','rc');

delete(myPool);
exit

