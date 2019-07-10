function J=crankNicholsonParallelTH2(psi1,psi0,ch,Kr,QC,dT,dZ,nL)

       th0=clappHornTh(psi0,ch);
       th1=clappHornTh(psi1,ch);
       
       psiC0=repmat((sum(Kr.*psi0)-QC')./sum(Kr),[nL 1]);
       psiC1=repmat((sum(Kr.*psi1)-QC')./sum(Kr),[nL 1]);
       
       A0=diag(Kr);
       b0=Kr/2.*(psiC0+psi0);
       A1=diag(Kr);
       b1=Kr/2.*(psiC1+psi1);
              
       QR0=-(Kr*2).*(psiC0-A0\b0);
       QR1=-(Kr*2).*(psiC1-A1\b1);
       QR=(QR0+QR1)/2;
       
       dTH=dT*-QR./(dZ);
       J=(1e6*((th1-th0)-dTH)).^2;
       
end