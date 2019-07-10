function J=crankNicholsonSeriesTH2(psi1,psi0,ch,A,Kr,QC,sc,dT,dZ)

       th0=clappHornTh(psi0,ch);
       th1=clappHornTh(psi1,ch);
       
       B0=Kr.*psi0;
       B0(1)=B0(1)-QC*sc(3);
       B1=Kr.*psi1;
       B1(1)=B1(1)-QC*sc(3);
       
       QR0=-Kr.*((A\B0)-psi0)/sc(3); 
       QR1=-Kr.*((A\B1)-psi1)/sc(3); 
       QR=(QR0+QR1)/2;
       
       dTH=dT*-QR./(dZ);
       J=(1e4*((th1-th0)-dTH)).^2;
       
end