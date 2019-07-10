function J=crankNicholsonFreeTH(psi1,psi0,ch,A,coeffs,KR,GC,dT,dZ)

       th0=clappHornTh(psi0,ch);
       th1=clappHornTh(psi1,ch);
       
       B0=formPentBGC(coeffs,psi0,GC(1));
       B1=formPentBGC(coeffs,psi1,GC(2));
       
       QR0=KR.*((A\B0)-psi0); 
       QR1=KR.*((A\B1)-psi1); 
       QR=(QR0+QR1)/2;
       dTH=dT*QR./(dZ);
       J=(1e4*((th1-th0)-dTH)).^2;
       
end