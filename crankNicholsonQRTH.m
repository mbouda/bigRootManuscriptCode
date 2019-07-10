function J=crankNicholsonQRTH(psi1,psi0,ch,A,L,kr,b2,c1,b12,b21,c12,c21,k12,k21,CnA,CpA,nL,GC,dT,dZ)

       th0=clappHornTh(psi0,ch);
       th1=clappHornTh(psi1,ch);
       B0=formPentAnBKG(b2,c1,b12,b21,c12,c21,k12,k21,CnA,CpA,nL,psi0,GC(1));
       B1=formPentAnBKG(b2,c1,b12,b21,c12,c21,k12,k21,CnA,CpA,nL,psi1,GC(2));
       QR0=L.*kr.*((A\B0)-psi0); 
       QR1=L.*kr.*((A\B1)-psi1); 
       QR=(QR0+QR1)/2;
       dTH=dT*QR./(dZ);
       J=(1e4*((th1-th0)-dTH)).^2;
       
end