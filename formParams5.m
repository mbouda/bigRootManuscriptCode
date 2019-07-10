function [b2,c1,c2,c3,c4,c5,b12,b21,c12,c21,k12,k21,CnA,CpA,CnB,CpB]=formParams5(kr,Kx,L)

    a2=kr./Kx;
    a=sqrt(a2);
    b=a.*L;
    b2=a2.*L;
    
    c1=sinh(b)./(b);
    c2=(1-cosh(b))./(b2);
    c3=tanh(b)./(b);
    c4=(sech(b(1))-1)./b2(1);
    c5=tanh(b/2)./(b);

    b21=b2(2:end)./b2(1:end-1);
    b12=b2(1:end-1)./b2(2:end);
    c21=c1(2:end)./c1(1:end-1);
    c12=c1(1:end-1)./c1(2:end);
    
    k21=Kx(2:end)./Kx(1:end-1);
    k12=Kx(1:end-1)./Kx(2:end);

    CnA=c2(2:end)+k21.*c21.*c2(1:end-1);
    CpA=c2(1:end-1)+k12.*c12.*c2(2:end);
    
    CpB=(c1(1:end-1)./c2(1:end-1)+c1(2:end)./(c2(2:end).*k12));
    CnB=(c1(2:end)./c2(2:end)+c1(1:end-1)./(c2(1:end-1).*k21));

end