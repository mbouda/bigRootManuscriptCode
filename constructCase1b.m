function [L,Kr,Kx,a,c1,c2,arch,inLayer,nL,rc]=constructCase1b()

    parents=[0 1 2 3 4 5 6]';
    M=ones(7,1);
    nL=size(parents,1);
    L=repmat(0.05,[nL 1]);

    kxA=1e-4;
    krA=1e-9;
    rA=0.5e-3;
    bA=0.1e-3;

    Kr=repmat(pi*(2*rA-bA)*krA/bA,[nL 1]);
    Kx=repmat(kxA*pi*(rA-bA)^2,[nL 1]);

    rc=rand(nL,3)+0.5; %random coefficients, to avoid degeneracy
    L=L.*rc(:,1);
    Kr=Kr.*rc(:,2);
    Kx=Kx.*rc(:,3);
    
    arch=[parents L M];
    a=sqrt(Kr./Kx);
    c1=sinh(a.*L)./(a.*L);
    c2=(1-cosh(a.*L))./(a.^2.*L);

    inLayer=(1:7)';
    
end
