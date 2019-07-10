function [c,ceq]=bigRootQRGCSubSumCon(X,sc,Kr,nL)

    X=abs(X);

    kr=X(1:nL)/sc(1);
    L=X(2*nL+(1:nL))/sc(3);    

    KrS=kr.*L;
    
    ceq=1e8*(sum(KrS)-sum(Kr));
    c=[];   
end