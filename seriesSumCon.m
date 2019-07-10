function [c,ceq]=seriesSumCon(X,KR,sc,nL)

    X=abs(X);
    Kr=X(1:nL)*sc(1);
    ceq=1e8*(sum(Kr)-KR);
    c=[];
end 