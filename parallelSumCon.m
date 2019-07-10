function [c,ceq]=parallelSumCon(X,KR,sc)

    Kr=abs(X)*sc(1);
    ceq=1e8*(sum(Kr)-KR);
    c=[];
end