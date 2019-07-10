function b=formPentBGC(params,psiS,GC)

    a=params(:,1);
    c=params(:,3);
    d=params(:,4);
    g=params(:,7);
    h=params(:,8);
    
    back=cat(1,zeros(2,1),psiS(1:end-2).*g(3:end));
    prev=cat(1,0,psiS(1:end-1).*c(2:end));
    this=-psiS;
    next=cat(1,psiS(2:end).*d(1:end-1),0);
    forw=cat(1,psiS(3:end).*h(1:end-2),zeros(2,1));
    
    b=back+prev+this+next+forw;
    b(1)=b(1)-a(1)*GC;

end