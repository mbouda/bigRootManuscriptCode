function A=formPentAGC(params)

    a=params(:,1);
    b=params(:,2);
    c=params(:,3);
    d=params(:,4);
    e=params(:,5);
    f=params(:,6);
    g=params(:,7);
    h=params(:,8);
    
    lld=g(3:end)+e(3:end);
    ld=a(2:end)+c(2:end)-e(2:end);
    md=b-a-1; md(1)=b(1)-1;
    ud=d(1:end-1)+f(1:end-1)-b(1:end-1);
    uud=h(1:end-2)-f(1:end-2);

    A=diag(lld,-2)+diag(ld,-1)+diag(md)+diag(ud,1)+diag(uud,2);


end
