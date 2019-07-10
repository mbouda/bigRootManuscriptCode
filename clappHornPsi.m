function P=clappHornPsi(th,ch)

    szT=size(th);
    if szT(2)>1
        
        THS=repmat(ch.ths',[szT(1) 1]);
        PSAT=repmat(ch.pSat',[szT(1) 1]);
        B=repmat(ch.B',[szT(1) 1]);
        P=PSAT.*(th./THS).^(-B);
    else
        P=ch.pSat.*(th./ch.ths).^(-ch.B);
    end

end