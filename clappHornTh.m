function th=clappHornTh(psi,ch)

    th=ch.ths.*(psi./ch.pSat).^(-1./ch.B);

end