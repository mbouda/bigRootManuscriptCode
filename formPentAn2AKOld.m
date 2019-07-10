function A=formPentAn2AKOld(c2,c3,c5,Cn,Cp,k12,k21,nL)

    a1=cat(1,(c5(1:end-2)./(c2(2:end-1).*Cp(2:end))-c5(1:end-2)./c5(2:end-1)),...
        -c5(end-1)/c3(end));
    an1=cat(1,-c5(2)/c5(1),...
        c5(3:end-1)./(c2(2:end-2).*Cn(1:end-2))-c5(3:end-1)./c5(2:end-2),...
        c3(end)./(c2(end-1).*Cn(end-1))-c3(end)./c5(end-1));
    a2=c5(1:end-2)./(c2(3:end).*k12(2:end).*Cp(2:end));
    an2=cat(1,c5(3:end-1)./(c2(1:end-3).*k21(1:end-2).*Cn(1:end-2)),...
        c3(end)./(c2(end-2).*k21(end-1).*Cn(end-1)));
    
    A=diag(ones(nL,1))+diag(a1,1)+diag(an1,-1)+diag(a2,2)+diag(an2,-2);

end