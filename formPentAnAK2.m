function A=formPentAnAK2(b2,c4,b12,b21,c12,c21,k12,k21,Cn,Cp,nL)

    a1=cat(1,-c4(1)*k21(1)*(b2(2)-1/Cp(2)),...
        -k21(2:end-1).*(1./(b2(2:end-2).*Cp(3:end))-b21(2:end-1)),...
        k21(end)*b21(end));
    an1=cat(1,-k12(1)/(b2(2)*c4(1)),...
        -k12(2:end).*(1./(b2(3:end).*Cn(1:end-1))-b12(2:end)));
    a2=cat(1,-c4(1)*k21(1)*c12(2)/Cp(2),...
        k21(2:end-1).*c12(3:end)./(b2(2:end-2).*Cp(3:end)));
    an2=k12(2:end).*c21(1:end-1)./(b2(3:end).*Cn(1:end-1));
   
    A=diag(ones(nL,1))+diag(a1,1)+diag(an1,-1)+diag(a2,2)+diag(an2,-2);
end


