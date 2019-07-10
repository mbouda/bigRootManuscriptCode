function [G0,G1]=gradJG(b2,c1,c12,c21,CnA,CpA,psiX,psiS,GC,nSol)
                 
    C1=repmat(c1,[1 nSol]);
    C12=repmat(c12,[1 nSol]);
    C21=repmat(c21,[1 nSol]);
    CNA=repmat(CnA,[1 nSol]);
    CPA=repmat(CpA,[1 nSol]);

    G1=cat(1,GC,...
        -(C21(1:end-1,:).*(psiX(1:end-2,:)-(1-C1(1:end-2,:)).*psiS(1:end-2,:))+(1-C1(2:end-1,:)).*psiS(2:end-1,:)-psiX(2:end-1,:))./CNA(1:end-1,:),...
        b2(end)*(psiX(end,:)-psiS(end,:)));
        
    G0=cat(1,GC-b2(1).*(psiX(1,:)-psiS(1,:)),...
        (C12(2:end,:).*(psiX(3:end,:)-(1-C1(3:end,:)).*psiS(3:end,:))+(1-C1(2:end-1,:)).*psiS(2:end-1,:)-psiX(2:end-1,:))./CPA(2:end,:),...
        zeros(1,nSol));

end
