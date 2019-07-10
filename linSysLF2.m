function C=linSysLF2(arch,hydro,psiC)
%linSysLF2(arch,hydro,psiC)
    % inputs: 
    %   arch=[parents L M]
    %   hydro=[Kr Kx psiS]
%This function solves the Landsberg-Fowkes problem on a network domain
%Allowing for both junctions (bifurcations of tree, where three links meet)
%and joints (nodes, where 2 links meet) 

    a=sqrt(hydro(:,1)./hydro(:,2));
    psiS=hydro(:,3);
    L=arch(:,2);
    M=arch(:,3);
    nL=size(arch,1);
    parents=arch(:,1);
    
    %% Kx ratios
%     ratioKx=hydro(:,2)*(1./hydro(:,2))';
    
    
    %% make base coefficient parts
    A=hydro(:,2).*(a.*tanh(a.*L)); 
    B=((1./a-exp(-a.*L)./(a.*cosh(a.*L))))./hydro(:,2);
    C=exp(-a.*L).*(tanh(a.*L)+1);
    D=-1./cosh(a.*L);

    %% combine them into coefficients according to topology
    
    nJunc=M(1)-1;
    nNode=length(unique(parents))-1; %need to use sparseUnique!
    nJoint=nNode-nJunc;
    
    [oP,nodeType]=occur(parents);
    junctP=unique(parents(oP==2)); %need to use sparseUnique!
    %can we just do: oP(nodeType==2)? Will scramble order? if so, will fail?
    jointP=parents(oP==1);
    jointP(1)=[];
    
    juncs=zeros(nJunc,3);
    for i=1:nJunc
        juncs(i,:)=cat(2,junctP(i),find(parents==junctP(i))');
    end
    
    joints=zeros(nJoint,2);
    for i=1:nJoint
        joints(i,:)=cat(2,jointP(i),find(parents==jointP(i)));
    end
    
    %G holds parents information
    G=zeros(nL);
    for i=nL:-1:2
        G(i,parents(i))=1;
        distal=G(:,i)>0;
        G(distal,parents(i))=G(distal,i)+1;
    end
    DL=G>0;
    DL=DL(:,any(DL));
    
    selfC=ones(nNode,1);
    proxC=zeros(nNode);
    distC=zeros(nNode);
    distRH=zeros(nNode,1);
    for i=1:nNode
        
        proxC(i,nodeType(:,1)==parents(nodeType(i,1)))=D(nodeType(i,1));

        if nodeType(i,2)==1 %joint
            selfC(i)=selfC(i)+B(nodeType(i,1))*A(joints(joints(:,1)==nodeType(i,1),2));
        elseif nodeType(i,2)==2 %junction
            selfC(i)=selfC(i)+B(nodeType(i,1))*sum(A(juncs(juncs(:,1)==nodeType(i,1),2:3)'));
        else
            throw(MException('linSysLF2:badNode','Links should only occur once or twice in parents vector'));
        end
        
        for j=find(DL(:,i))'; %DL holds link info, not node info...
            k=j;
            pathIJ=[];
            while k>nodeType(i,1)
                pathIJ=cat(1,k,pathIJ);
                k=parents(k);
            end
            pathIJ(end)=[]; %remove, C not applicable on first distal link
            E=B(nodeType(i,1))*A(j)*prod(C(pathIJ)); %we actually want C_2 and have C_4
            %use pathIJ to pick Kx ratios... and use those, also...
            if any(pathIJ)
                distC(i,nodeType(:,1)==parents(j))=distC(i,nodeType(:,1)==parents(j))+E;
            end
            distRH(i)=distRH(i)+E*psiS(j); %execute once for each distal link to ith node
        end
        
    end
    proxRH=psiS(nodeType(:,1)).*(1+D(nodeType(:,1)));

    %% make and solve system for all psiJ
    X=proxC+diag(selfC)+distC;
    b=cat(1,-D(1)*psiC,zeros(nNode-1,1))+proxRH+distRH;
    psiJ=X\b;

    %% Use psiJ to find all Gradients
    %distribute psiJ to links
    
    pJL=zeros(nL,1);
    for i=1:nNode
        if nodeType(i,2)==2
            pJL(juncs(juncs(:,1)==nodeType(i,1),2:3))=psiJ(i);
        else
            pJL(joints(joints(:,1)==nodeType(i,1),2))=psiJ(i);
        end
    end
    pJL(1)=psiC;
    
    %re-formulate in Q...
    QJ=A.*(pJL-psiS);
    QB=zeros(nL,1);
    for i=nL:-1:1
        j=parents==i;
        if any(j)
            QB(i)=sum(QJ(j));
            QJ(i)=QJ(i)+C(i)*QB(i);
        end
    end
    %is the order wrong?
    GB=QB./hydro(:,2);
    %% use psiJ and G to find all C
    C=zeros(nL,2);
    C(:,1)=(pJL-psiS+exp(-a.*L).*GB./a)./(exp(a.*L)+exp(-a.*L));
    C(:,2)=C(:,1)-GB./a;
    
end