for n=1:50
    H{n}=hilb(n);
    x=ones(n,1);
    b{n}=H{n}*x;
    
    %Gaussian
    yGauss(1:n,n)=H{n}\b{n};
    errG(1:n,n)=abs(yGauss(1:n,n)-[1]);
    
    %Cholesky
    if (n<=12)
        %Cholesky only applies to n<=12 normally.
        L=chol(H{n});
        %t=L'\b{n};
        yChol(1:n,n)=L\(L'\b{n});
        errC(1:n,n)=abs(yChol(1:n,n)-[1]);
    end
    
    %Condition Number
    Condition(n)=cond(H{n},2);
    
    %Balanced
%    for i=1:n
%        S(i)=max(H{n}(i,:));
%        D(i,i)=S(i)^-1;
%    end
%    DA{n}=D*H{n};
%    ConditionB(n)=cond(DA{n},2);
%    yBal(1:n,n)=DA{n}\(D*b{n});
%    errB(1:n,n)=abs(yBal{n}-[1]);
    
    %Tikhonov Regularization
%    alpha=0.001;
%    for i=1:20
%        alpha=alpha/10;
%        yTR{i+3}(1:n,n)=(alpha*eye(n)+H{n}*(H{n})')\((H{n})'*b{n});
%        errT{i+3}(1:n,n)=abs(yTR{i+3}(1:n,n)-[1]);
%    end
    alpha=1e-12;
    yTR(1:n,n)=(alpha*eye(n)+H{n}*(H{n})')\((H{n})'*b{n});
    errT(1:n,n)=abs(yTR(1:n,n)-[1]);
    
    %SOR
    ySOR(1:n,n)=0;
    SORD=diag(diag(H{n}));
    SORL=tril(H{n},-1);
    SORU=triu(H{n},1);
    ohm=0.1; %abs(2/(1+sqrt(1-vrho((SORD-SORL)\SORU))));
    for t=1:2000
        if n>=2
            ySOR(1:n,n)=(SORD+ohm*SORL)\(((1.0-ohm)*SORD-ohm*SORU)*ySOR(1:n,n)+ohm*b{n});
        else
            ySOR(1,1)=1;
        end
    end
    errSOR(1:n,n)=abs(ySOR(1:n,n)-[1]);
    
    %Conjugate Grandient
    yCG(1:n,n)=0;
    r=b{n}-H{n}*yCG(1:n,n);
    p=r;
    for k=0:100
        step=dot(r,p)/dot(H{n}*p,p);
        if dot(H{n}*p,p)==0 break
        end
        yCG(1:n,n)=yCG(1:n,n)+step*p;
        r=b{n}-H{n}*yCG(1:n,n);
        beta=-dot(r,H{n}*p)/dot(p,H{n}*p);
        p=r+beta*p;
    end
    errCG(1:n,n)=abs(yCG(1:n,n)-[1]);
    
end