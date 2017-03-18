step=2000;
range=[-1 1];

for n=6:5:101
    for j=1:n
            xm(j,(n-1)/5)=-1+2*(j-1)/(n-1);
            fref(j,(n-1)/5)=1/(1+25*xm(j,(n-1)/5)^2);
    end
    
    %Lagrangian
    for m=1:step
    x=range(1)+(range(2)-range(1))/step*m;
        for j=1:n
            for i=1:n
                lelement(j,i)=(x-xm(i,(n-1)/5))/(xm(j,(n-1)/5)-xm(i,(n-1)/5));
            end
            l(j,(n-1)/5)=prod(lelement(j,1:j-1))*prod(lelement(j,j+1:n));
        end
        L(m,(n-1)/5)=dot(fref(1:n,(n-1)/5),l(1:n,(n-1)/5));
    end
        
    %Linear
    for m=1:step
    x=range(1)+(range(2)-range(1))/step*m;
        for j=2:n-1
            if (or(x<=xm(j-1,(n-1)/5),x>=xm(j+1,(n-1)/5)))
                fai(j,(n-1)/5)=0;
            else
                if ((x>=xm(j-1,(n-1)/5)) && (x<=xm(j,(n-1)/5)))
                    fai(j,(n-1)/5)=(x-xm(j-1,(n-1)/5))/(xm(j,(n-1)/5)-xm(j-1,(n-1)/5));
                else
                    fai(j,(n-1)/5)=(x-xm(j+1,(n-1)/5))/(xm(j,(n-1)/5)-xm(j+1,(n-1)/5));
                end
            end
        end
        if ((x>=xm(1,(n-1)/5))&&(x<=xm(2,(n-1)/5)))
            fai(1,(n-1)/5)=(x-xm(2,(n-1)/5))/(xm(1,(n-1)/5)-xm(2,(n-1)/5));
        end
        if ((x>=xm(n-1,(n-1)/5))&&(x<=xm(n,(n-1)/5)))
            fai(n,(n-1)/5)=(x-xm(n-1,(n-1)/5))/(xm(n,(n-1)/5)-xm(n-1,(n-1)/5));
        else
            fai(n,(n-1)/5)=0;
        end
        I(m,(n-1)/5)=dot(fref(1:n,(n-1)/5),fai(1:n,(n-1)/5));
    end
        
    %Spline
    for j=3:n-2
        g(j)=3*n/4*(-fref(j-1,(n-1)/5)+fref(j+1,(n-1)/5));
    end
    g(2)=3*n/4*(-fref(1,(n-1)/5)+fref(3,(n-1)/5))-50/676/2;
    g(n-1)=3*n/4*(-fref(n-2,(n-1)/5)+fref(n,(n-1)/5))+50/676/2;
    M(2:n-1,(n-1)/5)=(diag(diag(2*ones(n-2)),0)+diag(diag(0.5*ones(n-3)),1)+diag(diag(0.5*ones(n-3)),-1))\g(2:n-1)';
    M(1,(n-1)/5)=50/676; M(n,(n-1)/5)=-50/676;
    
    for m=1:step
    x=range(1)+(range(2)-range(1))/step*m;
        for j=2:n-1
            if (x>=xm(j-1,(n-1)/5))&&(x<=xm(j,(n-1)/5))
                alpha(j,(n-1)/5)=(1+2*(x-xm(j,(n-1)/5))/(xm(j-1,(n-1)/5)-xm(j,(n-1)/5)))*((x-xm(j-1,(n-1)/5))/(xm(j,(n-1)/5)-xm(j-1,(n-1)/5)))^2;
                beta(j,(n-1)/5)=(x-xm(j,(n-1)/5))*((x-xm(j-1,(n-1)/5))/(xm(j,(n-1)/5)-xm(j-1,(n-1)/5)))^2;
            else
                if ((x>=xm(j,(n-1)/5)) && (x<=xm(j+1,(n-1)/5)))
                    alpha(j,(n-1)/5)=(1+2*(x-xm(j,(n-1)/5))/(xm(j+1,(n-1)/5)-xm(j,(n-1)/5)))*((x-xm(j+1,(n-1)/5))/(xm(j,(n-1)/5)-xm(j+1,(n-1)/5)))^2;
                    beta(j,(n-1)/5)=(x-xm(j,(n-1)/5))*((x-xm(j+1,(n-1)/5))/(xm(j,(n-1)/5)-xm(j+1,(n-1)/5)))^2;
                else
                    alpha(j,(n-1)/5)=0;
                    beta(j,(n-1)/5)=0;
                end
            end
        end
        if ((x>=xm(n-1,(n-1)/5))&&(x<=xm(n,(n-1)/5)))
            alpha(n,(n-1)/5)=(1+2*(x-xm(n,(n-1)/5))/(xm(n-1,(n-1)/5)-xm(n,(n-1)/5)))*((x-xm(n-1,(n-1)/5))/(xm(n,(n-1)/5)-xm(n-1,(n-1)/5)))^2;
            beta(n,(n-1)/5)=(x-xm(n,(n-1)/5))*((x-xm(n-1,(n-1)/5))/(xm(n,(n-1)/5)-xm(n-1,(n-1)/5)))^2;
        else
            alpha(n,(n-1)/5)=0; beta(n,(n-1)/5)=0;
        end
        if ((x>=xm(1,(n-1)/5)) && (x<=xm(2,(n-1)/5)))
            alpha(1,(n-1)/5)=(1+2*(x-xm(1,(n-1)/5))/(xm(2,(n-1)/5)-xm(1,(n-1)/5)))*((x-xm(2,(n-1)/5))/(xm(1,(n-1)/5)-xm(2,(n-1)/5)))^2;
            beta(1,(n-1)/5)=(x-xm(1,(n-1)/5))*((x-xm(2,(n-1)/5))/(xm(1,(n-1)/5)-xm(2,(n-1)/5)))^2;
        else
            alpha(1,(n-1)/5)=0;
            beta(1,(n-1)/5)=0;
        end
        S(m,(n-1)/5)=dot(fref(1:n,(n-1)/5),alpha(1:n,(n-1)/5))+dot(M(1:n,(n-1)/5),beta(1:n,(n-1)/5));
    end
    
end
for m=1:step
    x=range(1)+(range(2)-range(1))/step*m;
    f(m)=1/(1+25*x^2);
end
for j=1:(n-1)/5
    errorL(:,j)=L(:,j)-f';
    errorI(:,j)=I(:,j)-f';
    errorS(:,j)=S(:,j)-f';
end