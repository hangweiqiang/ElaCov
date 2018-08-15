function  [SS, S0] = elaCORE(X, C)

% C: a cx1 vector of thresholds 
% SS: pxpxc matrix.
% S0: pxp

[n,p] = size(X);
X = X - repmat(mean(X), n, 1);
X = X./repmat(std(X), n, 1)*sqrt(n/(n-1));

sd = eye(p); S0 = eye(p)-100;
for i = 1:p
    for j = i+1:p
        e = X(:,i).*X(:,j);
        mij = mean(e);
        S0(i,j) = mij;
        sd(i,j) = sqrt(sum((e-mij).^2)/(n-1));
%        S0(j,i) = S0(i,j);
%        sd(j,i) = sd(i,j);
    end
end
sd = sd/sqrt(n);


stmp = S0(S0>-10);
SS = ones(p,p,length(C));

%process stmp
stmp=sort(stmp);
csum=cumsum(stmp);
rsum=cumsum(stmp(end:-1:1));
rsum=rsum(end:-1:1);
 for ic = 1:length(C)
     Csd = C(ic)*sd; 
     for i = 1:p-1
         for j = i+1:p
             if(Csd(i,j)<0)
                 SS(i,j,ic)=0;
                 SS(j,i,ic)=0;
             else
                First=-Csd(i,j)+S0(i,j);
                Second=Csd(i,j)+S0(i,j);
                pl=0;pr=length(stmp);
                while(pl<pr)
                    mid=floor((pl+pr+1)/2);
                    if(stmp(mid)<=Second) pl=mid;
                    else pr=mid-1;
                    end
                end
                SecondPos=pl;
                pl=1;pr=length(stmp)+1;
                while(pl<pr)
                    mid=floor((pl+pr)/2);
                    if(stmp(mid)>=First) pr=mid;
                    else pl=mid+1;
                    end
                end
                FirstPos=pl;
                if(FirstPos~=length(stmp)+1&&SecondPos~=0&&FirstPos<=SecondPos)
                    SumIij=SecondPos-FirstPos+1;
                    SumStmpIij=rsum(FirstPos)+csum(SecondPos)-rsum(1);
                    SS(i,j, ic)=SumStmpIij/SumIij;
                    SS(j,i,ic)=SS(i,j,ic);
                else
                    SS(i,j,ic)=0;
                    SS(j,i,ic)=0;
                end
             end             
         end
     end
 end
end
