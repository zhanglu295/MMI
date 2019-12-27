function [CovMatrix,InvMatrix]=MDIin(expressdata,clusterassign,clusternum,shrink)
numbergene=size(clusternum,1);
CovMatrix=[];
aggregate(1)=0;
for i=2:1:numbergene
    aggregate(i)=aggregate(i-1)+clusternum(i-1,1);
end
for i=1:1:numbergene
X=expressdata(:,i);
clusterlabelX=clusterassign(i,:);
for j=i:1:numbergene
Y=expressdata(:,j);
clusterlabelY=clusterassign(j,:);
for m=1:1:clusternum(i,1)
    row=aggregate(i)+m;
    indexm=find(clusterlabelX==m);
    index_nom=find(clusterlabelX~=m);
    mean_m=mean(X(indexm));
    datam(indexm)=X(indexm);
    datam(index_nom)=mean_m;
    varm=var(datam);
for n=1:1:clusternum(j,1)
    col=aggregate(j)+n;
    indexn=find(clusterlabelY==n);
    index_non=find(clusterlabelY~=n);
    mean_n=mean(Y(indexn)); 
    datan=[];
    datan(indexn)=Y(indexn);
    datan(index_non)=mean_n;
    varn=var(datan);
    if row==col
    CovMatrix(row,col)=varm;
    end
    if row~=col    
     [findrow1,findcol1]=find(clusterlabelX==m&clusterlabelY==n);                    covdataA=X(findcol1,1);
     covdataB=Y(findcol1,1);
     share=size(covdataA,1);
     covcell=cov(datam,datan);
     if share==0
     CovMatrix(row,col)=0;
     else
     invSigma= pinv(covcell);
     bprime=fminbnd(@(bprime) myfuncov(bprime,invSigma(1,2),varm,varn),-sqrt(varm*varn),sqrt(varm*varn));
	CovMatrix(row,col)=bprime;
   end
   end
end
end
end
end
[aa,bb]=size(CovMatrix);
for i=1:1:aa
    for j=i+1:1:bb
        CovMatrix(j,i)= CovMatrix(i,j);
    end
end
InvMatrix=inv(CovMatrix'*CovMatrix+shrink*eye(aa))*CovMatrix';
dlmwrite('MDI_covin.txt',CovMatrix,'delimiter','\t');
dlmwrite('MDI_invcovin.txt',InvMatrix,'delimiter','\t');
end
