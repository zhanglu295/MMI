function MDIresult=MDIcombine(data,clusterassign,clusternum,MDIin_invMatrix,shink)
Data_ori=data;
clusterlabel=clusterassign';
numrow=1;
[x1,y1]=size(Data_ori);
MDIout_Matrix=[];
alpha=0.1;
lambda=x1/(1-alpha)*alpha;
totalcolumn=sum(clusternum);

aggregate(1)=0;
for i=2:1:y1
	aggregate(i)=aggregate(i-1)+clusternum(i-1,1);
end
aggregate_all(1)=0;
for i=2:1:y1
    aggregate_all(i)=aggregate_all(i-1)+clusternum(i-1,1);
end
Databin=zeros(x1,totalcolumn);
Wi=[];
for i=1:1:y1
    for j=1:1:clusternum(i,1)
     index=find(clusterlabel(:,i)==j);
     Databin(index,aggregate(i)+j)=1;
	 Wi(1,aggregate(i)+j)=(1/(x1+lambda))*(lambda/(clusternum(i,1))+size(index,1));
end
end

for i=1:1:y1
       numi=clusternum(i,1);
	  for j=i:1:y1
       numj=clusternum(j,1);
       for m=1:1:numi
           row=aggregate(i)+m;
         for n=1:1:numj
           col=aggregate(j)+n;
           [findrow,findcol]=find(Databin(:,row)==1&Databin(:,col)==1);
           if i==j & m~=n
           MDIout_Matrix(row,col)=-Wi(1,row)*Wi(1,col);
           else
		   MDIout_Matrix(row,col)=(1/(x1+lambda))*(lambda/((clusternum(i,1))*(clusternum(j,1)))+size(findrow,1))-Wi(1,row)*Wi(1,col);
           end
        end
       end
     end
end
[a,b]=size(MDIout_Matrix);
for i=1:1:a
    for j=i:1:b
        MDIout_Matrix(j,i)= MDIout_Matrix(i,j);
    end
end
incov=inv(MDIout_Matrix'*MDIout_Matrix+shink*eye(a))*MDIout_Matrix';
combineout=zeros(y1,y1);
combinein=zeros(y1,y1);
for i=1:1:y1
     numi=clusternum(i,1);
     for j=i:1:y1
      numj=clusternum(j,1);
       for m=1:1:numi
           row=aggregate(i)+m;
         for n=1:1:numj
           col=aggregate(j)+n;
           combineout(i,j)=combineout(i,j)+abs(incov(row,col));
         end
       end
    end
end
for i=1:1:y1
    for j=i+1:1:y1
        combineout(j,i)=combineout(i,j);
    end
end
PCallout=sum(sum(combineout));
for i=1:1:y1
     numi=clusternum(i,1);
     for j=i:1:y1
      numj=clusternum(j,1);
       for m=1:1:numi
           row=aggregate(i)+m;
         for n=1:1:numj
           col=aggregate(j)+n;
           combinein(i,j)=combinein(i,j)+abs(MDIin_invMatrix(row,col));
         end
       end
    end
end
for i=1:1:y1
    for j=i+1:1:y1
        combinein(j,i)=combinein(i,j);
    end
end
PCallin=sum(sum(combinein));
for i=1:1:y1
     numi=clusternum(i,1)-1;
     for j=i+1:1:y1
      numj=clusternum(j,1)-1;
    X=Data_ori(:,i);
    Y=Data_ori(:,j);
    clusterlabelX=clusterlabel(:,i);
    clusterlabelY=clusterlabel(:,j); 
    MDIout=0;
    MDIin=0;
	MDIout=combineout(i,j)-combineout(i,:)*combineout(:,j)/PCallout;
	MDIin=combinein(i,j)-combinein(i,:)*combinein(:,j)/PCallin;
	maxX=size(unique(clusterlabelX),1);
	maxY=size(unique(clusterlabelY),1);
	Xassign=unique(clusterlabelX);
	Yassign=unique(clusterlabelY);
	entropyx=0;
	entropyy=0;
	for m=1:1:maxX
	mm=Xassign(m);
	number=size(X(clusterlabelX==mm),1);
	varm=var((X(clusterlabelX==mm)));
	proportion=number/x1;
	newentropy=proportion*0.5*log2(2*pi*exp(1)*varm)-proportion*log2(proportion);
	if newentropy<0
	newentropy=0;
	end
	entropyx=entropyx+newentropy;
	end
	for n=1:1:maxY
	nn=Yassign(n);
	number=size(Y(clusterlabelY==nn),1);
	varn=var((Y(clusterlabelY==nn)));
	proportion=number/x1;
	newentropy=+proportion*0.5*log2(2*pi*exp(1)*varn)-proportion*log2(proportion);
	if newentropy<0
	newentropy=0;
	end
	entropyy=entropyy+newentropy;
	end
	MDIresult(numrow,1)=i;
	MDIresult(numrow,2)=j;
	MDIresult(numrow,3)=(MDIin+MDIout)/max(entropyx,entropyy);
        numrow=numrow+1;
    end
end
dlmwrite('MDI.result',MDIresult,'delimiter','\t');
