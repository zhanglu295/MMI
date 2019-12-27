function [clusternum,clusterassign]=dynamicprogramming(Data,maxallow)
[row col]=size(Data);
clusterassign=[];
for r=1:1:col
BICarray=[];
T=[];
for i=1:1:maxallow
obj = gmdistribution.fit(Data(:,r),i,'Replicates',10);
BICarray(i)=obj.BIC;
end
[BICvalue,I]=sort(BICarray);
clusternum(r,1)=I(1)
[data,IX]=sort(Data(:,r));
%Create Variance Matrix
%Square
datasquare=data.^2;
SquareMatrix=[];
for i=1:1:row
    for j=i+1:1:row
        SquareMatrix(i,j)=sum(datasquare(i:j));
    end
end 
MeanMatrix=[];
for i=1:1:row
    for j=i+1:1:row
        MeanMatrix(i,j)=mean(data(i:j));
    end 
end 
Variance=zeros(row,row);
for i=1:1:row
    for j=i+1:1:row
       sigma=SquareMatrix(i,j)-(j-i+1)* MeanMatrix(i,j)^2;
       if sigma<=0
           sigma=0.00001;
       end
       Variance(i,j)=-(j-i+1)*log(sigma)/2;
    end
end


for i=1:1:row
    for j=1:1:clusternum(r,1)-1
        T(i,j)=Variance(j,i);
    end
end
    T(:,clusternum(r,1))=-10000;


for k=2:clusternum(r,1)
    for i=2:1:row
        for j=1:1:i-1
            T(i,k)=max(T(j,k-1)+Variance(j+1,i),T(i,k));
            if T(i,k)==T(j,k-1)+Variance(j+1,i)
                B(i,k)=j;
            end
        end
    end
end
maxvalue=max(T(row,:));
[x,y]=find(T(row,:)==maxvalue);
Pos=[];
Pos(1,1)=B(row,y);
  k=y;
  index=2;
   while Pos(1,index-1)~=0
   Pos(1,index)=B(Pos(1,index-1),k-1);
   k=k-1;
   index=index+1;
   end
   Pos=[row Pos];
   for t=size(Pos,2):-1:2
	clusterassign(IX(Pos(1,t)+1:Pos(1,t-1),1),r)=t-1;
   end
end
