function MMIMatrix=MMI(X,Y,clusterlabelX,clusterlabelY)
totalnum=size(X,1);
MMIout=0;
MMIin=0;
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
proportion=number/totalnum;
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
proportion=number/totalnum;
newentropy=+proportion*0.5*log2(2*pi*exp(1)*varn)-proportion*log2(proportion);
if newentropy<0
newentropy=0;
end
entropyy=entropyy+newentropy;
end
for m=1:1:maxX
        mm=Xassign(m);
		dataX=X(clusterlabelX==mm);
		varm=var((X(clusterlabelX==mm)));
        for n=1:1:maxY
          nn=Yassign(n);
          dataY=Y(clusterlabelY==nn);
	      varn=var((Y(clusterlabelY==nn)));
          covdataA=[];
          covdataB=[];
          index=1;       
          [findrow1,findcol1]=find(clusterlabelX==mm&clusterlabelY==nn);                   
          covdataA=X(findrow1,1);
          covdataB=Y(findrow1,1);
          share=size(covdataA,1);
          if share>=1
			cluster1=sum(clusterlabelX==mm);
			cluster2=sum(clusterlabelY==nn);   
			MMIout=MMIout+share/totalnum*log2((share/(cluster1*cluster2/totalnum)));
          if var(covdataA)~=0 && var(covdataB)~=0
			covcell=cov(covdataA,covdataB);
			invSigma= pinv(covcell);
			bprime=fminbnd(@(bprime) myfuncov(bprime,invSigma(1,2),varm,varn),-sqrt(varm*varn),sqrt(varm*varn));
			pairmatrix=[varm bprime;bprime varn];
			MMIin=MMIin+0.5*log2(varm*varn/det(pairmatrix))*share/totalnum;
		  end
         end
       end
end
MMIMatrix=(MMIout+MMIin)/max(entropyx,entropyy);
