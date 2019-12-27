clc;
warning off;
%clear all;
data=load('exp.txt');

[clusternum,clusterassign]=dynamicprogramming(data,3);
[CovMatrix,InvMatrix]=MDIin(data,clusterassign,clusternum,0.3);
MDIresult=MDIcombine(data,clusterassign,clusternum,InvMatrix,0.1);
X=data(:,1);
Y=data(:,2);
clusterlabelX=clusterassign(:,1);
clusterlabelY=clusterassign(:,2);
MMIMatrix=MMI(X,Y,clusterlabelX,clusterlabelY);
