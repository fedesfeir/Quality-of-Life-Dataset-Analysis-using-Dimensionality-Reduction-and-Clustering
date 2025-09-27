QLCit24
rng(09042002)
rowp = randperm(178);
Xp=X(rowp, :);
Xr=Xp(1:60,:);

%standarize data
n=size(Xr,1)
J=size(Xr,2)
Jc=eye(n,n)-(1/n)*ones(n,n)

%coveriance matrix
Sc=(1/60)*Xr'*Jc*Xr;
D=diag(diag(Sc).^0.5);
Z=Jc*Xr*D^(-1);

%correlation matrix
Sz=(1/60)*Z'*Z;

%2. Cdistv=pdist(Z,"euclidean");
distv=pdist(Z,"euclidean");
Dist= squareform(distv);

%3. Compute the WSPP;
[UOtt2,bOtt2, aOtt2, fOtt2,iterOtt2]=WSPP(Dist, 2, 50);
[pf2,Dw2,Db2] = psF(Z,UOtt2)
[UOtt3,bOtt3, aOtt3, fOtt3,iterOtt3]=WSPP(Dist, 3, 50);
[pf3,Dw3,Db3] = psF(Z,UOtt3)
sum(UOtt2)
sum(UOtt3)
[UOtt4,bOtt4, aOtt4, fOtt4,iterOtt4]=WSPP(Dist, 4, 50);
[pf4,Dw4,Db4] = psF(Z,UOtt4)
sum(UOtt4)

%4. Compute the WSP;
[UOtt2,DbOtt2, DwOtt2, fOtt2,iterOtt2]=WSP(Dist, 2, 50);
[pf2,Dw2,Db2] = psF(Z,UOtt2)
[UOtt3,DbOtt3, DwOtt3, fOtt3,iterOtt3]=WSP(Dist, 3, 50);
[pf3,Dw3,Db3] = psF(Z,UOtt3)
[UOtt4,DbOtt4, DwOtt4, fOtt4,iterOtt4]=WSP(Dist, 4, 50);
[pf4,Dw4,Db4] = psF(Z,UOtt4)

%5. Compute the PD;
[UOtt2,DbOtt2, DwOtt2, fOtt2,iterOtt2]=PD(Dist, 2, 50);
[UOtt3,DbOtt3, DwOtt3, fOtt3,iterOtt3]=PD(Dist, 3, 50);
[UOtt4,DbOtt4, DwOtt4, fOtt4,iterOtt4]=PD(Dist, 4, 10);
[pf2,Dw2,Db2] = psF(Z,UOtt2);
[pf3,Dw3,Db3] = psF(Z,UOtt3);
[pf4,Dw4,Db4] = psF(Z,UOtt4);
Q2 = UOtt2*DbOtt2*UOtt2'+UOtt2*DwOtt2*UOtt2'-diag(diag(UOtt2*DwOtt2*UOtt2'));
Q3 = UOtt3*DbOtt3*UOtt3'+UOtt3*DwOtt3*UOtt3'-diag(diag(UOtt3*DwOtt3*UOtt3'));
Q4 = UOtt4*DbOtt4*UOtt4'+UOtt4*DwOtt4*UOtt4'-diag(diag(UOtt4*DwOtt4*UOtt4'));

%verify ultramtric
verultrametrica(Q2)
verultrametrica(Q3)
verultrametrica(Q4)
Z3 = linkage(Q3);
figure();
dendrogram(Z3,0)

%2. Compute PCA to determine the number of Principal Components;
[A,L]=eigs(Sz,9);
A2=A(:,1:2);
A2r=rotatefactors(A2);

%3. Compute K-means on the number of component identified in the step 1;
Y=Z*A2r;
[vc2,Ym2]=kmeans(Y,2,'Replicates',100);
[vc3,Ym3]=kmeans(Y,3,'Replicates',100);

%4. Identify the best K with pF;
I2=eye(2)
U2=I2(vc2,:);
I3=eye(3)
U3=I3(vc3,:);
[pf, Dw, Db]=psF(Y, U2);
psF2=psF(Y,U2)
[pf, Dw, Db]=psF(Y, U3);
psF3=psF(Y,U3)
silhouette(Y,vc2);
s2=silhouette(Y,vc2)
mean(s2)
silhouette(Y,vc3)
s3=silhouette(Y,vc3)
mean(s3)
figure;
gscatter(Y(:,1), Y(:,2), vc2, 'rgb', 'o^s', 8);
hold on;
plot(Ym2(:,1), Ym2(:,2), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('PC1');
ylabel('PC2');
title('Tandem Analysis');
legend('Cluster 1', 'Cluster 2', 'Centroids');
grid on;
hold off;

%5. Compute Reduced K-mean;
[Urkm,Arkm, Yrkm,frkm,inrkm]=REDKM(Z, 2, 2, 20);
pf2=psF(Yrkm,Urkm)
[Urkm3,Arkm3, Yrkm3,frkm3,inrkm3]=REDKM(Z, 2, 3, 20);
pf3=psF(Yrkm3,Urkm3)
U2'*Urkm  %
U3'*Urkm3
figure;
gscatter(Y(:,1), Y(:,2), Urkm,'rgb', 'o^s', 8);
hold on;
title('Reduced K-Means Clustering');
legend('Cluster 1', 'Cluster 2');
grid on;
hold off;

%6. Compute Factorial K-means;
[Ufkm,Afkm, Yfkm,ffkm,infkm]=FKM(Z, 2, 2, 20);
U2'*Ufkm
[Ufkm3,Afkm3, Yfkm3,ffkm3,infkm3]=FKM(Z, 2, 3, 20);
U3'*Ufkm
pf2FK=psF(Yfkm,Ufkm)
pf3FK=psF(Yfkm3,Ufkm3)
[~, cluster_idx] = max(Ufkm, [], 2);
figure;
hold on;
colors = {'r', 'b', 'g'};
markers = {'o', 's', 'd'};
for k = 1:max(cluster_idx)
    idx = cluster_idx == k;
    scatter(Yfkm(idx, 1), Yfkm(idx, 2), 50, colors{k}, markers{k}, 'filled', 'DisplayName',
    text(Yfkm(idx, 1), Yfkm(idx, 2), num2str(k), 'HorizontalAlignment', 'center', 'VerticalA
end
xlabel('First Dimension');
ylabel('Second Dimension');
title('Classification of Objects on First Two Dimensions of Factorial K-means');
legend show;
grid on;
hold off;

%7. Compute Clustering and Disjoint PCA;
[Vcdpca,Ucdpca,Acdpca, Ycdpca,fcdpca,incdpca]=CDPCA(Z, 2, 2, 20);
[Vcdpca3,Ucdpca3,Acdpca3, Ycdpca3,fcdpca3,incdpca3]=CDPCA(Z, 2, 3, 20);
pf2CDP=psF(Ycdpca,Ucdpca)
pf3CDP=psF(Ycdpca3,Ucdpa3)
[Vdpca,Adpca, Ydpca, fdpca, indpca] = DPCA(Z,2,20)

% Create a heatmap of the loadings matrix
figure;
h = heatmap(Acdpca);
h.Title = 'Loadings Matrix';
h.XLabel = 'Components';
h.YLabel = 'Variables';
h.Colormap = winter;
h.ColorLimits = [-1, 1];
%%%%%%%%
explainedVarianceDim1 = 62.19;
explainedVarianceDim2 = 14.75;
Dim1 = Ycdpca(:, 1);
Dim2 = Ycdpca(:, 2);
figure;
hold on;
scatter(Dim1(group1), Dim2(group1), 50, 'b', 'o', 'filled');
scatter(Dim1(group2), Dim2(group2), 50, 'r', '^', 'filled');
ellipse(mean(Dim1(group1)), mean(Dim2(group1)), std(Dim1(group1)), std(Dim2(group1)), 'b');
ellipse(mean(Dim1(group2)), mean(Dim2(group2)), std(Dim1(group2)), std(Dim2(group2)), 'r');
xlabel(['Dim 1 (' num2str(explainedVarianceDim1) '% Variance)']);
ylabel(['Dim 2 (' num2str(explainedVarianceDim2) '% Variance)']);
title('Clustering and Disjoint PCA');
legend('Cluster 1', 'Cluster 2');
grid on;
axis equal;
hold off;

%8. Compute the Double K-means.
[Vdkm,Udkm,Ymdkm, fdkm,indkm]=DKM(Z, 2, 2, 20);
[Vdkm3,Udkm3,Ymdkm3, fdkm3,indkm3]=DKM(Z, 2, 3, 20);

%9. Compare the results by using the confusion matrix (contingency table);
%confusion matrix between 4 and 5
U2'*Urkm
%confusion matrix between 4 and 6
U2'*Ufkm
%confusion matrix between 4 and 7
U2'*Ucdpca
