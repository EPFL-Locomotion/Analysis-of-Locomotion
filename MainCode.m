

%% Matrix with all the kinematic features
[FeaturesSCINoFloat]=SCIPartNoFloat(100);
[FeaturesSCIFloat]=SCIPartFloat(100);
[FeaturesHealthyFloat]=HealthyPartFloat(100);
[FeaturesHealthyNoFloat]=HealthyPartNoFloat(100);
NumbFeatures=(fieldnames(FeaturesSCINoFloat));

for i=1:size(fieldnames(FeaturesSCINoFloat),1)-6

    PCAKinMatrix(:,i)=[FeaturesSCINoFloat.(NumbFeatures{i});FeaturesSCIFloat.(NumbFeatures{i});FeaturesHealthyFloat.(NumbFeatures{i});FeaturesHealthyNoFloat.(NumbFeatures{i})];

end

% Pca on kinematic features and plot of the explain variance

[coeff,score,variance,~,explain]=pca(zscore(PCAKinMatrix));
figure;
bar(explain);
title('Principal Component explaning variances');
xlabel('Principal Components');
ylabel('Explain Variance');
grid on

%plot of the loading of the first principal component   
[sorted_coeff1P, sorting_value1P] = sort(abs(coeff(:,1)),'descend');
SortedFeatures = NumbFeatures(sorting_value1P);
figure
bar(sorted_coeff1P);
xticks(1:40);
xticklabels(SortedFeatures(1));
title({'Loadings of the 1st principal component (sorted)'});
xlabel('Features');


% rappresentation of coefficient of 1PC
coeff(:,1)=coeff(:,1)/max(coeff(:,1));%normalization of the loadings

figure;
imagesc(coeff(:,1));
colorbar;
caxis([-1 1]);
set(gca,'YDir','normal')
title('PC1');
set(gca,'TickLength',[0 0]);
set(gca,'xtick',[]);
for i=1:40
    hold on
    line([0.5 1.5], [i-0.5 i-0.5],'Color','k');
end   


% rappresentation of coefficient of 2PC

coeff(:,2)=coeff(:,2)/max(coeff(:,2));%normalization of the loadings

figure;
imagesc(coeff(:,2));
colorbar;
caxis([-1 1]);
set(gca,'YDir','normal')
title('PC2');
set(gca,'TickLength',[0 0]);
set(gca,'xtick',[]);
for i=1:40
    hold on
    line([0.5 1.5], [i-0.5 i-0.5],'Color','k');
end   



%rappresention in 3D of the kineamtic features
i=1;

Size1=size(FeaturesSCINoFloat.(NumbFeatures{i}),1);
Size2=size(FeaturesSCIFloat.(NumbFeatures{i}),1);
Size3=size(FeaturesHealthyFloat.(NumbFeatures{i}),1);
Size4=size(FeaturesHealthyNoFloat.(NumbFeatures{i}),1);

figure;
scatter3(score(1:Size1,1),score(1:Size1,2),score(1:Size1,3),'filled');
hold on
scatter3(score(Size1:Size1+Size2,1),score(Size1:Size1+Size2,2),score(Size1:Size1+Size2,3),'filled');
hold on
scatter3(score(Size1+Size2:Size1+Size2+Size3,1),score(Size1+Size2:Size1+Size2+Size3,2),score(Size1+Size2:Size1+Size2+Size3,3),'filled');
hold on
scatter3(score(Size1+Size2+Size3:Size1+Size2+Size3+Size4,1),score(Size1+Size2+Size3:Size1+Size2+Size3+Size4,2),score(Size1+Size2+Size3:Size1+Size2+Size3+Size4,3),'filled');
xlabel('1PC');
ylabel('2PC');
zlabel('3PC');
legend('SCINoFloat','SCIFloat','HealthyFloat','HealthyNoFloat');


% rappresenation in 2D of the kinematic features

figure;
scatter(score(1:Size1,1),score(1:Size1,2),'filled');
hold on
scatter(score(Size1:Size1+Size2,1),score(Size1:Size1+Size2,2),'filled');
hold on
scatter(score(Size1+Size2:Size1+Size2+Size3,1),score(Size1+Size2:Size1+Size2+Size3,2),'filled');
hold on
scatter(score(Size1+Size2+Size3:Size1+Size2+Size3+Size4,1),score(Size1+Size2+Size3:Size1+Size2+Size3+Size4,2),'filled');
xlabel('1PC');
ylabel('2PC');
legend('SCINoFloat','SCIFloat','HealthyFloat','HealthyNoFloat');

%% PCA on Kinematics and emg features


for i=1:size(fieldnames(FeaturesSCINoFloat),1)

    PCAEMGKinMatrix(:,i)=[FeaturesSCINoFloat.(NumbFeatures{i});FeaturesSCIFloat.(NumbFeatures{i});FeaturesHealthyFloat.(NumbFeatures{i});FeaturesHealthyNoFloat.(NumbFeatures{i})];

end

 

[coeff_EMG,score_EMG,variance_EMG,~,explainEMG]=pca(zscore(PCAEMGKinMatrix));
figure;
bar(explainEMG);
title('Principal Component explaning variances');
xlabel('Principal Components');
ylabel('Explain Variance');
grid on


%plot of the loading of the first principal component  
[sorted_coeff1PEMG, sorting_value1PEMG] = sort(abs(coeff_EMG(:,1)),'descend');
SortedFeaturesEMG = NumbFeatures(sorting_value1PEMG);
figure
bar(sorted_coeff1PEMG);
xticks(1:35);
xticklabels(SortedFeaturesEMG(1));
title({'Loadings of the 1st principal component (sorted)'});
xlabel('Features');

% rappresentation of coefficient of 1PC
coeff_EMG(:,1)=coeff_EMG(:,1)/max(coeff_EMG(:,1));%normalization of the loadings
figure;
imagesc(coeff_EMG(:,1));
cb=colorbar;
caxis([-1 1]);
set(gca,'YDir','normal')
title('PC1');
set(gca,'TickLength',[0 0]);
set(gca,'xtick',[]);
for i=1:46
    hold on
    line([0.5 1.5], [i-0.5 i-0.5],'Color','k');
end    

% rappresentation of coefficient of 2PC
coeff_EMG(:,2)=coeff_EMG(:,2)/max(coeff_EMG(:,2));%normalization of the loadings
figure;
imagesc(coeff_EMG(:,2));
colorbar;
caxis([-1 1]);
set(gca,'YDir','normal')
title('PC2');
set(gca,'TickLength',[0 0]);
set(gca,'xtick',[]);
for i=1:46
    hold on
    line([0.5 1.5], [i-0.5 i-0.5],'Color','k');
end  

% rappresenation in 3D of the kinematic features 
figure;
scatter3(score_EMG(1:Size1,1),score_EMG(1:Size1,2),score_EMG(1:Size1,3),'filled');
hold on
scatter3(score_EMG(Size1:Size1+Size2,1),score_EMG(Size1:Size1+Size2,2),score_EMG(Size1:Size1+Size2,3),'filled');
hold on
scatter3(score_EMG(Size1+Size2:Size1+Size2+Size3,1),score_EMG(Size1+Size2:Size1+Size2+Size3,2),score_EMG(Size1+Size2:Size1+Size2+Size3,3),'filled');
hold on
scatter3(score_EMG(Size1+Size2+Size3:Size1+Size2+Size3+Size4,1),score_EMG(Size1+Size2+Size3:Size1+Size2+Size3+Size4,2),score_EMG(Size1+Size2+Size3:Size1+Size2+Size3+Size4,3),'filled');
xlabel('1PC');
ylabel('2PC');
zlabel('3PC');
legend('SCINoFloat','SCIFloat','HealthyFloat','HealthyNoFloat');


% rappresenation in 2D of the kinematic features

figure;
scatter(score_EMG(1:Size1,1),score_EMG(1:Size1,2),'filled');
hold on
scatter(score(Size1:Size1+Size2,1),score_EMG(Size1:Size1+Size2,2),'filled');
hold on
scatter(score_EMG(Size1+Size2:Size1+Size2+Size3,1),score_EMG(Size1+Size2:Size1+Size2+Size3,2),'filled');
hold on
scatter(score_EMG(Size1+Size2+Size3:Size1+Size2+Size3+Size4,1),score(Size1+Size2+Size3:Size1+Size2+Size3+Size4,2),'filled');
xlabel('1PC');
ylabel('2PC');
legend('SCINoFloat','SCIFloat','HealthyFloat','HealthyNoFloat');