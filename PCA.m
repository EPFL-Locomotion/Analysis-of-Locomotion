%% PCA-Kinematics
%the total matrix is the following:

%Healthy_Float subject: matrix(9*29);
%Healthy_NO_Float subject: matrix(10*29);
%SCI_FLOAT subject: matrix(14*29);
%NO_SCI_FLOAT subject: matrix(21*29);


%% PCA-Kinematics (54*29)
[FeaturesSCINoFloat]=SCIPartNoFloat(100,1000);
[FeaturesSCIFloat]=SCIPartFloat(100,1000);
[FeaturesHealthyFloat]=HealthyPartFloat(100,1000);
[FeaturesHealthyNoFloat]=HealthyPartNoFloat(100,1000);
NumbFeatures=(fieldnames(FeaturesSCINoFloat));

for i=1:size(fieldnames(FeaturesSCINoFloat),1)-6

    PCAKinMatrix(:,i)=[FeaturesSCINoFloat.(NumbFeatures{i});FeaturesSCIFloat.(NumbFeatures{i});FeaturesHealthyFloat.(NumbFeatures{i});FeaturesHealthyNoFloat.(NumbFeatures{i})];

end


[coeff,score,variance]=pca(zscore(PCAKinMatrix));
[sorted_coeff1P, sorting_value1P] = sort(abs(coeff(:,1)),'descend');
SortedFeatures = NumbFeatures(sorting_value1P);
figure
bar(sorted_coeff1P);
xticks(1:29);
xticklabels(SortedFeatures(1));
title({'Loadings of the 1st principal component (sorted)'});
xlabel('Features');



%% PCA Kinematics plus emg (54*33)



for i=1:size(fieldnames(FeaturesSCINoFloat),1)

    PCAEMGKinMatrix(:,i)=[FeaturesSCINoFloat.(NumbFeatures{i});FeaturesSCIFloat.(NumbFeatures{i});FeaturesHealthyFloat.(NumbFeatures{i});FeaturesHealthyNoFloat.(NumbFeatures{i})];

end


[coeff_EMG,score_EMG,variance_EMG]=pca(zscore(PCAEMGKinMatrix));
[sorted_coeff1PEMG, sorting_value1PEMG] = sort(abs(coeff_EMG(:,1)),'descend');
SortedFeaturesEMG = NumbFeatures(sorting_value1PEMG);
figure
bar(sorted_coeff1PEMG);
xticks(1:29);
xticklabels(SortedFeaturesEMG(1));
title({'Loadings of the 1st principal component (sorted)'});
xlabel('Features');