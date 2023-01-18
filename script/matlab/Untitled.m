clear
clc
close all
%% cambiare in base al computer in uso

cd 'C:\Users\letcle\Documents\WORK\sfPCA-for-Bipolar-Disorder-and-Schizophrenia\script\matlab'
path_file = 'D:\ds000030_fmriprep'; %aggiornare indirizzo LaCie in base a computer

%% set paths
path_atlas = fullfile('..', '..','data', 'Hammers_mith_atlas_n30r83_delivery_Dec16');
path_data = fullfile('..', 'data');

path_roi = fullfile(path_atlas, 'n30r83_names_Excel.xlsx');
path_maps = fullfile('..', '..', 'maps');
path_participants = fullfile('..', '..', 'participants');

%% import relevant data

mass = spconvert(readmatrix(fullfile('..', '..', 'data', 'mass.txt')));
stiff = spconvert(readmatrix(fullfile('..', '..', 'data', 'stiff.txt')));

%% list maps available

maps_avail = dir(path_maps);
[~,region_names]  = xlsread(path_roi,'A3:A85');
%selected_rois = region_names(3:4); %si potrebbe prendere da maps avail e migliorare

% Define parameters sequence
Kfolds = 5; 
niter = 20; 
loglambdaseq = -3; 
N_PC = 5; 
%%
selected_rois = erase({maps_avail(3:size(maps_avail,1)).name}, '.csv');
selected_rois = erase(selected_rois, 'z_maps_')

%% meglio iterare su ROI per univocità
for r=3:size(selected_rois,2)
    tic
    fprintf(strcat('\n%%%%%%%%', selected_rois{r}, '%%%%%%%%\n'))
    %import map
    map = readmatrix(fullfile(path_maps,['z_maps_', selected_rois{r}, '.csv'])); % problema permessi matlab soliti
    participants = readtable(fullfile(path_participants,['part_ofinterest_', selected_rois{r}, '.csv']));
    
    % drop null rows 
    M = mean(map, 2, 'omitnan');
    X_lacking = find(M == 0);
    
    part_lacking = participants(X_lacking,:);
    writetable(part_lacking, fullfile(path_participants, strcat('part_lacking_', selected_rois{r}, '.csv')));

    map(X_lacking,:) = [];

    mean_map = mean(map);
    map_0 = map - ones(size(map,1),1)*mean_map;
    
    % Computation of the first N_PC PC function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Computation of first n PC %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    index_NA = isnan(map_0(1,:)); 
    
    F_not_normalized = zeros(N_PC, size(map_0,2)); 
    U_normalized = zeros(size(map_0,1), N_PC); 

    optimal_lambdas_indices = zeros(N_PC,1); 
    CV_tokeep = zeros(N_PC,size(loglambdaseq,2));

    X_residuals = map_0;
    for k = 1:N_PC
        disp(['Computing PC function' num2str(k) '...'])
        CVseq = fPCAManifold_Kfold(X_residuals, mass, stiff, loglambdaseq, Kfolds, niter, index_NA);
        CV_tokeep(k,:) = CVseq;
        [CV_min,CV_min_index] = min(CVseq);
        disp(['    Index chosen: ' num2str(CV_min_index) '.'])
        optimal_lambdas_indices(k) = CV_min_index;
        [F_not_normalized(k,:),U_normalized(:,k)] = fPCAManifold(X_residuals, mass, stiff, loglambdaseq(CV_min_index), niter, index_NA);
        %remove from data matrix PC computed, so we can repply the algorithm to
        %the residuals and get the next PC function
        X_residuals = X_residuals - (U_normalized(:,k) * F_not_normalized(k,:));
    end
    
    % Computation scores
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Computation scores %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    F_normalized = zeros(size(F_not_normalized));
    for k = 1:N_PC
        F_normalized(k,:) = F_not_normalized(k,:)./getL2Norm(F_not_normalized(k,:),mass);
    end

    %Un-normlized PC vectors (put back variance to each PC vector)
    U_not_normalized = zeros(size(U_normalized));
    for k = 1:N_PC
        U_not_normalized(:,k) = U_normalized(:,k).*getL2Norm(F_not_normalized(k,:),mass);
    end
    
    [Q,R] = qr(U_not_normalized);
    explained_var = diag(R).^2/(size(map_0,1));
    trace(U_not_normalized' * U_not_normalized)
    fig = plot(explained_var,'b--o');
    saveas(fig, fullfile('..', '..', 'plots', 'cumsum_variance', ['cumsumvariance_', selected_rois{r}, '.pdf']))
    
    X_0_clean = map_0;
    X_0_clean(isnan(map_0))=0;
    [U,S,V] = svd(X_0_clean,'econ');
    U_ALL = zeros(size(U_not_normalized));
    for i = 1:size(S,1)
        U_ALL(:,i) = U(:,i)*S(i,i)*getL2Norm(V(:,i),mass);
    end
    total_var = trace(U_ALL' * U_ALL)/(size(X_0_clean,1));

    
    writematrix(F_normalized, fullfile('..','..', 'results', [selected_rois{r}, '_F_normalized.csv']));
    writematrix(U_normalized, fullfile('..','..', 'results', [selected_rois{r}, '_U_normalized.csv']));
    writematrix(F_not_normalized, fullfile('..','..', 'results', [selected_rois{r}, '_F_not_normalized.csv']));
    writematrix(U_not_normalized, fullfile('..','..', 'results', [selected_rois{r}, '_U_not_normalized.csv']));
    writematrix(CV_tokeep, fullfile('..','..', 'results', [selected_rois{r}, '_CV.csv']));
    writematrix(cumsum(explained_var)./total_var, fullfile('..','..', 'results', [selected_rois{r}, '_exp_var.csv']));
    

    toc
end














