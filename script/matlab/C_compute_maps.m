%% script to create multiple maps, all for the BIPOLAR disorder data 
% will then create a script to compare and evaluate which are significant


%%
clear all
close all
clc

%%
cd '~/OneDrive - Politecnico di Milano/PhD.LAVORO/3zo_anno/sfPCA_BIP_SCHZ/script/matlab'
path_file = 'D:\ds000030_fmriprep'; %aggiornare indirizzo LaCie in base a computer

%%
path_atlas = fullfile('..', '..','data', 'Hammers_mith_atlas_n30r83_delivery_Dec16');
path_data = fullfile('..', 'data');

path_roi = fullfile(path_atlas, 'n30r83_names_Excel.xlsx');
path_maps = fullfile('..','..', 'maps');

directorycontent = dir(path_file);

%% import mesh
nodes = readmatrix(fullfile('..', '..', 'data','mesh', 'node_spaziodati_onlygm.csv'));
elem = readmatrix(fullfile('..', '..', 'data','mesh', 'elem_onlygm.csv')); %la 5ta colonna in realtà nonn serve, sono solo 4
labels = (readmatrix(fullfile('..', '..', 'data','mesh', 'labels_onlygm.csv')))';
participants = readtable(fullfile('..', '..', 'data','participants.csv'));
%%
[~,region_names]  = xlsread(path_roi,'A3:A85');

d = dir(path_file);
dfolders = d([d(:).isdir]) ;
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..', 'phenotype'}))

% filtering parameters
TR = 2; %seconds, according to dataset documentation
Fs = 1/TR;
F_cut = 0.01;
band = [F_cut Fs];


%% temporary regions selected

for r=61:63
    
    ROI_name = region_names{r}
    ROI_id = find(strcmp(region_names,'PL_postce_G_l_sum' ));
    ROI_nodes = find(labels == ROI_id); %vediamo
    
    fprintf(strcat('\n%%%%%%%%', ROI_name, '%%%%%%%%\n'))
    
    desinenza_task = '_task-rest_bold_space-MNI152NLin2009cAsym_preproc.nii.gz';
    part_ofinterest = participants( (strcmp(participants.diagnosis, 'CONTROL')|strcmp(participants.diagnosis, 'SCHZ')|strcmp(participants.diagnosis, 'BIPOLAR'))&(participants.rest == 1),:);
    
    
    r_maps_task = zeros(size(part_ofinterest,1), size(nodes,1));
    nsbj = size(part_ofinterest,1);
    for i=1:nsbj
        name_task = strcat(part_ofinterest.participant_id(i),desinenza_task);
        
        try
            task = niftiread(char(strcat(path_file,'\', part_ofinterest.participant_id(i),'\func\', name_task)));
            
            task_unrolled = zeros(size(nodes,1), size(task,4));
            for j=1:size(nodes,1)
                task_unrolled(j,:) = task(round(nodes(j,1)), round(nodes(j,2)), round(nodes(j,3)),:);
            end
            %fprintf('succesfully unrolled\n');
            
            nsample_task = size(task,4);
            t_task = (0:nsample_task-1)*TR;
            % task_filtered = zeros(size(task_unrolled));
            % filtraggio
            %task_filtered = filterTS(task_unrolled, t_task, band); %funz mod per non usare t
            task_filtered = task_unrolled;
            
            %fprintf('succesfully filtered\n');
            
            % serie media della ROI per BART
            ROI_ts = task_filtered(ROI_nodes,:);
            mean_ROI = mean(ROI_ts);
            
            for j = 1:size(nodes,1)
                r_maps_task(i,j)=corr(mean_ROI',task_filtered(j,:)');
            end
            
        catch
            warning('file did not exist or something')
        end
        
        % print per verifica
        fprintf("completed: %s, n %d/%d \n", char(part_ofinterest.participant_id(i)), i, nsbj)
    end
    %TIME : ~30 sec each
    
    z_maps_task = 0.5 * (log(1+r_maps_task) - log(1-r_maps_task));
    writematrix(z_maps_task, fullfile(path_maps, strcat('z_maps_', ROI_name, '.csv')));
    writetable(part_ofinterest, fullfile(path_maps, strcat('part_ofinterest_', ROI_name, '.csv')));
end