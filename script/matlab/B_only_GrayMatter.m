clear all
close all
clc

%% Per visualizzare immagini settare a 1
viz_fig = 1;

%%
cd 'C:\Users\letcle\Documents\WORK\UCLA_WORK\script'

%%
path_atlas = strcat('../data/Hammers_mith_atlas_n30r83_delivery_Dec16/');
path_template = strcat('../data/mni_icbm152_nlin_asym_09c_nifti/mni_icbm152_nlin_asym_09c');
path_data = '/Volumes/LaCie/ds000030_fmriprep/'

%%
atlas_hammers = niftiread(strcat(path_atlas,'/Hammers_mith_atlas_n30r83_SPM5.nii'));
%%
labels = readmatrix('../data/labels_refined.csv');
elem = readmatrix('../data/elem_refined.csv');
nodes = readmatrix('../data/node_refined_spaziodati.csv'); %questi sono nello spazio DEI DATI
node_transf = nodes;
node_atlas = readmatrix('../data/node_atlas_refined.csv');
%%
% Tissue ID for the outputs are as follow:
%   0-Air/background, 1-Scalp, 2-Skull, 3-CSF, 4-GM, 5-WM, 6-air pockets

elem_csf = elem(elem(:,5)==3,:);
elem_gm = elem(elem(:,5)==4,:);
elem_wm = elem(elem(:,5)==5,:);

if (viz_fig == 1)
    figure, title('CSF'), plotmesh(nodes,elem_csf, 'x<100| y<100 | z >100')
    figure, title('GM'), plotmesh(nodes,elem_gm, 'x<100| y<100 | z >100') %mi sembra ok
    figure, title('WM'), plotmesh(nodes,elem_wm, 'x<100| y<100 | z >100')
end
%%

elem_new = vertcat(elem_gm(:,1), elem_gm(:,2), elem_gm(:,3));
elem_new = sort(elem_new);
elem_new = unique(elem_new);
node_gm = node_transf(elem_new, :);

%%
elem_unrolled = vertcat(elem(:,1), elem(:,2), elem(:,3), elem(:,4));

gm_unrolled = vertcat(elem_gm(:,1), elem_gm(:,2), elem_gm(:,3), elem_gm(:,4));
gm_unique = sort(unique(gm_unrolled));

notgm_unrolled = elem_unrolled(~ismember(elem_unrolled, gm_unique));
notgm_unique = sort(unique(notgm_unrolled));

node_onlygm = node_transf(gm_unique,:);
elem_onlygm = elem_gm;

%%
for j=1:size(elem_onlygm,1)
    for k = 1:4
        elem_onlygm(j,k)=elem_onlygm(j,k) - length(find(notgm_unique <= elem_onlygm(j,k)));
    end
end

%%
plotmesh(node_onlygm,elem_onlygm,'x>60 | z<30','FaceColor',[0.5 0.5 0.5],'EdgeAlpha',0.6)
% in ../script_intermedi/preprocess_cerebra.m c'è anche estrazione solo wm

%% labels only GM

size(node_onlygm)
node_atlas_onlygm = node_atlas(gm_unique,:);
size(node_atlas_onlygm)

%%
for i=1:size(node_atlas_onlygm,1)
       labels_onlygm(i)=atlas_hammers(round(node_atlas_onlygm(i,1)),round(node_atlas_onlygm(i,2)),round(node_atlas_onlygm(i,3)));
end

%% save
% data viz inspection emerge che label assegnate ok, gm isolata
% correttamente
write_vtk_el(node_atlas_onlygm, elem_onlygm(:,1:4), labels_onlygm', '../results_preproc/2201_13_refined/onlygm/labels_onlygm.vtk')

plotmesh(node_atlas_onlygm,elem_onlygm,'x>60 | z<30','FaceColor',[1 1 1],'EdgeAlpha',0.6)

%%
writematrix(labels_onlygm, '../results_preproc_1701/labels_onlygm.csv');
writematrix(elem_onlygm, '../results_preproc_1701/elem_onlygm.csv');
writematrix(node_onlygm, '../results_preproc_1701/node_spaziodati_onlygm.csv');
writematrix(node_atlas_onlygm, '../results_preproc_1701/node_atlas_onlygm.csv');

