clear all
close all
clc

%% To visualize images set to 1
viz_fig = 1;

%%
cd 'C:\Users\letcle\Documents\WORK\UCLA_WORK\script'
%%
root_path = '/Users/letizia/OneDrive - Politecnico di Milano/PhD.LAVORO/2ndo_anno'

%%
path_atlas = strcat(root_path, '../data/Hammers_mith_atlas_n30r83_delivery_Dec16/');
path_template = strcat(root_path, '../data/mni_icbm152_nlin_asym_09c_nifti/mni_icbm152_nlin_asym_09c');
path_data = '/Volumes/LaCie/ds000030_fmriprep/'

%% import template and atlas
atlas_hammers = niftiread('../Hammers_mith_atlas_n30r83_delivery_Dec16/Hammers_mith_atlas_n30r83_SPM5.nii');
template_mni9c_csf = niftiread(strcat(path_template,'/mni_icbm152_csf_tal_nlin_asym_09c.nii'));
template_mni9c_gm = niftiread(strcat(path_template,'/mni_icbm152_gm_tal_nlin_asym_09c.nii'));
template_mni9c_wm = niftiread(strcat(path_template,'/mni_icbm152_wm_tal_nlin_asym_09c.nii'));
template_mni9c_t1 = niftiread(strcat(path_template,'/mni_icbm152_t1_tal_nlin_asym_09c.nii'));

%% visualizee template (per tissue)

if (viz_fig == 1)
    figure,
    for i=1:size(template_mni9c_csf,3)
        subplot(1,3,1)
        imagesc(template_mni9c_csf(:,:,i)), title('csf')
        subplot(1,3,2)
        imagesc(template_mni9c_gm(:,:,i)), title('gm')
        subplot(1,3,3)
        imagesc(template_mni9c_wm(:,:,i)), title('wm')
        pause(0.1)
    end
end
%% creazione struttura adeguata per brain2mesh
mni_tal_9c.wm=template_mni9c_wm;
mni_tal_9c.gm=template_mni9c_gm;
mni_tal_9c.csf=template_mni9c_csf;

%% impostazioni per mesh

cfg.radbound.wm=1.7;
cfg.radbound.gm=1.7;
cfg.radbound.csf=1.7;

cfg.maxnode = 13000;
cfg.maxvol = 700;

%% creazione mesh
[node,elem,face] = brain2mesh(mni_tal_9c,cfg);

%%
%writematrix(node, '../data/mesh/node_original.csv')
%writematrix(elem, '../data/mesh/elem_original.csv')

%%
% viz mesh
if(viz_fig == 1)
    figure,
    plotmesh(node,elem(elem(:,5)==5,:),'FaceColor',[1 1 1],'EdgeAlpha',0.6) %%wm
    hold on;
    plotmesh(node,elem(elem(:,5)==4,:),'x>40| y<20','FaceColor',[0.35 0.35 0.35],'EdgeAlpha',0.6) %%pial
    plotmesh(node,elem(elem(:,5)==3,:),'x>60 | z<30','FaceColor',[0.2 0.6 1],'EdgeAlpha',0.6) %%csf
end


%% import data for alignment (BART task)

d = dir(path_data);
dfolders = d([d(:).isdir]) ;
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..', 'phenotype'}))
desinenza_bart = '_task-bart_bold_space-MNI152NLin2009cAsym_preproc.nii.gz';

i = 1;
name_bart = [dfolders(i).name desinenza_bart]
bart = niftiread([path_data dfolders(i).name '/func/' name_bart]);

img_ref_bart = bart(:,:,:,1);
if(viz_fig == 1)
    figure,
    imagesc(img_ref_bart(:,:,10))
end

%% if rerun

bart = niftiread(['../sub-10159' '/func/' 'sub-10159_task-bart_bold_space-MNI152NLin2009cAsym_preproc.nii.gz']);
img_ref_bart = bart(:,:,:,1);

%% nodi are transfered in data's space

K = size(template_mni9c_csf)./size(img_ref_bart)
node_transf(:, 1) = node(:, 1)./K(1);
node_transf(:, 2) = node(:, 2)./K(2);
node_transf(:, 3) = node(:, 3)./K(3); 

%%
if(viz_fig == 1)
    figure,
    plotmesh(node_transf,elem(elem(:,5)==5,:),'FaceColor',[1 1 1],'EdgeAlpha',0.6) %%wm
    hold on;
    plotmesh(node_transf,elem(elem(:,5)==4,:),'x>40| y<20','FaceColor',[0.35 0.35 0.35],'EdgeAlpha',0.6) %%pial
    plotmesh(node_transf,elem(elem(:,5)==3,:),'x>60 | z<30','FaceColor',[0.2 0.6 1],'EdgeAlpha',0.6) %%csf
end
%% naif attempt

K = size(atlas_hammers)./size(template_mni9c_csf)

% proporzione tra template e atlante (nodi atlas in spazio atlante)
% node_atlas(:, 1) = node(:, 1)./K(1);
% node_atlas(:, 2) = node(:, 2)./K(2);
% node_atlas(:, 3) = node(:, 3)./K(3); 

%% more precise 
% i nodi originali sono a partire dal template
% per i dati --> template/dati
% per atlante --> template/atlante (non serve fare un passaggio in +)

node_atlas = zeros(size(node));

node_atlas(:, 1) = 18 + ((162 - 18)/(172 - 23))*(node(:, 1) - 23);
node_atlas(:, 2) = 23 + ((198 - 23)/(210 - 23))*(node(:, 2) - 23);
node_atlas(:, 3) = 7 + ((152 - 7)/(164 - 6))*(node(:, 3) - 6); 

% ATLAS
% Y (dim1):
% prima riga: 18
% ultima: 162
% X (dim2):
% prima riga: 23
% ultima: 198
% Z (dim3):
% prima vera slice: la 7
% ultima: 152


% TEMPLATE
% Y (dim1):
% prima riga: 23
% ultima: 172
% X (dim2):
% prima riga: 23
% ultima: 210
% Z (dim3):
% prima vera slice: la 6
% ultima: 164


%% assign a label to each node
for i=1:size(node_atlas,1)
       labels_hammers(i)=atlas_hammers(round(node_atlas(i,1)),round(node_atlas(i,2)),round(node_atlas(i,3)));
end


% % per viz:
% write_vtk_el(node_atlas, elem(:,1:4), labels_hammers, '../data/mesh/labels_hammers2.vtk')
% write_vtk_el(node_transf, elem(:,1:4), labels_hammers', '../results_preproc_1701/labels_dati.vtk')

%%
writematrix(labels_hammers, '../data/mesh/labels_refined.csv');
writematrix(elem, '../data/mesh/elem_refined.csv');
writematrix(node_transf, '../results_preproc_1701/node_refined_spaziodati.csv');
%%
writematrix(node_atlas,  '../data/mesh/node_atlas_refined.csv');