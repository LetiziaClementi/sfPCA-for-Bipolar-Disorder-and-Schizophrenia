clear all
close all
clc

%%
cd '~/OneDrive - Politecnico di Milano/PhD.LAVORO/3zo_anno/sfPCA_BIP_SCHZ/script/matlab'
addpath(genpath('FELICITY_Ver1.3.1'))

%%
output_path = fullfile('..', '..', 'results', 'sfPCA')
maps_path = fullfile('..', '..', 'maps');

%%

% z_maps = readmatrix(data_path);
X = z_maps;

% questp solo perché in script di eardi
N = isnan(X);
S = sum(isnan(X));
length(unique(S)) % capire bene come fare

%%

meanX = mean(X); %la media spalma i NaN su tutta la colonna (questo aiuta eliminazione)
X_nomean = X - ones(size(X,1),1)*meanX;
nan_nodes = find(isnan(X_nomean(1,:)));
% disp('Dataset Size ');
% size(X)

%% create manifold
manifold = [];
manifold.faces =  elem(:,1:4); %inversione nome ma sono sempre gli elementi
manifold.vertices = nodes;

%% FARE FEM!!!!
[R0,R1] = computeFEM(manifold); %questa  versione è specifica per superfici 2d immerse in uno spazio 3d
%%
mass = R0;
stiff = R1;

% mass = FEM(1).MAT;
% stiff = FEM(2).MAT;

%% Define parameters sequence
Kfolds = 5;
niter = 10;
loglambdaseq = 0:0.25:5;

%% importo i risultati fPCA da R

write_vtk_el(nodes, elem(:,1:4), fpcakfoldcoeff', '../plots/fpcakfoldcoeff.vtk', 'PC')
write_vtk_el(nodes, elem(:,1:4), fpcalambda1coeff', '../plots/fpcalambda1coeff.vtk', 'PC')





