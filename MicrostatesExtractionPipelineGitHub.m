%% SCALP MICROSTATE AND INVERSE SOLUTION

% Author: Sara Goldoni, Giulia Demartis
% Collaborators:Edoardo Paolini
%% Clear
clear
close all
clc

%% Add FieldTrip to the Path
addpath  % addpath fieldtrip
ft_defaults;

%% Loading Data
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
% Specify the path and the header of the BrainVision files
D = %directory were data are;
hdr = %name of the data file;
% load Brain Vision Data Exchange format dataset and return EEGLAB EEG structure
EEG = pop_loadbv(D, hdr);


%% Cut of 'bad' epochs 
%we had to remove bad epochs, skip this step if not needed
times = [1:10238, 10752:22655, 23019:23425, 23764:115270, 116486:161127];
EEG.data = EEG.data(:, times);
EEG.times = EEG.times(:, times);
% update channels and data
load('new_chanlocs_data.mat'); %spike montage
EEG.chanlocs = new_chanlocs;
EEG.data = EEG.data(new_data,:);
EEG.nbchan = length(EEG.chanlocs);
% Transpose it to FieldTrip
eeg_data = eeglab2fieldtrip(EEG, 'preprocessing');
ti = 1;
tf = min(158698, length(EEG.times)); % edit the first number
eeg_data.trial{1,1} = eeg_data.trial{1,1}(:,ti:tf);
eeg_data.time{1,1} = eeg_data.time{1,1}(:,ti:tf);

%% Create Layout for Scalp microstates using Fieldtrip
cfg =[];
cfg.elec = eeg_data.elec; 
cfg.rotate = 0;
cfg.outline = 'convex';
cfg.skipcomnt   = 'yes';
cfg.skipscale = 'yes';
layout_scalp = ft_prepare_layout(cfg, eeg_data); 

% Plot Layout
figure('Name', 'Layout');
ft_plot_layout(layout);

%% Scalp Microstates

% Initialize an empty microstate object
ms_scalp = microstate.individual;

% Export eeg_data to a +microstate object
ms_scalp = ms_scalp.import_fieldtrip(eeg_data);

% Estimate maps
ms_scalp = ms_scalp.cluster_estimatemaps(5);

% Plot Maps
figure('Name', 'Scalp Microstates');
ms_scalp.plot('maps',layout_scalp);

%% Cutting time for inverse solution
% Cut the time to run the program for inverse solution 
ti = 1;
tf = min(2000,length(EEG.times)); % edit the first number to suit your RAM capabilities
eeg_data.trial{1,1} = eeg_data.trial{1,1}(:,ti:tf);
eeg_data.time{1,1} = eeg_data.time{1,1}(:,ti:tf);

%% Loading HeadModel and SourceModel

% Headmodel: geometrical and conductive properties of the head
% Solution Points: location of the potential sources

% Load the Headmodel and Solution Points
load('headmodel.mat');
solution_points = ft_read_headshape('cortex_8196.surf.gii');
%% Create Layout for inverse solution

% Loads or creates a 2-D layout of the channel locations. 
% This layout is required for plotting the topographical distribution of the potential or
% field distribution, or for plotting timecourses in a topographical arrangement.
cfg =[];
cfg.elec = eeg_data.elec;  % info about electrodes
cfg.rotate = 0 ;
cfg.outline= 'convex';
layout = ft_prepare_layout(cfg, eeg_data);   
%% Electrodes alignment with Headmodel and Plotting

% Electrode alignment with the headmodel
cfg = [];
cfg.method = 'project'; % projects electrodes onto scalp surface
cfg.headshape = headmodel.bnd(1);
eeg_data.elec = ft_electroderealign(cfg, eeg_data.elec);

% Adjusting the units (mm -> cm)
headmodel = ft_convert_units(headmodel,eeg_data.elec.unit);
solution_points = ft_convert_units(solution_points,eeg_data.elec.unit);

% Plot the Headmodel
% Visualise the surface of the three layers of the headmodel
figure('Name', 'HeadModel');
ft_plot_mesh(headmodel.bnd(3),'facecolor','r','edgecolor','none');
ft_plot_mesh(headmodel.bnd(2),'facecolor','g','edgecolor','none');
ft_plot_mesh(headmodel.bnd(1),'facecolor','b','edgecolor','none');
view([45 0 0]); alpha 0.3; rotate3d;
% Plot the electrodes on the scalp surface
figure('Name', 'Electrodes on the Scalp Surface');
hold on
ft_plot_sens(eeg_data.elec, 'elecsize', 40, 'label', 'label');
ft_plot_headmodel(headmodel, 'facealpha', 0.5);
view(90, 0)

% Plot the Solution Points
% Visualise the surface of the sourcemodel
figure('Name', 'Surface of the sourcemodel');
ft_plot_mesh(solution_points,'facecolor','brain','edgecolor','none')
camlight('infinite'); lighting gouraud; rotate3d;
% Visualise the solution points of the sourcemodel
figure('Name', 'Solution Points');
plot3(solution_points.pos(:,1), solution_points.pos(:,2), solution_points.pos(:,3), 'r.')
axis equal; axis vis3d;
grid on;

%% Calculation of Timelocked Signal and GMFP, and Plotting

% This function compute the timelocked (over trials) average potentials.
% In our case, having only one trial, we would not need it, however, 
% we modify the structure of the variable so that we can perform the subsequent analyses

cfg = [];
cfg.covariance = 'yes';
eeg_data = ft_timelockanalysis(cfg, eeg_data);

% TEP Plotting 

% 1. Plot the EEG projected on each electrode

figure('Name', 'Multiplot ER')
cfg = [];
cfg.fontsize = 6;
cfg.layout = layout;
cfg.baseline = 'yes';
cfg.showlabels = 'yes';
ft_multiplotER(cfg, eeg_data);

% 2. Plot a topographic map of an EEG field as a 2-D circular view

figure('Name', 'Topoplot ER');
cfg = [];
cfg.xlim = [0:10:100];
cfg.zlim = [-20 20];
cfg.layout = layout;
cfg.colorbar = 'yes';
cfg.fontsize = 3;
ft_topoplotER(cfg, eeg_data);

% Global Mean Field Power (GMFP) Calculation
cfg = [];
cfg.method = 'amplitude';
eeg_gmfp = ft_globalmeanfield(cfg, eeg_data);

% 3. Plot of all the potentials, including the GMFP
figure('Name', 'Potentials and GMFP Plot');
plot(eeg_data.time, eeg_data.avg');
hold on;
plot(eeg_gmfp.time, eeg_gmfp.avg', 'LineWidth', 3, 'Color', 'k');
xlabel('Time (ms)');
ylabel('Amplitude (mV)');


%% LeadField Calculation
% Computes the forward model for many dipole locations on a sourcemodel
% and stores it for efficient inverse modelling
% leadfield = ft_prepare_leadfield(cfg);

cfg = [];
cfg.elec = eeg_data.elec;            
cfg.channel = eeg_data.label;  
cfg.headmodel = headmodel; 
cfg.grid = solution_points;
leadfield = ft_prepare_leadfield(cfg);

%% INVERSE SOLUTION
%% Inverse Solution computation using eLORETA

% Perform the inverse solution 
cfg = [];
cfg.method = 'eloreta';       
cfg.sourcemodel = leadfield;       
cfg.headmodel = headmodel;       
cfg.eloreta.lambda = 0.05;
cfg.eloreta.keepmom = 'yes'; 
cfg.eloreta.keepfilter = 'yes'; 
% mantenere il campo pow
eeg_source = ft_sourceanalysis(cfg, eeg_data);

% Computes descriptive parameters of the source analysis results.
% The interest is on 'mom' field which contains the time courses of 
% the event-related field at the source level. Colloquially, these 
% time courses are known as ‘virtual channels’, reflecting the signal
% that would be picked up if it could directly be recorded by a channel 
% at that location. Correnti in sorgente.
cfg = [];
cfg.projectmom = 'yes';
cfg.keeppow = 'yes';
eeg_source = ft_sourcedescriptives(cfg,eeg_source);

% Rearrange mom storage
mom = zeros(length(eeg_source.avg.mom), length(eeg_source.time));
for i=1:size(mom,1)
    if isempty(eeg_source.avg.mom{1,i})
        continue;
    end
    mom(i,:) = eeg_source.avg.mom{1,i};
end
eeg_source.avg.mom = mom;

% Displays the source reconstruction on a cortical mesh
cfg = [];
cfg.funparameter = 'mom';
ft_sourcemovie(cfg,eeg_source);

%% Source Interpolation with MRI Template

% Load MRI template
mri = ft_read_mri('MNI.nii');

% Interpolation of the sources with the MRI template 
cfg = [];
cfg.keeppow='yes';
cfg.parameter = 'mom';
source_int = ft_sourceinterpolate(cfg,eeg_source,mri);

cfg = [];
cfg.parameter = 'pow';
source_int_pow = ft_sourceinterpolate(cfg,eeg_source,mri);
%% Parcelization
% The output is a channel-based representation with the combined
% (e.g.averaged) representation of the source parameters per parcel.
% A region or cluster of voxels or vertices in the brain that is treated 
% as a single unit for analysis. 
% Parcellation involves dividing the brain into distinct regions, 
% often based on anatomical or functional criteria.

% Load the atlas for the parcellization
atlas = ft_read_atlas('schaefer.nii'); 
load('schaefer.mat');
atlas.tissuelabel=tissuelabel';
% Convert units
atlas = ft_convert_units(atlas, source_int.unit);
% Combines the source-reconstruction parameters over the parcels
cfg = [];
source_parcel = ft_sourceparcellate(cfg, source_int, atlas);

cfg = [];
source_parcel_pow = ft_sourceparcellate(cfg, source_int_pow, atlas);

%% Source microstate extraction

% Conversion into +microstate format
ms = microstate.individual();
ms= ms.import_fieldtrip(source_parcel,'modality','source');

% Assigning missing values
ms.data = source_parcel.mom'; 

% Microstates map estimation
ms= ms.cluster_estimatemaps(5);

% layout creation with +microstate
template= 'schaefer';
labels= atlas.tissuelabel;
layout_source=microstate.functions.layout_creator(template,labels);

% plotting source microstates
figure;
ms.plot('maps', layout_source, 'cscale', [0.5, 1]);


%%
for i=1:5
    source_parcel_pow.pow=ms.maps(:,i);
    figure('Name', strcat('Microstate:', i));
    cfg = [];
    cfg.method = 'ortho';
    cfg.funparameter = 'pow';
    cfg.funcolorlim= [0.1 0.37];
    cfg.funcolormap = 'hot';
    ft_sourceplot(cfg,source_parcel_pow,mri);
end


%% Plot EEG sources after parcellization
% Plot EEG sources after parcellization

figure('Name', 'EEG Sources after Parcellization');
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'mom';
cfg.funcolormap = 'jet';
ft_sourceplot(cfg,source_parcel);


