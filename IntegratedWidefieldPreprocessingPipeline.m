%% widefield preprocessing pipeline
%integrating Yixiang's with Xinxin's
%updated 03/2022
clear all; close all; clc

fileInfo = readtext('D:\Test_data\scripts\widefield_preprocessing_inputs.txt', ' ');
%1 is path, 2 is file, 3 is
%roiSet, 4 is spont or not (spont = 1), 5 is rig

%spatial factor downsample
spatialFactor = 2; %downsample factor
%Whether discard frames that contain large motions, recommend setting to 0
moveAssessFlag = 0; 
%flag for to use as F for the dF/F calculation, 1 = bottom 5th and 2 = mean
%values
flag = 1;

%kk = 1;
for kk = 1:size(fileInfo, 1)
    datapath = fileInfo{kk, 1};
    cd(datapath)
    filename = fileInfo{kk, 2};
    outputFolder = strcat('output', '_', filename(1:end-4));
    mkdir(outputFolder);
    loadtag = filename(1:3); %name of reference file for motion correct
    spont = fileInfo{kk, 4};
    rig = fileInfo{kk, 5};
    %% Read in movie
    %This function takes advantage of the openMovie function
    %from the wholeBrainDX-master and uses it to input frames in
    %a specific movie. It allows concatenation of several movies
    %with a common prefix. From Yixiang's script.
    A = openMovie(filename);
    tmp = dir([filename(1:length(filename)-4) '@00*.tif']);
    %Make a combined matrix for each recording
    if ~isempty(tmp)
        for j = 1:numel(tmp)
            fn = tmp(j).name;
            B = openMovie(fn);
            A = cat(3, A, B);
            clear B
        end
    end
    
    %% Downsample movie
    A1 = imresize(A, 1/spatialFactor, 'bilinear');
    
    %% Motion correction
    % NormCore is default, dft is less sensitive
    tic;
    
    % try NormCore method first then go to dft if failed
    disp('Trying NoRMCorre algo to register movies')
    
    [A2, movTag, tform_all, NormTform_all, movIdx_saved] = ...
       movAssess_NoRMCorre(A1, moveAssessFlag, outputFolder, filename, spatialFactor);
                        
    if isempty(A2)
        %If failed, try the dft method
        disp('The NoRMCorre algo did not work')
        disp('Switch to fast discrete fourier transformation')
        [A2, movTag, tform_all, NormTform_all, movIdx_saved] = ...
            movAssess(A1, moveAssessFlag, outputFolder, filename, spatialFactor, loadtag);
    end

    disp('Movement assessment finished...Time cost = ')
    toc;
    
    clear A1
    
    %Save movemet assessment results
    checkname = [filename(1:length(filename)-4) '_moveAssess' movTag '.mat'];
    save(fullfile(outputFolder,checkname),'tform_all','NormTform_all','movIdx_saved');
    
    %% Apply ROI mask(s)
    % Apply 3D mask to matrix A. First call the ROIMask function to generate 2D
    % mask, then extend 2D mask to 3D mask.
    sz = size(A2);
    ROIname = fileInfo{kk, 3};
    ROIData = ReadImageJROI(ROIname);
    
    A3 = ROI.ApplyMask(A2, ROIData, spatialFactor, rig);
    % Get 2D ROI mask from ROIData. If there are more than 1 ROI, merge
    % the ROIs to a single mask. ROI is considered a polygon are
    % defined in Image J. Inputs: ROIData, size of matrix, frame size.
    % Outputs: indices of vertices of polygons, 2d mask containing all
    % ROIs
    
    disp('Successfully apply ROI')
    
    %% df/f outside of ROI
    % used for partial seedbased correlation after preprocessing
    % spillover fluorescence from nearby areas regressed out
    % anything outside ROI correlated with the activity within could be
    % regressed out - warning
    
    A2 = reshape(A2,[sz(1)*sz(2),sz(3)]);
    A3 = reshape(A3, [sz(1)*sz(2),sz(3)]);
    outInd = A3(:,1) == 0;
    A_out = A2(outInd, :);
    A_out(A_out == 0) = nan;
    Avg_out = nanmean(A_out,1); Avg_out = reshape(Avg_out, [1, 1, sz(3)]);
    
    sz2 = size(Avg_out);
    A_re = reshape(Avg_out, [sz2(1)*sz2(2),sz2(3)]);
    
    if flag == 1
        %Use bottom 5th percentile as F0
        A_F0 = prctile(A_re,5,2);
        A_F0 = repmat(A_F0,[1,sz2(3)]);
        disp('Compute dF/F0 based on 5th percentile')
        save('A_F0_5th.mat', 'A_F0');
    else
        % Use mean values as F0
        A_F0 = repmat(mean(A_re,2),[1,sz2(3)]);
        disp('Compute dF/F0 based on 50th percentile')
        save('A_F0_50th_background.mat', 'A_F0');
    end
    Avg_out_dFoF = reshape(A_re./A_F0 - 1,sz2);
    
    if any(isnan(Avg_out_dFoF))
        Avg_out_dFoF = [];
        disp('Background not defined! Check if .roi file is provided!')
    end
    checkname = [filename(1:length(filename)-4) '_out_dFoF.mat'];
    save(fullfile(outputFolder,checkname),'Avg_out_dFoF');
    clear A2 Avg_out_dFoF A_F0
    
    %% focus on ROI and downsample
    
    %Focusing on just the ROI part of the movie
    A3 = reshape(A3, sz);
    A4 = focusOnroi(A3);
    clear A3
    
    %Get the downsampled roi mask
    sz = size(A4);
    ds_Mask = repmat((A4(:,:,1) ~= 0),[1,1,sz(3)]);
    smallMask = ds_Mask(:,:,1);
    
    %Save mat file with data
    checkname = [filename(1:length(filename)-4) '_instance.mat'];
    A_save = reshape(A4, [sz(1)*sz(2),sz(3)]);
    save(fullfile(outputFolder,checkname),'filename','flag', 'smallMask', ...
        'A_save', 'ROIData','-v7.3');
    disp('Instance.mat file saved, consider deleting raw .tif movies!')
    clear A_save
    
    %% photobleach correction
    tic;
    A5 = bleachCorrection(A4);
    toc;
    clear A4
    disp('Photobleaching corrected');
    
    %% Gaussian smoothing
    A6 = GauSmoo(A5,1); 
    A6(A5 == 0) = nan;
    disp('Gaussian smoothing is done');
    clear A5
    disp(' ')
    
    %% top hat filter 
    if spont == 1
        A7 = A6;
        disp('Not doing top-hat filtering here')
    else
        A7 = TopHatFiltering(A6, 300);
        disp('Top hat filtering complete')
        disp(' ')
    end
    clear A6
    
    %% gross dF/F
    
    %Centered data around origins (gross dFoF)
    A_mean = nanmean(A7,3);
    A7 = A7./A_mean - 1;
    
    %% SVD denoising
    % common to do
    % throws away trivial components that should mostly be containing noise (ie
    % keep just the pricinpal components), so now you also have the principle
    % components if you want them
    % time consuming and can be skipped
    % keep iniDim = 1; keep dimension 1 through (at least) 120
    % first component is usually noise
    
    iniDim = 1;
    disp(['Initial Dimention = ' num2str(iniDim)]);
    [A8,U,S,V,iniDim,PC_exp] = roiSVD(A7, iniDim);
    %Reaply downsampled roi mask
    if ~exist('ds_Mask','var')
        ds_Mask = smallMask;
        ds_Mask = repmat(ds_Mask, [1,1,size(A8,3)]);
    end
    A8 = A8.*ds_Mask;
    A8(A8 == 0) = nan;
    checkname = [filename(1:length(filename)-4) '_SVD.mat'];
    save(fullfile(outputFolder,checkname),'U','S','V','PC_exp','iniDim');
    disp('SVD denosing is done')
    
    clear A7
    
    %% df/F of downsampled matrix
    %Recover the data, this is important for later dFoF
    A8 = A8.*A_mean + A_mean;

    %Impose dFOverF to downsampled matrix
    if spont == 1
        sz3 = size(A8);
        A8 = reshape(A8,[sz3(1)*sz3(2),sz3(3)]);
        
        if flag == 1
            %Use bottom 5th percentile as F0
            %load A_F0_new.mat
            %A_F0 = A_F0_after;
            A_F0 = prctile(A8,5,2);
            A_F0_fnm = ['A_F0_5th_', outputFolder, '.mat'];
            save(A_F0_fnm, 'A_F0')
            A_F0 = repmat(A_F0,[1,sz3(3)]);
            disp('Compute dF/F0 based on 5th percentile')
        else
            % Use mean values as F0
            A_F0 = repmat(mean(A8,2),[1,sz3(3)]);
            disp('Compute dF/F0 based on 50th percentile')
            save('A_F0_50th.mat', 'A_F0');
        end
        A_dFoF = reshape(A8./A_F0 - 1,sz3);
        
        A_dFoF = A_dFoF.*ds_Mask;
        disp('Gloabal dFoverF is done')
        disp(' ')
    else
        sz3 = size(A8);
        A_dFoF = A8;
        disp('Do not do gross dFoF for stimulation experiment!')
    end
    clear A8;
    
    % make movie
    tmp = reshape(A_dFoF, [sz3(1)*sz3(2),sz3(3)]);
    s0 = nanstd(tmp(:));
    m0 = nanmean(tmp(:));
    I0=mat2gray(A_dFoF, [-2*s0+m0, 5*s0+m0]);   %scale the whole array so that min = 0, max = 1
    
    [Iarr0, ~] = gray2ind(I0, 256);
    [~, fn, ~] = fileparts(filename);
    Iarr2avi(Iarr0, 1, sz3(3), ['dff_', fn])
    
    movieRange = [-2*s0+m0, 5*s0+m0];
    
    %Save filtered matrix
    A_dFoF = single(A_dFoF);
    checkname = [filename(1:length(filename)-4) '_preprocessed' movTag '.mat'];
    save(fullfile(outputFolder,checkname),'A_dFoF','movieRange', '-v7.3');
    
    %% generate connected components
    
    renewCC(A_dFoF, 0.3, outputFolder, filename)
    
    disp(['Preprocessing done: ' filename]);
    clearvars -except kk fileInfo spatialFactor flag moveAssessFlag
end
%% go to WAIP
