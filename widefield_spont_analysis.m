% widefield spontaneous processing pipeline, integrating Xinxin and Yixiang's
% pipelines

clear all; close all; clc
animalInfo = readtext('preprocessing_list.txt', ' ');
sigma = 1; % for gaussian smoothing
dura_th = 8; % duration threshold in frames
dia_th = 10; % diameter threshold in pixels, Yixiang used 10
opticFlowFlag = 1; %whether or not to measure wave directionality


for ff = 1:2%:size(animalInfo, 1)
    datapath = animalInfo{ff, 1};
    cd(datapath)
    need_test_th = animalInfo{ff, 3};
    rig = animalInfo{ff, 5};
    fd_list = dir(fullfile('*output*'));
    
    
    %% test threshold on first movie
    if need_test_th == 1
        cd(fd_list(1).name)
        % load dA file, already downsampled
        flist = dir(fullfile('*preprocessed*.mat'));
        load(flist(1).name)
        imgall = A_dFoF(:,:,301:end);
        
        %Get disconnected rois from the loaded movie
        A1 = imgall(:, :, 1);
        sz2 = size(A1);
        if any(isnan(A1(:)))
            B1 = ~isnan(A1);
        else
            B1 = ~(A1 == 0);
        end
        
        %Get connected componets from frame1
        C1 = bwconncomp(B1);
        
        %Construct the cell array of different rois
        n_roi = animalInfo{ff, 7};
        if animalInfo{ff, 8} > 1
            cur_roi = zeros(sz2);
            cur_roi(C1.PixelIdxList{2}) = 1;
            roi{1} = cur_roi;
        else
            for i = 1 : n_roi
                cur_roi = zeros(sz2);
                cur_roi(C1.PixelIdxList{i}) = 1;
                roi{i} = cur_roi;
            end
        end
        
        % Get the total mask combining all rois
        totalMask = zeros(size(roi{1}));
        for r = 1:length(roi)
            totalMask = totalMask + roi{r};
        end
        
        % Zscore
        sz = size(imgall);
        A_mean = nanmean(imgall(:));
        nan_id = isnan(imgall);
        imgall((nan_id)) = A_mean;
        A_z = zscore(imgall(:));
        A_z = reshape(A_z, sz);
        A_z(nan_id) = nan;
        
        test_th = [1.5 2 3 5];
        for t = 1:length(test_th)
            thresh = test_th(t);
            for r = 1:length(roi)
                if r == length(roi) + 1
                    maskMatrix = totalMask;
                else
                    maskMatrix = roi{r};
                end
                
                maskId = find(maskMatrix > 0);
                maskMatrix = reshape(maskMatrix, sz(1)*sz(2), 1);
                maskMatrix = repmat(maskMatrix, 1, sz(3));
                maskMatrix = reshape(maskMatrix, sz(1), sz(2), sz(3));
                
                % apply the mask
                subMov = imgall.*maskMatrix;
                
                %Zscore
                subMov = reshape(subMov, sz(1) * sz(2), sz(3));
                temp_subMov = subMov(maskId, :);
                
                sz3 = size(temp_subMov);
                mean_subMov = nanmean(temp_subMov(:));
                nan_id = isnan(temp_subMov);
                temp_subMov((nan_id)) = mean_subMov;
                temp_subMov = zscore(temp_subMov(:));
                temp_subMov = reshape(temp_subMov, sz3);
                temp_subMov(nan_id) = nan;
                
                subMov = zeros(size(subMov));
                subMov(maskId, :) = temp_subMov;
                
                %Binarization
                activeMov = subMov > thresh;
                subMov = reshape(subMov, sz(1), sz(2), sz(3));
                activeMov = reshape(activeMov,sz(1), sz(2), sz(3));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % detect connected components
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Get connected components
                CC = bwconncomp(activeMov);
                
                %Get regional properties
                STATS = regionprops(CC, activeMov, 'Area', 'BoundingBox', 'Centroid', 'PixelList');
                
                %Compute duration, diameter, area and centroid
                roiBoundingBox = zeros(length(STATS),6);
                
                for i = 1:length(STATS)
                    roiBoundingBox(i, :) = STATS(i).BoundingBox;
                end
                
                dura = roiBoundingBox(:,6);
                dia = mean([roiBoundingBox(:,4) roiBoundingBox(:,5)], 2);
                area = vertcat(STATS.Area);
                
                %Save all properties
                validId{r, 1} = (dura > dura_th & dia > dia_th);
                durations{r, 1} = dura(validId{r, 1});
                diameters{r, 1} = dia(validId{r, 1});
                roiCentr{r, 1} = vertcat(STATS(validId{r, 1}).Centroid);
                roiArea{r, 1} = [STATS(validId{r, 1}).Area];
                boundBox{r, 1} = [STATS(validId{r, 1}).BoundingBox];
                valid_tmp = find(validId{r, 1} > 0);
                
                for v = 1:length(valid_tmp)
                    pixel{r, 1}{v} = STATS(valid_tmp(v)).PixelList; % pixel index of all valid events
                    %each cell is pixel x, pixel y, and frame num
                end
                
                valid{r, 1} = valid_tmp; %all valid events
                
                
                %Reconstruct binary mov based on valid connected components
                valid_mov = subMov;
                valid_activeMov{r} = activeMov;
                badId = find(validId{r, 1} == 0);
                for id = 1:length(badId)
                    removePixel = CC.PixelIdxList{badId(id)};
                    valid_mov(removePixel) = 0;
                    valid_activeMov{r}(removePixel) = 0;
                end
                
                %valid_activeMov_down = imresize(valid_activeMov{r}, 0.25, 'bilinear');
                %tmpsz = size(valid_activeMov_down);
                %valid_activeMov_down = reshape(valid_activeMov_down, tmpsz(1)*tmpsz(2), tmpsz(3));
                
                if r == 1
                    total_ActiveMovie = valid_activeMov{r};
                else
                    total_ActiveMovie = total_ActiveMovie + valid_activeMov{r};
                end
                clear activeMov subMov temp_subMov maskId valid_mov badId removePixel
            end
            
            % create segmented movie
            total_ActiveMovie = total_ActiveMovie > 0;
            for fr = 1:sz(3)
                [I2, map2] = gray2ind(total_ActiveMovie(:,:,fr), 8); %figure; imshow(I2,map)
                F(fr) = im2frame(I2,map2);  %setup the binary segmented mask movie
            end
            mov_fnm = [flist(1).name(1:end-4), 'mask_th', num2str(thresh), '.avi'];
            writeMovie_xx(F, mov_fnm, 0);
        end
        test_th = input('Please input the threshold level you pick for this animal:');
    else
       test_th = animalInfo{ff, 6};
       %test_th = [1.5 2 3];
    end
    
    clearvars -except animalInfo dura_th dia_th sigma ff ...
        datapath fd_list th opticFlowFlag need_test_th test_th
    
    cd (datapath)
    
    %% after picking threshold, do property analysis with same threshold
    for t = 1:length(test_th)
        th = test_th(t);
        for f_id = 1:length(fd_list)
            cd(fd_list(f_id).name)
            flist = dir(fullfile('*preprocessed*.mat'));
            load(flist(1).name)
            save_subfd_name = flist(1).name(1:end-4);
            
            %Get disconnected rois from the loaded movie
            A1 = A_dFoF(:, :, 1);
            sz2 = size(A1);
            if any(isnan(A1(:)))
                B1 = ~isnan(A1);
            else
                B1 = ~(A1 == 0);
            end
            
            %Get connected componets from frame1
            C1 = bwconncomp(B1);
            
            %Construct the cell array of different rois
            %n_roi = size(C1.PixelIdxList, 2);
            n_roi = animalInfo{ff, 7};
            if animalInfo{ff, 8} > 1
                cur_roi = zeros(sz2);
                cur_roi(C1.PixelIdxList{2}) = 1;
                roi{1} = cur_roi;
            else  
                for i = 1 : n_roi
                    cur_roi = zeros(sz2);
                    cur_roi(C1.PixelIdxList{i}) = 1;
                    roi{i} = cur_roi;
                end
            end
            
            % Z-score the current movie
            imgall2 = A_dFoF(:,:,501:end);
            sz4 = size(imgall2);
            A_mean = nanmean(imgall2(:));
            nan_id = isnan(imgall2);
            imgall2((nan_id)) = A_mean;
            A_z = zscore(imgall2(:));
            A_z = reshape(A_z, sz4);
            A_z(nan_id) = nan;
            
            sz = size(A_z);
            
            smallSize = 0.5; % downSample the flowfield size %%%%%%%%%%% normally 0.5
            imgall = reshape(A_z, sz(1), sz(2), sz(3));
            imgall(isnan(imgall)) = 0;
            
            if opticFlowFlag == 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % compute flow field
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Computes optic flow field (with normalized vectors) for imgall,
                %Xinxin's version
                [normVx, normVy] = computeFlowField_normalized_xx(imgall, sz, smallSize); %LK method
                totalMask = ~isnan(A_z(:,:,1));
                smallMask = imresize(totalMask, smallSize, 'bilinear');
                clear A_z;
                
                %Plot the optic flow field
                h = figure; imagesc(smallMask); hold on
                quiver(mean(normVx, 3).*smallMask, mean(normVy, 3).*smallMask); axis image;
                title(save_subfd_name)
                set(h, 'Position', [0, 0, 1200, 900]);
                h.PaperPositionMode = 'auto';
                saveas(gcf, [save_subfd_name, '_quiver.png'])
                
                %Plot the optic flow field (further downsampled)
                h = figure; imagesc(imresize(smallMask, .5, 'bilinear')); hold on
                quiver(imresize(mean(normVx, 3), .5, 'bilinear') .* imresize(smallMask, .5, 'bilinear'), ...
                    imresize(mean(normVy, 3), .5, 'bilinear') .* imresize(smallMask, .5, 'bilinear'));
                axis image;
                title(save_subfd_name)
                set(h, 'Position', [0, 0, 1200, 900]);
                h.PaperPositionMode = 'auto';
                saveas(gcf, [save_subfd_name, '_quiver_s.png'])
            end
            
            %Binarize and filter
            % Initialize
            mkdir(['th_', num2str(th)])
            cd(['th_', num2str(th)])
            total_ActiveMovie = [];
            angle =[]; n = 1; frameRate = 10;
            RHO = []; p_Interval1 = []; m_Interval1 = []; m_p_Duration = [];
            sz = size(imgall);
            
            % Get the total mask combining all rois
            totalMask = zeros(size(roi{1}));
            for r = 1:length(roi)
                totalMask = totalMask + roi{r};
            end
            
            for r = 1:length(roi)
                if r == length(roi) + 1
                    maskMatrix = totalMask;
                else
                    maskMatrix = roi{r};
                end
                maskId = find(maskMatrix > 0);
                maskMatrix = reshape(maskMatrix, sz(1)*sz(2), 1);
                maskMatrix = repmat(maskMatrix, 1, sz(3));
                maskMatrix = reshape(maskMatrix, sz(1), sz(2), sz(3));
                
                % apply the mask
                imgall = reshape(imgall, sz(1), sz(2), sz(3));
                subMov = imgall.*maskMatrix;
                
                %Zscore
                subMov = reshape(subMov, sz(1) * sz(2), sz(3));
                temp_subMov = subMov(maskId, :);
                
                sz3 = size(temp_subMov);
                mean_subMov = nanmean(temp_subMov(:));
                nan_id = isnan(temp_subMov);
                temp_subMov((nan_id)) = mean_subMov;
                temp_subMov = zscore(temp_subMov(:));
                temp_subMov = reshape(temp_subMov, sz3);
                temp_subMov(nan_id) = nan;
                
                subMov = zeros(size(subMov));
                subMov(maskId, :) = temp_subMov;
                
                %Binarization
                activeMov = subMov > th;
                subMov = reshape(subMov, sz(1), sz(2), sz(3));
                activeMov = reshape(activeMov,sz(1), sz(2), sz(3));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % detect connected components
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Get connected components
                CC = bwconncomp(activeMov);
                
                %Get regional properties
                STATS = regionprops(CC, activeMov, 'Area', 'BoundingBox', 'Centroid', 'PixelList');
                
                %Compute duration, diameter, area and centroid
                roiBoundingBox = zeros(length(STATS),6);
                
                for i = 1:length(STATS)
                    roiBoundingBox(i, :) = STATS(i).BoundingBox;
                end
                
                dura = roiBoundingBox(:,6);
                dia = mean([roiBoundingBox(:,4) roiBoundingBox(:,5)], 2);
                area = vertcat(STATS.Area);
                
                %Save all properties
                validId{r, n} = (dura > dura_th & dia > dia_th);
                durations{r, n} = dura(validId{r, n});
                diameters{r, n} = dia(validId{r, n});
                roiCentr{r, n} = vertcat(STATS(validId{r, n}).Centroid);
                roiArea{r, n} = [STATS(validId{r, n}).Area];
                boundBox{r, n} = [STATS(validId{r, n}).BoundingBox];
                valid_tmp = find(validId{r, n} > 0);
                
                for v = 1:length(valid_tmp)
                    pixel{r, n}{v} = STATS(valid_tmp(v)).PixelList;
                end
                
                valid{r, n} = valid_tmp;
                
                
                %Reconstruct binary mov based on valid connected components
                valid_mov = subMov;
                valid_activeMov{r} = activeMov;
                badId = find(validId{r, n} == 0);
                for id = 1:length(badId)
                    removePixel = CC.PixelIdxList{badId(id)};
                    valid_mov(removePixel) = 0;
                    valid_activeMov{r}(removePixel) = 0;
                end
                
                if r == 1
                    total_ActiveMovie = valid_activeMov{r};
                else
                    total_ActiveMovie = total_ActiveMovie + valid_activeMov{r};
                end
                
                % do wave property analysis
                if opticFlowFlag == 1
                    % compute flow field without normalization
                    % Xinxin's version
                    [AVx, AVy] = computeFlowField_xx(imgall, sz); %LK method
                    nDomains = sum(validId{r, n});
                    validDomains = CC.PixelIdxList(validId{r, n});
                    theta = []; rho = [];
                    %Go through each domain and transform Cartesian coordinates to polar or cylindrical
                    for k = 1:nDomains
                        p_id = intersect(1 : sz(1)*sz(2)*size(AVx, 3), validDomains{k});
                        [theta(k), rho(k)]= cart2pol(sum(AVx(p_id)), sum(AVy(p_id))); %transformation
                        if mod(k, 100) == 0
                            disp(['loop ', num2str(k), ' of ', num2str(nDomains), ' complete']);
                        end
                    end
                    
                    clear AVx AVy
                    
                    angle{r, n} = theta;
                    RHO{r, n} = rho;
                    
                    h = figure;
                    if ~isempty(angle{r,n})
                        rose(angle{r, n});
                        set(gca,'YDir','reverse'); % 90 degree is moving downwards
                        title([save_subfd_name, '_th', num2str(th)]);
                        saveas(h, [save_subfd_name, 'roi', num2str(r), 'th', num2str(th), '_rosePlot.png'])
                    end
                end
                
                %Plot durations and diameters of detected components
                savefn2 = [save_subfd_name, 'roi', num2str(r), '_th', num2str(th)];
                h(1) = figure; hist(durations{r, n}, 50); xlabel('durations (frames)'); title(['thresh=', num2str(th)])
                saveas(h(1), [savefn2, '_durations.png'])
                h(2) = figure; hist(diameters{r, n}, 50); xlabel('diameters (pixels)'); title(['thresh=', num2str(th)])
                saveas(h(2), [savefn2, '_diameters.png'])
                h(3) = figure; scatter(durations{r, n}, diameters{r, n}); xlabel('durations'); ylabel('diameters'); title(['thresh=', num2str(th)])
                saveas(h(3), [savefn2, '_duraVSdia.png'])
                
                %Reconstruct binary mov based on valid connected components
                valid_mov = subMov;
                valid_activeMov{r} = activeMov;
                clear activeMov subMov
                
                %Find noise
                badId = find(validId{r, n} == 0);
                for id = 1:length(badId)
                    removePixel = CC.PixelIdxList{badId(id)};
                    valid_mov(removePixel) = 0;
                    valid_activeMov{r}(removePixel) = 0;
                end
                
                valid_activeMov_down{r} = imresize(valid_activeMov{r}, 0.5, 'bilinear');
                tmpsz = size(valid_activeMov_down{r});
                valid_activeMov_down{r} = reshape(valid_activeMov_down{r}, tmpsz(1)*tmpsz(2), tmpsz(3));
                
                %Combine differ rois
                if isempty(total_ActiveMovie)
                    total_ActiveMovie = valid_activeMov{r};
                else
                    total_ActiveMovie = total_ActiveMovie + valid_activeMov{r};
                end
                
                
                % compute pixel active duration and pixel event interval
                for p = 1:tmpsz(1)*tmpsz(2)
                    activeOn{p} = find(valid_activeMov_down{r}(p, 2:end) - valid_activeMov_down{r}(p, 1:end-1) > 0) + 1;
                    activeOff{p} = find(valid_activeMov_down{r}(p, 2:end) - valid_activeMov_down{r}(p, 1:end-1) < 0);
                    
                    if (isempty(activeOn{p}) + isempty(activeOff{p})) == 1
                        
                        activeOn{p} = [];
                        activeOff{p} = [];
                        
                    elseif (isempty(activeOn{p}) + isempty(activeOff{p})) == 0
                        
                        if activeOn{p}(1) > activeOff{p}(1)
                            activeOff{p} = activeOff{p}(2:end);
                            if isempty(activeOff{p})
                                activeOn{p} = [];
                                activeOff{p} = [];
                            end
                        end
                        
                        if ~isempty(activeOff{p})
                            if activeOn{p}(end) > activeOff{p}(end)
                                activeOn{p} = activeOn{p}(1:end-1);
                                if isempty(activeOn{p})
                                    activeOn{p} = [];
                                    activeOff{p} = [];
                                end
                            end
                        end
                    end
                    
                    pixelDuration{p} = activeOff{p} - activeOn{p};
                    
                    goodId = pixelDuration{p} > 2;
                    activeOn{p} = activeOn{p}(goodId);
                    activeOff{p} = activeOff{p}(goodId);
                    pixelDuration{p} = pixelDuration{p}(goodId);
                    meanDuration(p) = sum(pixelDuration{p}) / sz(3);
                    
                    %Interval: from end of one event to the beginning of the following event
                    pixelInterval1{p} = activeOn{p}(2:end) - activeOff{p}(1:end-1);
                    meanInterval1(p) = mean(pixelInterval1{p});
                    
                    %Interval: between the center of each event
                    eventCenterTime{p} = activeOff{p} - activeOn{p};
                    pixelInterval2{p} = eventCenterTime{p}(2 : end) - eventCenterTime{p}(1 : end-1);
                    meanInterval2(p) = mean(pixelInterval2{p});
                end
                
                %Plot event intervals
                p_Interval1{r, n} = pixelInterval1;
                meanInterval1(isnan(meanInterval1)) = 2000/frameRate;
                meanInterval1 = reshape(meanInterval1, tmpsz(1), tmpsz(2))/frameRate;
                h = figure; imagesc(meanInterval1); colorbar; colormap jet
                caxis([0, 50]); axis image
                title(savefn2);
                saveas(h, [savefn2, '_interval.png']);
                
                %Plot mean durations
                meanDuration = reshape(meanDuration, tmpsz(1), tmpsz(2));
                h = figure; imagesc(meanDuration); colorbar; colormap jet
                caxis([0, 0.3]); axis image
                title(savefn2);
                saveas(h, [savefn2, '_duration.png']);
                
                m_Interval1{r, n} = meanInterval1;
                m_p_Duration{r, n} = meanDuration;
                
            end
            
            clear valid_activeMov imgall
            
            %Record all properteis as a struct
            rp.validId = validId;
            rp.durations = durations;
            rp.diameters = diameters;
            rp.roiCentr = roiCentr;
            rp.roiArea = roiArea;
            rp.boundBox = boundBox;
            rp.valid = valid;
            if exist('pixel', 'var')
            rp.pixel = pixel;
            end
            rp.angle = angle;
            rp.RHO = RHO;
            rp.p_Interval1 = p_Interval1;
            rp.m_Interval1 = m_Interval1;
            rp.m_p_Duration = m_p_Duration;
            
            % create segmented movie
            total_ActiveMovie = total_ActiveMovie > 0;
            for fr = 1:sz(3)
                [I2, map2] = gray2ind(total_ActiveMovie(:,:,fr), 8); %figure; imshow(I2,map)
                F(fr) = im2frame(I2,map2);  %setup the binary segmented mask movie
            end
            mov_fnm = [flist(1).name(1:end-4), 'mask_th', num2str(th), '.avi'];
            writeMovie_xx(F, mov_fnm, 0);
            
            %Downsample
             total_ActiveMovie = imresize(total_ActiveMovie, .5, 'bilinear');
            
            if opticFlowFlag == 1
                %Get wave properties (angle, rho, vector matrices)
                angle = rp.angle;
                RHO = rp.RHO;
                
                %Save opticflow properties for this movie
                save([savefn2, '_opticFlow.mat'], 'angle', 'normVx', 'normVy', 'RHO', ...
                    'total_ActiveMovie', '-v7.3');
                
                clear normVx normVy
            end
            
            clear total_ActiveMovie
            
            % Get more wave properties
            validId = rp.validId;
            durations = rp.durations ;
            diameters = rp.diameters ;
            roiCentr = rp.roiCentr ;
            roiArea = rp.roiArea;
            boundBox = rp.boundBox;
            valid = rp.valid;
            if exist('pixel', 'var')
            pixel = rp.pixel;
            end
            
            p_Interval1 = rp.p_Interval1;
            m_Interval1 = rp.m_Interval1;
            m_p_Duration = rp.m_p_Duration;
            
            
            % get mean dF/F of all events
            % initialize
            if exist('pixel', 'var')
            for r = 1:length(pixel)
                num_events(r) = length(pixel{r});
                event_dFF{r} = zeros(1, num_events(r));
            end
            
            for r = 1:length(pixel)
                events = pixel{r};
                for ee = 1:num_events(r)
                    event_coords = events{ee};
                    %                 if length(event_coords) > 100000
                    %                     disp('large event')
                    %                 end
                    total_dFF = zeros(1,length(event_coords));
                    for zz = 1:length(event_coords)
                        total_dFF(zz) = A_dFoF(event_coords(zz,2), event_coords(zz,1), event_coords(zz, 3));
                        %                     if length(event_coords) > 100000
                        %                         if mod(zz, 10000) == 0
                        %                             disp(['event coordinates ', num2str(zz), ' of ', num2str(length(event_coords)), ' complete']);
                        %                         end
                        %                     end
                    end
                    
                    total_dFF = total_dFF(total_dFF > 0);
                    event_dFF{r}(ee) = mean(total_dFF);
                    clear event_coords total_dFF
                end
                clear events
            end
            end
            
            %Save all wave properties as a data summary file
            if exist('event_dFF', 'var') && exist('pixel', 'var')
            save([savefn2, '_dataSummary.mat'], 'p_Interval1', 'm_Interval1', 'event_dFF', ...
                'validId', 'durations', 'diameters', 'roiCentr', 'roiArea', 'angle', 'RHO', 'm_p_Duration', 'boundBox', 'pixel', '-v7.3');
            else
                save([savefn2, '_dataSummary.mat'], 'p_Interval1', 'm_Interval1', ...
                'validId', 'durations', 'diameters', 'roiCentr', 'roiArea', 'angle', 'RHO', 'm_p_Duration', 'boundBox', '-v7.3');
            end
            
            %Store the regionprop struct
            rp_total{f_id} = rp;
            if exist('event_dFF', 'var')
            event_dFF_total{f_id} = event_dFF;
            end
            
            %Go back to root directory
            cd(datapath)
            close all
            clearvars -except animalInfo sz dura_th dia_th sigma ff ...
                datapath flist th opticFlowFlag f_id fd_list rp_total ...
                need_test_th event_dFF_total t test_th
        end
        
        % save summary data for this animal
        mkdir(['th_', num2str(th)])
        cd(['th_', num2str(th)])
        if contains(datapath, 'dhbe')
            save_fnm = datapath(9:end);
        elseif contains(datapath, 'het')
            save_fnm = datapath(18:end);
        elseif contains(datapath, 'homo')
            save_fnm = datapath(19:end);
        end
        save([save_fnm, '_th', num2str(th), '_dataSummary.mat'], 'rp_total', 'event_dFF_total', '-v7.3');
        disp(['Finished threshold ', num2str(th)])
        cd(datapath)
    end
    disp(['Finished animal ', num2str(ff) ' of ', num2str(size(animalInfo, 1))])
    
    clearvars -except animalInfo dura_th dia_th sigma ff opticFlowFlag need_test_th fd_list t test_th
    
end