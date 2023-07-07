% widefield spontaneous analysis
% integrate summary data and compare before/after injection

clear all; close all; clc
animalInfo = readtext('preprocessing_list.txt', ' ');
%rig = 20 * 20.4 * 2; %upstairs rig
%rig = 40 * 23.36 * 2; %3rd rig downstairs
%rig = 40 * 16.92 * 2; %downstairs rig

%tag = {'ctl', 'inj'};
tag = {'spont'};
frameRate = 10;
%all_th = [1 1.5 2 3 5];
all_th = 0;

for ff = 1:2%size(animalInfo, 1)
    datapath = animalInfo{ff, 1};
    cd(datapath)
    fd_list = dir(fullfile('*output*'));
    
    for t = 1:length(all_th)
        th = animalInfo{ff, 6};
        cd(['th_', num2str(th)])
        fname = dir(fullfile('*dataSummary.mat*'));
        load(fname(1).name)
        num_rois = length(rp_total{1}.roiCentr);
        %num_rois = animalInfo{ff, 7};
        
        rig = animalInfo{ff, 5};
        if rig == 4
            mag = (20 * 20.4 * 2)/animalInfo{ff,4};
        elseif rig == 3
            mag = (40 * 23.36 * 22)/animalInfo{ff,4};
            mag = mag/10;
        elseif rig == 2
            mag = (40 * 16.92 * 2)/animalInfo{ff, 4};
            if mag >= 180
                warning('Assuming using the old objective!')
                mag = 180/2.3;
            end
            mag = 100/mag;
        end
        
        for g = 1:length(tag)
            if g == 1
                %idx = [7:12];
                idx = [1:size(fd_list,1)];
            else
                idx = [1:6];
            end
            
            for ii = 1:length(idx)
                
                validId = rp_total{idx(ii)}.validId;
                boundBox = rp_total{idx(ii)}.boundBox;
                diameters = rp_total{idx(ii)}.diameters;
                durations = rp_total{idx(ii)}.durations;
                event_dFF = event_dFF_total{idx(ii)};
                m_Interval1 = rp_total{idx(ii)}.m_Interval1;
                m_p_Duration = rp_total{idx(ii)}. m_p_Duration;
                p_Interval1 = rp_total{idx(ii)}.p_Interval1;
                if isfield(rp_total{idx(ii)}, 'pixel')
                    pixel = rp_total{idx(ii)}.pixel;
                else
                    pixel = {0};
                end
                angle = rp_total{idx(ii)}.angle;
                RHO = rp_total{idx(ii)}.RHO;
                roiArea = rp_total{idx(ii)}.roiArea;
                roiCentr = rp_total{idx(ii)}.roiCentr;
                img_sz = 2 * size(m_p_Duration{1, 1});
                sz = size(diameters);
                
                for r = 1:num_rois
                    
                    total_diameter{ii}{r} = [];
                    total_duration{ii}{r} = [];
                    total_area{ii}{r} = [];
                    total_roiArea{ii}{r} = [];
                    total_center_dist{ii}{r} = [];
                    total_center_stepDist_sum{ii}{r} = [];
                    total_center_speed{ii}{r} = [];
                    total_center_stepSpeed{ii}{r} = [];
                    total_max_area{ii}{r} = [];
                    total_center_on{ii}{r} = [];
                    total_event_dFF{ii}{r} = [];
                    total_angle{ii}{r} = [];
                    
                    regionFreq{ii}{r} = length(diameters{r});
                    total_diameter{ii}{r} = [total_diameter{ii}{r}; diameters{r}*mag];
                    total_duration{ii}{r} = [total_duration{ii}{r}; durations{r}/frameRate];
                    total_area{ii}{r} = [total_area{ii}{r}; roiArea{r}' * (mag^2) ./durations{r}];
                    total_roiArea{ii}{r} = [total_roiArea{ii}{r}; roiArea{r}' * (mag^2)];
                    if ~isempty(event_dFF)
                        total_event_dFF{ii}{r} = [total_event_dFF{ii}{r}; event_dFF{r}];
                    else
                        total_event_dFF{ii}{r} = [total_event_dFF{ii}{r}; 0];
                    end
                    total_angle{ii}{r} = [total_angle{ii}{r}; angle{r}];
                    
                    if size(pixel{r}, 2) > 1
                        for w = 1:size(pixel{r}, 2)
                            minT = min(pixel{r}{w}(:, 3));
                            maxT = max(pixel{r}{w}(:, 3));
                            timeSteps = minT : maxT;
                            for t = 1 : length(timeSteps)
                                p_id{t} = find(pixel{r}{w}(:, 3) == timeSteps(t));
                                length_p{ii}{r}(w, t) = length(p_id{t});
                                center_time{ii}{r}{w}(t, :) = [mean(pixel{r}{w}(p_id{t}, 1)), mean(pixel{r}{w}(p_id{t}, 2))];
                                if t > 1
                                    center_stepDist{ii}{r}{w}(t-1) = pdist([center_time{ii}{r}{w}(t-1, :); center_time{ii}{r}{w}(t, :)],'euclidean')...
                                        * mag;
                                end
                            end
                            
                            center_on{ii}{r}(w, :) = center_time{ii}{r}{w}(1, :);
                            
                            % centroid travel distance
                            center_dist{ii}{r}(w) = pdist([center_time{ii}{r}{w}(1, :); center_time{ii}{r}{w}(end, :)],'euclidean') * mag;
                            center_stepDist_sum{ii}{r}(w) = sum(center_stepDist{ii}{r}{w});
                            
                            % centroid moving speed
                            center_speed{ii}{r}(w) = center_dist{ii}{r}(w) / length(timeSteps) * frameRate; % unit: um/s
                            center_stepSpeed{ii}{r}(w) = mean(center_stepDist{ii}{r}{w}) * frameRate; % averaged from each time point, unit: um/s
                            
                            % maximum area during each wave
                            max_area{ii}{r}(w) = max(length_p{ii}{r}(w, :)) * (mag^2);
                        end
                    else
                        center_dist{ii}{r} = 0;
                        center_stepDist_sum{ii}{r} = 0;
                        center_speed{ii}{r} = 0;
                        center_stepSpeed{ii}{r} = 0;
                        max_area{ii}{r} = 0;
                        center_on{ii}{r} = 0;
                    end
                    
                    total_center_dist{ii}{r} = [total_center_dist{ii}{r}, center_dist{ii}{r}];
                    total_center_stepDist_sum{ii}{r} = [total_center_stepDist_sum{ii}{r}, center_stepDist_sum{ii}{r}];
                    total_center_speed{ii}{r} = [total_center_speed{ii}{r}, center_speed{ii}{r}];
                    total_center_stepSpeed{ii}{r} = [total_center_stepSpeed{ii}{r}, center_stepSpeed{ii}{r}];
                    total_max_area{ii}{r} = [total_max_area{ii}{r}, max_area{ii}{r}];
                    total_center_on{ii}{r} = [total_center_on{ii}{r}; center_on{ii}{r}];
                    
                end
                
                % pixel wave interval
                totalInterval1{ii}{r} = zeros(1, img_sz(1)*img_sz(2)/4);
                for p = 1:length(p_Interval1{1, 1})
                    for r = 1:num_rois
                        interval = [];
                        interval = [interval, p_Interval1{r}{p}];
                        totalInterval1{ii}{r}(p) = mean(interval);
                    end
                end
                
                % pixel active duration
                for r = 1:size(durations, 1)
                    d = zeros(size(m_p_Duration{r, 1}));
                    if ~isempty(m_p_Duration{r})
                        d = d + m_p_Duration{r};
                    end
                    total_p_duration{ii}{r} = d;
                end
                
            end
            cd ..
            cd(['th_', num2str(th)])
            
            save([tag{g}, '_th', num2str(th), '_dataSummary.mat'], 'totalInterval1', 'total_p_duration', 'total_diameter', 'total_duration', 'total_area', 'regionFreq', ...
                'total_center_dist', 'total_center_on', 'total_center_stepDist_sum', 'total_center_speed', 'total_center_stepSpeed', 'total_max_area', ...
                'center_dist', 'center_on', 'center_stepDist_sum', 'center_speed', 'center_stepSpeed', 'max_area', 'total_event_dFF', 'total_angle', ...
                'total_roiArea')
        end
        
        clearvars -except animalInfo rig tag frameRate mag num_rois ff all_th th t fd_list
        
        % summary data for between animal comparison
        
        for g = 1:length(tag)
            fname = [tag{g},'_th', num2str(th), '_dataSummary.mat'];
            load(fname)
            sz = size(total_duration, 2);
            num_rois = length(total_duration{1});
            
            for r = 1:num_rois
                Areas{g, r} = [];
                Areas2{g, r} = [];
                Durations{g, r} = [];
                Intervals{g, r} = [];
                NumEvents{g, r} = [];
                MeanEventDff{g, r} = [];
                DomAngles{g, r} = [];
                OpAngles{g, r} = [];
                DirInd{g, r} = [];
                Angles{g, r} = [];
                Diameters{g, r} = [];
                for n = 1:sz
                    %total_area{n}{r} = total_area{n}{r} ./ total_duration{n}{r};
                    medianRoiArea(g, n, r) = median(total_roiArea{n}{r});
                    medianArea(g, n, r) = median(total_area{n}{r});
                    medianArea_max(g, n, r) = median(total_max_area{n}{r});
                    Areas{g, r} = [Areas{g, r}; total_area{n}{r}];
                    Areas2{g, r} = [Areas2{g, r}; total_roiArea{n}{r}];
                    Diameters{g, r} = [Diameters{g, r}; total_diameter{n}{r}];
                    
                    medianDuration(g, n, r) = median(total_duration{n}{r});
                    Durations{g, r} = [Durations{g, r}; total_duration{n}{r}];
                    
                    validP = (totalInterval1{n}{r} > 10) & (totalInterval1{n}{r} < 2000);
                    intervals1{n, r} = totalInterval1{n}{r}(validP) / frameRate;
                    medianInterval(g, n, r) = median(intervals1{n, r});
                    Intervals{g, r} = [Intervals{g, r}, intervals1{n, r}];
                    
                    median_center_dist(g, n, r) = median(total_center_dist{n}{r});
                    median_center_stepDist_sum(g, n, r) = median(total_center_stepDist_sum{n}{r});
                    median_center_speed(g, n, r) = median(total_center_speed{n}{r});
                    median_center_stepSpeed(g, n, r) = median(total_center_stepSpeed{n}{r});
                    
                    NumEvents{g, r} = [NumEvents{g, r}; length(total_duration{n}{r})];
                    MeanEventDff{g, r} = [MeanEventDff{g,r}; total_event_dFF{n}{r}'];
                    median_mean_event_dff(g, n, r) = median(total_event_dFF{n}{r});
                    
                    Angles{g, r} = [Angles{g, r}; total_angle{n}{r}'];
                    find_angles = rad2deg(wrapTo2Pi(total_angle{n}{r}));
                    tmp = find(find_angles > 240);
                    tmp2 = find(find_angles < 300);
                    tmp3 = find(find_angles < 120);
                    tmp4 = find(find_angles > 60);
                    tmp5 = intersect(tmp, tmp2);
                    tmp6 = intersect(tmp3, tmp4);
                    % fraction of waves in dominant direction
                    DomAngles{g, r} = [DomAngles{g, r}; length(tmp5)/length(find_angles)];
                    OpAngles{g, r} = [OpAngles{g, r}; length(tmp6)/length(find_angles)];
                    clear find_angles tmp tmp2 tmp3 tmp4 tmp5 tmp6
                end
                DirInd_tmp = (DomAngles{g, r} - OpAngles{g,r})./(DomAngles{g, r} + OpAngles{g, r});
                if sum(isnan(DirInd_tmp)) >= 1
                    idx = find(isnan(DirInd_tmp));
                    DirInd_tmp(idx) = 0;
                end
                DirInd{g, r} = [DirInd{g, r}; DirInd_tmp'];
                clear DirInd_tmp idx
            end
        end
        
        %data structs are such that each row is a condition and each column is
        %an roi
        
        % get mean across rois and events, each row is a condition
        medianRoiArea2 = nanmean(nanmean(medianRoiArea,3),2);
        % median area
        medianArea2 = nanmean(nanmean(medianArea, 3),2);
        % median area max
        medianArea_max2 = nanmean(nanmean(medianArea_max, 3),2);
        % median center dist
        median_center_dist2 = nanmean(nanmean(median_center_dist, 3), 2);
        % median center step speed
        median_center_stepSpeed2 = nanmean(nanmean(median_center_stepSpeed, 3), 2);
        % median center step dist
        median_center_stepDist_sum2 = nanmean(nanmean(median_center_stepDist_sum, 3), 2);
        % median duration
        medianDuration2 = nanmean(nanmean(medianDuration, 3), 2);
        % median interval
        medianInterval2 = nanmean(nanmean(medianInterval, 3), 2);
        % num events
        NumEvents2 = zeros(size(NumEvents));
        for g = 1:length(tag)
            for r = 1:num_rois
                NumEvents2(g, r) = mean(NumEvents{g,r});
            end
        end
        NumEvents2 = mean(NumEvents2, 2);
        % mean event dff
        median_mean_event_dff2 = nanmean(nanmean(median_mean_event_dff, 3), 2);
        % median fraction of waves in dominant direction
        median_dom_angles = zeros(size(DomAngles));
        for g = 1:length(tag)
            for r = 1:num_rois
                median_dom_angles(g, r) = nanmedian(DomAngles{g, r});
            end
        end
        median_dom_angles = median(median_dom_angles, 2);
        % median directionality index
        median_dir_ind = zeros(size(DirInd));
        for g = 1:length(tag)
            for r = 1:num_rois
                median_dir_ind(g, r) = nanmedian(DirInd{g, r});
            end
        end
        median_dir_ind = median(median_dir_ind, 2);
        
        save(['within_animal_summ_th', num2str(th), '.mat'], 'medianArea2', 'medianArea_max2', 'median_center_dist2', ...
            'median_center_stepSpeed2', 'median_center_stepDist_sum2', 'medianDuration2', 'medianInterval2', ...
            'NumEvents2', 'median_mean_event_dff2', 'Angles', 'median_dom_angles', 'medianRoiArea2', 'Diameters',...
            'Durations', 'median_dir_ind');
        
        if length(tag) > 1
            % within animal line plots
            colorVector = [0 0 255; 255 0 0] / 255;
            colorVector_line = [105 105 105; 105 105 105] / 255;
            colors = colorVector([1, 2], :);
            % areas
            figure;
            line_cords_x = [1 2];
            for g = 1:length(tag)
                scatter(g, medianArea2(g), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
                hold on
                scatter(g, medianArea_max2(g), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
            end
            line_cords_y = [medianArea2 medianArea_max2];
            plot(line_cords_x, line_cords_y, 'k')
            
            txt = 'median area (um^2)   '; text(1, medianArea2(1,1), txt, 'FontSize', 8, 'HorizontalAlignment', 'right')
            txt = 'peak area (um^2)   '; text(1, medianArea_max2(1, 1), txt, 'FontSize', 8, 'HorizontalAlignment', 'right')
            xlim([0 3]); h = gca; h.XAxis.TickLength = [0 0]; ylabel('area (um^2)')
            set(h, 'XTickLabel', {' ', ' ', 'ctl', ' ', 'inj', ' ', ' '})
            saveas(gcf, 'areas.png')
            
            % distances
            figure;
            line_cords_x = [1 2];
            for g = 1:length(tag)
                scatter(g, median_center_dist2(g), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
                hold on
                scatter(g, median_center_stepDist_sum2(g), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
            end
            line_cords_y = [median_center_dist2 median_center_stepDist_sum2];
            plot(line_cords_x, line_cords_y, 'k')
            
            txt = 'dist (um)   '; text(1, median_center_dist2(1,1), txt, 'FontSize', 8, 'HorizontalAlignment', 'right')
            txt = 'step dist (um)  '; text(1, median_center_stepDist_sum2(1, 1), txt, 'FontSize', 8, 'HorizontalAlignment', 'right')
            xlim([0 3]); h = gca; h.XAxis.TickLength = [0 0]; ylabel('dist (um)')
            set(h, 'XTickLabel', {' ', ' ', 'ctl', ' ', 'inj', ' ', ' '})
            saveas(gcf, 'distances.png')
            
            % speeds, intervals, durations
            figure;
            line_cords_x = [1 2];
            for g = 1:length(tag)
                scatter(g, median_center_stepSpeed2(g), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
                hold on
                scatter(g, medianDuration2(g), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
                scatter(g, medianInterval2(g), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
            end
            line_cords_y = [median_center_stepSpeed2 medianDuration2 medianInterval2];
            plot(line_cords_x, line_cords_y, 'k')
            
            txt = 'speed (s)   '; text(1, median_center_stepSpeed2(1,1), txt, 'FontSize', 8, 'HorizontalAlignment', 'right')
            txt = 'duration (s)  '; text(1, medianDuration2(1, 1), txt, 'FontSize', 8, 'HorizontalAlignment', 'right')
            txt = 'interval (s)  '; text(1, medianInterval2(1, 1), txt, 'FontSize', 8, 'HorizontalAlignment', 'right')
            xlim([0 3]); h = gca; h.XAxis.TickLength = [0 0]; ylabel('seconds')
            set(h, 'XTickLabel', {' ', ' ', 'ctl', ' ', 'inj', ' ', ' '})
            saveas(gcf, 'speeds_intervals_durations.png')
            
            %num events and mean dFF
            figure; h1 = subplot(1, 2, 1);
            line_cords_x = [1 2];
            for g = 1:length(tag)
                scatter(g, NumEvents2(g), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
                hold on
            end
            line_cords_y = [NumEvents2];
            plot(line_cords_x, line_cords_y, 'k')
            
            txt = 'number events   '; text(1, NumEvents2(1,1), txt, 'FontSize', 8, 'HorizontalAlignment', 'right')
            xlim([0 3]); h1 = gca; h1.XAxis.TickLength = [0 0];
            set(h1, 'XTickLabel', {' ', 'ctl', 'inj', ' '})
            
            h2 = subplot(1, 2, 2);
            for g = 1:length(tag)
                scatter(g, median_mean_event_dff2(g), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
                hold on
            end
            line_cords_y = [median_mean_event_dff2];
            plot(line_cords_x, line_cords_y, 'k')
            
            txt = 'mean event dFF   '; text(1, median_mean_event_dff2(1,1), txt, 'FontSize', 8, 'HorizontalAlignment', 'right')
            xlim([0 3]); h2 = gca; h2.XAxis.TickLength = [0 0]; ylabel('dFF')
            set(h2, 'XTickLabel', {' ', 'ctl', 'inj', ' '})
            saveas(gcf, 'numevents_meandFF.png')
        end
        disp(['finished threshold ', num2str(th)])
        cd ..
    end
    disp(['finished animal ', num2str(ff), ' of ', num2str(size(animalInfo,1))])
    clearvars -except animalInfo tag frameRate ff all_th t
    
end




