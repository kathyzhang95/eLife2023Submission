% widefield spontaneous plot comparisons across animals

clear all; close all; clc
tag = {'het', 'homo'};
same_animal = 0;
animalInfo = readtext('preprocessing_list.txt', ' ');
%all_th = [1 1.5 2 3 5];
all_th = 0;
%all_th = 1.5;

for t = 1:length(all_th)
    %th = all_th(t);
    median_center_dist = [];
    median_center_stepDist_sum = [];    
    median_center_stepSpeed = [];
    median_mean_event_dff = [];
    medianArea = [];
    medianRoiArea = [];
    medianArea_max = [];
    medianDuration = [];
    medianInterval = [];
    numEvents = [];
    median_dom_angles2 = [];
    median_dir_ind2 = [];
    Angles2 = cell(1,2);
    Diameters2 = cell(1,2);
    Durations2 = cell(1,2);
    
    
    %load all data
    
    for ff = 1:5%size(animalInfo, 1)
        datapath = animalInfo{ff, 1};
        th = animalInfo{ff, 6};
        cd(datapath)
        cd(['th_', num2str(th)])
        f_list = dir(fullfile('*within_animal_summ*'));
        load(f_list(1).name)
        
        if contains(datapath, 'het')
            median_center_dist = [median_center_dist, median_center_dist2/100];
            median_center_stepDist_sum = [median_center_stepDist_sum, median_center_stepDist_sum2/100];
            median_center_stepSpeed = [median_center_stepSpeed, median_center_stepSpeed2/100];
            median_mean_event_dff = [median_mean_event_dff, median_mean_event_dff2];
            medianArea = [medianArea, medianArea2];
            medianArea_max = [medianArea_max, medianArea_max2];
            medianRoiArea = [medianRoiArea, medianRoiArea2];
            medianDuration = [medianDuration, medianDuration2];
            medianInterval = [medianInterval, medianInterval2];
            numEvents = [numEvents, NumEvents2];
            for r = 1:size(Angles,2)
                Angles2{1} = [Angles2{1}; Angles{1,r}];
                Diameters2{1} = [Diameters2{1}; Diameters{1,r}];
                Durations2{1} = [Durations2{1}; Durations{1,r}];
                %             if same_animal
                %                 Angles2{2} = [Angles2{2}; Angles{2,r}];
                %                 Diameters2{2} = [Diameters2{2}; Diameters{2,r}];
                %                 Durations2{2} = [Durations2{2}; Durations{2,r}];
                %             end
            end
            median_dom_angles2 = [median_dom_angles2, median_dom_angles];
            median_dir_ind2 = [median_dir_ind2, median_dir_ind];
        end
    end
    sz2 = length(median_center_dist);
    
    for ff = 1:5%size(animalInfo, 1)
        datapath = animalInfo{ff, 1};
        th = animalInfo{ff, 6};
        cd(datapath)
        cd(['th_', num2str(th)])
        f_list = dir(fullfile('*within_animal_summ*'));
        load(f_list(1).name)
        
        if contains(datapath, 'homo')
            median_center_dist(2,ff-sz2) = median_center_dist2/100;
            median_center_stepDist_sum(2,ff-sz2) = median_center_stepDist_sum2/100;
            median_center_stepSpeed(2, ff-sz2) = median_center_stepSpeed2/100;
            median_mean_event_dff(2, ff-sz2) = median_mean_event_dff2;
            medianArea(2, ff-sz2) = medianArea2;
            medianRoiArea(2, ff-sz2) = medianRoiArea2;
            medianArea_max(2, ff-sz2) =  medianArea_max2;
            medianDuration(2, ff-sz2) =  medianDuration2;
            medianInterval(2, ff-sz2) = medianInterval2;
            numEvents(2, ff-sz2) = NumEvents2;
            for r = 1:length(Angles)
                Angles2{2} = [Angles2{2}; Angles{r}];
                Diameters2{2} = [Diameters2{2}; Diameters{r}];
                Durations2{2} = [Durations2{2}; Durations{r}];
            end
            median_dom_angles2(2, ff-sz2) = median_dom_angles;
            median_dir_ind2(2, ff-sz2) = median_dir_ind;
        end
        
    end
    
    median_center_dist(median_center_dist == 0) = NaN;
    median_center_stepDist_sum(median_center_stepDist_sum == 0) = NaN;
    median_center_stepSpeed(median_center_stepSpeed == 0) = NaN;
    median_mean_event_dff(median_mean_event_dff == 0) = NaN;
    medianArea(medianArea == 0) = NaN;
    medianRoiArea(medianRoiArea == 0) = NaN;
    medianArea_max(medianArea_max == 0) = NaN;
    medianDuration(medianDuration == 0) = NaN;
    medianInterval(medianInterval == 0) = NaN;
    numEvents(numEvents == 0) = NaN;
    median_dom_angles2(median_dom_angles2 == 0) = NaN;
    median_dir_ind2(median_dir_ind2 == 0) = NaN;

%median_dir_ind2 = abs(median_dir_ind2);    
% 
% idx_tmp = isnan(median_center_dist);
% if sum(sum(idx_tmp)) > 0
%     median_center_dist(idx_tmp) = 0;
% end
% idx_tmp = isnan(median_center_stepDist_sum);
% if sum(sum(idx_tmp)) > 0
%     median_center_stepDist_sum(idx_tmp) = 0;
% end
% idx_tmp = isnan(median_center_stepSpeed);
% if sum(sum(idx_tmp)) > 0
%     median_center_stepSpeed(idx_tmp) = 0;
% end
% idx_tmp = isnan(median_mean_event_dff);
% if sum(sum(idx_tmp)) > 0
%     median_mean_event_dff(idx_tmp) = 0;
% end
% idx_tmp = isnan(medianArea);
% if sum(sum(idx_tmp)) > 0
%     medianArea(idx_tmp) = 0;
% end
% idx_tmp = isnan(medianRoiArea);
% if sum(sum(idx_tmp)) > 0
%     medianRoiArea(idx_tmp) = 0;
% end
% idx_tmp = isnan(medianArea_max);
% if sum(sum(idx_tmp)) > 0
%     medianArea_max(idx_tmp) = 0;
% end
% idx_tmp = isnan(medianDuration);
% if sum(sum(idx_tmp)) > 0
%     medianDuration(idx_tmp) = 0;
% end
% idx_tmp = isnan(medianInterval);
% if sum(sum(idx_tmp)) > 0
%     medianInterval(idx_tmp) = 0;
% end
% idx_tmp = isnan(numEvents);
% if sum(sum(idx_tmp)) > 0
%     numEvents(idx_tmp) = 0;
% end
% idx_tmp = isnan(median_dom_angles2);
% if sum(sum(idx_tmp)) > 0
%     median_dom_angles2(idx_tmp) = 0;
% end
% idx_tmp = isnan(median_dir_ind2);
% if sum(sum(idx_tmp)) > 0
%     median_dir_ind2(idx_tmp) = 0;
% end
%     
    
    %% plots
    cd ..
    cd ..
    cd ..
    %mkdir(['th_', num2str(th)])
    %cd(['th_', num2str(th)])
    %colorVector = [128 128 128; 0 202 100] / 255;
    colorVector = [0 202 100; 128 128 128] / 255;
    %colorVector = [0 202 100; 0 0 0] / 255;
    %colorVector = [0,0,0.6; 0.8,0,0];
    colorVector_line = [105 105 105; 105 105 105] / 255;
    colors = colorVector([1, 2], :);
    
    % area boxplots
    line_cords_x = [1 2];
    c = figure; set(c, 'position', [0 0 1500 450])
    num_animals = size(medianArea, 2);
    
    subplot(1, 3, 1)
    for g = 1:length(tag)
        data(g, :) = medianArea(g, :);
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
        hold on
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
    end
    if same_animal
        for kk = 1:num_animals
            line_cords_y = [data(1,kk) data(2,kk)];
            plot(line_cords_x, line_cords_y, 'k')
        end
    end
    clear tmp
    tmp{1} = data(1,:);
    tmp{2} = data(2,:);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', tag, 'whisker', 1000)
    box off
    
    ylabel('median area (um^2)')
    if same_animal
        [~, p1] = ttest(data(1,:), data(2,:));
    else
            [~, p1] = ttest2(data(1,:), data(2,:));
    end
    title(['p = ', num2str(p1)])
    clear data
    
    subplot(1, 3, 2)
    for g = 1:length(tag)
        data(g, :) = medianArea_max(g, :);
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
        hold on
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
    end
    if same_animal
        for kk = 1:num_animals
            line_cords_y = [data(1,kk) data(2,kk)];
            plot(line_cords_x, line_cords_y, 'k')
        end
    end
    clear tmp
    tmp{1} = data(1,:);
    tmp{2} = data(2,:);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', tag, 'whisker', 1000)
    box off
    
    ylabel('peak area (um^2)')
    if same_animal
        [~, p2] = ttest(data(1,:), data(2,:));
    else
            [~, p2] = ttest2(data(1,:), data(2,:));
    end
    title(['p = ', num2str(p2)])
    clear data
    
    subplot(1, 3, 3)
    for g = 1:length(tag)
        data(g, :) = medianRoiArea(g, :);
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
        hold on
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
    end
    if same_animal
        for kk = 1:num_animals
            line_cords_y = [data(1,kk) data(2,kk)];
            plot(line_cords_x, line_cords_y, 'k')
        end
    end
    clear tmp
    tmp{1} = data(1,:);
    tmp{2} = data(2,:);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', tag, 'whisker', 1000)
    box off
    
    ylabel('median Roi area (um^2)')
    if same_animal
        [~, p11] = ttest(data(1,:), data(2,:));
    else
            [~, p11] = ttest2(data(1,:), data(2,:));
    end
    title(['p = ', num2str(p11)])
    clear data
    saveas(c, ['areas_th', num2str(th), '.png'])
    
    %% speed and distance plots
    line_cords_x = [1 2];
    c = figure; set(c, 'position', [0 0 1500 450])
    num_animals = size(medianArea, 2);
    subplot(1, 3, 1)
    for g = 1:length(tag)
        data(g, :) = median_center_dist(g, :);
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
        hold on
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
    end
    if same_animal
        for kk = 1:num_animals
            line_cords_y = [data(1,kk) data(2,kk)];
            plot(line_cords_x, line_cords_y, 'k')
        end
    end
    clear tmp
    tmp{1} = data(1,:);
    tmp{2} = data(2,:);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', tag, 'whisker', 1000)
    box off
    
    ylabel('distance (mm)')
    if same_animal
        [~, p3] = ttest(data(1,:), data(2,:));
    else
           [~, p3] = ttest2(data(1,:), data(2,:));
    end
    title(['p = ', num2str(p3)])
    clear data
    
    subplot(1, 3, 2)
    for g = 1:length(tag)
        data(g, :) = median_center_stepSpeed(g, :);
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
        hold on
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
    end
    if same_animal
        for kk = 1:num_animals
            line_cords_y = [data(1,kk) data(2,kk)];
            plot(line_cords_x, line_cords_y, 'k')
        end
    end
    clear tmp
    tmp{1} = data(1,:);
    tmp{2} = data(2,:);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', tag, 'whisker', 1000)
    box off
    
    ylabel('step speed (mm/s)')
    if same_animal
        [~, p4] = ttest(data(1,:), data(2,:));
    else
           [~, p4] = ttest2(data(1,:), data(2,:));
    end
    title(['p = ', num2str(p4)])
    clear data
    
    subplot(1, 3, 3)
    for g = 1:length(tag)
        data(g, :) = median_center_stepDist_sum(g, :);
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
        hold on
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
    end
    if same_animal
        for kk = 1:num_animals
            line_cords_y = [data(1,kk) data(2,kk)];
            plot(line_cords_x, line_cords_y, 'k')
        end
    end
    clear tmp
    tmp{1} = data(1,:);
    tmp{2} = data(2,:);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', tag, 'whisker', 1000)
    box off
    
    ylabel('step distance sum (mm)')
    if same_animal
        [~, p5] = ttest(data(1,:), data(2,:));
    else
           [~, p5] = ttest2(data(1,:), data(2,:));
    end
    title(['p = ', num2str(p5)])
    clear data
    
    saveas(c, ['speeds and distances_th', num2str(th), '.png'])
    
    %% durations, intervals, num events plots
    line_cords_x = [1 2];
    c = figure; set(c, 'position', [0 0 1500 450])
    num_animals = size(medianArea, 2);
    
    subplot(1, 3, 1)
    for g = 1:length(tag)
        data(g, :) = medianDuration(g, :);
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
        hold on
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
    end
    if same_animal
        for kk = 1:num_animals
            line_cords_y = [data(1,kk) data(2,kk)];
            plot(line_cords_x, line_cords_y, 'k')
        end
    end
    clear tmp
    tmp{1} = data(1,:);
    tmp{2} = data(2,:);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', tag, 'whisker', 1000)
    box off
    
    ylabel('median duration (s)')
    if same_animal
        [~, p6] = ttest(data(1,:), data(2,:));
    else
           [~, p6] = ttest2(data(1,:), data(2,:));
    end
    title(['p = ', num2str(p6)])
    clear data
    
    subplot(1, 3, 2)
    for g = 1:length(tag)
        data(g, :) = medianInterval(g, :);
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
        hold on
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
    end
    if same_animal
        for kk = 1:num_animals
            line_cords_y = [data(1,kk) data(2,kk)];
            plot(line_cords_x, line_cords_y, 'k')
        end
    end
    clear tmp
    tmp{1} = data(1,:);
    tmp{2} = data(2,:);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', tag, 'whisker', 1000)
    box off
    
    ylabel('interval (s)')
    if same_animal
        [~, p7] = ttest(data(1,:), data(2,:));
    else
            [~, p7] = ttest2(data(1,:), data(2,:));
    end
    title(['p = ', num2str(p7)])
    clear data
    
    subplot(1, 3, 3)
    freqEvents = numEvents./10;
    for g = 1:length(tag)
        data(g, :) = freqEvents(g, :);
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
        hold on
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
    end
    if same_animal
        for kk = 1:num_animals
            line_cords_y = [data(1,kk) data(2,kk)];
            plot(line_cords_x, line_cords_y, 'k')
        end
    end
    clear tmp
    tmp{1} = data(1,:);
    tmp{2} = data(2,:);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', tag, 'whisker', 1000)
    box off
    
    ylabel('number of events/min')
    %ylim([0 max(ylim)])
    if same_animal
        [~, p8] = ttest(data(1,:), data(2,:));
    else
           [~, p8] = ttest2(data(1,:), data(2,:));
    end
    title(['p = ', num2str(p8)])
    clear data
    
    saveas(c, ['durations intervals numevents th', num2str(th), '.png'])
    %% mean event dFF
    line_cords_x = [1 2];
    num_animals = size(medianArea, 2);
    figure;
    for g = 1:length(tag)
        data(g, :) = median_mean_event_dff(g, :);
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
        hold on
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
    end
    if same_animal
        for kk = 1:num_animals
            line_cords_y = [data(1,kk) data(2,kk)];
            plot(line_cords_x, line_cords_y, 'k')
        end
    end
    clear tmp
    tmp{1} = data(1,:);
    tmp{2} = data(2,:);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', tag, 'whisker', 1000)
    box off
    
    ylabel('mean event dFF')
    if same_animal
        [~, p9] = ttest(data(1,:), data(2,:));
    else
            [~, p9] = ttest2(data(1,:), data(2,:));
    end
    title(['p = ', num2str(p9)])
    clear data
    
    saveas(gcf, ['mean event dFF th', num2str(th), '.png'])
    
    %% angles plots
    line_cords_x = [1 2];
    num_animals = size(medianArea, 2);
    subplot(1, 2, 1)
    for g = 1:length(tag)
        data(g, :) = median_dom_angles2(g, :);
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
        hold on
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
    end
    if same_animal
        for kk = 1:num_animals
            line_cords_y = [data(1,kk) data(2,kk)];
            plot(line_cords_x, line_cords_y, 'k')
        end
    end
    clear tmp
    tmp{1} = data(1,:);
    tmp{2} = data(2,:);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', tag, 'whisker', 1000)
    box off
    
    ylabel('median fraction of waves in dominant direction')
    if same_animal
        [~, p10] = ttest(data(1,:), data(2,:));
    else
           [~, p10] = ttest2(data(1,:), data(2,:));
    end
    title(['p = ', num2str(p10)])
    clear data
    
    subplot(1, 2, 2)
    for g = 1:length(tag)
        data(g, :) = median_dir_ind2(g, :);
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
        hold on
        scatter(g*ones(1, num_animals), data(g,:), 'MarkerEdgeColor',  colorVector(g, :)', 'LineWidth', 2)
    end
    if same_animal
        for kk = 1:num_animals
            line_cords_y = [data(1,kk) data(2,kk)];
            plot(line_cords_x, line_cords_y, 'k')
        end
    end
    clear tmp
    tmp{1} = data(1,:);
    tmp{2} = data(2,:);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', tag, 'whisker', 1000)
    box off
    
    ylabel('median wave DSI')
    if same_animal
        [~, p15] = ttest(data(1,:), data(2,:));
    else
           [~, p15] = ttest2(data(1,:), data(2,:));
    end
    title(['p = ', num2str(p15)])
    save('waveDSI', 'data');
    clear data
    
    saveas(gcf, ['wave directionality th', num2str(th), '.png'])
    
    %rose plots
    for g = 1:length(tag)
        figure;
        % polar histogram
        polarhistogram(Angles2{g}, 20, 'normalization', 'probability', 'LineWidth', 1, 'FaceColor', colors(g, :));
        thetaticks(0:45:315);
        rticks([0.05, 0.1, 0.15])
        rticklabels({0.05, 0.1, 0.15})
        ax = gca;
        ax.LineWidth = 1;
        title(tag{g});
        %rlim([0 0.08])
        save_fnm = [tag{g}, '_th', num2str(th), '_rosePlot.png'];
        saveas(gcf, save_fnm)
    end
    
    clear counts centers proportion
    %% diameters histograms
    data_all1 = Diameters2{1};
    data_all2 = Diameters2{2};
    data_all = [data_all1; data_all2];
    figure;
    for g = 1:length(tag)
        data = Diameters2{g};
        median_data = median(data, 'omitnan');
        [counts(g, :), centers(g, :)] = hist(data,  0:1:max(max(data_all)) );
        proportion(g, :) = counts(g, :) / length(data);
        if g == 1
            facevalue = 0.3;
        else
            facevalue = 0.2;
        end
        bar(centers(g, :), proportion(g, :), 'barWidth', 10, 'FaceColor', colorVector(g, :)); hold on
        line(median_data * [1 1], [0 max(ylim)], 'lineStyle', '- -', 'color', colorVector(g, :), 'lineWidth', 3); hold on
        clear data
    end
    xlabel('diameters');
    [~, p12] = ttest2(data_all1, data_all2);
    title(['p = ', num2str(p12)])
    saveas(gcf, ['diameters', num2str(th), '.png'])
    
    clear data_all data_all1 data_all2 data counts centers proportion
    %% durations histograms
    data_all1 = Durations2{1};
    data_all2 = Durations2{2};
    data_all = [data_all1; data_all2];
    figure;
    for g = 1:length(tag)
        data = Durations2{g};
        median_data = median(data, 'omitnan');
        [counts(g, :), centers(g, :)] = hist(data,  0:1:max(max(data_all)) );
        proportion(g, :) = counts(g, :) / length(data);
        if g == 1
            facevalue = 0.3;
        else
            facevalue = 0.2;
        end
        bar(centers(g, :), proportion(g, :), 'barWidth', 0.5, 'FaceColor', colorVector(g, :)); hold on
        line(median_data * [1 1], [0 max(ylim)], 'lineStyle', '- -', 'color', colorVector(g, :), 'lineWidth', 3); hold on
        clear data
    end
    xlabel('durations');
    [~, p13] = ttest2(data_all1, data_all2);
    title(['p = ', num2str(p13)])
    saveas(gcf, ['durations', num2str(th), '.png'])
    
    close all
    clearvars -except tag same_animal animalInfo all_th t
end
