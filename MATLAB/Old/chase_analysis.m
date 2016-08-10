%% Pick run folder
clear
mfile = mfilename;
scriptPath = mfilename('fullpath');
scriptPath = scriptPath(1:end-length(mfilename));
runPaths  = uipickfiles('FilterSpec',scriptPath);

subj_fly_heard_vols = [];
subj_rep_heard_vols = [];
subj_fly_missed_vols = [];
subj_rep_missed_vols = [];

%% Get file structures
for runpath = runPaths
    xplane_dir = [runpath{:},'\'];
    xplane_fn = dir([xplane_dir,'*data.txt']);
    xplane_fn = [xplane_dir, xplane_fn.name];
    eeg_fn = dir([xplane_dir,'*.cog']);
    eeg_fn = eeg_fn.name;
    cd(xplane_dir)

    %% Import xplane and eeg data
    if isempty(dir('quickload.mat'))

        %% Import EEG data
        % If you get an error, go into text file and fix line with error. Seems to
        % be problem with newline character being put in at certain times. I went 
        % in and separated the two lines, then took the average of the last column 
        % value of the before and after lines because it was missing. Also, delete 
        % the last line because it will not have all of the columns
        eeg_data = cog_load(eeg_fn,64,0,3,1);

        %% Xplane
        % Parse xml file for var names
        tree = xml2struct([xplane_dir,'xplane_udp.xml']);
        nVars = length(tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item);
        varNames = cell(1,nVars);
        for itemN = 1:nVars
           varNames{itemN} = tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item{itemN}.Attributes.ref;
        end
        xplane_data = load(xplane_fn);

        % Map variable names, performance hard coded into data file
        xplane_time = xplane_data(:,1);
        for varN = 1:nVars
            varName = strrep(varNames{varN}, '/', '_');
            eval([varName,' = xplane_data(:,1+varN);'])
        end
        performance = xplane_data(:,end);

        % Extract distance
        x1 = sim_flightmodel_position_local_x;
        y1 = sim_flightmodel_position_local_y;
        z1 = sim_flightmodel_position_local_z;
        x2 = sim_multiplayer_position_plane1_x;
        y2 = sim_multiplayer_position_plane1_y;
        z2 = sim_multiplayer_position_plane1_z;
        distance = sqrt( (x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2);

        %% Clean distance and perf data
        nXPPoints = length(performance);

        % Distance
        badDcounter = 0;
        interpDist = distance;
        interpDist(1:15) = NaN;
        for sample = 15:nXPPoints
            if interpDist(sample) > 200
               badDcounter = badDcounter + 1;
               interpDist((sample-15):sample+15) = NaN;
               continue
            end
        end

        % Performance
        badPcounter = 0;
        interpPerf = performance;
        interpPerf(1:15) = NaN;
        for sample = 15:nXPPoints
            if interpPerf(sample) < 1.5 %special performance threshold for rob, generally can't be less than 1
               badPcounter = badPcounter + 1;
               interpPerf((sample-15):sample+15) = NaN;
               continue
            end
        end

        interpPerf = inpaint_nans(interpPerf, 1);
        interpDist = inpaint_nans(interpDist, 1);

        xplane_markers = sim_cockpit_radios_nav2_freq_hz;

%         close all
        figure; plot(performance, 'b'); 
        hold on; plot(interpPerf, 'r.')
        title('Performance over time');
        figure; plot(distance, 'b'); 
        hold on; plot(interpDist, 'r.')
        title('Distance over time');
        clearvars -except xplane_time performance distance xplane_dir ...
            xplane_fn eeg_fn interpPerf interpDist eeg_data xplane_markers runpath...
            subj_fly_heard_vols subj_rep_heard_vols subj_fly_missed_vols subj_rep_missed_vols...
            runPaths
        save 'quickload.mat'
    else
        load('quickload.mat')
        close all
        figure; plot(performance, 'b'); hold on; plot(interpPerf, 'r.')
        figure; plot(distance, 'b'); hold on; plot(interpDist, 'r.')
    end

    %% Import text data
    cond_files = dir([xplane_dir,'*y*.txt']);
    text_beep_res = [];
    text_beep_onsets = [];
    text_heard_beep_onsets = [];
    text_missed_beep_onsets = [];
    text_beep_vols = [];
    text_heard_beep_vols = [];
    text_missed_beep_vols = [];

    for cond_file_in = 1:length(cond_files)
        cond_file = cond_files(cond_file_in).name;
        vol_data = read_mixed_csv(cond_file,' ');

        % Line by line, identify beep onset, result, and volume
        for line_ind = 1:length(vol_data)
            if strcmp(vol_data(line_ind,2),'on')
                text_beep_res = [text_beep_res vol_data(line_ind+2,2)];
                text_beep_vols = [text_beep_vols str2double(vol_data(line_ind,4))];
                text_beep_onsets = [text_beep_onsets str2double(vol_data(line_ind,1))];
            end
        end

    end
    for beep = 1:length(text_beep_res)
        if strcmp(text_beep_res(beep),'heard')
            text_heard_beep_onsets = [text_heard_beep_onsets text_beep_onsets(beep)];
            text_heard_beep_vols = [text_heard_beep_vols text_beep_vols(beep)];
        end
        if strcmp(text_beep_res(beep),'missed')
            text_missed_beep_onsets = [text_missed_beep_onsets text_beep_onsets(beep)];
            text_missed_beep_vols = [text_missed_beep_vols text_beep_vols(beep)];
        end
    end

    %% Extract EEG event info
    % Different event codes for different data streams EEG and XPlane
    eeg_markers = eeg_data(:,end);

    eeg_beep_onsets = [];
    eeg_beep_offsets = [];
    eeg_heard_onsets = [];
    eeg_missed_onsets = [];
    eeg_heard_beep_onsets = [];
    eeg_missed_beep_onsets = [];
    eeg_fly_beep_onsets =[];
    eeg_rep_beep_onsets = [];
    eeg_fly_no_onsets =[];
    eeg_rep_no_onsets = [];

    for s = 2:length(eeg_markers)
        % Heard
        if eeg_markers(s)==768 && eeg_markers(s-1) == 0
            eeg_heard_onsets = [eeg_heard_onsets s]; % Heard indicator

            % Look backwards for beep offset
            back_s = 1;
            while eeg_markers(s-back_s) == 0 
                back_s = back_s+1;
            end
            eeg_beep_offsets = [eeg_beep_offsets s-back_s];

            % Look backwards for beep onset, knowing it was heard
            back_s = 1;
            while eeg_markers(s-back_s) == 0 || ( eeg_markers(s-back_s) == 2304 && eeg_markers(s-back_s-1) == 2304 )
                back_s = back_s+1;
            end
            eeg_heard_beep_onsets = [eeg_heard_beep_onsets s-back_s];
        end

        % Missed
        if eeg_markers(s)==512 && eeg_markers(s-1) == 0
            eeg_missed_onsets = [eeg_missed_onsets s]; % Missed indicator

            % Look backwards for beep offset
            back_s = 1;
            while eeg_markers(s-back_s) == 0 
                back_s = back_s+1;
            end
            eeg_beep_offsets = [eeg_beep_offsets s-back_s];

            % Look backwards for beep onset, knowing it was heard
            back_s = 1;
            while eeg_markers(s-back_s) == 0 || ( eeg_markers(s-back_s) == 2304 && eeg_markers(s-back_s-1) ~= 0 )
                back_s = back_s+1;
            end
            eeg_missed_beep_onsets = [eeg_missed_beep_onsets s-back_s+1];
        end

        % Beep onset regardless of heard or missed
        if eeg_markers(s)==2304 && eeg_markers(s-1) ~= 2304
            eeg_beep_onsets = [eeg_beep_onsets s];
        end
        if eeg_markers(s)==256 && eeg_markers(s-1) ~= 256
            eeg_fly_beep_onsets = [eeg_fly_beep_onsets s];
        end
        if eeg_markers(s)==1280 && eeg_markers(s-1) ~= 1280
            eeg_rep_beep_onsets = [eeg_rep_beep_onsets s];
        end
        if eeg_markers(s)==2560 && eeg_markers(s-1) ~= 2560
            eeg_fly_no_onsets = [eeg_fly_no_onsets s];
        end
        if eeg_markers(s)==12800 && eeg_markers(s-1) ~= 12800
            eeg_rep_no_onsets = [eeg_rep_no_onsets s];
        end
    end
    eeg_times = (1:length(eeg_data))/300;
    close all
    figure; plot(eeg_times,eeg_markers); hold on;
    ax = gca;
%     ax.XAxis.Exponent = 0;
%     ax.XAxis.TickLabelFormat = '%,d';
    % plot(eeg_fly_beep_onsets, eeg_markers(eeg_fly_beep_onsets), 'rx')
    % plot(eeg_rep_beep_onsets, eeg_markers(eeg_rep_beep_onsets), 'kx')
    % plot(eeg_fly_no_onsets, eeg_markers(eeg_fly_no_onsets), 'bx')
    % plot(eeg_rep_no_onsets, eeg_markers(eeg_rep_no_onsets), 'gx')
    % plot(eeg_beep_onsets, eeg_markers(eeg_beep_onsets), 'kx')
    % plot(eeg_beep_offsets, eeg_markers(eeg_beep_offsets), 'kx')
    plot(eeg_times(eeg_heard_beep_onsets), eeg_markers(eeg_heard_beep_onsets), 'k+')
    plot(eeg_times(eeg_heard_onsets), eeg_markers(eeg_heard_onsets), 'k+') % Indicates heard beep is before
    plot(eeg_times(eeg_missed_beep_onsets), eeg_markers(eeg_missed_beep_onsets), 'r*')
    plot(eeg_times(eeg_missed_onsets), eeg_markers(eeg_missed_onsets), 'r*')% Indicates missed beep is before


    %% Extract XPlane event codes
    diff_nav2 = diff(xplane_markers);
    startIndex = find(diff_nav2,1)+1;
    startTime = xplane_time(startIndex);

    xplane_rest_inds = find(xplane_markers==250);
    xplane_beep_offsets_inds = [];
    xplane_heard_beep_offsets_inds = [];
    xplane_missed_beep_offsets_inds = [];
    xplane_heard_onsets_inds = [];
    xplane_missed_onsets_inds = [];
    xplane_heard_beep_onsets_inds = [];
    xplane_missed_beep_onsets_inds = [];
    xplane_1000_onsets_inds = [];
    xplane_rep_beep_onsets_inds = [];

    for s = 1:length(xplane_markers)
        % Heard
        if xplane_markers(s)==3000 && xplane_markers(s-1) == 0 % Indicator always changes to 3000 from 0 to indicate a recent beep was heard
            xplane_heard_onsets_inds = [xplane_heard_onsets_inds s]; % Heard indicator

            % Look backwards from indicator for beep offset, marked by changing from 0
            % (to 1000)
            back_s = 1;
            while xplane_markers(s-back_s) == 0 
                back_s = back_s+1;
            end
            xplane_beep_offsets_inds = [xplane_beep_offsets_inds s-back_s];
            xplane_heard_beep_offsets_inds = [xplane_heard_beep_offsets_inds s-back_s];

            % Look backwards for beep onset, knowing it was heard
            back_s = 1;
            while xplane_markers(s-back_s) == 0 || ( xplane_markers(s-back_s) == 1000 && xplane_markers(s-back_s-1) == 1000 )
                back_s = back_s+1;
            end
            xplane_heard_beep_onsets_inds = [xplane_heard_beep_onsets_inds s-back_s];
        end

        % Missed
        if xplane_markers(s)==2000 && xplane_markers(s-1) == 0
            xplane_missed_onsets_inds = [xplane_missed_onsets_inds s]; % Missed indicator

            % Look backwards for beep offset
            back_s = 1;
            while xplane_markers(s-back_s) == 0 
                back_s = back_s+1;
            end
            xplane_beep_offsets_inds = [xplane_beep_offsets_inds s-back_s]; % carried over from heard beeps
            xplane_missed_beep_offsets_inds = [xplane_missed_beep_offsets_inds s-back_s];

            % Look backwards for beep onset, knowing it was missed
            back_s = 1;
            while xplane_markers(s-back_s) == 0 || ( xplane_markers(s-back_s) == 1000 && xplane_markers(s-back_s-1) == 1000 )
                back_s = back_s+1;
            end
            xplane_missed_beep_onsets_inds = [xplane_missed_beep_onsets_inds s-back_s];
        end

        % Beep onset regardless of heard or missed
        if xplane_markers(s)==1000 && xplane_markers(s-1) ~= 1000
            xplane_1000_onsets_inds = [xplane_1000_onsets_inds s]; % mixed up heard onsets with fly with beeps indicator
        end
        if xplane_markers(s)==500 && xplane_markers(s-1) ~= 500
            xplane_rep_beep_onsets_inds = [xplane_rep_beep_onsets_inds s];
        end
    end

    xplane_beep_onsets_inds = sort([xplane_heard_beep_onsets_inds, xplane_missed_beep_onsets_inds]);
    xplane_fly_with_beep_cond_onsets_inds = setdiff(xplane_1000_onsets_inds, xplane_beep_onsets_inds);

    close all
    figure; plot(xplane_markers); hold on; 
    plot(xplane_beep_onsets_inds, xplane_markers(xplane_beep_onsets_inds), 'kx')
		%%
		xplane_beep_onsets = round(xplane_time(xplane_beep_onsets_inds)/1000);
%     plot(xplane_beep_offsets, xplane_markers(xplane_beep_offsets), 'kx')
%     plot(xplane_heard_beep_onsets, xplane_markers(xplane_heard_beep_onsets), 'kx')
%     plot(xplane_heard_onsets, xplane_markers(xplane_heard_onsets), 'k+')
%     plot(xplane_missed_beep_onsets, xplane_markers(xplane_missed_beep_onsets), 'rx')
%     plot(xplane_missed_onsets, xplane_markers(xplane_missed_onsets), 'r+')
%     plot(xplane_rest_inds, 250, 'r*')
%     plot(xplane_fly_with_beep_cond_onsets, xplane_markers(xplane_fly_with_beep_cond_onsets), 'g*')

    %%
    % Don't need the text time, except to make sure it isn't too diff from more
    % precise xplane_time
    % EEG.data = eeg_data(1:64)';
    % EEG.times = xplane_time/1000;
    % z_xplane_time = xplane_time-xplane_time(1);
    % 
    % diff([text_beep_onsets; round(ms_time(xplane_beep_onsets))'],1)

    % the function dsearchn also works
    % xplane_beep_ind = dsearchn(frex',freqIwant);

    %%REACTION TIME
    %% find large diff vols
    large_log = (abs(diff(text_beep_vols))>0.15);
    large_ind = 1:length(text_beep_vols);
    large_ind = [large_ind(large_log), large_ind(large_log)];
    text_beep_vols_fixed = text_beep_vols;
    text_beep_vols_fixed(large_ind)=0;

    close all
    figure;
    plot(text_beep_vols); hold on;
    plot(abs(diff(text_beep_vols)));
    plot(text_beep_vols_fixed);


    %% find which of the xplane_heard_beep fall in fly vs replay
    cond_tab = sortrows([[xplane_fly_with_beep_cond_onsets_inds; ones(1,length(xplane_fly_with_beep_cond_onsets_inds))]' ;[xplane_rep_beep_onsets_inds; zeros(1,length(xplane_rep_beep_onsets_inds))]']); % FLY = 1, REP = 0

    fly_beep_inds = [];
    rep_beep_inds = [];
    for row = 2:length(cond_tab)
        for heard_beep_onset_ind = 1:length(xplane_heard_beep_onsets_inds)
            if xplane_heard_beep_onsets_inds(heard_beep_onset_ind) > cond_tab(row-1,1) && xplane_heard_beep_onsets_inds(heard_beep_onset_ind) < cond_tab(row,1)
                if cond_tab(row-1,2) %Fly
                   fly_beep_inds = [fly_beep_inds heard_beep_onset_ind]; 
                else % Replay
                   rep_beep_inds = [rep_beep_inds heard_beep_onset_ind]; 
                end
            end
            if xplane_heard_beep_onsets_inds(heard_beep_onset_ind) > cond_tab(end,1) && row == length(cond_tab)
                if cond_tab(end,2) %Fly
                   fly_beep_inds = [fly_beep_inds heard_beep_onset_ind]; 
                else % Replay
                   rep_beep_inds = [rep_beep_inds heard_beep_onset_ind]; 
                end
            end
        end
    end

    fly_beep_inds = setdiff(fly_beep_inds,large_ind);
    rep_beep_inds = setdiff(rep_beep_inds,large_ind);
    fly_heard_vols = text_heard_beep_vols(fly_beep_inds);
    rep_heard_vols = text_heard_beep_vols(rep_beep_inds);
    subj_fly_heard_vols = [subj_fly_heard_vols fly_heard_vols];
    subj_rep_heard_vols = [subj_rep_heard_vols rep_heard_vols];

    % figure;
    % plot(xplane_markers); hold on;
    % plot(xplane_heard_beep_onsets,xplane_markers(xplane_heard_beep_onsets),'gx');
    % plot(xplane_heard_beep_onsets(fly_beep_inds),xplane_markers(xplane_heard_beep_onsets(fly_beep_inds)),'bo');
    % plot(xplane_heard_beep_onsets(rep_beep_inds),xplane_markers(xplane_heard_beep_onsets(rep_beep_inds)),'ro');

    % std_fly_heard_vols = std(fly_heard_vols)
    % std_rep_heard_vols = std(rep_heard_vols)
    % mean_fly_heard_vols = mean(fly_heard_vols)
    % mean_rep_heard_vols = mean(rep_heard_vols)
    % figure; hist(fly_heard_vols)
    % figure; hist(rep_heard_vols)
    % [h,p] = ttest2(fly_heard_vols,rep_heard_vols, 'vartype', 'unequal')
    % [h,p] = ttest2(fly_heard_vols,rep_heard_vols)
    % [p,h] = ranksum(fly_heard_vols,rep_heard_vols)

    % out_heard_fly_vols = fly_heard_vols((mean_fly_heard_vols-2*std_fly_heard_vols)<fly_heard_vols<(mean_fly_heard_vols+2*std_fly_heard_vols));
    % out_heard_rep_vols = rep_heard_vols((mean_rep_heard_vols-2*std_rep_heard_vols)<rep_heard_vols<(mean_rep_heard_vols+2*std_rep_heard_vols));
    % std_out_fly_heard_vols = std(out_heard_fly_vols)
    % std_out_rep_heard_vols = std(out_heard_rep_vols)
    % mean_out_fly_heard_vols = mean(out_heard_fly_vols)
    % mean_out_rep_heard_vols = mean(out_heard_rep_vols)
    % figure; hist(out_heard_fly_vols)
    % figure; hist(out_heard_rep_vols)
    % [h,p] = ttest2(out_heard_fly_vols,out_heard_rep_vols, 'vartype', 'unequal')
    % [h,p] = ttest2(out_heard_fly_vols,out_heard_rep_vols)


    %% find which of the xplane_missed_beep fall in fly vs replay
    cond_tab = sortrows([[xplane_fly_with_beep_cond_onsets_inds; ones(1,length(xplane_fly_with_beep_cond_onsets_inds))]'... % FLY = 1
    ;[xplane_rep_beep_onsets_inds; zeros(1,length(xplane_rep_beep_onsets_inds))]']); % REP = 0

    fly_beep_inds = [];
    rep_beep_inds = [];
    for row = 2:length(cond_tab)
        for missed_beep_onset_ind = 1:length(xplane_missed_beep_onsets_inds)
            if xplane_missed_beep_onsets_inds(missed_beep_onset_ind) > cond_tab(row-1,1) && xplane_missed_beep_onsets_inds(missed_beep_onset_ind) < cond_tab(row,1)
                if cond_tab(row-1,2) %Fly
                   fly_beep_inds = [fly_beep_inds missed_beep_onset_ind]; 
                else % Replay
                   rep_beep_inds = [rep_beep_inds missed_beep_onset_ind]; 
                end
            end
            if xplane_missed_beep_onsets_inds(missed_beep_onset_ind) > cond_tab(end,1) && row == length(cond_tab)
                if cond_tab(end,2) %Fly
                   fly_beep_inds = [fly_beep_inds missed_beep_onset_ind]; 
                else % Replay
                   rep_beep_inds = [rep_beep_inds missed_beep_onset_ind]; 
                end
            end
        end
    end

    fly_beep_inds = setdiff(fly_beep_inds,large_ind);
    rep_beep_inds = setdiff(rep_beep_inds,large_ind);
    fly_missed_vols = text_missed_beep_vols(fly_beep_inds);
    rep_missed_vols = text_missed_beep_vols(rep_beep_inds);
    subj_fly_missed_vols = [subj_fly_missed_vols fly_missed_vols];
    subj_rep_missed_vols = [subj_rep_missed_vols rep_missed_vols];

    % figure;
    % plot(xplane_markers); hold on;
    % plot(xplane_missed_beep_onsets,xplane_markers(xplane_missed_beep_onsets),'gx');
    % plot(xplane_missed_beep_onsets(fly_beep_inds),xplane_markers(xplane_missed_beep_onsets(fly_beep_inds)),'bo');
    % plot(xplane_missed_beep_onsets(rep_beep_inds),xplane_markers(xplane_missed_beep_onsets(rep_beep_inds)),'ro');

    % std_fly_missed_vols = std(fly_missed_vols)
    % std_rep_missed_vols = std(rep_missed_vols)
    % mean_fly_missed_vols = mean(fly_missed_vols)
    % mean_rep_missed_vols = mean(rep_missed_vols)
    % figure; hist(fly_missed_vols)
    % figure; hist(rep_missed_vols)
    % [h,p] = ttest2(fly_missed_vols,rep_missed_vols, 'vartype', 'unequal')
    % [h,p] = ttest2(fly_missed_vols,rep_missed_vols)
    % [p,h] = ranksum(fly_missed_vols,rep_missed_vols)

    % out_missed_fly_vols = fly_missed_vols((mean_fly_missed_vols-2*std_fly_missed_vols)<fly_missed_vols<(mean_fly_missed_vols+2*std_fly_missed_vols));
    % out_missed_rep_vols = rep_missed_vols((mean_rep_missed_vols-2*std_rep_missed_vols)<rep_missed_vols<(mean_rep_missed_vols+2*std_rep_missed_vols));
    % std_out_fly_missed_vols = std(out_missed_fly_vols)
    % std_out_rep_missed_vols = std(out_missed_rep_vols)
    % mean_out_fly_missed_vols = mean(out_missed_fly_vols)
    % mean_out_rep_missed_vols = mean(out_missed_rep_vols)
    % figure; hist(out_missed_fly_vols)
    % figure; hist(out_missed_rep_vols)
    % [h,p] = ttest2(out_missed_fly_vols,out_missed_rep_vols, 'vartype', 'unequal')
    % [h,p] = ttest2(out_missed_fly_vols,out_missed_rep_vols)

end
%%
% close all;
% figure; plot(interpPerf,interpDist);
% [rho,p] = corr(interpPerf,interpDist, 'type', 'Spearman')
% figure; plot(abs(diff(interpPerf)),interpDist(1:end-1));
% [rho,p] = corr(abs(diff(interpPerf)),interpDist(1:end-1), 'type', 'Spearman')

%%
% heard_perf = interpPerf(xplane_heard_beep_onsets);
% missed_perf = interpPerf(xplane_missed_beep_onsets);
% 
% heard_dist = interpDist(xplane_heard_beep_onsets);
% missed_dist = interpDist(xplane_missed_beep_onsets);
% 
% [p,h] = ranksum(heard_perf,missed_perf)
% mean(heard_perf)
% mean(missed_perf)
% 
% [p,h] = ranksum(heard_dist,missed_dist)
% mean(heard_dist)
% mean(missed_dist)

%%
% figure;
% plot(interpPerf); hold on;
% plot(xplane_heard_beep_onsets,interpPerf(xplane_heard_beep_onsets),'k+')
% plot(xplane_missed_beep_onsets,interpPerf(xplane_missed_beep_onsets),'rx')

%%
% hist(subj_fly_heard_vols)
% hist(subj_fly_heard_vols)
figure; hold on;
plot(subj_fly_heard_vols);
plot(subj_rep_heard_vols);
plot(subj_fly_missed_vols);
plot(subj_rep_missed_vols);
legend({'Heard while Flying', 'Heard while Passive', 'Missed while Flying', 'Missed while Passive'}, 'Position', [0.6 0.2 0.1 0.1]);
legend boxoff

%%
% close all
% Fit
% [muhat,sigmahat] = normfit(x);
% % Bin
% xmin = -5; % lower bound
% xmax = 5; % upper bound
% h = 0.5; % bin width
% edges = xmin:h:xmax;
% % Histogram
% n = histc(x,edges);
% figure;
% bar(edges,n,'histc');
% Overlay the distribution
% hold;
% f = @(z) N*h*normpdf(z,muhat,sigmahat);
% fplot(f,[xmin xmax],'color','r');

bins = 20:0.1:30;
nbins = length(bins);
figure; hold on;

subplot(4,1,1);
N = hist(subj_fly_heard_vols,bins);
hist(subj_fly_heard_vols,bins);
line([mean(subj_fly_heard_vols) mean(subj_fly_heard_vols)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
line([median(subj_fly_heard_vols) median(subj_fly_heard_vols)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
title('Heard while Flying')
xlim([22,26])
ylim([0,max(N)])

subplot(4,1,2);
N = hist(subj_rep_heard_vols,bins);
hist(subj_rep_heard_vols,bins);
line([mean(subj_rep_heard_vols) mean(subj_rep_heard_vols)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
line([median(subj_rep_heard_vols) median(subj_rep_heard_vols)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
title('Heard while Passive')
xlim([22,26])
ylim([0,max(N)])

subplot(4,1,3);
N = hist(subj_fly_missed_vols,bins);
hist(subj_fly_missed_vols,bins)
line([mean(subj_fly_missed_vols) mean(subj_fly_missed_vols)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
line([median(subj_fly_missed_vols) median(subj_fly_missed_vols)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
title('Missed while Flying')
xlim([22,26])
ylim([0,max(N)])

subplot(4,1,4);
N = hist(subj_rep_missed_vols,bins);
hist(subj_rep_missed_vols,bins);
l1 = line([mean(subj_rep_missed_vols) mean(subj_rep_missed_vols)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
l2 = line([median(subj_rep_missed_vols) median(subj_rep_missed_vols)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
title('Missed while Passive')
xlim([22,26])
ylim([0,max(N)])
legend([l1,l2],{'Mean','Median'},'Position',[0.8 0.1675 0.025 0.025])
% legend boxoff
% [h,p] = ttest2(subj_fly_missed_vols,subj_rep_missed_vols, 'vartype', 'unequal')
% [h,p] = ttest2(subj_fly_missed_vols,subj_rep_missed_vols)
[p,h] = ranksum(subj_fly_missed_vols,subj_rep_missed_vols)

% [h,p] = ttest2(subj_fly_heard_vols,subj_rep_heard_vols, 'vartype', 'unequal')
% [h,p] = ttest2(subj_fly_heard_vols,subj_rep_heard_vols)
[p,h] = ranksum(subj_fly_heard_vols,subj_rep_heard_vols)