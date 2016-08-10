%% Pick subject folders to analyze
clear
mfile = mfilename;
scriptPath = mfilename('fullpath');
scriptPath = scriptPath(1:end-length(mfilename));
% subjPaths  = uipickfiles('FilterSpec',scriptPath);
% subjPaths  = uipickfiles('FilterSpec','H:\Sensorimotor_RT_8-14-2015');
subjPaths = {'C:\Users\Rob\Desktop\Dropbox\EAPSI\Analysis\Behavioral\MEG_Subject28'};

%% Initialize containers for group
fly_heard_vols = [];
rep_heard_vols = [];
fly_missed_vols = [];
rep_missed_vols = [];

fly_heard_rts = [];
rep_heard_rts = [];
fly_missed_rts = [];
rep_missed_rts = [];

%% Plotting options
plot_subj_rts = 0;
plot_subj_vols = 0;
plot_subj_hists = 0;

plot_rts = 0;
plot_vols = 0;
plot_rts_hists = 0;
plot_vols_hists = 0;

%% Loop through subject folders, gather vols
for subjPath = subjPaths 
	
	%% Handle path
	subjPath = subjPath{:};
	cd(subjPath)
	subjPath
	blockPaths = dir;
	
	%% Initialize subject containers for each subject
	subj_beeps = [];
	subj_fly_heard_vols = [];
	subj_rep_heard_vols = [];
	subj_fly_missed_vols = [];
	subj_rep_missed_vols = [];

	subj_fly_heard_onsets = [];
	subj_rep_heard_onsets = [];
	subj_fly_missed_onsets = [];
	subj_rep_missed_onsets = [];
	
	subj_fly_heard_rts = [];
	subj_rep_heard_rts = [];
	subj_fly_missed_rts = [];
	subj_rep_missed_rts = [];
	
	%% Loop through block folders in subject folder
	for blockPathInd = 1:length(blockPaths) 
		%% Make sure it's a folder not a file
		runPath = blockPaths(blockPathInd);
		if runPath.isdir == 0 || strcmp(runPath.name,'..') || strcmp(runPath.name,'.')
			continue
		end
		cd(subjPath)
		cd(runPath.name);

		%% Import replay text run data
		replay_cond_files = dir('*replay*.txt');
		replay_text_res = [];
		replay_text_onsets = [];
		replay_text_heard_onsets = [];
		replay_text_missed_onsets = [];
		replay_text_vols = [];
		replay_text_heard_vols = [];
		replay_text_missed_vols = [];
		replay_text_rts = [];
		replay_text_heard_rts = [];
		replay_text_missed_rts = [];
		
		
		if ~isempty(replay_cond_files)
			for cond_file_in = 1:length(replay_cond_files)
					cond_file = replay_cond_files(cond_file_in).name;
					display(['Reading in ', cond_file]);
					
					vol_data = read_mixed_csv(cond_file,' ');
					
					% Line by line, identify beep onset, result, and volume
					try
						for line_ind = 1:length(vol_data)
							if strcmp(vol_data(line_ind,2),'on')
								replay_text_res = [replay_text_res vol_data(line_ind+2,2)];
								replay_text_vols = [replay_text_vols str2double(vol_data(line_ind,4))];
								replay_text_onsets = [replay_text_onsets str2double(vol_data(line_ind,1))];
								replay_text_rts = [replay_text_rts str2double(vol_data(line_ind+2,1))-str2double(vol_data(line_ind,1))];
							end
						end
					catch ME
						display(ME.message)
						display('Couldnt read last trial, most likely')
					end
			end
			for beep = 1:length(replay_text_res)
				if strcmp(replay_text_res(beep),'heard')
					replay_text_heard_onsets = [replay_text_heard_onsets replay_text_onsets(beep)];
					replay_text_heard_vols = [replay_text_heard_vols replay_text_vols(beep)];
					replay_text_heard_rts = [replay_text_heard_rts replay_text_rts(beep)];
				end
				if strcmp(replay_text_res(beep),'missed')
					replay_text_missed_onsets = [replay_text_missed_onsets replay_text_onsets(beep)];
					replay_text_missed_vols = [replay_text_missed_vols replay_text_vols(beep)];
					replay_text_missed_rts = [replay_text_missed_rts replay_text_rts(beep)];
				end
			end
		end
		
		%% Import fly text run data
		fly_cond_files = dir('*fly.txt*');
		fly_text_res = [];
		fly_text_onsets = [];
		fly_text_heard_onsets = [];
		fly_text_missed_onsets = [];
		fly_text_vols = [];
		fly_text_heard_vols = [];
		fly_text_missed_vols = [];
		fly_text_rts = [];
		fly_text_heard_rts = [];
		fly_text_missed_rts = [];
		
		if ~isempty(fly_cond_files)
			for cond_file_in = 1:length(fly_cond_files)
					cond_file = fly_cond_files(cond_file_in).name;
					display(['Reading in ', cond_file]);
					vol_data = read_mixed_csv(cond_file,' ');
					
					% Line by line, identify beep onset, result, and volume
					try
						for line_ind = 1:length(vol_data)
							if strcmp(vol_data(line_ind,2),'on')
								fly_text_res = [fly_text_res vol_data(line_ind+2,2)];
								fly_text_vols = [fly_text_vols str2double(vol_data(line_ind,4))];
								fly_text_onsets = [fly_text_onsets str2double(vol_data(line_ind,1))];
								fly_text_rts = [fly_text_rts str2double(vol_data(line_ind+2,1))-str2double(vol_data(line_ind,1))];
							end
						end
					catch ME
						display(ME.message)
						display('Couldnt read last trial, most likely')
					end
			end
			if ~isempty(fly_text_res)
				for beep = 1:length(fly_text_res)
					if strcmp(fly_text_res(beep),'heard')
						fly_text_heard_onsets = [fly_text_heard_onsets fly_text_onsets(beep)];
						fly_text_heard_vols = [fly_text_heard_vols fly_text_vols(beep)];
						fly_text_heard_rts = [fly_text_heard_rts fly_text_rts(beep)];
					end
					if strcmp(fly_text_res(beep),'missed')
						fly_text_missed_onsets = [fly_text_missed_onsets fly_text_onsets(beep)];
						fly_text_missed_vols = [fly_text_missed_vols fly_text_vols(beep)];
						fly_text_missed_rts = [fly_text_missed_rts fly_text_rts(beep)];
					end
				end
			end
		end
		
		%% Make beep matrix for the block
		% [1/0(heard/missed) 1/0(fly/replay) vol onset]
		run_beeps =  [ones(1,length(fly_text_heard_vols)), ones(1,length(replay_text_heard_vols)), zeros(1,length(fly_text_missed_vols)), zeros(1,length(replay_text_missed_vols));...
			ones(1,length(fly_text_heard_vols)), zeros(1,length(replay_text_heard_vols)), ones(1,length(fly_text_missed_vols)), zeros(1,length(replay_text_missed_vols));...
			fly_text_heard_vols, replay_text_heard_vols, fly_text_missed_vols, replay_text_missed_vols;...
			fly_text_heard_onsets, replay_text_heard_onsets, fly_text_missed_onsets, replay_text_missed_onsets;...
			fly_text_heard_rts, replay_text_heard_rts, fly_text_missed_rts, replay_text_missed_rts]';
		
		%% Append run beeps to subj_beeps in order
		if ~isempty(run_beeps)
			sorted_beeps = sortrows(run_beeps,4);
			subj_beeps = [subj_beeps; sorted_beeps];
		end
		
	end
	
	%% Catch empty subjects
	if isempty(subj_beeps)
		display('Subject had no behavioral data!')
		continue
	end
	
	%% For the subject, remove large beepdiffs, accounting for different conditions
	keepInds = zeros(length(subj_beeps),1);
	for beep = 2:length(subj_beeps)
		beepdiff = abs(subj_beeps(beep-1,3)-subj_beeps(beep,3));
		conddiff = abs(subj_beeps(beep-1,2)-subj_beeps(beep,2));
		if beepdiff < 0.11
			keepInds(beep-1) = 1;
		else
			if conddiff == 1
				keepInds(beep-1) = 1;
			else
				keepInds(beep-1) = 0;
			end
		end
	end
% 	subj_beeps = [subj_beeps, keepInds];
	
	%% Sort vols and onsets into conditions
	keepBeeps = subj_beeps(keepInds==1,:);
	heardBeeps = keepBeeps(keepBeeps(:,1)==1,:);
	missedBeeps = keepBeeps(keepBeeps(:,1)==0,:);
	
	heardFlyBeeps = heardBeeps(heardBeeps(:,2)==1,:);
	heardReplayBeeps = heardBeeps(heardBeeps(:,2)==0,:);
	missedFlyBeeps = missedBeeps(missedBeeps(:,2)==1,:);
	missedReplayBeeps = missedBeeps(missedBeeps(:,2)==0,:);
	
	subj_fly_heard_vols = heardFlyBeeps(:,3);
	subj_rep_heard_vols = heardReplayBeeps(:,3);
	subj_fly_missed_vols = missedFlyBeeps(:,3);
	subj_rep_missed_vols = missedReplayBeeps(:,3);
	
	subj_fly_heard_rts = heardFlyBeeps(:,5);
	subj_rep_heard_rts = heardReplayBeeps(:,5);
	subj_fly_missed_rts = missedFlyBeeps(:,5);
	subj_rep_missed_rts = missedReplayBeeps(:,5);
	
% 	subj_fly_heard_onsets = heardFlyBeeps(:,4);
% 	subj_rep_heard_onsets = heardReplayBeeps(:,4);
% 	subj_fly_missed_onsets = missedFlyBeeps(:,4);
% 	subj_rep_missed_onsets = missedReplayBeeps(:,4);

	%% Save subject vols for group stats
	fly_heard_vols = [fly_heard_vols; subj_fly_heard_vols];
	rep_heard_vols = [rep_heard_vols; subj_rep_heard_vols];
	fly_missed_vols = [fly_missed_vols; subj_fly_missed_vols];
	rep_missed_vols = [rep_missed_vols; subj_rep_missed_vols];
	
	fly_heard_rts = [fly_heard_rts; subj_fly_heard_rts];
	rep_heard_rts = [rep_heard_rts; subj_rep_heard_rts];
	fly_missed_rts = [fly_missed_rts; subj_fly_missed_rts];
	rep_missed_rts = [rep_missed_rts; subj_rep_missed_rts];

	%% Plotting subj beep vols over time
	if plot_subj_vols == 1
		close all

		%% Plot heard vs missed beep vols
		figure; hold on;
		for beep = 1:length(keepBeeps)
				if keepBeeps(beep,1) == 1
						plot(beep,keepBeeps(beep,3),'r*')
				else
						plot(beep,keepBeeps(beep,3),'b*')
				end
		end
		title('Heard vs. Missed');

		%% Plot fly vs replay beep vols
		figure; hold on;
		for beep = 1:length(keepBeeps)
				if keepBeeps(beep,2) == 1
						plot(beep,keepBeeps(beep,3),'g*')
				else
						plot(beep,keepBeeps(beep,3),'k*')
				end
		end
		title('Fly vs. Replay');
	end
	
	%% Plot subj histograms of vols by condition
	if plot_subj_hists == 1
		bins = 20:0.1:60;
		nbins = length(bins);
		close all
		set(0,'DefaultFigureWindowStyle','docked')
		figure; hold on;
		ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
			1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

		text(0.5, 1,subjPath(end-12:end),'HorizontalAlignment'...
			,'center','VerticalAlignment', 'top')

		subplot(4,1,1);
		N = hist(subj_fly_heard_vols,bins);
		hist(subj_fly_heard_vols,bins);
		line([mean(subj_fly_heard_vols) mean(subj_fly_heard_vols)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
		line([median(subj_fly_heard_vols) median(subj_fly_heard_vols)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
		title('Heard while Flying')
		% xlim([22,26])
		% ylim([0,max(N)])

		subplot(4,1,2);
		N = hist(subj_rep_heard_vols,bins);
		hist(subj_rep_heard_vols,bins);
		line([mean(subj_rep_heard_vols) mean(subj_rep_heard_vols)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
		line([median(subj_rep_heard_vols) median(subj_rep_heard_vols)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
		title('Heard while Passive')
		% xlim([22,26])
		% ylim([0,max(N)])

		subplot(4,1,3);
		N = hist(subj_fly_missed_vols,bins);
		hist(subj_fly_missed_vols,bins)
		line([mean(subj_fly_missed_vols) mean(subj_fly_missed_vols)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
		line([median(subj_fly_missed_vols) median(subj_fly_missed_vols)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
		title('Missed while Flying')
		% xlim([22,26])
		% ylim([0,max(N)])

		subplot(4,1,4);
		N = hist(subj_rep_missed_vols,bins);
		hist(subj_rep_missed_vols,bins);
		l1 = line([mean(subj_rep_missed_vols) mean(subj_rep_missed_vols)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
		l2 = line([median(subj_rep_missed_vols) median(subj_rep_missed_vols)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
		title('Missed while Passive')
		% xlim([22,26])
		% ylim([0,max(N)])
		legend([l1,l2],{'Mean','Median'},'Position',[0.8 0.1675 0.025 0.025])
	end
	
	%% Plot histograms of rts by condition 
	if plot_subj_rts == 1
		bins = 0:20:1000;
		nbins = length(bins);
		close all
		set(0,'DefaultFigureWindowStyle','docked')
		figure; hold on;
		ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
			1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

		text(0.5, 1,subjPath(end-12:end),'HorizontalAlignment'...
			,'center','VerticalAlignment', 'top')

		subplot(4,1,1);
		N = hist(subj_fly_heard_rts,bins);
		hist(subj_fly_heard_rts,bins);
		line([mean(subj_fly_heard_rts) mean(subj_fly_heard_rts)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
		line([median(subj_fly_heard_rts) median(subj_fly_heard_rts)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
		title('Heard while Flying')
		% xlim([22,26])
		% ylim([0,max(N)])

		subplot(4,1,2);
		N = hist(subj_rep_heard_rts,bins);
		hist(subj_rep_heard_rts,bins);
		line([mean(subj_rep_heard_rts) mean(subj_rep_heard_rts)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
		line([median(subj_rep_heard_rts) median(subj_rep_heard_rts)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
		title('Heard while Passive')
		% xlim([22,26])
		% ylim([0,max(N)])

		subplot(4,1,3);
		N = hist(subj_fly_missed_rts,bins);
		hist(subj_fly_missed_rts,bins)
		line([mean(subj_fly_missed_rts) mean(subj_fly_missed_rts)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
		line([median(subj_fly_missed_rts) median(subj_fly_missed_rts)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
		title('Missed while Flying')
		% xlim([22,26])
		% ylim([0,max(N)])

		subplot(4,1,4);
		N = hist(subj_rep_missed_rts,bins);
		hist(subj_rep_missed_rts,bins);
		l1 = line([mean(subj_rep_missed_rts) mean(subj_rep_missed_rts)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
		l2 = line([median(subj_rep_missed_rts) median(subj_rep_missed_rts)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
		title('Missed while Passive')
		% xlim([22,26])
		% ylim([0,max(N)])
		legend([l1,l2],{'Mean','Median'},'Position',[0.8 0.1675 0.025 0.025])
	end
	
end

%% Plot group vol histograms
if plot_vols_hists == 1
	close all

	bins = 20:0.1:60;
	nbins = length(bins);
	figure; hold on;

	subplot(4,1,1);
	N = hist(fly_heard_vols,bins);
	hist(fly_heard_vols,bins);
	line([mean(fly_heard_vols) mean(fly_heard_vols)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
	line([median(fly_heard_vols) median(fly_heard_vols)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
	title('Heard while Flying Vol')
	xlim([20,45])
	% ylim([0,max(N)])

	subplot(4,1,2);
	N = hist(rep_heard_vols,bins);
	hist(rep_heard_vols,bins);
	line([mean(rep_heard_vols) mean(rep_heard_vols)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
	line([median(rep_heard_vols) median(rep_heard_vols)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
	title('Heard while Passive Vol')
	xlim([20,45])
	% ylim([0,max(N)])

	subplot(4,1,3);
	N = hist(fly_missed_vols,bins);
	hist(fly_missed_vols,bins)
	line([mean(fly_missed_vols) mean(fly_missed_vols)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
	line([median(fly_missed_vols) median(fly_missed_vols)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
	title('Missed while Flying Vol')
	xlim([20,45])
	% ylim([0,max(N)])

	subplot(4,1,4);
	N = hist(rep_missed_vols,bins);
	hist(rep_missed_vols,bins);
	l1 = line([mean(rep_missed_vols) mean(rep_missed_vols)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
	l2 = line([median(rep_missed_vols) median(rep_missed_vols)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
	title('Missed while Passive Vol')
	xlim([20,45])
	% ylim([0,max(N)])
	legend([l1,l2],{'Mean','Median'},'Position',[0.8 0.1675 0.025 0.025])
	% legend boxoff
end

%% Plot group rts histograms
if plot_rts_hists == 1
	close all

	bins = 0:5:2000;
	nbins = length(bins);
	figure; hold on;

	subplot(2,1,1);
	N = hist(fly_heard_rts,bins);
	hist(fly_heard_rts,bins);
	line([mean(fly_heard_rts) mean(fly_heard_rts)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
	line([median(fly_heard_rts) median(fly_heard_rts)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
	title('Heard while Flying RT')
	xlim([0,1000])
	ylim([0,max(N)])

	subplot(2,1,2);
	N = hist(rep_heard_rts,bins);
	hist(rep_heard_rts,bins);
	l1 = line([mean(rep_heard_rts) mean(rep_heard_rts)], [0 max(N)], 'Color', 'r', 'Linewidth', 3);
	l2 = line([median(rep_heard_rts) median(rep_heard_rts)], [0 max(N)], 'Color', 'g', 'Linewidth', 3);
	title('Heard while Passive RT')
	xlim([0,1000])
	ylim([0,max(N)])
	
	legend([l1,l2],{'Mean','Median'},'Position',[0.8 0.1675 0.025 0.025])
	% legend boxoff
	
end

%% Compute vol statistics
% [h,p] = ttest2(fly_missed_vols,rep_missed_vols)
% [h,p] = ttest2(fly_heard_vols,rep_heard_vols)

[h,p] = ttest2(fly_missed_vols,rep_missed_vols, 'vartype', 'unequal')
[h,p] = ttest2(fly_heard_vols,rep_heard_vols, 'vartype', 'unequal')

[p,h] = ranksum(fly_missed_vols,rep_missed_vols)
[p,h] = ranksum(fly_heard_vols,rep_heard_vols)

%% Compute rts statistics
% [h,p] = ttest2(fly_missed_rts,rep_missed_rts)
% [h,p] = ttest2(fly_heard_rts,rep_heard_rts)

% [h,p] = ttest2(fly_missed_rts,rep_missed_rts, 'vartype', 'unequal')
% [h,p] = ttest2(fly_heard_rts,rep_heard_rts, 'vartype', 'unequal')

[p,h] = ranksum(fly_missed_rts,rep_missed_rts)
[p,h] = ranksum(fly_heard_rts,rep_heard_rts)

%% Trash
%{
	
    subj_fly_heard_vols = [subj_fly_heard_vols fly_text_heard_vols];
    subj_rep_heard_vols = [subj_rep_heard_vols replay_text_heard_vols];
    subj_fly_missed_vols = [subj_fly_missed_vols fly_text_missed_vols];
    subj_rep_missed_vols = [subj_rep_missed_vols replay_text_missed_vols];

    subj_fly_heard_onsets = [subj_fly_heard_onsets fly_text_heard_onsets];
    subj_rep_heard_onsets = [subj_rep_heard_onsets replay_text_heard_onsets];
    subj_fly_missed_onsets = [subj_fly_missed_onsets fly_text_missed_onsets];
    subj_rep_missed_onsets = [subj_rep_missed_onsets replay_text_missed_onsets];

    run_beeps =  [ones(1,length(subj_fly_heard_vols)), ones(1,length(subj_rep_heard_vols)), zeros(1,length(subj_fly_missed_vols)), zeros(1,length(subj_rep_missed_vols));...
    ones(1,length(subj_fly_heard_vols)), zeros(1,length(subj_rep_heard_vols)), ones(1,length(subj_fly_missed_vols)), zeros(1,length(subj_rep_missed_vols));...
    subj_fly_heard_vols, subj_rep_heard_vols, subj_fly_missed_vols, subj_rep_missed_vols;...
     subj_fly_heard_onsets, subj_rep_heard_onsets, subj_fly_missed_onsets, subj_rep_missed_onsets]';
%}