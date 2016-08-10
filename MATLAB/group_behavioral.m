%% Pick subject folders to analyze
clear
mfile = mfilename;
scriptPath = mfilename('fullpath');
scriptPath = scriptPath(1:end-length(mfilename));
addpath(genpath('C:\Users\Rob\Desktop\Dropbox\EAPSI\Analysis\EAPSI_Analysis_Code\MATLAB'));
% subjPaths  = uipickfiles('FilterSpec',scriptPath);
% subjPaths  = uipickfiles('FilterSpec','H:\Sensorimotor_RT_8-14-2015');
subjPaths = {'C:\Users\Rob\Desktop\Dropbox\EAPSI\Analysis\Data\MEG_Subject28'};


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

rhoPDs = [];
rhodPdDs = [];

%% Loop through subject folders, gather vols
for subjPath = subjPaths 
	
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
	
	subj_fly_heard_onsets = heardFlyBeeps(:,4);
	subj_rep_heard_onsets = heardReplayBeeps(:,4);
	subj_fly_missed_onsets = missedFlyBeeps(:,4);
	subj_rep_missed_onsets = missedReplayBeeps(:,4);
	
	plot(sim_cockpit_radios_nav2_freq_hz)
	mriStartIndex
	xplane_markers = sim_cockpit_radios_nav2_freq_hz;
	xplane_markers = round(xplane_markers,0);
	for sampleInd = 1:length(xplane_markers)
		sample = xplane_markers(sampleInd);
		if sample == 3000 || sample == 2000
			back_s = 1;
			while xplane_markers(sampleInd-back_s) == 0 || ( xplane_markers(sampleInd-back_s) == 1000 && xplane_markers(sampleInd-back_s-1) == 1000 )
					back_s = back_s+1;
					if back_s > 1000 % If beep doesn't occur within 1 second before heard/miss
						error('First beep not found');
					end
			end
			break
		end
	end
	earliest_onset_xplane_ind = sampleInd-back_s;
	earliest_onset_mri_ind = earliest_onset_xplane_ind;
% 	earliest_onset_text_ind = ;
	earliest_onset_text_time = min([subj_fly_heard_onsets;subj_rep_heard_onsets;subj_fly_missed_onsets;subj_rep_missed_onsets]);
	earliest_onset_xplane_time = xplane_time(earliest_onset_xplane_ind)/1000;
	earliest_onset_mri_time = mri_time(first_beep_onset_ind);
	text_aligned_xplane_time = xplane_time
	mri_time(first_beep_onset_ind)-mri_time(mriStartIndex)

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
% p = vartestn([fly_heard_vols;rep_heard_vols],'TestType','LeveneAbsolute')
% [h,p] = ttest2(fly_missed_vols,rep_missed_vols)
% [h,p] = ttest2(fly_heard_vols,rep_heard_vols)

% [h,p,ci,stats] = ttest2(fly_missed_vols,rep_missed_vols, 'vartype', 'unequal')
% [h,p,ci,stats] = ttest2(fly_heard_vols,rep_heard_vols, 'vartype', 'unequal')

% [p,h,stats] = ranksum(fly_missed_vols,rep_missed_vols)
% [p,h,stats] = ranksum(fly_heard_vols,rep_heard_vols)

% stats = mwwtest(fly_missed_vols,rep_missed_vols)
% stats = mwwtest(fly_heard_rts,rep_heard_rts)
%% Compute rts statistics
% p = vartestn(MPG,Model_Year,'TestType','LeveneAbsolute')

% [h,p] = ttest2(fly_missed_rts,rep_missed_rts)
% [h,p] = ttest2(fly_heard_rts,rep_heard_rts)


% [h,p] = ttest2(fly_missed_rts,rep_missed_rts, 'vartype', 'unequal')
% [h,p] = ttest2(fly_heard_rts,rep_heard_rts, 'vartype', 'unequal')

% [p,h] = ranksum(fly_missed_rts,rep_missed_rts)
% [p,h] = ranksum(fly_heard_rts,rep_heard_rts)

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