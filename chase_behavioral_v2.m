%% Pick subject folders to analyze
clear
mfile = mfilename;
scriptPath = mfilename('fullpath');
scriptPath = scriptPath(1:end-length(mfilename));
% subjPaths  = uipickfiles('FilterSpec',scriptPath);
% subjPaths  = uipickfiles('FilterSpec','H:\Sensorimotor_RT_8-14-2015');
% subjPaths = {'C:\Users\Rob\Desktop\Dropbox\EAPSI\Analysis\Behavioral\MEG_Subject28'};
subjPaths = {'C:\Users\Rob\Desktop\Dropbox\EAPSI\Analysis\EAPSI_MRI\MRI_Subject12'};
quickload = 1;
xplane = 1;

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
		blockPath = blockPaths(blockPathInd);
		if blockPath.isdir == 0 || strcmp(blockPath.name,'..') || strcmp(blockPath.name,'.')
			continue
		end
		cd(subjPath)
		cd(blockPath.name);

		%% Xplane
		if xplane == 1 
			if quickload == 0 || isempty(dir('quickload.mat'))
				% Parse xml file for var names
				tree = xml2struct('xplane_udp.xml');
				nVars = length(tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item);
				varNames = cell(1,nVars);
				for itemN = 1:nVars
					varNames{itemN} = tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item{itemN}.Attributes.ref;
				end

				xplane_fn = dir('*data.txt');
				try
					xplane_fn = xplane_fn.name;
					fid= fopen(xplane_fn);
					n = 0;
					tline = fgetl(fid);
					while ischar(tline)
						tline = fgetl(fid);
						n = n+1;
					end
					fclose(fid);

					xplane_data = dlmread(xplane_fn,' ',[0 0 n-100 89]);
					if sum(xplane_data(:,89)) == 0
						xplane_data(:,89) = [];
					end
				catch ME
					display(ME.message);
					display('Couldnt read text file, most likely empty');
					display(n);
				end

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

				save 'quickload.mat'
			else
				load('quickload.mat');
			end
			
			
			%% Create MRI time vector
			findIndices= find(diff(sim_cockpit_radios_nav1_freq_hz(1:1000)))+1;
			mriStartIndex = findIndices(1);
			mri_time = (xplane_time-xplane_time(mriStartIndex))/1000000;
			
	% 		close all
	% 		figure; plot(performance, 'b');
	% 		hold on; plot(interpPerf, 'r.')
	% 		title('Performance over time');
	% 		ylim([1,2])
	% 		figure; plot(distance, 'b');
	% 		hold on; plot(interpDist, 'r.')
	% 		title('Distance over time');
	% 		ylim([0,200])
% 			eeg_fn = dir('*.cog');
% 			eeg_fn = eeg_fn.name;
% 			eeg_data = cog_load(eeg_fn,64,0,3,1);
% 			
% 			[rhoPD, pval] = corr(performance,distance,'type','Spearman')
% 			[rhodPdD, pval] = corr(diff(performance),diff(distance),'type','Spearman')
% 			rhoPDs = [rhoPDs rhoPD];
% 			rhodPdDs = [rhodPdDs rhodPdD];
% 			continue
		end
		
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