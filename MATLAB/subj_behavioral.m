%% Initialize
clear
quickload = 1;
xplane = 1;
interp = 1;
plot_block = 0;
block_stats = 0;
eeg = 0;
mri = 0;
subjPath = 'C:\\Users\\Rob\\Desktop\\Dropbox\\EAPSI\\Analysis\\Data\\MEG_Subject28';
addpath(genpath(subjPath));

%% Handle path
cd(subjPath)
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
	cd(subjPath);
	cd(blockPath.name);
	
	%% Import Xplane block data
	if xplane == 1
		%% Gets xplane data, saves if quickload is disabled or quickload file is absent
		if ~quickload || isempty(dir('quickload.mat'))
			%% Parse xml file for var names
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
				
				xplane_block_data = dlmread(xplane_fn,' ',[0 0 n-100 89]);
				if sum(xplane_block_data(:,89)) == 0
					xplane_block_data(:,89) = [];
				end
			catch ME
				display(ME.message);
				display('Couldnt read text file, most likely empty');
				display(n);
			end
			
			%% Map variable names
			xplane_time = xplane_block_data(:,1);
			for varN = 1:nVars
				varName = strrep(varNames{varN}, '/', '_');
				eval([varName,' = xplane_block_data(:,1+varN);'])
			end
			
			%% Extract distance and performance and markers
			performance = xplane_block_data(:,end); % hard coded into data file
			xplane_markers = sim_cockpit_radios_nav2_freq_hz;
			x1 = sim_flightmodel_position_local_x;
			y1 = sim_flightmodel_position_local_y;
			z1 = sim_flightmodel_position_local_z;
			x2 = sim_multiplayer_position_plane1_x;
			y2 = sim_multiplayer_position_plane1_y;
			z2 = sim_multiplayer_position_plane1_z;
			distance = sqrt( (x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2);
			
			mri_markers = sim_cockpit_radios_nav1_freq_hz;
			
			%% Quicksave
			save('quickload.mat', 'distance', 'performance', 'xplane_block_data', 'xplane_markers', 'xplane_time');
			
		else
			%% Quickload
			load('quickload.mat');
		end
	end
	
	%% Create block MRI time vector
	if mri
	findIndices = find(diff(mri_markers))+1;
	mriStartIndex = findIndices(1);
	mri_time = (xplane_time-xplane_time(mriStartIndex))/1000000;
	end
	
	%% Clean distance and perf data
	if interp
		
		%% Interpolate distance
		nXPPoints = length(performance);
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
		interpDist = inpaint_nans(interpDist, 1);
		
		%% Interpolate performance
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
		
	end
	
	%% Plot block distance and perf
	if plot_block
		close all
		figure; plot(performance, 'b');
		hold on; plot(interpPerf, 'r.')
		title('Performance over time');
		ylim([1,2])
		figure; plot(distance, 'b');
		hold on; plot(interpDist, 'r.')
		title('Distance over time');
		ylim([0,200])
	end
	
	%% Import block EEG
	if eeg
		eeg_block_fn = dir('*.cog');
		eeg_block_fn = eeg_block_fn.name;
		eeg_block_data = cog_load(eeg_block_fn,64,0,3,1);
	end
	
	%% Initialize replay containers, find run files
	text_run_rep_files = dir('*replay*.txt');
	text_rep_res = [];
	text_rep_onsets = [];
	text_rep_heard_onsets = [];
	text_rep_missed_onsets = [];
	text_rep_vols = [];
	text_rep_heard_vols = [];
	text_rep_missed_vols = [];
	text_rep_rts = [];
	text_rep_heard_rts = [];
	text_rep_missed_rts = [];
	
	%% Gather replay beeps
	if ~isempty(text_run_rep_files)
		
		%% Loop through replay run files
		for text_run_rep_file_in = 1:length(text_run_rep_files)
			%% Read in file
			text_run_rep_file = text_run_rep_files(text_run_rep_file_in).name;
			display(['Reading in ', text_run_rep_file]);
			vol_data = read_mixed_csv(text_run_rep_file,' ');
			
			%% Line by line, identify beep onset, result, and volume
			try
				for line_ind = 1:length(vol_data)
					if strcmp(vol_data(line_ind,2),'on')
						text_rep_res = [text_rep_res vol_data(line_ind+2,2)];
						text_rep_vols = [text_rep_vols str2double(vol_data(line_ind,4))];
						text_rep_onsets = [text_rep_onsets str2double(vol_data(line_ind,1))];
						text_rep_rts = [text_rep_rts str2double(vol_data(line_ind+2,1))-str2double(vol_data(line_ind,1))];
					end
				end
			catch ME
				display(ME.message)
				display('Couldnt read last trial, most likely')
			end
		end
		
		%% Identify replay heard vs. missed
		for beep = 1:length(text_rep_res)
			if strcmp(text_rep_res(beep),'heard')
				text_rep_heard_onsets = [text_rep_heard_onsets text_rep_onsets(beep)];
				text_rep_heard_vols = [text_rep_heard_vols text_rep_vols(beep)];
				text_rep_heard_rts = [text_rep_heard_rts text_rep_rts(beep)];
			end
			if strcmp(text_rep_res(beep),'missed')
				text_rep_missed_onsets = [text_rep_missed_onsets text_rep_onsets(beep)];
				text_rep_missed_vols = [text_rep_missed_vols text_rep_vols(beep)];
				text_rep_missed_rts = [text_rep_missed_rts text_rep_rts(beep)];
			end
		end
	end
	
	%% Initialize fly containers, find run files
	text_run_fly_files = dir('*fly.txt*');
	text_fly_res = [];
	text_fly_onsets = [];
	text_fly_heard_onsets = [];
	text_fly_missed_onsets = [];
	text_fly_vols = [];
	text_fly_heard_vols = [];
	text_fly_missed_vols = [];
	text_fly_rts = [];
	text_fly_heard_rts = [];
	text_fly_missed_rts = [];
	
	%% Gather fly beeps
	if ~isempty(text_run_fly_files)
		%% Loop through fly run files
		for text_run_fly_file_in = 1:length(text_run_fly_files)
			text_run_fly_file = text_run_fly_files(text_run_fly_file_in).name;
			display(['Reading in ', text_run_fly_file]);
			vol_data = read_mixed_csv(text_run_fly_file,' ');
			
			%% Line by line, identify beep onset, result, and volume
			try
				for line_ind = 1:length(vol_data)
					if strcmp(vol_data(line_ind,2),'on')
						text_fly_res = [text_fly_res vol_data(line_ind+2,2)];
						text_fly_vols = [text_fly_vols str2double(vol_data(line_ind,4))];
						text_fly_onsets = [text_fly_onsets str2double(vol_data(line_ind,1))];
						text_fly_rts = [text_fly_rts str2double(vol_data(line_ind+2,1))-str2double(vol_data(line_ind,1))];
					end
				end
			catch ME
				display(ME.message)
				display('Couldnt read last trial, most likely')
			end
		end
		
		%% Identify replay heard vs. missed
		if ~isempty(text_fly_res)
			for beep = 1:length(text_fly_res)
				if strcmp(text_fly_res(beep),'heard')
					text_fly_heard_onsets = [text_fly_heard_onsets text_fly_onsets(beep)];
					text_fly_heard_vols = [text_fly_heard_vols text_fly_vols(beep)];
					text_fly_heard_rts = [text_fly_heard_rts text_fly_rts(beep)];
				end
				if strcmp(text_fly_res(beep),'missed')
					text_fly_missed_onsets = [text_fly_missed_onsets text_fly_onsets(beep)];
					text_fly_missed_vols = [text_fly_missed_vols text_fly_vols(beep)];
					text_fly_missed_rts = [text_fly_missed_rts text_fly_rts(beep)];
				end
			end
		end
	end
	
	%% Fill missing xplane_data with NaN from text data
	% 	xplane_block_data
	% 	text_rep_onsets
	% 	text_fly_onsets
	text_fly_onsets_ones = [text_fly_onsets; ones(1,length(text_fly_onsets))];
	text_rep_onsets_zeros = [text_rep_onsets; zeros(1,length(text_rep_onsets))];
	text_onsets = sortrows([text_fly_onsets, text_rep_onsets]');
	text_onsets_ones_zeros = [text_fly_onsets_ones, text_rep_onsets_zeros]';
	sorted = sortrows(text_onsets_ones_zeros);
	
	
	%%
	xplane_time_diff = diff(xplane_time);
% 	plot(xplane_time_diff);
	xplane_srate = 1000000/median(diff(xplane_time));
	xplane_period_inUs = median(diff(xplane_time));
% 	plot(sorted);
	% Loop through text events/beeps
	% Loop through xplane samples (events)
	xplane_time_zMs = (xplane_time-xplane_time(1))/1000;
	trans = [];
	for xplane_marker_ind = 2:length(xplane_markers)
		if xplane_markers(xplane_marker_ind-1) ~= xplane_markers(xplane_marker_ind)
			if xplane_markers(xplane_marker_ind) == 0 % this is either first start flying/replay/
				trans = [trans xplane_markers(xplane_marker_ind-1) xplane_markers(xplane_marker_ind)];
			end
		end

% 		if isFirstChange == 1
% 			[xplane_markers(xplane_marker_ind-1), xplane_markers(xplane_marker_ind)]
% 			isFirstChange = 0;
% 		end
		
% 		xplane_time_zMs(xplane_marker_ind)/1000
	end
	
	% If xplane_sample_time == text_event_time, check if xplane marker
	% is nearby, if not, add in text event and fill in missing data
	
	%% Make beep matrix for the block
	% [1/0(heard/missed) 1/0(fly/replay) vol onset]
	block_beeps =  [ones(1,length(text_fly_heard_vols)), ones(1,length(text_rep_heard_vols)), zeros(1,length(text_fly_missed_vols)), zeros(1,length(text_rep_missed_vols));...
		ones(1,length(text_fly_heard_vols)), zeros(1,length(text_rep_heard_vols)), ones(1,length(text_fly_missed_vols)), zeros(1,length(text_rep_missed_vols));...
		text_fly_heard_vols, text_rep_heard_vols, text_fly_missed_vols, text_rep_missed_vols;...
		text_fly_heard_onsets, text_rep_heard_onsets, text_fly_missed_onsets, text_rep_missed_onsets;...
		text_fly_heard_rts, text_rep_heard_rts, text_fly_missed_rts, text_rep_missed_rts]';
	
	%% Append block beeps to subj_beeps in order
	if ~isempty(block_beeps)
		sorted_beeps = sortrows(block_beeps,4);
		subj_beeps = [subj_beeps; sorted_beeps];
	end
	
	%% Run block stats
	if block_stats
		[rhoPD, pval] = corr(performance,distance,'type','Spearman')
		[rhodPdD, pval] = corr(diff(performance),diff(distance),'type','Spearman')
		rhoPDs = [rhoPDs rhoPD];
		rhodPdDs = [rhodPdDs rhodPdD];
	end
end






