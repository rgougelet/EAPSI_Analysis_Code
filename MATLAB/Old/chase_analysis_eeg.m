xplane_dir = [cd,'\'];
xplane_fn = dir([xplane_dir,'*data.txt']);
xplane_fp = [xplane_dir, xplane_fn.name];
eeg_fn = dir([xplane_dir,'*.cog']);
eeg_fn = eeg_fn.name;

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG = eeg_checkset( EEG );
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
eeg_data = cog_load(eeg_fn,64,0,3,1);
eeg_data = eeg_data(:,[1:64,end])';
EEG = pop_importdata('dataformat','array','nbchan',0,'data','eeg_data','setname',xplane_fn.name(1:end-4),'srate',300,'pnts',0,'xmin',0);% EEG.data = eeg_data(:,[1:64,end])';
EEG = eeg_checkset( EEG );
EEG = pop_chanevent(EEG, EEG.nbchan,'edge','leading','edgelen',0);
EEG = eeg_checkset( EEG );
EEG = pop_chanedit(EEG, 'load',{'C:\\Users\\Rob\\Desktop\\Dropbox\\EAPSI\\Analysis\\Behavioral\\cognionics_channels64.ced' 'filetype' 'autodetect'});
EEG = eeg_checkset( EEG );
EEG = pop_eegfiltnew(EEG, 1, 50, 990, 0, [], 0);
EEG = eeg_checkset( EEG );
EEG = pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
EEG = eeg_checkset( EEG );
EEG = pop_resample( EEG, 32);
EEG.setname = xplane_fn.name(1:end-4);
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',[EEG.setname,'.set'],'filepath',xplane_dir);

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

% EEG = clean_rawdata(EEG, 5, [-1], 0.8, 4, 20, 0.5);
% EEG = eeg_checkset( EEG );
    
% EEG = pop_epoch( EEG, {  '2304'  }, [-0.25           1], 'newname', [EEG.setname,'_E'], 'epochinfo', 'yes');
eeglab redraw