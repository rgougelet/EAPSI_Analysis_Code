function [out] = cog_load(fname, EEG_CHS, EXT_CHS, ACC_CHS, notch); 

%System parameters to convert raw ADC units to physical units
VREF_EEG = 2.5;
GAIN_EEG = 3;
SCALE_EEG = 2^32;
ISTIM = 24e-9;
EEG_TO_VOLTS = 2*VREF_EEG/(GAIN_EEG*SCALE_EEG);

%Paramters of Aux Box v1
VREF_EXT = 4.75;
GAIN_EXT = 0.5;
SCALE_EXT = 2^32;
EXT_TO_VOLTS = VREF_EXT/GAIN_EXT/SCALE_EXT;

%Accelerometer based on ADXL327 at 2.5V supply
VREF_ACC = 2.5;
SCALE_ACC = 2^24;
ACC_G_PER_V = 0.362;

%compute total number of channels in data file
%add two additional for packet counter and trigger
CHS_TOTAL = EEG_CHS+EXT_CHS+ACC_CHS+2;

%open file
f1 = fopen(fname);
out = fread(f1,'int32');
fclose(f1);
%reshape raw data into matrix
out = out(1:(floor(length(out)/CHS_TOTAL))*CHS_TOTAL);
out = reshape(out, CHS_TOTAL, length(out)/CHS_TOTAL)';

%strip impedance check data on EEG channels if requested
if(notch==1)
    %fs/4 notch
    B = [0.85 0 0.85];
    A = [1 0 0.7];
    out(:,1:EEG_CHS) = filter(B,A, out(:,1:EEG_CHS));

    %fs/2 notch
    B = [0.8 0.8];
    A = [1 0.6];
    out(:,1:EEG_CHS) = filter(B,A, out(:,1:EEG_CHS));
end

%scale EEG channels into volts
out(:,1:EEG_CHS) = EEG_TO_VOLTS*out(:,1:EEG_CHS);

%scale EXT channels if present
if(EXT_CHS>0)
    out(:,EEG_CHS+1:EEG_CHS+EXT_CHS) = EXT_TO_VOLTS*out(:,EEG_CHS+1:EEG_CHS+EXT_CHS);
end

%scale ACC channels if present, converts to units of g (9.8 m/s^2)
if(ACC_CHS>0)
    out(:,EEG_CHS+EXT_CHS+1:EEG_CHS+EXT_CHS+ACC_CHS) = VREF_ACC*(out(:,EEG_CHS+EXT_CHS+1:EEG_CHS+EXT_CHS+ACC_CHS)/SCALE_ACC-0.5)/ACC_G_PER_V;
end

%No scaling for Packet Counter or Trigger necessary
end
