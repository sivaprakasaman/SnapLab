% This is a full transparent template script that can be used to quickly
% playout many trials of a given stimulus at a chosen level and throw
% synchrounous (within 1 sample) triggers for the EEG.
% Note that the stimulus is re-written to the TDT buffer in each trial.
% This allows for stims for each trial to be distinct if needed.

% Copyright 2016-2019 Hari Bharadwaj. All rights reserved.
% hbharadwaj@purdue.edu

clear all; close all hidden; clc; %#ok<CLALL>

try
    % Loading random generator seed and state so that anything generated
    % randomly will be regenerated with the same realization everytime
    load('s.mat');
    rng(s);
    
    % Initialize play circuit
    fig_num=99;
    USB_ch=1;
    IAC = -1;
    FS_tag = 3;
    Fs = 48828.125;
    [f1RZ,RZ,FS]=load_play_circuit(FS_tag,fig_num,USB_ch,0,IAC);
    
    
    % Experiment parameters
    ntrials = 5; % This is per polarity. We typically want about 3-500.
    levels = [60, 70]; % Could be a list the way code is written
    isi = 0.2; % Average interstimulus interval (only approx guarenteed)
    
    % Some cushion for time taken to write each trial to buffer
    % Stim has to be short enough to be written in this time, typically 50
    % ms should be enough for most stims.
    bufferwritetime = 0.05; 
    jitter = 0.1; % Maximum value of uniformly distributed jitter in ISI
    
    % Send trigger to EEG to start saving
    invoke(RZ, 'SetTagVal', 'trigval',253);
    invoke(RZ, 'SoftTrg', 6);
    
    % Wait at least 3 seconds before presenting any stims (e.g., for
    % filter transients, etc.)
    WaitSecs(3.0);
    
    % Load or generate a stimulus with the intent of playing it in both
    % polarities. Replace with other stims here.
    % Must be mono (i.e., size Nsamples-by-1) with a variable called 'y'.
    % This can also be generated within the main loop if each trial has a
    % different stimulus, but depending on stim generation code, in some
    % cases that might slowdown the overall timing (but not affect the
    % synchrony of triggers).
    % The sampling rate has to Fs = 48828.125 Hz based on earlier hardcoded
    % properties of this template script.
    
    
    % Make a 1s-long, 4000 Hz tone modulated at 223 Hz, just as an example
    t = 0:(1/Fs):(1 - 1/Fs);
    t = t(:); % Note that we want size of Nsample-by-1 rather
    risetime = 0.025;
    carrier = sin(2*pi*4000*t);
    env = 0.5 + 0.5*sin(2*pi*223*t);
    y = scaleSound(rampsound(carrier .* env,Fs,risetime));
    
    % Make a different variable for each polarity
    yplus = y;
    yminus = -1 * y;
    
    
    % Using jitter to make sure that background noise averages out across
    % trials. We use jitter that is random between 0 and 'jitter' seconds.
    % Average duration added by the jitter is jitter/2 seconds
    jitlist = rand(ntrials, numel(levels), 2)*jitter;
    isi = isi - jitter/2;
    
    if isi < 0.05
        error('Interstimulus interval too short, cannot continue');
    end
    
    
    % Get the duration of the stimulus
    dur = numel(y)/Fs;
    
    % Keep track of time if needed
    tstart = tic;
    
    
    for j = 1:ntrials
        for L = 1:numel(levels)
            
            for p = [1, 2]
                if(p==1)
                    y = yplus;
                else
                    y = yminus;
                end
                
                stimrms = rms(y);
                
                
                
                chanL = y;
                chanR = y;
                
                stimTrigger = L + (p-1)*numel(levels);
                
                
                jit = jitlist(j, L, 1);
                
                stimlength = numel(y); % Recalculate, just in case
                
                %---------------------------------------------------------
                % Stimulus calibration calculations based on known
                % hardware and known .rcx circuit properties
                %---------------------------------------------------------
                % ER-2s give 100dB SPL for a 1kHz tone with a 1V-rms drive.
                % BasicPlay.rcx converts +/- 1 in MATLAB to +/- 5V at TDT
                % output, so you get 114 dB SPL for MATLAB rms of 1.
                %
                % So if we want a signal with rms of "stimrms" to be X dB,
                % then you have to attenuate in harwdware by below amount:
                
                dropL = 114 - levels(L) + db(stimrms);
                dropR = 114 - levels(L) + db(stimrms); % Assumes diotic
                
                
                invoke(RZ, 'SetTagVal', 'trigval', stimTrigger);
                invoke(RZ, 'SetTagVal', 'nsamps', stimlength);
                %write to buffer left ear
                invoke(RZ, 'WriteTagVEX', 'datainL', 0, 'F32', chanL);
                %write to buffer right ear
                invoke(RZ, 'WriteTagVEX', 'datainR', 0, 'F32', chanR);
                
                %setting analog attenuation L
                invoke(RZ, 'SetTagVal', 'attA', dropL);
                %setting analog attenuation R
                invoke(RZ, 'SetTagVal', 'attB', dropR);
                
                % Just giving time for data to be written into buffer
                WaitSecs(bufferwritetime);
                
                %Start playing from the buffer:
                invoke(RZ, 'SoftTrg', 1); %Playback trigger
                fprintf(1,' Trial Number %d/%d\n', j, ntrials);
                WaitSecs(dur + isi + jit);
            end
            
        end
    end
    toc(tstart); % Just to help get a sense of how long things really take
    
    %Clearing I/O memory buffers:
    invoke(RZ,'ZeroTag','datainL');
    invoke(RZ,'ZeroTag','datainR');
    
    % Have at least 3 seconds of no stim EEG data at the end
    WaitSecs(3.0);
    
    % Send trigger to EEG computer asking it to stop saving
    invoke(RZ, 'SetTagVal', 'trigval', 254);
    invoke(RZ, 'SoftTrg', 6);
    
    close_play_circuit(f1RZ,RZ);
    fprintf(1,'\n Done with data collection!\n');
    
catch me
    close_play_circuit(f1RZ,RZ);
    rethrow(me);
end
