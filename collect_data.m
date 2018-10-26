%% set path
restoredefaultpath; % start from scratch
addpath(genpath('C:\Users\mvdm\Documents\GitHub\vandermeerlab\code-matlab\shared')); % clone from GitHub

%% params
cfg_master = [];
cfg_master.debug = 0;
cfg_master.base_fp = 'C:\data\ros-data';
cfg_master.theta = [5.5 8.5];
cfg_master.lowgamma = [30 45]; % avoid notch filter
cfg_master.highgamma = [55 80];

what = {'theta', 'lowgamma', 'highgamma'}; % freqs to process, note must have freq ranges defined in cfg_master

%% get the list of folders we need to process
cd(cfg_master.base_fp);
fc = FindFiles('*.set');
nSessions = length(fc);

%% initialize variables
ALL_data = [];

%% main loop
for iS = 1:nSessions
    
    fprintf('\n*** Session %d/%d... ***\n', iS, nSessions);
    
    sd = ParseSetFileName(fc{iS}); % get this session's info: subject, condition, session, time
    sid = cat(2, sd.condition, num2str(sd.session)); % session ID
    cd(sd.fd);
    
    egf_list = FindFiles('*.egf*');
    for iF = 1:length(egf_list)
        
        % load this file
        this_egf = egf_list{iF}; [~, fn, fe] = fileparts(this_egf); this_egf = cat(2, fn, fe); fe = fe(2:end);
        cfg_csc = []; cfg_csc.fc = {this_egf};
        csc = LoadCSC_Axona(cfg_csc);
        
        % power spectrum
        [P, F] = pwelch(csc.data, csc.cfg.Fs, csc.cfg.Fs/2, 2^16, csc.cfg.Fs); P = 10*log10(P);
        
        if cfg_master.debug
            figure;
            subplot(221); 
            plot(F, 10*log10(P));
            set(gca, 'FontSize', 18, 'XLim', [0 100], 'XTick', 0:10:10); grid on; box off;
            title(this_egf);
           
        end
        
        % extract power and peak frequency in requested bands
        this_pow = [];
        for iW = 1:length(what)
            
            [this_pf, this_pv, this_pa] = FindPeak(P, F, cfg_master.(what{iW}));
            
            ALL_data.(sd.subject).(sid).(what{iW}).(fe).pf(str2num(sd.time)) = this_pf; % should place in loop somehow
            ALL_data.(sd.subject).(sid).(what{iW}).(fe).pv(str2num(sd.time)) = this_pv;
            ALL_data.(sd.subject).(sid).(what{iW}).(fe).pa(str2num(sd.time)) = this_pa;
            
        end
        
        % filter signal in requested bands
        
        
        % theta-phase to gamma power coupling
        
        
        % low-gamma power to high-gamma power cross-correlation
        
        
    end % of egfs
    
end % of sessions


%
function sd = ParseSetFileName(in)

[fp fn fe] = fileparts(in);

subj_idx = regexp(fp, 'NN\d');
sd.subject = fp(subj_idx:subj_idx + 2);

if isempty(regexp(fp, 'Drug'))
    sd.condition = 'vehicle';
else
    sd.condition = 'drug';
end

digit_idx = regexp(fn, '\d');
sd.session = fn(digit_idx(end - 1)); % assume session number is the second to last digit
sd.time = fn(digit_idx(end)); % assume time is the last digit
sd.fd = fp;

end

%
function [pf, pv, pa] = FindPeak(P, F, freq_range)

keep_idx = F > freq_range(1) & F <= freq_range(2);
this_F = F(keep_idx); this_P = P(keep_idx);
[pv, p_idx] = max(this_P); % value at peak
pf = this_F(p_idx); % frequency at peak in freq range
pa = nanmean(this_P); % average in freq range

end