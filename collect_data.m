%% set path
restoredefaultpath; % start from scratch
addpath(genpath('C:\Users\mvdm\Documents\GitHub\vandermeerlab\code-matlab\shared')); % clone from GitHub
addpath(genpath('C:\Users\mvdm\Documents\GitHub\htt-lfp'));
%% params
cfg_master = [];
cfg_master.debug = 0;
cfg_master.base_fp = 'C:\data\ros-data';
cfg_master.theta = [5.5 8.5];
cfg_master.lowgamma = [30 45]; % avoid notch filter
cfg_master.highgamma = [55 80];
cfg_master.pbins = -pi:pi/18:pi; % phase bins for cross-freq coupling

what = {'theta', 'lowgamma', 'highgamma'}; % freqs to process, note must have freq ranges defined in cfg_master

cfg_filter = [];
cfg_filter.type = 'fdesign';
cfg_filter.order = 10;
cfg_filter.f = cfg_master.theta;

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
        for iW = 1:length(what) % loop over freq bands
            
            [this_pf, this_pv, this_pa] = FindPeak(P, F, cfg_master.(what{iW}));
            
            ALL_data.(sd.subject).(sid).(what{iW}).(fe).pf(str2num(sd.time)) = this_pf; % should place in loop somehow
            ALL_data.(sd.subject).(sid).(what{iW}).(fe).pv(str2num(sd.time)) = this_pv;
            ALL_data.(sd.subject).(sid).(what{iW}).(fe).pa(str2num(sd.time)) = this_pa;
            
        end % of freq band loop
        
        % get theta phase
        cfg_filter.f = cfg_master.theta;
        csc_theta = FilterLFP(cfg_filter, csc);
        csc_theta.data = angle(hilbert(csc_theta.data));
        
        % gamma powers
        cfg_filter.f = cfg_master.lowgamma; csc_lowgamma = FilterLFP(cfg_filter, csc);
        csc_lowgamma.data = abs(hilbert(csc_lowgamma.data));
        
        cfg_filter.f = cfg_master.highgamma; csc_highgamma = FilterLFP(cfg_filter, csc);
        csc_highgamma.data = abs(hilbert(csc_highgamma.data));
        
        % theta-phase to gamma power coupling
        cfg_tc = []; cfg_tc.bins = cfg_master.pbins; cfg_tc.interp = 'nearest';
        [~, t_lg, ~] = MakeTC_1D(cfg_tc, csc_theta.tvec, csc_theta.data, csc_lowgamma.tvec, csc_lowgamma.data);
        [~, t_hg, ~] = MakeTC_1D(cfg_tc, csc_theta.tvec, csc_theta.data, csc_highgamma.tvec, csc_highgamma.data);
        
        % low-gamma power to high-gamma power cross-correlation
        [lg_x_hg, xbin] = xcorr(csc_lowgamma.data, csc_highgamma.data, 300, 'coeff');
        
        if cfg_master.debug
            subplot(222); clear h;
            h(1) = plot(cfg_tc.bins(1:end-1), t_lg, 'b', 'LineWidth', 1); hold on;
            h(2) = plot(cfg_tc.bins(1:end-1), t_hg, 'r', 'LineWidth', 1);
            set(gca, 'FontSize', 18, 'XTick', -pi:pi/2:pi);
            legend(h, {'lg', 'hg'}); legend boxoff;
           
            subplot(223);
            plot(xbin, lg_x_hg); set(gca, 'FontSize', 18, 'XTick', -300:150:300);
            
            drawnow; pause; close all
        end
        
        % store summary stats in big struct
        ALL_data.(sd.subject).(sid).cfc.(fe).tlg(str2num(sd.time)) = (max(t_lg) - min(t_lg)) ./ (max(t_lg) + min(t_lg)); % cfc magnitude
        ALL_data.(sd.subject).(sid).cfc.(fe).thg(str2num(sd.time)) = (max(t_hg) - min(t_hg)) ./ (max(t_hg) + min(t_hg)); % cfc magnitude
        
        [~, max_p] = max(t_lg);
        ALL_data.(sd.subject).(sid).cfc.(fe).tlgp(str2num(sd.time)) = cfg_master.pbins(max_p);
        [~, max_p] = max(t_hg);
        ALL_data.(sd.subject).(sid).cfc.(fe).thgp(str2num(sd.time)) = cfg_master.pbins(max_p);
        
        
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