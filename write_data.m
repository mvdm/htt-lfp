%%
addpath('C:\Users\mvdm\Documents\GitHub\htt-lfp\util');

%%
cfg_master = [];
cfg_master.outfile = 'C:\data\ros-data\out.xlsx';


% subject condition freqband file

this_row = 2;

subj = fieldnames(ALL_data);
for iS = 1:length(subj)
    
    this_subj = subj{iS};
    
    cond = fieldnames(ALL_data.(this_subj));
    for iC = 1:length(cond)
        
        this_cond = cond{iC};
        
        freq = fieldnames(ALL_data.(this_subj).(this_cond));
        for iF = 1:length(freq)
            
           this_freq = freq{iF};
           
           fl = fieldnames(ALL_data.(this_subj).(this_cond).(this_freq));
           for iFl = 1:length(fl)
               
               this_fl = fl{iFl};
               this_data = ALL_data.(this_subj).(this_cond).(this_freq).(this_fl);
               this_data = cat(2, this_data.pf, this_data.pv, this_data.pa); % note this determines labels
               
               this_label = {this_subj, this_cond, this_freq, this_fl};
               
               this_cell = conv2xlscell(this_row, 2);
               xlswrite(cfg_master.outfile, this_label, 1, this_cell);
               this_cell = conv2xlscell(this_row, 2 + length(this_label));
               xlswrite(cfg_master.outfile, this_data, 1, this_cell);
               
               this_row = this_row + 1;
               
           end % of files
           
        end % of freqs
        
    end % of conditions
    
end % of subjects