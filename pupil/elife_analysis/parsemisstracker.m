function missing_data=parsemisstracker(sub_data)


% this function is passed the output from tracker chop and calcuates the
% missing data for each trial

%sub_data.tracker_data=out;
sample_rate=2;
baseline_size=1000;
% size of time window collected (in ms) after each event
assessment_size=6000;
time_window=length(-baseline_size:sample_rate:assessment_size-1);

trial_number=size(sub_data.tracker_data.rew_outcome,1);
data=sub_data.tracker_data.data;
%initialise data matrices
options_missa=zeros(trial_number,time_window);
rew_outcome_missa=zeros(trial_number,time_window);
pun_outcome_missa=zeros(trial_number,time_window);
options_missb=zeros(trial_number,time_window);
rew_outcome_missb=zeros(trial_number,time_window);
pun_outcome_missb=zeros(trial_number,time_window);
last_time_point=data(end,1);

% missing data
l_missing=sub_data.tracker_data.removed_data(:,1) | sub_data.tracker_data.removed_data(:,3) ; %missing or blinked
r_missing=sub_data.tracker_data.removed_data(:,2) | sub_data.tracker_data.removed_data(:,4);
both_missing=l_missing & r_missing;
any_missing=l_missing | r_missing;


options_time=sub_data.tracker_data.timing.options;
rew_outcome_time=sub_data.tracker_data.timing.rew_outcome;
pun_outcome_time=sub_data.tracker_data.timing.pun_outcome;

for i=1:trial_number
    % first get the data for the options period
    opt_start_time=options_time(i,1)-baseline_size;
    opt_stop_time=options_time(i,1)+assessment_size;
    
    % if tracker data starts after trial onset or ends before trial finish
    % then disregard
    if opt_stop_time<last_time_point && opt_start_time>data(1,1)
        options_missa(i,:)=any_missing(data(:,1)> opt_start_time & data(:,1)<=opt_stop_time,1);
        options_missb(i,:)=both_missing(data(:,1)> opt_start_time & data(:,1)<=opt_stop_time,1);
    else
        options_missa(i,:)=1;
        options_missb(i,:)=1;
    end
    
    %get data for reward outcome period
    rew_outcome_start_time=rew_outcome_time(i,1)-baseline_size;
    rew_outcome_stop_time=rew_outcome_time(i,1)+assessment_size;
    
    if rew_outcome_stop_time<last_time_point && rew_outcome_start_time>data(1,1)
        rew_outcome_missa(i,:)=any_missing(data(:,1)> rew_outcome_start_time & data(:,1)<=rew_outcome_stop_time,1);
        rew_outcome_missb(i,:)=both_missing(data(:,1)> rew_outcome_start_time & data(:,1)<=rew_outcome_stop_time,1);
    else
        rew_outcome_missa(i,:)=1;
        rew_outcome_missb(i,:)=1;
    end
    % get data from punishment outcome period
    pun_outcome_start_time=pun_outcome_time(i,1)-baseline_size;
    pun_outcome_stop_time=pun_outcome_time(i,1)+assessment_size;
    
    
    if pun_outcome_stop_time<last_time_point && pun_outcome_start_time>data(1,1)
        pun_outcome_missa(i,:)=any_missing(data(:,1)> pun_outcome_start_time & data(:,1)<=pun_outcome_stop_time,1);
        pun_outcome_missb(i,:)=both_missing(data(:,1)> pun_outcome_start_time & data(:,1)<=pun_outcome_stop_time,1);
    else
        pun_outcome_missa(i,:)=1;
        pun_outcome_missb(i,:)=1;
    end
    
end

missing_data=struct;

missing_data.options_missa=options_missa;
missing_data.options_missb=options_missb;
missing_data.rew_outcome_missa=rew_outcome_missa;
missing_data.rew_outcome_missb=rew_outcome_missb;
missing_data.pun_outcome_missa=pun_outcome_missa;
missing_data.pun_outcome_missb=pun_outcome_missb;
