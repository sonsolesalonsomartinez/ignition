%% calculate_dynamics
%
% - Loads all the timeseries saved in data_timeseries.mat
% - Bandpass filters the time series
% - Calculates ignition based meassures of integration, hierarchy and
%   metastability; local meassures of ignition and ignition variability;
%   and phase-based global metastability.
%
% Input:  Time series data
%             - data_ts is a matrix of subjects x regions x volumes
% Output: Local meassures
%             - ignition  
%             - ignition_variability 
%         Global meassures
%             - global_integration 
%             - spatial_diversity 
%             - temporal_variability 
%             - phase_metastability
%             - phase_metastability_surro
%
% Code from Gustavo Deco
% Eddited by Sonsoles Alonso Martinez
% alonsomartinez.sonsoles@gmail.com
% Last edited Feb 2019
%

clearvars

% set path
mypath = '/your/path/to/data/';
% load time series
load([mypath, 'data_timeseries.mat'],'data_ts');
 
% Basic parameters  
TR= 2; % sampling interval(TR)  
[n_Subjects, N_areas, Tmax] =size(data_ts);
Isubdiag = find(tril(ones(N_areas),-1));
T = 1:Tmax; 
nTRs = 3; % nTRs size window
n_surrogates =10000;

% Filter parameters
flp = .02;                    % lowpass frequency of filter
fhi = .1;                     % highpass
delt = TR;                    % sampling interval
k = 9;                        % 9th order butterworth filter
fnq = 1/(2*delt);             % Nyquist frequency
Wn = [flp/fnq fhi/fnq];       % butterworth bandpass non-dimensional frequency
[bfilt,afilt] = butter(k,Wn); % construct the filter  

% Realocate variables
ignition = cell(n_Subjects,1);              % local ignition
ignition_variability = cell(n_Subjects,1);  % local ignition variability  
global_integration = cell(n_Subjects,1);    % ignition-based global integration
spatial_diversity = cell(n_Subjects,1);     % ignition-based global hierarchy
temporal_variability = cell(n_Subjects,1);  % ignition-based global metastability
phase_metastability = cell(n_Subjects,1);   % phase-based metastability
phase_metastability_surro = cell(n_surrogates,n_Subjects);% phase-based metastability surrogates


% Get meassures
disp(['Getting dynamics for ' num2str(n_Subjects) ' subjects'])

for nsub = 1:n_Subjects 
    
    BOLD = squeeze(data_ts(nsub,:,:));  %obtain regions x volumes matrix; 
    
    % --------------------- events & phases ----------------------------
    phases=zeros(N_areas,Tmax);
    events=zeros(N_areas,Tmax);   
    for seed = 1:N_areas
      % demean time series
      ts = demean(detrend(BOLD(seed,:)));
      % filter time series
      signal_filt = filtfilt(bfilt,afilt,ts); 
      % get phases
      phases(seed,:) = angle(hilbert(signal_filt));
      % get events
      tise = detrend(demean(signal_filt(T)));
      ev1 = tise>std(tise)+mean(tise);
      ev2 = [0 ev1(1:end-1)];
      events(seed,:) = (ev1-ev2)>0;          
    end

    % ----------------------- integration --------------------------------
    integ=zeros(1,T(end));
    for t = T
        % For each time point a matrix of co-occurring events is calculated
        events_matrix=zeros(N_areas,N_areas);
        for i = 1:N_areas 
            for j = 1:N_areas 
              events_matrix(i,j) = events(i,t)*events(j,t); 
            end
        end 
        A = events_matrix; 
        A = A-eye(N_areas);     
        [~, csize]=get_components(abs(A));% number of regions within comp.
        cs = max(csize); % size of the largest component
        integ(t)=cs/N_areas;  
    end 

    % --------------------- ignition-driven integration -----------------
    nevents2 = zeros(N_areas,1);
    IntegStim2 = zeros(N_areas,N_areas,1);
    % save events and integration values for nTRs after the event
    for seed = 1:N_areas
          flag = 0; % flag 0 to start counting the events
          for t = T
                %detect first event (nevents = matrix with 1xnode and
                %number of events in each cell)
                if events(seed,t) == 1 && flag == 0 
                  flag = 1;
                  nevents2(seed) = nevents2(seed)+1;  
                end
                %save integration value for nTRs after the first event
                %(nodesx(nTR-1)xevents)
                if flag > 0
                  IntegStim2(seed,flag,nevents2(seed)) = integ(t);
                  flag = flag+1;
                end
                %after nTRs, set flag to 0 and wait for the next event
                %(then, integ saved for (nTRs-1) events)
                if flag == nTRs
                  flag = 0;
                end
          end
    end

    %---------------------------------------------------------------
    %                       LOCAL MEASSURES
    %---------------------------------------------------------------
    % mean and std of the max ignition within nTRs for each node
    mevokedinteg2 = zeros(N_areas,1);
    stdevokedinteg2 = zeros(N_areas,1);
    for seed = 1:N_areas
      % nodal ignition
      mevokedinteg2(seed) = mean(max(squeeze(IntegStim2(seed,:,1:nevents2(seed))))); 
      % nodal ignition variability
      stdevokedinteg2(seed) = std(max(squeeze(IntegStim2(seed,:,1:nevents2(seed))))); 
    end
    ignition{nsub} = mevokedinteg2;  
    ignition_variability{nsub} = stdevokedinteg2;

    %---------------------------------------------------------------
    %                       GLOBAL MEASSURES
    %---------------------------------------------------------------    
    % ignition-based integration
    global_integration{nsub} = mean(mevokedinteg2);   
    % ignition-based hierarchy
    spatial_diversity{nsub} = std(mevokedinteg2); 
    % ignition-based metastability
    temporal_variability{nsub} = mean(stdevokedinteg2);

    
    %------------------- Phase-based Metastability  ------------------
    pattern=zeros(T(end),length(Isubdiag));
    syncdata=zeros(1,T(end));
    for t=T
        patt=zeros(N_areas,N_areas-1);
        for i=1:N_areas
            for j=1:i-1
             patt(i,j)=cos(adif(phases(i,t),phases(j,t)));
            end
        end 
        pattern(t,:)=patt(Isubdiag);
        % kuramoto order parameter
        kudata=sum(complex(cos(phases(:,t)),sin(phases(:,t))))/N_areas;
        syncdata(t)=abs(kudata); 
    end
    phase_metastability{nsub} = std(syncdata);

    %--------------- Phase-based Metastability surrogates --------------
    for ns=1:n_surrogates
        for i=1:N_areas
            phases_s(i,:)=phases(i,randperm(Tmax));
        end
        for t=T
            kudata=sum(complex(cos(phases_s(:,t)),sin(phases_s(:,t))))/N_areas;
            syncdata(t)=abs(kudata);
        end
        phase_metastability_surro{ns,nsub} = std(syncdata);
    end


end
save([mypath, 'results.mat'],'ignition', 'ignition_variability',...
                   'global_integration', 'spatial_diversity',...
                   'temporal_variability','phase_metastability',...
                   'Tmax','T','N_areas','n_Subjects','phase_metastability_surro');
