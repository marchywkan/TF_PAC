%% This program computes MP based dynamic PAC for a simulated data
warning off;
addpath("eeglab2019_0");
eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% load the data
Fs=2000;
EEG = pop_importdata('dataformat','matlab','nbchan',0,'data','restEEG_81_84.mat','srate',Fs,'pnts',0,'xmin',0);

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
% EEG.data=syn_data2;

% set the normalized reconstruction energy threshold and number of minimum MP atoms. 
nEnergy = 0.98; 
min_atoms = 30;

%% Extract MP atoms and features

% Edited this to capture higher frequecies
[BOOK , LASTCOM] = pop_mp_calc(EEG , 1 , 1 , 40 , 20000 , 0.01 , nEnergy , min_atoms , 2000 , 1, 1);
atomall = real(squeeze(BOOK.reconstruction(1,1,:,:)));

atom_freq= round(BOOK.parameters.frequencies);
atom_widths = BOOK.parameters.widths;
atmtag=size(atom_widths,1);
latency=zeros(1,atmtag);
tstrt=zeros(1,atmtag);

for i=1:atmtag
    [~ , ind_max1] = max(real(BOOK.reconstruction(1,1,i,:)));
    latency(i) = ind_max1 / Fs;
    ind_poc1 = find(real(BOOK.reconstruction(1,1,i,:)) > 0.2 , 1);
    if isempty(ind_poc1)
    ind_poc1=0;
    end
    tstrt(i) = ind_poc1 / Fs;
end
zindx = find(atom_freq<=1)';
atom_widths(zindx)=[];
latency(zindx)=[];
tstrt(zindx)=[];
atom_freq1=atom_freq;
atom_freq1(zindx)=[];
tstp=tstrt+atom_widths';
allftr=[tstrt', tstp', latency',atom_freq1];
atoms=atomall;
atoms(zindx,:)=[];

zindx1 = find(mean(atoms,2)==0)';
atoms(zindx1,:)=[];

% Added this to see what frequencies were extracted in MP
disp('Extracted frequencies:');
disp(unique(atom_freq1));
%}
%% Use the features for clustering
X=allftr(:,1:3);
[IDX,C,SUMD,K]=kmeans_opt(X);
idx=IDX; C1= idx;
gr= 70*ones(size(idx'));
figure; scatter3(X(:,1),X(:,2),X(:,3),gr,C1,'filled');
xlabel('tstart');
ylabel('tstop');
zlabel('Latency');

cluster_atom = cell([1,K]);
cluster_MVL = cell([1,K]);
cluster_tf_MVL = cell([1,K]);

for jk=1:K
    cluster_atom = atoms((idx==jk),:);
    cluster_freq = atom_freq1((idx==jk));

     % Add these lines here to inspect what's going on
    disp(['Cluster ', num2str(jk), ' frequencies:']);
    disp(unique(cluster_freq));

     hstep=2;
     high=[15 50]; % set the required amplitude frequency range
     low=[2 12]; % set the required phase frequency range
     % Changed this, can move this around to see if different coupling is
     % present wiht lower "high" freq
     
     highfreq=high(1):2:high(2);
     amp_length=length(highfreq);
     lowfreq=low(1):1:low(2);
     phase_length=length(lowfreq);
     MP_MVL_all = zeros(amp_length,phase_length);
     tf_mvl=zeros(amp_length,phase_length,size(atoms,2));
      for i1=1:phase_length
           for j=1:amp_length
             l_freq= lowfreq(i1);
             h_freq= highfreq(j);
             [MP_MVL_all(j,i1), tf_mvl(j,i1,:)]=mpMVL(cluster_atom, h_freq, l_freq,cluster_freq);
           end
      end
      
      cluster_MVL{1,jk} = MP_MVL_all;
      cluster_tf_MVL{1,jk} = tf_mvl;
end

MP_PAC = cat(3, cluster_MVL{:});
MP_PAC_all = sum(MP_PAC,3);
plot_comodulogram(MP_PAC_all,high,low,hstep);
 
T = 1/Fs;             % Sampling period
L = length(EEG.data);             % Length of signal
timeaxis = (0:L-1)*T;
All_MVL1= cat(4, cluster_tf_MVL{:});
All_MVL = squeeze(sum(All_MVL1,4));
lf_PAC= squeeze(mean(All_MVL, 1));
figure; pcolor( timeaxis,low(1):1:low(2),lf_PAC); shading(gca,'interp');
colormap(jet);
 
hf_PAC= squeeze(mean(All_MVL, 2));
figure; pcolor( timeaxis,high(1):hstep:high(2),hf_PAC); shading(gca,'interp');
colormap(jet);

function [MVL, MVL_all] = mpMVL(atoms, h_freq, l_freq, atom_freq1)
        % The three lines after this were edited by me in an attempt to fix
        % an error
        tolerance = 0.5;
        amp_sig_in=find(abs(atom_freq1 - h_freq) < tolerance, 1);
        phs_sig_in=find(abs(atom_freq1 - l_freq) < tolerance, 1);
      if isempty(amp_sig_in) || isempty(phs_sig_in)
          MVL=0;
          MVL_all= zeros(1,size(atoms,2));
      else
          Amp_sig= atoms(amp_sig_in,:);
          Ph_sig= atoms( phs_sig_in,:);
           [m_MVL,MVL_a] = calc_MVL(Ph_sig,Amp_sig);
           [~,~,th]= MVL_surrogate_new(Amp_sig,Ph_sig);
           m_MVL(m_MVL<th)=0;
           MVL= (m_MVL); 
           MVL_all= (MVL_a); 
      end
end

function [mnsurr, stdsurr,threshold]= MVL_surrogate_new(Amp,Phase)
phase=Phase;
amplitude=Amp;
numpoints=length(Amp);   %% number of sample points in raw signal 
numsurrogate=50;   %% number of surrogate values to compare to actual value 
minskip=floor(numpoints/3);   %% time lag must be at least this big 
maxskip=numpoints-minskip; %% time lag must be smaller than this 
c1=ceil(minskip + (maxskip-minskip) .* rand(numsurrogate, 1));
skip=c1;
surrogate_m=zeros(numsurrogate,1); 

for s=1:numsurrogate
    surrogate_amplitude=[amplitude(skip(s):end) amplitude(1:skip(s)-1)];
    surrogate_m(s)=calc_MVL(phase,surrogate_amplitude);     
end 
   M_sur = surrogate_m;
   mnsurr= mean(surrogate_m);
   
   stdsurr = std( M_sur);
   pd = fitdist(M_sur,'Normal');
   ci = paramci(pd,'Alpha',.05);
   threshold=ci(2,1);
   
   function [m_MVL,MVL] = calc_MVL(Phase,Amp)
           Amp=abs(detrend(hilbert(Amp))); %figure; plot(Amp);
           Phase=angle(detrend(hilbert((Phase)))); 
           z1=(exp(1i*Phase));
           z=Amp.*(z1);% Get complex valued signal
           MVL = abs(z);
           m_MVL = abs(mean((z)));
   end
end

 function plot_comodulogram(PAC_mat,high,low,hstep)
% This function plot the comodulogram
%   INPUTS:
%   PAC_mat   : PAC values for required frequency
%   high      : High frequency range for comodulogram
%   low       : Low frequency range for comodulogram

%% plot comodulogram
figure; 
xticks = [low(1):1:low(2)]; yticks=[high(1):5:high(2)];
pcolor([low(1):1:low(2)],[high(1):hstep:high(2)],PAC_mat); shading(gca,'interp'); 
colormap(jet);
set(gca,'FontSize',10);
xlabel('Phase Frequency (Hz)','FontSize',10);ylabel('Amplitude Frequency (Hz)','FontSize',10);
set(gca,'FontName','Arial');
set(gca,'XTick',xticks);
set(gca,'YTick',yticks);
end
%}

function [m_MVL,MVL] = calc_MVL(Phase,Amp)
           Amp=abs(detrend(hilbert(Amp))); %figure; plot(Amp);
           Phase=angle(detrend(hilbert((Phase)))); 
           z1=(exp(1i*Phase));
           z=Amp.*(z1);% Get complex valued signal
           MVL = abs(z);
           m_MVL = abs(mean(nnz(z)))*0.001;
end
