%
%   contrec1anal_ELS.m
%       For analyzing data obtained with contrec1
%

%% EOD Detection
recdata_rect = sqrt(recdata.^2);
recdata_smooth = medfilt1(recdata_rect',5)'; %% Uses a filter to smooth the data to find peaks
[a,b] = size(recdata_smooth);
eod_binned = zeros(a,b); %% Bins times in the recording where the threshold is crossed, = 0 or 1
ipi = [];
peod=[];
figure
for i=1:a,
    disp(' ')
    disp(['Running rep #' sprintf('%s', num2str(i))])
    disp(' ')
    eod = [];
    amp = [];
    a1 = subplot(211);
    plot(rectime,recdata(i,:),'k-')
    ylabel('Potential (V)')
    a2 = subplot(212);
    plot(rectime,recdata_smooth(i,:),'k-')
    hold on
    xlabel('Time (s)')
    ylabel('Rectified, Smoothed Potential (V)')			%% the smoothed EOD is only positive
    linkaxes([a1 a2],'x')
    zoom on
    if i==1,
        thold = input('set threshold... ');                     % input threshold
        snip_width = input('set snippet width in seconds (1 or .02)... ');             % input snippet width
    end
    j=1;
    [maxtab, mintab] = peakdet(recdata_smooth(i,:), thold);     % detect peaks above threshold
    if size(maxtab,1)==0
    	continue
    end
    eod = rectime(maxtab(:,1));
    amp = recdata_smooth(i,maxtab(:,1));
    plot(eod,amp,'g.')                                          % plot detected peaks
    if i==1,
        input('Hit return to continue.');
    end
    hold off
    zoom on
    eod_bin_tmp = bin(eod,0.05);
    eod_binned(i,1:length(eod_bin_tmp)) = eod_bin_tmp;
    ipi = [ipi diff(eod).*1000];
    peod_tmp=[];
	if snip_width==1,
	    for j=1:size(maxtab,1),				%% creating an array of the potential eods in the recording
    		if maxtab(j,1)<48828
    			peod_tmp(j,:)=recdata(i,1:48828);
    		elseif maxtab(j,1)>439455
    			peod_tmp(j,:)=recdata(i,439455:488282);
    		else
    			peod_tmp(j,:)=recdata(i,(maxtab(j,1)-24414):(maxtab(j,1)+24413));
    		end
    	end
    elseif snip_width==.02
    	for j=1:size(maxtab,1),				%% creating an array of the potential eods in the recording
    		if maxtab(j,1)<976
    			peod_tmp(j,:)=recdata(i,1:976);
    		elseif maxtab(j,1)>487307
    			peod_tmp(j,:)=recdata(i,487307:488282);
    		else
    			peod_tmp(j,:)=recdata(i,(maxtab(j,1)-488):(maxtab(j,1)+487));
    		end
    	end
    end
    peod=[peod;peod_tmp];
end
close(gcf)

% %% Power Spectrum
% nfft = 2^nextpow2(b);
% freq = srate/2*linspace(0,1,nfft/2);
% for i=1:a,
%     p = fft(recdata(i,:)-mean(recdata(i,:)),nfft)/length(recdata);
%     praw = abs(p(1:nfft/2));
%     power(i,:) = praw/max(praw);
% end
% powerDB = 10*log10(power);
% figure
% subplot(221)
% semilogx(freq,powerDB,'k-',freq,mean(powerDB),'r-')
% title('Power Spectrum');
% xlabel('frequency [Hz]');
% ylabel('power (dB)');

% Power Spectrum
for i=1:a,
    recanal = recdata(i,:) - mean(recdata(i,:));
    [power(i,:),freq] = pwelch(recanal,length(recanal),[],[],srate);
    power(i,:) = power(i,:)./max(power(i,:));
end
powerDB = 10*log10(power);
figure
subplot(221)
semilogx(freq,powerDB,'k-',freq,mean(powerDB),'r-')
title('Power Spectrum');
xlabel('frequency [Hz]');
ylabel('power (dB)');

%% ISI Distribution
for i=1:a,
    [ipi_t, ipi_pd(i,:), ipi_mean(i), ipi_cv(i)] = isidist(eod_binned(i,:),0.05,30,0.05);
end
subplot(222)
plot(ipi_t,ipi_pd,'k-',ipi_t,mean(ipi_pd),'r-')
title('Interspike interval distribution');
xlabel('interspike interval [msec]');
ylabel('probability per bin');
xlim([0 30])

%% Power Spectral Density of Spike Train
for i=1:a,
    [ipi_f,ipi_Pxxn(i,:),ipi_tvect,ipi_Cxx(i,:)] = psautospk(eod_binned(i,:),0.05);
end
subplot(223)
plot(ipi_f,ipi_Pxxn,'k-',ipi_f,mean(ipi_Pxxn),'r-')
title('Power spectrum of the spike train');
xlabel('Frequency [Hz]');
ylabel('Power spectral density [(spk/sec)^2/Hz]');
subplot(224)
plot(ipi_tvect,ipi_Cxx,'k-',ipi_tvect,mean(ipi_Cxx),'r-')
title('Autocorrelation function');
xlabel('time lag [msec]');
ylabel('Autocorrelation [(spk/sec)^2]')