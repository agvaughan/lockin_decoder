function [smoothed_output, demodulation_freq, demodulation_phase] = decode_lockin_fn( ...
    lockin_sample_input, ...
    chopper_reference_input, ...
    gold_standard_input , ...
    chopper_frequency, ...
    sample_rate,...
    do_plot)

% Serves as a software lockin amplifier, by decoding the lockin_sample_input

dbstop if error

%% Parameters
enforce_chopper_frequency  = 1; % FORCE use of exactly chopper frequency?
enforce_phase_frequency    = 1; % FORCE use of zero phase.
%freq_detection_method      = 'pwelch'; %'pwelch','interpolation','smoothing' % Attempt to find peak of FFT more precisely.

use_quadrature = 1;

%LowPassCutoff = round(chopper_frequency/10);% t_1/2 scalar in ms
LowPassCutoff =  5;% t_1/2 scalar in ms

%% Parse inputs

n_seconds = length(lockin_sample_input)/sample_rate;

if nargin < 6,
    do_plot = 0;
end

if nargin < 5,
    chopper_frequency = 470;
end

if nargin < 4,
    sample_rate = 10000;
end

if nargin < 3,
    gold_standard_input = NaN;
end

% If chopper_reference_input is not included, we can make do with the sample
% vector

if nargin < 2 || isempty(chopper_reference_input) || isnan(chopper_reference_input)
    chopper_reference_input = NaN;
end


if do_plot,
    close all
end


%% Pre-process inputs.
% Ensure column vector
lockin_sample_input     = lockin_sample_input(:);
chopper_reference_input = chopper_reference_input(:);
gold_standard_input     = gold_standard_input(:);
n_samples               = length(lockin_sample_input);

% Zscore incoming signals
lockin_sample_input = zscore(lockin_sample_input);
chopper_reference_input = zscore(chopper_reference_input);


%% Take FFT of 

% FFT of Output Light
NFFT = 2^(nextpow2(n_samples)); % Next power of 2 from length of y

% Frequencies for FFT
fft_freqs = sample_rate/2*linspace(0,1,NFFT/2);

% FFT of sample signal and reference signal
% We always work on the chopper_reference_input, but if it is not provided then we simply copy the lockin_sample_input.
sample_fft = fft(lockin_sample_input,NFFT);
if  isempty(chopper_reference_input) || all(isnan(chopper_reference_input))
    % Defaults to using sample (via reference)
    reference_fft   = sample_fft;     % Work on this from now on.
    chopper_reference_input = lockin_sample_input; % Work on this from now on.
else
    reference_fft = fft(chopper_reference_input,NFFT);
end


%% Find center frequency of chopper reference input.

if enforce_chopper_frequency
    % Just use nominal frequency.
    demodulation_freq = chopper_frequency;
else
    % Phases of the reference and sample ffts.
    reference_phases = unwrap(angle(reference_fft(1:end/2)));
    
    % Calculate power of reference, excitation, and emission power
    reference_power = 2*abs(reference_fft(1:NFFT/2+1));
    
    % Normalize FFT of excitation and emission
    freq_ind_near_chopper  = find(fft_freqs > (chopper_frequency*0.9),1) : find(fft_freqs > (chopper_frequency*1.1),1);
    power_near_chopper  = reference_power(freq_ind_near_chopper);
    reference_power     = reference_power / max(power_near_chopper);
    [demodulation_freq,freqs_pwelch,power_pwelch,f_ind_pwelch] = find_peak_freq_pwelch( chopper_reference_input, chopper_frequency, sample_rate );
    
    if do_plot
        [~,freq_ind_near_chopper] = max(reference_power(freq_ind_near_chopper));
        %demodulation_phase_raw = phase_near_chopper(freq_ind_near_chopper);
        %demodulation_freq_raw  = freqs_near_chopper(freq_ind_near_chopper);
        
        figure(1); clf; hold on; clear h
        h(1) = plot(freqs_pwelch,power_pwelch,'g'); axis tight;
        h(2) = plot([1 1]*chopper_frequency,ylim,'k');
        h(3) = plot([1 1]*freqs_pwelch(f_ind_pwelch),ylim,'r');
        h(4) = plot([1 1]*demodulation_freq,ylim,'b:');
        set(h,'LineWidth',1.5)
        legend({'Power (Welch)','Chopper actual','Freq Pwelch','Freq FFT'})
        title('Pwelch power and identified frequencies')
        xlabel('Frequency (hz)')
        ylabel('PSD')
    end
    
end


%% Find phase via demodulation

if enforce_phase_frequency
    % Just set phase to zero.
    demodulation_phase = 0;
    phase_hedge           = demodulation_freq / sample_rate *pi;
    ref_cycles            = (n_samples / sample_rate) * demodulation_freq;
    reference_sinusoid    = cos(      - phase_hedge + demodulation_phase  + linspace( 0, 2*pi*ref_cycles, n_samples))';
    reference_sinusoid_90 = cos( pi/2 - phase_hedge + demodulation_phase  + linspace( 0, 2*pi*ref_cycles, n_samples))';
    
else
    % Find phase via demodulation.
    phase_hedge = demodulation_freq / sample_rate * pi;
    % Demodulate the reference input to get the phase
    demod_phase = demod(chopper_reference_input,demodulation_freq*1.00001,sample_rate,'pm');
    % F
    demod_phase_medfilt = unwrap(medfilt1( demod_phase, sample_rate/2)  + phase_hedge);
    reference_sinusoid = modulate( demod_phase_medfilt,demodulation_freq,sample_rate,'pm',1);
    reference_sinusoid_90 = modulate( demod_phase_medfilt+pi/2,demodulation_freq,sample_rate,'pm',1);
    
    if do_plot
        figure(9); clf
        subplot(2,1,1); hold on
        plot(demod_phase)
        plot(demod_phase_medfilt)
        if use_synthetic_data
            hline(synthetic_phase,'r:')
            % disp(synthetic_phase - mean(demod_phase))
            % disp(synthetic_phase - mean(demod_phase_medfilt))
        end
        legend('Demodulation phase','medfilt(demodulation phase)')
        %
        subplot(2,1,2); cla; hold on
        strips(zscore([chopper_reference_input' reference_sinusoid']),sample_rate*2,sample_rate)
        xlim([0 sample_rate/50])
        legend({'Input','Demodulated'})
    end
    demodulation_phase = NaN;
    
end


%% Multiply the output_lockin by the reference signal

output_lockin_0deg = lockin_sample_input .* reference_sinusoid;

if use_quadrature,
    % Use the vector strength of quadrature decoding.
    output_lockin_90deg  = lockin_sample_input .* reference_sinusoid_90;
    output_lockin_all = sqrt((output_lockin_0deg).^2 + (output_lockin_90deg).^2);
    %output_lockin_all = abs( output_lockin_0deg + i * output_lockin_90deg );
    output_lockin_to_use = output_lockin_all;
else
    % Just use the best estimate of phase.
    output_lockin_to_use = output_lockin_0deg;
end

%% Low-pass filter lockin output

% Low-pass time-constant
% For conv approach, this is the t1/2 of the exponential kernel
% For FiltfiltM with a butterworth filter, this is the shoulder frequency.
% In almost all cases, it is applied twice.

% Can use FilterM / FilterX for speed.
% Uses the FilterM / FilterX functions from
% Jan Simon, Heidelberg, (C) 2011 matlab.THISYEAR(a)nMINUSsimon.de
% Single dimension, so can just use FilterX as the direct .mex file
% call.
% This is about 20% faster than filtfilt.

% Reload filters ; if not, generate.
try
    d = load('butterworth_filter_saved.mat');
    if d.StopBandFrequency ~= LowPassCutoff
        error('*** Filter has changed - rebuilding. ***')
    end
    
catch err
    fprintf('Error loading filter :: %s\n',err.message')    
    % Lowpass workes fine.
    %%
    %LowPassCutoff = 20
    d = designfilt('lowpassiir','FilterOrder',10,... % Can't be 20
        'PassbandFrequency', LowPassCutoff ,...
        'PassbandRipple',0.1,...
        'SampleRate',sample_rate);    
    %freqz(d,logspace(-1,2,200),sample_rate)
    %keyboard
    %%
    save butterworth_filter_saved;

end


%% Compiled equivalent of FilterX

% Subtract mean
offset = mean(output_lockin_to_use);
% Add a ~1s buffer on either size to minimize start/stop artifacts.
buffer_size = 1 * sample_rate;
pre_buffer = ones(buffer_size,1) * (output_lockin_to_use(1)-offset);
post_buffer = ones(buffer_size,1) * (output_lockin_to_use(end)-offset);
% Filter.
output_lockin_to_use_smoothed = filtfilt(d, [ pre_buffer; output_lockin_to_use-offset; post_buffer]);
output_lockin_to_use_smoothed = output_lockin_to_use_smoothed([buffer_size+1:end-buffer_size]) + offset;
assert( ~any(isnan(output_lockin_to_use_smoothed)), 'FOUND NANS IN output_lockin_to_use_smoothed')

if do_plot
    figure(10); clf; hold on
    h = plot(zscore([output_lockin_to_use output_lockin_to_use_smoothed gold_standard_input]))
    set(h,'LineWidth',2)
    legend({'Unsmoothed','Smoothed','Gold Standard'})
    axis tight;
    title(sprintf('Buffer size :: %.1fs',buffer_size/sample_rate))
    vline(sample_rate),vline(n_samples+sample_rate)
end

smoothed_output  = output_lockin_to_use_smoothed;

% DONE 



%% Plotting for sanity-check

if do_plot
    if use_synthetic_data
        fprintf('Reference :: %.3fhz (vs. %.3f) @ %.2f radians ( vs %.2f  - error %.2f/%.2fdeg )\n',...
            demodulation_freq, chopper_frequency,...
            demodulation_phase,synthetic_phase, ...
            demodulation_phase-synthetic_phase,(demodulation_phase-synthetic_phase)/pi*180)
        %fprintf('Sample    :: %.2fhz (vs. %.2f) @ %.2f radians (%.2f deg) [ expect :: %.2f ]\n',sample_freq,chopper_frequency,sample_phase,sample_phase*180/pi,synthetic_phase)
        
    else
        fprintf('Reference :: %.3fhz (vs. %.3f) @ %.2f radians (%.2f deg)\n',demodulation_freq,chopper_frequency,demodulation_phase,demodulation_phase*180/pi)
        %fprintf('Sample    :: %.2fhz (vs. %.2f) @ %.2f radians (%.2f deg)\n',sample_freq,chopper_frequency,sample_phase,sample_phase*180/pi)
    end
    
    figure(4); clf;
    t_inds = (1:100) + 5000;
    % Reference power.
    subplot(3,1,1); cla; hold on
    plot([fft_freqs 0],reference_power)
    %[pxx,fs] = pwelch(lockin_sample_input,2048,0,2048,sample_rate);
    set(gca,'YScale','log')
    xlabel('Frequency hz'); ylabel('Power / hz');
    xlim(chopper_frequency * [0.9999 1.0001])
    h = vline(demodulation_freq,'r:'); set(h,'LineWidth',2)
    %
    subplot(3,2,3); cla; hold on
    plot(fft_freqs,reference_phases,'.','MarkerSize',20)
    h = vline(demodulation_freq,'--'); set(h,'LineWidth',2)
    h = hline(demodulation_phase,'r--'); set(h,'LineWidth',2)
    ylim([-2*pi 2*pi]); xlim(demodulation_freq * [0.5 2])
    ylabel({'Phase','Radians'})
    %
    subplot(3,2,4); cla; hold on
    %%
    t_inds_multi = ceil(bsxfun(@plus,1:(2*sample_rate/chopper_frequency),linspace(0,n_samples*0.9,11)'));
    for ind = 1:10,
        plot(      lockin_sample_input(t_inds_multi(ind,:)),reference_sinusoid(t_inds_multi(ind,:)),'LineWidth',1.5)
    end
    %%
    xlabel('Sample Signal'); ylabel({'Reference','sinusoid'})
    minmax = @(x) [min(x) max(x)];
    plot( minmax(lockin_sample_input(t_inds)), minmax(reference_sinusoid(t_inds)),'r')
    xlim(minmax(lockin_sample_input(t_inds))); ylim(minmax(reference_sinusoid(t_inds)))
    axis square
    %
    subplot(3,1,3); cla; hold on
    t = linspace(0,n_seconds,n_seconds*sample_rate);
    plot(t(t_inds_multi(10,:)),zscore(lockin_sample_input(t_inds_multi(10,:))),'k')
    plot(t(t_inds_multi(10,:)),zscore(reference_sinusoid(t_inds_multi(10,:))),'b','LineWidth',2)
    axis tight;xlim([0 5/chopper_frequency])
    if ~isnan(chopper_reference_input)
        plot(t(t_inds_multi(10,:)),zscore(chopper_reference_input(t_inds_multi(10,:))),'r')
        legend({'Sample input','Recovered sinusoid','Reference input'})
    else
        legend({'Sample signal','Recovered sinusoid'},'Location','EastOutside')
    end
    axis tight; ylabel({'Sinusoids','(zscored)'}); xlabel('Time (s)')
    
    
    figure(5);clf
    t = linspace(0,n_seconds,n_seconds*sample_rate);
    
    subplot(3,1,1); hold on;
    % Time frame to plot
    %t_inds = (1:5*ceil(sample_rate/1)) + ceil(n_samples/2);
    %t_inds = sample_rate:(n_samples-sample_rate);
    t_inds = 1:n_samples;
    
    
    % Plot pwelch of all interesting signals, except those that are NaN;
    signals_to_plot = {lockin_sample_input, reference_sinusoid, output_lockin_to_use, output_lockin_to_use_smoothed,gold_standard_input};
    signal_legend = {'lockin_sample_input', 'reference_sinusoid', 'output_lockin_to_use', 'output_lockin_to_use_smoothed','gold_standard_input'};
    notnan_pwelch_to_plot = cellfun(@(x) ~all(isnan(x)), signals_to_plot);
    %
    this_signals_to_plot = cell2mat(signals_to_plot(notnan_pwelch_to_plot));
    signal_legend = signal_legend(notnan_pwelch_to_plot);
    pwelch(this_signals_to_plot,ceil(sample_rate/2),[],ceil(sample_rate/4),sample_rate,'power'),
    legend(signal_legend,'Location','EastOutside','Interpreter','none')
    set(get(gca,'children'),'LineWidth',1.5)
    xlim([1 chopper_frequency* 2.5]/1000); ylim([-150 0])
    vline([20 33 60],'b');
    vline(LowPassCutoff,'r')
    xlabel('Hz'); ylabel('Power/hz')
    %
    
    subplot(3,1,2); cla; hold on
    this_signals_to_plot = cell2mat(signals_to_plot(notnan_pwelch_to_plot));
    this_signals_to_plot(:,4:end) = NaN;
    set(get(gca,'children'),'LineWidth',1.5)
    h = plot(t(t_inds),this_signals_to_plot(t_inds,:));
    set(h(3),'LineWidth',2)
    title('Input Signals')
    hline(0)
    xlims = xlim;
    
    subplot(3,1,3); cla; hold on
    % Plots output_lockin_to_use_smoothed,gold_standard_input
    this_signals_to_plot = cell2mat(signals_to_plot(notnan_pwelch_to_plot));
    this_signals_to_plot = this_signals_to_plot(:,3:end);
    plot(t(t_inds),zscore(this_signals_to_plot(t_inds,:)))
    set(get(gca,'Children'),'LineWidth',2)
    title('Output Signals')
    legend(signal_legend(3:end),'Interpreter','none')
    xlim(xlims)
    
    hline(0,'r')
    
    
    %%
    figure(6); clf
    subplot(2,2,1);
    strips(zscore([lockin_sample_input reference_sinusoid output_lockin_to_use]),1,sample_rate)
    xlim([0 1/100])
    legend({'Raw Input','Lockin Reference'})
    
    %%
    subplot(2,2,2); cla; hold on
    n_subregions = 10;
    samples_per_subregion = 10*sample_rate/chopper_frequency;
    t_inds_multi = ceil(bsxfun(@plus            ,...
        1:samples_per_subregion   ,...
        linspace(n_samples*0.1,n_samples*0.9,n_subregions)'));
    for i = 1:size(t_inds_multi,1)
        plot( ...
            gold_standard_input(t_inds_multi(i,:)),...
            output_lockin_to_use_smoothed(t_inds_multi(i,:)),...
            '.','LineWidth',1.5)
    end
    xlabel('Gold standard')
    ylabel('Output')
    xlims = xlim; xlim(xlims + [-0.1 0.1]*range(xlims));
    ylims = ylim; ylim(ylims + [-0.1 0.1]*range(ylims));
    %%
    subplot(2,2,3);
    strips(zscore([gold_standard_input(1:10:end),output_lockin_to_use_smoothed(1:10:end)]),1/2,sample_rate)
    set(get(gca,'Children'),'LineWidth',1.5)
    legend({'Gold Standard','Lockin Output'})
    xlim([0 1/10])
    
    subplot(2,2,4);
    strips(zscore([gold_standard_input(1:10:end),output_lockin_to_use_smoothed(1:10:end)]),1/2,sample_rate)
    set(get(gca,'Children'),'LineWidth',1.5)
    legend({'Gold Standard','Lockin Output'})
    xlim([0 1/100])
    %xlim([0 sample_rate/100])
    
    
    %%
    tilefigs
    keyboard
end

%
% %% Envelope method?
% figure(11)
% tmp = output_lockin_to_use;
% n = 5;
% subplot(3,1,1)
% envelope(tmp,n*sample_rate/chopper_frequency,'rms');
% axis tight
% subplot(3,1,2); cla; hold on
% [env_up,env_low] =  envelope(tmp,n*sample_rate/chopper_frequency,'rms');
% plot(filtfilt(d_incoming,env_up - env_low)); axis tight
% plot(filtfilt(d_incoming,gold_standard_input)); axis tight
%
% subplot(3,1,3)
% plot(output_lockin_to_use_smoothed)
%
% smoothed_output = (env_up - env_low);

end % END DECODE_LOCKIN MAIN FUNCTION

%% %%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%% %%

function [demodulation_freq,freqs_pwelch,power_pwelch,f_ind_pwelch] = find_peak_freq_pwelch( signal, chopper_frequency, sample_rate )
freqs_to_sample = linspace(chopper_frequency*0.999,chopper_frequency*1.001,1000);
[power_pwelch,freqs_pwelch] = pwelch(signal,sample_rate/2,[],freqs_to_sample,sample_rate);
[~,f_ind_pwelch] = max(power_pwelch);
demodulation_freq  = freqs_pwelch(f_ind_pwelch);
end

function [demodulation_freq,demodulation_phase] = find_peak_freq_interpolated( freqs, powers,phases ,interpolation_scale)
% Find peak power and phase across a frequency range interpolated to "interpolation_scale"-times resolution
freqs_interp = linspace(min(freqs),max(freqs),length(freqs)*interpolation_scale);
power_interp = interp1( ...
    freqs, ... % x
    powers, ... % y
    freqs_interp, ... % new x
    'spline' );
phase_interp = interp1( ...
    freqs, ... % x
    phases, ... % y
    freqs_interp, ... % new x
    'spline' );

[~,freq_ind_interp] = max(power_interp);
demodulation_freq   = freqs_interp(freq_ind_interp);
demodulation_phase  = mod(phase_interp(freq_ind_interp),2*pi);

end

