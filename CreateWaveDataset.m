% Script used to create wave dataset
clc; clear;
global GRAVACC FLUIDDENSITY;      
GRAVACC = 9.80655; FLUIDDENSITY = 1025;
% First step: Define wave parameters and create a wave class

% Enter all wave input parameters
wavepara.heading_direction   = pi/2;         % angle in radians
wavepara.spectrum_name       = 'PiersonMoskowitz'; % 'PiersonMoskowitz'; % JONSWAP % [] regular wave;
wavepara.significant_height  = 3;
wavepara.peak_period         = 11;
wavepara.current             = struct('velocity', [0.609 0.609 0.609],'coordinate', [0 -80 -320]); %struct('velocity', [0 0], 'coordinate', [0 -320]);     % Linear current struct('velocity', [0 0], 'coordinate', [0 -320]);
wavepara.is_interacting      = 1;           % 1 for TRUE, 0 for FALSE
wavepara.wave_file           = [];

wave = Wave(wavepara);

% figure(1); hold on;
% plot(wave.angular_frequency, wave.spectrum,'linewidth',2)
% xlabel('Wave frequency (rad/s)')
% ylabel('Wave spectrum (m^2s/rad)')
% xlim([0 2.5]); box on;

% Second step: Create the dataset and save it 

TotDraft             = 120;
Platform.nStrips     = 240;
Platform.StripDepths = linspace(0,TotDraft,Platform.nStrips);

t = 0:0.0125:20;

% Create the dataset at heights = strip depths of the platform
wave.create_dataset(t, [0 0], -Platform.StripDepths' )

% the output file has the format WaveHsXX_TpXX_DirXX_UniCurXXXX_TfXXX.mat 