classdef Wave < handle
%WAVE defines the wave class, having the following functionality.
%-------------------------------------------------------------------------%   
% # Get wave velocity and acceleration fields;
% # Able to consider the wave direction and current;
% # Able to consider the wave direction and create wave fields in 
%   in three dimension
% Updated by Lin Chen on 19 Jan 2018.
% Copyright (c) 2016-2018 Lin Chen <l.chen.tj@gmail.com>
%-------------------------------------------------------------------------%  

properties
    
    % Basic parameters.
    heading_direction;    % Anlge with X-axis, rad.
    significant_height;   % [m]
    peak_period;          % [s]
    idx                   % index of tower coordinates where loads are to be calculated
    spectrum_name;        % If empty, representing regular wave.
    water_depth;          % [m]      
    
    % Derived parameter
    amplitude;            % A vector saving amplitude of wave components.
    wavenumber;           
    phase;
    angular_frequency;    % [rad]
    spectrum;             % Matrix of two column and the first column is 
                          % frequency and the second is spectral amplitude.
    
    % Current: presently consider current is the same or opposite direction
    % of the waves.
    current;
    is_interacting;     % Whether to consider the wave current interaction.
    surface_current_velocity;
    surface_current_slope;
end
    
properties 
    flow_velocity;       
    flow_acceleration;   
    surface_elevation;   
    dynamic_pressure;    
end

% When wave velocities and accelearations have been solved before 
% the simulation and then used by interpolation. At each time step the
% wave velocity and acceleration will be read into flow_velocity and
% flow_acceleration. If the flow_velocity and flow_acceleration are 
% calculated in real-time, the dataset is used to save the time 
% histories.
properties             
    dataset;             %
end

methods
   
    function obj = Wave(opt)

        global GRAVACC;

        if isempty(opt.wave_file)

            % Default settings
%             obj.heading_direction = 0;
%             obj.spectrum_name = [];
% %             obj.spectrum_name      = 'JONSWAP';
%             obj.spectrum_name      = 'PiersonMoskowitz';
%             obj.significant_height = 6; %3.66; %9.14; %3.0;
%             obj.peak_period = 10; %9.7;  %13.6; %10.0;
%             
% %             One component indicates uniform current and the first value
% %             of the current must be the surface current velocity.
%             obj.current = struct('velocity', [0.609 0.457 0.305 0.152 0],...
%             'coordinate', [0 -80 -160 -240 -320]);      % Linear current 
% %             obj.current = struct('velocity', [0.609	0.609], 'coordinate', [0 -320]); % Uniform current
% %             obj.current = struct('velocity', [0 0], 'coordinate', [0 -320]);
%             
%             obj.is_interacting = 1;
            
            
            % My settings
            obj.heading_direction  = opt.heading_direction;
            obj.spectrum_name      = opt.spectrum_name;
            obj.significant_height = opt.significant_height;
            obj.peak_period        = opt.peak_period;
            obj.current            = opt.current;           
            obj.is_interacting     = opt.is_interacting;
            
            % =====X==X====
            
            if obj.is_interacting
                obj.surface_current_velocity = obj.current.velocity(1);
                obj.surface_current_slope = (obj.current.velocity(1) - obj.current.velocity(2))/(obj.current.coordinate(1) - obj.current.coordinate(2));
            else
                obj.surface_current_velocity = 0;
                obj.surface_current_slope    = 0;
            end
            
            obj.water_depth = 320;
            
            n_coordinate = length(obj.current.coordinate);
            if n_coordinate > 4 
                n_coordinate = 4;
            end
            
            obj.current.polynomials = polyfit(obj.current.coordinate,...
                obj.current.velocity,n_coordinate-1);
            obj.current.poly_order = n_coordinate-1;

            % For the case of a regular wave.
            if isempty(obj.spectrum_name)
                
                obj.amplitude = obj.significant_height/2.0;
                obj.phase = 0.0;
                obj.angular_frequency = 2*pi/obj.peak_period;
                obj.solve_wavenumbers();

            % For the case of irregular waves.
            else
                % Number of bins for discretization of the spectrum.
                n_frequency = 512;   %8192
                frequency_limit = 200;

                frequency_min = 0.0001/2/pi; % [Hz]
                frequency_max = 1/0.01/2;
                if frequency_max > frequency_limit
                    frequency_max = frequency_limit;
                end

                if frequency_max > 1/obj.peak_period * 4
                    frequency_max = 1/obj.peak_period * 4;
                end

                obj.angular_frequency  = 2*pi*linspace(frequency_min, ...
                    frequency_max, n_frequency)';

                % Computation spectrum.
                temporal_constant = obj.angular_frequency * ...
                    obj.peak_period/2/pi;

                switch obj.spectrum_name
                    case 'JONSWAP'                    

                        if obj.peak_period/sqrt(obj.significant_height) <= 3.6
                            gama = 5;
                        elseif obj.peak_period/sqrt(obj.significant_height) > 5
                            gama = 1;
                        else
                            gama = exp(5.75-1.15*obj.peak_period/...
                                sqrt(obj.significant_height));
                        end

                        sigma = ones(size(obj.angular_frequency)) * 0.09;
                        index = obj.angular_frequency<=2*pi/obj.peak_period;
                        sigma(index) = 0.07;

                        obj.spectrum = 1/(2*pi)*5/16 * obj.significant_height^2 ...
                            * (1-0.287*log(gama)) * obj.peak_period *...
                            (temporal_constant.^(-5)) .* ...
                            (exp(-1.25 * temporal_constant.^(-4))) .* ...
                            (gama.^(exp(-0.5*((temporal_constant - 1) ./ sigma) .^2)));

                    case 'PiersonMoskowitz'

                        obj.spectrum = 1/(2*pi)*5/16 * obj.significant_height^2 * ...
                            obj.peak_period * (temporal_constant.^(-5)) .* ...
                            (exp(-1.25 * temporal_constant.^(-4))); 

                    otherwise 

                        error('Wave spectrum unknown.');

                end

                obj.solve_wavenumbers();
                
                % Effect of the current.
                if abs(obj.surface_current_velocity)>10E-6

                    obj.spectrum = obj.spectrum * 4 ./...
                        ((1+sqrt(1+4*obj.angular_frequency*obj.surface_current_velocity/GRAVACC)).^2 .* ...
                        sqrt(1+4*obj.angular_frequency*obj.surface_current_velocity/GRAVACC) ); 
                    
                    % Equilibrium range for the case of waves on adverse 
                    % current.
                    spectrum_equil_range = 0.015 * GRAVACC^2 ./ ...
                        (obj.angular_frequency - obj.wavenumber * ...
                        obj.surface_current_velocity).^5 ./ (1 + 2*obj.surface_current_velocity * ...
                        (obj.angular_frequency - obj.wavenumber * ...
                        obj.surface_current_velocity) / GRAVACC);

                    if obj.current.velocity < 0
                        for i = 1:n_frequency
                            if obj.spectrum(i) > spectrum_equil_range(i)
                                obj.spectrum(i) = spectrum_equil_range(i);
                            end
                        end
                    end

                end

                angular_frequency_interval = obj.angular_frequency(2) - obj.angular_frequency(1);
                obj.amplitude= sqrt(obj.spectrum) * sqrt(2*angular_frequency_interval);

                % Get random Phase or load predefined random set.
%                 obj.phase             = 2*pi*rand(obj.nFrequency,1);
                loaded = load('phase.mat'); % For comparison.
                obj.phase = loaded.phase;
            end

        elseif strcmp(opt.wave_profile, 'load')

%             if strcmp(opt.wave_profile, 'load')
                if isempty(opt.wave_file)
                    error('Wave data file unknown.');
                end
                loaded = load(opt.wave_file);
                obj = loaded.wave;
                clear loaded;
%             end

        end
    end % Constructor ends.
end % Method definition ends.
end % Wave class definition ends.