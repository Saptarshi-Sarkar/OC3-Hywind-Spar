function compute_airy_wave_flow(wave, time, xy, zs)
%COMPUTE_AIRY_WAVE_FLOW computes wave velocities and accelerations using 
% linear Airy wave theory
% # This function evaluates the flow velocity and acceleartions along a 
%   a vertical line.
% # xy - horizontal coordinates along X and Y
% # zs - mutiple depth along the vertical directions
%-------------------------------------------------------------------------%
% Updated by Lin Chen on 16 Dec 2017
% Copyright (c) 2016-2017 Lin Chen <l.chen.tj@gmail.com>
%-------------------------------------------------------------------------%
global GRAVACC FLUIDDENSITY;

% Allocation
n_z = length(zs);
wave.flow_velocity = zeros(3, n_z);
wave.flow_acceleration = zeros(3, n_z);
wave.dynamic_pressure = zeros(1,n_z);
current_velocity = polyval(wave.current.polynomials,zs);

% Computation
horizon_coord =  xy(1) * cos(wave.heading_direction) ...
                +xy(2) * sin(wave.heading_direction);

% Note: surface current velocity: wave.current.velocity(1)
z_independent = (wave.angular_frequency - wave.surface_current_velocity * ...
    wave.wavenumber).* wave.amplitude;

instant_angle = wave.angular_frequency * time ...
                - wave.wavenumber * horizon_coord + wave.phase;
cos_instant_angle = cos(instant_angle);
sin_instant_angle = sin(instant_angle);
wave.surface_elevation = wave.amplitude' * sin_instant_angle;

n_frequency = length(wave.angular_frequency);
cosh_z_dependent = zeros(n_frequency, n_z);
sinh_z_dependent = zeros(n_frequency, n_z);
cosh_z_dependent_in_pressure = zeros(n_frequency, n_z);

for i = 1:1:n_z
    z_dependent = wave.wavenumber * (zs(i) + wave.water_depth);
    temp_var = cosh(z_dependent)./sinh(wave.wavenumber * ...
        wave.water_depth);
    nan_index = isnan(temp_var)==1;
    temp_var(nan_index) = 0.0;
    cosh_z_dependent(:,i) = temp_var;
    
    temp_var = cosh(z_dependent)./cosh(wave.wavenumber * ...
        wave.water_depth);
    nan_index = isnan(temp_var)==1;
    temp_var(nan_index) = 0.0;
    cosh_z_dependent_in_pressure(:,i) = temp_var;
    
    temp_var = sinh(z_dependent)./sinh(wave.wavenumber * ...
        wave.water_depth);
    nan_index = isnan(temp_var)==1;
    temp_var(nan_index) = 0.0;
    sinh_z_dependent(:,i) = temp_var;
end

% The fluid density rho and gravitational values are not included here 
% which will be 
wave.dynamic_pressure = FLUIDDENSITY * GRAVACC * ...
    (wave.amplitude .* sin_instant_angle)' * cosh_z_dependent_in_pressure;

horizon_velocity = (z_independent .* sin_instant_angle)' * ...
    cosh_z_dependent + current_velocity'; 

% Velocity in X,Y,Z directions.
wave.flow_velocity(1,:) = horizon_velocity * cos(wave.heading_direction);
wave.flow_velocity(2,:) = horizon_velocity * sin(wave.heading_direction);
wave.flow_velocity(3,:) = ((z_independent .* cos_instant_angle)' ...
                          * sinh_z_dependent);

horizon_acc = (wave.angular_frequency .* z_independent .* ...
    cos_instant_angle)'* cosh_z_dependent;

wave.flow_acceleration(1,:) =  horizon_acc * cos(wave.heading_direction);
wave.flow_acceleration(2,:) =  horizon_acc * sin(wave.heading_direction);
wave.flow_acceleration(3,:) = -(wave.angular_frequency .* z_independent ...
                                .* sin_instant_angle)'* sinh_z_dependent;
                            
end 