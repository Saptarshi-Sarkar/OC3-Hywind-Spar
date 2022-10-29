function create_dataset( wave, times, xy, zs )
%CREATE_DATASET creates dataset to be used in subsequent analysis
%-------------------------------------------------------------------------%
% Currently, the wave velocity and acceleration are considered at x=0 & y=0
% Updated by Lin Chen on 16 Dec 2017, TCD, Dublin 2, Ireland
% Email: L.Chen.tj@gmail.com
% Copyright (c) 2016-2017 Lin Chen.
%-------------------------------------------------------------------------%
% Solve wave profile for time and depth

n_time = length(times);
n_z = length(zs);

velocity = struct('x', zeros(n_time, n_z), ...
                  'y', zeros(n_time, n_z), ...
                  'z', zeros(n_time, n_z));
                  
acceleration = struct('x', zeros(n_time, n_z), ...
                      'y', zeros(n_time, n_z), ...
                      'z', zeros(n_time, n_z));

wave.dataset = struct('time', times, 'planar_location', xy,...
                      'vertical_location', zs, 'velocity', ...
                       velocity, 'acceleration', acceleration, ...
                       'surface_elevation', zeros(n_time,1), ...
                       'dynamic_pressure', zeros(n_time, n_z));
clear velocity acceleration;

for i = 1:n_time
    disp(times(i));
    wave.compute_airy_wave_flow(times(i), [0,0], zs);
    wave.dataset.surface_elevation(i)= wave.surface_elevation;
    wave.dataset.dynamic_pressure(i,:) = wave.dynamic_pressure;
    wave.dataset.velocity.x(i,:)     = wave.flow_velocity(1,:);
    wave.dataset.velocity.y(i,:)     = wave.flow_velocity(2,:);
    wave.dataset.velocity.z(i,:)     = wave.flow_velocity(3,:);
    wave.dataset.acceleration.x(i,:) = wave.flow_acceleration(1,:);
    wave.dataset.acceleration.y(i,:) = wave.flow_acceleration(2,:);
    wave.dataset.acceleration.z(i,:) = wave.flow_acceleration(3,:);
end

if wave.current.velocity(1) == 0
    Cur = 'NoCur';
elseif wave.current.velocity(1) == wave.current.velocity(2)
    Cur = 'UniCur';
elseif wave.current.velocity(1) ~= wave.current.velocity(2)
    Cur = 'LinCur';
end

% Save the dataset
filename = strcat('WaveHs', strrep(num2str(wave.significant_height),...
    '.', '_'), '_Tp', strrep(num2str(wave.peak_period),...
    '.', '_'), '_Dir', strrep(num2str(rad2deg(wave.heading_direction)),...
    '.', '_'), '_', Cur , strrep(num2str(wave.current.velocity(1)),...
    '.', '_'), '.mat');
save(filename, 'wave');

end