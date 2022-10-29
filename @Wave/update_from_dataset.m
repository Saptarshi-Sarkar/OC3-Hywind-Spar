function update_from_dataset( wave, time )
%UPDATE_FROM_DATASET reads instant velocities and accelerations from the 
% pre-calculated dataset
%-------------------------------------------------------------------------%
% Modified by Lin Chen on 17 Jan 2018
% Copyright (c) 2016-2018 Lin Chen <l.chen.tj@gmail.com>
%-------------------------------------------------------------------------%
if isempty(wave.dataset)
    error('Wave dataset unavailable.')
end
if time > wave.dataset.time(end)
    error('Simulation time beyond the time scope of wave dataset.')
end
[~,index]              = min(abs(wave.dataset.time-time));
wave.dynamic_pressure  = wave.dataset.dynamic_pressure(index,:);
wave.flow_velocity(1,:)     = wave.dataset.velocity.x(index,:);
wave.flow_velocity(2,:)     = wave.dataset.velocity.y(index,:);
wave.flow_velocity(3,:)     = wave.dataset.velocity.z(index,:);
wave.flow_acceleration(1,:) = wave.dataset.acceleration.x(index,:);
wave.flow_acceleration(2,:) = wave.dataset.acceleration.y(index,:);
wave.flow_acceleration(3,:) = wave.dataset.acceleration.y(index,:);
end
