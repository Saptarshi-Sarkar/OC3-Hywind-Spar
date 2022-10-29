function wavenum = solve_wavenumber(wave, angular_frequency, tolerance)
%SOLVE_WAVENUMBER solves the wave number from the implicit linear dispersion
%  relation using secant method 
%-------------------------------------------------------------------------%
%
% Created by Lin Chen on 21/10/2016, TCD, Dublin 2, Ireland
% Email: L.Chen.tj@gmail.com
% Copyright (c) 2016 Lin Chen.
%-------------------------------------------------------------------------%
global GRAVACC;
% Linear dispersion relation function
fun_dispersion = @(x) (GRAVACC*x - (angular_frequency - x*wave.surface_current_velocity)*wave.surface_current_slope)*tanh(x*wave.water_depth)...
                      - (angular_frequency-x*wave.surface_current_velocity)^2;
% Note that wave.surface_current_velocity is the current velocity at Z=0

% Secand method for solving: initialization 
fun0 = fun_dispersion(0);
x0 = 0;
x1 = angular_frequency; 
fun1 = fun_dispersion(x1);
while fun1 * fun0 > 0
    x1 = x1 + angular_frequency;
    fun1 = fun_dispersion(x1);
end

% Secant method for solving the nonlinear equation, and positive solution
% is pursued only
while abs(fun1) > tolerance
    x2 = x1 - fun1 *(x1 - x0)/(fun1 - fun0);
    fun2 = fun_dispersion(x2);
    if fun2 * fun1 < 0
        x0 = x1;
        x1 = x2;
        fun0 = fun1;
        fun1 = fun2;
    else 
        x1 = x2;
        fun1 = fun2;
    end
end
wavenum = x2; % Return the solution
end