function solve_wavenumbers( obj )
%GETWAVENUMBER solves wave number for given wave frequency vector
%-------------------------------------------------------------------------%
% Modified on 17 Dec 2017
% Copyright (c) 2016-2017 Lin Chen <l.chen.tj@gmail.com>
%-------------------------------------------------------------------------%
% disp('- Now solving wave numbers which takes a moment.');
obj.wavenumber = zeros(size(obj.angular_frequency));
for i  = 1:length(obj.wavenumber)
    obj.wavenumber(i) = solve_wavenumber(obj, obj.angular_frequency(i),...
        1E-12);
end
% disp('- done.');
end