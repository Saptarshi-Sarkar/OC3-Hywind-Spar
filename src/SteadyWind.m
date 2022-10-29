function [ t, velocity ] = SteadyWind(y, z, HHWindSpeed, ShearCoeff )

t = 0:100:400;
velocity = zeros(length(t), 3, length(y), length(z));
for i = 1:length(t)
    for j = 1:length(y)
        velocity(i,1,j,:) = HHWindSpeed*(z/90).^ShearCoeff;
    end
end



