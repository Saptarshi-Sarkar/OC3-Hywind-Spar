function [ O1, O2, O3, Q1, Q2, Q3, dO1, dO2, dO3, dQ1, dQ2, dQ3, s11_B1, s22_B1, s33_B1, s12_B1, s13_B1, s23_B1 ] = BladeModeShapes( FlapMode1, FlapMode2, EdgeMode1, Twist, BldSec )
BldSec(end) = 61.5;
ddFlapMode1 = polyval(polyder(polyder(FlapMode1)),BldSec);
ddFlapMode2 = polyval(polyder(polyder(FlapMode2)),BldSec);
ddEdgeMode1 = polyval(polyder(polyder(EdgeMode1)),BldSec);

ddO1 = ddFlapMode1.*cos(Twist);
ddO2 = ddFlapMode2.*cos(Twist);
ddO3 = ddEdgeMode1.*sin(Twist);

ddQ1 = -ddFlapMode1.*sin(Twist);
ddQ2 = -ddFlapMode2.*sin(Twist);
ddQ3 = ddEdgeMode1.*cos(Twist);

%% Find O1
prevsum = 0;
dO1 = zeros(length(ddO1),1);
for i = 2:length(ddO1)
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = interp1([BldSec(i-1), BldSec(i)], [ddO1(i-1), ddO1(i)], x, 'pchip');    
    dO1(i) = prevsum + trapz(x,y);    
    prevsum = dO1(i);    
end

prevsum = 0;
O1 = zeros(length(ddO1),1);
for i = 2:length(dO1)    
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = interp1([BldSec(i-1), BldSec(i)], [dO1(i-1), dO1(i)], x, 'pchip');    
    O1(i) = prevsum + trapz(x,y);    
    prevsum = O1(i);    
end

%% Find O2
prevsum = 0;
dO2 = zeros(length(ddO2),1);
for i = 2:length(ddO2)    
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(ddO2(i-1), ddO2(i), 5);    
    dO2(i) = prevsum + trapz(x,y);    
    prevsum = dO2(i);   
end

prevsum = 0;
O2 = zeros(length(ddO2),1);
for i = 2:length(dO2)   
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(dO2(i-1), dO2(i), 5);    
    O2(i) = prevsum + trapz(x,y);    
    prevsum = O2(i);    
end

%% Find O3
prevsum = 0;
dO3 = zeros(length(ddO3),1);
for i = 2:length(ddO3)    
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(ddO3(i-1), ddO3(i), 5);    
    dO3(i) = prevsum + trapz(x,y);   
    prevsum = dO3(i);    
end

prevsum = 0;
O3 = zeros(length(ddO3),1);
for i = 2:length(dO3)   
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(dO3(i-1), dO3(i), 5);    
    O3(i) = prevsum + trapz(x,y);    
    prevsum = O3(i);    
end

% hold on
% plot(BldSec,O1,BldSec,O2,BldSec,O3)

%% Find Q1
prevsum = 0;
dQ1 = zeros(length(ddQ1),1);
for i = 2:length(ddQ1)
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(ddQ1(i-1), ddQ1(i), 5);    
    dQ1(i) = prevsum + trapz(x,y);    
    prevsum = dQ1(i);    
end

prevsum = 0;
Q1 = zeros(length(ddQ1),1);
for i = 2:length(dQ1)    
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(dQ1(i-1), dQ1(i), 5);    
    Q1(i) = prevsum + trapz(x,y);    
    prevsum = Q1(i);    
end

%% Find Q2
prevsum = 0;
dQ2 = zeros(length(ddQ2),1);
for i = 2:length(ddQ2)    
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(ddQ2(i-1), ddQ2(i), 5);    
    dQ2(i) = prevsum + trapz(x,y);    
    prevsum = dQ2(i);   
end

prevsum = 0;
Q2 = zeros(length(ddQ2),1);
for i = 2:length(dQ2)   
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(dQ2(i-1), dQ2(i), 5);    
    Q2(i) = prevsum + trapz(x,y);    
    prevsum = Q2(i);    
end

%% Find Q3
prevsum = 0;
dQ3 = zeros(length(ddQ3),1);
for i = 2:length(ddQ3)    
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(ddQ3(i-1), ddQ3(i), 5);    
    dQ3(i) = prevsum + trapz(x,y);   
    prevsum = dQ3(i);    
end

prevsum = 0;
Q3 = zeros(length(ddQ3),1);
for i = 2:length(dQ3)   
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(dQ3(i-1), dQ3(i), 5);    
    Q3(i) = prevsum + trapz(x,y);    
    prevsum = Q3(i);    
end

% hold on
% plot(BldSec,Q1,BldSec,Q2,BldSec,Q3)

%% Find s11_B1
prevsum = 0;
s11_B1 = zeros(length(ddQ3),1);
q = dO1.*dO1 + dQ1.*dQ1;
for i = 2:length(s11_B1)   
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(q(i-1), q(i), 5);    
    s11_B1(i) = prevsum + trapz(x,y);    
    prevsum = s11_B1(i);    
end

%% Find s22_B2
prevsum = 0;
s22_B1 = zeros(length(ddQ3),1);
q = dO2.*dO2 + dQ2.*dQ2;
for i = 2:length(s22_B1)   
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(q(i-1), q(i), 5);    
    s22_B1(i) = prevsum + trapz(x,y);    
    prevsum = s22_B1(i);    
end

%% Find s33_B2
prevsum = 0;
s33_B1 = zeros(length(ddQ3),1);
q = dO3.*dO3 + dQ3.*dQ3;
for i = 2:length(s33_B1)   
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(q(i-1), q(i), 5);    
    s33_B1(i) = prevsum + trapz(x,y);    
    prevsum = s33_B1(i);    
end

%% Find s12_B1
prevsum = 0;
s12_B1 = zeros(length(ddQ3),1);
q = dO1.*dO2 + dQ1.*dQ2;
for i = 2:length(s12_B1)   
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(q(i-1), q(i), 5);    
    s12_B1(i) = prevsum + trapz(x,y);    
    prevsum = s12_B1(i);    
end

%% Find s13_B1
prevsum = 0;
s13_B1 = zeros(length(ddQ3),1);
q = dO1.*dO3 + dQ1.*dQ3;
for i = 2:length(s13_B1)   
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(q(i-1), q(i), 5);    
    s13_B1(i) = prevsum + trapz(x,y);    
    prevsum = s13_B1(i);    
end

%% Find s23_B1
prevsum = 0;
s23_B1 = zeros(length(ddQ3),1);
q = dO2.*dO3 + dQ2.*dQ3;
for i = 2:length(s23_B1)   
    x = linspace(BldSec(i-1), BldSec(i),5);
    y = linspace(q(i-1), q(i), 5);    
    s23_B1(i) = prevsum + trapz(x,y);    
    prevsum = s23_B1(i);    
end
% hold on
% plot(BldSec,s12_B1,BldSec,s13_B1,BldSec,s23_B1)
