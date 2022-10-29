function [yout, Cons] = ode4(F,t0,h,tfinal,y0,Cons0)
% ODE4  Classical Runge-Kutta ODE solver.
%   yout = ODE4(F,t0,h,tfinal,y0) uses the classical
%   Runge-Kutta method with fixed step size h on the interval
%      t0 <= t <= tfinal
%   to solve
%      dy/dt = F(t,y)
%   with y(t0) = y0.

%   Copyright 2014 - 2015 The MathWorks, Inc.
t = t0 : h : tfinal-h; LT = 0; Ltrem = 0;
n = length(t);
yout = zeros(n+1,length(y0));   yout(1,:) = y0;
Cons  = zeros(n+1,length(Cons0)); Cons(1,:) = Cons0;
y = y0;
cons = Cons0;

tic;
for i = 1:n
  [s1,Cons1] = F(t(i),y, cons);
  [s2,Cons2] = F(t(i)+h/2, y+h*s1/2, Cons1);
  [s3,Cons3] = F(t(i)+h/2, y+h*s2/2, Cons2);
  [s4,Cons4] = F(t(i)+h, y+h*s3, Cons3);
  y = y + h*(s1 + 2*s2 + 2*s3 + s4)/6;
  cons  = (Cons1 + 2*Cons2 + 2*Cons3 + Cons4)/6;
  yout(i+1,:) = y;  
  Cons(i+1,:)  = cons;  
  
  %% Calculating remaining time
  if t(i) - LT > 5.0
      trem = (n-i)*toc/i/60;
      if LT == 0
          fprintf('Simulation time elapsed: %.0f secs. Computation time left: %.0f minutes.',t(i),trem);
      else
          fprintf([repmat('\b',1,numel(num2str(LT))+ numel(num2str(Ltrem)) +39),'%.0f secs. Computation time left: %.0f minutes.'],t(i),trem);
      end
      LT = round(t(i),0);
      Ltrem = round(trem,0);
  end
end
fprintf('\nFinished\n');
