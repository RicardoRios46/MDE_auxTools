function [B, b_val] = get_Btensor(g, dt)
% g: gradient waveform, n*3, in T/m
% d: step time for each point in gradient waveformt, in s

gamma = 42.6e6*2*pi;

% integrating from 0 to t each time
x = 0:dt:dt*size(g,1);
for i = 1:size(g,1)-1
    x_interval = x(1:i+1);
    g_interval = g(1:i+1,:);
    q(i,:) = trapz(x_interval,g_interval);
end

q = gamma*q; %SI-units

% ToDo why do I have to do the integral this way?
B = dt*(q'*q); %s/m^2

b_val = trace(B)*1e-6; % s/mm^2

end