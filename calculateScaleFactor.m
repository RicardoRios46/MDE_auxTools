function scaleFactor =  calculateScaleFactor(g, dt, desiredBvalue, tol)
% calculateScaleFactor - Calculate a Scale factor for a waveform to scale it
% to a desired B value
step = .001;
scaleFactor = 1;
[~, bValue] = calculateBTensor(g, dt);

while not(desiredBvalue-tol < bValue && bValue < desiredBvalue+tol)
    if bValue < desiredBvalue
        scaleFactor = scaleFactor+step;
        gScaled = g*scaleFactor;
        [~, bValue] = calculateBTensor(gScaled, dt);
    else %b_val > desiredB
        scaleFactor=scaleFactor-step;
        gScaled = g*scaleFactor;
        [~, bValue] = calculateBTensor(gScaled, dt);
    end
end

end