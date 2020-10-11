function checkGradientUse(g, scaleFactor, gMaxBruker)
% checkGradeintUse - checks percentage gradient strength used, trows a
% warning if 100% is exceded
g = g*scaleFactor./gMaxBruker;
gVec = g(:);

maxGradientUsed = max(abs(gVec))*100;
fprintf('Max gradient strength used: %f %% \n',maxGradientUsed)
if maxGradientUsed > 100
    warning("Excedding 100% gradient strength for the desired B_value")                
end
end