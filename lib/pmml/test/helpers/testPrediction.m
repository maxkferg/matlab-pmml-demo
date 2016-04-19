function testPrediction(muPredict,sdPredict)
% testPrediction. Assert the predicted values (scores) match documented values
    % Throws an error on failure. Prints total error on success
    % See <documentationUrl> for an example case which generates these scores
    muTol = 0.01; % 1 percent tolerance
    sdTol = 0.10; % 10 percent tolerance (example has high round-off errors)
    muActual = 1.0095;
    sdActual = 0.0226;
    
    muError = sum(abs(muActual-muPredict)./muActual);
    sdError = sum(abs(sdActual-sdPredict)./sdActual);
    
    if muError > muTol
        error('Predicted mu value does not match expected value');
    end
    if sdError > sdTol
        error('Predicted sd value does not match expected value');
    end
    fprintf('Predicted values are within %.1f percent\n',100*max(muError,sdError));
end