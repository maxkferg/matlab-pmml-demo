function testHyparameters(hypOpt)
% testHyparameters. Assert the hyperparameters match documented values
    % Throws an error on failure. Prints total error on success
    % See <documentationUrl> for an example case which generates these hyperparams
    gammaStar = 2.4890;
    lambda1 = 1.5164;
    lambda2 = 59.3113;
    
    tol = 0.01; % 1 percent tolerance
    expected = [log(lambda1) log(lambda2) log(sqrt(gammaStar))];
    error = sum(abs(hypOpt.cov-expected')./expected');
    
    if error > tol
        throw 'Hyperparameters do not match expected values'
    end
    fprintf('Hyperparameters are within %.1f percent\n',100*error);
end