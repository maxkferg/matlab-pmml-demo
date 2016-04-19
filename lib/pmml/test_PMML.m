function test_PMML()
% test_PMML Some simple tests for the PMML class
%   Test that the PMML class can create valid PMML
%   Test that the PMML class can read pmml and score new values
%   Test that that the PMML class can score without read/writing
    clear all; clc;
    addpath(genpath('lib'));
    addpath(genpath('test'));
     
    testWriting();
    testReadThenScore();
    testScoring();
    testOptimize();
end



function testWriting()
% Test that the PMML class can create valid PMML
% Would be nice to have a xsd to test this
    expected = 'test/fixtures/expected.pmml';
    filename = 'test/fixtures/output.pmml';

    % Define valid function inputs matching the documentation example
    % The hyperparameters are defined in the same way that gpml returns them
    % This make the PMML package easier to use with gpml, but requires the
    % PMML package to make conversions internally
    sn = 0.1051; lambda1=1.5164; lambda2=59.3113; gamma=sqrt(2.4890);
    hyp.lik = log(sn);
    hyp.mean = [];
    hyp.cov = log([lambda1; lambda2; gamma]);
    meanfunc = 'MeanZero';
    covfunc = 'ARDSquaredExponentialKernel';
    likfunc = 'Gaussian';
    inffunc = 'Exact';
    xTrain = [1,3; 2,6];
    yTrain = [1; 2];

    p = pmml.GaussianProcess(hyp, inffunc, meanfunc, covfunc, likfunc, xTrain, yTrain);
    p.toPMML(filename);

    % Compare output to expected output
    expect = fopen(expected);
    actual = fopen(filename);
    while 1
        eline = fgetl(expect);
        aline = fgetl(actual);
        if ~ischar(eline)
            break
        end
        assert(strcmp(eline,aline),sprintf('%s does not match %s',eline,aline));
    end
    fclose(expect); fclose(actual);
    fprintf('GP Test: testWriting passed\n');
end



function testReadThenScore()
% Test that the PMML class can read pmml and score new values
    xNew = [1,4];
    filename = 'test/fixtures/expected.pmml';

    % Load model from PMML file
    model = pmml.GaussianProcess(filename);

    % Score the example values
    [mu,s] = model.score(xNew);

    testPrediction(mu,s);
    fprintf('GP Test: testScoring passed\n');
end



function testScoring()
% Test that that the PMML class can score without read/writing
% Complete the nist example using the gpml package
    xTrain = [1,3; 2,6];
    yTrain = [1; 2];
    xNew = [1,4];

    % Define mean and cov function
    likfunc = @likGauss;
    meanfunc = @meanZero;  % Zero mean function
    covfunc = @covSEard;   % ARD Squared exponential cov function
    gamma = sqrt(3);       % Realistic starting value for gamma
    lambda1 = 2;           % Realistic starting value for lambda1
    lambda2 = 60;          % Realistic starting value for lambda2

    sn = 0.1; % sigma

    hyp.lik = log(sn);
    hyp.mean = [];
    hyp.cov = log([lambda1; lambda2; gamma]);

    % Optimize hyperparameters
    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, xTrain, yTrain);

    % Represent the trained model as a PMML.GaussianProcess object
    model = pmml.GaussianProcess(hyp, 'Exact', 'MeanZero', 'ARDSquaredExponentialKernel', 'Gaussian', xTrain, yTrain);

    % Test that we are predicting the right values
    [mu,s] = model.score(xNew);
    testPrediction(mu,s);
    fprintf('GP Test: testScoring passed\n');
end


function testOptimize()
% Test that the GaussianProcess object can optimize it's own
% Hyperparameters.
    xTrain = [1,3; 2,6];
    yTrain = [1; 2];
    xNew = [1,4];

    % Define mean and cov function
    likfunc = @likGauss;
    meanfunc = @meanZero;  % Zero mean function
    covfunc = @covSEard;   % ARD Squared exponential cov function
    gamma = sqrt(3);       % Realistic starting value for gamma
    lambda1 = 2;           % Realistic starting value for lambda1
    lambda2 = 60;          % Realistic starting value for lambda2

    sn = 0.1; % sigma

    hyp.lik = log(sn);
    hyp.mean = [];
    hyp.cov = log([lambda1; lambda2; gamma]);

    % Represent the trained model as a PMML.GaussianProcess object
    model = pmml.GaussianProcess(hyp, 'Exact', 'MeanZero', 'ARDSquaredExponentialKernel', 'Gaussian', xTrain, yTrain);

    % Optimize the hyperparameters
    hyp = model.optimize();
    
    % Test that the hyperparameters are correct
    testHyperparameters(hyp);
    
    % Test that we are predicting the right values
    [mu,s] = model.score(xNew);
    testPrediction(mu,s);
    fprintf('GP Test: testOptimize passed\n');
end







