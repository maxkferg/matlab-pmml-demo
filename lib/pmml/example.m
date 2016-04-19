function example()
% Train a GP Regression classifier, write it to PMML, and then use it to
% score values
    clear all; clc;
    addpath(genpath('lib'));
    addpath(genpath('test'));
    filename = 'test/fixtures/output.pmml';

    % Define valid function inputs matching the documentation example
    % The hyperparameters are defined in the same way that gpml returns them
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

    % Train and optimize a GP model
    p = pmml.GaussianProcess(hyp, inffunc, meanfunc, covfunc, likfunc, xTrain, yTrain);
    p.optimize(-100);
    
    % Write the model to a PMML file
    p.toPMML(filename);

    % Load model from PMML file
    model = pmml.GaussianProcess(filename);

    % Score some new values
    xNew = [1,4];

    % Score the example values
    [mu,s] = model.score(xNew);

    testPrediction(mu,s);
    fprintf('GP Test: testScoring passed\n');
end

