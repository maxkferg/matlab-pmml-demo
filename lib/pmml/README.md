# Gaussian Process Regression PMML Support for Matlab

Save and load a trained Gaussian Process regression (GPR) model to/from PMML. This package exposes the
`pmml.GaussianProcess` class which is used to represent a trained GPR model. The model hyperparameters
can be optimized using the GPML package, or directly on any `GaussianProcess` object.
`GaussianProcess` objects can be used to generate scores for new x values, regardless of whether they
were initialized from a PMML file, or a trained GPML model.

## Creating GaussianProcess objects

### pmml.GaussianProcess(hyperparameters, infFunc, meanFunc, covFunc, likFunc, xTrain, yTrain)
Create a new GaussianProcess object from either a PMML file or input parameters.

Where:
* hyp        struct of column vectors of mean/cov/lik hyperparameters
* infFunc    function specifying the inference method
* meanFunc   mean kernel function from ['MeanZero']
* covFunc    covariance kernel function
* likFunc    likelihood function from ['Gaussian']
* x          n by D matrix of training inputs
* y          column vector of length n of training targets

The covariance kernel function can be any valid function name,
as specified in the PMML standard e.g. ['RadialBasisKernel','ARDSquaredExponentialKernel','AbsoluteExponentialKernel','GeneralizedExponentialKernel']

The hyp parameter should take the same form as used by the GPML package.
* hyp.lik - The log of the noise variance. log(sigma_n)
* hyp.mean - parameters for the mean kernel (empty matrix)
* hyp.cov - parameters for the cov kernel. log([lambda1; lambda2; gamma])


### pmml.GaussianProcess(filename)
Create a new GaussianProcess object from an existing PMML file.
This method of creating GaussianProcess objects is used to load trained models from PMML.

Where:
* filename - the path to a valid PMML filename

## Object methods
Once a GaussianProcess object has been created it can be used to score new
x values or it can be saved to a PMML file. For this section, assume that
`p` is a valid GaussianProcess object.

### p.optimize()
Optimize the the hyperparameters of this model, for the training
values kernel type passed at initialization.

### p.score(xNew)
Return scores for the new x values. xNew should be an m x n matrix of values
where each row represents a test point. The method will return an m x 1
column vector of y values (scores).

### p.toPMML(filename)
Return the trained GPR model as valid PMML. If the optional filename
parameter is provided, the PMML will be saved to file.




## Example

```matlab
    % Define valid function inputs matching the documentation example
    % The hyperparameters are defined in the same way that gpml returns them
    % This make the PMML package easier to use with gpml, but requires the
    % PMML package to make conversions internally
    sn = 0.1051;
    lambda1=1.5164;
    lambda2=59.3113;
    gamma=sqrt(2.4890);
    hyp.lik = log(sn);
    hyp.mean = [];
    hyp.cov = log([lambda1; lambda2; gamma]);

    meanfunc = 'MeanZero';
    covfunc = 'ARDSquaredExponentialKernel';
    likfunc = 'Gaussian';
    inffunc = 'Exact';
    xTrain = [1,3; 2,6];
    yTrain = [1; 2];

    % Create a GPR model
    p = pmml.GaussianProcess(hyp, inffunc, meanfunc, covfunc, likfunc, xTrain, yTrain);

    % Optimize the hyperparameters
    p.optimize();

    % Access the hyperparameters
    p.hyp

    % Score some new values
    p.score([1,4])

    % Save the pmml model
    p.toPMML('output.pmml');
```
The GPR model and training points are now saved in the PMML format.
The model can be loaded and used to score some new values.

```matlab
	% Load the GPR model from file
	p = pmml.GaussianProcess('output.pmml');

	% Score the new values
    p.score([1,4])

    % Score multiple trianing points
    p.score([1,4; 2,3; 0,3])
```

## GPML Support
This package is designed to work flawlessly with the GPML package. GPML objects can easily be converted
to PMML by padding it's hyperparameters to the `pmml.GaussianProcess` class.

## TODO
- Support for 'RadialBasisKernel' and 'GeneralizedExponentialKernel' in Matlab package
- Test with a large number of inputs
- Support column naming (other than 'x1','x2',...)

## License
MIT

