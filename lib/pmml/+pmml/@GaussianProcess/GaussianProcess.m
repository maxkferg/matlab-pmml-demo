classdef GaussianProcess < handle
% GaussianProcess. Gaussian Process Regression model for PMML package
%   This object represents a trained GPR model. It can be saved to a PMML
%   file using the toPMML() method. A similar object can be created from
%   a PMML file by passing the filename to the constructor.
%
%   Internally this class depends on the GPML matlab package for scoring.
%   The hyperparameters are always passed/stored in the form defined by the
%   GPML package to maintain compatibility.
    
   
    properties (Access=public)
        hyp;  % Struct containing hyperparameters
    end
    
    properties (Access=private)
        inferenceFunc;  % Function object from gpml
        meanFunc;       % Function object from gpml
        covFunc;        % Function object from gpml
        likFunc;        % Function object from gpml
        xTrain;         % Array of values
        yTrain;         % Array of values
    end
    
    methods
        
        function self = GaussianProcess(hyp, inf, mean, cov, lik, x, y)
            % Create a new GP object from either a PMML file or input parameters
            % This class can be initialized in two ways:
            %    GaussianProccess(hyp, inf, mean, cov, lik, x, y)
            %    GaussianProccess(filename)
            %
            % Where:
            %    hyp       struct of column vectors of mean/cov/lik hyperparameters
            %    inf       function specifying the inference method 
            %    mean      mean kernel function from ['MeanZero']
            %    cov       covariance kernel function from ['RadialBasisKernel','ARDSquaredExponentialKernel','AbsoluteExponentialKernel','GeneralizedExponentialKernel']
            %    lik       likelihood function from ['Gaussian']
            %    x         n by D matrix of training inputs
            %    y         column vector of length n of training targets
            %    filename  The path to a valid pmml file
            
            % Type checking for first parameter
            if (nargin==1) && ~isa(hyp,'char')
                throw 'Filename should be a string'
            elseif (nargin>1) && isa(hyp,'char')
                throw 'Hyperparameters should not be a string'
            end
            
            if (nargin==1)
                % There is only one input parameters so the object is 
                % being loaded from file.
                [hyp, inf, mean, cov, lik, x, y] = self.fromFile(hyp);
            end
            
            % Store the important parameters
            self.hyp = hyp;
            self.inferenceFunc = inf;
            self.meanFunc = mean;
            self.covFunc = cov;
            self.likFunc = lik;
            self.xTrain = x;
            self.yTrain = y;
            
            % Do some naive input validation
            self.getMeanFunc();
            self.getCovFunc();
            self.getLikFunc();
            self.getInferenceFunc();
        end 
            
          
        function pmml = toPMML(self,filename)
            % Translate a trained gpml object into PMML
            % Only supports gpml objects trained with ARD squared exponential
            % cov function and zero mean function. 
            %
            % @param(String) filename. Filename for output PMML
            % @output(String) pmml. Return pmml as a string
            infk = self.inferenceFunc;
            meank = self.meanFunc;
            covk = self.covFunc;
            lik = self.likFunc;
            x = self.xTrain;
            y = self.yTrain;
         
            pmml = self.translate();
            if (nargin>0)
                % Save string to file
                fid = fopen(filename,'wt');
                fprintf(fid, pmml);
                fclose(fid);
            end
        end
        
        
        function hyp = optimize(self,n)
            % Optimize the hyperparameters of this model
            % Uses the training values that were supplied at initialization
            % @param{int} n. The maximum number of optimization steps
            if (nargin<2)
                n = -100;
            end     
            meanfunc = self.getMeanFunc();  % Zero mean function
            covfunc = self.getCovFunc();   % ARD Squared exponential cov function
            likfunc = self.getLikFunc();
            infer = self.getInferenceFunc();
            x = self.xTrain;
            y = self.yTrain;
            hyp = minimize(self.hyp, @gp, n, infer, meanfunc, covfunc, likfunc, x, y);
            self.hyp = hyp;
        end
       
        
        function [yPredict,sPredict] = score(self,xNew)
        % Score. Score new x values using GaussianProcess Regression
        %   @param(array)  xNew. New x values to score where each row is a test point
        %   @output(array) yPredict. The predicted y values
        %   @output(array) sPredict. The predicted standard deviation

            % Check that xNew has the same number of columns as the training values
            % xNew should be stored as rows of training points
            if (size(xNew,2)~=size(self.xTrain,2))
                error('xNew must have %i columns',size(self.xTrain,2));
            end

            meanfunc = self.getMeanFunc();  % Zero mean function
            covfunc = self.getCovFunc();   % ARD Squared exponential cov function
            likfunc = self.getLikFunc();
            infer = self.getInferenceFunc();
            x = self.xTrain;
            y = self.yTrain;
            [yPredict,sPredict] = gp(self.hyp, infer, meanfunc, covfunc, likfunc, x, y, xNew);
        end

        
        function gridSearch(self,hyp,xParam,xRange,yParam,yRange)
        % Perform a grid search over the covariance hyperparameters and 
        % plot the results
        % @param{Int} xParam. The index of the hyp.cov parameter to change
        % @param{Array} xRange. The range of values to try
        % @param{Int} yParam. The index of the hyp.cov parameter to change
        % @param{Array} yRange. The range of values to try     
            meanfunc = self.getMeanFunc();
            covfunc = self.getCovFunc();
            likfunc = self.getLikFunc();
            infer = self.getInferenceFunc();
            x = linspace(xRange(1),xRange(2),100);
            y = linspace(yRange(1),yRange(2),100);
            r = zeros(length(x),length(y));
            for i=1:length(x)
                for j=1:length(y)
                    hyp.cov(xParam) = x(i);
                    hyp.cov(yParam) = y(j);
                    r(i,j) = gp(hyp, infer, meanfunc, covfunc, likfunc, self.xTrain, self.yTrain);
                end
            end
            contour(x,y,r);
            colorbar;
        end
    end
    
    
    
    methods (Access=private)
           
        function [hyper,inf,mean,cov,lik,x,y] = fromFile(self,filename)
        % Setup the object by loading parameters from a pmml file
        % Also validates the parameters were loaded correctly
            [hyper,inf,mean,cov,lik,x,y] = self.parse(filename);
        end
        
        function func = getMeanFunc(self)
        % Return the mean kernal function. We only support the MeanZero
        % function according to the PMML specification
            functionName = self.meanFunc;
            if strcmp(functionName,'MeanZero')
                func = @meanZero;
            else 
                throw(['Unknown mean function ', functionName]);
            end
        end
        
        
        function func = getCovFunc(self)
        % Return the covariance kernal function.
        % The PMML documentation refers to this function as the 
        % 'kernal function' although we will adopt the name covariance
        % function to distinguish between and mean and cov kernel
            functionName = self.covFunc;
            if strcmp(functionName,'RadialBasisKernel')
                func = @unsure;
            elseif strcmp(functionName,'ARDSquaredExponentialKernel')
                func = @covSEard; 
            elseif strcmp(functionName,'AbsoluteExponentialKernel')
                func = @unsure;
            elseif strcmp(functionName,'GeneralizedExponentialKernel')
                func = @unsure;
            else 
                throw('Unknown covariance function ' + functionName);
            end
        end
        
        
        function func = getLikFunc(self)
        % Return the likilihood function. We only support the Gaussian
        % likilihood function according to the PMML specification
            functionName = self.likFunc;
            if strcmp(functionName,'Gaussian')
                func = @likGauss;
            else 
                throw(['Unsupported likilihood function ', functionName]);
            end
        end
        
        
        function func = getInferenceFunc(self)
        % Return the inference function. We only support the 'Exact'
        % inference function.
            functionName = self.inferenceFunc;
            if strcmp(functionName,'Exact')
                func = @infExact;
            else 
                throw(['Unsupported inference function ', functionName]);
            end
        end
    end
end

