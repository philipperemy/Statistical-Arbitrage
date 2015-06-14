function [h,pValue,stat,cValue,reg] = adftest(y,varargin)
%ADFTEST Augmented Dickey-Fuller test for a unit root
%
% Syntax:
%
%   [h,pValue,stat,cValue,reg] = adftest(y)
%   [h,pValue,stat,cValue,reg] = adftest(y,param1,val1,param2,val2,...)
%
% Description:
%
%   Dickey-Fuller tests assess the null hypothesis of a unit root in a
%   univariate time series y. All tests use the model
%
%       y(t) = c + d*t + a*y(t-1) + b1*(1-L)y(t-1)
%                                 + b2*(1-L)y(t-2)
%                                 + ... 
%                                 + bp*(1-L)y(t-p)
%                                 + e(t),
%
%   where L is the lag operator Ly(t) = y(t-1). The null hypothesis
%   restricts a = 1. Variants of the test, appropriate for series with
%   different growth characteristics, restrict the drift and deterministic
%   trend coefficients, c and d, respectively, to be 0. Lagged differences
%   bk*(1-L)y(t-k), k = 1, ..., p, "augment" the test to account for serial
%   correlations in the innovations process e(t).
%
% Input Arguments:
%
%   y - Vector of time-series data. The last element is the most recent
%       observation. NaNs indicating missing values are removed.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME        VALUE
%
%   'lags'      Scalar or vector of nonnegative integers indicating the
%               number p of lagged changes of y to include in the model.
%               The default value is 0.
%
%   'model'     String or cell vector of strings indicating the model
%               variant. Values are 'AR' (autoregressive), 'ARD'
%               (autoregressive with drift), or 'TS' (trend stationary).
%               The default value is 'AR'.
%
%               o When the value is 'AR', the null model
%
%                   y(t) = y(t-1) + b1*(1-L)y(t-1)
%                                 + b2*(1-L)y(t-2)
%                                 + ... 
%                                 + bp*(1-L)y(t-p)
%                                 + e(t)                   
%
%                 is tested against the alternative model
%
%                   y(t) = a*y(t-1) + b1*(1-L)y(t-1)
%                                   + b2*(1-L)y(t-2)
%                                   + ... 
%                                   + bp*(1-L)y(t-p)
%                                   + e(t)
%
%                 with AR(1) coefficient a < 1.
%
%               o When the value is 'ARD', the null model
%
%                   y(t) = y(t-1) + b1*(1-L)y(t-1)
%                                 + b2*(1-L)y(t-2)
%                                 + ... 
%                                 + bp*(1-L)y(t-p)
%                                 + e(t)                   
%
%                 is tested against the alternative model
%
%                   y(t) = c + a*y(t-1) + b1*(1-L)y(t-1)
%                                       + b2*(1-L)y(t-2)
%                                       + ... 
%                                       + bp*(1-L)y(t-p)
%                                       + e(t)
%
%                 with drift coefficient c and AR(1) coefficient a < 1.
%
%               o When the value is 'TS', the null model
%
%                   y(t) = c + y(t-1) + b1*(1-L)y(t-1)
%                                     + b2*(1-L)y(t-2)
%                                     + ... 
%                                     + bp*(1-L)y(t-p)
%                                     + e(t)                   
%
%                 is tested against the alternative model
%
%                   y(t) = c + d*t + a*y(t-1) + b1*(1-L)y(t-1)
%                                             + b2*(1-L)y(t-2)
%                                             + ... 
%                                             + bp*(1-L)y(t-p)
%                                             + e(t)
%
%                 with drift coefficient c, deterministic trend coefficient
%                 d, and AR(1) coefficient a < 1. 
%
%   'test'      String or cell vector of strings indicating the type of
%               test statistic. Values are 't1', 't2', or 'F'. The default
%               value is 't1'. 
%
%               o When the value is 't1', a standard t statistic
%
%                   t1 = (a-l)/se
%
%                 is computed from OLS estimates of the AR(1) coefficient a
%                 and its standard error se in the alternative model. The
%                 test assesses the significance of the restriction a = 1.
%
%               o When the value is 't2', a lag-adjusted, "unstudentized" t
%                 statistic
%
%                   t2 = T*(a-1)/(1-b1-...-bp)
%
%                 is computed from OLS estimates of the AR(1) coefficient a
%                 and the stationary coefficients b1, ..., bp in the
%                 alternative model. T is the effective sample size,
%                 adjusted for lags and missing values. The test assesses
%                 the significance of the restriction a = 1.
%
%               o When the value is 'F', an F statistic is computed to
%                 assess the significance of a joint restriction on the
%                 alternative model. If the value of 'model' is 'ARD', the
%                 restriction is a = 1 and c = 0. If the value of 'model'
%                 is 'TS', the restriction is a = 1 and d = 0. An F test is
%                 invalid when the value of 'model' is 'AR'.
%
%   'alpha'     Scalar or vector of nominal significance levels for the
%               tests. Values must be between 0.001 and 0.999. The default
%               value is 0.05.
%
%   Scalar or single string parameter values are expanded to the length of
%   any vector value (the number of tests). Vector values must have equal
%   length. If any value is a row vector, all outputs are row vectors.
%
% Output Arguments:
%
%   h - Vector of Boolean decisions for the tests, with length equal to the
%       number of tests. Values of h equal to 1 indicate rejection of the
%       unit-root null in favor of the alternative model. Values of h equal
%       to 0 indicate a failure to reject the unit-root null.
%
%   pValue - Vector of p-values of the test statistics, with length equal
%       to the number of tests. When the value of 'test' is 't1' or 't2',
%       p-values are left-tail probabilities. When the value of 'test' is
%       'F', p-values are right-tail probabilities.
%
%   stat - Vector of test statistics, with length equal to the number of
%       tests. Statistics are computed using OLS estimates of the
%       coefficients in the alternative model.
%
%   cValue - Vector of critical values for the tests, with length equal to
%       the number of tests. When the value of 'test' is 't1' or 't2',
%       critical values are for left-tail probabilities. When the value of
%       'test' is 'F', critical values are for right-tail probabilities.
%
%   reg - Structure of regression statistics from the OLS estimation of 
%       coefficients in the alternative model. The number of records is
%       equal to the number of tests. Each record has the following fields:
%
%       num         Length of the input series y, with NaNs removed
%       size        Effective sample size, adjusted for lags, difference*
%       names       Regression coefficient names			
%       coeff       Estimated coefficient values
%       se          Estimated coefficient standard errors
%       Cov         Estimated coefficient covariance matrix
%       tStats      t statistics of coefficients and p-values
%       FStat       F statistic and p-value
%       yMu         Mean of y, adjusted for lags, difference*
%       ySigma      Standard deviation of y, adjusted for lags, difference*
%       yHat        Fitted values of y, adjusted for lags, difference*
%       res         Regression residuals
%       DWStat      Durbin-Watson statistic
%       SSR         Regression sum of squares
%       SSE         Error sum of squares
%       SST         Total sum of squares
%       MSE         Mean squared error
%       RMSE        Standard error of the regression
%       RSq         R^2 statistic
%       aRSq        Adjusted R^2 statistic
%       LL          Loglikelihood of data under Gaussian innovations
%       AIC         Akaike information criterion
%       BIC         Bayesian (Schwarz) information criterion
%       HQC         Hannan-Quinn information criterion
%
%       *Lagging and differencing a time series reduce the sample size.
%       Absent any presample values, if y(t) is defined for t = 1:N, then
%       the lagged series y(t-k) is defined for t = k+1:N. Differencing
%       reduces the time base to k+2:N. With p lagged differences, the
%       common time base is p+2:N and the effective sample size is N-(p+1).
%
% Notes:
%
%   o A suitable value for 'lags' must be determined in order to draw valid
%     inferences from the test. One method is to begin with a maximum lag,
%     such as the one recommended by Schwert [7], and then test down by
%     assessing the significance of the coefficient of the largest lagged
%     change in y, bp. The usual t statistic is appropriate, as reported in
%     the reg output structure. Another method is to combine a measure of
%     fit, such as SSR, with information criteria such as AIC, BIC, and
%     HQC. These statistics are also reported in the reg output structure.
%     Ng and Perron [6] provide further guidelines.
%
%   o The value of 'model' is determined by the growth characteristics of
%     the time series being tested, and should be chosen with a specific
%     testing strategy in mind. As discussed in Elder & Kennedy [4],
%     including too many regressors results in lost power, while including
%     too few biases the test in favor of the null. In general, if a series
%     is growing, the 'TS' model provides a reasonable trend-stationary
%     alternative to a unit-root process with drift. If a series is not
%     growing, 'AR' and 'ARD' models provide reasonable stationary
%     alternatives to a unit-root process without drift. The 'ARD'
%     alternative has mean c/(1-a); the 'AR' alternative has mean 0.
%
%   o Dickey-Fuller statistics follow nonstandard distributions under the
%     null, even asymptotically. Critical values for a range of sample
%     sizes and significance levels have been tabulated using Monte Carlo
%     simulations of the null model with Gaussian innovations and five
%     million replications per sample size. For small samples, values are
%     valid only for Gaussian innovations; for large samples, values are
%     also valid for non-Gaussian innovations. Critical values and p-values
%     are interpolated from the tables. Tables for tests of type 't1' and
%     't2' are identical to those for PPTEST.
%
% Example:
%
%   % Test GDP data for a unit root using a trend-stationary alternative
%   % with 0, 1, and 2 lagged differences:
%
%   load Data_GDP
%   y = log(Data);
%   h = adftest(y,'model','TS','lags',0:2)
%
%   % The test fails to reject the unit-root null with each alternative.
%
% References:
%
%   [1] Davidson, R. and J. G. MacKinnon. Econometric Theory and Methods.
%       Oxford, UK: Oxford University Press, 2004.
%  
%   [2] Dickey, D. A., and W. A. Fuller. "Distribution of the Estimators
%       for Autoregressive Time Series with a Unit Root." Journal of the
%       American Statistical Association. Vol. 74, 1979, pp. 427-431.
% 
%   [3] Dickey, D. A., and W. A. Fuller. "Likelihood Ratio Statistics for
%       Autoregressive Time Series with a Unit Root." Econometrica. Vol.
%       49, 1981, pp. 1057-1072.
%
%   [4] Elder, J., and P. E. Kennedy. "Testing for Unit Roots: What Should
%       Students Be Taught?" Journal of Economic Education. Vol. 32, 2001, 
%       pp. 137-146.
% 
%   [5] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
%   [6] Ng, S., and P. Perron. "Unit Root Tests in ARMA Models with
%       Data-Dependent Methods for the Selection of the Truncation Lag."
%       Journal of the American Statistical Association. Vol. 90, 1995, 
%       pp. 268-281.
%
%   [7] Schwert, W. "Tests for Unit Roots: A Monte Carlo Investigation."
%       Journal of Business and Economic Statistics. Vol. 7, 1989,
%       pp. 147-159.
%
% See also PPTEST, LMCTEST, KPSSTEST, VRATIOTEST. 

% Copyright 2009-2010 The MathWorks, Inc.

% Parse inputs and set defaults:

parseObj = inputParser;
parseObj.addRequired('y',@yCheck);
parseObj.addParamValue('lags',0,@lagsCheck);
parseObj.addParamValue('model','AR',@modelCheck);
parseObj.addParamValue('test','t1',@testCheck);
parseObj.addParamValue('alpha',0.05,@alphaCheck);

parseObj.parse(y,varargin{:});

y = parseObj.Results.y;
lags = parseObj.Results.lags;
model = parseObj.Results.model;
test = parseObj.Results.test;
alpha = parseObj.Results.alpha;

% Check parameter values for commensurate lengths, expand scalars and
% single strings, and convert all variables to columns:

[numTests,rowOutput,lags,model,test,alpha] = sizeCheck(lags,model,test,alpha);
y = y(:);

% Adjust y for missing, lagged, and differenced values:

y(isnan(y)) = []; % Remove missing values
N = length(y);    % Size of observed data
T = N-(lags+1);   % Effective sample size 

% Set row/column grid for tables of critical values:

sampSizes = [10 15 20 25 30 40 50 75 100 150 200 300 500 1000 10000];
minT = min(sampSizes);              % Minimum effective sample size
maxT = max(max(sampSizes),max(T));  % Maximum effective sample size
sampSizes(end) = maxT;              % Force maxT into table

sigLevels = [0.001 (0.005:0.005:0.10) (0.125:0.025:0.20) ...
                   (0.80:0.025:0.875) (0.90:0.005:0.995) 0.999];
minAlpha = min(sigLevels); % Minimum significance level
maxAlpha = max(sigLevels); % Maximum significance level

% Check if any tests are outside of tables:

if any(T < minT)
    
    error(message('econ:adftest:EffectiveSampleSizeLessThanTabulatedValues', minT))
      
end

if any(alpha < minAlpha) || any(alpha > maxAlpha)

    error(message('econ:adftest:AlphaOutOfRange', sprintf( '%5.3f', minAlpha ), sprintf( '%5.3f', maxAlpha )))
      
end

% Preallocate output variables:

switch nargout
    
    case {0,1}
        
        needPValue = false;
        needRegOut = false;
        
        h = false(numTests,1);
        
    case 2
        
        needPValue = true;
        needRegOut = false;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        
    case 3
        
        needPValue = true;
        needRegOut = false;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        stat = pValue;
        
    case 4
        
        needPValue = true;
        needRegOut = false;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        stat = pValue;
        cValue = pValue;
        
    case 5
        
        needPValue = true;
        needRegOut = true;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        stat = pValue;
        cValue = pValue;
        regFields = {'num','size','names','coeff','se','Cov','tStats',...
                     'FStat','yMu','ySigma','yHat','res','DWStat','SSR',...
                     'SSE','SST','MSE','RMSE','RSq','aRSq','LL','AIC',...
                     'BIC','HQC'};
        reg = cell2struct(cell(length(regFields),numTests),regFields,1);

end

% Initialize loop variables:

lastLags = [];
lastModel = [];
lastType = [];

% Run the tests:

for i = 1:numTests
    
    testT = T(i);
    testLags = lags(i);
    testModel = model{i};
    testType = test{i};
    testAlpha = alpha(i);
    
    if strcmpi(testModel,'AR') && strcmpi(testType,'F')
        
        error(message('econ:adftest:InvalidFTest'))
          
    end
    
    % Check to see if statistic is unchanged from last test:
    
    sameStat = isequal(testLags,lastLags) && ...
               isequal(testModel,lastModel) && ...
               isequal(testType,lastType);
    
    if ~sameStat % Recompute the statistic only if it changes
    
        % Perform the regression:
    
        testReg = runReg(i,y,testT,testLags,testModel,needRegOut);
    
        % Get appropriate table of critical values:
    
        CVTable = getCVTable(testModel,testType);
    
        % Compute the statistic:
    
        [testStat,testPValue] = getStat(i,testT,testLags,testModel,testType,testReg,sigLevels,sampSizes,CVTable,needPValue);
    
    end
    
    % Test the statistic:
    
    [testCValue,testH] = runTest(sigLevels,sampSizes,testT,testType,testAlpha,testStat,CVTable);
    
    % Add the test results to the outputs:
    
    switch nargout
    
        case {0,1}
        
            h(i) = testH;
        
        case 2
        
            h(i) = testH;
            pValue(i) = testPValue;
        
        case 3
        
            h(i) = testH;
            pValue(i) = testPValue;
            stat(i) = testStat;
        
        case 4
        
            h(i) = testH;
            pValue(i) = testPValue;
            stat(i) = testStat;
            cValue(i) = testCValue;
        
        case 5
        
            h(i) = testH;
            pValue(i) = testPValue;
            stat(i) = testStat;
            cValue(i) = testCValue;
            reg(i) = testReg;

    end
    
    % Save values to check for changes, next loop:
    
    lastLags = testLags;
    lastModel = testModel;
    lastType = testType;

end

% Display outputs as row vectors if any parameter value is a row vector:

if rowOutput
    
    switch nargout
        
        case {0,1}

            h = h';

        case 2
            
            h = h';
            pValue = pValue';
            
        case 3
            
            h = h';
            pValue = pValue';
            stat = stat';
            
        case 4
            
            h = h';
            pValue = pValue';
            stat = stat';
            cValue = cValue';
            
        case 5
            
            h = h';
            pValue = pValue';
            stat = stat';
            cValue = cValue';
            reg = reg';
        
    end
    
end

%-------------------------------------------------------------------------
% Check input y
function OK = yCheck(y)
            
    if isempty(y)
        
        error(message('econ:adftest:DataUnspecified'))
          
    elseif ~isnumeric(y)
        
        error(message('econ:adftest:DataNonNumeric'))
          
    elseif ~isvector(y)
        
        error(message('econ:adftest:DataNonVector'))
          
    else
        
        OK = true;
        
    end

%-------------------------------------------------------------------------
% Check value of 'lags' parameter
function OK = lagsCheck(lags)
    
    if ~isnumeric(lags)
        
        error(message('econ:adftest:LagsNonNumeric'))
          
    elseif ~isvector(lags)
        
        error(message('econ:adftest:LagsNonVector'))
          
    elseif any(mod(lags,1) ~= 0) || any(lags < 0)
        
        error(message('econ:adftest:LagsOutOfRange'))
          
    else
        
        OK = true;
        
    end

%-------------------------------------------------------------------------
% Check value of 'model' parameter
function OK = modelCheck(model)
    
    if ~isvector(model)
        
        error(message('econ:adftest:ModelNonVector'))
          
    elseif isnumeric(model) || (iscell(model) && any(cellfun(@isnumeric,model)))
        
        error(message('econ:adftest:ModelNumeric'))
    
    elseif ~all(ismember(upper(model),{'AR','ARD','TS'}))
        
        error(message('econ:adftest:ModelInvalid'))
          
    else
        
        OK = true;
        
    end

%-------------------------------------------------------------------------
% Check value of 'test' parameter
function OK = testCheck(test)
    
    if ~isvector(test)
        
        error(message('econ:adftest:TypeOfTestNonVector'))
          
    elseif isnumeric(test) || (iscell(test) && any(cellfun(@isnumeric,test)))
        
        error(message('econ:adftest:TypeOfTestNumeric'))
    
    elseif ~all(ismember(upper(test),{'T1','T2','F'}))
        
        error(message('econ:adftest:TypeOfTestInvalid'))
          
    else
        
        OK = true;
        
    end

%-------------------------------------------------------------------------
% Check value of 'alpha' parameter
function OK = alphaCheck(alpha)
    
    if ~isnumeric(alpha)
        
        error(message('econ:adftest:AlphaNonNumeric'))
          
    elseif ~isvector(alpha)
        
        error(message('econ:adftest:AlphaNonVector'))
          
    else
        
        OK = true;
        
    end
 
%-------------------------------------------------------------------------
% Check parameter values for commensurate lengths, expand scalars and
% single strings, and convert all variables to columns
function [numTests,rowOutput,varargout] = sizeCheck(varargin)

% Initialize outputs:

numTests = 1;
rowOutput = false;

% Determine vector lengths, number of tests, row output flag:

for i = 1:nargin
    
    ivar = varargin{i};
    iname = inputname(i);
    
    if isnumeric(ivar) || iscell(ivar)
        paramLength.(iname) = length(ivar);
        if ~isscalar(ivar)
        	rowOutput = rowOutput || (size(ivar,1) == 1);
        end    
    else        
        paramLength.(iname) = 1;   % Single string
        varargin{i} = varargin(i); % Convert to cell        
    end
    
    numTests = max(numTests,paramLength.(iname));
    
end

% Check for commensurate vector lengths:

for i = 1:(nargin-1)
    iname = inputname(i);
    for j = (i+1):nargin
        jname = inputname(j);
        if (paramLength.(iname) > 1) && (paramLength.(jname) > 1) ...
            && (paramLength.(iname) ~= paramLength.(jname))
        
            error(message('econ:adftest:ParameterSizeMismatch', iname, jname))
              
        end        
    end
end

% Expand scalars and single strings:

for i = 1:nargin
    
    ivar = varargin{i};
    if paramLength.(inputname(i)) == 1
        varargout{i} = ivar(ones(numTests,1));
    else
        varargout{i} = ivar(:);  % Column output
    end
    
end

%-------------------------------------------------------------------------
% Perform the test regression
function testReg = runReg(i,y,testT,testLags,testModel,needRegOut)

% Construct lagged ys and changes:

N = length(y);
yLags = lagmatrix(y,0:(testLags+1));
testY = yLags((testLags+2):end,1);
deltaYLags = -diff(yLags,1,2);

% Set up the regression:

switch upper(testModel)
    
    case 'AR'
        
    % Model: y(t) = a*y(t-1) + b1*(1-L)y(t-1) + b2*(1-L)y(t-2) + ... + bp*(1-L)y(t-p) + e(t)
        
        numNonDelta = 1; % Number of non-stationary coefficients
        
        % Design matrix:
        
        X = [yLags((testLags+2):end,2),deltaYLags((testLags+2):end,2:end)];
        
        % Coefficient names:
        
        if testLags > 0
            names = cat(1,'a',strcat({'b'},num2str((1:testLags)','%-d')));
        else
            names = {'a'};
        end
        
    case 'ARD'
        
    % Model: y(t) = c + a*y(t-1) + b1*(1-L)y(t-1) + b2*(1-L)y(t-2) + ... + bp*(1-L)y(t-p) + e(t)
        
        numNonDelta = 2; % Number of non-stationary coefficients
        
        % Design matrix:
        
        X = [ones(N-(testLags+1),1),yLags((testLags+2):end,2),deltaYLags((testLags+2):end,2:end)];
        
        % Coefficient names:
        
        if testLags > 0
            names = cat(1,'c','a',strcat({'b'},num2str((1:testLags)','%-d')));
        else
            names = {'c';'a'};
        end
        
    case 'TS'
        
    % Model: y(t) = c + d*t + a*y(t-1) + b1*(1-L)y(t-1) + b2*(1-L)y(t-2) + ... + bp*(1-L)y(t-p) + e(t)
        
        numNonDelta = 3; % Number of non-stationary coefficients
        
        % Design matrix:
        
        X = [ones(N-(testLags+1),1),(1:testT)',yLags((testLags+2):end,2),deltaYLags((testLags+2):end,2:end)];
        
        % Coefficient names:
        
        if testLags > 0
            names = cat(1,'c','d','a',strcat({'b'},num2str((1:testLags)','%-d')));
        else
            names = {'c';'d';'a'};
        end

end

if size(X,1) < size(X,2)
    
    error(message('econ:adftest:RankDeficientDesignMatrix'));
      
end

% Run the regression:

[Q,R] = qr(X,0);
coeff = R\(Q'*testY);
numParams = length(coeff);
yHat = X*coeff;
res = testY-yHat;
SSE = res'*res;
dfe = testT-numParams;
MSE = SSE/dfe;
S = R\eye(numParams);
Cov = S*S'*MSE;

% Write results to the regression record:

if needRegOut % Compute all statistics

    yBar = mean(testY);
    regRes = yHat-yBar;
    SSR = regRes'*regRes;
    dfr = numParams-1;
    dft = testT-1;
    diffRes = diff(res);
    
    testReg.num = N;        
    testReg.size = testT;
    testReg.names = names;
    testReg.coeff = coeff;
    testReg.se = sqrt(diag(Cov));
    testReg.Cov = Cov;
    testReg.tStats.t = coeff./testReg.se;
    testReg.tStats.pVal = 2*(tcdf(-abs(testReg.tStats.t),dfe));    
    testReg.FStat.F = (SSR/dfr)/(SSE/dfe);
    testReg.FStat.pVal = 1-fcdf(testReg.FStat.F,dfr,dfe);
    testReg.yMu = yBar;
    testReg.ySigma = std(testY);
    testReg.yHat = yHat;
    testReg.res = res;
    testReg.DWStat = (diffRes'*diffRes)/SSE;
    testReg.SSR = SSR;
    testReg.SSE = SSE;
    testReg.SST = SSR+SSE;
    testReg.MSE = MSE;
    testReg.RMSE = sqrt(MSE);
    testReg.RSq = 1-SSE/testReg.SST;
    testReg.aRSq = 1-(SSE/testReg.SST)*(dft/dfe);
    testReg.LL = -normlike([0,sqrt(MSE)],res);
    testReg.AIC = 2*numParams-2*testReg.LL;
    testReg.BIC = numParams*log(testT)-2*testReg.LL;
    testReg.HQC = 2*numParams*log(log(testT))-2*testReg.LL;
                 
else % Compute only the statistics needed to run the test 
    
    testReg.names = names;
    testReg.coeff = coeff;
    testReg.Cov = Cov;
    
end

% Warn if the roots of 1 - b1*z - ... - bp*z^p are inside the unit circle,
% invalidating the test (Hamilton [5] pp. 517-518 and notes to Table 17.3):

lagPolyRoots = roots([-coeff(end:-1:numNonDelta+1);1]);

    if any(abs(lagPolyRoots) <= 1)
        
        warning(message('econ:adftest:InvalidStatistic', i))
            
    end

%-------------------------------------------------------------------------
% Get table of critical values
function CVTable = getCVTable(testModel,testType)

switch upper(testModel)
    
    case 'AR'
        
        switch upper(testType)

            case 'T1'

              % Row headers (effective sample sizes) are at the right of the table.
              % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
              % Alpha       0.001    0.005    0.010    0.015    0.020    0.025    0.030    0.035    0.040    0.045    0.050    0.055    0.060    0.065    0.070    0.075    0.080    0.085    0.090    0.095    0.100    0.125    0.150    0.175    0.200    0.800    0.825    0.850    0.875    0.900    0.905    0.910    0.915    0.920    0.925    0.930    0.935    0.940    0.945    0.950    0.955    0.960    0.965    0.970    0.975    0.980    0.985    0.990    0.995    0.999      T
              % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
                CVTable = [-3.9976  -3.1423  -2.7891  -2.5834  -2.4351  -2.3201  -2.2260  -2.1460  -2.0766  -2.0148  -1.9594  -1.9095  -1.8629  -1.8201  -1.7799  -1.7427  -1.7075  -1.6744  -1.6427  -1.6126  -1.5839  -1.4569  -1.3501  -1.2567  -1.1730   0.4848   0.5897   0.7056   0.8373   0.9913   1.0258   1.0621   1.0996   1.1392   1.1810   1.2254   1.2722   1.3230   1.3775   1.4361   1.5004   1.5720   1.6520   1.7437   1.8514   1.9818   2.1485   2.3836   2.7840   3.7494   % 10
                           -3.7444  -3.0262  -2.7138  -2.5289  -2.3954  -2.2909  -2.2044  -2.1305  -2.0659  -2.0083  -1.9560  -1.9084  -1.8648  -1.8239  -1.7858  -1.7501  -1.7161  -1.6841  -1.6537  -1.6248  -1.5970  -1.4738  -1.3690  -1.2773  -1.1950   0.4548   0.5589   0.6733   0.8024   0.9521   0.9856   1.0203   1.0568   1.0946   1.1352   1.1775   1.2225   1.2703   1.3218   1.3771   1.4378   1.5044   1.5790   1.6643   1.7627   1.8807   2.0310   2.2387   2.5818   3.3667   % 15
                           -3.6182  -2.9694  -2.6777  -2.5029  -2.3760  -2.2757  -2.1928  -2.1217  -2.0590  -2.0034  -1.9524  -1.9061  -1.8635  -1.8237  -1.7866  -1.7517  -1.7187  -1.6872  -1.6578  -1.6292  -1.6020  -1.4806  -1.3773  -1.2862  -1.2043   0.4428   0.5465   0.6604   0.7883   0.9363   0.9691   1.0033   1.0390   1.0768   1.1162   1.1577   1.2017   1.2485   1.2988   1.3529   1.4115   1.4756   1.5478   1.6301   1.7244   1.8376   1.9796   2.1742   2.4956   3.2078   % 20
                           -3.5489  -2.9338  -2.6541  -2.4850  -2.3620  -2.2652  -2.1845  -2.1150  -2.0538  -1.9992  -1.9496  -1.9043  -1.8626  -1.8236  -1.7870  -1.7523  -1.7200  -1.6892  -1.6598  -1.6319  -1.6052  -1.4850  -1.3824  -1.2917  -1.2102   0.4344   0.5376   0.6513   0.7788   0.9264   0.9589   0.9929   1.0284   1.0654   1.1042   1.1451   1.1885   1.2352   1.2846   1.3380   1.3957   1.4600   1.5303   1.6107   1.7031   1.8140   1.9527   2.1397   2.4475   3.1150   % 25
                           -3.5090  -2.9118  -2.6402  -2.4762  -2.3565  -2.2604  -2.1810  -2.1130  -2.0532  -1.9992  -1.9501  -1.9054  -1.8639  -1.8250  -1.7888  -1.7544  -1.7221  -1.6913  -1.6620  -1.6340  -1.6074  -1.4881  -1.3863  -1.2960  -1.2147   0.4294   0.5324   0.6457   0.7725   0.9192   0.9515   0.9853   1.0204   1.0575   1.0962   1.1370   1.1801   1.2262   1.2756   1.3288   1.3863   1.4495   1.5186   1.5973   1.6889   1.7965   1.9325   2.1164   2.4142   3.0575   % 30
                           -3.4538  -2.8850  -2.6232  -2.4629  -2.3456  -2.2530  -2.1756  -2.1085  -2.0488  -1.9957  -1.9476  -1.9034  -1.8624  -1.8244  -1.7889  -1.7552  -1.7233  -1.6932  -1.6645  -1.6368  -1.6106  -1.4925  -1.3910  -1.3014  -1.2202   0.4217   0.5246   0.6375   0.7640   0.9097   0.9422   0.9756   1.0109   1.0475   1.0859   1.1262   1.1690   1.2142   1.2631   1.3153   1.3718   1.4342   1.5024   1.5796   1.6692   1.7762   1.9088   2.0896   2.3777   2.9899   % 40
                           -3.4129  -2.8678  -2.6112  -2.4535  -2.3384  -2.2473  -2.1709  -2.1051  -2.0467  -1.9940  -1.9466  -1.9029  -1.8623  -1.8244  -1.7890  -1.7556  -1.7240  -1.6940  -1.6653  -1.6381  -1.6117  -1.4942  -1.3932  -1.3037  -1.2230   0.4188   0.5217   0.6346   0.7603   0.9061   0.9383   0.9716   1.0066   1.0428   1.0813   1.1214   1.1638   1.2092   1.2573   1.3092   1.3655   1.4270   1.4955   1.5726   1.6606   1.7655   1.8970   2.0726   2.3564   2.9504   % 50
                           -3.3676  -2.8431  -2.5968  -2.4432  -2.3307  -2.2404  -2.1650  -2.1000  -2.0425  -1.9913  -1.9446  -1.9014  -1.8608  -1.8238  -1.7890  -1.7557  -1.7244  -1.6945  -1.6664  -1.6393  -1.6132  -1.4961  -1.3956  -1.3068  -1.2260   0.4143   0.5171   0.6296   0.7555   0.9004   0.9323   0.9656   1.0001   1.0361   1.0743   1.1143   1.1568   1.2020   1.2501   1.3016   1.3576   1.4188   1.4862   1.5618   1.6495   1.7530   1.8808   2.0535   2.3300   2.9151   % 75
                           -3.3462  -2.8336  -2.5892  -2.4383  -2.3265  -2.2375  -2.1634  -2.0988  -2.0421  -1.9910  -1.9444  -1.9012  -1.8611  -1.8239  -1.7892  -1.7565  -1.7254  -1.6959  -1.6676  -1.6405  -1.6143  -1.4978  -1.3975  -1.3088  -1.2284   0.4114   0.5144   0.6268   0.7524   0.8970   0.9289   0.9621   0.9966   1.0326   1.0705   1.1106   1.1528   1.1974   1.2452   1.2962   1.3517   1.4121   1.4793   1.5554   1.6418   1.7454   1.8726   2.0433   2.3199   2.8880   % 100
                           -3.3259  -2.8200  -2.5810  -2.4314  -2.3210  -2.2333  -2.1593  -2.0960  -2.0396  -1.9888  -1.9425  -1.9001  -1.8606  -1.8239  -1.7890  -1.7565  -1.7256  -1.6962  -1.6682  -1.6414  -1.6153  -1.4990  -1.3992  -1.3104  -1.2304   0.4089   0.5113   0.6237   0.7492   0.8933   0.9255   0.9586   0.9930   1.0291   1.0671   1.1067   1.1486   1.1930   1.2407   1.2918   1.3470   1.4075   1.4745   1.5493   1.6354   1.7375   1.8630   2.0324   2.3030   2.8613   % 150
                           -3.3187  -2.8136  -2.5753  -2.4281  -2.3189  -2.2319  -2.1588  -2.0952  -2.0390  -1.9883  -1.9423  -1.9000  -1.8605  -1.8233  -1.7890  -1.7563  -1.7256  -1.6962  -1.6684  -1.6416  -1.6158  -1.4998  -1.3999  -1.3113  -1.2313   0.4073   0.5095   0.6218   0.7477   0.8919   0.9236   0.9569   0.9914   1.0273   1.0649   1.1047   1.1467   1.1910   1.2384   1.2892   1.3441   1.4042   1.4704   1.5455   1.6312   1.7329   1.8586   2.0276   2.2977   2.8536   % 200
                           -3.3051  -2.8137  -2.5759  -2.4289  -2.3199  -2.2322  -2.1583  -2.0949  -2.0382  -1.9875  -1.9417  -1.8995  -1.8607  -1.8242  -1.7897  -1.7573  -1.7265  -1.6970  -1.6690  -1.6422  -1.6165  -1.5005  -1.4009  -1.3125  -1.2326   0.4065   0.5091   0.6214   0.7467   0.8903   0.9224   0.9552   0.9895   1.0253   1.0627   1.1025   1.1443   1.1887   1.2363   1.2870   1.3418   1.4019   1.4681   1.5423   1.6276   1.7291   1.8553   2.0237   2.2924   2.8395   % 300
                           -3.2997  -2.8061  -2.5701  -2.4239  -2.3157  -2.2291  -2.1562  -2.0930  -2.0371  -1.9870  -1.9411  -1.8991  -1.8600  -1.8233  -1.7890  -1.7565  -1.7258  -1.6966  -1.6687  -1.6416  -1.6160  -1.5005  -1.4008  -1.3127  -1.2328   0.4054   0.5082   0.6202   0.7454   0.8897   0.9215   0.9546   0.9889   1.0247   1.0624   1.1019   1.1435   1.1878   1.2347   1.2855   1.3398   1.3999   1.4655   1.5402   1.6257   1.7273   1.8518   2.0192   2.2854   2.8322   % 500
                           -3.2900  -2.8030  -2.5694  -2.4235  -2.3156  -2.2298  -2.1571  -2.0938  -2.0381  -1.9872  -1.9416  -1.8994  -1.8603  -1.8238  -1.7894  -1.7568  -1.7262  -1.6970  -1.6691  -1.6425  -1.6167  -1.5015  -1.4018  -1.3138  -1.2337   0.4040   0.5069   0.6193   0.7441   0.8883   0.9202   0.9531   0.9875   1.0232   1.0609   1.1006   1.1426   1.1868   1.2346   1.2856   1.3405   1.4000   1.4665   1.5405   1.6259   1.7268   1.8514   2.0181   2.2840   2.8320   % 1000
                           -3.2864  -2.7999  -2.5662  -2.4229  -2.3155  -2.2289  -2.1562  -2.0936  -2.0376  -1.9873  -1.9416  -1.8994  -1.8604  -1.8241  -1.7898  -1.7575  -1.7268  -1.6977  -1.6697  -1.6430  -1.6175  -1.5022  -1.4025  -1.3144  -1.2343   0.4048   0.5074   0.6198   0.7448   0.8889   0.9208   0.9535   0.9878   1.0239   1.0614   1.1009   1.1424   1.1866   1.2336   1.2843   1.3388   1.3989   1.4646   1.5387   1.6239   1.7250   1.8489   2.0156   2.2816   2.8240]; % 10000
            case 'T2'

              % Row headers (effective sample sizes) are at the right of the table.
              % ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   ----- 
              % Alpha        0.001    0.005    0.010    0.015    0.020    0.025    0.030    0.035    0.040    0.045    0.050    0.055    0.060    0.065    0.070    0.075    0.080    0.085    0.090    0.095    0.100    0.125    0.150    0.175    0.200    0.800    0.825    0.850    0.875    0.900    0.905    0.910    0.915    0.920    0.925    0.930    0.935    0.940    0.945    0.950    0.955    0.960    0.965    0.970    0.975    0.980    0.985    0.990    0.995    0.999      T
              % ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
                CVTable = [-13.7043 -11.2023  -9.9515  -9.1596  -8.5743  -8.0998  -7.7040  -7.3643  -7.0628  -6.7942  -6.5512  -6.3285  -6.1240  -5.9350  -5.7580  -5.5912  -5.4350  -5.2865  -5.1462  -5.0140  -4.8874  -4.3259  -3.8636  -3.4691  -3.1252   0.6182   0.7355   0.8615   1.0017   1.1647   1.2013   1.2394   1.2796   1.3220   1.3672   1.4149   1.4660   1.5213   1.5807   1.6457   1.7183   1.8000   1.8921   2.0001   2.1294   2.2913   2.5052   2.8174   3.3865   4.8865   % 10
                           -15.6866 -12.5024 -10.9510  -9.9967  -9.3032  -8.7507  -8.2910  -7.9013  -7.5598  -7.2554  -6.9804  -6.7329  -6.5033  -6.2914  -6.0931  -5.9081  -5.7352  -5.5720  -5.4182  -5.2720  -5.1328  -4.5263  -4.0281  -3.6065  -3.2418   0.5626   0.6743   0.7933   0.9241   1.0745   1.1077   1.1425   1.1789   1.2172   1.2577   1.3005   1.3461   1.3950   1.4483   1.5059   1.5699   1.6412   1.7216   1.8152   1.9269   2.0641   2.2451   2.5042   2.9673   4.1421   % 15
                           -16.9638 -13.2686 -11.5387 -10.4844  -9.7183  -9.1210  -8.6266  -8.2051  -7.8386  -7.5121  -7.2172  -6.9507  -6.7053  -6.4804  -6.2727  -6.0774  -5.8956  -5.7240  -5.5614  -5.4086  -5.2632  -4.6299  -4.1123  -3.6753  -3.2992   0.5393   0.6482   0.7642   0.8907   1.0350   1.0669   1.1004   1.1354   1.1722   1.2109   1.2515   1.2951   1.3413   1.3910   1.4455   1.5056   1.5722   1.6482   1.7352   1.8394   1.9662   2.1304   2.3660   2.7835   3.8138   % 20
                           -17.8344 -13.7769 -11.9166 -10.7873  -9.9798  -9.3436  -8.8267  -8.3861  -8.0015  -7.6640  -7.3604  -7.0825  -6.8318  -6.5987  -6.3836  -6.1831  -5.9962  -5.8192  -5.6524  -5.4942  -5.3453  -4.6967  -4.1659  -3.7203  -3.3369   0.5243   0.6319   0.7461   0.8706   1.0120   1.0433   1.0759   1.1100   1.1456   1.1830   1.2227   1.2648   1.3100   1.3583   1.4113   1.4691   1.5338   1.6068   1.6913   1.7904   1.9119   2.0708   2.2958   2.6887   3.6556   % 25
                           -18.4608 -14.1486 -12.1881 -11.0179 -10.1840  -9.5299  -8.9906  -8.5353  -8.1398  -7.7917  -7.4775  -7.1922  -6.9319  -6.6939  -6.4705  -6.2641  -6.0706  -5.8889  -5.7190  -5.5579  -5.4056  -4.7438  -4.2073  -3.7534  -3.3641   0.5153   0.6221   0.7348   0.8581   0.9972   1.0281   1.0602   1.0933   1.1281   1.1650   1.2040   1.2453   1.2892   1.3367   1.3886   1.4452   1.5083   1.5795   1.6617   1.7589   1.8765   2.0285   2.2442   2.6248   3.5546   % 30
                           -19.3220 -14.6373 -12.5462 -11.3101 -10.4276  -9.7487  -9.1884  -8.7152  -8.2987  -7.9346  -7.6085  -7.3153  -7.0479  -6.8025  -6.5737  -6.3616  -6.1635  -5.9772  -5.8035  -5.6380  -5.4823  -4.8053  -4.2564  -3.7962  -3.4005   0.5024   0.6080   0.7197   0.8408   0.9780   1.0080   1.0393   1.0721   1.1065   1.1422   1.1804   1.2207   1.2636   1.3101   1.3604   1.4153   1.4764   1.5457   1.6249   1.7182   1.8323   1.9794   2.1895   2.5524   3.4318   % 40
                           -19.7788 -14.9237 -12.7582 -11.4840 -10.5730  -9.8708  -9.3014  -8.8142  -8.3967  -8.0246  -7.6940  -7.3949  -7.1212  -6.8688  -6.6368  -6.4208  -6.2197  -6.0315  -5.8532  -5.6866  -5.5274  -4.8395  -4.2849  -3.8192  -3.4189   0.4970   0.6018   0.7126   0.8331   0.9689   0.9988   1.0298   1.0621   1.0956   1.1311   1.1684   1.2083   1.2510   1.2965   1.3462   1.4003   1.4600   1.5281   1.6056   1.6962   1.8085   1.9533   2.1574   2.5108   3.3661   % 50
                           -20.5058 -15.3097 -13.0662 -11.7298 -10.7869 -10.0515  -9.4594  -8.9536  -8.5229  -8.1401  -7.8010  -7.4953  -7.2155  -6.9558  -6.7183  -6.4990  -6.2919  -6.0992  -5.9176  -5.7458  -5.5838  -4.8859  -4.3206  -3.8477  -3.4431   0.4884   0.5926   0.7025   0.8213   0.9550   0.9842   1.0144   1.0462   1.0794   1.1144   1.1511   1.1901   1.2319   1.2764   1.3248   1.3778   1.4361   1.5023   1.5780   1.6669   1.7747   1.9151   2.1132   2.4540   3.2737   % 75
                           -20.8574 -15.5513 -13.2236 -11.8579 -10.8964 -10.1506  -9.5474  -9.0381  -8.5940  -8.2067  -7.8643  -7.5503  -7.2675  -7.0061  -6.7663  -6.5419  -6.3346  -6.1391  -5.9551  -5.7831  -5.6191  -4.9118  -4.3422  -3.8661  -3.4585   0.4839   0.5878   0.6973   0.8154   0.9483   0.9777   1.0078   1.0392   1.0721   1.1068   1.1432   1.1819   1.2231   1.2673   1.3150   1.3672   1.4254   1.4903   1.5647   1.6524   1.7602   1.8979   2.0926   2.4265   3.2180   % 100
                           -21.2716 -15.7405 -13.3679 -11.9752 -10.9886 -10.2347  -9.6216  -9.1022  -8.6578  -8.2649  -7.9178  -7.6010  -7.3148  -7.0517  -6.8088  -6.5835  -6.3719  -6.1745  -5.9903  -5.8165  -5.6506  -4.9371  -4.3617  -3.8814  -3.4718   0.4792   0.5829   0.6918   0.8088   0.9405   0.9694   0.9995   1.0310   1.0635   1.0982   1.1345   1.1725   1.2133   1.2572   1.3049   1.3567   1.4140   1.4783   1.5520   1.6391   1.7449   1.8804   2.0736   2.4061   3.1864   % 150
                           -21.5166 -15.8389 -13.4301 -12.0334 -11.0503 -10.2860  -9.6697  -9.1485  -8.6978  -8.3015  -7.9487  -7.6305  -7.3420  -7.0757  -6.8290  -6.6020  -6.3895  -6.1919  -6.0075  -5.8322  -5.6663  -4.9506  -4.3715  -3.8898  -3.4791   0.4769   0.5802   0.6886   0.8059   0.9377   0.9664   0.9966   1.0281   1.0606   1.0947   1.1308   1.1689   1.2097   1.2531   1.3002   1.3520   1.4088   1.4726   1.5459   1.6322   1.7375   1.8717   2.0619   2.3872   3.1667   % 200
                           -21.6927 -16.0114 -13.5560 -12.1393 -11.1316 -10.3523  -9.7233  -9.1932  -8.7357  -8.3359  -7.9756  -7.6587  -7.3682  -7.1021  -6.8570  -6.6282  -6.4149  -6.2140  -6.0273  -5.8492  -5.6821  -4.9624  -4.3840  -3.9001  -3.4867   0.4753   0.5786   0.6867   0.8033   0.9345   0.9632   0.9928   1.0239   1.0562   1.0903   1.1262   1.1642   1.2046   1.2479   1.2948   1.3461   1.4030   1.4668   1.5393   1.6255   1.7297   1.8646   2.0530   2.3742   3.1301   % 300
                           -21.8988 -16.0669 -13.5975 -12.1626 -11.1600 -10.3850  -9.7497  -9.2170  -8.7593  -8.3565  -8.0015  -7.6793  -7.3852  -7.1170  -6.8701  -6.6412  -6.4259  -6.2253  -6.0373  -5.8609  -5.6939  -4.9690  -4.3879  -3.9021  -3.4909   0.4735   0.5765   0.6846   0.8013   0.9323   0.9609   0.9903   1.0211   1.0531   1.0869   1.1227   1.1605   1.2009   1.2441   1.2908   1.3418   1.3982   1.4620   1.5348   1.6199   1.7232   1.8559   2.0425   2.3618   3.1180   % 500
                           -21.9844 -16.1445 -13.6637 -12.2158 -11.1999 -10.4207  -9.7848  -9.2469  -8.7868  -8.3843  -8.0244  -7.7004  -7.4055  -7.1357  -6.8875  -6.6576  -6.4426  -6.2409  -6.0511  -5.8732  -5.7040  -4.9813  -4.3981  -3.9106  -3.4967   0.4717   0.5750   0.6830   0.7997   0.9301   0.9587   0.9882   1.0193   1.0516   1.0855   1.1210   1.1587   1.1987   1.2419   1.2885   1.3398   1.3962   1.4598   1.5323   1.6175   1.7208   1.8529   2.0402   2.3594   3.1100   % 1000
                           -22.1154 -16.1974 -13.6952 -12.2627 -11.2401 -10.4532  -9.8133  -9.2735  -8.8112  -8.4057  -8.0443  -7.7189  -7.4208  -7.1504  -6.9002  -6.6687  -6.4540  -6.2529  -6.0643  -5.8857  -5.7170  -4.9897  -4.4043  -3.9158  -3.5011   0.4717   0.5745   0.6824   0.7988   0.9288   0.9572   0.9868   1.0177   1.0498   1.0836   1.1191   1.1566   1.1967   1.2401   1.2867   1.3374   1.3938   1.4564   1.5287   1.6138   1.7169   1.8491   2.0341   2.3518   3.1011]; % 10000
  
            case 'F'
                
                error(message('econ:adftest:FTestForARModel'))
                
        end
        
    case 'ARD'
        
        switch upper(testType)
            
            case 'T1'

              % Row headers (effective sample sizes) are at the right of the table.
              % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
              % Alpha       0.001    0.005    0.010    0.015    0.020    0.025    0.030    0.035    0.040    0.045    0.050    0.055    0.060    0.065    0.070    0.075    0.080    0.085    0.090    0.095    0.100    0.125    0.150    0.175    0.200    0.800    0.825    0.850    0.875    0.900    0.905    0.910    0.915    0.920    0.925    0.930    0.935    0.940    0.945    0.950    0.955    0.960    0.965    0.970    0.975    0.980    0.985    0.990    0.995    0.999      T
              % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
                CVTable = [-6.1728  -4.8563  -4.3427  -4.0541  -3.8538  -3.6994  -3.5751  -3.4691  -3.3782  -3.2976  -3.2260  -3.1609  -3.1022  -3.0484  -2.9978  -2.9510  -2.9071  -2.8657  -2.8268  -2.7900  -2.7552  -2.6001  -2.4713  -2.3608  -2.2633  -0.7119  -0.6179  -0.5134  -0.3960  -0.2595  -0.2292  -0.1974  -0.1644  -0.1299  -0.0930  -0.0546  -0.0139   0.0303   0.0768   0.1279   0.1836   0.2433   0.3110   0.3881   0.4772   0.5851   0.7205   0.9090   1.2323   1.9914   % 10
                           -5.2336  -4.3409  -3.9642  -3.7444  -3.5882  -3.4673  -3.3673  -3.2817  -3.2073  -3.1422  -3.0831  -3.0287  -2.9793  -2.9341  -2.8917  -2.8519  -2.8148  -2.7791  -2.7453  -2.7132  -2.6830  -2.5474  -2.4334  -2.3345  -2.2459  -0.7656  -0.6733  -0.5711  -0.4560  -0.3226  -0.2929  -0.2618  -0.2295  -0.1954  -0.1593  -0.1218  -0.0820  -0.0395   0.0057   0.0547   0.1079   0.1666   0.2323   0.3054   0.3906   0.4916   0.6203   0.7972   1.0849   1.7257   % 15
                           -4.8765  -4.1350  -3.8106  -3.6175  -3.4789  -3.3702  -3.2803  -3.2034  -3.1356  -3.0750  -3.0215  -2.9724  -2.9270  -2.8848  -2.8452  -2.8083  -2.7736  -2.7405  -2.7093  -2.6793  -2.6507  -2.5243  -2.4173  -2.3231  -2.2388  -0.7911  -0.6996  -0.5989  -0.4851  -0.3515  -0.3221  -0.2915  -0.2594  -0.2259  -0.1907  -0.1531  -0.1139  -0.0723  -0.0276   0.0211   0.0730   0.1306   0.1943   0.2669   0.3512   0.4514   0.5753   0.7454   1.0240   1.6287   % 20
                           -4.6989  -4.0165  -3.7208  -3.5416  -3.4119  -3.3118  -3.2278  -3.1547  -3.0912  -3.0349  -2.9846  -2.9376  -2.8946  -2.8546  -2.8171  -2.7819  -2.7482  -2.7167  -2.6865  -2.6577  -2.6307  -2.5087  -2.4056  -2.3144  -2.2327  -0.8056  -0.7154  -0.6153  -0.5026  -0.3705  -0.3411  -0.3107  -0.2787  -0.2454  -0.2101  -0.1731  -0.1338  -0.0922  -0.0477   0.0003   0.0523   0.1094   0.1730   0.2453   0.3274   0.4254   0.5484   0.7128   0.9842   1.5729   % 25
                           -4.5788  -3.9518  -3.6686  -3.4995  -3.3770  -3.2792  -3.1977  -3.1285  -3.0677  -3.0129  -2.9634  -2.9182  -2.8765  -2.8379  -2.8015  -2.7675  -2.7354  -2.7048  -2.6756  -2.6477  -2.6209  -2.5020  -2.4004  -2.3110  -2.2310  -0.8151  -0.7256  -0.6264  -0.5132  -0.3814  -0.3519  -0.3219  -0.2901  -0.2571  -0.2222  -0.1855  -0.1469  -0.1056  -0.0613  -0.0135   0.0381   0.0949   0.1579   0.2288   0.3112   0.4100   0.5302   0.6966   0.9632   1.5302   % 30
                           -4.4417  -3.8676  -3.6060  -3.4472  -3.3307  -3.2382  -3.1603  -3.0944  -3.0365  -2.9845  -2.9376  -2.8941  -2.8544  -2.8163  -2.7812  -2.7484  -2.7171  -2.6878  -2.6594  -2.6327  -2.6069  -2.4910  -2.3929  -2.3065  -2.2278  -0.8288  -0.7399  -0.6411  -0.5286  -0.3981  -0.3693  -0.3391  -0.3077  -0.2745  -0.2395  -0.2032  -0.1646  -0.1232  -0.0789  -0.0317   0.0195   0.0757   0.1386   0.2094   0.2900   0.3865   0.5064   0.6694   0.9330   1.4861   % 40
                           -4.3679  -3.8235  -3.5693  -3.4173  -3.3048  -3.2155  -3.1405  -3.0758  -3.0194  -2.9683  -2.9219  -2.8790  -2.8400  -2.8036  -2.7694  -2.7369  -2.7064  -2.6776  -2.6503  -2.6243  -2.5990  -2.4857  -2.3887  -2.3031  -2.2260  -0.8364  -0.7476  -0.6488  -0.5372  -0.4063  -0.3773  -0.3473  -0.3156  -0.2825  -0.2477  -0.2109  -0.1722  -0.1307  -0.0865  -0.0390   0.0119   0.0685   0.1309   0.2019   0.2832   0.3788   0.4972   0.6595   0.9215   1.4687   % 50
                           -4.2748  -3.7604  -3.5208  -3.3747  -3.2665  -3.1812  -3.1099  -3.0479  -2.9931  -2.9443  -2.9002  -2.8592  -2.8215  -2.7860  -2.7528  -2.7217  -2.6925  -2.6644  -2.6377  -2.6120  -2.5875  -2.4774  -2.3830  -2.2992  -2.2236  -0.8450  -0.7573  -0.6588  -0.5477  -0.4181  -0.3893  -0.3593  -0.3283  -0.2951  -0.2602  -0.2240  -0.1854  -0.1451  -0.1011  -0.0544  -0.0032   0.0532   0.1146   0.1844   0.2649   0.3600   0.4805   0.6393   0.8977   1.4358   % 75
                           -4.2214  -3.7298  -3.4980  -3.3565  -3.2505  -3.1667  -3.0962  -3.0355  -2.9816  -2.9334  -2.8897  -2.8497  -2.8123  -2.7776  -2.7450  -2.7141  -2.6846  -2.6571  -2.6305  -2.6054  -2.5811  -2.4721  -2.3782  -2.2956  -2.2205  -0.8493  -0.7611  -0.6633  -0.5517  -0.4222  -0.3936  -0.3634  -0.3319  -0.2997  -0.2652  -0.2287  -0.1900  -0.1487  -0.1047  -0.0575  -0.0063   0.0498   0.1111   0.1809   0.2608   0.3551   0.4735   0.6323   0.8853   1.4156   % 100
                           -4.1862  -3.7032  -3.4774  -3.3387  -3.2348  -3.1526  -3.0834  -3.0241  -2.9716  -2.9242  -2.8812  -2.8416  -2.8049  -2.7705  -2.7382  -2.7077  -2.6791  -2.6519  -2.6256  -2.6008  -2.5769  -2.4694  -2.3770  -2.2952  -2.2209  -0.8553  -0.7679  -0.6700  -0.5592  -0.4306  -0.4022  -0.3720  -0.3407  -0.3078  -0.2734  -0.2372  -0.1986  -0.1576  -0.1141  -0.0671  -0.0164   0.0392   0.1012   0.1707   0.2508   0.3451   0.4634   0.6217   0.8745   1.4054   % 150
                           -4.1587  -3.6859  -3.4616  -3.3261  -3.2248  -3.1442  -3.0760  -3.0177  -2.9651  -2.9186  -2.8760  -2.8371  -2.8006  -2.7668  -2.7351  -2.7047  -2.6761  -2.6489  -2.6229  -2.5982  -2.5744  -2.4671  -2.3750  -2.2935  -2.2194  -0.8562  -0.7686  -0.6711  -0.5604  -0.4311  -0.4026  -0.3727  -0.3417  -0.3091  -0.2748  -0.2384  -0.2003  -0.1599  -0.1163  -0.0698  -0.0191   0.0361   0.0972   0.1669   0.2475   0.3430   0.4603   0.6176   0.8710   1.3933   % 200
                           -4.1353  -3.6751  -3.4561  -3.3198  -3.2183  -3.1371  -3.0700  -3.0116  -2.9595  -2.9135  -2.8708  -2.8322  -2.7961  -2.7625  -2.7308  -2.7010  -2.6729  -2.6461  -2.6207  -2.5959  -2.5722  -2.4658  -2.3743  -2.2929  -2.2193  -0.8597  -0.7726  -0.6748  -0.5639  -0.4350  -0.4061  -0.3758  -0.3446  -0.3122  -0.2776  -0.2413  -0.2030  -0.1622  -0.1188  -0.0721  -0.0210   0.0351   0.0967   0.1661   0.2454   0.3405   0.4578   0.6164   0.8683   1.3894   % 300
                           -4.1168  -3.6611  -3.4423  -3.3084  -3.2097  -3.1314  -3.0651  -3.0073  -2.9562  -2.9099  -2.8680  -2.8293  -2.7931  -2.7599  -2.7286  -2.6988  -2.6708  -2.6440  -2.6183  -2.5937  -2.5702  -2.4647  -2.3736  -2.2928  -2.2195  -0.8623  -0.7747  -0.6769  -0.5663  -0.4371  -0.4080  -0.3779  -0.3468  -0.3141  -0.2797  -0.2434  -0.2050  -0.1647  -0.1213  -0.0744  -0.0236   0.0319   0.0930   0.1624   0.2417   0.3363   0.4534   0.6111   0.8627   1.3812   % 500
                           -4.0982  -3.6539  -3.4378  -3.3035  -3.2051  -3.1266  -3.0606  -3.0035  -2.9529  -2.9066  -2.8648  -2.8258  -2.7902  -2.7569  -2.7259  -2.6963  -2.6686  -2.6420  -2.6166  -2.5923  -2.5686  -2.4625  -2.3718  -2.2912  -2.2182  -0.8630  -0.7754  -0.6780  -0.5669  -0.4380  -0.4090  -0.3793  -0.3482  -0.3158  -0.2813  -0.2452  -0.2068  -0.1662  -0.1226  -0.0757  -0.0250   0.0302   0.0906   0.1606   0.2396   0.3339   0.4504   0.6087   0.8628   1.3790   % 1000
                           -4.0932  -3.6446  -3.4307  -3.2987  -3.2005  -3.1217  -3.0563  -2.9990  -2.9483  -2.9028  -2.8610  -2.8227  -2.7875  -2.7544  -2.7231  -2.6940  -2.6663  -2.6395  -2.6142  -2.5899  -2.5666  -2.4615  -2.3713  -2.2904  -2.2174  -0.8630  -0.7756  -0.6782  -0.5669  -0.4378  -0.4094  -0.3796  -0.3481  -0.3153  -0.2814  -0.2455  -0.2073  -0.1672  -0.1238  -0.0778  -0.0275   0.0276   0.0898   0.1593   0.2395   0.3342   0.4515   0.6093   0.8588   1.3761]; % 10000
   
            case 'T2'
                
              % Row headers (effective sample sizes) are at the right of the table.
              % ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
              % Alpha        0.001    0.005    0.010    0.015    0.020    0.025    0.030    0.035    0.040    0.045    0.050    0.055    0.060    0.065    0.070    0.075    0.080    0.085    0.090    0.095    0.100    0.125    0.150    0.175    0.200    0.800    0.825    0.850    0.875    0.900    0.905    0.910    0.915    0.920    0.925    0.930    0.935    0.940    0.945    0.950    0.955    0.960    0.965    0.970    0.975    0.980    0.985    0.990    0.995    0.999      T
              % ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
                CVTable = [-15.8375 -13.7705 -12.7397 -12.0828 -11.5980 -11.1995 -10.8662 -10.5720 -10.3111 -10.0777  -9.8618  -9.6626  -9.4784  -9.3066  -9.1462  -8.9946  -8.8518  -8.7142  -8.5843  -8.4584  -8.3390  -7.8052  -7.3499  -6.9495  -6.5887  -1.5192  -1.3019  -1.0689  -0.8156  -0.5283  -0.4651  -0.3996  -0.3319  -0.2617  -0.1880  -0.1098  -0.0275   0.0607   0.1548   0.2564   0.3674   0.4891   0.6241   0.7778   0.9582   1.1790   1.4607   1.8588   2.5438   4.2812   % 10
                           -18.8053 -16.0373 -14.6542 -13.8012 -13.1676 -12.6568 -12.2275 -11.8592 -11.5339 -11.2439 -10.9800 -10.7371 -10.5117 -10.3003 -10.1019  -9.9204  -9.7488  -9.5856  -9.4291  -9.2799  -9.1367  -8.5029  -7.9693  -7.5036  -7.0924  -1.6248  -1.4069  -1.1747  -0.9231  -0.6434  -0.5817  -0.5188  -0.4533  -0.3847  -0.3132  -0.2386  -0.1603  -0.0769   0.0110   0.1058   0.2091   0.3223   0.4475   0.5884   0.7517   0.9453   1.1939   1.5380   2.1208   3.5134   % 15
                           -20.8206 -17.4581 -15.8350 -14.8357 -14.1025 -13.5192 -13.0399 -12.6239 -12.2604 -11.9329 -11.6377 -11.3649 -11.1139 -10.8823 -10.6637 -10.4612 -10.2706 -10.0891  -9.9162  -9.7523  -9.5952  -8.9024  -8.3198  -7.8178  -7.3742  -1.6760  -1.4569  -1.2242  -0.9732  -0.6935  -0.6328  -0.5708  -0.5054  -0.4388  -0.3687  -0.2961  -0.2195  -0.1386  -0.0523   0.0401   0.1384   0.2470   0.3677   0.5046   0.6597   0.8456   1.0816   1.4023   1.9352   3.1866   % 20
                           -22.2099 -18.3844 -16.6279 -15.5255 -14.7259 -14.0928 -13.5676 -13.1137 -12.7212 -12.3690 -12.0515 -11.7582 -11.4884 -11.2426 -11.0116 -10.7946 -10.5929 -10.4014 -10.2174 -10.0442  -9.8781  -9.1463  -8.5356  -8.0101  -7.5473  -1.7037  -1.4848  -1.2538  -1.0045  -0.7265  -0.6661  -0.6038  -0.5398  -0.4734  -0.4039  -0.3314  -0.2556  -0.1750  -0.0904   0.0007   0.0981   0.2055   0.3228   0.4565   0.6092   0.7906   1.0160   1.3235   1.8344   2.9974   % 25
                           -23.2847 -19.1040 -17.1857 -16.0273 -15.1852 -14.5147 -13.9589 -13.4857 -13.0641 -12.6927 -12.3570 -12.0519 -11.7731 -11.5147 -11.2735 -11.0465 -10.8347 -10.6328 -10.4430 -10.2632 -10.0906  -9.3310  -8.7032  -8.1608  -7.6822  -1.7247  -1.5056  -1.2724  -1.0221  -0.7451  -0.6851  -0.6232  -0.5596  -0.4937  -0.4252  -0.3532  -0.2791  -0.2001  -0.1157  -0.0252   0.0714   0.1760   0.2929   0.4239   0.5727   0.7500   0.9711   1.2773   1.7756   2.8900   % 30
                           -24.6748 -20.0868 -17.9665 -16.6939 -15.7773 -15.0576 -14.4607 -13.9571 -13.5137 -13.1165 -12.7609 -12.4371 -12.1411 -11.8637 -11.6073 -11.3680 -11.1430 -10.9300 -10.7313 -10.5399 -10.3601  -9.5657  -8.9029  -8.3360  -7.8401  -1.7522  -1.5328  -1.3006  -1.0506  -0.7749  -0.7157  -0.6540  -0.5915  -0.5255  -0.4569  -0.3855  -0.3105  -0.2307  -0.1480  -0.0592   0.0359   0.1395   0.2545   0.3835   0.5289   0.7009   0.9161   1.2099   1.6878   2.7574   % 40
                           -25.6040 -20.6964 -18.4622 -17.1337 -16.1678 -15.4068 -14.7833 -14.2531 -13.7864 -13.3725 -13.0052 -12.6709 -12.3638 -12.0791 -11.8147 -11.5671 -11.3359 -11.1193 -10.9142 -10.7177 -10.5309  -9.7097  -9.0296  -8.4479  -7.9380  -1.7675  -1.5475  -1.3158  -1.0652  -0.7881  -0.7284  -0.6672  -0.6029  -0.5374  -0.4692  -0.3982  -0.3240  -0.2447  -0.1604  -0.0717   0.0221   0.1254   0.2398   0.3668   0.5109   0.6824   0.8959   1.1834   1.6475   2.6896   % 50
                           -27.0188 -21.5221 -19.1134 -17.6888 -16.6679 -15.8719 -15.2117 -14.6490 -14.1565 -13.7246 -13.3386 -12.9885 -12.6668 -12.3707 -12.0958 -11.8382 -11.5977 -11.3698 -11.1571 -10.9529 -10.7595  -9.9049  -9.2039  -8.6024  -8.0767  -1.7857  -1.5655  -1.3322  -1.0831  -0.8096  -0.7497  -0.6880  -0.6246  -0.5594  -0.4916  -0.4211  -0.3471  -0.2695  -0.1875  -0.0999  -0.0058   0.0968   0.2071   0.3321   0.4763   0.6455   0.8523   1.1376   1.5959   2.5937   % 75
                           -27.6550 -22.0112 -19.4956 -18.0122 -16.9461 -16.1210 -15.4505 -14.8675 -14.3645 -13.9188 -13.5204 -13.1588 -12.8290 -12.5270 -12.2465 -11.9797 -11.7315 -11.4954 -11.2768 -11.0673 -10.8678 -10.0032  -9.2866  -8.6779  -8.1445  -1.7926  -1.5715  -1.3392  -1.0894  -0.8147  -0.7557  -0.6946  -0.6320  -0.5665  -0.4983  -0.4279  -0.3540  -0.2762  -0.1931  -0.1054  -0.0113   0.0902   0.2007   0.3253   0.4668   0.6334   0.8387   1.1193   1.5693   2.5426   % 100
                           -28.4397 -22.4737 -19.8593 -18.3296 -17.2385 -16.3874 -15.6923 -15.0991 -14.5785 -14.1241 -13.7150 -13.3460 -13.0047 -12.6958 -12.4052 -12.1347 -11.8809 -11.6434 -11.4179 -11.2052 -11.0028 -10.1195  -9.3893  -8.7662  -8.2247  -1.8075  -1.5868  -1.3540  -1.1042  -0.8305  -0.7711  -0.7106  -0.6471  -0.5816  -0.5141  -0.4432  -0.3698  -0.2923  -0.2105  -0.1228  -0.0300   0.0708   0.1831   0.3064   0.4475   0.6151   0.8180   1.0931   1.5384   2.5077   % 150
                           -28.8138 -22.6913 -20.0186 -18.4666 -17.3621 -16.5055 -15.7911 -15.1899 -14.6710 -14.2126 -13.8005 -13.4226 -13.0797 -12.7654 -12.4747 -12.2033 -11.9461 -11.7063 -11.4758 -11.2573 -11.0520 -10.1595  -9.4292  -8.8026  -8.2539  -1.8095  -1.5888  -1.3557  -1.1056  -0.8299  -0.7710  -0.7098  -0.6469  -0.5821  -0.5149  -0.4450  -0.3720  -0.2963  -0.2144  -0.1279  -0.0346   0.0652   0.1754   0.2990   0.4396   0.6061   0.8086   1.0814   1.5247   2.4841   % 200
                           -29.1931 -22.9641 -20.2618 -18.6653 -17.5354 -16.6522 -15.9297 -15.3172 -14.7866 -14.3169 -13.9009 -13.5212 -13.1763 -12.8533 -12.5564 -12.2766 -12.0181 -11.7763 -11.5468 -11.3274 -11.1231 -10.2194  -9.4778  -8.8443  -8.2945  -1.8171  -1.5954  -1.3617  -1.1121  -0.8360  -0.7764  -0.7159  -0.6526  -0.5876  -0.5203  -0.4503  -0.3767  -0.2991  -0.2180  -0.1315  -0.0382   0.0635   0.1741   0.2962   0.4356   0.6017   0.8068   1.0804   1.5171   2.4536   % 300
                           -29.4773 -23.1410 -20.4013 -18.8081 -17.6688 -16.7770 -16.0449 -15.4259 -14.8850 -14.4128 -13.9926 -13.6025 -13.2532 -12.9305 -12.6320 -12.3549 -12.0932 -11.8460 -11.6127 -11.3921 -11.1842 -10.2743  -9.5244  -8.8835  -8.3260  -1.8212  -1.5997  -1.3657  -1.1151  -0.8392  -0.7800  -0.7191  -0.6562  -0.5907  -0.5235  -0.4541  -0.3808  -0.3031  -0.2229  -0.1355  -0.0431   0.0573   0.1669   0.2900   0.4282   0.5925   0.7956   1.0671   1.5009   2.4355   % 500
                           -29.7109 -23.3043 -20.5073 -18.8833 -17.7392 -16.8359 -16.0975 -15.4767 -14.9354 -14.4605 -14.0384 -13.6547 -13.3049 -12.9775 -12.6716 -12.3915 -12.1274 -11.8790 -11.6460 -11.4241 -11.2148 -10.2969  -9.5453  -8.9024  -8.3420  -1.8221  -1.6002  -1.3672  -1.1167  -0.8409  -0.7816  -0.7206  -0.6579  -0.5932  -0.5263  -0.4564  -0.3832  -0.3057  -0.2239  -0.1379  -0.0452   0.0541   0.1629   0.2853   0.4235   0.5866   0.7876   1.0589   1.4964   2.4235   % 1000
                           -29.9070 -23.4318 -20.6163 -18.9767 -17.8146 -16.9008 -16.1627 -15.5440 -14.9996 -14.5175 -14.0903 -13.7007 -13.3462 -13.0192 -12.7130 -12.4276 -12.1641 -11.9147 -11.6807 -11.4598 -11.2464 -10.3215  -9.5627  -8.9202  -8.3587  -1.8220  -1.6005  -1.3660  -1.1160  -0.8401  -0.7808  -0.7203  -0.6574  -0.5930  -0.5259  -0.4558  -0.3833  -0.3073  -0.2264  -0.1410  -0.0495   0.0497   0.1607   0.2838   0.4227   0.5869   0.7898   1.0579   1.4882   2.4173]; % 10000
    
            case 'F'
                
              % Row headers (effective sample sizes) are at the right of the table.
              % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
              % Alpha      0.001    0.005    0.010    0.015    0.020    0.025    0.030    0.035    0.040    0.045    0.050    0.055    0.060    0.065    0.070    0.075    0.080    0.085    0.090    0.095    0.100    0.125    0.150    0.175    0.200    0.800    0.825    0.850    0.875    0.900    0.905    0.910    0.915    0.920    0.925    0.930    0.935    0.940    0.945    0.950    0.955    0.960    0.965    0.970    0.975    0.980    0.985    0.990    0.995    0.999      T
              % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
               CVTable = [23.5197  14.5283  11.6684  10.1602   9.1959   8.4921   7.9435   7.4979   7.1210   6.8015   6.5231   6.2746   6.0526   5.8534   5.6690   5.5019   5.3486   5.2068   5.0743   4.9504   4.8357   4.3474   3.9645   3.6502   3.3862   0.8663   0.8027   0.7385   0.6730   0.6045   0.5902   0.5757   0.5610   0.5462   0.5309   0.5154   0.4994   0.4830   0.4663   0.4492   0.4314   0.4129   0.3933   0.3728   0.3511   0.3272   0.3004   0.2696   0.2300   0.1751   % 10
                          15.8122  10.9448   9.1854   8.2326   7.5829   7.1015   6.7186   6.3999   6.1309   5.8965   5.6890   5.5044   5.3379   5.1861   5.0501   4.9226   4.8032   4.6920   4.5885   4.4918   4.4005   4.0082   3.6938   3.4327   3.2091   0.8923   0.8285   0.7636   0.6968   0.6273   0.6129   0.5982   0.5832   0.5681   0.5523   0.5363   0.5201   0.5035   0.4864   0.4684   0.4500   0.4309   0.4106   0.3890   0.3658   0.3405   0.3122   0.2790   0.2350   0.1740   % 15
                          13.3913   9.7048   8.2968   7.5123   6.9762   6.5703   6.2443   5.9752   5.7412   5.5399   5.3593   5.1984   5.0533   4.9211   4.7986   4.6854   4.5801   4.4817   4.3895   4.3023   4.2202   3.8655   3.5782   3.3396   3.1338   0.9059   0.8420   0.7768   0.7098   0.6396   0.6248   0.6098   0.5948   0.5795   0.5637   0.5476   0.5310   0.5139   0.4963   0.4780   0.4591   0.4394   0.4187   0.3968   0.3733   0.3477   0.3183   0.2840   0.2385   0.1731   % 20
                          12.2474   9.0658   7.8140   7.1158   6.6378   6.2709   5.9789   5.7315   5.5204   5.3346   5.1727   5.0248   4.8910   4.7676   4.6543   4.5489   4.4516   4.3598   4.2738   4.1926   4.1159   3.7813   3.5109   3.2837   3.0876   0.9132   0.8491   0.7841   0.7165   0.6459   0.6310   0.6160   0.6008   0.5852   0.5694   0.5531   0.5366   0.5196   0.5019   0.4833   0.4645   0.4445   0.4236   0.4012   0.3775   0.3511   0.3214   0.2860   0.2399   0.1729   % 25
                          11.5062   8.6789   7.5322   6.8938   6.4519   6.1058   5.8296   5.5975   5.3951   5.2198   5.0650   4.9236   4.7943   4.6767   4.5690   4.4699   4.3769   4.2894   4.2067   4.1286   4.0546   3.7357   3.4747   3.2556   3.0632   0.9187   0.8547   0.7890   0.7214   0.6502   0.6355   0.6205   0.6054   0.5897   0.5738   0.5577   0.5410   0.5238   0.5058   0.4874   0.4684   0.4483   0.4273   0.4047   0.3803   0.3537   0.3236   0.2878   0.2407   0.1717   % 30
                          10.7257   8.2446   7.2104   6.6276   6.2207   5.9047   5.6449   5.4311   5.2437   5.0778   4.9321   4.8005   4.6806   4.5713   4.4685   4.3733   4.2847   4.2020   4.1239   4.0500   3.9792   3.6751   3.4254   3.2133   3.0302   0.9246   0.8607   0.7952   0.7273   0.6558   0.6409   0.6257   0.6102   0.5946   0.5787   0.5622   0.5454   0.5281   0.5102   0.4919   0.4725   0.4525   0.4311   0.4085   0.3838   0.3568   0.3263   0.2901   0.2421   0.1712   % 40
                          10.3059   8.0110   7.0424   6.4834   6.0947   5.7937   5.5468   5.3401   5.1601   5.0028   4.8626   4.7358   4.6183   4.5118   4.4125   4.3202   4.2340   4.1536   4.0779   4.0062   3.9386   3.6421   3.3988   3.1922   3.0126   0.9288   0.8648   0.7991   0.7311   0.6595   0.6445   0.6294   0.6140   0.5981   0.5821   0.5655   0.5486   0.5313   0.5134   0.4948   0.4754   0.4551   0.4336   0.4108   0.3861   0.3591   0.3286   0.2920   0.2431   0.1727   % 50
                           9.8372   7.6930   6.8006   6.2836   5.9214   5.6380   5.4072   5.2130   5.0438   4.8944   4.7609   4.6402   4.5307   4.4289   4.3358   4.2473   4.1660   4.0889   4.0170   3.9479   3.8821   3.5965   3.3624   3.1625   2.9885   0.9334   0.8693   0.8036   0.7355   0.6635   0.6485   0.6332   0.6177   0.6019   0.5856   0.5692   0.5523   0.5348   0.5168   0.4981   0.4787   0.4581   0.4365   0.4133   0.3887   0.3613   0.3299   0.2931   0.2436   0.1712   % 75
                           9.5491   7.5378   6.6916   6.1938   5.8436   5.5674   5.3432   5.1535   4.9904   4.8450   4.7138   4.5946   4.4872   4.3880   4.2960   4.2110   4.1316   4.0563   3.9852   3.9175   3.8536   3.5729   3.3422   3.1463   2.9736   0.9365   0.8721   0.8066   0.7384   0.6664   0.6512   0.6359   0.6202   0.6044   0.5882   0.5713   0.5544   0.5371   0.5189   0.5000   0.4805   0.4598   0.4380   0.4149   0.3901   0.3631   0.3319   0.2943   0.2448   0.1726   % 100
                           9.3412   7.4086   6.5894   6.1078   5.7693   5.5002   5.2814   5.0955   4.9368   4.7953   4.6684   4.5546   4.4500   4.3530   4.2636   4.1800   4.1011   4.0272   3.9583   3.8919   3.8288   3.5546   3.3267   3.1328   2.9642   0.9388   0.8747   0.8088   0.7407   0.6687   0.6536   0.6381   0.6225   0.6067   0.5905   0.5737   0.5567   0.5390   0.5209   0.5021   0.4825   0.4619   0.4400   0.4167   0.3914   0.3640   0.3325   0.2952   0.2455   0.1720   % 150
                           9.2227   7.3388   6.5317   6.0617   5.7268   5.4650   5.2525   5.0704   4.9115   4.7729   4.6476   4.5351   4.4312   4.3359   4.2473   4.1645   4.0869   4.0144   3.9451   3.8792   3.8163   3.5421   3.3176   3.1250   2.9573   0.9388   0.8749   0.8094   0.7410   0.6690   0.6538   0.6386   0.6230   0.6071   0.5908   0.5742   0.5570   0.5394   0.5211   0.5023   0.4825   0.4620   0.4403   0.4170   0.3920   0.3642   0.3329   0.2951   0.2451   0.1719   % 200
                           9.0996   7.2794   6.4896   6.0261   5.6933   5.4359   5.2227   5.0407   4.8857   4.7482   4.6247   4.5121   4.4093   4.3156   4.2287   4.1468   4.0698   3.9980   3.9290   3.8644   3.8025   3.5332   3.3100   3.1188   2.9520   0.9414   0.8771   0.8112   0.7427   0.6706   0.6556   0.6403   0.6247   0.6087   0.5924   0.5757   0.5585   0.5409   0.5229   0.5040   0.4842   0.4634   0.4414   0.4181   0.3928   0.3651   0.3333   0.2958   0.2453   0.1712   % 300
                           9.0079   7.2151   6.4347   5.9784   5.6536   5.4002   5.1963   5.0203   4.8669   4.7310   4.6090   4.4980   4.3966   4.3024   4.2156   4.1341   4.0583   3.9854   3.9186   3.8544   3.7930   3.5246   3.3038   3.1148   2.9482   0.9435   0.8792   0.8131   0.7445   0.6723   0.6571   0.6418   0.6261   0.6101   0.5939   0.5772   0.5602   0.5427   0.5243   0.5053   0.4855   0.4647   0.4427   0.4192   0.3938   0.3659   0.3341   0.2963   0.2454   0.1716   % 500
                           8.9331   7.1857   6.4194   5.9624   5.6368   5.3849   5.1795   5.0050   4.8531   4.7178   4.5970   4.4873   4.3868   4.2932   4.2060   4.1255   4.0490   3.9782   3.9101   3.8464   3.7860   3.5196   3.2985   3.1091   2.9441   0.9428   0.8787   0.8131   0.7447   0.6725   0.6574   0.6420   0.6265   0.6105   0.5942   0.5772   0.5601   0.5423   0.5240   0.5052   0.4855   0.4646   0.4426   0.4192   0.3942   0.3664   0.3346   0.2966   0.2458   0.1714   % 1000
                           8.8985   7.1439   6.3854   5.9387   5.6169   5.3682   5.1632   4.9871   4.8369   4.7035   4.5818   4.4734   4.3733   4.2809   4.1954   4.1152   4.0404   3.9692   3.9025   3.8393   3.7784   3.5142   3.2951   3.1069   2.9427   0.9439   0.8797   0.8139   0.7454   0.6730   0.6578   0.6423   0.6267   0.6107   0.5943   0.5774   0.5602   0.5425   0.5244   0.5054   0.4855   0.4650   0.4430   0.4197   0.3945   0.3666   0.3349   0.2970   0.2461   0.1711]; % 10000
                       
        end
        
    case 'TS'
        
        switch upper(testType)
            
            case 'T1'
               
              % Row headers (effective sample sizes) are at the right of the table.
              % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
              % Alpha       0.001    0.005    0.010    0.015    0.020    0.025    0.030    0.035    0.040    0.045    0.050    0.055    0.060    0.065    0.070    0.075    0.080    0.085    0.090    0.095    0.100    0.125    0.150    0.175    0.200    0.800    0.825    0.850    0.875    0.900    0.905    0.910    0.915    0.920    0.925    0.930    0.935    0.940    0.945    0.950    0.955    0.960    0.965    0.970    0.975    0.980    0.985    0.990    0.995    0.999      T
              % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
                CVTable = [-7.6228  -5.9606  -5.3370  -4.9888  -4.7452  -4.5600  -4.4087  -4.2833  -4.1754  -4.0816  -3.9972  -3.9219  -3.8528  -3.7898  -3.7315  -3.6773  -3.6269  -3.5787  -3.5336  -3.4911  -3.4510  -3.2760  -3.1308  -3.0072  -2.8991  -1.3728  -1.2895  -1.1951  -1.0881  -0.9631  -0.9351  -0.9059  -0.8755  -0.8437  -0.8101  -0.7745  -0.7368  -0.6968  -0.6537  -0.6074  -0.5566  -0.5014  -0.4394  -0.3690  -0.2868  -0.1879  -0.0648   0.1060   0.3967   1.0702   % 10
                           -6.1981  -5.1511  -4.7279  -4.4842  -4.3123  -4.1772  -4.0672  -3.9743  -3.8934  -3.8218  -3.7567  -3.6979  -3.6444  -3.5944  -3.5487  -3.5057  -3.4652  -3.4272  -3.3911  -3.3564  -3.3238  -3.1788  -3.0572  -2.9520  -2.8588  -1.4512  -1.3750  -1.2892  -1.1890  -1.0703  -1.0437  -1.0160  -0.9869  -0.9564  -0.9241  -0.8903  -0.8540  -0.8150  -0.7737  -0.7294  -0.6809  -0.6277  -0.5686  -0.5019  -0.4246  -0.3317  -0.2155  -0.0573   0.2011   0.7686   % 15
                           -5.6626  -4.8453  -4.4979  -4.2923  -4.1434  -4.0277  -3.9310  -3.8486  -3.7773  -3.7139  -3.6569  -3.6050  -3.5570  -3.5130  -3.4716  -3.4328  -3.3965  -3.3624  -3.3296  -3.2984  -3.2686  -3.1366  -3.0250  -2.9281  -2.8416  -1.4854  -1.4118  -1.3285  -1.2325  -1.1181  -1.0925  -1.0650  -1.0368  -1.0072  -0.9760  -0.9428  -0.9077  -0.8702  -0.8301  -0.7867  -0.7391  -0.6874  -0.6293  -0.5636  -0.4883  -0.3980  -0.2837  -0.1298   0.1207   0.6610   % 20
                           -5.4014  -4.6925  -4.3755  -4.1869  -4.0526  -3.9455  -3.8576  -3.7830  -3.7168  -3.6572  -3.6040  -3.5557  -3.5109  -3.4694  -3.4305  -3.3942  -3.3597  -3.3270  -3.2957  -3.2663  -3.2381  -3.1133  -3.0073  -2.9150  -2.8323  -1.5067  -1.4339  -1.3528  -1.2596  -1.1474  -1.1218  -1.0954  -1.0672  -1.0380  -1.0071  -0.9745  -0.9398  -0.9026  -0.8626  -0.8201  -0.7740  -0.7223  -0.6651  -0.6002  -0.5259  -0.4364  -0.3252  -0.1732   0.0740   0.6055   % 25
                           -5.2503  -4.5873  -4.2943  -4.1198  -3.9927  -3.8912  -3.8077  -3.7369  -3.6740  -3.6181  -3.5670  -3.5206  -3.4778  -3.4376  -3.4007  -3.3660  -3.3329  -3.3016  -3.2718  -3.2437  -3.2167  -3.0959  -2.9939  -2.9043  -2.8236  -1.5197  -1.4482  -1.3678  -1.2751  -1.1642  -1.1389  -1.1128  -1.0852  -1.0564  -1.0257  -0.9932  -0.9590  -0.9226  -0.8832  -0.8409  -0.7948  -0.7439  -0.6870  -0.6231  -0.5482  -0.4590  -0.3485  -0.1986   0.0463   0.5630   % 30
                           -5.0682  -4.4754  -4.2049  -4.0416  -3.9233  -3.8299  -3.7523  -3.6851  -3.6264  -3.5739  -3.5262  -3.4822  -3.4420  -3.4044  -3.3695  -3.3365  -3.3053  -3.2755  -3.2472  -3.2202  -3.1943  -3.0785  -2.9804  -2.8940  -2.8163  -1.5365  -1.4655  -1.3863  -1.2946  -1.1856  -1.1611  -1.1353  -1.1084  -1.0800  -1.0498  -1.0183  -0.9841  -0.9477  -0.9087  -0.8667  -0.8209  -0.7709  -0.7148  -0.6513  -0.5782  -0.4901  -0.3806  -0.2330   0.0091   0.5156   % 40
                           -4.9607  -4.4031  -4.1492  -3.9946  -3.8827  -3.7933  -3.7195  -3.6548  -3.5983  -3.5480  -3.5022  -3.4595  -3.4200  -3.3835  -3.3497  -3.3175  -3.2875  -3.2591  -3.2320  -3.2061  -3.1812  -3.0692  -2.9731  -2.8883  -2.8119  -1.5458  -1.4756  -1.3966  -1.3060  -1.1988  -1.1745  -1.1490  -1.1222  -1.0940  -1.0641  -1.0322  -0.9988  -0.9630  -0.9245  -0.8823  -0.8370  -0.7876  -0.7324  -0.6686  -0.5955  -0.5079  -0.3995  -0.2536  -0.0160   0.4885   % 50
                           -4.8378  -4.3225  -4.0857  -3.9406  -3.8339  -3.7497  -3.6787  -3.6175  -3.5630  -3.5146  -3.4707  -3.4302  -3.3930  -3.3583  -3.3254  -3.2947  -3.2658  -3.2377  -3.2113  -3.1861  -3.1618  -3.0538  -2.9610  -2.8794  -2.8055  -1.5584  -1.4886  -1.4109  -1.3214  -1.2146  -1.1909  -1.1657  -1.1390  -1.1112  -1.0817  -1.0504  -1.0170  -0.9817  -0.9430  -0.9017  -0.8565  -0.8074  -0.7515  -0.6889  -0.6160  -0.5297  -0.4225  -0.2760  -0.0407   0.4522   % 75
                           -4.7673  -4.2841  -4.0541  -3.9146  -3.8101  -3.7275  -3.6579  -3.5983  -3.5456  -3.4982  -3.4555  -3.4163  -3.3796  -3.3455  -3.3136  -3.2833  -3.2550  -3.2279  -3.2021  -3.1774  -3.1536  -3.0475  -2.9564  -2.8756  -2.8025  -1.5644  -1.4951  -1.4176  -1.3281  -1.2222  -1.1982  -1.1735  -1.1471  -1.1192  -1.0900  -1.0587  -1.0260  -0.9906  -0.9524  -0.9110  -0.8662  -0.8170  -0.7619  -0.6990  -0.6274  -0.5413  -0.4343  -0.2896  -0.0562   0.4339   % 100
                           -4.7130  -4.2426  -4.0230  -3.8854  -3.7857  -3.7059  -3.6384  -3.5799  -3.5288  -3.4827  -3.4407  -3.4023  -3.3665  -3.3332  -3.3019  -3.2723  -3.2444  -3.2176  -3.1922  -3.1680  -3.1447  -3.0397  -2.9499  -2.8703  -2.7982  -1.5707  -1.5014  -1.4245  -1.3358  -1.2307  -1.2069  -1.1821  -1.1554  -1.1276  -1.0987  -1.0678  -1.0351  -0.9999  -0.9618  -0.9206  -0.8764  -0.8267  -0.7723  -0.7101  -0.6381  -0.5519  -0.4453  -0.3023  -0.0675   0.4149   % 150
                           -4.6894  -4.2217  -4.0055  -3.8719  -3.7735  -3.6947  -3.6291  -3.5715  -3.5209  -3.4758  -3.4342  -3.3956  -3.3604  -3.3273  -3.2963  -3.2672  -3.2397  -3.2133  -3.1882  -3.1640  -3.1406  -3.0374  -2.9482  -2.8694  -2.7980  -1.5743  -1.5057  -1.4284  -1.3400  -1.2352  -1.2114  -1.1868  -1.1608  -1.1331  -1.1041  -1.0732  -1.0403  -1.0056  -0.9676  -0.9265  -0.8818  -0.8332  -0.7784  -0.7166  -0.6452  -0.5599  -0.4531  -0.3078  -0.0775   0.4119   % 200
                           -4.6604  -4.2049  -3.9903  -3.8583  -3.7612  -3.6831  -3.6180  -3.5613  -3.5118  -3.4668  -3.4259  -3.3881  -3.3534  -3.3209  -3.2903  -3.2613  -3.2338  -3.2078  -3.1831  -3.1591  -3.1359  -3.0332  -2.9449  -2.8671  -2.7965  -1.5775  -1.5088  -1.4320  -1.3442  -1.2396  -1.2162  -1.1914  -1.1654  -1.1385  -1.1093  -1.0785  -1.0460  -1.0108  -0.9732  -0.9324  -0.8879  -0.8382  -0.7838  -0.7222  -0.6513  -0.5661  -0.4602  -0.3163  -0.0848   0.3919   % 300
                           -4.6251  -4.1859  -3.9778  -3.8478  -3.7519  -3.6744  -3.6097  -3.5542  -3.5043  -3.4600  -3.4193  -3.3821  -3.3475  -3.3154  -3.2849  -3.2561  -3.2288  -3.2032  -3.1784  -3.1548  -3.1325  -3.0306  -2.9428  -2.8650  -2.7947  -1.5799  -1.5112  -1.4346  -1.3469  -1.2421  -1.2183  -1.1936  -1.1677  -1.1402  -1.1110  -1.0801  -1.0471  -1.0121  -0.9745  -0.9340  -0.8895  -0.8408  -0.7865  -0.7249  -0.6531  -0.5680  -0.4608  -0.3176  -0.0871   0.3998   % 500
                           -4.6152  -4.1778  -3.9692  -3.8399  -3.7446  -3.6689  -3.6040  -3.5483  -3.4992  -3.4549  -3.4145  -3.3778  -3.3434  -3.3108  -3.2804  -3.2518  -3.2248  -3.1990  -3.1746  -3.1510  -3.1284  -3.0276  -2.9406  -2.8631  -2.7929  -1.5806  -1.5122  -1.4354  -1.3478  -1.2429  -1.2195  -1.1946  -1.1686  -1.1413  -1.1125  -1.0821  -1.0490  -1.0140  -0.9762  -0.9355  -0.8914  -0.8424  -0.7887  -0.7271  -0.6556  -0.5713  -0.4656  -0.3213  -0.0912   0.3952   % 1000
                           -4.5959  -4.1652  -3.9592  -3.8340  -3.7389  -3.6632  -3.5991  -3.5441  -3.4961  -3.4525  -3.4123  -3.3750  -3.3409  -3.3089  -3.2787  -3.2504  -3.2235  -3.1982  -3.1739  -3.1504  -3.1277  -3.0274  -2.9405  -2.8634  -2.7934  -1.5826  -1.5141  -1.4378  -1.3502  -1.2466  -1.2229  -1.1982  -1.1726  -1.1452  -1.1162  -1.0856  -1.0532  -1.0185  -0.9811  -0.9405  -0.8964  -0.8470  -0.7930  -0.7314  -0.6591  -0.5737  -0.4678  -0.3240  -0.0957   0.3836]; % 10000
                 
            case 'T2'
                
              % Row headers (effective sample sizes) are at the right of the table.
              % ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
              % Alpha        0.001    0.005    0.010    0.015    0.020    0.025    0.030    0.035    0.040    0.045    0.050    0.055    0.060    0.065    0.070    0.075    0.080    0.085    0.090    0.095    0.100    0.125    0.150    0.175    0.200    0.800    0.825    0.850    0.875    0.900    0.905    0.910    0.915    0.920    0.925    0.930    0.935    0.940    0.945    0.950    0.955    0.960    0.965    0.970    0.975    0.980    0.985    0.990    0.995    0.999      T
              % ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
                CVTable = [-18.5945 -16.5780 -15.6998 -15.1440 -14.7287 -14.3914 -14.1039 -13.8514 -13.6277 -13.4232 -13.2361 -13.0611 -12.8982 -12.7449 -12.6009 -12.4647 -12.3350 -12.2100 -12.0917 -11.9780 -11.8680 -11.3698 -10.9370 -10.5518 -10.2007  -4.4678  -4.1798  -3.8584  -3.4979  -3.0833  -2.9917  -2.8963  -2.7976  -2.6947  -2.5875  -2.4742  -2.3542  -2.2245  -2.0875  -1.9389  -1.7786  -1.6017  -1.4051  -1.1825  -0.9215  -0.6055  -0.2080   0.3485   1.3098   3.6765   % 10
                           -22.2637 -19.8543 -18.6433 -17.8746 -17.2993 -16.8339 -16.4386 -16.0986 -15.7942 -15.5190 -15.2679 -15.0352 -14.8201 -14.6177 -14.4258 -14.2459 -14.0755 -13.9138 -13.7596 -13.6104 -13.4668 -12.8261 -12.2756 -11.7891 -11.3511  -4.7183  -4.4230  -4.1028  -3.7476  -3.3398  -3.2503  -3.1566  -3.0601  -2.9580  -2.8529  -2.7417  -2.6247  -2.5021  -2.3692  -2.2276  -2.0740  -1.9082  -1.7230  -1.5171  -1.2800  -0.9995  -0.6483  -0.1722   0.6109   2.3703   % 15
                           -25.1513 -22.0941 -20.5859 -19.6474 -18.9513 -18.3898 -17.9150 -17.5057 -17.1449 -16.8212 -16.5259 -16.2522 -16.0007 -15.7667 -15.5453 -15.3363 -15.1381 -14.9510 -14.7718 -14.5989 -14.4326 -13.6967 -13.0717 -12.5219 -12.0279  -4.8461  -4.5420  -4.2140  -3.8564  -3.4520  -3.3622  -3.2708  -3.1740  -3.0742  -2.9694  -2.8608  -2.7459  -2.6244  -2.4936  -2.3535  -2.2030  -2.0392  -1.8605  -1.6601  -1.4330  -1.1623  -0.8267  -0.3768   0.3499   1.9516   % 20
                           -27.2584 -23.7057 -21.9768 -20.9009 -20.1125 -19.4778 -18.9466 -18.4848 -18.0806 -17.7165 -17.3878 -17.0821 -16.8030 -16.5429 -16.2992 -16.0663 -15.8495 -15.6409 -15.4447 -15.2573 -15.0768 -14.2736 -13.5967 -13.0044 -12.4707  -4.9313  -4.6207  -4.2884  -3.9279  -3.5213  -3.4311  -3.3391  -3.2428  -3.1428  -3.0391  -2.9306  -2.8160  -2.6940  -2.5661  -2.4294  -2.2812  -2.1193  -1.9416  -1.7435  -1.5194  -1.2555  -0.9286  -0.4925   0.2106   1.7301   % 25
                           -28.8695 -24.8459 -22.9626 -21.7970 -20.9337 -20.2474 -19.6713 -19.1707 -18.7353 -18.3427 -17.9883 -17.6640 -17.3634 -17.0829 -16.8209 -16.5738 -16.3421 -16.1210 -15.9140 -15.7129 -15.5209 -14.6740 -13.9569 -13.3296 -12.7704  -4.9871  -4.6719  -4.3360  -3.9710  -3.5609  -3.4710  -3.3783  -3.2838  -3.1840  -3.0797  -2.9709  -2.8554  -2.7357  -2.6071  -2.4692  -2.3230  -2.1637  -1.9894  -1.7910  -1.5670  -1.3066  -0.9853  -0.5571   0.1288   1.5797   % 30
                           -31.1160 -26.5053 -24.3460 -23.0288 -22.0699 -21.3026 -20.6667 -20.1090 -19.6277 -19.1994 -18.8064 -18.4518 -18.1254 -17.8243 -17.5405 -17.2759 -17.0258 -16.7866 -16.5620 -16.3441 -16.1383 -15.2255 -14.4554 -13.7878 -13.1925  -5.0574  -4.7359  -4.3937  -4.0234  -3.6116  -3.5218  -3.4291  -3.3338  -3.2335  -3.1297  -3.0220  -2.9084  -2.7878  -2.6602  -2.5244  -2.3783  -2.2212  -2.0465  -1.8532  -1.6344  -1.3758  -1.0615  -0.6404   0.0252   1.4097   % 40
                           -32.5622 -27.5182 -25.1961 -23.7865 -22.7760 -21.9667 -21.2873 -20.7068 -20.1939 -19.7414 -19.3314 -18.9623 -18.6191 -18.3004 -18.0025 -17.7219 -17.4590 -17.2090 -16.9742 -16.7494 -16.5317 -15.5775 -14.7753 -14.0772 -13.4562  -5.1025  -4.7776  -4.4338  -4.0591  -3.6437  -3.5533  -3.4609  -3.3652  -3.2656  -3.1623  -3.0548  -2.9415  -2.8214  -2.6944  -2.5595  -2.4130  -2.2554  -2.0825  -1.8905  -1.6718  -1.4149  -1.1006  -0.6919  -0.0436   1.3099   % 50
                           -34.7768 -29.0978 -26.5157 -24.9552 -23.8237 -22.9284 -22.1973 -21.5775 -21.0202 -20.5283 -20.0814 -19.6781 -19.3045 -18.9618 -18.6392 -18.3369 -18.0552 -17.7893 -17.5395 -17.2993 -17.0682 -16.0545 -15.2019 -14.4636 -13.8124  -5.1609  -4.8290  -4.4801  -4.1046  -3.6851  -3.5949  -3.5024  -3.4049  -3.3040  -3.1995  -3.0903  -2.9768  -2.8594  -2.7319  -2.5954  -2.4497  -2.2925  -2.1211  -1.9277  -1.7121  -1.4611  -1.1516  -0.7439  -0.1094   1.1898   % 75
                           -35.8705 -29.9028 -27.1762 -25.5647 -24.3787 -23.4373 -22.6708 -22.0106 -21.4404 -20.9321 -20.4737 -20.0546 -19.6684 -19.3168 -18.9829 -18.6713 -18.3791 -18.1053 -17.8441 -17.5954 -17.3577 -16.3123 -15.4341 -14.6723 -14.0024  -5.1892  -4.8559  -4.5047  -4.1259  -3.7034  -3.6125  -3.5192  -3.4228  -3.3225  -3.2183  -3.1097  -2.9961  -2.8761  -2.7516  -2.6168  -2.4712  -2.3136  -2.1407  -1.9485  -1.7339  -1.4827  -1.1776  -0.7757  -0.1486   1.1302   % 100
                           -37.1431 -30.7593 -27.8713 -26.1711 -24.9387 -23.9680 -23.1674 -22.4859 -21.8910 -21.3599 -20.8840 -20.4448 -20.0427 -19.6754 -19.3313 -19.0079 -18.7070 -18.4159 -18.1449 -17.8856 -17.6387 -16.5588 -15.6547 -14.8758 -14.1881  -5.2197  -4.8821  -4.5271  -4.1460  -3.7236  -3.6322  -3.5380  -3.4411  -3.3402  -3.2365  -3.1279  -3.0143  -2.8949  -2.7698  -2.6353  -2.4905  -2.3351  -2.1627  -1.9718  -1.7550  -1.5052  -1.2004  -0.8025  -0.1763   1.0617   % 150
                           -37.8724 -31.1869 -28.2442 -26.4984 -25.2245 -24.2437 -23.4297 -22.7340 -22.1183 -21.5776 -21.0938 -20.6523 -20.2412 -19.8655 -19.5156 -19.1868 -18.8778 -18.5854 -18.3114 -18.0489 -17.8022 -16.6940 -15.7799 -14.9905 -14.2939  -5.2395  -4.9011  -4.5453  -4.1631  -3.7395  -3.6483  -3.5550  -3.4576  -3.3566  -3.2528  -3.1437  -3.0301  -2.9090  -2.7823  -2.6468  -2.5026  -2.3463  -2.1756  -1.9867  -1.7719  -1.5214  -1.2180  -0.8164  -0.2013   1.0505   % 200
                           -38.5337 -31.6601 -28.6382 -26.8079 -25.5187 -24.5046 -23.6752 -22.9618 -22.3384 -21.7871 -21.2925 -20.8396 -20.4247 -20.0411 -19.6854 -19.3522 -19.0372 -18.7406 -18.4615 -18.1949 -17.9439 -16.8297 -15.8984 -15.0947 -14.3895  -5.2546  -4.9151  -4.5578  -4.1748  -3.7494  -3.6572  -3.5643  -3.4658  -3.3658  -3.2614  -3.1519  -3.0374  -2.9181  -2.7912  -2.6584  -2.5133  -2.3587  -2.1877  -1.9982  -1.7840  -1.5351  -1.2335  -0.8379  -0.2195   0.9992   % 300
                           -38.9555 -31.9802 -28.8935 -27.0725 -25.7621 -24.7377 -23.8800 -23.1514 -22.5164 -21.9553 -21.4547 -20.9929 -20.5701 -20.1878 -19.8245 -19.4878 -19.1687 -18.8685 -18.5858 -18.3164 -18.0590 -16.9285 -15.9907 -15.1796 -14.4662  -5.2666  -4.9261  -4.5668  -4.1809  -3.7563  -3.6638  -3.5701  -3.4722  -3.3709  -3.2662  -3.1563  -3.0418  -2.9232  -2.7957  -2.6601  -2.5155  -2.3596  -2.1871  -1.9980  -1.7840  -1.5361  -1.2336  -0.8390  -0.2249   1.0135   % 500
                           -39.4974 -32.2665 -29.1380 -27.2701 -25.9420 -24.8840 -24.0208 -23.2914 -22.6493 -22.0786 -21.5733 -21.1094 -20.6872 -20.2966 -19.9295 -19.5875 -19.2651 -18.9645 -18.6743 -18.4045 -18.1475 -17.0134 -16.0612 -15.2362 -14.5177  -5.2707  -4.9275  -4.5686  -4.1835  -3.7594  -3.6677  -3.5734  -3.4761  -3.3754  -3.2709  -3.1605  -3.0448  -2.9256  -2.7976  -2.6624  -2.5181  -2.3614  -2.1927  -2.0030  -1.7914  -1.5444  -1.2417  -0.8453  -0.2359   1.0024   % 1000
                           -39.7292 -32.5214 -29.3358 -27.4684 -26.1172 -25.0678 -24.1936 -23.4529 -22.8026 -22.2219 -21.7066 -21.2426 -20.8101 -20.4155 -20.0461 -19.6993 -19.3748 -19.0695 -18.7835 -18.5064 -18.2441 -17.0924 -16.1340 -15.3135 -14.5868  -5.2802  -4.9392  -4.5794  -4.1924  -3.7667  -3.6759  -3.5805  -3.4836  -3.3833  -3.2787  -3.1696  -3.0554  -2.9349  -2.8092  -2.6747  -2.5301  -2.3725  -2.2014  -2.0122  -1.7980  -1.5529  -1.2478  -0.8524  -0.2459   0.9738]; % 10000
                               
            case 'F'
                        
              % Row headers (effective sample sizes) are at the right of the table.
              % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
              % Alpha      0.001    0.005    0.010    0.015    0.020    0.025    0.030    0.035    0.040    0.045    0.050    0.055    0.060    0.065    0.070    0.075    0.080    0.085    0.090    0.095    0.100    0.125    0.150    0.175    0.200    0.800    0.825    0.850    0.875    0.900    0.905    0.910    0.915    0.920    0.925    0.930    0.935    0.940    0.945    0.950    0.955    0.960    0.965    0.970    0.975    0.980    0.985    0.990    0.995    0.999      T
              % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   -----
               CVTable = [35.6551  21.5349  17.1540  14.9466  13.5352  12.4893  11.6811  11.0307  10.4865  10.0206   9.6144   9.2592   8.9445   8.6624   8.4016   8.1652   7.9462   7.7457   7.5593   7.3844   7.2192   6.5287   5.9898   5.5496   5.1792   1.6196   1.5270   1.4332   1.3363   1.2347   1.2137   1.1924   1.1705   1.1484   1.1256   1.1022   1.0783   1.0536   1.0279   1.0013   0.9736   0.9444   0.9132   0.8797   0.8432   0.8030   0.7571   0.7023   0.6289   0.5224   % 10
                          21.6618  15.0130  12.6637  11.4058  10.5597   9.9305   9.4285   9.0127   8.6581   8.3548   8.0844   7.8437   7.6275   7.4303   7.2472   7.0785   6.9227   6.7780   6.6439   6.5164   6.3983   5.8860   5.4739   5.1311   4.8351   1.6846   1.5912   1.4961   1.3970   1.2924   1.2705   1.2482   1.2254   1.2020   1.1779   1.1535   1.1281   1.1019   1.0750   1.0466   1.0169   0.9860   0.9531   0.9175   0.8787   0.8354   0.7860   0.7251   0.6426   0.5148   % 15
                          17.5668  12.9524  11.2049  10.2269   9.5535   9.0407   8.6276   8.2866   7.9983   7.7421   7.5147   7.3116   7.1292   6.9626   6.8099   6.6675   6.5339   6.4102   6.2927   6.1826   6.0782   5.6292   5.2645   4.9582   4.6925   1.7132   1.6205   1.5247   1.4245   1.3179   1.2958   1.2732   1.2501   1.2266   1.2024   1.1773   1.1516   1.1249   1.0971   1.0683   1.0379   1.0057   0.9716   0.9351   0.8949   0.8501   0.7988   0.7364   0.6484   0.5121   % 20
                          15.7964  11.9921  10.4725   9.6256   9.0302   8.5824   8.2196   7.9133   7.6525   7.4274   7.2242   7.0410   6.8739   6.7207   6.5807   6.4498   6.3276   6.2139   6.1062   6.0055   5.9091   5.4924   5.1541   4.8678   4.6186   1.7313   1.6384   1.5425   1.4418   1.3351   1.3127   1.2898   1.2662   1.2421   1.2174   1.1921   1.1658   1.1385   1.1102   1.0810   1.0501   1.0174   0.9825   0.9458   0.9053   0.8596   0.8064   0.7413   0.6513   0.5092   % 25
                          14.7731  11.3892  10.0170   9.2406   8.7051   8.2932   7.9602   7.6739   7.4309   7.2215   7.0326   6.8609   6.7049   6.5624   6.4294   6.3053   6.1919   6.0840   5.9822   5.8874   5.7968   5.4023   5.0792   4.8041   4.5654   1.7416   1.6489   1.5532   1.4526   1.3453   1.3227   1.2998   1.2762   1.2519   1.2271   1.2017   1.1754   1.1484   1.1200   1.0904   1.0592   1.0263   0.9914   0.9531   0.9119   0.8654   0.8116   0.7462   0.6547   0.5088   % 30
                          13.6827  10.7326   9.5301   8.8400   8.3519   7.9749   7.6690   7.4131   7.1914   6.9952   6.8215   6.6634   6.5185   6.3859   6.2629   6.1486   6.0413   5.9410   5.8451   5.7564   5.6707   5.2990   4.9921   4.7329   4.5068   1.7557   1.6634   1.5674   1.4670   1.3583   1.3357   1.3127   1.2890   1.2646   1.2394   1.2137   1.1868   1.1594   1.1310   1.1009   1.0692   1.0360   1.0006   0.9623   0.9206   0.8737   0.8192   0.7525   0.6589   0.5088   % 40
                          13.0368  10.3560   9.2462   8.5976   8.1406   7.7862   7.5020   7.2614   7.0513   6.8653   6.6970   6.5479   6.4116   6.2859   6.1675   6.0567   5.9544   5.8592   5.7683   5.6819   5.6004   5.2450   4.9502   4.6980   4.4768   1.7638   1.6719   1.5760   1.4750   1.3667   1.3438   1.3204   1.2964   1.2720   1.2472   1.2214   1.1946   1.1672   1.1382   1.1082   1.0765   1.0432   1.0069   0.9683   0.9258   0.8783   0.8235   0.7560   0.6604   0.5081   % 50
                          12.3174   9.9229   8.9004   8.3099   7.8902   7.5651   7.2985   7.0716   6.8764   6.7019   6.5451   6.4045   6.2748   6.1549   6.0442   5.9415   5.8447   5.7541   5.6673   5.5853   5.5075   5.1658   4.8837   4.6429   4.4308   1.7737   1.6813   1.5856   1.4848   1.3764   1.3536   1.3300   1.3062   1.2815   1.2562   1.2303   1.2033   1.1753   1.1463   1.1157   1.0836   1.0499   1.0138   0.9746   0.9319   0.8838   0.8285   0.7599   0.6628   0.5080   % 75
                          11.9200   9.7120   8.7429   8.1780   7.7731   7.4547   7.1983   6.9794   6.7895   6.6215   6.4702   6.3337   6.2080   6.0928   5.9851   5.8844   5.7905   5.7014   5.6182   5.5384   5.4625   5.1314   4.8561   4.6195   4.4114   1.7796   1.6871   1.5913   1.4901   1.3809   1.3581   1.3345   1.3103   1.2854   1.2601   1.2339   1.2068   1.1790   1.1497   1.1196   1.0869   1.0528   1.0165   0.9768   0.9335   0.8848   0.8282   0.7589   0.6625   0.5063   % 100
                          11.6112   9.5005   8.5826   8.0400   7.6515   7.3507   7.1040   6.8927   6.7099   6.5456   6.3982   6.2655   6.1438   6.0302   5.9256   5.8292   5.7373   5.6509   5.5691   5.4920   5.4178   5.0926   4.8235   4.5921   4.3870   1.7854   1.6933   1.5971   1.4959   1.3869   1.3639   1.3405   1.3164   1.2919   1.2664   1.2403   1.2131   1.1848   1.1556   1.1245   1.0920   1.0578   1.0213   0.9821   0.9386   0.8900   0.8334   0.7633   0.6650   0.5080   % 150
                          11.4893   9.3975   8.5066   7.9775   7.5960   7.3008   7.0562   6.8510   6.6699   6.5092   6.3655   6.2358   6.1162   6.0050   5.9007   5.8035   5.7137   5.6287   5.5478   5.4713   5.3991   5.0775   4.8110   4.5815   4.3798   1.7884   1.6962   1.6004   1.4995   1.3901   1.3668   1.3435   1.3195   1.2945   1.2691   1.2425   1.2153   1.1870   1.1578   1.1271   1.0952   1.0609   1.0240   0.9844   0.9411   0.8924   0.8351   0.7643   0.6658   0.5049   % 200
                          11.3273   9.3053   8.4303   7.9077   7.5350   7.2463   7.0071   6.8033   6.6286   6.4720   6.3310   6.2025   6.0835   5.9736   5.8721   5.7773   5.6886   5.6045   5.5255   5.4495   5.3764   5.0593   4.7950   4.5682   4.3685   1.7915   1.6992   1.6029   1.5014   1.3919   1.3687   1.3449   1.3205   1.2954   1.2700   1.2437   1.2164   1.1883   1.1588   1.1281   1.0958   1.0613   1.0246   0.9849   0.9416   0.8923   0.8350   0.7643   0.6655   0.5064   % 300
                          11.1903   9.2247   8.3667   7.8579   7.4898   7.2063   6.9699   6.7686   6.5935   6.4402   6.2991   6.1718   6.0552   5.9472   5.8458   5.7524   5.6650   5.5833   5.5043   5.4294   5.3577   5.0435   4.7825   4.5575   4.3590   1.7927   1.7006   1.6044   1.5033   1.3942   1.3710   1.3474   1.3233   1.2986   1.2728   1.2464   1.2191   1.1905   1.1609   1.1299   1.0974   1.0631   1.0264   0.9863   0.9428   0.8938   0.8366   0.7660   0.6663   0.5051   % 500
                          11.0965   9.1788   8.3253   7.8192   7.4601   7.1763   6.9434   6.7433   6.5700   6.4163   6.2775   6.1518   6.0370   5.9301   5.8304   5.7369   5.6493   5.5658   5.4873   5.4130   5.3422   5.0313   4.7718   4.5487   4.3515   1.7938   1.7012   1.6048   1.5037   1.3946   1.3715   1.3478   1.3235   1.2985   1.2728   1.2466   1.2194   1.1909   1.1615   1.1305   1.0980   1.0637   1.0266   0.9868   0.9427   0.8929   0.8355   0.7658   0.6661   0.5047   % 1000
                          10.9918   9.1146   8.2838   7.7823   7.4278   7.1492   6.9191   6.7227   6.5519   6.3992   6.2634   6.1385   6.0244   5.9198   5.8203   5.7285   5.6417   5.5594   5.4812   5.4068   5.3362   5.0265   4.7686   4.5454   4.3489   1.7960   1.7037   1.6076   1.5063   1.3970   1.3737   1.3503   1.3258   1.3008   1.2753   1.2487   1.2215   1.1933   1.1637   1.1327   1.1003   1.0656   1.0290   0.9892   0.9457   0.8961   0.8382   0.7677   0.6669   0.5051]; % 10000
                       
        end
        
end

%-------------------------------------------------------------------------
% Compute the statistic:
function [testStat,testPValue] = getStat(i,testT,testLags,testModel,testType,testReg,sigLevels,sampSizes,CVTable,needPValue)

beta = testReg.coeff;
Cov = testReg.Cov;
aId = strncmp('a',testReg.names,1);
a = beta(aId);
bId = strncmp('b',testReg.names,1);
b = beta(bId);
se = sqrt(diag(Cov));
se_a = se(aId);
        
switch upper(testType)
    
    case 'T1'

        testStat = (a-1)/se_a;
         
    case 'T2'

        testStat = testT*(a-1)/(1-sum(b));
        
    case 'F'
        
        % Follows Hamilton [5] pp.521-527.
        %        
        % Restrictions of the form:
        % R*(beta-beta0) = R*beta - R*beta0
        %                = R*beta - r
        %                = 0
        switch upper(testModel)
                       
            case 'ARD'
                
                R = [eye(2),zeros(2,testLags)];
                beta0 = [0;1;zeros(testLags,1)];
                
            case 'TS'
                
                R = [zeros(2,1),eye(2),zeros(2,testLags)];
                beta0 = [0;0;1;zeros(testLags,1)];
                
        end
        
        % The inverse in the Wald form of the test statistic is replaced
        % with a more efficient computation using Cholesky factorization
        % and triangular backsolve:
        %
        % U = chol(R*Cov*R')
        % u = R*(beta-beta0)
        % testStat = (u'*inv(R*Cov*R')*u)/2
        %          = (u'*inv(U'*U)*u)/2
        %          = (u'*inv(U)*inv(U')*u)/2
        %          = ((u'/U)*(U'\u))/2
        %          = (v'*v)/2
        % where v = U'\u.
        
        U = chol(R*Cov*R');
        u = R*(beta-beta0);
        v = U'\u;
        testStat = (v'*v)/2;
        
end

if needPValue
    testPValue = getPValue(i,sigLevels,sampSizes,CVTable,testT,testType,testStat);
else
    testPValue = NaN;
end

%-------------------------------------------------------------------------
% Get p-values
function testPValue = getPValue(i,sigLevels,sampSizes,CVTable,testT,testType,testStat)

% P-values are estimated using two successive 1D interpolations:
%
% 1. Find all critical values associated with the effective sample size.
% 2. Find the cumulative probability associated with the test statistic.

CVTableRowT = interp2(sigLevels,sampSizes,CVTable,sigLevels,testT,'linear');

if strcmpi(testType,'F') % Right-tailed test
    
    if testStat >= CVTableRowT(1)
        
        testPValue = sigLevels(1);
        warning(message('econ:adftest:RightTailStatTooBig', i, sprintf( '%5.3f', testPValue )))
            
    elseif testStat <= CVTableRowT(end)
        
        testPValue = sigLevels(end);
        warning(message('econ:adftest:RightTailStatTooSmall', i, sprintf( '%5.3f', testPValue )))
            
    else 
        
        testPValue = interp1(CVTableRowT,sigLevels,testStat,'linear');
            
    end        
    
else % Left-tailed test
    
    if testStat <= CVTableRowT(1)
        
        testPValue = sigLevels(1);
%        warning(message('econ:adftest:LeftTailStatTooSmall', i, sprintf( '%5.3f', testPValue )))
            
    elseif testStat >= CVTableRowT(end)
        
        testPValue = sigLevels(end);
%        warning(message('econ:adftest:LeftTailStatTooBig', i, sprintf( '%5.3f', testPValue )))
            
    else
        
        testPValue = interp1(CVTableRowT,sigLevels,testStat,'linear');
            
    end
    
end

%-------------------------------------------------------------------------
% Test the statistic:
function [testCValue,testH] = runTest(sigLevels,sampSizes,testT,testType,testAlpha,testStat,CVTable)

testCValue = interp2(sigLevels,sampSizes,CVTable,testAlpha,testT,'linear');

switch upper(testType)
    
    case {'T1','T2'}        

        testH = (testStat < testCValue); % Left-tailed test
        
    case 'F'

        testH = (testStat > testCValue); % Right-tailed test
    
end