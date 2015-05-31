function [h,pValue,stat,cValue,mles] = jcitest(Y,varargin)
%JCITEST Johansen cointegration test
%
% Syntax:
%
%   [h,pValue,stat,cValue,mles] = jcitest(Y)
%   [h,pValue,stat,cValue,mles] = jcitest(Y,param,val,...)
%
% Description:
%
%   Johansen tests assess the null hypothesis H(r) of cointegration rank
%   less than or equal to r among the numDims-dimensional time series in Y
%   against alternatives H(numDims) (trace test) or H(r+1) (maxeig test).
%   The tests also produce maximum likelihood estimates of the parameters
%   in a vector error-correction (VEC) model of the cointegrated series.
%
% Input Arguments:
%
%   Y - numObs-by-numDims matrix representing numObs observations of a
%       numDims-dimensional time series y(t), with the last observation the
%       most recent. Observations containing NaN values are removed.
%       Initial values for lagged variables in VEC model estimation are
%       taken from the beginning of the data.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME        VALUE
%
%   'model' 	String or cell vector of strings specifying the form of the
%               deterministic components of the VEC(q) model of y(t):
%
%                 (1-L)y(t) = C*y(t-1)
%                             + B1*(1-L)y(t-1) + ... + Bq*(1-L)y(t-q)
%                             + D*x
%                             + e(t)
%
%               where L is the lag operator Ly(t) = y(t-1). If r < numDims
%               is the cointegration rank, then C = A*B' where A is a
%               numDims-by-r matrix of error-correction speeds and B is a
%               numDims-by-r matrix of basis vectors for the space of
%               cointegrating relations. x contains any exogenous terms
%               representing deterministic trends in the data. For maximum
%               likelihood estimation, it is assumed that e(t) ~ NID(0,Q),
%               where Q is the innovations covariance matrix.
%
%               Values of 'model' are those considered by Johansen [2]:
%
%               VALUE       FORM OF A*B'*y(t-1) + D*x
%
%               o 'H2'      A*B'*y(t-1). There are no intercepts or trends
%                           in the cointegrating relations and there are no
%                           trends in the data. This model is only
%                           appropriate if all series have zero mean.
%
%               o 'H1*'     A*(B'*y(t-1) + c0). There are intercepts in the
%                         	cointegrating relations and there are no
%                           trends in the data. This model is appropriate
%                           for nontrending data with nonzero mean.
%
%               o 'H1'      A*(B'*y(t-1) + c0) + c1. There are intercepts
%                           in the cointegrating relations and there are
%                           linear trends in the data. This is a model of
%                           "deterministic cointegration," where the
%                           cointegrating relations eliminate both
%                           stochastic and deterministic trends in the
%                           data. This is the default value.
%
%               o 'H*'      A*(B'*y(t-1) + c0 + d0*t) + c1. There are
%                           intercepts and linear trends in the
%                           cointegrating relations and there are linear
%                           trends in the data. This is a model of
%                           "stochastic cointegration," where the
%                           cointegrating relations eliminate stochastic
%                           but not deterministic trends in the data.
%
%               o 'H'       A*(B'*y(t-1) + c0 + d0*t) + c1 + d1*t. There
%                           are intercepts and linear trends in the
%                           cointegrating relations and there are quadratic
%                           trends in the data. Unless quadratic trends are
%                           actually present in the data, this model may
%                           produce good in-sample fits but poor out-of-
%                           sample forecasts.
%
%               Deterministic terms outside of the cointegrating relations,
%               c1 and d1, are identified by projecting constant and linear
%               regression coefficients, respectively, onto the orthogonal
%               complement of A, as in [2].
%
%   'lags'      Scalar or vector of nonnegative integers indicating the
%               number q of lagged differences in the VEC(q) model of y(t).
%               The default value is 0.
%
%               Lagging and differencing a time series reduce the sample
%               size. Absent any presample values, if y(t) is defined for
%               t = 1:N, then the lagged series y(t-k) is defined for
%               t = k+1:N. Differencing reduces the time base to k+2:N.
%               With q lagged differences, the common time base is q+2:N
%               and the effective sample size is T = N-(q+1).
%
%   'test'      String or cell vector of strings indicating the type of
%               test to be performed. Values are 'trace' or 'maxeig'. The
%               default value is 'trace'. Both tests assess the null
%               hypothesis H(r) of cointegration rank less than or equal to
%               r.  Statistics are computed using the effective sample size
%               T and ordered estimates of the eigenvalues of C = A*B',
%               lambda(1) > ... > lambda(numDims).
%
%               o When the value is 'trace', the alternative hypothesis is
%                 H(numDims). Statistics are:
%
%                   -T*[log(1-lambda(r+1)) + ... + log(1-lambda(numDims))]
%
%               o When the value is 'maxeig', the alternative hypothesis is
%                 H(r+1). Statistics are:
%
%                   -T*log(1-lambda(r+1))
%
%   'alpha'     Scalar or vector of nominal significance levels for the
%               tests. Values must be between 0.001 and 0.999. The default
%               value is 0.05.
%
%	'display'   String or cell vector of strings indicating whether or not
%               to display a summary of test results and parameter
%               estimates in the command window. Values are:
%
%               VALUE           DISPLAY
%
%               o 'off'         No display to the command window. This is
%                               the default if JCITEST is called with only
%                               one output argument (h).
%               
%               o 'summary'     Display a summary of test results. Null
%                               ranks r = 0:numDims-1 are displayed in the
%                               first column of each summary. Multiple
%                               tests are displayed in separate summaries.
%                               This is the default if JCITEST is called
%                               with more than one output argument (that
%                               is, if pValue is computed), and is
%                               unavailable if JCITEST is called with only
%                               one output argument (h).
%
%               o 'params'      Display maximum likelihood estimates of the
%                               parameter values associated with the
%                               reduced-rank VEC(q) model of y(t). This
%                               display is only available if JCITEST is
%                               called with five output arguments (that is,
%                               if mles is computed). Displayed parameter
%                               values are returned in mles.rn(m).paramVals
%                               for null rank r = n and test m.
%
%               o 'full'        Display both 'summary' and 'params'.
%
%   Scalar or single string values are expanded to the length of any vector
%   value (the number of tests). Vector values must have equal length.
%
% Output Arguments:
%
%   All outputs arguments are numTests-by-numDims dataset arrays. Rows
%   correspond to different tests specified by the input parameters,
%   labeled t1, t2, etc. Columns correspond to different maintained
%   cointegration ranks r = 0, ..., numDims-1, labeled r0, r1, etc. To
%   access information for the model associated with null rank r = n and
%   test m, use <argName>.rn(m), where <argName> is the argument name.
%
%   h - Boolean decisions for the tests. Values of h equal to 1 (true)
%       indicate rejection of the null of cointegration rank r in favor of
%       the alternative. Values of h equal to 0 (false) indicate a failure
%       to reject the null.
%
%   pValue - Right-tail probabilities of the test statistics.
%
%   stat - Test statistics, determined by the 'test' parameter.
%
%   cValue - Critical values for right-tail probabilities, determined by
%       the 'alpha' parameter. JCITEST loads tables of critical values from
%       the file Data_JCITest.mat, then linearly interpolates test critical
%       values from the tables. Tabulated values were computed using
%       methods described in [3].
%
%   mles - Structures of maximum likelihood estimates associated with the
%       VEC(q) model of y(t). Each structure has the following fields:
%
%       paramNames  Cell vector of parameter names, of the form:
%
%                   	{'A','B','B1',...,'Bq','c0','d0','c1','d1'}'.
%
%                   Elements depend on the values of 'lags' and 'model'.
%
%       paramVals   Structure of parameter estimates with field names
%                   corresponding to the parameter names in paramNames.
%
%       res         T-by-numDims matrix of residuals, where T is the
%                   effective sample size, obtained by fitting the VEC(q)
%                   model of y(t) to the input data.
%
%       EstCov      Estimated covariance Q of the innovations process e(t).
%
%       eigVal      Eigenvalue associated with H(r+1).
%
%       eigVec      Eigenvector associated with the eigenvalue in eigVal.
%                   Eigenvectors v are normalized so that v'*S11*v = 1,
%                   where S11 is defined as in [2].
%
%       rLL         Restricted loglikelihood of Y under the null.
%
%       uLL         Unrestricted loglikelihood of Y under the alternative.
%
% Notes:
%
%   o If JCITEST fails to reject the null of cointegration rank r = 0, the
%     inference is that the error-correction coefficient C is zero, and the
%     VEC(q) model reduces to a standard VAR(q) model in first differences.
%     If JCITEST rejects all cointegration ranks r less than numDims, the
%     inference is that C has full rank, and y(t) is stationary in levels.
%
%   o The parameters A and B in the reduced-rank VEC(q) model are not
%     uniquely identified, though their product C = A*B' is. JCITEST
%     constructs B = V(:,1:r) using the eigenvector matrix V, normalized so
%     that V'*S11*V = I, as in [2]. This normalization produces a simpler
%     formula for the MLE of A.
%
%   o To test linear constraints on the error-correction speeds A and the
%     space of cointegrating relations spanned by B, use JCONTEST.
%
%   o Stationary series, associated with standard unit vectors in the space
%     of cointegrating relations, can be removed from cointegration
%     analysis. Individual series can be pretested using ADFTEST, PPTEST,
%     KPSSTEST, and LMCTEST. As an alternative, JCONTEST can test for
%     standard unit vectors in the context of the full model.
%
%   o To convert VEC(q) model parameters in the mles output to VAR(q+1)
%     model parameters, use VECTOVAR.
%
% Example:
%
%   % Data on the term structure of interest rates in Canada:
% 
%   load Data_Canada
%   Y = Data(:,3:end);
%   names = series(3:end);
%   plot(dates,Y)
%   legend(names,'location','NW')
%   grid on
% 
%   % Test for cointegration:
% 
%   [h,pValue,stat,cValue,mles] = jcitest(Y,'model','H1');
% 
%   % Plot the estimated cointegrating relations B'*y(t-1)+c0:
%   
%   YLag = Y(2:end,:);
%   T = size(YLag,1);
%   B = mles.r2.paramVals.B;
%   c0 = mles.r2.paramVals.c0;  
%   plot(dates(2:end),YLag*B+repmat(c0',T,1))
%   grid on
%
% References:
% 
%   [1] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
%   [2] Johansen, S. Likelihood-Based Inference in Cointegrated Vector
%       Autoregressive Models. Oxford: Oxford University Press, 1995.
%
%   [3] MacKinnon, J. G., A. A. Haug, and L. Michelis. â€œNumerical
%       Distribution Functions of Likelihood Ratio Tests for
%       Cointegration.â€? Journal of Applied Econometrics. v. 14, 1999,
%       pp. 563-577.
%
%   [4] Turner, P. M. "Testing for Cointegration Using the Johansen
%       Approach: Are We Using the Correct Critical Values?" Journal of
%       Applied Econometrics. v. 24, 2009, pp. 825-831.
%
% See also EGCITEST, JCONTEST, VECTOVAR.

% Copyright 2011 The MathWorks, Inc.

% Parse inputs and set defaults:

parseObj = inputParser;
parseObj.addRequired('Y',@YCheck);
parseObj.addParamValue('model','H1',@modelCheck);
parseObj.addParamValue('lags',0,@lagsCheck);
parseObj.addParamValue('test','trace',@testCheck);
parseObj.addParamValue('alpha',0.05,@alphaCheck);
parseObj.addParamValue('display','summary',@displayCheck);

parseObj.parse(Y,varargin{:});

Y = parseObj.Results.Y;
model = upper(parseObj.Results.model);
lags = parseObj.Results.lags;
test = lower(parseObj.Results.test);
alpha = parseObj.Results.alpha;
display = lower(parseObj.Results.display);

% Check parameter values for commensurate lengths, expand single-element
% values, and convert all variables to columns:

[numTests,model,lags,test,alpha,display] = sizeCheck(model,lags,test,alpha,display);

% Capture data name for display:

dataName = inputname(1);

% Turn off display if no p-value is computed:

if nargout < 2
    
   display = repmat({'off'},numTests,1); 
    
end

% Remove rows of Y with missing values:

Y(any(isnan(Y),2),:) = [];
[numObs,numDims] = size(Y);

% Check if any lags are too big:

bigLags = find(numObs < lags+2);
if ~isempty(bigLags)

    error(message('econ:jcitest:ZeroEffectiveSampleSize', bigLags))

end

% Set row/column grid for tables of asymptotic critical values:

sigLevels = [0.001 (0.005:0.005:0.10) (0.125:0.025:0.20) ...
                   (0.80:0.025:0.875) (0.90:0.005:0.995) 0.999];
minAlpha = sigLevels(1);
maxAlpha = sigLevels(end);

maxDims = 12; % Table rows 1, ..., maxDims (null ranks maxDims-1, ..., 0)

% Check if any alpha is outside of the tables:

if any(alpha < minAlpha) || any(alpha > maxAlpha)

    error(message('econ:jcitest:AlphaOutOfRange', sprintf( '%5.3f', minAlpha ), sprintf( '%5.3f', maxAlpha )))
      
end

% Check if the number of dimensions is too large:

if numDims > maxDims
    
	error(message('econ:jcitest:numDimsExceedsTabulatedValues', maxDims))
      
end

% Preallocate output variables:

switch nargout
    
    case {0,1}
        
        needPValue = false;
        needMLEs = false;
        
        h = false(numTests,numDims);
        
    case 2
        
        needPValue = true;
        needMLEs = false;
        
        h = false(numTests,numDims);
        pValue = NaN(numTests,numDims);
        
    case 3
        
        needPValue = true;
        needMLEs = false;
        
        h = false(numTests,numDims);
        pValue = NaN(numTests,numDims);
        stat = pValue;
        
    case 4
        
        needPValue = true;
        needMLEs = false;
        
        h = false(numTests,numDims);
        pValue = NaN(numTests,numDims);
        stat = pValue;
        cValue = pValue;
        
    case 5
        
        needPValue = true;
        needMLEs = true;
        
        h = false(numTests,numDims);
        pValue = NaN(numTests,numDims);
        stat = pValue;
        cValue = pValue;
        MLEFields = {'paramNames','paramVals','res','EstCov',...
                     'eigVal','eigVec','rLL','uLL'};
        mles0 = cell2struct(cell(size(MLEFields)),MLEFields,2);
        mles = repmat(mles0,numTests,numDims);

end

% Load critical values:

load Data_JCITest JCV

% Initialize loop variables:

lastModel = [];
lastLags = [];
lastType = [];

% Run the tests:

for testNum = 1:numTests
    
    testModel = model{testNum};
    testLags = lags(testNum);
    testType = test{testNum};
    testAlpha = alpha(testNum);
    testDisplay = display{testNum};
    
    tBase = (testLags+2):numObs; % Commensurate time base, all lags
    testT = length(tBase); % Effective sample size
    
    % Check to see if statistics are unchanged from last test:
    
    sameStat = isequal(testModel,lastModel) && ...
               isequal(testLags,lastLags) && ...
               isequal(testType,lastType);
    
    if ~sameStat % Recompute the statistic only if it changes
        
        % Concentrate the Bi and D parameters out of the VEC(q) model to
        % form the reduced-rank regression R0 = A*B'*R1 + ehat. Return LY,
        % DY, and DLY for use in computing MLEs:
        
        [R0,R1,LY,DY,DLY] = conVEC;
            
        % Compute the product moment matrices of the residuals:
        
        S00 = (R0'*R0)/testT;
        S01 = (R0'*R1)/testT;
        S10 = (R1'*R0)/testT;
        S11 = (R1'*R1)/testT;
        
        % Estimate the eigenvectors and eigenvalues of C = A*B':
        
        [V,D] = eig(S10*(S00\S01),S11);
        [lambda,sortIdx] = sort(diag(D),'descend'); % Ordered eigenvalues
        V = V(:,sortIdx); % Corresponding eigenvectors
                
        % Normalize the eigenvectors so that V'*S11*V = I, as in [2]:
        
        U = V'*S11*V; % V diagonalizes S11
        V = bsxfun(@rdivide,V,sqrt(diag(U))');
        
        % Get the appropriate table of critical values:
    
        CVTable = getCVTable;
    
        % Compute the statistic:
        
        [testStat,logLambda] = getStat;
        
        % Compute the p-value, if requested:
        
        if needPValue
            
            testPValue = getPValue;
            
        end
               
    end
    
    % Test the statistic:
    
    [testCValue,testH] = runTest;
    
    % Compute mles, if requested:
    
    if needMLEs
        
        testMLEs = getMLEs;
        
    end

    % Add the test results to the outputs:
    
    switch nargout
    
        case {0,1}
        
            h(testNum,:) = testH;
        
        case 2
        
            h(testNum,:) = testH;
            pValue(testNum,:) = testPValue;
        
        case 3
        
            h(testNum,:) = testH;
            pValue(testNum,:) = testPValue;
            stat(testNum,:) = testStat;
        
        case 4
        
            h(testNum,:) = testH;
            pValue(testNum,:) = testPValue;
            stat(testNum,:) = testStat;
            cValue(testNum,:) = testCValue;
        
        case 5
        
            h(testNum,:) = testH;
            pValue(testNum,:) = testPValue;
            stat(testNum,:) = testStat;
            cValue(testNum,:) = testCValue;
            mles(testNum,:) = testMLEs;
            
    end
    
    % Save values to check for changes, next loop:
    
    lastModel = testModel;
    lastLags = testLags;
    lastType = testType;
    
    % Display results in the command window:
    
    if ~strcmp(testDisplay,'off')
        displayResults
    end

end

% Convert outputs to dataset arrays with rows labeled by test numbers t1,
% t2, etc. and columns labeled by null rank r0, r1, etc.:

rankNames = strcat({'r'},num2str((0:numDims-1)','%-d'));
testNames = strcat({'t'},num2str((1:numTests)','%-d'));

switch nargout

    case {0,1}

        h = dataset([h,rankNames(:)'],'ObsNames',testNames);

    case 2

        h = dataset([h,rankNames(:)'],'ObsNames',testNames);
        pValue = dataset([pValue,rankNames(:)'],'ObsNames',testNames);

    case 3

        h = dataset([h,rankNames(:)'],'ObsNames',testNames);
        pValue = dataset([pValue,rankNames(:)'],'ObsNames',testNames);
        stat = dataset([stat,rankNames(:)'],'ObsNames',testNames);

    case 4

        h = dataset([h,rankNames(:)'],'ObsNames',testNames);
        pValue = dataset([pValue,rankNames(:)'],'ObsNames',testNames);
        stat = dataset([stat,rankNames(:)'],'ObsNames',testNames);
        cValue = dataset([cValue,rankNames(:)'],'ObsNames',testNames);

    case 5

        h = dataset([h,rankNames(:)'],'ObsNames',testNames);
        pValue = dataset([pValue,rankNames(:)'],'ObsNames',testNames);
        stat = dataset([stat,rankNames(:)'],'ObsNames',testNames);
        cValue = dataset([cValue,rankNames(:)'],'ObsNames',testNames);
        mles = dataset({mles,rankNames{:}},'ObsNames',testNames); %#ok

end

%-------------------------------------------------------------------------
% Check input Y
function OK = YCheck(Y)
            
    if isempty(Y)
        
        error(message('econ:jcitest:DataUnspecified'))
          
    elseif isvector(Y)
        
        error(message('econ:jcitest:DataIsVector'))
          
                    
    elseif ~isnumeric(Y)
        
        error(message('econ:jcitest:DataNonNumeric'))
          
    else
        
        OK = true;
        
    end
    
end % YCHECK
    
%-------------------------------------------------------------------------
% Check value of 'model' parameter
function OK = modelCheck(model)
    
    if ~isvector(model)
        
        error(message('econ:jcitest:ModelNonVector'))
          
    elseif isnumeric(model) || (iscell(model) && any(cellfun(@isnumeric,model)))
        
        error(message('econ:jcitest:ModelNumeric'))
    
    elseif ~all(ismember(upper(model),{'H2','H1*','H1','H*','H'}))
        
        error(message('econ:jcitest:ModelInvalid'))
          
    else
        
        OK = true;
        
    end
    
 end % MODELCHECK

%-------------------------------------------------------------------------
% Check value of 'lags' parameter
function OK = lagsCheck(lags)
    
    if ~isvector(lags)
        
        error(message('econ:jcitest:LagsNonVector'))
          
    elseif ~isnumeric(lags)
        
        error(message('econ:jcitest:LagsNonNumeric'))
          
    elseif any(mod(lags,1) ~= 0) || any(lags < 0)
        
        error(message('econ:jcitest:LagsOutOfRange'))
          
    else
        
        OK = true;
        
    end
    
end % LAGSCHECK

%-------------------------------------------------------------------------
% Check value of 'test' parameter
function OK = testCheck(test)
    
    if ~isvector(test)
        
        error(message('econ:jcitest:TypeOfTestNonVector'))
          
    elseif isnumeric(test) || (iscell(test) && any(cellfun(@isnumeric,test)))
        
        error(message('econ:jcitest:TypeOfTestNumeric'))
    
    elseif ~all(ismember(lower(test),{'trace','maxeig'}))
        
        error(message('econ:jcitest:TypeOfTestInvalid'))
          
    else
        
        OK = true;
        
    end
    
end % TESTCHECK

%-------------------------------------------------------------------------
% Check value of 'alpha' parameter
function OK = alphaCheck(alpha)
    
    if ~isvector(alpha)
        
        error(message('econ:jcitest:AlphaNonVector'))
          
    elseif ~isnumeric(alpha)
        
        error(message('econ:jcitest:AlphaNonNumeric'))
          
    else
        
        OK = true;
        
    end
    
end % ALPHACHECK

%-------------------------------------------------------------------------
% Check value of 'display' parameter
function OK = displayCheck(display)
    
    if ~isvector(display)
        
        error(message('econ:jcitest:DisplayFlagNonVector'))
          
    elseif isnumeric(display) || (iscell(display) && any(cellfun(@isnumeric,display)))
        
        error(message('econ:jcitest:DisplayFlagNumeric'))
    
    elseif ~all(ismember(lower(display),{'off','summary','params','full'}))
        
        error(message('econ:jcitest:DisplayFlagtInvalid'))
          
    else
        
        OK = true;
        
    end
    
end % DISPLAYCHECK

%-------------------------------------------------------------------------
% Check parameter values for commensurate lengths, expand single-element
% values, and convert all variables to columns
function [numTests,varargout] = sizeCheck(varargin)

% Initialize outputs:

numTests = 1;

% Determine vector lengths, number of tests:

for i = 1:nargin
    
    ivar = varargin{i};
    iname = inputname(i);
    
    if isnumeric(ivar) || iscell(ivar)
        paramLength.(iname) = length(ivar);  
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
        
            error(message('econ:jcitest:ParameterSizeMismatch', iname, jname))
              
        end        
    end
end

% Expand scalars and single strings:

for i = 1:nargin
    
    ivar = varargin{i};
    if paramLength.(inputname(i)) == 1
        varargout{i} = ivar(ones(numTests,1)); % Expand to column output
    else
        varargout{i} = ivar(:);  % Column output
    end
    
end

end % SIZECHECK

%-------------------------------------------------------------------------
% Concentrate the VEC model
function [R0,R1,LY,DY,DLY] = conVEC
    
YLags = lagmatrix(Y,0:(testLags+1)); % Y(t-k) on observed time base
LY = YLags(tBase,(numDims+1):2*numDims); % Y(t-1) on commensurate time base

% Form multidimensional differences so that the kth numDims-wide block of
% columns in DelatYLags contains (1-L)Y(t-k+1):

DeltaYLags = zeros(testT,(testLags+1)*numDims);
for k = 1:(testLags+1)
    DeltaYLags(:,((k-1)*numDims+1):k*numDims) = ...
               YLags(tBase,((k-1)*numDims+1):k*numDims) ...
             - YLags(tBase,(k*numDims+1):(k+1)*numDims);
end

DY = DeltaYLags(:,1:numDims); % (1-L)Y(t)
DLY = DeltaYLags(:,(numDims+1):end); % [(1-L)Y(t-1), ..., (1-L)Y(t-q)]

% Express the VEC(q) model in the form Z0 = A*B'*Z1 + PSI*Z2 + e(t):

switch testModel

    case 'H2'
        
        Z0 = DY;
        Z1 = LY;
        Z2 = DLY;
        
    case 'H1*'
        
        Z0 = DY;
        Z1 = [LY,ones(testT,1)]; % Restricted constant
        Z2 = DLY;
        
    case 'H1'
        
        Z0 = DY;
        Z1 = LY;
        Z2 = [DLY,ones(testT,1)]; % Unrestricted constant
        
    case 'H*'
        
        Z0 = DY;
        Z1 = [LY,(1:testT)']; % Restricted linear term
        Z2 = [DLY,ones(testT,1)]; % Unrestricted constant
        
    case 'H'
        
        Z0 = DY;
        Z1 = LY;
        Z2 = [DLY,ones(testT,1),(1:testT)']; % Unrestricted constant, linear term
        
end

% Perform the concentrating regressions:

R0 = Z0-Z2*(Z2\Z0); % Z0 on Z2 residuals
R1 = Z1-Z2*(Z2\Z1); % Z1 on Z2 residuals
    
end % CONVEC

%-------------------------------------------------------------------------
% Get table of critical values
function CVTable = getCVTable

switch testModel

    case 'H2'

        switch testType

            case 'trace'

                CVTable = JCV(:,:,1,1);

            case 'maxeig'

                CVTable = JCV(:,:,1,2);

        end

    case 'H1*'

        switch testType

            case 'trace'

                CVTable = JCV(:,:,2,1);

            case 'maxeig'

                CVTable = JCV(:,:,2,2);

        end

    case 'H1'

        switch testType

            case 'trace'

                CVTable = JCV(:,:,3,1);

            case 'maxeig'

                CVTable = JCV(:,:,3,2);

        end

    case 'H*'

        switch testType

            case 'trace'

                CVTable = JCV(:,:,4,1);

            case 'maxeig'

                CVTable = JCV(:,:,4,2);

        end
        
    case 'H'
        
        switch testType

            case 'trace'

                CVTable = JCV(:,:,5,1);

            case 'maxeig'

                CVTable = JCV(:,:,5,2);

        end        

end

end % GETCVTABLE

%-------------------------------------------------------------------------
% Compute the statistic
function [testStat,logLambda] = getStat

% Form summands:

logLambda = log(abs(1-lambda(1:numDims))); % ABS for rounding lambda ~ 1

% Compute statistics for cointegration ranks r = 0, ..., numDims-1:

switch testType
    
    case 'trace' % -T*[log(1-lambda(r+1)) + ... + log(1-lambda(numDims))]
        
        testStat = -testT*flipud(cumsum(flipud(logLambda)));
        
    case 'maxeig' % -T*log(1-lambda(r+1))
        
        testStat = -testT*logLambda;
        
end

end % GETSTAT

%-------------------------------------------------------------------------
% Get p-values
function testPValue = getPValue

testPValue = NaN(numDims,1);

for d = numDims:-1:1

    CVTableRow = CVTable(d,:); % Row for rank difference d = numDims-r
    idx = numDims-d+1; % Index of corresponding testStat, pValue

    % Right-tailed test
    
        if testStat(idx) >= CVTableRow(1)
        
            testPValue(idx) = sigLevels(1);
           % warning(message('econ:jcitest:RightTailStatTooBig', idx, sprintf( '%5.3f', testPValue( idx ) )))
            
        elseif testStat(idx) <= CVTableRow(end)
        
            testPValue(idx) = sigLevels(end);
            warning(message('econ:jcitest:RightTailStatTooSmall', idx, sprintf( '%5.3f', testPValue( idx ) )))
            
        else
        
            testPValue(idx) = interp1(CVTableRow,sigLevels,testStat(idx),'linear');
            
        end
    
end

end % GETPVALUE

%-------------------------------------------------------------------------
% Test the statistic
function [testCValue,testH] = runTest
    
testCValue = NaN(numDims,1);

for d = numDims:-1:1

    CVTableRow = CVTable(d,:); % Row for rank difference d = numDims-r
    idx = numDims-d+1; % Index of corresponding testStat, pValue
    
    testCValue(idx) = interp1(sigLevels,CVTableRow,testAlpha,'linear');
    
end

testH = (testStat > testCValue); % Right-tailed test

end % RUNTEST

%-------------------------------------------------------------------------
% Compute MLEs of the VEC model
function testMLEs = getMLEs
    
testMLEs = repmat(mles0,numDims,1);
    
for r = 0:(numDims-1)
    
    idx = r+1;
    
    A = S01*V(:,1:r);
    B = V(1:numDims,1:r);   
    
    switch testModel
        
        case 'H2'
            
            detNames = cell(0,1);
            W = DY-LY*B*A';
            X = DLY;
            P = (X\W)'; % [B1, ..., Bq]
            detVals = cell(0,1);
            [paramNames,paramVals] = assembleParams;
            
        case 'H1*'
            
            detNames = {'c0'};
            c0 = V(end,1:r)';
            W = DY-(LY*B+repmat(c0',testT,1))*A';
            X = DLY;
            P = (X\W)'; % [B1, ..., Bq]
            detVals = {c0};
            [paramNames,paramVals] = assembleParams;
            
        case 'H1'
            
            detNames = {'c0';'c1'};
            W = DY-LY*B*A';
            X = [DLY,ones(testT,1)];
            P = (X\W)'; % [B1, ..., Bq, c]
            c = P(:,end);
            c0 = A\c; % Project c on A
            c1 = c-A*c0; % Orthogonal component of c
            detVals = {c0;c1};
            [paramNames,paramVals] = assembleParams;
            
        case 'H*'
            
            detNames = {'c0';'d0';'c1'};
            d0 = V(end,1:r)';
            W = DY-(LY*B+(1:testT)'*d0')*A';
            X = [DLY,ones(testT,1)];
            P = (X\W)'; % [B1, ..., Bq, c]
            c = P(:,end);
            c0 = A\c; % Project c on A
            c1 = c-A*c0; % Orthogonal component of c
            detVals = {c0;d0;c1};
            [paramNames,paramVals] = assembleParams;
            
        case 'H'
            
            detNames = {'c0';'d0';'c1';'d1'};
            W = DY-LY*B*A';
            X = [DLY,ones(testT,1),(1:testT)'];
            P = (X\W)'; % [B1, ..., Bq, c, d]
            c = P(:,end-1);
            c0 = A\c; % Project c on A
            c1 = c-A*c0; % Orthogonal component of c
            d = P(:,end);
            d0 = A\d; % Project d on A
            d1 = d-A*d0; % Orthogonal component of d
            detVals = {c0;d0;c1;d1};
            [paramNames,paramVals] = assembleParams;
            
    end
    
    testMLEs(idx).paramNames = paramNames;    
    testMLEs(idx).paramVals = paramVals;
    testMLEs(idx).res = W-X*P';
    testMLEs(idx).EstCov = S00-A*A';
    testMLEs(idx).eigVal = lambda(idx);
    testMLEs(idx).eigVec = V(:,idx);
    testMLEs(idx).rLL = getRLL;
    testMLEs(idx).uLL = getULL;
            
end

    %----------------------------------------------------------------------
    % Assemble parameter names/values
    function [paramNames,paramVals] = assembleParams

    if testLags > 0
        lagNames = strcat({'B'},num2str((1:(testLags))','%-d'));
    else
        lagNames = cell(0,1);
    end

    paramNames = cat(1,'A','B',lagNames,detNames);

    lagVals = cell(testLags,1);
    for lag = 1:testLags        
        lagVals{lag} = P(:,(((lag-1)*numDims)+1):lag*numDims);        
    end

    paramVals0 = cat(1,{A},{B},lagVals,detVals);
    
    % Convert to structure:
    
    for param = 1:length(paramNames)
        paramVals.(paramNames{param}) = paramVals0{param};
    end

    end % ASSEMBLEPARAMS
    
    %----------------------------------------------------------------------
    % Compute restricted loglikelihood under null H(r)
    function rLL = getRLL
        
        rLL = -(testT*numDims/2)*(log(2*pi)+1)...
              -(testT/2)*(log(det(S00))+sum(logLambda(1:r)));

    end % GETRLL

    %----------------------------------------------------------------------
    % Compute unrestricted loglikelihood
    function uLL = getULL
        
        switch testType
            
            case 'trace' % Under alternative H(numDims)
                
               uLL = -(testT*numDims/2)*(log(2*pi)+1)...
                     -(testT/2)*(log(det(S00))+sum(logLambda)); 
                
            case 'maxeig' % Under alternative H(r+1)
                
               uLL = -(testT*numDims/2)*(log(2*pi)+1)...
                     -(testT/2)*(log(det(S00))+sum(logLambda(1:(r+1)))); 
                
        end

    end % GETULL

end % GETMLES

%--------------------------------------------------------------------------
% Display results in the command window
function displayResults
    
    if ~strcmp(testDisplay,'params')
        
    %fprintf('\n************************')
    %fprintf('\nResults Summary (Test %d)\n',testNum)

    %fprintf('\nData: %s',dataName)
    %fprintf('\nEffective sample size: %d',testT)
    %fprintf('\nModel: %s',testModel)
    %fprintf('\nLags: %d',testLags)
    %fprintf('\nStatistic: %s',testType)
    %fprintf('\nSignificance level: %3.2f\n\n',testAlpha)

    %fprintf('\n%-3s%-3s%-10s%-9s%-9s%-9s',...
           % 'r','h','stat','cValue','pValue','eigVal')
    %fprintf('\n----------------------------------------')
    
    for r = 0:numDims-1

        idx = r+1;

        %fprintf('\n%-3d%-3d%-9.4f% -9.4f% -9.4f% -9.4f',...
                %r,testH(idx),testStat(idx),testCValue(idx),testPValue(idx),lambda(idx))

    end

    %fprintf('\n')
    
    end
    
    if any(strcmp(testDisplay,{'params','full'}))
        
        %fprintf('\n****************************')
        %fprintf('\nParameter Estimates (Test %d)\n',testNum)
        
        for r = 0:numDims-1
            
            idx = r+1;
            paramNames = testMLEs(idx).paramNames;
            paramVals = testMLEs(idx).paramVals;
            %fprintf('\nr = %d',r)
            %fprintf('\n------\n')
            
            if r == 0 % Remove empty parameters from display
                paramNames(ismember(paramNames,{'A','B','c0','d0'})) = [];
            end
            
            for param = 1:length(paramNames)
                
                %fprintf('%s =\n',paramNames{param})
                disp(paramVals.(paramNames{param}))
                
            end
            
        end
        
    end

end % DISPLAYRESULTS

%-------------------------------------------------------------------------

end % JCITEST