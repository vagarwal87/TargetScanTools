 
function [estimates, model]=fitcurve(xdata, ydata) 
model = @expfun; 
start_point = [-0.5 1.5 3 2]; 
estimates = fminsearch(model, start_point); 
    function [sse, FittedCurve]=expfun(params) 
        b_0 = params(1); 
        b_1 = params(2); 
        b_2 = params(3); 
        b_3 = params(4); 
        FittedCurve = max(0,b_0 + b_1.* (1 ./ (1 + exp(-b_2.*xdata+ b_3)))); 
        ErrorVector = FittedCurve - ydata; 
        sse = sum(ErrorVector .^ 2) + 10 .*  (FittedCurve(1)+FittedCurve(2)); 
    end 
end 
