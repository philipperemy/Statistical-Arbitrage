function [  ] = SpreadProcessor( spread )
   spr = spread.px;
   mu = mean(spr);
   sigma = std(spr);
   spr_norm = (spr - mu)/sigma;
   f = ReversionFrequencyCalculator(spr_norm, 50);
   mse_m = mse(spr);
   fprintf('mean = %f, std = %f, rev = %f, mse = %f\n', mu, sigma, 1/f, mse_m);
   %plot(spr_norm);
end

