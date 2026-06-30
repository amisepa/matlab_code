function fm = my_specparam_fit_wrapper(freqs, psd, freq_range, ~)
% SPECPARAM_FIT  Thin MATLAB wrapper for Python specparam.
% 4th argument (settings) accepted but ignored — specparam defaults used.
% 
% % Test several known exponents
% for true_exp = [0.5, 1.0, 1.5, 2.0]
%     psd = 1 ./ (1:30).^true_exp;
%     fm  = my_specparam_fit_wrapper(1:30, psd, [1 30], struct());
%     fprintf('true=%.1f  fit=%.3f  R2=%.4f\n', true_exp, fm.aperiodic_params(end), fm.r_squared);
% end
% % --> Fitted exponents should match exactly (R²≈1) since there are no peaks to confuse the fit.
%
% % Compare with python specparam:
% sp = py.specparam.SpectralModel(pyargs('verbose', false));
% sp.fit(py.numpy.array(1:30), py.numpy.array(1./(1:30)), py.list({1, 30}));
% fm = my_specparam_fit_wrapper(1:30, 1./(1:30), [1 30], struct());
% py_exp = double(sp.results.params.aperiodic.params);
% fprintf('Python exp: %.4f\n', py_exp(end));
% fprintf('Wrapper exp: %.4f\n', fm.aperiodic_params(end));
% 
% 
% Cedric Cannard, March 2026


sp = py.specparam.SpectralModel(pyargs('verbose', false));
sp.fit(py.numpy.array(freqs), py.numpy.array(psd), ...
       py.list({freq_range(1), freq_range(2)}));

ap      = double(sp.results.params.aperiodic.params);
metrics = sp.results.metrics.results;

fm.aperiodic_params = ap;
fm.r_squared        = double(metrics{'gof_rsquared'});
fm.error            = double(metrics{'error_mae'});
fm.freqs            = freqs;

end

