# PBI model and validation

## Model

1. **Extract a template-informed time series.** Within the extraction mask, regress the fMRI spatial pattern at each time point on a fixed biological map (for example, a PET receptor template) and an intercept. Retain the template coefficient as a function of time. Record the template, support mask, normalization, preprocessing, and time grid. Treat a single-template fit and a fit conditional on other templates as different estimands.
2. **Estimate task response (planned).** Fit the scan-level BOLD component with HRF-convolved condition regressors and confounds. Retain task-related variation during extraction. Return condition estimates and contrasts such as large-reward versus no-reward anticipation.
3. **Estimate task-modulated coupling.** For the deconvolution-based gPPI model, multiply the microtime physiological signal by each condition vector before HRF convolution and sampling at scan times. Fit target BOLD with condition main effects, a scan-level physiological main effect, the condition-specific interaction columns, and confounds. Contrast interaction coefficients to estimate the difference in coupling. Preserve the physiological signal and named design columns so this model can be audited.

## Implementation and validation

- Use the same inputs and named outputs to compare the MATLAB/SPM and Python workflows. Compare template coefficients, sampled design columns, and contrast estimates before comparing statistic maps.
- Document the deconvolution, temporal filtering, and autocorrelation model for each implementation. Specify numerical tolerances for operations intended to match and quantify differences when the algorithms differ.
- Start with synthetic volume data with known task responses and coupling. Include a case with nonzero common coupling and zero condition difference, and a case with a known interaction. Then validate with a shareable task fMRI example.
- Keep component extraction and target geometry separate. Verify alignment of the biological map and fMRI within the extraction mask, and preserve the target grid or vertex mapping in the output. For a volume-defined component and surface targets, the same Bio time series can be used with each target geometry. A cortical-only map projection changes the extraction domain and should be documented accordingly.
