## Resubmission

This is a resubmission. In response to the previous feedback, we have made the following changes:

* Expanded the Description field to a one-paragraph overview of the package's functionality
* Added \value sections to the Rd files for the exported functions, as follows:
  - eBa0.Rd: Describes return as no return value, called for side effects.
  - elbo_fixed_diagonal.Rd: Describes return as no return value, called for side effects.
  - generate_log_prob.Rd: Describes return as no return value, called for side effects.
  - plot.CVIoutput.Rd: Describes return as a ggplot object representing visualisation
* In R/cvi_singlerun.R, replaced cat() with message().

## R CMD check results

0 errors | 0 warnings | 0 note
