# Changelog

## frpp 0.0.0.9000

- Add [`fsim()`](https://ksublee.github.io/frpp/reference/operations.md)
  as the primary simulation function. Existing
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) methods call it,
  preserving path and RNG behavior and old serialized classes.

- Assemble numeric state-backend likelihoods using the experimental
  public `emhawkes::exp_kernel_replay()` interface, without Hawkes
  likelihood evaluation or full component histories when inference is
  not requested. This requires the emhawkes development interface; the
  bridge remains available with 1.0.0.

- Use [`fspec()`](https://ksublee.github.io/frpp/reference/frpp_spec.md)
  as the primary composition constructor. Keep
  [`frpp_spec()`](https://ksublee.github.io/frpp/reference/frpp_spec.md)
  and old serialized classes compatible; add
  [`as_fspec()`](https://ksublee.github.io/frpp/reference/frpp_spec.md).

- Add [`ffit()`](https://ksublee.github.io/frpp/reference/ffit.md) with
  maxLik and stats optimizers, optional joint residual estimation,
  shared/fixed parameters, constraints and a domain validator.

- Add fitted-model restoration, initialization replay, coefficient,
  likelihood, covariance, summary and diagnostic methods. Residual
  parameters stay fixed by default; invalid Wald covariance is withheld
  rather than repaired.
