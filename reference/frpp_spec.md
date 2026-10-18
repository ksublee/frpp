# Compose a Kernel Specification and a Residual Law

Create an experimental FRPP object containing an emhawkes kernel and a
separate residual-law object. The result does not inherit from `hspec`.

## Usage

``` r
fspec(kernel, residual)
frpp_spec(kernel, residual)
as_fspec(object, name = "legacy")
as_frpp(object, name = "legacy")
as_legacy_hspec(object)
```

## Arguments

- kernel:

  An `hspec` without residual-distribution slots.

- residual:

  A law created by `residual_law`.

- object:

  For conversion, an existing FRPP `hspec` or a `frpp_spec`. For the
  legacy adapter, a `frpp_spec`.

- name:

  A descriptive name for a converted legacy law.

## Details

The current model refreshes every dimension's candidate residual after
every event in the combined process. It does not retain separate
residual clocks. Conversion preserves callback functions and their
default bindings without executing them. Captured environments remain
references, as in emhawkes. The primary constructor is `fspec`, matching
the `hspec` naming convention. `frpp_spec` and `as_frpp` remain
compatibility names. New objects have classes `c("fspec", "frpp_spec")`;
saved objects from the first prototype can still use the old simulation
and likelihood methods. The adapter is transitional: it restores the old
FRPP representation to use existing public emhawkes operations. It is
not a pure Hawkes specification.

## Value

`frpp_spec` and `as_frpp` return a composition object. `as_legacy_hspec`
returns an `hspec` carrying the residual callbacks.

## See also

[`residual_law`](https://ksublee.github.io/frpp/reference/residual_law.md),
[`simulate.frpp_spec`](https://ksublee.github.io/frpp/reference/operations.md)

## Examples

``` r
kernel <- methods::new("hspec", mu = .3, alpha = .2, beta = 1)
law <- residual_law(c(shape = 2, rate = 2), stats::dgamma,
                    stats::pgamma, stats::qgamma, stats::rgamma,
                    name = "Gamma(2,2)")
model <- fspec(kernel, law)
model
#> FRPP specification (fspec)
#> Event types: 1
#> Residual law: Gamma(2,2)
#> Residual candidates refreshed after every global event.
stopifnot(!methods::is(model, "hspec"))
```
