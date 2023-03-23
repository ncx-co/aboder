# aboder <img src="man/figures/logo.png" align="right" height="139" />

*Approximate Bayesian Optimal Design of Experiments in R*

This package aids in Bayesian Optimal Experimental Design (BOED) utilizing the `brms` package as the base modeling framework.

The main purpose is to estimate the expected information gain (EIG) for a potential new data point (design) to collect.
This can be used to compare multiple designs to determine which design will be most useful to collect next.

All that is needed is a fitted `brms` model and a design data frame.
The fitted model can be as simple as a fitted prior, no data is needed to start with.
The design data frame is composed of independent variables that will mirror what will be used to fit your model.

## Example Usage

```r
dummy_data <- data.frame(x = 0, y = 0, weight = 0)
model <- brms::brm(
  y | weights(weight) ~ x, # use weights to fit a prior with no/dummy data
  data = dummy_data,
  sample_prior = "only",
  prior = brms::prior(normal(0, 1), class = "b")
)

designs <- data.frame(x = 1:10)
designs$eig <- aboder::batch_nmc_eig(designs, model)
designs
```

## References

- [Variational Bayesian Optimal Experimental Design](https://arxiv.org/pdf/1903.05480.pdf)
- [Simulation-based optimal Bayesian experimental design for nonlinear systems](https://arxiv.org/pdf/1108.4146.pdf)
- [On Nesting Monte Carlo Estimators](https://arxiv.org/pdf/1709.06181.pdf)
- [Bayesian Experimental Design: A Review](https://projecteuclid.org/journals/statistical-science/volume-10/issue-3/Bayesian-Experimental-Design-A-Review/10.1214/ss/1177009939.full)
- [Pyro: Designing Adaptive Experiments to Study Working Memory](http://pyro.ai/examples/working_memory.html)
