# `tuts`: time-uncertain time series in R

`tuts` is an R package that enables time-uncertain time series analysis in R using an errors-in-variables approach. A stable version is not yet on CRAN but the GitHub version can be installed with:

```
library(devtools) # Install if not available
install_github('peterfranke/tuts')
library(tuts)
```

You then might like to look at the examples in `help(tubfs)` for time uncertain Bayesian frequency selection, or `help(tuar1)` for time uncertain Bayesian autoregressive modelling (currently only order 1), amongst others. A more complete vignette will be available shortly.

If you find `tuts` useful please cite it as:
Franke, P. M., Huntley, B., & Parnell, A. C. (2018). Frequency selection in paleoclimate time series: A model‐based approach incorporating possible time uncertainty. Environmetrics, 29(2), e2492.

Or using bibles:

```
@article{franke2018frequency,
  title={Frequency selection in paleoclimate time series: A model-based approach incorporating possible time uncertainty},
  author={Franke, Peter M and Huntley, Brian and Parnell, Andrew C},
  journal={Environmetrics},
  volume={29},
  number={2},
  pages={e2492},
  year={2018},
  publisher={Wiley Online Library}
}
```