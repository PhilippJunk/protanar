# protanar: PROTemics ANAlysis in R

A R package for the tidy analysis of proteomics data.

## Installation

To install the package without building the vignettes (currently recommended), please run

```R
remotes::install_github(
   'philippjunk/protanar'
)
```

If you want to install it with the vignette, run this
```R
remotes::install_github(
    'philippjunk/protanar',
    build_vignettes=TRUE,
    dependencies=TRUE
)
```

Since the vignette currently depends on the download of several files, which appears to not be working only sometimes, it is recommended to install the package without the vignette at the moment.

## Examples

Check out the vignettes for an example workflow.

## License

protanar is freely available under the [MIT](https://choosealicense.com/licenses/mit/) license.
