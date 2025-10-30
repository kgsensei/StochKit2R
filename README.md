# StochKit2R

## Development

Rebuild c++ exported functions for Rcpp:

```r
Rcpp::compileAttributes()
```

Rebuild vignettes:

```r
devtools::build_vignettes()
# or to see all of stdout/stderr
devtools::build_vignettes(quiet = FALSE)
```
