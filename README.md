# StochKit2R

## Notes About Changes

- `CXX_STD = CXX17` in `Makevars` will likely require users to be running R4.5+
  - It might not be required to have StochKit2R work but it's a good safety
    measure either way.

## Development

Rebuild c++ Exported Functions for Rcpp:

```r
Rcpp::compileAttributes()
```

Rebuild Documentation & Exported Functions:

```r
devtools::document()
```

Rebuild Vignettes:

```r
devtools::build_vignettes()
# or to see all of stdout/stderr
devtools::build_vignettes(quiet = FALSE)
```

Install Package from Local Directory:

```r
install.packages("/path/", repos = NULL, type = "source")
```
