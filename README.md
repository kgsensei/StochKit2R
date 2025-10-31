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

## Todo

- Move `inst/dimer_decay.xml` and `inst/schlogl.xml` into a `data/` directory
  and update documentation and `DESCRIPTION` accordingly. This seems to be a
  modern standard for packages that come with datasets.
- Verify current solution works on Linux (who cares about Mac anyway) and other
  potentially different environments.
- Attempting to browse vignettes currently displays douplicate vignettes for
  StochKit2R, I'll have to look into this later and see if it's a me-only
  problem.
