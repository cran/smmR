# smmR 1.0.1

## Minor Improvements

* Added a `NEWS.md` file to track changes to the package.

* Remove `Depends` field in DESCRIPTION file.

## Bug fixes

* For initial distribution, transition matrices and conditional sojourn time 
distributions, we allow the fact that the sum of the probabilities may not be 
equal to 1 exactly due to round-off errors (we introduce `.Machine$double.eps`).