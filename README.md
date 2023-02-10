# threemc

<!-- badges: start -->
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![R build status](https://github.com/mrc-ide/threemc/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/threemc/actions)
<!-- badges: end -->

## Installation

To install `threemc`:

```r
remotes::install_github("mrc-ide/threemc", upgrade = FALSE)
```
## threemc ##

To see how `threemc` works for a simple example, please see the relevant vignette.

## Development steps ##

- Make changes in a new branch,
- Run checks (requires `goodpractice` package). Also documents functions with
  `devtools::check()`. These checks can be run directly from the script
  `scripts/goodpractice.R` or from the command line with `make goodpractice`
  (use `make goodpractice_vignette` if making changes to vignettes, which are
  very slow to build),
- When branch is ready for merging create a PR and add a reviewer,
- Ensure that the version number has been updated according to semantic 
  versioning and add a news item describing the change, and
- Reviewer should check code and ensure the build passes on Buildkite before 
  merging.

## License

MIT © Imperial College of Science, Technology and Medicine
