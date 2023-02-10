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

(largely copied from [naomi](https://github.com/mrc-ide/naomi))
- Make changes in a new branch,
- Run checks (requires `goodpractice` package). Also documents functions with
  `devtools::check()`. These checks can be run directly from the script
  `scripts/goodpractice.R` or from the command line with `make goodpractice`
  (use `make goodpractice_vignette` if making changes to vignettes, which are
  very slow to build),
- When branch is ready for merging create a PR and add a reviewer,
- Ensure that the version number has been updated according to 
  [semantic versioning](https://semver.org/) and add a news item describing the 
  change, and
- Reviewer should check code and ensure the build passes on Buildkite before 
  merging.

Additionally: 
- Wherever possible, add arguments for optional new functionality, with old 
  functionality as default,
- If changing `threemc_aggregate` or `threemc_ppc` (the slowest "purely R" 
  functions in this package), a before vs after profiling would be useful,
  using e.g. `profvis`,
- Follow the [tidverse styleguide](https://style.tidyverse.org) as much as
  possible, but at least have your commits pass the call to `lintr` from 
  `goodpractice`,
- Ensure that git commits (and R comments) are in the imperative case (see 
  [here](https://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html)
  for Tim Pope's rational behind this),
- To avoid notes from `devtools::check()` when using non-standard evaluation (NSE), 
  - in `dplyr`, refer to all columns with ``` rlang::`.data` ```, e.g. `.data$col`,
  - in `data.table`, assign columns referenced in NSE to the value NULL at the 
  beginning of a function (see `threemc_ppc` for an example). 

## License

MIT © Imperial College of Science, Technology and Medicine
