## Test environments
* local Manjaro Linux 24.1.2 install, R 4.4.1
* win-builder (release, oldrelease, devel)
* mac-builder
* R-hub (configurations: linux, macos, macos-arm64, windows, valgrind)


## R CMD check results

### local
0 errors | 0 warnings | 0 notes

### win-builder
0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Florian Junge <diso.fbmn@h-da.de>'

New submission

Possibly misspelled words in DESCRIPTION:
  Bonferroni (15:52)
  FWER (3:8)
  Guo (18:73)
  Hochberg (15:70)
  Holm (15:64)
  Zhu (18:67)
  Döhler (17:15)

This is a new package. The "misspelled" words are names and an acronym.

### mac-builder
0 errors | 0 warnings | 0 notes

### R-hub
0 errors | 0 warnings | 0 notes
