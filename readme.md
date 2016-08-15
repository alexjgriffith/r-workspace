# R Workspace

## Development work done for CCCA

## Alexander Griffith MASc, UOttawa

### Version 0.1

* Functions are labled *month*-function.r
* Variables are labled *month*-variables.r

### Notes
Specific tasks should have their own names. e.g. those that provide examples to gene assosiation or motif denovo. The examples should be self contained and when used transcribed into the specific *month*-*day*.r file appropriate for the day. There should be a refference included to any examples and as functions become codified add Roxygen documentation and integrate them into the CCCA package.


Following the reorganize commit things may be broken.

To fix switch:

```R
source("~/r-workspace/nov-functions.r")
source("~/r-workspace/nov-variables.r")
```

To:

```R
source("~/r-workspace/functions/nov-functions.r")
source("~/r-workspace/variables/nov-variables.r")
```
