# GenelandHandleR
An R package to support analyses with Geneland

At this stage, this package has functions to run Geneland in parallel (this should work on Windows, unix and Mac machines), 
run MCMC diagnostics and plot results (bar plots and on a map).

## Quick start
Install the package from version control from within R:
```
library(devtools)
install_github("carlopacioni/GenelandHandleR")
```

If you are on Windows and have not used devtools before, then you have to download the Rtools executable file from CRAN webpage and run it. devtools can be installed from within R with

Then you can use the following to install it:
```
install.packages("devtools")
```

If you want to access the package tutorials, you can use the option:
```
install_github("carlopacioni/HexSimR",  build_vignettes=TRUE)
```
but this requires some setting up. 
