# Introduction

This package implements the analyses from Bilinski and Hatfield (2025) Statistics in Medicine "Nothing to See Here: A non-inferiority approach to parallel trends. 
The key function is `run_NI_test`, which fits reduced and expanded models and tests for the difference between linear combinations of their parameters.
To reproduce the results shown in the manuscript for the analysis of the dependent coverage provision on any health insurance, use the following code:

```{r run_models}
library(NonInfParTren)
data("depcov")
base.model <- as.formula(paste("anyhi","~",paste(c("fipstate","factor(trend)",
                                                      "enact.trt.month","impl.trt.month", # these are treat.by.post
                                                      "female","hispanic",
                                                      "white", "asian", "other",
                                                      "mar", "student","fpl_ratio","fpl_ratio_2",
                                                      paste0("age",c(17:25,27:29))),collapse="+")))

grp.lin.yr <- run_NI_test(data=dat.for.reg,reduced = base.model,
                           # Add treatment-group specific linear trends
                           expanded = update.formula(base.model,"~ . + fedelig:trend"),
                           lincom_var = 'impl.trt',cluster = 'fipstate',weight = 'weight')
```

