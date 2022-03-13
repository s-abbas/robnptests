---
title: "robnptests -- An R package for robust two-sample location and variability tests"
authors:
- affiliation: 1
  name: Sermad Abbas
  orcid: 0000-0001-9162-9792
- affiliation: 2
  name: Barbara Brune
- affiliation: 1
  name: Roland Fried
date: 09 January 2021
output:
  pdf_document: default
  html_document: default
tags:
- R
- robust statistics
- nonparametric statistics
affiliations:
- index: 1
  name: TU Dortmund University
- index: 2
  name: Technical University of Vienna
bibliography: REFERENCES.bib
csl: ../vignettes/csda.csl
---

# Summary
The package `robnptests` is a compilation of two-sample tests  based on two criteria: The tests are (i) robust against outliers and (ii) (approximately) distribution free. 
Regarding the latter aspect, we implemented tests that keep an intended significance level and provide a reasonably high efficiency, both  under a variety of continuous distributions.
Robustness is achieved by using test statistics that are based on robust location and scale measures. 

In what follows, we give a brief description of the package's contents.
More details can be found in the introductory vignette of the package, which can be called by `vignette("robnptests-vignette")`, and the cited papers.

# Data situation
We consider two samples of independent and identically distributed (i.i.d.) random variables $X_1, ..., X_m$ and $Y_1, ..., Y_n$, respectively. 
The underlying distributions are assumed to be continuous with cumulative distribution functions $F_X$ and $F_Y$.

The tests can be used for either one of the following cases:

1. Assuming that the parameters of both distributions are equal while the location parameters may be unequal, i.e. $F_X(x) = F_Y(x + \Delta)$, $\Delta \in \mathbb{R}$, the tests can be used to detect a location difference between the samples.
2. For equal location parameters with possibly unequal scale parameters, i.e. $F_X(x) = F_Y(x/\theta)$, $\theta > 0$, a transformation of the observations enables to identify differing scale parameters.


# Statement of need
A popular test for the location setting is the two-sample $t$-test.
It is considered to be robust against deviations from the normality assumption by keeping the specified significance level due to the central limit theorem.
However, non-normality can result in a loss of power [@Wil03appl].
In addition, particularly for small samples, the $t$-test is prone to outliers [@FriDeh11robu].
Distribution-free tests, like the two-sample Wilcoxon rank-sum test, can be nearly as powerful as the $t$-test under normality and may have higher power under non-normality.
Still, they also can be vulnerable to outliers [@FriGat07rank].
The two-sample tests in `robnptests` combine (approximate) distribution independence and robustness against outliers.
Thus, they are better suited for outlier-corrupted samples from unknown data-generating distributions.
At the same time, such tests can be nearly as powerful as popular procedures like the aforementioned $t$-test or the Wilcoxon test on uncontaminated samples.

Figure 1 compares the power of the $t$-test, the Wilcoxon test and two robust tests - one based on the one-sample Hodges-Lehmann estimator (HL1 test) [@HodLeh63esti] and one based on Huber's M-estimator [@Hub64robu] - for a fixed location difference between the samples and a single outlier of increasing size.
The power of the $t$-test decreases to zero, while the loss in power of the Wilcoxon test and both robust tests is bounded. 
The robust tests provide a higher power than the Wilcoxon test.
The differences between the Wilcoxon test and robust tests like those in the example can become more apparent when more outliers are involved [@FriDeh11robu].
  
![Power of the two-sample $t$-test, the Wilcoxon rank sum test, and two robust tests - one based on the one-sample Hodges-Lehmann estimator and one based on Huber's M-estimator - on two samples of size $m = n = 10$ from two normal distributions with a location difference of $\Delta = 2$ and a single outlier of increasing size.](img/fig1_-_power_under_outliers.pdf){width=4in height=4in}

Like the $t$-test for the location problem, the performance of parametric tests for the scale problem may deteriorate when their distributional assumption is not fulfilled.
Moreover, nonparametric tests may be robust against violations of the distributional assumptions or outliers, but many of them do not cope well with asymmetry [@Fri12onli].
Applying the robust location tests to transformed observations as proposed by @Fri12onli, yields good results in terms of power and size under asymmetry and outlier corruption.
However, the tests may be less efficient under symmetry. 

# Implemented two-sample tests
Each test statistic consists of a robust estimator for the location difference between the two populations that should be compared.
This difference is divided by a robust estimator for the within-sample variance. 

To obtain a distribution-free test decision, the $p$-value can be computed by using the permutation principle, the randomization principle, or a normal approximation.
With the permutation principle, the tests hold the desired significance level exactly at the cost of large computing times even for quite small samples such as $m = n = 10$.
The time can be reduced by using a randomization distribution and, even more, by taking advantage of the asymptotic normality of the location-difference estimators.
The latter approach, however, is only advisable for large sample sizes $m, n > 30$.

The tests based on the following location estimators are described in @FriDeh11robu:

* The _difference of the sample medians_ helps to achieve high robustness. However, this estimator is not very efficient under the normal distribution or distributions that do not deviate too much from it.
* To improve the efficiency, one can use the difference of the _one-sample Hodges-Lehmann estimators_ [@HodLeh63esti] at the cost of losing robustness due to the lower breakdown point.
* Similarly, the _two-sample Hodges-Lehmann estimator_ leads to a robust test with a higher power under normality than the tests based on the sample median.

For scaling, we use different estimators based on medians and pairwise differences, see @FriDeh11robu for a detailed description.

In addition, we implemented tests based on _M-estimators_. This approach to robust location estimation allows for flexibility in how outliers are treated through the specification of the parameters of the corresponding $\rho$-function. 
We focus on Huber's $\rho$-function, the bisquare function and the Hampel $\rho$-function in a similar manner as described in @Abo92robu.
Moreover, the package contains Yuen's $t$-test which uses the difference of _trimmed means_ to estimate the location difference [@YueDix73appr].

Some of the robust scale estimators may become zero when ties occur in the data.
This can happen in real-world applications when the measurements are rounded or stem from discrete distributions. 
Discretization can also lead to a loss in power or conservative tests.
To cope with this, we add random noise from a uniform distribution with a small variance to each observation. This procedure is called `wobbling`. 
The objective to enable the computation of the test statistic without distorting the observations too much.

A more detailed overview of the implemented tests and corresponding test statistics can be found in the vignettes of the package.

## Applications
Besides conventional two-sample problems, the tests can be used for the online detection of structural breaks in outlier-contaminated time series.
@AbbFriGat16dete describe how intensity changes in image sequences generated by a virus sensor are automatically detected by applying the tests to the individual pixel time series of the sequence.
Moreover, the test statistics can be used as control statistics for robust, (approximately) distribution-free control charts for time series with a time-varying signal [@AbbFri17cont; @AbbFri20robu].
In @AbbFriHei19dete, the tests are applied to detect unusual sequences in time series of crack widths in concrete beams by searching for sudden scale changes.

## Other packages with robust two-sample tests
The CRAN Task View for robust statistical methods currently lists the three packages `WRS2` [@MaiWil19wrs2], `walrus` [@LovMai18walr], and `robeth` [@Mar20robe] that explicitly deal with robust hypothesis tests for the two-sample problem.
`WRS2` contains a large collection of robust procedures which are presented in the book *Introduction to Robust Estimation and Hypothesis Testing* by [@Wil17intr]. `walrus` provides a different user interface for the functions in `WRS2`. 
The package `robeth` [@Mar20robe] contains some robust tests for linear hypotheses.

The functions in `WRS2` concentrate on the heteroscedastic setting, whereas our focus lies on the homoscedastic case.
The reason is that especially for small samples, estimating the within-sample variance separately for both samples, as is the case under heteroscedasticity, may lead to unreliable estimates.
Moreover, choosing equal sample sizes $m = n$, can protect against a deteriorating test performance of our implemented tests under heteroscedasticity in terms of size and power [@StaShe90robu, p. 179].

The R package `nptest` [@Hel21npte] contains nonparametric versions of the two-sample $t$-test, which is realized by using the permutation and randomization principles described above on the $t$-statistic.
This principle has also been studied in @AbbFri17cont, and, while achieving distribution independence, the test statistic lacks robustness against outliers.

# References
