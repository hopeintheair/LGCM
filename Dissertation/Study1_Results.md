# Results

All estimators are evaluated against a single, analytically defined population
target. Under the generating second-order growth model, effects-coding (EC)
scale alignment maps the generating growth parameters
$\boldsymbol{\theta}_{\mathrm{gen}}=(\boldsymbol{\kappa}_{\mathrm{gen}},
\boldsymbol{\Phi}_{\mathrm{gen}})$ to the aligned target

$$
\boldsymbol{\Phi}_{\mathrm{ec}} = \mathbf{M}\,\boldsymbol{\Phi}_{\mathrm{gen}}\,\mathbf{M}',
\qquad
\boldsymbol{\kappa}_{\mathrm{ec}} = \mathbf{M}\,\boldsymbol{\kappa}_{\mathrm{gen}}
  + (\mathbf{B}'\mathbf{B})^{-1}\mathbf{B}'\mathbf{c},
\qquad
\mathbf{M}=(\mathbf{B}'\mathbf{B})^{-1}\mathbf{B}'\mathbf{D}\mathbf{B},
$$

where $\mathbf{B}$ is the growth basis, $\mathbf{D}=\operatorname{diag}(s_t)$
collects the per-wave loading sums $s_t=\sum_j\lambda_{jt}$, and
$\mathbf{c}$ the per-wave intercept sums. Because both the sum-score
1-LGCM and the EC-fitted 2-LGCM depend on the measurement model only through
$(s_t,c_t)$, $\boldsymbol{\Phi}_{\mathrm{ec}}$ is the *common* estimand for
both models, and "bias" throughout denotes departure from
$\boldsymbol{\Phi}_{\mathrm{ec}}$ rather than from the raw latent metric.
In Study 1 loadings are time-invariant, so $s_t\equiv s$, $\mathbf{M}=s\,\mathbf{I}$,
and the alignment reduces to a uniform rescaling that preserves correlations
(the nontrivial, time-varying $\mathbf{M}$ is exercised in Study 2).

We report the results in the order that isolates the two sources of error.
We first establish the *asymptotic* behavior of each estimator through its
probability limit (plim), computed analytically by fitting each model to the
exact population moments of the generating mechanism; these point estimates are
the values to which the Monte Carlo estimates converge as $N\to\infty$. We
then characterize *finite-sample* behavior, stratified throughout by the number
of measurement occasions $T$, which—as shown below—governs which of the two
error sources dominates. We never pool across $T$.

## Specification conditions

The design crosses two manipulations that bear on model specification:
correlated uniqueness (CU) in the data ($\rho\in\{0,.15\}$) and the number of
occasions ($T\in\{3,5\}$). Table 0 maps these onto specification status. When
$\rho=0$ both fitted models are correctly specified. When $\rho=.15$ the
data contain a lag-1 CU that the fitted model should accommodate; this is
possible at $T=5$ (the sum-score 1-LGCM has $df=6$) but **not** at $T=3$,
where adding lag-1 CU to the sum-score model would require 10 parameters for
only 9 nonredundant moments ($df=-1$). For comparability the CU terms are
therefore omitted from *both* fitted models at $T=3$, making the
$\{T=3,\rho=.15\}$ cells the only misspecified condition in the design.

**Table 0. Specification map.**

| Data CU | Fitted CU | Design cells | Identification | Specification |
|---|---|---|---|---|
| No ($\rho=0$) | — | $\rho=0,\,T=3$ | identified | correct |
| No ($\rho=0$) | — | $\rho=0,\,T=5$ | identified | correct |
| Yes ($\rho=.15$) | omitted | $\rho=.15,\,T=3$ | CU not identified in 1-LGCM ($df=-1$) | **misspecified** |
| Yes ($\rho=.15$) | modeled | $\rho=.15,\,T=5$ | identified ($df=6$) | correct |

The loading-shape factor (tau-equivalent vs. congeneric) is not a separate
specification condition: since both estimators depend on the loadings only
through $s_t$, and the congeneric loadings were constructed to hold $s_t$
fixed, the two loading shapes are observationally equivalent for the growth
moments. Results are therefore averaged over loading shape (and over the number
of indicators $J$, which only rescales $s$); the full per-cell results
appear in the Supplement.

## Asymptotic behavior: probability limits

**Consistency under correct specification.** Across every correctly specified
cell—both models, both reliabilities, all $J$, $T\in\{3,5\}$,
$\rho\in\{0,.15\}$ with CU modeled where identified—the plim of each estimator
equals the aligned target to numerical precision (Table 1a; maximum absolute
relative bias $<1\times10^{-6}$ for every parameter). Two conclusions follow.
First, the longitudinal EC alignment is internally validated: the 2-LGCM
recovers exactly the $\boldsymbol{\Phi}_{\mathrm{ec}}$ implied by the
derivation. Second, the sum-score 1-LGCM is *also* consistent for the same
estimand whenever it is correctly specified—including the $\{T=5,\rho=.15\}$
condition once lag-1 CU is modeled. Differences between the two models under
correct specification are therefore purely finite-sample, not asymptotic.

**Table 1a. Maximum $|$asymptotic relative bias$|$ across all correctly
specified cells.**

| | $\phi_{11}$ | $\phi_{22}$ | $\phi_{12}$ | $\rho_{IS}$ |
|---|---|---|---|---|
| max $|$arbias$|$ | $2.2\times10^{-7}$ | $4.2\times10^{-7}$ | $4.5\times10^{-7}$ | $4.8\times10^{-7}$ |

**Asymptotic bias under CU omission.** Bias survives at $N\to\infty$ only in
the misspecified $\{T=3,\rho=.15\}$ cells, and there it is severe (Table 1b).
For the EC 2-LGCM the slope-variance plim exceeds the target by 53% under good
reliability ($R^2=.45$) and by 229% under poor reliability ($R^2=.16$); the
intercept–slope association is attenuated by 70% (good) and *reverses sign*
under poor reliability (population $\rho_{IS}$ of $+.40$ projecting to
$\approx-.24$). Critically, the sum-score 1-LGCM exhibits essentially the same
asymptotic bias ($+54\%$ and $+234\%$; $\rho_{IS}$ bias $-.71$ and
$-1.61$). When CU is present but unmodeled, the second-order machinery of the
EC model confers **no asymptotic protection**: both estimators are projecting
the same omitted covariance onto the growth factors.

**Table 1b. Asymptotic relative bias in the misspecified condition
($T=3,\rho=.15$).**

| Model | Reliability | $\phi_{11}$ | $\phi_{22}$ | $\phi_{12}$ | $\rho_{IS}$ |
|---|---|---|---|---|---|
| EC 2-LGCM | $R^2=.45$ | $+.22$ | $+.53$ | $-.60$ | $-.70$ |
| EC 2-LGCM | $R^2=.16$ | $+.93$ | $+2.29$ | $-2.57$ | $-1.60$ |
| Sum 1-LGCM | $R^2=.45$ | $+.22$ | $+.54$ | $-.61$ | $-.71$ |
| Sum 1-LGCM | $R^2=.16$ | $+.95$ | $+2.34$ | $-2.63$ | $-1.61$ |

The growth means $(\kappa_I,\kappa_S)$ were recovered without appreciable
asymptotic bias in every condition, misspecified cells included
($|$relative bias$|<.01$; not tabled). Misspecification of the CU structure
distorts the second-order moments of growth while leaving the first moments
intact.

## Finite-sample behavior

Figure 1 superimposes the Monte Carlo relative bias at $N\in\{200,500,800\}$
on the plim ($N=\infty$), separately for the slope variance $\phi_{22}$ and
the intercept–slope correlation $\rho_{IS}$. The four columns make the
governing role of $T$ explicit.

![Figure 1. Finite-sample bias against sample size, with the analytic
probability limit as the asymptotic anchor.](figure1_bias_convergence.png)

*Figure 1. Relative bias of $\phi_{22}$ (top) and $\rho_{IS}$ (bottom) as a
function of sample size; the rightmost point ($\infty$) is the analytic plim.
Color distinguishes the two models, line type the two reliability levels.
Bias decays to zero with $N$ only where the plim is itself zero.*

**$T=5$: the estimand is recovered (Table 2, Panel B).** With five occasions
the EC 2-LGCM is essentially unbiased at all sample sizes—$\phi_{22}$ relative
bias $\le.004$ under good reliability and $\le.066$ under poor reliability
even at $N=200$—with confidence-interval coverage at the nominal level
($.94$–$.98$). The sum-score model matches this when $\rho=0$ and shows
only a small, reliability-dependent finite-sample bias under unmodeled-free CU
($\rho=.15$, poor reliability: $\phi_{22}$ bias $.26$ at $N=200$ shrinking
to $.07$ at $N=800$), consistent with its zero plim. At $T=5$ the practical
behavior tracks the asymptotic theory closely for both models.

**$T=3$: finite-sample bias dominates (Table 2, Panel A).** With only three
occasions a large *upward* finite-sample bias in the slope variance appears even
under correct specification ($\rho=0$). For the EC model under good
reliability it falls from $+.47$ at $N=200$ to $+.16$ at $N=800$; under
poor reliability from $+.86$ to $+.47$. Because the corresponding plim is
exactly zero (Table 1a), this is unambiguously a small-$T$ finite-sample
artifact—the slope variance is estimated from a single increment—rather than a
metric or specification problem. It nonetheless renders the $T=3$ growth-
variance estimates untrustworthy at the sample sizes common in applied work.

**Convergence behavior and admissibility (Table 2).** The $T=3$ condition also
exposes a selection problem that bears directly on how these models are used in
practice. The EC 2-LGCM frequently yields inadmissible solutions at small
samples—down to 37% admissible under poor reliability at $N=200$—whereas the
sum-score model converges almost always (85–100% admissible). Convergence,
however, is not a signal of trustworthiness. In the misspecified
$\{T=3,\rho=.15\}$ cell the sum-score model converges in 92–100% of replications
and its admissible-sample mean *converges to its badly biased plim*
($\phi_{22}$ bias $+.55$ reliable, $+2.31$ poor at $N=800$), with
coverage collapsing accordingly ($\phi_{22}$ coverage $.38$; $\rho_{IS}$
coverage $.16$ under poor reliability at $N=800$). By contrast, the EC
model's admissible-sample mean in the same cell *moves away* from its (equally
biased) plim as $N$ grows, because the inadmissible solutions discarded at
small $N$ are disproportionately those with extreme slope-variance estimates;
the retained subset is a selected, non-representative sample. A naive comparison
of admissible-only results would thus make the EC model appear less biased at
$T=3$ than the sum-score model, even though the two share the same probability
limit. The cleaner convergence of the sum-score model is therefore not evidence
of better recovery—if anything it removes the very warning (inadmissibility)
that flags the underlying misspecification.

**Table 2. Finite-sample performance, stratified by $T$ (averaged over loading
shape and $J$); selected rows.**

*Panel A — $T=3$*

| Model | $\rho$ | Reliability | $N$ | $\phi_{22}$ rbias | $\rho_{IS}$ rbias | $\phi_{22}$ cov | $\rho_{IS}$ cov | % admiss |
|---|---|---|---|---|---|---|---|---|
| EC | .00 | $R^2=.45$ | 200 | .47 | $-.44$ | .97 | .85 | .65 |
| EC | .00 | $R^2=.45$ | 800 | .16 | $-.12$ | .97 | .90 | .84 |
| EC | .00 | $R^2=.16$ | 200 | .86 | $-.88$ | .99 | .79 | .36 |
| Sum | .15 | $R^2=.45$ | 800 | .55 | $-.61$ | .82 | .64 | .98 |
| Sum | .15 | $R^2=.16$ | 800 | 2.31 | $-1.52$ | .38 | .16 | 1.00 |

*Panel B — $T=5$*

| Model | $\rho$ | Reliability | $N$ | $\phi_{22}$ rbias | $\rho_{IS}$ rbias | $\phi_{22}$ cov | $\rho_{IS}$ cov | % admiss |
|---|---|---|---|---|---|---|---|---|
| EC | .15 | $R^2=.45$ | 200 | .003 | .065 | .95 | .95 | .97 |
| EC | .15 | $R^2=.45$ | 800 | $-.000$ | .017 | .95 | .95 | 1.00 |
| EC | .15 | $R^2=.16$ | 800 | .011 | .044 | .96 | .94 | .95 |
| Sum | .15 | $R^2=.16$ | 200 | .26 | $-.49$ | .97 | .83 | .68 |
| Sum | .15 | $R^2=.16$ | 800 | .065 | $-.075$ | .97 | .91 | .85 |

## The asymptotic fit-index paradox

A final result concerns whether global fit would alert an analyst to the
misspecification responsible for the $T=3$ asymptotic bias. It would not.
Table 3 reports the *population* fit of the misspecified models—fit computed at
the exact generating moments, free of sampling noise. The EC 2-LGCM, despite a
53–229% bias in the slope variance, attains a population RMSEA of $.044$ and a
CFI of $.97$ under good reliability—values that pass every conventional
close-fit guideline. The sum-score 1-LGCM is more striking still: at $T=3$ it
is nearly saturated ($df=1$) and reproduces the sum-score moments *exactly*
(population RMSEA $=0$, CFI $=1.00$), yet its growth parameters carry the
identical bias. Good fit does not imply accurate recovery of growth parameters,
and the parsimony of the sum-score model actively conceals the misspecification:
the projection of the omitted covariance onto the growth factors is absorbed
without residual misfit. For contrast, the correctly specified
$\{T=5,\rho=.15\}$ models show both perfect population fit and zero asymptotic
bias.

**Table 3. Population (asymptotic) fit versus parameter bias.**

| Cell | Model | Reliability | RMSEA | SRMR | CFI | $df$ | $\phi_{22}$ arbias | $\rho_{IS}$ arbias |
|---|---|---|---|---|---|---|---|---|
| **MISSPEC** $T=3,\rho=.15$ | EC | $R^2=.45$ | .044 | .020 | .968 | 76.5 | $+.53$ | $-.70$ |
| **MISSPEC** $T=3,\rho=.15$ | EC | $R^2=.16$ | .044 | .031 | .861 | 76.5 | $+2.29$ | $-1.60$ |
| **MISSPEC** $T=3,\rho=.15$ | Sum | $R^2=.45$ | .000 | .000 | 1.000 | 1.0 | $+.54$ | $-.71$ |
| **MISSPEC** $T=3,\rho=.15$ | Sum | $R^2=.16$ | .000 | .000 | 1.000 | 1.0 | $+2.34$ | $-1.61$ |
| Correct $T=5,\rho=.15$ | EC | $R^2=.45$ | .000 | .000 | 1.000 | 206.5 | $.00$ | $.00$ |
| Correct $T=5,\rho=.15$ | Sum | $R^2=.45$ | .000 | .000 | 1.000 | 6.0 | $.00$ | $.00$ |

## Summary of findings

(1) After longitudinal effects-coding alignment, both the EC-fitted 2-LGCM and
the sum-score 1-LGCM are consistent for the same population growth structure
$\boldsymbol{\Phi}_{\mathrm{ec}}$ whenever the model is correctly specified;
the alignment derivation is validated to numerical precision. (2) With
$T=5$ occasions this consistency is realized in practice, the EC model being
essentially unbiased with nominal coverage across sample sizes. (3) With
$T=3$ occasions, growth-variance estimates are dominated by a large
finite-sample bias that is present even under correct specification and that
vanishes asymptotically. (4) When correlated uniqueness is present but cannot be
modeled ($T=3$), both estimators incur the same severe *asymptotic* bias,
attenuating or reversing the intercept–slope association; the second-order
EC model offers no protection. (5) This misspecification is invisible to global
fit indices—indeed the parsimonious sum-score model conceals it most
completely—and clean convergence is not a marker of recovery: the sum-score
model converges reliably to biased estimates, while the EC model's
inadmissibility, though inconvenient, is the more honest signal.
