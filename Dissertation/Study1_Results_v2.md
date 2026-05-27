# Results

A simulation that reports only finite-sample bias cannot distinguish an
estimator that is *biased because the sample is small* from one that is
*biased because it is the wrong model*. The two carry opposite practical
implications—the first is remedied by collecting more data, the second is not—
yet at any single sample size they are indistinguishable. The central analytic
device of this study is therefore a decomposition of an estimator's total error
into components with distinct asymptotic fates, made possible by pairing the
Monte Carlo estimates with the estimator's analytic probability limit (plim).

## The decomposition

All estimators are benchmarked against the effects-coding–aligned population
target $\boldsymbol{\Phi}_{\mathrm{ec}}=\mathbf{M}\boldsymbol{\Phi}_{\mathrm{gen}}\mathbf{M}'$,
the *common* estimand on which the sum-score 1-LGCM and the EC-fitted 2-LGCM are
made commensurable (see Method). With a common target in place, the relative
bias of a parameter at sample size $N$ decomposes additively:

$$
\underbrace{\frac{E[\hat\theta_N]-\theta_{\mathrm{ec}}}{\theta_{\mathrm{ec}}}}_{\text{total relative bias}}
=
\underbrace{\frac{E[\hat\theta_N]-\operatorname{plim}\hat\theta}{\theta_{\mathrm{ec}}}}_{\text{finite-sample component}}
+
\underbrace{\frac{\operatorname{plim}\hat\theta-\theta_{\mathrm{ec}}}{\theta_{\mathrm{ec}}}}_{\text{asymptotic component}} .
$$

The **asymptotic component** is the discrepancy that survives as $N\to\infty$:
it is zero if and only if the estimator is consistent for the target, and
nonzero precisely when the fitted model is misspecified in a way it cannot
represent. The **finite-sample component** is the remaining gap between the
expected estimate at $N$ and the plim; for a consistent estimator it is the
entire bias and it vanishes with $N$. We obtain $\operatorname{plim}\hat\theta$
analytically—by fitting each model to the *exact* population moments of the
generating mechanism—so the decomposition is exact and free of Monte Carlo
noise. Because the three quantities share the denominator $\theta_{\mathrm{ec}}$,
they add exactly.

Once the target is the aligned $\boldsymbol{\Phi}_{\mathrm{ec}}$, a *scale*
discrepancy between composite and latent metrics cannot enter the bias at all;
it has been absorbed by the alignment. What remains for the decomposition to
adjudicate are three substantively different sources: a vanishing
**small-sample** effect (nonzero finite-sample component, zero asymptotic
component), a genuine **inconsistency** (nonzero asymptotic component), and—
as shown below—a **selection** artifact that arises when the admissibility
screening applied to Monte Carlo replications removes part of the sampling
distribution.

## Asymptotic component: consistency holds except under unrepresentable misspecification

The asymptotic component is zero, to numerical precision, in every correctly
specified cell—both models, both reliabilities, all $J$, $T\in\{3,5\}$, and
$\rho=.15$ once lag-1 correlated uniqueness (CU) is modeled where identified
(maximum $|$asymptotic relative bias$|<1\times10^{-6}$; Table 1). The
longitudinal EC alignment is thereby validated, and the sum-score 1-LGCM is
shown to be consistent for the *same* estimand whenever correctly specified.

A nonzero asymptotic component appears only in the $\{T=3,\rho=.15\}$ cells,
where the data carry a lag-1 CU that cannot be modeled: adding CU to the
sum-score 1-LGCM at $T=3$ would require 10 parameters for 9 nonredundant moments
($df=-1$), so CU is omitted from both fitted models for comparability. There the
asymptotic bias is large—the slope variance is inflated by 53% (good
reliability) to 234% (poor), and the intercept–slope correlation is attenuated
or sign-reversed (Table 1). Crucially the EC and sum-score models share this
asymptotic bias almost exactly; when CU is present but unrepresentable, the
second-order structure confers no asymptotic protection.

**Table 1. Asymptotic component (plim $-$ target, relative).**

| Condition | Model | Reliability | $\phi_{22}$ | $\rho_{IS}$ |
|---|---|---|---|---|
| All correctly specified cells | both | both | $<10^{-6}$ | $<10^{-6}$ |
| $T=3,\rho=.15$ (CU omitted) | EC | $R^2=.45$ | $+.53$ | $-.70$ |
| $T=3,\rho=.15$ | EC | $R^2=.16$ | $+2.29$ | $-1.60$ |
| $T=3,\rho=.15$ | Sum | $R^2=.45$ | $+.54$ | $-.71$ |
| $T=3,\rho=.15$ | Sum | $R^2=.16$ | $+2.34$ | $-1.61$ |

## The decomposition in practice: same symptom, different disease

Table 2 applies the decomposition to the slope variance $\phi_{22}$ in three
conditions whose *total* bias is similar in magnitude but whose *composition*—
and therefore whose diagnosis—is entirely different. Figure 1 shows the same
information as convergence toward the $N=\infty$ plim.

**Table 2. Decomposition of total relative bias in $\phi_{22}$
(total $=$ finite-sample $+$ asymptotic).**

| Case | Model, condition | $N$ | total | finite-sample | asymptotic |
|---|---|---|---|---|---|
| **A** | EC, $T=3,\rho=0$ ($R^2=.45$) | 200 | $.47$ | $.47$ | $.00$ |
| | | 500 | $.24$ | $.24$ | $.00$ |
| | | 800 | $.16$ | $.16$ | $.00$ |
| **B** | Sum, $T=3,\rho=.15$ ($R^2=.45$) | 200 | $.78$ | $.24$ | $.54$ |
| | | 500 | $.59$ | $.04$ | $.54$ |
| | | 800 | $.55$ | $.01$ | $.54$ |
| **C** | EC, $T=3,\rho=.15$ ($R^2=.45$) | 200 | $.45$ | $-.08$ | $.53$ |
| | | 500 | $.22$ | $-.32$ | $.53$ |
| | | 800 | $.13$ | $-.40$ | $.53$ |

**Case A — small sample (a benign, vanishing artifact).** Under correct
specification at $T=3$, the EC model's slope-variance bias is large at small
samples ($+47\%$ at $N=200$) but is *entirely* finite-sample: the asymptotic
component is zero, and the bias falls monotonically toward it ($+16\%$ at
$N=800$). The slope variance is simply hard to estimate from a single increment;
the estimator is consistent and the problem is remedied by more data (or more
occasions). Reporting the $N=200$ figure alone would wrongly indict a correct
model.

**Case B — genuine inconsistency (the bias that more data cannot fix).** With
unrepresentable CU, the sum-score model's bias of the same apparent size
($+55\%$ at $N=800$) decomposes almost entirely into the asymptotic component
($+.54$), the finite-sample part having already shrunk to $\approx.01$. This
estimator converges—cleanly, and to the *wrong* value. No sample size repairs
it. The decomposition draws a sharp line between Case A and Case B even though
their total biases are nearly identical: same symptom, opposite disease.

**Case C — selection (a screening artifact, not consistency).** The EC model in
the *identical* misspecified condition appears, on admissible replications, to
have a *smaller* and *shrinking* total bias ($+45\%\to+13\%$). Taken at face
value this would suggest the EC model is the better choice at $T=3$. The
decomposition shows otherwise: its asymptotic component ($+.53$) equals the
sum-score model's, but its finite-sample component is large and *negative*, and
grows more negative with $N$ ($-.08\to-.40$)—the admissible-sample mean is
moving *away* from the plim. The cause is that the EC model's unconstrained
probability limit in this cell is itself **inadmissible**: the population
solution places a negative disturbance on the final wave ($s_3=-.092$), a
Heywood case at $N=\infty$. The admissibility screen applied to each replication
therefore removes the upper tail of the slope-variance distribution, so the
admissible-only mean targets a constrained boundary rather than the plim. The
apparently milder bias of the EC model at $T=3$ is thus a selection artifact of
conditioning on admissibility, not evidence of better recovery.

![Figure 1. Total relative bias against sample size, with the analytic plim as
the $N=\infty$ anchor.](figure1_bias_convergence.png)

*Figure 1. Total relative bias of $\phi_{22}$ (top) and $\rho_{IS}$ (bottom)
versus $N$; the rightmost point ($\infty$) is the analytic plim. A line that
descends to zero (e.g., $T=3,\rho=0$) is all finite-sample; a line that flattens
onto a nonzero plim ($T=5$ vs.\ $T=3,\rho=.15$ comparison) carries an asymptotic
component; a line that departs from its own plim (EC at $T=3,\rho=.15$) is under
selection.*

That Case A and Case B are genuinely different is confirmed at $T=5$, where the
finite-sample component is the *only* component for both models. There the
sum-score slope-variance bias under $\rho=.15$ is finite-sample throughout
(poor reliability: $+.26$ at $N=200$, $+.07$ at $N=800$, asymptotic component
$0$) and the EC model is essentially unbiased at every $N$. With enough
occasions to identify and model CU, both estimators are consistent and behave
accordingly; the decomposition attributes their finite-$N$ differences to
sampling, not to misspecification.

## Diagnostic signatures

The three components leave distinct, observable signatures that an applied
analyst can read without access to the plim (Table 3). A consistent estimator
under small samples (Case A) retains near-nominal coverage—its intervals are
merely wide—and its bias declines with $N$. A genuinely inconsistent estimator
(Case B) shows coverage that *collapses* as $N$ grows and the intervals tighten
around the wrong value (e.g., $\phi_{22}$ coverage $.38$, $\rho_{IS}$ coverage
$.16$ at $N=800$ under poor reliability), even while it converges without
complaint (admissibility $\ge .92$). A selection regime (Case C) announces
itself through a high *inadmissibility* rate (down to 37% admissible at $N=200$)
rather than through poor coverage. The practical lesson is that clean
convergence is not a marker of recovery: the parsimonious sum-score model
converges reliably to biased estimates, whereas the EC model's inadmissibility,
though inconvenient, is the more honest warning that the growth structure is not
recoverable at $T=3$ under CU.

**Table 3. Observable signatures of each component (selected rows).**

| Case | Model, condition | $N$ | total bias | $\phi_{22}$ cov | $\rho_{IS}$ cov | % admissible |
|---|---|---|---|---|---|---|
| A (small sample) | EC, $T=3,\rho=0$, poor | 200 | $.86$ | $.99$ | $.79$ | $.36$ |
| B (inconsistency) | Sum, $T=3,\rho=.15$, poor | 200 | $2.62$ | $.76$ | $.41$ | $.92$ |
| B (inconsistency) | Sum, $T=3,\rho=.15$, poor | 800 | $2.31$ | $.38$ | $.16$ | $1.00$ |
| recovered | EC, $T=5,\rho=.15$, good | 800 | $-.00$ | $.95$ | $.95$ | $1.00$ |

## Summary of findings

The decomposition resolves the apparent failures of these growth estimators into
their true causes. (1) After alignment to the common estimand
$\boldsymbol{\Phi}_{\mathrm{ec}}$, scale discrepancy is removed by construction.
(2) Under correct specification both the EC-fitted 2-LGCM and the sum-score
1-LGCM are consistent; with $T=5$ occasions this is realized in practice, the EC
model being essentially unbiased with nominal coverage. (3) At $T=3$, growth-
variance bias under correct specification is *entirely finite-sample*—large but
vanishing—and must not be mistaken for a defect of the model. (4) When
correlated uniqueness is present but unrepresentable ($T=3$), both estimators
incur the *same* asymptotic bias, which no sample size removes; the second-order
EC model offers no asymptotic protection. (5) On admissible replications the
EC model nonetheless *appears* less biased at $T=3$, but this is a selection
artifact: its unconstrained plim is an inadmissible Heywood solution, so
conditioning on admissibility distorts the comparison. Distinguishing these four
sources—scale, small sample, inconsistency, and selection—requires the plim;
finite-sample bias alone conflates them.
