# Appendix A. Analytic Population Targets and Probability Limits

This appendix documents how the analytic benchmarks used throughout Study 1 are
obtained. Two quantities are needed for every design cell: (i) the population
moments implied by the data-generating model, from which (ii) the probability
limit (plim) of each fitted estimator is computed by minimizing the maximum-
likelihood (ML) discrepancy at the population. Pairing each plim with the
effects-coding–aligned target $\boldsymbol\Phi_{\mathrm{ec}}$ permits the
decomposition of total error into finite-sample and asymptotic components
reported in the Results. All quantities are exact (free of Monte Carlo error).

## A.1 Generating model

Let $T$ be the number of measurement occasions and $J$ the number of indicators
per occasion. The model has three layers:

$$
\begin{aligned}
\text{growth:}\quad & (\alpha,\beta)' \sim N\!\big(\mu_{\mathrm{gen}},\,\boldsymbol\Phi_{\mathrm{gen}}\big),\\[2pt]
\text{latent:}\quad & \eta_t=\alpha+\beta\,c_t+\zeta_t,\qquad \zeta_t\sim N(0,\psi_t)\ \text{independent},\quad c_t=t-1,\\[2pt]
\text{measurement:}\quad & y_{jt}=\tau_{jt}+\lambda_{jt}\,\eta_t+e_{jt},\qquad
\mathbf e_j\sim N(0,\boldsymbol\Sigma_{e_j}),
\end{aligned}
$$

where $\mathbf e_j=(e_{j1},\dots,e_{jT})'$ is the residual vector of indicator
$j$ across time, with within-indicator lag-1 correlated uniqueness (CU) and
zero cross-indicator dependence. The growth basis is the $T\times2$ matrix
$\mathbf B$ with rows $(1,c_t)$.

## A.2 Population moments in closed form

**Proposition A.1 (item-level moments).** Under the model of A.1, the population
mean vector and covariance matrix of the $JT$ observed indicators are

$$
E[\mathbf y]=\boldsymbol\tau+\boldsymbol\Lambda\,\mathbf B\,\mu_{\mathrm{gen}},
\qquad
\operatorname{Cov}(\mathbf y)=\boldsymbol\Lambda\,
\big(\mathbf B\,\boldsymbol\Phi_{\mathrm{gen}}\,\mathbf B'+\boldsymbol\Psi\big)\,
\boldsymbol\Lambda'+\boldsymbol\Theta,
$$

where indicators are stacked in time-major order ($\text{index}=(t-1)J+j$),
$\boldsymbol\Psi=\operatorname{diag}(\psi_1,\dots,\psi_T)$, the loading matrix
$\boldsymbol\Lambda$ is $JT\times T$ with
$\boldsymbol\Lambda[(t-1)J+j,\,u]=\lambda_{jt}\,\mathbb 1\{u=t\}$, and the
residual covariance $\boldsymbol\Theta$ is block-structured,
$\boldsymbol\Theta[(t-1)J+j,\,(u-1)J+k]=\mathbb 1\{j=k\}\,\boldsymbol\Sigma_{e_j}[t,u]$.

*Proof.* Write $\boldsymbol\eta=\mathbf B(\alpha,\beta)'+\boldsymbol\zeta$. Since
$(\alpha,\beta)\perp\boldsymbol\zeta$,
$E[\boldsymbol\eta]=\mathbf B\mu_{\mathrm{gen}}$ and
$\operatorname{Cov}(\boldsymbol\eta)=\mathbf B\boldsymbol\Phi_{\mathrm{gen}}\mathbf B'+\boldsymbol\Psi$.
Stacking the measurement equations gives
$\mathbf y=\boldsymbol\tau+\boldsymbol\Lambda\boldsymbol\eta+\mathbf e$ with
$\boldsymbol\eta\perp\mathbf e$, so $E[\mathbf y]=\boldsymbol\tau+\boldsymbol\Lambda E[\boldsymbol\eta]$
and $\operatorname{Cov}(\mathbf y)=\boldsymbol\Lambda\operatorname{Cov}(\boldsymbol\eta)\boldsymbol\Lambda'+\operatorname{Cov}(\mathbf e)$.
Cross-indicator independence makes $\operatorname{Cov}(\mathbf e)=\boldsymbol\Theta$ block-diagonal in $j$. $\square$

Equivalently, in scalar form,

$$
\operatorname{Cov}(\eta_t,\eta_u)=\operatorname{Var}(\alpha)+(c_t+c_u)\operatorname{Cov}(\alpha,\beta)+c_tc_u\operatorname{Var}(\beta)+\mathbb 1\{t=u\}\psi_t,
$$
$$
\operatorname{Cov}(y_{jt},y_{ku})=\lambda_{jt}\lambda_{ku}\operatorname{Cov}(\eta_t,\eta_u)+\mathbb 1\{j=k\}\,\boldsymbol\Sigma_{e_j}[t,u],
$$

and the CU block is $\boldsymbol\Sigma_{e_j}=\mathbf D_j\mathbf R\mathbf D_j$ with
$\mathbf D_j=\operatorname{diag}(\sqrt{\theta_{j1}},\dots,\sqrt{\theta_{jT}})$ and
$\mathbf R[t,u]=\mathbb 1\{t=u\}+\rho\,\mathbb 1\{|t-u|=1\}$.

**Proposition A.2 (sum-score moments).** Let $S_t=\sum_{j}y_{jt}$ be the wave-$t$
sum score and $\mathbf A$ the $T\times JT$ aggregation matrix,
$\mathbf A[t,(t-1)J+j]=1$. Then $E[\mathbf S]=\mathbf A\,E[\mathbf y]$ and
$\operatorname{Cov}(\mathbf S)=\mathbf A\operatorname{Cov}(\mathbf y)\mathbf A'$,
which in scalar form is

$$
\operatorname{Cov}(S_t,S_u)=s_t\,s_u\,\operatorname{Cov}(\eta_t,\eta_u)+\sum_{j}\boldsymbol\Sigma_{e_j}[t,u],
\qquad s_t=\sum_{j}\lambda_{jt}.
$$

**Remark A.1 (loadings enter only through their sum).** Proposition A.2 shows
that both the sum-score moments and—because the wave-$t$ first-order factor in
the EC model carries total loading $s_t$—the aligned target depend on the
measurement model only through the per-wave loading sums $s_t$ (and intercept
sums $c_t=\sum_j\tau_{jt}$). Two loading configurations with the same $s_t$ are
therefore observationally equivalent for the growth moments. This is why the
tau-equivalent and congeneric conditions, constructed to hold $s_t$ fixed,
produce identical bias, and why the longitudinal alignment matrix reduces to a
scalar $\mathbf M=s\,\mathbf I$ when loadings are time-invariant ($s_t\equiv s$).

## A.3 Aligned target

Effects-coding scale alignment (derived in the Method) maps the generating
growth parameters to the common estimand against which both models are
benchmarked,

$$
\boldsymbol\Phi_{\mathrm{ec}}=\mathbf M\,\boldsymbol\Phi_{\mathrm{gen}}\,\mathbf M',
\qquad
\boldsymbol\kappa_{\mathrm{ec}}=\mathbf M\,\mu_{\mathrm{gen}}+(\mathbf B'\mathbf B)^{-1}\mathbf B'\mathbf c,
\qquad
\mathbf M=(\mathbf B'\mathbf B)^{-1}\mathbf B'\mathbf D\mathbf B,
$$

with $\mathbf D=\operatorname{diag}(s_t)$ and $\mathbf c=(c_1,\dots,c_T)'$ the
intercept sums. $\boldsymbol\Phi_{\mathrm{ec}}$ places the composite and latent
estimators on a single metric, so that any discrepancy from it reflects
estimation or projection error rather than a difference of scale.

## A.4 Probability limits via population-moment fitting

For a fitted model with implied moments $\Sigma(\theta),\mu(\theta)$, the ML
estimator minimizes the discrepancy

$$
F_{\mathrm{ML}}(\theta;\mathbf S,\bar{\mathbf x})=\log|\Sigma(\theta)|+\operatorname{tr}\!\big(\mathbf S\,\Sigma(\theta)^{-1}\big)-\log|\mathbf S|-p+(\bar{\mathbf x}-\mu(\theta))'\Sigma(\theta)^{-1}(\bar{\mathbf x}-\mu(\theta)),
$$

where $\mathbf S,\bar{\mathbf x}$ are the sample covariance and mean and $p$ the
number of observed variables. Because $\mathbf S\xrightarrow{p}\Sigma_0$ and
$\bar{\mathbf x}\xrightarrow{p}\mu_0$, the estimator converges to the
pseudo-true value (White, 1982)

$$
\hat\theta_N\xrightarrow{p}\theta^\*=\arg\min_\theta F_{\mathrm{ML}}(\theta;\Sigma_0,\mu_0),
$$

i.e., the **plim**. We therefore obtain $\theta^\*$ directly by minimizing
$F_{\mathrm{ML}}$ at the exact population moments of Propositions A.1–A.2, rather
than by large-sample simulation. Operationally we supply $\Sigma_0$ and $\mu_0$
(item-level for the EC 2-LGCM, sum-level for the 1-LGCM) to the SEM fitting
routine as summary statistics, fitting the *same* model syntax used in the Monte
Carlo study so that the analytic and empirical analyses are exactly paired. Two
points warrant note:

1. **Sample size is immaterial to the plim.** The point estimates minimizing
   $F_{\mathrm{ML}}$ do not depend on the nominal $N$ assigned to the population
   moments; $N$ rescales only the $\chi^2$ statistic ($\chi^2=(N-1)F_{\mathrm{ML}}$
   at the optimum) and the standard errors, not the location of the minimum.
   We use $N=10^5$ and disable the default covariance rescaling so that the fed
   moments are used exactly.

2. **The minimizer is unconstrained.** No inequality (admissibility)
   constraints are imposed, so $\theta^\*$ is the unconstrained ML probability
   limit. Under severe misspecification $\theta^\*$ may be inadmissible (e.g., a
   negative wave-specific disturbance variance—a population Heywood case). When
   this occurs, the admissibility screening applied in the Monte Carlo study
   makes the admissible-conditional estimator converge to a *constrained*
   boundary value distinct from $\theta^\*$; the resulting divergence between the
   admissible-sample mean and the unconstrained plim is itself diagnostic of a
   selection regime (see Results, Case C).

**Proposition A.3 (consistency under correct specification).** If the fitted
model is correctly specified—its implied moment structure can reproduce
$(\Sigma_0,\mu_0)$ at the aligned parameters—then $F_{\mathrm{ML}}$ attains its
minimum of $0$ there and $\theta^\*=\boldsymbol\Phi_{\mathrm{ec}}$ (and
$\boldsymbol\kappa_{\mathrm{ec}}$); the estimator is consistent for the target.
Otherwise $\min_\theta F_{\mathrm{ML}}>0$ and $\theta^\*\neq$ target, the gap
being the asymptotic bias; the minimized $F_{\mathrm{ML}}$ maps to the population
fit indices (RMSEA, CFI) reported for the misspecified cells.

## A.5 Bias decomposition

With the target $\theta_{\mathrm{ec}}$ and plim $\theta^\*=\operatorname{plim}\hat\theta$
in hand, the relative bias at sample size $N$ decomposes additively (common
denominator $\theta_{\mathrm{ec}}$):

$$
\underbrace{\frac{E[\hat\theta_N]-\theta_{\mathrm{ec}}}{\theta_{\mathrm{ec}}}}_{\text{total}}
=
\underbrace{\frac{E[\hat\theta_N]-\theta^\*}{\theta_{\mathrm{ec}}}}_{\text{finite-sample}}
+
\underbrace{\frac{\theta^\*-\theta_{\mathrm{ec}}}{\theta_{\mathrm{ec}}}}_{\text{asymptotic}} .
$$

The asymptotic component is nonzero iff the estimator is inconsistent for the
target (Proposition A.3); the finite-sample component vanishes as $N\to\infty$
for a consistent estimator. This decomposition is what distinguishes a
small-sample artifact (asymptotic component $=0$) from a genuine inconsistency
(asymptotic component $\neq 0$), and—when the finite-sample component fails to
vanish because $\theta^\*$ is inadmissible—a selection artifact.

## A.6 Numerical verification

The closed-form moments of Propositions A.1–A.2 were verified against a
brute-force simulation of $N=2\times10^6$ observations from the generating
mechanism: the maximum absolute discrepancy between the analytic and empirical
moments was $\le 3\times10^{-3}$ (Monte Carlo error at this $N$) for both
item-level and sum-score covariances and means. The plim computation was
verified against a large-sample Monte Carlo ($N=10^5$): the simulated mean
estimate matched the analytic plim to three decimals (e.g., slope variance
$\hat\phi_{22}=.0731$ vs.\ plim $.0729$). R code implementing
Propositions A.1–A.2 (`pop_moments`) and the plim computation appears in the
replication materials (`Study1_calculation.R`, `Study1_plim.R`).

---

**Reference.** White, H. (1982). Maximum likelihood estimation of misspecified
models. *Econometrica, 50*(1), 1–25.
