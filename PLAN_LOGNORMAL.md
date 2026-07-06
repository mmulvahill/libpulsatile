# Plan: configurable pulse parameterization (papers-default)

Support multiple published model specifications as **configurable options**,
defaulting to the published papers while keeping the current natural-scale
"research" behavior available as opt-ins.

## Reference specifications
- **Mulvahill MS thesis** (single-subject): `ln(alpha_k), ln(omega_k) ~ N(mu, sigma^2)`,
  `sigma_alpha, sigma_omega ~ Uniform(0, .)`, Strauss/hard-core/order-statistic
  location priors, `sigma_eps^2 ~ IG`.
- **Horton et al. 2017** (population): same log-normal pulse mass/width with
  `Uniform` SD priors; hierarchical subject + population levels; order-statistic
  locations; `sigma^-2 ~ Gamma`.

Both parameterize pulse mass/width on the **log scale** with **Uniform** SD
priors and **no** per-pulse t-scale. The current code uses natural-scale
truncated-normal + a Student-t (kappa) scale-mixture + half-Cauchy SD prior --
i.e. an additional research direction that matches neither paper. The PDFs are in
`/root/.claude/uploads/6fb472a5-.../` (read locally with `pdftotext -layout` or
the Read tool `pages=`; publisher hosts are egress-blocked).

## Design: three configurable axes, papers by default

| Axis | Default (papers) | Opt-in (current research code) |
|---|---|---|
| Pulse RE scale | **log-normal** (`log theta ~ N(mu, sigma^2/kappa)`) | natural-scale truncated-normal |
| Pulse SD prior | **Uniform(0, .)** | half-Cauchy |
| Per-pulse kappa | **Gaussian** (`student_t_pulses = FALSE`) | Student-t mixture |

`{lognormal, uniform, gaussian}` reproduces the thesis (single) and Horton
(population) exactly. All three current behaviors survive as opt-ins.

## Status
- **master**: PR #18 merged (kappa fix, un-frozen/estimated width SDs,
  `student_t_pulses` toggle, `identifiability_check()`, population `sv` bound guard).
- **Branch `claude/hormone-mcmc-diagnostics-7wxxy5`** @ `07fd780`: log-normal RE
  draw done. `Patient::lognormal_pulses` flag + `SS_DrawRandomEffects` log-normal
  branch (prior `log theta ~ N(mu, sigma^2/kappa)`; symmetric natural-scale RW
  proposal => `1/theta` Jacobian in the ratio; no truncation constant). Catch2
  test verifies the closed form. Compiles + passes on the installed toolchain.

## Remaining tasks (implement -> compile -> Catch2/testthat -> commit, each)

**12. Log-scale mean / SD / tvarscale / birth-death (single-subject).**
- `ss_draw_fixedeffects.h`: log-normal branch uses `(log theta - mu)`; **drop**
  the truncation normalizing constant (`oldint/newint` pnorm terms). Thread a
  `lognormal` bool via the **constructor** (its container type is `bool`, not
  `Patient`, so `parameter_support` can't read the flag otherwise); mu support =
  all reals under log-normal, `>0` on natural.
- `ss_draw_sdrandomeffects.h`: log residuals `(log theta - mu)`; drop truncation.
- `ss_draw_tvarscale.h`: log residuals when kappa is active (skipped under the
  gaussian default).
- `birthdeath.h::add_new_pulse` + `joint_birthdeath.h::add_new_response_pulse`:
  log-normal -> `new_width = exp(rnorm(mu_omega, sigma_omega/sqrt(kappa)))`, no
  reject loop; same for mass.

**13. Uniform SD-prior option (single-subject).**
- `Patient::uniform_sd_prior` flag + `mass_sd_max`/`width_sd_max` in
  `PatientPriors`. Branch `ss_draw_sdrandomeffects.h` (uniform: support
  `0 < sigma < max`, prior ratio cancels; half-Cauchy: current `fourth_part`).
  Default uniform.

**14. Thread axes through R specs + paper-aligned defaults + population/joint C++.**
- `pulse_spec`/`population_spec`/`joint_spec`: add
  `pulse_distribution = c("lognormal","truncnorm")`,
  `sd_prior = c("uniform","half_cauchy")`; flip `student_t_pulses` default to
  **FALSE**. Carry all three via the priors list (same `Rf_asReal` pattern).
- **Recompute log-scale defaults from the papers -- do NOT guess.** Read thesis
  Tables III.1/III.2 and Horton section 4 (approx: mu_log-mass ~ 1,
  mu_log-width ~ 3, sigma ~ Uniform(0, .)). e.g. `width_mean` 42 -> ~3.0 (log),
  `width_sd` -> a log-scale SD with a uniform bound ~2-3.
- Population: log-residual branches in `pop_draw_pulse_sds.h`/`pop_draw_means.h`
  (pop already uses uniform SD prior -- matches). Read + set flags in
  `population.cpp`/`jointsinglesubject.cpp`.
- Docs, `.Rd`, `NAMESPACE`, validation, spec tests. **Breaking-change callout**:
  flips default results and the meaning of `width_mean`/`width_sd`; add a
  regression test that the `truncnorm` option reproduces pre-change behavior.

**15. Sim-study validation.** Simulate with a known log-scale `sigma_omega`; fit
papers-default; confirm `width_sd` recovers the truth and the trace is
well-behaved (vs the old space-filling). `simulate.R` has a log-normal path --
verify it is consistent with the new fit parameterization.

## Definition of done
`devtools::check()` = **0 errors / 0 warnings**; full testthat + Catch2 green;
sim-study shows width-SD recovery. Open a **new PR** off `master` (PR #18 is
merged -- do not reuse it).

## Environment (toolchain not preinstalled; CRAN egress-blocked, apt mirror OK)
```
apt-get update && apt-get install -y --no-install-recommends \
  r-base-dev gfortran libopenblas-dev poppler-utils locales \
  r-cran-rcpp r-cran-rcpparmadillo r-cran-testthat r-cran-devtools \
  r-cran-remotes r-cran-dplyr r-cran-ggplot2 r-cran-tibble \
  r-cran-tidyr r-cran-rlang r-cran-magrittr r-cran-rinside
locale-gen en_US.UTF-8
```
Build/test: `cd R-package && Rscript -e 'devtools::load_all()'` (compiles C++
incl. Catch2), `devtools::test()`, `LANG=en_US.UTF-8 Rscript -e 'devtools::check()'`.
`.Rd`/`NAMESPACE` are hand-maintained; only regenerate if an exported C++
signature changes (this work does not). Branch pushes: `git push --force-with-lease`
(history diverged from merged PR #18 by design).

## Orchestration
**Standing default for this project:** all C++/R implementation, analysis, and
test-execution is carried out by **Opus 4.8 subagents** (Agent tool,
`model: opus`). The lead session (Claude Fable 5) owns planning, sequencing, task
decomposition, and review of returned work -- it does not do the implementation
itself. Pinning the domain/statistics work to a single fixed executor also keeps
a model change from landing partway through a long implementation step.

Each subagent gets a self-contained slice: the specific files, the exact MH math
from the task above, and the compile/test bar. It returns its diff and test
results; the lead reviews and commits. Keep slices small enough to review in one
pass.

---

## Handoff prompt (paste into the new session)

> **Context.** `libpulsatile` -- C++/Rcpp Bayesian MCMC deconvolution of pulsatile
> hormone concentration time-series data (biostatistics). Continue on branch
> `claude/hormone-mcmc-diagnostics-7wxxy5` (@ `07fd780`, based on `master` with
> merged PR #18). Read `CLAUDE.md` and `PLAN_LOGNORMAL.md` first -- the plan has
> the full task breakdown, the exact MH math, the environment setup, and the
> definition of done.
>
> **Goal.** Implement the configurable pulse parameterization: three axes
> (RE scale lognormal/truncnorm, SD prior uniform/half_cauchy, kappa
> gaussian/student-t), defaulting to the papers `{lognormal, uniform, gaussian}`
> which reproduce the Mulvahill thesis and Horton 2017. Done through task 11
> (log-normal RE draw); do tasks 12-15 in order, compiling + Catch2/testthat
> testing each step, to a bar of `devtools::check()` 0 errors / 0 warnings, then
> open a new PR off `master`.
>
> **Orchestration (standing default).** Do the planning/sequencing/review
> yourself; run *all* the C++/R implementation, analysis, and test execution
> through Opus 4.8 subagents (Agent tool, `model: opus`) rather than doing it in
> the lead session. Give each subagent a self-contained slice (files + math from
> `PLAN_LOGNORMAL.md` + the compile/test bar); it reports its diff and test
> results back for your review before you commit.
>
> **Reference PDFs** (publisher hosts egress-blocked; read locally): thesis and
> Horton 2017 in `/root/.claude/uploads/6fb472a5-.../`. Get the log-scale default
> values from the thesis Tables III.1/III.2 and Horton section 4.
