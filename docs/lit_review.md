# Literature review: alternatives and extensions to BD-MCMC deconvolution

**Date:** 2026-06-07
**Status:** Background research. Nothing here is scheduled work; it is the
evidence base behind the future-appendix entries (I-O) in
`2026-06-07-vignette-roadmap-design.md`.

---

Several substantive methodological developments since the ~2015-2017 work
underlying libpulsatile are worth investigating.

## 1. HormoneBayes: SMC replaces BD-MCMC entirely (2024)

The most directly relevant paper. Voliotis, Abbara, Prague, Veldhuis, et al.
sidestep the transdimensional sampling problem by reformulating the generative
model [1]:

- Instead of estimating the number/location of pulses (requiring BD-MCMC), they
  model the underlying hypothalamic ON/OFF switching process as a two-state
  Markov chain with transition rates tau_ON and tau_OFF.
- The observed LH is driven by both a pulsatile signal P_t and a basal signal
  B_t.
- Inference uses a Gibbs sampler where:
  - Sequential Monte Carlo (SMC) with ancestral sampling handles the latent
    state trajectories (H_t, B_t)
  - Simplified manifold MALA (sMMALA) samples the continuous parameters
    (k, d, f) using gradient information
  - Adaptive Metropolis-Hastings handles (tau_ON, tau_OFF)

**Relevance:** No variable-dimension sampling at all. The "number of pulses"
emerges naturally from the ON/OFF process rather than being a discrete
parameter to sample. They report reliable convergence across healthy men,
pre/post-menopausal women, PCOS, and HA. Open-source C++ implementation
available. They cite Johnson et al. (the BD-MCMC approach libpulsatile
implements) as the prior art they build on [1].

**Trade-off:** Their model is more parsimonious (5 parameters + latent states
vs. libpulsatile's richer per-pulse parameter structure). Whether that is a
feature or a limitation depends on the scientific questions.

## 2. Compressed sensing / sparse optimization (2014, 2022)

Faghih, Dahleh, Brown et al. reformulate pulse deconvolution as a sparse
recovery problem, avoiding MCMC entirely [2]:

- Cortisol secretion is modeled as a 2nd-order linear ODE with sparse pulsatile
  inputs.
- Solved via coordinate descent: the FOCUSS algorithm (iteratively reweighted
  L1 minimization) plus generalized cross-validation for regularization.
- Assumes 15-22 secretory events over 24 hours; GCV balances sparsity against
  residual error.
- Achieves R^2 > 0.92 on real 24-hour cortisol data sampled every 10 minutes.

Amin, Faghih et al. (2022) extended this to multi-hormone sparse system
identification (leptin-cortisol dynamics) using state-space models with sparse
recovery [3].

**Relevance:** Orders of magnitude faster than MCMC, and the compressed-sensing
formulation naturally regularizes the ill-posedness seen in this model class.
The downside is point estimates only -- no full posterior uncertainty. But if
convergence trouble stems from a genuinely multi-modal or weakly identified
posterior, L1 regularization can give a cleaner answer about where the
identifiability problems lie.

## 3. Variational inference for point-process deconvolution (2020)

Shibue & Komaki use marked point processes as latent variables for calcium
imaging deconvolution (structurally analogous to hormone-pulse deconvolution)
and solve with variational inference rather than MCMC [4]. This simultaneously
estimates event times, amplitudes, and shapes. VI yields approximate posteriors
(uncertainty quantification) at MCMC-competitive quality but much faster.

## 4. HMC for birth-death-sampling models (2024)

Shao, Magee, Suchard developed scalable gradients enabling Hamiltonian Monte
Carlo for episodic birth-death-sampling models, achieving 10-200x efficiency
gains over standard MCMC [5]. The key insight is a linear-time gradient
computation that makes HMC feasible for these models.

This is phylodynamics, not endocrinology, but the mathematical structure
(birth-death process with observations) is closely related. If efficient
gradients can be derived for this model's likelihood, HMC/NUTS could
dramatically improve mixing without changing the model.

## 5. Normalizing flows to accelerate MCMC (2022)

Gabrie, Rotskoff, Vanden-Eijnden formalize Monte Carlo augmented with
normalizing flows [6]: learn a normalizing flow as an adaptive proposal
distribution for MCMC, targeting the posterior's modes. Requires limited prior
data to train. This addresses slow mixing directly -- the flow learns to
propose jumps between modes that standard random-walk proposals miss.

## 6. Posterior-based proposals (2019)

Pooley et al. introduce posterior-based proposals (PBPs) for accelerating MCMC,
finding them "significantly faster than or competitive with existing methods"
across various model types [7]. Uses information from the posterior to inform
the proposal distribution.

---

## Practical recommendations

Given libpulsatile's known pain points (slow BD-MCMC, weak identifiability,
priors constrained to achieve convergence):

1. **Highest-value investigation:** HormoneBayes's approach of modeling the
   generator process rather than the pulses themselves. This eliminates the
   transdimensional problem and its convergence difficulties. A hybrid is
   conceivable: their SMC-based latent-state approach with libpulsatile's
   richer observation model.
2. **Quick diagnostic win:** Implement the Faghih compressed-sensing approach
   as a fast initialization for BD-MCMC. Even keeping full Bayesian inference,
   good starting values from L1-regularized optimization could dramatically
   reduce burn-in.
3. **Medium-term:** Investigate whether the model admits efficient gradient
   computation (autodiff through the likelihood). If so, HMC/NUTS via Stan
   could replace the custom sampler with minimal model changes and a 10-100x
   speedup.
4. **On weak identifiability specifically:** HormoneBayes encodes physiological
   knowledge directly in the generative model (LH half-life priors, known
   ON/OFF dynamics) rather than as constraints on pulse parameters. The
   needing-constrained-priors-to-converge problem may be a symptom of the
   pulse-centric parameterization being inherently weakly identified, not an
   implementation bug.

---

## References

[1] Voliotis M, Abbara A, Prague JK, Veldhuis JD, Dhillo WS,
Tsaneva-Atanasova K. "HormoneBayes: A novel Bayesian framework for the
analysis of pulsatile hormone dynamics." PLOS Computational Biology 20,
e1011928 (2024). doi:10.1371/journal.pcbi.1011928

[2] Faghih RT, Dahleh MA, Adler GK, Klerman EB, Brown EN. "Deconvolution of
Serum Cortisol Levels by Using Compressed Sensing." PLoS ONE 9, e85204 (2014).
doi:10.1371/journal.pone.0085204

[3] Amin MR, Pednekar DD, Azgomi HF, van Wietmarschen H, Aschbacher K,
Faghih RT. "Sparse System Identification of Leptin Dynamics in Women With
Obesity." Frontiers in Endocrinology 13, 769951 (2022).
doi:10.3389/fendo.2022.769951

[4] Shibue R, Komaki F. "Deconvolution of calcium imaging data using marked
point processes." PLoS Computational Biology 16, e1007650 (2020).
doi:10.1371/journal.pcbi.1007650

[5] Shao Y, Magee AF, Vasylyeva TI, Suchard MA. "Scalable gradients enable
Hamiltonian Monte Carlo sampling for phylodynamic inference under episodic
birth-death-sampling models." PLOS Computational Biology 20, e1011640 (2024).
doi:10.1371/journal.pcbi.1011640

[6] Gabrie M, Rotskoff GM, Vanden-Eijnden E. "Adaptive Monte Carlo augmented
with normalizing flows." PNAS 119, e2109420119 (2022).
doi:10.1073/pnas.2109420119

[7] Pooley CM, Bishop SC, Doeschl-Wilson A, Marion G. "Posterior-based
proposals for speeding up Markov chain Monte Carlo." Royal Society Open
Science 6, 190619 (2019). doi:10.1098/rsos.190619
