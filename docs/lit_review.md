# Literature review: alternatives and extensions to BD-MCMC deconvolution

**Date:** 2026-06-07  
**Status:** Background research. Nothing here is scheduled work; it is the
evidence base behind the future-appendix entries (I-O) in
`2026-06-07-vignette-roadmap-design.md`.

---

Here's a synthesis of what I found across 8M+ biomedical papers. Several substantive developments since your
  ~2015-2017 work are worth investigating.

  ---
  1. HormoneBayes: SMC replaces BD-MCMC entirely (2024)

  This is the most directly relevant paper. Voliotis, Abbara, Prague, Veldhuis, et al. sidestep the
  transdimensional sampling problem by reformulating the generative model [1]:

  - Instead of estimating the number/location of pulses (requiring BD-MCMC), they model the underlying
  hypothalamic ON/OFF switching process as a two-state Markov chain with transition rates tau_ON and tau_OFF
  - The observed LH is driven by both a pulsatile signal P_t and a basal signal B_t
  - Inference uses a Gibbs sampler where:
    - Sequential Monte Carlo (SMC) with ancestral sampling handles the latent state trajectories (H_t, B_t)
    - Simplified manifold MALA (sMMALA) samples the continuous parameters (k, d, f) using gradient information
    - Adaptive Metropolis-Hastings for (tau_ON, tau_OFF)

  Why this matters for you: No variable-dimension sampling at all. The "number of pulses" emerges naturally from
   the ON/OFF process rather than being a discrete parameter to sample. They report reliable convergence across
  healthy men, pre/post-menopausal women, PCOS, and HA. Open-source C++ implementation available. They cite
  Johnson et al. (your BD-MCMC approach) as the prior art they're building on [1].

  The tradeoff: their model is more parsimonious (5 parameters + latent states vs. your richer per-pulse
  parameter structure). Whether that's a feature or a bug depends on your scientific questions.

  ---
  2. Compressed Sensing / Sparse Optimization (2014, 2022)

  Faghih, Dahleh, Brown et al. reformulate pulse deconvolution as a sparse recovery problem, avoiding MCMC
  entirely [2]:

  - Model cortisol secretion as a 2nd-order linear ODE with sparse pulsatile inputs
  - Solve via coordinate descent: FOCUSS algorithm (iteratively reweighted L1 minimization) + generalized
  cross-validation for regularization
  - Assumes 15-22 secretory events over 24 hours; GCV balances sparsity vs. residual error
  - Achieves R^2 > 0.92 on real 24-hour cortisol data sampled every 10 minutes

  Amin, Faghih et al. (2022) extended this to multi-hormone sparse system identification (leptin-cortisol
  dynamics) using state-space models with sparse recovery [3].

  Why this matters for you: Orders of magnitude faster than MCMC. The compressed sensing formulation naturally
  regularizes the ill-posedness you're experiencing. The downside: point estimates only -- no full posterior
  uncertainty. But if your convergence issues stem from the posterior being genuinely multi-modal or poorly
  identified, the L1 regularization might actually give you a cleaner answer about where the identifiability
  problems lie.

  ---
  3. Variational Inference for Point Process Deconvolution (2020)

  Shibue & Komaki use marked point processes as latent variables for calcium imaging deconvolution (structurally
   analogous to your problem) and solve with variational inference rather than MCMC [4]. This simultaneously
  estimates event times, amplitudes, and shapes. VI gives you approximate posteriors (uncertainty
  quantification) at MCMC-competitive quality but much faster.

  ---
  4. HMC for Birth-Death-Sampling Models (2024)

  Shao, Magee, Suchard developed scalable gradients enabling Hamiltonian Monte Carlo for episodic
  birth-death-sampling models, achieving 10-200x efficiency gains over standard MCMC [5]. Their key insight: a
  linear-time gradient computation algorithm makes HMC feasible for these models.

  This is for phylodynamics, not endocrinology, but the mathematical structure (birth-death process with
  observations) is closely related. If you could derive efficient gradients for your model's likelihood,
  HMC/NUTS could dramatically improve mixing without changing your model.

  ---
  5. Normalizing Flows to Accelerate MCMC (2022)

  Gabrié, Rotskoff, Vanden-Eijnden formalize Monte Carlo augmented with normalizing flows [6]. The idea: learn a
   normalizing flow as an adaptive proposal distribution for MCMC, targeting the modes of your posterior.
  Requires limited prior data to train. This could address slow mixing directly -- the flow learns to propose
  jumps between modes that standard random-walk proposals miss.

  ---
  6. Posterior-Based Proposals (2019)

  Pooley et al. introduce posterior-based proposals (PBPs) for accelerating MCMC, finding them "significantly
  faster than or competitive with existing methods" across various model types [7]. Uses information from the
  posterior to inform the proposal distribution.

  ---
  Practical Recommendations

  Given your specific pain points (slow BD-MCMC, under-specification, constrained priors needed):

  1. Highest-value investigation: HormoneBayes's approach of modeling the generator process rather than the
  pulses themselves. This eliminates the transdimensional problem and the associated convergence difficulties.
  You could potentially hybridize: use their SMC-based latent state approach while keeping your richer
  observation model.
  2. Quick diagnostic win: Implement the Faghih compressed sensing approach as a fast initialization for your
  BD-MCMC. Even if you want full Bayesian inference, getting good starting values from L1-regularized
  optimization could dramatically reduce burn-in.
  3. Medium-term: Investigate whether your model admits efficient gradient computation (autodiff through the
  likelihood). If so, HMC/NUTS via Stan could replace your custom sampler with minimal model changes and give
  you 10-100x speedup.
  4. For the under-specification issue specifically: The HormoneBayes paper's approach of encoding physiological
   knowledge directly in the generative model (LH half-life priors, known ON/OFF dynamics) rather than as
  constraints on pulse parameters may be the conceptual shift that helps. Your constrained-prior-to-converge
  problem might be a symptom of the pulse-centric parameterization being inherently weakly identified, not a
  coding bug.

  ---
  References

  [1] Voliotis M, Abbara A, Prague JK, Veldhuis JD, Dhillo WS, Tsaneva-Atanasova K. "HormoneBayes: A novel
  Bayesian framework for the analysis of pulsatile hormone dynamics." PLOS Computational Biology 20, e1011928
  (2024). doi:10.1371/journal.pcbi.1011928
      https://citations.gxl.ai/papers/PMC10931492#L39,L44,L45

  [2] Faghih RT, Dahleh MA, Adler GK, Klerman EB, Brown EN. "Deconvolution of Serum Cortisol Levels by Using
  Compressed Sensing." PLoS ONE 9, e85204 (2014). doi:10.1371/journal.pone.0085204
      https://citations.gxl.ai/papers/PMC3904842#L12,L16,L18

  [3] Amin MR, Pednekar DD, Azgomi HF, van Wietmarschen H, Aschbacher K, Faghih RT. "Sparse System
  Identification of Leptin Dynamics in Women With Obesity." Frontiers in Endocrinology 13, 769951 (2022).
  doi:10.3389/fendo.2022.769951
      https://citations.gxl.ai/papers/PMC9037068#L6,L8

  [4] Shibue R, Komaki F. "Deconvolution of calcium imaging data using marked point processes." PLoS
  Computational Biology 16, e1007650 (2020). doi:10.1371/journal.pcbi.1007650
      https://citations.gxl.ai/papers/PMC7093033#L0,L8

  [5] Shao Y, Magee AF, Vasylyeva TI, Suchard MA. "Scalable gradients enable Hamiltonian Monte Carlo sampling
  for phylodynamic inference under episodic birth-death-sampling models." PLOS Computational Biology 20,
  e1011640 (2024). doi:10.1371/journal.pcbi.1011640
      https://citations.gxl.ai/papers/PMC11006205#L0

  [6] Gabrié M, Rotskoff GM, Vanden-Eijnden E. "Adaptive Monte Carlo augmented with normalizing flows." PNAS
  119, e2109420119 (2022). doi:10.1073/pnas.2109420119
      https://citations.gxl.ai/papers/PMC8915891#L0,L7

  [7] Pooley CM, Bishop SC, Doeschl-Wilson A, Marion G. "Posterior-based proposals for speeding up Markov chain
  Monte Carlo." Royal Society Open Science 6, 190619 (2019). doi:10.1098/rsos.190619
      https://citations.gxl.ai/papers/PMC6894579#L0
