# Proof audits — November 26–27, 2026

The audit used the actual 2,152-line July 15 source recovered from the original chat download, not a reconstruction from summaries. The main C1–C4 proof remains intact after the corrections below. The canonical Poissonized Good–Turing companion is a new, shorter exposition based on a separately recovered draft and a direct independent derivation.

## Corrections to the July manuscript

1. **Quadratic-factor regularity.** The old lemma said that division by `(x-u)^2` preserves every extra differentiability order of F. This is false: `F(x)=x^2|x|^(3/2)` is C3, but its quotient is only C1. The revised lemma explicitly assumes jointly continuous `F_txx` and the lower mixed derivatives, giving continuous G and G_t by an integral formula. Analytic F gives analytic G. Every application in this paper is analytic, so the concrete phase theorems are unchanged.
2. **Actual nonzero frozen marginal.** The abstract theorem did not exclude an identically zero frozen equation. A constant objective makes every point critical, so no active-level conclusion follows in that case. The revised theorem requires the actual frozen function to be nonzero. The fixed-sample C3 proof already establishes this using its unique highest exponential rate and positive coefficient alpha_t.
3. **Stationary arcs versus the global boundary.** The implicit function theorem at a nondegenerate stationary Maxwell tie gives an analytic local tie curve. With potentially competing ties, it does not by itself identify that curve with the global phase boundary. The proposition now asserts the local curve and positive slope; identification with beta requires the stated unique global nonuniform tie and strict uniform stability. No global analyticity claim is used in C1/C2.
4. **Truncated formal series.** The positivity remark equated an exponential truncated logarithmic series with the full product. The corrected identity is a congruence modulo z^(t+1), which gives exactly the coefficients required in C3.

Two related editorial repairs distinguish pointwise local phase orientation from a stronger neighborhood-monotonicity assertion, and state precisely the fixed-r differentiability conclusion from definable stratification.

## Claim and dependency ledger

| Claim / TeX label | Dependencies actually checked | Audit conclusion |
|---|---|---|
| `lem:one-valley` | Interior improvement at zeros; one marginal root on each monotonicity branch; integral plateau perturbation | Valid, including flat second derivatives on a strictly increasing marginal branch |
| `lem:quadratic-factor`, `lem:first-return` | Integral Taylor factor; compact maximizing subsequence; strictly positive parameter derivative near a tie | Valid with repaired mixed regularity; handles isolated pinch zeros without a nondegenerate Hessian |
| `thm:expected-threshold` | q'(z), q''(z), negative endpoint, upward interior crossing, positive spinodal derivative | Valid; no unsupported branch-uniqueness assumption |
| `prop:mgf-one-heavy` | Explicit curvature numerator; one sign crossing; boundary marginal strictly smaller than h'(0) | Valid full-support one-heavy reduction |
| `lem:full-duality` | Exponentiate h(lambda,s)-h(s,lambda)=(s-lambda)z; affine simplex constraint | Valid pointwise excess and optimizer-set duality |
| `lem:upward-transversality-new` | Negative shared marginal; h(x)>h(u)>h(y); convex bracket at a negative argument | Valid and strictly positive |
| `lem:spinodal-strong` | Unique root after factored sign threshold; d<2; exact K_lambda and K_s signs | Valid; numerical roots independently checked |
| `thm:late-ray-new`, `thm:beta-cont-new` | First return, endpoint exclusion, transversality, uniform quadratic factor | Valid; no regular-arc assumption |
| `thm:comparison-main`, Appendix A | Exact S/N normalization and decomposition; three exponential monotonicities; positive telescoping increments | Independently recomputed identities and checked each elementary sign argument |
| `thm:early-down-new` | Integrated comparison excludes A q_mu>=0 at a tie; dual-kernel order; convex bracket | Valid strictly negative s derivative |
| `thm:local-orientation`, `cor:beta-strict-new` | Finite analytic tie set; negative s derivatives; uniform factor; compact negative margins; first-return topology | Valid; includes degenerate and multiple ties |
| `thm:beta-zero-new`, `prop:spinodal-bound` | C2-uniform small-s expansion and quadratic gap; exact root equation | Valid endpoint and large-s two-sided bounds |
| `thm:C1C2-new` | Late ray + duality + strict increase + small-s endpoint | Exact stated zero-set classification follows |
| `prop:regular-maxwell-arcs` | Analytic IFT and transversality | Corrected scope: local tie curve; global beta only under additional hypotheses |
| `thm:C3-new` | Multinomial EGF; Stirling log coefficients; real characteristic roots; nonzero highest exponential rate | Valid t(t+1) bound for s!=0; smaller-span dimension alone does not improve it |
| `thm:C4-support` | Upper semicontinuity; top-t mass criterion; uniform saturated construction | Valid existential support bound and sharp saturation; no universal shape theorem |
| GT companion structural theorem | Coordinatewise unbiasedness; disjoint error events; derivative cubic; Descartes sign rule; mass-preserving second variation | Valid canonical order-r structure, including ordinary full support |

The appendix proof does not rely on an unprinted positive-coefficient certificate. Its crucial induction is `J_(k+1)-J_k = p^(-j) E_0 - theta_j E_j > 0`, using both `theta_j<1` and `p^j E_j<E_0`. These inequalities are proved independently in the manuscript, rather than inferred from tests.

For the no-reentrance proof, an interior tie is handled by integrating a positive parameter derivative backwards from a global zero. A limiting uniform tie is treated through the quadratic factor G and the positive derivative G_lambda at the spinodal. This is precisely what the earlier component-endpoint argument failed to establish.

## Independent computations

`python scripts/audit.py` passed:

- 34 symbolic identities, including all three appendix identities for k=1,...,10 and the root-reduced spinodal derivative relation.
- 900 seeded 70-digit adversarial comparison cases across logarithmic parameter scales, checking the three monotonicity inequalities, positive remainder, and conditional comparison.
- 100 fixed-sample cases comparing direct multinomial enumeration, the EGF product, truncated logarithmic coefficients, positive frozen coefficients, and a directly differentiated tangent derivative.
- 1,968 exact rational tail cases, comparing direct multinomial sums against inclusion–exclusion and the top-t plateau criterion. The explicit counterexample gives 9/25 versus 1/3 exactly.
- 24 independently computed late boundaries for n in {2,3,5,10}, with 504 factorized-excess sign checks, spinodal orientation, strict sampled boundary ordering, and the concavity/spinodal sandwich.
- 240 full-simplex SLSQP optimizations compared with the independently searched one-heavy envelope.

`python scripts/audit_good_turing.py` passed:

- Exact symbolic first and second derivatives of the general risk kernel.
- Exact positive-root counts for curvature cubics at theta=2,...,50.
- 225 high-precision coordinate bias/risk identities.
- 288 independent full-simplex optimizer runs compared with the one-exceptional families across all support sizes.

Numerical phase values and optimizer comparisons are not interval certificates. They exercise alternative formulations and search for contradictions; the general conclusions rest on the analytic proofs. Results and method parameters are saved in `checks/`.

## Historical Good–Turing gaps identified on November 26

The earlier active-uniform v13 manuscript claims a complete ordinary Poissonized risk phase diagram and a golden-ratio threshold expansion. The November 26 audit withheld these claims for two reasons.

- The stationary-branch proof asserts that a substituted numerator for b'' is positive without giving a sign certificate. Strict concavity of log(-phi') alone does not supply that omitted inequality.
- The surrogate-to-true exponential displacement argument evaluates one candidate but does not fully localize the optimal branch and control the optimized value throughout the comparison interval.

The November 27 revision supplies the missing arguments in [the companion's phase sections](good-turing/phase-transition.tex).

## November 27 Good–Turing proof audit

| Claim | Decisive check |
|---|---|
| Paired heavy coordinate, `gtp:pair` | Average the logarithmic marginal derivative across the two coordinates. Every denominator is positive. The midpoint derivative has absolute value below one and its second derivative is a strictly positive integral of a square. The chain rule gives positive heavy-coordinate curvature. |
| Fold, stability and coexistence, `gtp:phase` | The branch slopes run from negative infinity to n−2. At the fold the scalar cubic is nonzero. Light-splitting directions have negative curvature, so the historical “local minimum” is a full-simplex saddle for n≥3. The excess strictly increases with intensity on both branches. Integration from the saddle to uniform gives a strictly negative fold excess. Uniform itself is a saddle at the upper endpoint for n≥3. |
| Degenerate cases | A singleton alphabet has one feasible law. The binary transition coalesces at twice the inflection intensity, equal to 6 at order zero. At the higher-order boundary t = n times the kernel peak, a nonuniform candidate would have negative marginal, exclude all zeros, and force every coordinate above the peak, contradicting the total. |
| Localization, `gtp:window` | Test the branch equation at the kernel peak plus 1/n. Its partner is below 2 log n for large n. This brackets the maximizing stationary point throughout the critical window and bounds the gain from optimizing the light coordinates by the square of the exponential tail times a polynomial. |
| Displacement and expansion, `gtp:shift`, `gtp:gold` | The window contains both roots. The mean-value denominator is of order n to the power −1/2 throughout it. A first exponential bound permits evaluation at the upper surrogate root with exponentially small relative error. Analytic inversion supplies the algebraic remainder. The archived coefficients agree through order n to the power −1/2, and the companion includes the next coefficient and the leading exponential constant. |
| Higher-order scope, `gtp:light` | Above the full-support peak the same positive denominators apply at every order. Below it, a rational inequality forces a light singleton for n=3 and r=10 on a right neighborhood of twice the peak. A general classification of the earlier phases remains open. |

The full `make check` passed with the pinned dependencies. The [new checker](scripts/audit_good_turing_phase.py) verifies the ordinary Taylor identities through degree five and the rational bound for the higher-order example. It compares the integral certificate with direct implicit differentiation at 20 paired roots for occupancy orders 0, 1, 2, 6, and 10. An explicit quadratic parametrization of the ordinary pair independently yields the fold and global crossing for n = 3, 4, 10, 100, 1,000, and 10,000. For n = 3 these are 8.62088142129 and 8.66460008083. For n = 10,000 the displacement is approximately 1.30185 × 10⁻⁷⁹, about 1.03116 times its leading term. [Recorded phase checks](checks/good-turing-phase.json).

The phase conclusions are analytic. The optimization concerns the canonical estimator with independent Poisson counts. General fixed-count extremizers, arbitrary biased estimators, and exact minimax estimation remain open here.

## Literature and release limits

Primary sources checked include Berend–Kontorovich, *The Missing Mass Problem* (arXiv:1111.2328); McAllester–Ortiz's JMLR concentration paper; and Skorski, *On Missing Mass Variance* (arXiv:2104.07028). They motivate careful separation of expected mass, positive MGFs, variance, and estimator risk. The current main manuscript explicitly cites the fixed-sample expected-mass antecedent. A targeted search did not identify a competing positive-MGF phase-diagram result, but an exhaustive priority check has not been established. This audit is not external refereeing or formal proof verification.

The Good–Turing companion additionally cites [Skorski's fixed-sample asymptotic MSE analysis](https://arxiv.org/abs/2104.07029), whose formula contains a negative squared coupling term and whose asymptotic extremizer is Dirac–uniform, and [Painsky's general-order convergence guarantees](https://www.jmlr.org/papers/v23/21-1528.html). The companion explicitly distinguishes its exact independent-Poisson risk from these antecedents.

On November 26 both PDFs compiled without LaTeX warnings after three passes. All 35 main-paper pages and the four companion pages of that version were rendered and visually inspected. On November 27 both manuscripts were rebuilt and the revised nine-page companion was rendered for inspection. The complete Python checks also passed in a second fresh CPython 3.13.18 environment using the pinned requirements.
