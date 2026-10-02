# Bounded source check — priority HOLD

2026-10-02. Это карта зависимостей и ограниченный поиск, не exhaustive audit.

- McGinnis–Li–Mori, [2606.15443v1](https://arxiv.org/html/2606.15443v1),
  §4 Eq. (23), §§5,7: weighted CL error, cubic remainder, active coercivity,
  local convergence near a regular solution. Read from previously retained
  exact version text; official HTML reopened this pass. These tools are baseline.
- VCS core v0.2.1, Thm 6.1: global passive trace, 2-Lipschitz bound in l1
  physical conductance coordinates, including hidden rank changes. Local exact
  dependency copy is unchanged. The dependency is a theorem, not an old ACCEPT.
- Alman–Lian–Tran, [1309.3011](https://arxiv.org/abs/1309.3011): circular
  planar response positivity. Search abstract only this pass. No planar
  characterization is assumed in the new proof.
- Two recent **discovery candidates only**, not read/validated as full baselines:
  [2609.01672](https://arxiv.org/abs/2609.01672), arbitrary superport minor
  formulas; [2609.14385](https://arxiv.org/abs/2609.14385), cohomological
  electrical-network response minors. Their surfaced abstracts concern minor
  representations, not a checked physical-learning lower bound. Do not infer
  novelty or absence of coverage from abstracts.

Queries: `resistor network learning Jacobian response matrix minors Dirichlet
form coercivity positive parts`; `electrical network response matrix positive
negative parts Dirichlet principle nonplanar minors`.

Potential independent remainder to audit: quantitative physical sensitivity
over all realizations of measured ports with bounded l2 norm, at collided
outputs, plus an explicit original projected CL corridor. No claim that
Dirichlet clipping, Cauchy–Schwarz, response minors or bootstrap are new.
