# Module 1: Tau Phosphorylation Dynamics — A Systems Biology ODE Model

**What this is:** a dynamical model (not a data pipeline) of the kinase–phosphatase
tug-of-war that sets tau's phosphorylation state. This is the actual distinction
between computational biology and bioinformatics — we are *simulating a mechanism*,
not processing a dataset.

## The biology being modeled
- 4 kinases ADD phosphate to tau: **GSK3B, CDK5, CSNK2A1, AKT1**
- 1 phosphatase REMOVES it: **PP2A** (gene *PPP2CA*)
- **PIN1** facilitates PP2A's action (isomerizes the phospho-Thr-Pro motif so PP2A can act)

## The equations
State variable: `Tp` = phosphorylated tau fraction (0–1). `Tu = T_TOTAL − Tp`.

```
R_phos(Tu)   = Σ_k  Vmax_k · Tu / (Km_k + Tu)              [sum over 4 kinases]
f_PIN1       = PIN1_activity / (PIN1_activity + K_pin1)
R_dephos(Tp) = Vmax_PP2A · f_PIN1 · Tp / (Km_PP2A + Tp)

dTp/dt = R_phos(Tu) − R_dephos(Tp)
```

## Scenarios simulated (all literature-motivated)
| Scenario | Change | Why |
|---|---|---|
| Healthy (baseline) | default params | reference state |
| PP2A reduced 50% | Vmax_PP2A × 0.5 | PP2A activity is reduced in AD brain (Gong et al., J Neurochem) |
| PP2A + PIN1 reduced | + PIN1_activity × 0.5 | PIN1 is depleted in AD-vulnerable neurons (Lu et al., 1999, *Nature*) |
| GSK3B hyperactive | Vmax_GSK3B × 1.5 | GSK3B hyperactivation is a long-standing AD hypothesis (Hooper et al., 2008) |

## Honest caveat on parameters
Vmax/Km values are **illustrative**, chosen for qualitatively sensible dynamics —
not extracted from kinetic assays. This is a proof-of-concept qualitative model.
State this plainly if you present it. Treat the *direction* of the result
(perturbations → higher steady-state phosphorylated tau) as the finding, not the
exact numbers.

## Files
- `tau_ode_model.py` — full model, runs standalone or in Colab (just paste and run)
- `tau_phosphorylation_dynamics.png` — time-course plot, 4 scenarios
- `tau_steady_state_comparison.png` — bar chart of steady-state Tp per scenario

## How to run
```
pip install scipy numpy matplotlib --break-system-packages   # if needed
python3 tau_ode_model.py
```
Or paste the whole script into a Google Colab cell and run.

## Results from this run
| Scenario | Steady-state Tp |
|---|---|
| Healthy (baseline) | 0.674 |
| PP2A reduced 50% | 0.880 |
| PP2A + PIN1 reduced | 0.917 |
| GSK3B hyperactive | 0.769 |

All three AD-like perturbations push the system toward a higher phosphorylated-tau
steady state than healthy baseline — qualitatively consistent with hyperphosphorylated
tau accumulation in AD.

## Next steps (Module 1b — optional upgrade)
Right now the kinase/phosphatase Vmax values are guesses. You can make this
*quantitatively grounded* by scaling each enzyme's Vmax by that gene's relative
expression level in your own GSE5281 / GSE44770 / GSE33000 / GSE36980 data —
this turns your existing omics analysis directly into model parameters, and is
a genuinely novel contribution (static network → dynamic, data-informed model).

## LinkedIn-postable framing
"Built a systems biology model simulating how literature-documented enzyme changes
in Alzheimer's (reduced PP2A, depleted PIN1, hyperactive GSK3B) shift tau toward a
hyperphosphorylated steady state — turning a static kinase-phosphatase network map
into a predictive dynamical model."
