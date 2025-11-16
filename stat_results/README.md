Statistical Analysis of Interaction Performance Metrics (INP/LoAF) — inp_stat

This repository contains the Python implementation of the statistical analysis pipeline used in the study “An Inter-Metric Approach to Interaction Performance: Combining INP, LoAF and Event Timing for Practical Frontend Engineering.”
It complements the experimental INP-first measurement environment (interaction scenarios, optimization variants, and raw telemetry) described in the main article.

The goal of this project is to validate whether specific frontend optimization techniques provide statistically significant reductions in Interaction to Next Paint (INP) and to quantify the magnitude of these effects across multiple scenarios and CPU-throttling levels.

The analysis follows a reproducible, framework-agnostic methodology, based on raw measurements exported from the INP-first experimental runner (INP, Event Timing phases, LoAF frame costs).

Overview

Modern web interfaces exhibit complex event-processing pipelines where responsiveness depends not only on JavaScript execution speed but also on queueing delays, DOM mutation patterns, and rendering-pipeline pressure.
The companion experiment (conducted separately from this repository) measures INP and LoAF across three representative UI scenarios:

Content — light DOM mutations, accordion interactions

Dashboard — large lists (≈3000 items), partial re-renders, sorting/filtering

Form — synchronous validation logic, heavy processing tasks

Each scenario is tested across multiple optimization strategies (I1–I5) and CPU-throttling levels (×2, ×4, ×6, ×8).
This repository performs post-hoc statistical analysis of those measurements.

Objectives

The statistical pipeline provides answers to the following research questions:

Do optimization variants statistically reduce INP relative to the baseline B0?

How strong is each effect, and is it stable across scenarios and CPU levels?

Do LoAF-related improvements support/corroborate the observed INP effects?

Which optimizations are reliable across conditions (architecture-level) vs. conditional/local (scenario-specific)?

Methods
1. Non-parametric Hypothesis Testing

INP values exhibit asymmetric, right-tailed distributions with occasional outliers 

stat

.
For each pair (Baseline B0 ↔ Variant V) within a fixed (scenario × CPU level) group, we test:

H₀: “INP distributions do not differ in the direction of reduction.”

H₁: “Variant V produces smaller INP than baseline B0.”

The pipeline uses:

scipy.stats.mannwhitneyu(sample_V, sample_B0, alternative="less")


Each group contains 15 independent INP observations per variant, following the experimental protocol.

2. Multiple-Testing Correction

Since each scenario has several optimization variants, the analysis applies:

Holm–Bonferroni correction
to control family-wise Type I error (α = 0.05).

3. Effect Size Estimation

For every B0–V comparison, two complementary effect size metrics are calculated:

AUC (Area Under ROC Curve)
Interpretable as the probability that a random INP sample from V is smaller than from B0.
Values ≈0.7–0.8 indicate medium–large effects.

Cliff’s Delta (δ)
Ranges from −1 to 1.
Positive values mean Variant V tends to reduce INP.
|δ| close to 1 indicates dominant, practically meaningful effects.

Both metrics are implemented in this repository according to the definitions in the study.

4. Output

The pipeline generates:

per-scenario statistical tables

adjusted p-values (Holm–Bonferroni)

effect size metrics (AUC, Cliff’s delta)

optional plots (ECDF, boxplots)

aggregation summaries similar to Tables 1–2 in the statistical report 

stat

Key Findings (Summary)

Statistical results strongly reinforce the qualitative patterns observed in the raw INP/LoAF telemetry:

1. Architectural optimizations containing I2 produce large, robust, and universal improvements.

Significant for all scenarios (Content, Dashboard, Form)

Significant for all CPU levels (except one outlier in Content at ×2)

Cliff’s delta ≈ 1.0 in most conditions — meaning nearly every measurement from the optimized variant is smaller than every measurement of B0

p-values after Holm–Bonferroni correction reach 10⁻⁵–10⁻⁶

2. Local optimizations without I2 show inconsistent or weak statistical effects.

I1, I3, I4, I3+I4:

occasionally statistically significant

but not stable across CPU conditions

effect sizes smaller and noisy

improvements often too small to be practically meaningful in UX terms (INP remains high)

3. I5 (CSS containment) is beneficial but scenario-dependent.

Meaningful improvements only in Dashboard

Best interpreted as a supporting optimization, not a primary remediation

4. Extreme improvements in the “form” scenario confirm the importance of reducing processing (D₂) load.

Median INP drops from ~2.2 s to ~24 ms under variants containing I2

Distributions barely overlap — Cliff’s delta = 1.0

Strongest and clearest statistical effect in the whole study

These patterns validate the broader inter-metric model from the main paper:

INP ≈ D₁ (input delay) + D₂ (processing) + D₃ (presentation) 

article


and confirm that optimizations targeting processing and queueing (I2/I3) are the most impactful across real-world interaction types.

Repository Structure
.
├── stats_pipeline.py        # Main analysis script
├── requirements.txt         # Python dependencies
├── results/                 # Raw INP/LoAF measurements from the experiment
├── stat_results/            # Generated statistical tables and plots
└── .gitignore

Installation

Clone and install:

git clone https://github.com/olehpapka28/inp_stat.git
cd inp_stat
pip install -r requirements.txt

Usage

Run the full statistical pipeline:

python stats_pipeline.py


This will:

load experiment results from results/

compute Mann–Whitney U tests + Holm–Bonferroni correction

calculate effect sizes

generate summary tables in stat_results/

optionally produce diagnostic plots

Data Requirements

The pipeline expects experiment outputs in the structure generated by the INP-first runner, including:

INP values per interaction

LoAF frame data

scenario labels

variant labels (B0, I1, I2, I2+I3+I4, …)

CPU throttling level

Each (scenario × variant × CPU) group must contain 15 replications, as defined in the experimental protocol.

Reproducibility Notes

The analysis uses only non-parametric tests due to non-Gaussian distribution shapes.

The environment is fully deterministic and does not rely on random seeds.

Holm–Bonferroni ensures valid family-wise error control.

All effect-size computations follow standard definitions reported in the paper.

Related Work

This repository is part of a larger research effort combining:

Interaction to Next Paint (INP)

Event Timing (input delay + processing phases)

Long Animation Frames (LoAF)
into a unified diagnostic framework for identifying the root causes of interaction latency in modern web interfaces.

For background and motivation, see the accompanying article and experimental methodology.
(The measurement environment and INP experiment runner are hosted in a different repository.)

License

MIT License (or specify if different).

Contact

For research questions or collaboration:

Authors: R. Oliinyk, O. Papka

Institution: Lviv University of Trade and Economics

Project: Inter-metric analysis of interaction responsiveness (INP/LoAF/Event Timing)
