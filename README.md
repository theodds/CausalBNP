# CausalBNP

 Code for reproducing examples in the book by Daniels, Linero, and Roy.

# Natural Planning

## Why?

__Goal:__ add code here so that folks can easily reproduce at-least-some of the
material from the case studies. Focus on examples from the case studies, for
simplicity.

## Done?

Ability to replicate all of the results from the various chapters, with some
description of how to do everything in the README. Possibly implemented in a
package. The methods we would need to be able to address are:

1. Chapter 8: BART + DP causal inference on quantiles
   - Package for doing this
   - Illustration: showing how to get the density
   - Intervals for quantiles, plus density plot
   
1. Chapter 9: Causal Inference with EDPM model [TODO]

1. Chapter 10: DDP + GP for causal inference with marginal structural models [TODO]

1. Chapter 11: DPMs for dropout:
   - Simulated data
   - Input: `Y`, with  missing data, as well as treatment vector `A`
   - Output: Intervals for each time point; goodness of fit checks

1. Chapter 12: DPMs for no dropout
   - Simulated data
   - Input: `Y` with missing data, as well a treatment vector `A`
   
1. Chapter 13: Causal Mediation Using DPMs [TODO]

1. Chapter 14: Causal Mediation Using BART

1. Chapter 15: Semicompeting Risks [TODO]


## Brainstorm

- How do we want to organize the code?

  - Folders for each chapter?
  - Single package with functions? Probably this is the way, and then illustrate
    the functions using different scripts in the folder
  - But what dependencies? And probably we should push this stuff up to CRAN?

- What packages do I need? Is there a way to avoid having to install too many
  packages up front? I guess we will need a package to do this stuff?

## Next Actions

- [ ] Run the MEPS example (QTE)
- [ ] Plot the MEPS example (QTE)
