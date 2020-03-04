# screening_outbreak_delay

Interventions targeting air travellers early in the pandemic may delay local outbreaks of SARS-CoV-2. 

| Authors |
| :-- |
| Dr. Sam Clifford |
| Dr. Carl AB Pearson |
| Assist. Prof. Petra Klepac |
| Mr Kevin van Zandvoort |
| Mr. Billy Quilty |
| Assoc. Prof. Rosalind Eggo |
| Assoc. Prof. Stefan Flasche |
| Other members of CMMID at LSHTM |

Interventions for infected travellers may delay or prevent an outbreak of SARS-CoV-2 in the destination country, such as

1. Exit screening: looking for symptoms (e.g. cough, fever) at the start airport 
2. Entry screening: looking for symptoms (e.g. cough, fever) at the destination airport
3. Sensitisation: e.g. instructions for travellers on what to do if they become ill

Given the scenario set out in the panel, this simulation calculates the likely delay to an outbreak. We present the minimum length of the delay that is seen in 97.5%, 75%, 50%, 25% and 2.5% of simulations. In cases where the plot is a line at or close to 0%, this indicates that the delays are less than one day in length. As traveller sensitisation becomes more effective, the simulate percentage of outbreaks completely delayed (i.e. averted) increases.

Here we make the simplifying assumption that the arrival rate is constant over time, which **may not be appropriate as an outbreak progresses**. As the outbreak continues, it may be more useful to account for the transmission of the disease not just at its source but in locations that it spreads to. We are developing new mathematical models to more accurately reflect that the number of infected travellers will likely increase over time in the absence of any travel restrictions.

These results are based on a stochastic simulation, which may take a few seconds to run, and there may be small variations in the calculated delay when the same parameters are used twice.

This repository contains the code for a [medR&chi;iv preprint](https://www.medrxiv.org/content/10.1101/2020.02.12.20022426v1) looking at whether or not syndromic screening and traveller sensitisation can delay an outbreak. The app is deployed at http://cmmid-lshtm.shinyapps.io/screening_outbreak_delay  The code is provided as-is, in order to demonstrate the application of the model. As the article is currently in pre-print status, it has not yet completed the peer review process and may be updated on the recommendations of reviewers.
