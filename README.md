# A Mathematical Framework for Community Testing Strategies 
Correspondance: casey.middleton@colorado.edu

## Purpose
This model combines viral kinetics and probability theory to predict how a given testing regimen reduces the effective reproductive number, R0. 
This repository contains the files needed to walk through an example of model capabilities, and replicate author results. 

## A worked example 
Examplevirus is currently taking the world by storm, pressuring public health officials to determine the best strategies for containment. 
There is no vaccine and no treatment for examplevirus, but perhaps *community screening programs* could be used to detect and quarantine infections early, preventing further transmission.
It is still early in the epidemic, but public health officials know a few things about examplevirus:
1. On average, viral growth and decay each take 5 days, giving a total infection duration of 10 days. 
2. Viral growth and decay is exponential, peaking on day 5 at approximately Log10(5) copies of examplevirus / mL.
3. An individual's infectiousness, or the likelihood of transmission given a contact, is proportional to their viral load at a given time. 

There are two viral detection tests currently available on the market:
1. High Sensitivity Viropath (HSV) can identify the virus immediately, at Log10(0) copies of examplevirus / mL, but has a high cost per test, requires invasive testing strategies, and may take 1-3 days to process results.
2. Low Sensitivity Viropath (LSV) is slightly less sensitive, requiring Log10(2) copies of examplevirus / mL to detect the virus, but is 1/3 the price of HSV, less invasive, and provides immediate results. 
Which test should public health officals recommend for community testing strategies? How often should these tests be administered, and what outcome is expected under this testing procedure?

Since examplevirus infectiousness is proportional to viral load, consider the area under the viral load curve a measure of total infectiousness over the course of a given infection. This is the area shown in gray below. We can then calculate the first detectable instance under an LSV testing program. If tests are admnistered once a week, under a uniform distribution of testing, the cumulative probability of detection grows linearly from the first detectable instance, eventually reaching guaranteed detection. Intuitively, the infectiousness that occurs *before* the first detectable instance has no chance of being eliminated through isolation of the infected individual, while the infectiousness that occurs *after* guaranteed detection is certainly eliminated through isolation (assuming perfect isolation given a positive test result). What happens between these two points, E(infectiousness | testing), can be predicted by multiplying infectiousness without testing by the cumulative probability of detection. The expected infectiousness with testing can then be calculated by the pink hatched area. Finally, the reduction in infectiousness can be linked to the reduction in transmission using the ratio of the pink area to the gray area. 
![alt text](https://github.com/CaseyMiddleton/TestingFramework/blob/main/example_figs/example_model_diagram_2.png)



