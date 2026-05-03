# Varying parameters
A subanalysis of the main article in which we run sensitivity analyses
to see how the predictive accuracy of the fALS criteria is altered by using
a different range of input parameters

## Setup
The sensitivity analyses consist of 12 analyses in which we keep all parameters 
the same as the main analysis, except for one specific input parameter which we
vary. 

In the main analysis we conducted one run of 1.000.000 pedigrees. In this subanalysis 
we conduct 10 runs of (at least) 100.000 pedigrees. If necessary, we increase 
the number of pedigrees simulated, to keep the number of index patients equal across 
the varying range of a parameter. For example, we run more pedigree simulations when the
disease allele pentrance is 5% then when it is 20%, to keep the same number of index
patients to analyse the fALS criteria. 

## Varying paramater values
### Demographic parameters
- age_onset_monogenic
    - We run the analyses as if the ALS-FTD common disease allele had a mean onset
    10 years earlier/later then in the literature estimates
    - In the main analysis, the monogenic ALS-FTD has an earlier onset than polygenic disease
- age_onset_polygenic
    - We run the analyses as if the polygenic ALS had a mean onset 10 years earlier/
    later then in the literature estimates
    - In the main analysis, the polygenic ALS has a later onset than polygenic disease
    - See: [historical_disease_onset_9x3.pdf](../rationale/disease_onset/historical_disease_onset_9x3.pdf)

- fert_rate
    - We run the analyses with a set fertility rate for all individuals, ranging 1-5 offspring per mother
    - In the main analysis, the fertility rate was year-of-offspring dependent. Historically, 
    people used to have more offspring in the past then in recent years, resulting in bigger pedigrees
    in the past
    - See: [simulated_fertility_rate.pdf](../rationale/fertility_rate/simulated_fertility_rate.pdf)

- life_exp
    - We run the analyses with a set life expectancy for all individuals, ranging 50 to 100 years for all individuals
    - In the main analysis, the life expectancy was year-of-birth dependent. Historically, 
    the life expectancy has increased in recent years, giving individuals who were to be disease-affected more 'chance' to actually develop the phenotype as they were to get older. 
    - See: [simulated_life_expectancy_comparison.pdf](../rationale/life_expectancy/simulated_life_expectancy_comparison.pdf)

## Rest of varying params will follow 