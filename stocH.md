# Quantifying the stochastic component of epigenetic aging  

## Introduction

Epigenetic clocks are multivariate linear models that can predict chronological and biological age, 
based on DNA methylation (DNAm) data {cite}`horvathDNAMethylationAge2013`. 
Age-associated DNAm changes are observed in two phenomena – first, the evidence of specific 
age-associated DNAm changes, and second, the evidence of "erosion" or increasing uniformity 
of DNAm landscape, associated with age. 
Thus, two phenomena can be considered deterministic and stochastic components of quasi-stochastic 
epigenetic changes associated with age. 
In this regard, understanding the contribution of stochastic component to the 
accuracy of linear models (clocks) prediction represents a fundamental interest. 

In 2024 H.Tong *et al* proposed a method for quantification of stochastic component of epigenetic 
aging, and estimated that approximately 66–75% of the accuracy underpinning Horvath’s clock could 
be driven by a stochastic process
{cite}`tongQuantifyingStochasticComponent2024`. The goal of this project was to check this statement by reproducing the analysis proprosed by authors. 


## Results
### Modelling single-cell one site DNAm state

For single site in single cell one can model methylation status as a Markov chain. 
Let $X\in \{0, 1\}$ be random variable, which is equal to 1, if site is methylated and 0, if it is not.
Considering also fixed $p$ as a probability of methylation of demethylated site and fixed $q$ 
as a probability of demethylation of methylated site, one can define a homogenous Markov chain of 
methylation dynamics ({numref}`1_dnam_dyn`, A). This process is homogenous in a sence of independence
of the transition probabilities from time. For systems of such kind, the stationary state can be reached, 
as shown on ({numref}`1_dnam_dyn`, B), and the independent probaility of each state does not depend on 
initial state {cite}`shiryaevProbability12016`. However, for more shorter times dynamics can be explored.

```{figure} figures/1_dnam_dyn_with_model.svg
:name: 1_dnam_dyn
Modelling single-cell methylation dynamics of one CpG.
**A**, Markov chain process and it's schematic representation.
**B**, Dynamics of methylation $(X=1)$ state probability within a model with different parameters.
```

### Constructing and implementing stochastic model for cell population

#### Procedure description

When switching from the one-cell to cell population, model is needed to generalize. 
To do so, I will provide several definitions, which authors used in the paper.

First, instead of one CpG site, several (353) ones are considered, which belong to Horvath's clock. 
For each of CpG $c$, independent procedure of estimaion of methylation change applied. 
Second, instead of one boolinan variable, $X \in \{1, 0\}$ defining methylation status in single cell, 
I operated with the fraction of methylated sites in particular position, named DNAm  ({numref}`dyn_pop`, 
A).

Thus, taking into account the described differences, a simulation of the stochastic will consist of the 
following steps.

##### 1. Initial methylation status distribution.  
The initial methylation status corresponds to fraction of methylated sites in particular position at the 
beginning of observation (DNAm). Following the paper, I defined the initial methylation status for each CpG 
$c$ as an average DNAm for young samples (avgDNAm(Young)). 

##### 2. Definition of probability of transition  
I modelled the change of methylation status at CpG $c$ occuring with probability $p_c$:

$$p_c = 1 - e^{-\gamma |\text{EffSize}_c|} $$
where $\gamma$ is global probability of a DNAm change (discussed below), and $\text{EffSize}_{c}$ is an 
effect size, which defined as difference in methylation fraction DNAm for betweem old and young  ({numref}`dyn_pop`, B). 
Mathematically, for each CpG $c$ in a CpG set:
$$ \text{EffSize}_c = \text{avgDNAm(Old)} - \text{avgDNAm(Young)}$$

##### 3. Definition of transition
3.1. Since switching to cell population, the transition process should be defined.
Since we have a population, the change in methylation status of CpG $c$ will change of DNAm 
from $\beta _c ^{(t)}$ at time step $t$ to $\beta _c ^{(t+1)}$. The DNAm values are defined 
as $\beta_c$ since they are beta-distributed ({numref}`dyn_pop`, left on C).  
3.2. The size of DNAm change is defined as a random deviation $r_{c}$ from truncated normal 
distribution $r_{c} \sim \mathcal{N}_{+}(0, \sigma)$. Sign of change is same as siogn of 
effect size: $\text{sign}(\text{EffSize}_{c})$.  
3.3. Change in DNAm. Since the DNAms $\beta _c ^{(t)}$ are beta-distributed and the $r$ belongs
to truncated normal distribution I first transform the DNAm to normal distribution, change the
value, and then convert it back to the beta-values. Thus,
$$\beta _c ^{(t)} \quad \rightarrow \quad x_c^{(t)} = iF(\beta_{c}^{(t)}),$$ 
where $iF$ is inverse of the normal cumulative distribution function. Next,
$$ x_c^{(t+1)} = x_c^{(t)} + \text{sign}( \text{EffSize}_{c} ) \cdot r_{c}, \ 
r_{c} \sim \mathcal{N}_{+}(0, \sigma)$$
and, finally,
$$\beta _c ^{(t+1)} = F(x_{c}^{(t+1)})$$ 
The equation above describes changes of DNAm for each CpG $c$ in a case of need to change.  
Strictly speaking, I am not sure, whether the defined process is Markov one, since the set of states is 
not countable anymore. But no properties of Makov process used in following analysis, so it does not matter.
Thus, the resulting model ({numref}`dyn_pop`, A) requires determination of experimental data for 
definition of initial states and establishment of 2 parameters ($\gamma$ and $\sigma$), which control model dynamics.  

```{figure} figures/2_simulation_agg.svg
:name: dyn_pop
Modelling methylation dynamics of independent CpGs in cell population.
**A**, Transition process scheme.
**B**, Age distribution of Multi-Ethnic Study of Atherosclerosis (MESA) monocytes methylome dataset, 
used to build stochastic model. Red area highlights 43 youngest samples with age less than 46, purple – 11 
oldest samples with age more than 80. 
**C**, (left) DNAm values are beta distribured. (center) The average DNAm for the youngest (AvgDNAm(Young)) 
and the oldest (AvgDNAm(Old)) samples for each CpG $c$ of the Horvath's clocks. (right) Distribution of an 
absolute value of an effect size.
**D** DNAm dynamics for three CpGs with different effect size in three independent replicates.
```


#### Multi-Ethnic Study of Atherosclerosis (MESA) sorted monocytes methylome dataset

To specify the model, authors used the MESA dataset of methylome of monocytes. 
From this dataset, DNAm for Horvath 353 CpGs were selected
(GEO:
[GSE56046](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56046])){cite}`reynoldsAgerelatedVariationsMethylome2014`. The rationality of this choice is 
determined by the subsequent comparison with Horvath's clocks.  
To simulate stochastic DNAm change, estimated effect size was used, coupled with parameters
($\gamma$, $\sigma$) provided in the paper. The implemented procedure, in principle, resembles
stochastic DNAm dynamics ({numref}`dyn_pop`, D): *(i)* DNAm changes with time monotonously, 
*(ii)* the increase of EffSize increases chance of DNAm change, while 
*(iii)* its' sign defines direction (increase or decrease) of DNAm change.

### Optimizing parameters of stochastic model
In the procedure described above, there are two parameters need to be optimized: global 
methylation change rate $\gamma$, and $\sigma$, defining the methylation magnitude of 
methylation change. So next step of analysis was dedicated to optimization of the parameters such that the obtained 
stochastic clocks would have the EffSize at the end of simulation similar to the observed ones 
({numref}`dyn_pop`, right at C). The parameters were evaluated with mean absolute error 
(MAE, so the less MAE, the better parameters combination).
To do so, the grid search for optimal parameters was used, where simulated EffSize was computed for each 
pair of the parameters ({numref}`optimal_parameters`).

```{figure} figures/3_optimizing_parameters_with_highlight.svg
:name: optimal_parameters
:width: 400
Search for optimal parameters of simulation. The optimal parameters combination is marked.
```

Obtained optimal parameters ($\gamma=9$, $\sigma=0.0004$, $\text{MAE}=0.003$) differ from the optimal parameters 
from paper ($\gamma=9.5$, $\sigma=0.00045$). One reason could be is that I ran only one simulation for 
each pair of parameters, while the authors ran 50 independent simulations per pair.
However, the trend of the MAE on a given grid clearly shows that MAE has a flat extremum 
({numref}`optimal_parameters`). At best, this needs to be checked on a wider grid of parameters. A good way 
to do this is to use gradient descent instead of using a grid approach.


### Constructing the stochastic Horvath clocks (StocH)

#### Simulate an artificial cohorts of samples

For each time step (or, age value), DNAm values can be simulated, using the 
optimized $\sigma$ and $\gamma$ parameters, similarly to the procedure described at previous 
section. For initian DNAm values (relative age 0, or age 45), I used the same average DNAm 
level for young samples.

Following the paper, I simulated three artificial cohorts of DNAm samples: train, model selection
and test sets. Each cohort consists of five samples per age value: from 45 to 83 (195 samples 
in total for each cohort).  
For the samples of year 45 I used the random number of time steps less than 35 (number of stem 
cells divsions per year {cite}`teschendorffComparisonEpigeneticMitoticlike2020`) to has some 
variety in DNAm values for this age.  So, 35 time steps were taken in one year of life for this simulation.  
For each age, the DNAm were simulated from the same DNAm values independently, reflecting 
the stochastic process in each individual sample. Obtained beta-values were that transformed by zeroing the 
mean and scaling to unit variance (with StandardScaler for each cohort independently).

Next, the linear model with regularization (ElasticNet) was train for a set of penalties $\lambda$ 
(from 0 for 1 with step 0.001) on the first cohort. Models were tested on the second cohort to find best 
penalty $\lambda$ ({numref}`optimal_model`, left). Finally, best model was validated on the 
third cohort ({numref}`optimal_model`, right). Low root mean square error (RMSE) reflects 
that obtained weights of the model are sensitive to pure stochastic change of DNAm level.  
Resulting model with best parameters trained on simulated is a stochastic Horvath's clocs (StocH).

```{figure} figures/4_optimizing_lambda.svg
:name: optimal_model
Fitting StocH on simulated aging samples. (Left), Search for best linear model regularization parameter.
(Right), Validation of model ability to predict age caused by stochastic process.
```

### Prediction of age in MESA dataset using StocH
For retained samples from MESA monocytes datased, chronological age was predicted with StocH ({numref}`stoch_component`, 
left upper at A). 
Thus, the predicted age values differ from those found in the original paper, 
mainly reflected in a smaller spread of values. Optimal parameter values imply 
a lower probability and magnitude of transition, and this can explain 
the smaller effect.

```{figure} figures/5_predicting_age_agg.svg
:name: stoch_component
Quantification of stochastic component of epigenetic aging. **A**, Chronological age, predicted with StoH and Horvath clocks for MESA datasets of monocytes and T-cells. **B**, estimated ratio of squared Pearson correlation coefficients. 
```


### Quantification of stochastic component of epigenetic aging

#### Init Horvath clocks

It's not clearly explained, how the Horvath clocks were trained in the paper. One can use the 
same dataset for training and test with cross-validation, without excluding validation subset 
from training data, which is not good approach. In the paper, the number of samples used to 
prediction of age is similar to the whole dataset size (1148, thus, samples in range of 46 to 80). 
This, it's unclear, whether some validation used or not.
The following code can be used to fit Horvath's clocks:
```python
horvath = ElasticNetCV(l1_ratio=0.5, alphas=np.arange(0, 1.00001, 0.001))
horvath.fit(X, y)
y_pred = horvath.predict(X)
```

Another option is to use weights and parameters from the Horvath's initial paper. It should 
be mentioned that the Horvath's datasets did not include the MESA dataset, so basically I 
do not expect data leakage. 

#### Comparing prediction of StocH and Horvath's clocks for MESA monocytes and T-cells
Once the predicted chronological age is obtained, one can try to estimate how much of the accuracy of the Horvath's clock prediction is due to a stochastic process. The ratio of the squared Pearson coefficients for the two predictions can be used as a measure of this component:
$$RR2 = \frac{R^2(\text{StocH})}{R^2(\text{Horvath})}$$

In general, $R^2$ indicates what fraction of an accuracy is explained by the model. For both datasets examined, the Horvath clock gave the best prediction, and based on the ratio above, an estimate of the stochastic component of 0.44 and 0.25 was obtained. Of note, both predictions were below claimed estimates.

## Discussion
As a result, the project achieved all of its objectives. Not all the values obtained during 
the following of the paper methods were reproduced. Particularly, the RR2 values obtained from 
analysis are lower for two used datasets (MESA Mono and CD4T). 
The analysis is available as Python code for review and evaluation at GitHub: 
[ombystoma-young/stochastic-component-epigen-proj](https://github.com/ombystoma-young/stochastic-component-epigen-proj), 
and can, in principle, serve as a setup for further investigations.


### Questions

1. Can resulted simulated DNAm exceed 1? There is no restricting conditions.
2. Which datasets authors used to build clocks? Did they used any validation data?
3. There is no information about $R^2$ definition in paper. Did authors use determination 
scores or just squared Pearson correlation? The $R^2$ for Stoc-clocks can be negative (I have an evidence).
4. Could such a difference in results be explained by change in parameters?


### Suggestions

1. Can we consider non-homogeneous process by define time-dependent global methylation rate $\gamma$?
2. Can we use different values of $\gamma$ for sites with positive and negative EffSize, reasoning different processes underlay the events of methylation and demethylation?
3. Can we use stochastic clocks as random matrix in mixed linear models?

## Credits

This project was performed by [Oksana Kotovskaya](https://x.com/ktv_xn).

## References

```{bibliography}
:style: plain
:filter: docname in docnames
```