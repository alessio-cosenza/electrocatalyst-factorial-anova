# Learning ANOVA through a catalyst experiment

This guide connects experimental reasoning, statistical theory, and the Python notebook. Read one section, run the corresponding notebook cells, and explain the output in your own words before moving on.

## 1. Define the question before choosing a test

The response is **specific activity (SA)**. The three factors are nominal Pt loading, carbon identity, and Pt locality. A factor is an input whose levels are compared. A level is one particular setting. A cell is one combination of all three factors.

| Factor | Levels | Example question |
|---|---|---|
| Nominal loading | 15, 40 wt% | Does average SA differ between loading groups? |
| Carbon | Disordered, Graphitized, FCX-B | Do average activities differ across the three tested supports? |
| Locality | IN, OUT | Does average SA differ between placement categories? |

All 2 × 3 × 2 = 12 combinations are present, with two electrode tests each. The layout is **balanced** at the electrode level: every cell has the same number of observations. It is a **full factorial** layout because every combination is observed.

Factorial analysis is the study of these factors together, including their interactions. ANOVA is a way to divide the response variation into model terms and compare each term with unexplained variation. Here we use a three-factor fixed-effects ANOVA. “Fixed” means the comparisons concern the specified levels, rather than a random sample of all possible carbon supports.

This is different from *factor analysis*, a dimensionality-reduction method that seeks latent variables behind many correlated measurements.

## 2. State what was actually replicated

The two replicates are separate electrode tests of the same catalyst batch. They provide information about electrode preparation/testing variation conditional on that batch. They do not measure variation between independently synthesized batches.

Consequently, the p-values and confidence intervals in this project compare the tested material samples against electrode-level variation. They do not justify population-level claims about the reproducible effects of synthesis loading, carbon, or locality. With one catalyst batch per factor combination, batch-specific properties and treatment combinations are inseparable. Calling the 24 rows independent syntheses would be pseudoreplication.

The errors of separate electrodes may reasonably be modeled as independent conditional on the sample, but shared ink, preparation sessions, instrument drift, or other clustering would need additional metadata. Independence is a property of how data were collected; a normality test cannot establish it.

For generalizable synthesis-level inference, prepare multiple independent catalyst batches per factor combination and test multiple electrodes within each batch. Then consider a hierarchical model with batch nested within the factor combination. Choose the number of batches using desired precision or power and a defensible estimate of batch variation; the current electrode variation alone is insufficient for that calculation.

## 3. Understand main effects

A main effect summarizes the difference after averaging over the other factors. Because this layout is balanced, ordinary averages give equal weight to every level combination.

For loading:

\[
\Delta_L = \overline{SA}_{40} - \overline{SA}_{15} = 30.266.
\]

The estimated equal-weight average increase is about 30.27 µA/cm²Pt. The conditional pointwise 95% confidence interval is [26.95, 33.58]. This is an average across the tested carbons and localities, not a universal increment for every catalyst.

For locality, the average OUT − IN difference is 10.551, with a conditional pointwise 95% interval [7.24, 13.86]. Carbon has three levels, so its main-effect F-test jointly asks whether all three marginal means are equal. It does not identify which pair differs or give a direction by itself.

## 4. Understand two-factor interactions

An interaction means that the difference associated with one factor changes when another factor changes. On the raw SA scale, the loading × locality interaction asks whether the OUT − IN **absolute difference** changes between loading levels.

For Disordered carbon:

| Nominal loading | IN mean | OUT mean | OUT − IN |
|---|---:|---:|---:|
| 15 wt% | 27.550 | 31.450 | 3.900 |
| 40 wt% | 55.700 | 78.150 | 22.450 |

The difference of differences is:

\[
(78.15-55.70)-(31.45-27.55)=22.45-3.90=18.55.
\]

Thus the observed OUT advantage grows by 18.55 on this carbon. In an interaction plot, parallel mean lines imply a constant difference; nonparallel lines show an observed interaction pattern. Noise can also produce nonparallel lines, so the graph must be interpreted with uncertainty.

“One comparison is significant and another is not” is not a test that the two effects differ. The difference of differences directly tests that question.

## 5. Understand the three-factor interaction

The loading × locality pattern can itself vary across carbons:

| Carbon | OUT − IN at 15 | OUT − IN at 40 | Change in the gap |
|---|---:|---:|---:|
| Disordered | 3.900 | 22.450 | +18.550 |
| Graphitized | 8.640 | 4.350 | −4.290 |
| FCX-B | 6.600 | 17.365 | +10.765 |

A three-factor interaction asks whether these changes are equal across the three carbons. It has two degrees of freedom because there are three carbon-specific changes to compare.

This is a useful scientific question even when the carbon main effect is small. Similar carbon averages can conceal different behavior at individual loading/locality combinations.

## 6. Write the model you are fitting

Let L denote loading, C carbon, P Pt locality, and r the electrode test. The model is:

\[
y_{ijkr}=\mu+L_i+C_j+P_k+(LC)_{ij}+(LP)_{ik}+(CP)_{jk}+(LCP)_{ijk}+\epsilon_{ijkr}.
\]

The terms represent the grand mean, three main effects, three two-factor interactions, one three-factor interaction, and an error. Sum-to-zero constraints identify the effect parameters. The parameters describe these tested samples conditional on their batches.

ANOVA and ordinary least squares regression are the same linear-model framework here. OLS chooses coefficients to minimize squared residuals. ANOVA organizes the fitted model into hypotheses about groups of coefficients. The categorical model does not require a straight-line relationship between SA and measured loading.

There are 12 freely estimated cell means: one intercept plus 11 effect degrees of freedom. The full model reproduces each sample's mean exactly. It is saturated for the 12 means, but still has 24 − 12 = 12 residual degrees of freedom from the electrode repeats. Its high R² describes separation of these observed means relative to electrode variation; it is not evidence of prediction quality on new batches.

## 7. Follow the variance calculation

Start with total variation around the grand mean:

\[
SS_T=\sum_i(y_i-\bar y)^2.
\]

For this balanced full factorial model, it decomposes into:

\[
SS_T=SS_L+SS_C+SS_P+SS_{LC}+SS_{LP}+SS_{CP}+SS_{LCP}+SS_E.
\]

Residual variation is the variation around the relevant cell mean:

\[
SS_E=\sum_{cells}\sum_r(y_{cell,r}-\bar y_{cell})^2=166.27265.
\]

For the first OUT/15/Disordered sample, the mean is (33.9 + 29.0)/2 = 31.45. Its contribution is (33.9 − 31.45)² + (29.0 − 31.45)² = 12.005. Add the analogous contributions from all 12 samples.

Degrees of freedom count independent pieces of information:

| Source | df |
|---|---:|
| Loading | 2 − 1 = 1 |
| Carbon | 3 − 1 = 2 |
| Locality | 2 − 1 = 1 |
| Loading × Carbon | 1 × 2 = 2 |
| Loading × Locality | 1 × 1 = 1 |
| Carbon × Locality | 2 × 1 = 2 |
| Three-factor interaction | 1 × 2 × 1 = 2 |
| Residual | 12 × (2 − 1) = 12 |
| Total | 24 − 1 = 23 |

Divide a sum of squares by its df to obtain its **mean square**. The pooled error mean square is 166.27265/12 = 13.8561. Its square root, 3.7224 µA/cm²Pt, estimates the conditional within-sample electrode standard deviation under the common-variance assumption.

For a term, calculate:

\[
F=\frac{MS_{term}}{MS_E}.
\]

For loading, F = (5496.123/1)/13.8561 ≈ 396.66. Under that term's null hypothesis and the error assumptions, F follows an F distribution with numerator df equal to the term df and denominator df 12. The p-value is the probability of an F at least as large under that null model.

## 8. Interpret p-values, uncertainty, and effect size separately

A p-value is not the probability that the null hypothesis is true, the probability the result happened “by chance,” or the size of a physical effect. A small p-value indicates incompatibility with the specified null/error model. A large p-value does not establish equivalence or absence of an effect.

Report an effect estimate in SA units, a confidence interval, and the analysis scope. A frequentist 95% CI is produced by a procedure that covers the fixed true quantity in 95% of repeated comparable experiments under the model; it is not a guarantee about this particular interval or a range containing 95% of individual electrodes.

For an OUT − IN comparison of two cells with two electrodes each:

\[
SE(\bar y_{OUT}-\bar y_{IN})=\sqrt{MS_E(1/2+1/2)}=3.7224.
\]

The interval is estimate ± t(0.975,12) × SE, giving a half-width of about 8.1104. These intervals pool error across all 12 cells and require a common variance. They are **pointwise**, not simultaneous intervals for all comparisons.

Eta squared is SS_term/SS_total. Loading accounts for about 80.9% of the observed centered SS in this balanced partition. Partial eta squared is SS_term/(SS_term + SS_error); for loading it is about 97.1%. These have different denominators. Partial eta squared values do not add to one, and neither quantity means that loading causally explains that percentage of catalyst physics.

## 9. Handle multiple questions explicitly

Looking at many tests increases the chance of at least one false positive. This project reports raw p-values and Holm-adjusted p-values, with the tested family named in each column:

- Seven omnibus factorial terms: a broad exploratory family.
- Six OUT − IN comparisons: one within each carbon/loading combination.
- Three carbon-specific loading × locality difference-of-differences: a separate exploratory family.

The adjustments control each named family, not every test across the entire notebook. The analyses are educational and exploratory, not preregistered. The log-scale analysis is a sensitivity analysis, not an additional search for significant results.

The raw-scale loading × carbon, loading × locality, and three-factor interaction p-values are about 0.0196, 0.0178, and 0.0284. Each becomes about 0.0889 after Holm adjustment across the seven terms. Treat their patterns as useful leads, not settled discoveries. Loading and locality remain clearly distinguishable relative to electrode error under that adjustment.

Among the six locality contrasts, the Disordered/40 and FCX-B/40 comparisons have Holm-adjusted p-values below 0.05. The Graphitized/15 comparison has raw p ≈ 0.0387 but adjusted p ≈ 0.1548. Its pointwise CI excludes zero even though the adjusted test does not reject: the CI and adjusted p-value address different multiplicity conventions.

## 10. Check the assumptions honestly

Classical F-tests require an appropriate mean model, independent conditional errors, common error variance, and normally distributed errors for exact small-sample inference. They do not require the pooled SA observations from different treatments to look normal.

Use residual-versus-fitted and normal Q–Q plots as descriptive checks. With only two tests per cell, each pair of residuals is mechanically +d/2 and −d/2. Symmetry in the Q–Q plot is therefore partly built into the model; it does not validate normality. Each cell SD has only one degree of freedom, so large apparent differences in SD may be unstable.

We do not use a normality-test p-value as permission to run ANOVA. We also avoid a median-centered Levene test across the 12 two-observation cells: both absolute deviations in a cell are identical, producing a degenerate variance-test setup. Robust covariance estimates or bootstrapping cannot create missing synthesis replication. Deleting observations simply to improve significance is inappropriate.

Replicate labels are not a random effect by default. If electrode 1 for every sample was run together on day 1 and electrode 2 on day 2, that would be a meaningful block to consider. The current labels alone do not establish such a structure.

## 11. Separate nominal loading from measured loading

The primary model compares the nominal 15 and 40 wt% preparation groups. It does not estimate an SA change per one measured wt% point. Exact measured loadings vary between localities: for low-loading FCX-B, OUT is 12.3 wt% and IN is 18.8 wt%. A locality comparison therefore compares materials that also differ somewhat in measured composition and may differ in other properties.

Do not add measured loading as an extra covariate to this full categorical model and claim it was independently controlled. It is constant within each cell and therefore already lies in the space spanned by the 12 fitted cell means. Adding its column leaves matrix rank at 12 while creating 13 columns: its additional coefficient is not separately identifiable.

A continuous-loading model asks a different question and imposes a functional form, such as linearity. With only two loading settings per carbon/locality pair, a straight line fits the two means exactly; this cannot validate linearity. A better future design includes intermediate loading values, closer IN/OUT matching, and independent synthesis batches.

## 12. Remember that interactions depend on scale

Raw SA compares absolute changes. Log(SA) compares proportional changes; a constant multiplicative ratio can create a raw-scale interaction as baseline activity increases.

The included log(SA) sensitivity analysis gives p ≈ 0.597 for loading × locality and p ≈ 0.166 for the three-factor term, versus about 0.0178 and 0.0284 on the raw scale. Thus those raw-scale patterns are not scale-invariant findings. Choose the scientific estimand—absolute activity gain or relative change—before interpreting an interaction, rather than choosing whichever scale produces a smaller p-value.

## 13. Read the Python model

```python
model = ols(
    "SA ~ C(Loading_level, Sum) * C(Carbon, Sum) * C(Locality, Sum)",
    data=df,
).fit()
table = anova_lm(model, typ=3)
```

`SA` is the response. `~` means “modeled as a function of.” `C(...)` treats a variable as categorical. `Sum` sets sum-to-zero coding so the main effects represent equal-weight marginal comparisons in this complete balanced design. `*` expands to all main effects and interactions; it is formula syntax, not multiplication of measured values. `A:B` denotes just an interaction. `.fit()` estimates the coefficients by OLS. `typ=3` requests each term's test in the full model; we omit the intercept test because whether the grand mean is zero is not our research question.

Type I SS are sequential and can depend on term order in an unbalanced dataset. Type II tests account for other terms except interactions containing the term being tested. Type III tests account for the full model, but main-effect hypotheses depend on coding. In this complete balanced design, Type III with sum contrasts gives the familiar orthogonal factorial partition. Do not transfer this recipe blindly to missing cells or an unbalanced dataset, or pair Type III with default treatment coding and assume the same hypotheses.

The ordinary coefficient p-values in `model.summary()` are not the same table as the joint ANOVA tests. For example, carbon occupies two coefficient dimensions; ANOVA tests them jointly. Use direct mean contrasts for interpretable SA differences instead of interpreting an arbitrary coded coefficient as the whole effect.

## 14. Explain what you did

**Methods example:** “We organized specific activity measurements in a 2 × 3 × 2 factorial layout comprising nominal Pt loading, carbon support, and Pt locality. We fitted an ordinary least squares model containing all main effects and interactions and calculated Type III ANOVA tests using sum contrasts. Differences were quantified with model-based contrasts, pointwise 95% confidence intervals, and Holm-adjusted p-values within stated comparison families. Each catalyst sample was measured using two separately prepared electrodes; inferential results therefore describe electrode-level variation conditional on the tested batches and do not quantify synthesis-to-synthesis reproducibility.”

**Results example:** “The average specific activity difference between nominal 40 and 15 wt% samples was 30.27 µA/cm²Pt (conditional 95% CI 26.95–33.58). The OUT − IN difference varied across loading/support combinations; for Disordered carbon it increased from 3.90 to 22.45 µA/cm²Pt. The omnibus interaction patterns were exploratory: their raw p-values were below 0.05, but their Holm-adjusted p-values across seven factorial terms were approximately 0.089, and some patterns changed on the logarithmic scale.”

**Portfolio example:** “Built a reproducible Python workflow to analyze a 24-observation catalyst experiment, estimate factorial effects and interactions, quantify uncertainty, and identify the limits imposed by electrode-level replication and measured-loading differences.”

Statistical analysis adds numerical effect estimates, uncertainty, explicit interaction tests, and transparent limits to a visual trend. It cannot recover a carbon-property mechanism, establish causation, or create independent replication. Carbon identity is a category here; separating surface area, porosity, or chemistry would require suitable measured predictors and a design that distinguishes their effects.

## 15. Practice before reading the answers

1. Calculate the mean of the two high-loading OUT/FCX-B SA values. Then subtract the IN mean.
2. Why does carbon have two main-effect degrees of freedom?
3. If every OUT value were exactly 10 units above its corresponding IN mean at both loadings, what would loading × locality represent?
4. Why does more precise electrode measurement not solve missing synthesis replication?
5. What changes if SA is modeled on the log scale?

Answers: (1) OUT mean 72.60; IN mean 55.235; difference 17.365. (2) Three means have two independent comparisons after accounting for the grand mean. (3) The interaction would be zero for those mean differences. (4) It reduces one source of error while leaving batch variation unidentified. (5) Comparisons concern ratios/proportional differences rather than absolute differences.

## References

- [NIST: multi-factor ANOVA](https://www.itl.nist.gov/div898/handbook/eda/section3/eda355.htm)
- [NIST: modeling designed experiments](https://www.itl.nist.gov/div898/handbook/pri/section4/pri43.htm)
- [NIST: residual assumptions](https://www.itl.nist.gov/div898/handbook/pri/section2/pri245.htm)
- [statsmodels: interactions, ANOVA, and sum contrasts](https://www.statsmodels.org/stable/examples/notebooks/generated/interactions_anova.html)
- [statsmodels: anova_lm](https://www.statsmodels.org/stable/generated/statsmodels.stats.anova.anova_lm.html)
- [statsmodels: multiple-testing adjustments](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html)

All numerical results above are calculated from the supplied dataset. The references support statistical methods, not claims about these particular catalysts.
