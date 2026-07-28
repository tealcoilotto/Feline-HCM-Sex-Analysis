# Results & Discussion: (MV)MR Analysis of Depression, Risk-Taking, Alcohol Use Disorder (AUD), and Non-Alcohol Substance Use Disorders

This document provides an in-depth discussion of the two-sample MR and MVMR results from [`MVMR_Analysis.R`](MVMR_Analysis.R). All figures referenced below are written to `mr_results/` when the script is run. All plots and tables produced by the script are also readily available in the [`Figures`](Figures) and [`Tables`](Tables) folders.

In addition to the (MV)MR pipeline, a bivariate [`LD score regression`](LDSC_Genetic_Correlation_Depression_Risk.R) (LDSC via GenomicSEM) was run on the major depression and risk-taking GWAS summary statistics to estimate the genome-wide correlation between the two exposures. This serves as a secondary check on the discussed MVMR attenuation findings.

---

## Background

### Why Mendelian Randomization?

Observational associations between depression, risk-taking, and substance-related outcomes are well established but are vulnerable to confounding (e.g., shared environmental adversity, socioeconomic factors) and reverse causation (e.g., AUD contributing to symptoms of depression rather than the reverse). Mendelian randomization (MR) addresses these limitations by using genetic variants associated with each exposure as instruments, serving as proxies for lifetime genetic liability to a trait. Because genotype is fixed at conception, these instruments are not impacted by later-life outcomes, environmental factors, and other confounding factors that arise when directly measuring exposures. However, MR relies on three core assumptions: 

1. **Relevance:** The instruments are robustly associated with the exposure (assessed with F-statistics in this analysis)
2. **Independence:** The instruments affect the outcome only through the exposure and not through other pathways (no horizontal pleiotropy,  assessed via MR-Egger, Cochran's Q, and MR-PRESSO)
3. **Exclusion Restriction:** The instruments are independent of other factors that may influence the outcome (e.g., ancestry, lifestyle). This cannot be directly tested, but using ancestry-matched and unrelated samples helps reduce the risk of this assumption being violated.
   
Two-sample MR draws on exposure and outcome SNP associations from two separate, non-overlapping (or minimally overlapping) genome-wide association study (GWAS) datasets. This increases statistical power and avoids many forms of confounding often present in single-sample analyses. The tradeoff is that both samples should come from comparable populations. In this analysis, both exposure and outcome GWAS datasets are sourced from European cohorts (mostly due to limited public availability of other datasets). Multivariable MR (MVMR) extends this further to estimate the independent effect of one exposure on an outcome while adjusting for additional exposures. In this analysis, MVMR is used to understand the distinct contributions of depression and risk-taking on the risk of AUD and non-alcohol substance use disorders.

### Motivating Question

This project aims to understand if genetically predicted liability to depression and risk-taking shows evidence of causal relationships with alcohol use disorder (AUD), and its complement, non-alcohol substance use disorders, and if each exposure's effect is independent of the other.

### Exposures

*Note: both exposure datasets reflect European populations, and both outcome datasets reflect Finnish populations, a key limitation of these analyses. Ideally, these analyses would be performed on larger, global datasets to reduce the risk of geographic and demographic confounding. However, that may cause issues with the MR assumptions.*

**Major depression** ([`ieu-b-102`](https://opengwas.io/datasets/ieu-b-102), UK Biobank and PGC - excluding 23andMe): Major depression is a broad, highly polygenic mood disorder with an estimated SNP heritability around ~8.9% [6]. It is genetically and phenotypically linked to substance use disorders through multiple potential pathways, including self-medication (alcohol/substance use to cope with symptoms), shared neurobiological vulnerability (e.g., disrupted reward pathways), and shared genetic risk factors [7]. *This dataset is binary for case/control, with log-odds units.*

**Risk-taking** ([`ukb-b-14147`](https://opengwas.io/datasets/ukb-b-14147), UK Biobank): In this dataset, risk-taking is a self-reported behavioral trait (*"would you describe yourself as someone who takes risks?"*), characterizing a general risk tolerance rather than a specific behavior. It has a slightly lower SNP heritability of around ~5% [8] and has been linked in prior literature to many externalizing behaviors, including substance use and substance use disorders [9]. *This dataset is binary self-report (yes/no), with standardized (standard deviation) units.*

### Outcomes

**Alcohol use disorder (AUD)** ([`finn-b-AUD`](https://opengwas.io/datasets/finn-b-AUD), FinnGen): Alcohol use disorder is a clinical diagnosis reflecting a problematic pattern of alcohol use, leading to impairment and/or distress [10]. In this sample, it has a ~4% prevalence, which is slightly lower than the global rate of ~7% among adults [11]. *This dataset is binary for case/control, with log-odds units.*

**Substance use, excluding alcohol** ([`finn-b-F5_SUBSNOALCO`](https://opengwas.io/datasets/finn-b-F5_SUBSNOALCO), FinnGen): Despite this dataset's name, it does not reflect general substance use. It captures non-alcohol substance use disorders, spanning opioids, cannabis, sedatives/hypnotics, cocaine, other stimulants, hallucinogens, tobacco, volatile solvents, and other psychoactive substances. Like AUD, substance use disorder is a clinical diagnosis reflecting a problematic pattern of substance use, leading to impairment and/or distress. In this sample, it has a ~1.8% prevalence (global comparison omitted - global rates in adults are highly variable / often include alcohol). *This dataset is binary for case/control, with log-odds units.*

---

## Methods

**⚠️ Note on Effect Scales:** Depression and the two outcomes (AUD, non-alcohol substance use disorders) were derived from logistic GWAS and are expressed in log-odds. Risk-taking, however, is expressed on a standardized (SD) scale. As a result, a single unit change in genetically predicted liability is not directly equivalent between depression and risk-taking. For depression, a single unit change represents a single log-odds unit, while for risk-taking, a single unit change represents a single standard deviation. This does not affect the validity of MR estimates for each exposure individually, or comparisons of significance and attenuation across models, but direct comparisons of effect-size *magnitude* between the two exposures should be interpreted with these scale differences in mind.

### Metric Definitions

* **β (beta):** The estimated causal effect of the genetically predicted exposure on the outcome. Since both outcomes are binary and derived from logistic GWAS, β represents the estimated change in outcome log-odds per unit increase in the exposure. Note that one unit is defined as log-odds for depression, but standard deviation for risk-taking.
* **SE (standard error):** The estimated uncertainty of the β estimate, expressed on the same scale as β.
* **Odds Ratio (OR):** A more interpretable form of the β estimate (OR = e^β), representing the multiplicative change in odds of the outcome per single unit increase of the genetically predicted exposure (again, one unit is defined differently for each exposure).
* **F-statistic:** A unitless measure of the genetic instrument strength, with larger values indicating stronger instruments.
* **R²:** The proportion of variance in the exposure explained by the evaluated genetic instruments, reported as a proportion or percentage.
* **Conditional F-statistic:** A measure of instrument strength for an exposure in MVMR after accounting for the other exposures. 
* **Cochran's Q:** A measure of heterogeneity among SNP-specific causal estimates. Q is unitless, with larger values indicating greater heterogeneity.
* **MR-Egger Intercept:** An estimate of the average directional pleiotropy across the evaluated genetic instruments. The intercept is expressed on the β scale.
* **MR-PRESSO Global Test:** A statistical test for the presence of horizontal pleiotropy due to outlier variants. The result is reported as a p-value.
* **Steiger Directionality Test:** A test assessing whether genetic instruments explain more variance in the exposure than the outcome, providing evidence for the direction of the hypothesized causal relationship. The result is reported as a p-value.
* **qhet_mvmr:** A heterogeneity-robust MVMR estimator that accounts for heterogeneity across genetic instruments when estimating direct causal effects. Estimates are expressed on the β scale.
* **SNP Heritability (h²):** The proportion of phenotypic differences that can be explained by common genetic variants, reported as a proportion or percentage.
* **LDSC Intercept:** A unitless check for inflation of GWAS results not attributable to genetic signal. A value close to 1 indicates limited inflation.
* **Genetic Covariance (g_cov):** The estimated covariance between genetic effects for two traits.
* **Genetic Correlation (r_g):** The standardized genetic covariance between two traits, ranging from −1 to 1.
* **h² Z-statistic:** The ratio of the estimated SNP heritability to its SE, measuring the statistical evidence for nonzero SNP heritability.
* **Genetic Covariance Z-statistic:** The ratio of the estimated genetic covariance to its SE, measuring evidence for nonzero genetic covariance.

### Instrument Selection
- Genome-wide significant SNPs (p < 5×10⁻⁸) were extracted for each exposure via `TwoSampleMR::extract_instruments()`
- LD clumping was used to retain independent instruments (r² < 0.001, 10,000 kb distance)
- Risk-taking: 29 independent instruments (mean F-statistic = 39.9)
- Depression: 48 independent instruments (mean F-statistic = 38.7)
- Proxy SNPs were used, when necessary, for exposure instruments during extraction

### Harmonization
Exposure and outcome datasets were harmonized using `TwoSampleMR::harmonise_data()` (action = 2), aligning effect alleles and inferring the positive strand orientation. Palindromic SNPs were removed when allele alignment could not be reliably determined.

### Steiger Directionality
For each exposure-outcome combination, `TwoSampleMR::directionality_test()` was used to assess whether genetic instruments explained more variance in the exposure than the outcome, supporting the assumed causal direction.

### Univariable MR
Primary estimates were computed using inverse-variance weighted (IVW) regression. Sensitivity analyses included MR-Egger, weighted median, and weighted mode methods. Associations were considered statistically significant at the p < 0.05 significance level.

### Sensitivity and Pleiotropy Diagnostics
For each exposure-outcome combination, the following were computed:
- Cochran's Q (heterogeneity across instruments)
- MR-Egger intercept test (directional pleiotropy)
- Leave-one-out analysis (influence of individual SNPs)
- MR-PRESSO (global test and outlier detection)
- Instrument strength (mean F-statistic)

### Multivariable MR (MVMR)
Risk-taking and depression were jointly assessed for each outcome using `TwoSampleMR::mv_extract_exposures()` and `TwoSampleMR::mv_multiple()`, with conditional instrument strength (F-statistics) assessed through `MVMR::strength_mvmr()`. Pleiotropy and heterogeneity were evaluated using the MVMR Q-statistic, MVMR-Egger, and the heterogeneity-robust `qhet_mvmr` estimator, with sensitivity analyses performed across a range of instrument correlation parameter (`pcor`) values.

### Software
Full R version and package details are available in [`session_info.txt`](session_info.txt) and the [`README.md`](README.md).

---

## 1. Univariable Mendelian Randomization Analysis of the Effect of Liability for Risk-Taking on Alcohol Use Disorder (AUD) Risk

### Instruments and Directionality

29 independent SNPs (r² < 0.001, 10,000 kb clumping window) were retained as genetic instruments for risk-taking, with one proxy substitution (`rs3943670` for `rs7829912`) and no SNPs lost during harmonization. Instrument strength was strong throughout (mean F ≈ 40, median F ≈ 35, all SNPs F > 10), and the instruments explained ~0.259% of variance in risk-taking liability, which is modest but more typical for polygenic behavioral traits. The Steiger directionality test supported the hypothesized causal direction (risk-taking → AUD) with high confidence (p = 2.54 * 10⁻⁴⁶), meaning the genetic instruments explain significantly more variance in risk-taking than in AUD, as expected if risk-taking is upstream.

### Main Estimate

<br>

<div>
  
| Method| β | β 95% CI | SE | OR | OR 95% CI | p-value |
|:-------:|:------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| IVW | 1.57 | (0.416, 2.73) | 0.590 | 4.81 | (1.52, 15.3) | 0.00770 |
| MR Egger | 0.966 | (-4.62, 6.55) | 2.85 | 2.63 | (0.00984, 701) | 0.737 |
| Weighted Median | 1.49 | (-0.0262, 3.00) | 0.77 | 4.42 | (0.974, 20.1) | 0.054 |
| Weighted Mode | 1.18 | (-1.48, 3.83) | 1.35 | 3.24 | (0.227, 46.2) | 0.393 |

</div>
<p><i>Table 1: Comparison of Univariable Methods for Risk-Taking → Alcohol Use Disorder</i></p>

<br>

The scatter plot in Figure 1 shows all four MR methods used (IVW, MR-Egger, weighted median, weighted mode) estimating a positive slope, meaning higher genetic liability for risk-taking is associated with a greater AUD risk. Only the IVW estimate reached significance at p < 0.05, though all methods agreed in direction, a sign of consistency among methods with varying assumptions and validity for this dataset. The majority of individual SNP estimates clustered together, with a few points exhibiting smaller effects on AUD or larger effects on risk-taking.

<br>

![Risk-Taking → AUD scatter](Figures/risk_aud_scatter_plot.png)
<p align="center"><i>Figure 1: Scatter Plot of Risk-Taking → Alcohol Use Disorder, Across Methods</i></p>

<br>

The method comparison forest plot in Figure 2 details the IVW confidence interval fully above the null, with the weighted median having borderline significance, and MR-Egger / weighted mode having wide and imprecise intervals. The large bounds of MR-Egger and weighted mode methods are consistent with their lower statistical power and do not necessarily suggest conflicting results.

<br>

![Risk-Taking → AUD odds ratio forest plot](Figures/risk_aud_methods_forest_plot.png)
<p align="center"><i>Figure 2: Odds Ratios across MR Methods for Risk-Taking Effects on AUD</i></p>

<br>

### Heterogeneity and Pleiotropy

Cochran's Q showed no significant heterogeneity, and the MR-Egger intercept was not significant, collectively suggesting the SNP-specific estimates are reasonably consistent with each other and that directional pleiotropy does not significantly distort the IVW estimate. Additionally, MR-PRESSO results did not detect any significant pleiotropic outliers.

Despite the reassuring global statistics, the SNP forest and funnel plots (Figures 3-4) reveal that a few SNPs exhibit directional discordance. The funnel plot is roughly symmetric around the IVW estimate, evidence against systematic directional pleiotropy, despite individual SNP estimates varying in magnitude and sign.

<br>

![Risk-Taking → AUD SNP forest plot](Figures/risk_aud_forest_plot.png)
<p align="center"><i>Figure 3: SNP Forest Plot for Risk-Taking → Alcohol Use Disorder</i></p>

<br>
<br>

![Risk-Taking → AUD funnel plot](Figures/risk_aud_funnel_plot.png)
<p align="center"><i>Figure 4: SNP Funnel Plot for Risk-Taking → Alcohol Use Disorder</i></p>

<br>

### SNP-Level Sensitivity: rs279846 (*GABRA2*)
One SNP, `rs279846`, stood out in the forest plot, potentially skewing the pooled estimate downward. This variant lies within *GABRA2*, a gene with a well-established link to alcohol and substance disorders [12]. It would be plausible that this SNP's effect on liability to AUD may reflect biological pathways beyond risk-taking, raising the risk of horizontal pleiotropy. Therefore, sensitivity checks for `rs279846` were conducted: 

- **IVW estimate with rs279846 included:** β = 1.57, Odds Ratio = 4.81, p = 0.00770
- **IVW estimate with rs279846 excluded:** β = 1.93, Odds Ratio = 6.89, p = 0.000508
- The SNP-specific F-statistic is 32.1, ruling out weak-instrument artifacts as the explanation for its divergence.

Once `rs279846` is removed, the IVW estimate direction holds, but with increased magnitude and significance. Because neither Cochran's Q or MR-PRESSO flagged significant heterogeneity or pleiotropy issues, `rs279846` was retained in the primary analysis. Notably, the removal of this potentially pleiotropic SNP strengthened the estimated effect (as opposed to weakening), suggesting that the significance of the risk-taking and AUD relationship is not the result of one extreme variant.

### Leave-One-Out

The leave-one-out plot in Figure 5 shows that no single SNP disproportionately drives the estimate, with all confidence intervals above the null. Removing `rs279846` increases the overall estimate more than removing any other SNP, consistent with the above sensitivity checks.

<br>

![Risk-Taking → AUD leave-one-out](Figures/risk_aud_loo_plot.png)
<p align="center"><i>Figure 5: Leave-One-Out Plot for Risk-Taking → Alcohol Use Disorder</i></p>

<br>

### Overall Interpretation

Taken together, these results provide consistent evidence of a positive causal effect of genetic liability for risk-taking on alcohol use disorder risk. The IVW estimate was significant, supported by strong genetic instruments and correct Steiger directionality. Sensitivity analyses showed no significant evidence of heterogeneity or directional pleiotropy, and all four methods were directionally consistent, regardless of significance. Additionally, the assessment of a potential pleiotropic SNP did not weaken, and if anything, strengthened the observed result.

The univariable IVW analysis estimated that each standard deviation increase in genetic liability to risk-taking was associated with a 4.81-fold higher odds of alcohol use disorder, consistent with a potential causal effect of risk-taking liability on risk of alcohol use disorder, under the assumptions of Mendelian randomization.

---

## 2. Univariable Mendelian Randomization Analysis of the Effect of Liability for Depression on Alcohol Use Disorder (AUD) Risk

### Instruments and Directionality

48 independent SNPs (r² < 0.001, 10,000 kb clumping window) were retained as genetic instruments for depression, with 3 proxy substitutions and 2 SNPs omitted during harmonization. Instrument strength was strong throughout (mean F ≈ 39, median F ≈ 35, all SNPs F >10), and the instruments explained ~0.371% of variance in depression liability. The Steiger directionality test supported the hypothesized causal direction (depression → AUD) with high confidence (p = 7.80 * 10⁻⁵⁸), meaning the genetic instruments explain significantly more variance in depression than in AUD, as expected if depression is upstream.

### Main Estimate

<br>

<div>
  
| Method| β | β 95% CI | SE | OR | OR 95% CI | p-value |
|:-------:|:------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| IVW | 0.411 | (0.188, 0.634) | 0.114 | 1.51 | (1.21, 1.89) | 0.000300 |
| MR Egger | 0.749 | (-0.610, 2.11) | 0.694 | 2.11 | (0.543, 8.23) | 0.286 |
| Weighted Median | 0.293 | (0.0144, 0.572) | 0.142 | 1.34 | (1.01, 1.77) | 0.0393 |
| Weighted Mode | 0.161 | (-0.455, 0.778) | 0.315 | 1.17 | (0.634, 2.18) | 0.611 |

</div>
<p><i>Table 2: Comparison of Univariable Methods for Depression → Alcohol Use Disorder</i></p>

<br>

Figure 6 depicts that all four MR methods used (IVW, MR-Egger, weighted median, weighted mode) estimated a positive slope, meaning higher genetic liability for depression is associated with greater AUD risk. The IVW and weighted median estimates both reached significance at p < 0.05, and all methods agreed in direction, a sign of consistency among methods with different assumptions and validity for this dataset. The majority of individual SNP estimates clustered together, with a few points having larger estimates.

<br>

![Depression → AUD scatter](Figures/dep_aud_scatter_plot.png)
<p align="center"><i>Figure 6: Scatter Plot of Depression → Alcohol Use Disorder, Across Methods</i></p>

<br>

The method comparisons in Figure 7 show the IVW and weighted median confidence intervals fully above the null, with MR-Egger and the weighted mode having wide and imprecise intervals. The large bounds of MR-Egger and the weighted mode methods are consistent with their lower statistical power and do not necessarily suggest conflicting results among the four effect estimates.

<br>

![Depression → AUD odds ratio forest plot](Figures/dep_aud_methods_forest_plot.png)
<p align="center"><i>Figure 7: Odds Ratios across MR Methods for Depression Effects on AUD</i></p>

<br>

### Heterogeneity and Pleiotropy

Unlike risk-taking, Cochran's Q was significant here, indicating heterogeneity among the 48 SNP-specific estimates, which is reasonable given the polygenic nature of depression. However, the MR-Egger intercept was not significant, and MR-PRESSO did not detect any significant pleiotropic outliers. Collectively, these results suggest that individual SNP effects may be noisy around a true shared effect, with little evidence of directional pleiotropy or confounding from other pathways.

The SNP forest and funnel plots in Figures 8 and 9 illustrate that a few SNPs exhibit directional discordance, as expected given the significant evidence of heterogeneity via Cochran's Q. The funnel plot is fairly symmetric around the IVW estimate, not indicating systematic directional pleiotropy, despite variability in sign and magnitude among select SNPs.

<br>

![Depression → AUD SNP forest plot](Figures/dep_aud_forest_plot.png)
<p align="center"><i>Figure 8: SNP Forest Plot for Depression → Alcohol Use Disorder</i></p>

<br>
<br>

![Depression → AUD funnel plot](Figures/dep_aud_funnel_plot.png)
<p align="center"><i>Figure 9: SNP Funnel Plot for Depression → Alcohol Use Disorder</i></p>

<br>

### Leave-One-Out

The leave-one-out plot in Figure 10 highlights that no single SNP is disproportionately driving the estimated effect, with all confidence intervals remaining above the null with similar ranges. 

<br>

![Depression → AUD leave-one-out](Figures/dep_aud_loo_plot.png)
<p align="center"><i>Figure 10: Leave-One-Out Plot for Depression → Alcohol Use Disorder</i></p>

<br>

### Overall Interpretation

Altogether, these univariable results provide consistent evidence of a positive causal effect of genetic liability for depression on alcohol use disorder risk. The IVW and weighted median estimates were both significant, supported by strong genetic instruments, correct Steiger directionality, and consistent direction across all estimates, regardless of significance. Sensitivity analyses showed no evidence of directional pleiotropy, despite the presence of SNP-level heterogeneity, which did not appear to distort the overall effect estimate.

The univariable IVW and weighted median analyses estimated that each log-odds unit increase in genetic liability to depression was associated with 1.51-fold and 1.34-fold higher odds of alcohol use disorder, respectively. This is consistent with a potential causal effect of depression liability on the risk of alcohol use disorder, under the assumptions of Mendelian randomization.

---

## 3. Multivariable Mendelian Randomization Analysis on the Effect of Liability to Risk-Taking and Depression on Alcohol Use Disorder (AUD) Risk

This analysis seeks to understand whether liability for risk-taking and depression have independent effects on AUD risk, or if one of the exposures is capturing their shared genetic liability.

### Instrument Overlap

67 unique SNPs were jointly retained as genetic instruments after harmonization, with zero overlapping SNPs between the risk-taking and depression instrument sets (confirmed in Part 7), suggesting that any observed attenuation reflects shared underlying genetic liability or correlated pathways, rather than the same variants directly appearing in both exposure sets.

### Conditional Instrument Strength

Both exposures maintained adequate conditional F-statistics (F > 10), with the conditional F-statistic for depression being noticeably larger than that of risk-taking. This indicates that once one accounts for the correlation between the two exposures, the risk-taking instruments retain comparatively less unique predictive power, foreshadowing an observed attenuation pattern.

### Main MVMR Estimate

Figure 11 details that after mutual adjustment, liability for depression retains a significant independent association with AUD, while the estimate for risk-taking lost significance, with a much wider confidence interval. However, both point estimates remained positive, in alignment with the univariable results. More details on the specific estimates are outlined in the next subsection.

<br>

![MVMR AUD forest plot (OR)](Figures/mv_aud_forest_OR_plot.png)
<p align="center"><i>Figure 11: MVMR Odds Ratios for Risk-Taking + Depression → Alcohol Use Disorder</i></p>

<br>

### Direct Comparison: Univariable vs. MVMR-Adjusted Estimates

As detailed in Table 3 and Figure 12, the point estimate for depression barely moves between the univariable and MVMR-adjusted models, while the estimate for risk-taking shifts substantially toward the null. 

<br>

<div>
  
| Analysis| β | β 95% CI | SE | OR | OR 95% CI | p-value |
|:-------:|:------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| Univariable Risk-Taking → AUD | 1.57 | (0.416, 2.73) | 0.590 | 4.81 | (1.52, 15.3) | 0.00770 |
| Multivariable Risk-Taking → AUD (adjusted for depression) | 0.750 | (-0.620, 2.12) | 0.700 | 2.12 | (0.538, 8.34) | 0.287 |
| Univariable Depression → AUD | 0.411 | (0.188, 0.634) | 0.114 | 1.51 | (1.21, 1.89) | 0.000300 |
| Multivariable Depression → AUD (adjusted for risk-taking) | 0.395 | (0.181, 0.610) | 0.110 | 1.48 | (1.20, 1.84) | 0.000599 |

</div>
<p><i>Table 3: Comparison of Univariable and Multivariable MR Results for Risk-Taking + Depression → Alcohol Use Disorder (all IVW)</i></p>

<br>
<br>

![Univariable vs MVMR — AUD (odds ratios)](Figures/univariable_vs_mv_aud_OR.png)
<p align="center"><i>Figure 12: Univariable vs. MVMR Odds Ratios for Risk-Taking + Depression → Alcohol Use Disorder</i></p>

<br>

### Heterogeneity, Pleiotropy, and Heterogeneity-Robust Estimates

The MVMR Q-statistic was significant, indicating heterogeneity among the 67 joint instruments. The MVMR-Egger intercept, however, was not significant, providing no evidence of directional pleiotropy across the joint set of genetic instruments. To determine if the primary IVW MVMR estimate was compromised by this heterogeneity, a heterogeneity-robust `qhet_mvmr()` estimation was run, both under the default assumption of no correlation between SNP-exposure effects and using an empirically estimated instrument correlation, included within a tested sensitivity range. Results were virtually unchanged across the entire sensitivity range, with depression remaining significant and risk-taking remaining non-significant. This suggests that the IVW MVMR results for AUD are not an artifact of the observed heterogeneity and are not sensitive to a plausible range of instrument correlations.

### Attenuation

The attenuation percent computed in Part 7 of [MVMR_Analysis.R](MVMR_Analysis.R) quantifies the change from each exposure's univariable to multivariable IVW estimate. Depression showed little attenuation, with the beta estimate decreasing by just 3.86% compared to risk-taking's more drastic 52.2% decrease after adjustment. Although the univariable analysis suggested that liability for risk-taking has a causal relationship with AUD risk, multivariable adjustment reveals that the observed effect of liability for risk-taking on AUD may be largely explained by the liability shared with depression.


### Overall Interpretation

The MVMR analysis for AUD risk supports an independent causal relationship of depression on AUD risk, while the univariable association for risk-taking likely reflects genetic liability correlated with depression, rather than an independent pathway. This conclusion is substantiated by maintained conditional instrument strength, no evidence of MVMR-Egger pleiotropy, consistency across heterogeneity-robust estimates, a stable and significant depression effect, and an extremely attenuated and non-significant risk-taking effect. 

The multivariable IVW analysis estimated that each log-odds unit increase in genetic liability to depression was associated with a 1.48-fold higher odds of alcohol use disorder, consistent with the univariable analysis. The analysis did not find a significant effect of liability to risk-taking on alcohol use disorder risk.

---

## 4. Univariable Mendelian Randomization Analysis of the Effect of Liability for Risk-Taking on Non-Alcohol Substance Use Disorder Risk

### Instruments and Directionality

The same 29 risk-taking SNPs were used, with one proxy substitution and no losses during harmonization. Instrument strength, like in the univariable analysis for AUD, was strong throughout. The Steiger directionality test supported the hypothesized causal direction (risk-taking → substance use disorder, excluding alcohol) with high confidence (p = 2.89 * 10⁻³⁷), meaning the genetic instruments explain significantly more variance in risk-taking than in non-alcohol substance use disorders, as expected if risk-taking is upstream. 

### Main Estimate

<br>

<div>
  
| Method| β | β 95% CI | SE | OR | OR 95% CI | p-value |
|:-------:|:------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| IVW | 2.17 | (-0.0334, 4.38) | 1.13 | 8.77 | (0.967, 79.6) | 0.0536 |
| MR Egger | 6.01 | (-4.55, 16.6) | 5.39 | 408 | (0.0106, 1.58 * 10⁷ ) | 0.274 |
| Weighted Median | 3.03 | (0.534, 5.54) | 1.28 | 20.9 | (1.706, 255) | 0.0173 |
| Weighted Mode | 3.07 | (-1.29, 7.43) | 2.22 | 21.6 | (0.277, 1.68 * 10³) | 0.178 |

</div>
<p><i>Table 4: Comparison of Univariable Methods for Risk-Taking → Non-Alcohol Substance Use Disorders</i></p>

<br>

The scatter plot in Figure 13 shows all four MR methods used (IVW, MR-Egger, weighted median, weighted mode) estimating a positive slope, meaning higher genetic liability for risk-taking is associated with a greater risk of substance use disorders (excluding alcohol). Only the weighted median estimate reached significance at p < 0.05, though all methods agreed in direction, a sign of consistency among methods with different assumptions and validity for this dataset. The majority of individual SNP estimates clustered together, with a few points exhibiting smaller effects on non-alcohol substance use disorders or larger effects on risk-taking.

<br>

![Risk-Taking → Substance Use scatter](Figures/risk_sub_scatter_plot.png)
<p align="center"><i>Figure 13: Scatter Plot of Risk-Taking → Non-Alcohol Substance Use Disorders, Across Methods</i></p>

<br>

The method comparison forest plot in Figure 14 shows the weighted median confidence interval fully above the null, with the IVW interval having borderline significance, and MR-Egger / weighted mode having wide and imprecise intervals. The large bounds of the weighted mode are consistent with its lower statistical power and do not necessarily suggest conflicting results with the weighted median and IVW. The immensely large bounds of the MR-Egger interval also reflect low power, as well as heterogeneity, and serve as evidence of the instability of these estimates, not evidence that the true effect is astronomically large.

<br>

![Risk-Taking → Substance Use odds ratio forest plot](Figures/risk_sub_methods_forest_plot.png)
<p align="center"><i>Figure 14: Odds Ratios across MR Methods for Risk-Taking Effects on Non-Alcohol Substance Use Disorders</i></p>

<br>

Although the weighted median was the only significant estimate, due to its low confidence and precision, the remainder of this univariable analysis is computed with IVW methods. The IVW estimate had a large interval as well, but to a smaller extent. Throughout the discussion of this section's results, it is important to note that the IVW estimate is both non-significant and imprecise.

### Heterogeneity and Pleiotropy

Both Cochran's Q and the MR-PRESSO global test were significant here, unlike in the univariable risk-taking-AUD analysis. However, MR-PRESSO did not identify any significant pleiotropic outliers, suggesting that the observed heterogeneity is distributed across the instruments, rather than being concentrated among select variants. This is corroborated by the MR-Egger intercept test, which was not significant, indicating no significant evidence of directional pleiotropy. Although there is not significant directional pleiotropy, the significant heterogeneity warrants cautious interpretation of the estimated effect.

The SNP forest and funnel plots in Figures 15 and 16 show that many SNPs exhibit directional discordance, expected given the significant evidence of heterogeneity via Cochran's Q and MR-PRESSO. The funnel plot is fairly symmetric around the IVW estimate, not indicating any directional pleiotropy, despite variability in the sign and magnitude of many SNPs.

<br>

![Risk-Taking → Substance Use SNP forest plot](Figures/risk_sub_forest_plot.png)
<p align="center"><i>Figure 15: SNP Forest Plot for Risk-Taking → Substance Use Disorders (Excluding Alcohol)</i></p>

<br>
<br>

![Risk-Taking → Substance Use funnel plot](Figures/risk_sub_funnel_plot.png)
<p align="center"><i>Figure 16: SNP Funnel Plot for Risk-Taking → Substance Use Disorders (Excluding Alcohol)</i></p>

<br>

### Leave-One-Out

In Figure 17, the leave-one-out plot indicates that no single SNP disproportionately drives the overall estimate. However, unlike in the univariable risk-taking-AUD analysis, the majority of leave-one-out intervals overlap the null, attesting to the imprecision and unreliability in the risk-taking-substance-use-disorder analysis. 

<br>

![Risk-Taking → Substance Use leave-one-out](Figures/risk_sub_loo_plot.png)
<p align="center"><i>Figure 17: Leave-One-Out Plot for Risk-Taking → Substance Use Disorder (Excluding Alcohol)</i></p>

<br>

### Overall Interpretation

Among all of the univariable analyses, this is the weakest result. In support of a risk-taking effect, the instrument strength and directionality are appropriate, and there is no significant directional pleiotropy. However, extreme heterogeneity, evidenced by Cochran's Q and MR-PRESSO, as well as instability in estimate magnitudes across all methods, are a large concern. Further, most leave-one-out intervals contain the null, indicating that the computed significance of the weighted median is not robust. At best, there may be a positive effect of liability for risk-taking on the risk of non-alcohol substance use disorders. Regardless of whether an effect is present, there is considerable uncertainty and insufficient robust evidence, requiring further investigation, potentially with a larger exposure sample size. 

---

## 5. Univariable Mendelian Randomization Analysis of the Effect of Liability for Depression on Risk for Non-Alcohol Substance Use Disorders

### Instruments and Directionality

The same 48 depression SNPs were used, with 3 proxy substitutions and 2 SNPs omitted during harmonization. Instrument strength, like in the univariable analysis for AUD, was strong throughout. The Steiger directionality test supported the hypothesized causal direction (depression → substance use disorders, excluding alcohol) with high confidence (p = 6.10 * 10⁻⁶⁵), meaning the genetic instruments explain significantly more variance in depression than in non-alcohol substance use disorders, as expected if depression is upstream. 

### Main Estimate

<br>

<div>
  
| Method| β | β 95% CI | SE | OR | OR 95% CI | p-value |
|:-------:|:------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| IVW | 0.406 | (0.123, 0.689) | 0.144 | 1.50 | (1.13, 1.99) | 0.00491 |
| MR Egger | 0.630 | (-1.09, 2.35) | 0.880 | 1.88 | (0.335, 10.5) | 0.478 |
| Weighted Median | 0.403 | (0.0161, 0.790) | 0.197 | 1.50 | (1.02, 2.20) | 0.0412 |
| Weighted Mode | 0.123 | (-0.707, 0.949) | 0.422 | 1.13 | (0.493, 2.58) | 0.776 |

</div>
<p><i>Table 5: Comparison of Univariable Methods for Depression → Non-Alcohol Substance Use Disorders</i></p>

<br>

In Figure 18, all four MR methods used (IVW, MR-Egger, weighted median, weighted mode) estimated a positive slope, meaning higher genetic liability for depression is associated with a greater risk of non-alcohol substance use disorders. The IVW and weighted median estimates both reached significance at p < 0.05, and all methods agreed in direction, a sign of consistency among methods with different assumptions and validity for this dataset. The majority of individual SNP estimates clustered together, with a few points having larger estimates for the effect on depression.

<br>

![Depression → Substance Use scatter](Figures/dep_sub_scatter_plot.png)
<p align="center"><i>Figure 18: Scatter Plot of Depression → Non-Alcohol Substance Use Disorders, Across Methods</i></p>

<br>

Figure 19 compares the methods used, showing that both IVW and the weighted median confidence intervals are above the null, with MR-Egger and the weighted mode having wider, imprecise intervals. The large bounds of both the MR-Egger and weighted mode methods are consistent with their lower statistical power and do not necessarily suggest conflicting results among the four effect estimates.

<br>

![Depression → Substance Use odds ratio forest plot](Figures/dep_sub_methods_forest_plot.png)
<p align="center"><i>Figure 19: Odds Ratios across MR Methods for Depression Effects on Substance Use Disorders (Excluding Alcohol)</i></p>

<br>

### Heterogeneity and Pleiotropy

Cochran's Q, the MR-Egger intercept, and MR-PRESSO were all non-significant, showing no evidence of heterogeneity, pleiotropy, or pleiotropic outliers. Collectively, these results imply that the SNP-specific effect estimates are reasonably consistent, with no notable distorting characteristics. 

Though the global statistics did not identify significant issues, the SNP forest and funnel plots (Figures 20-21) illustrate that a few SNPs exhibit directional discordance. The funnel plot is rather symmetric around the IVW estimate, not indicating any systematic directional pleiotropy, despite some individual SNP variance. 

<br>

![Depression → Substance Use SNP forest plot](Figures/dep_sub_forest_plot.png)
<p align="center"><i>Figure 20: SNP Forest Plot for Depression → Substance Use Disorders (Excluding Alcohol)</i></p>

<br>
<br>

![Depression → Substance Use funnel plot](Figures/dep_sub_funnel_plot.png)
<p align="center"><i>Figure 21: SNP Funnel Plot for Depression → Substance Use Disorders (Excluding Alcohol)</i></p>

<br>

### SNP-Level Sensitivity: rs843812

Similar to `rs279846` in risk-taking-AUD, one SNP, `rs843812`, noticeably deviated from other points in the forest plot, potentially skewing the overall estimate to be lower. A sensitivity analysis was conducted:

- **IVW estimate with rs843812 included:** β = 0.406, Odds Ratio = 1.50, p = 0.00491
- **IVW estimate with rs843812 excluded:** β = 0.460, Odds Ratio = 1.58, p = 0.00103
- The SNP-specific F-statistic is 31.8, ruling out weak-instrument artifacts as the explanation for its divergence.

Once `rs843812` was removed, the IVW estimate direction was unchanged, with a slightly increased magnitude and significance. Because neither Cochran's Q or MR-PRESSO flagged significant heterogeneity or pleiotropy issues, `rs843812` was retained in the primary analysis.

### Leave-One-Out

In Figure 22, the leave-one-out plot shows that no single SNP disproportionately drives the estimate, and all confidence intervals are above the null with high consistency among estimates. Removing `rs843812` increases the overall estimate more than removing any other SNP, consistent with its sensitivity checks.

<br>

![Depression → Substance Use leave-one-out](Figures/dep_sub_loo_plot.png)
<p align="center"><i>Figure 22: Leave-One-Out Plot for Depression → Substance Use Disorders (Excluding Alcohol)</i></p>

<br>

### Overall Interpretation

These results collectively provide consistent evidence of a positive causal effect of genetic liability for depression on substance use disorders (excluding alcohol). Both the IVW and weighted median estimates were significant, supported by strong genetic instruments and correct Steiger directionality. Sensitivity analyses showed no significant evidence of heterogeneity or directional pleiotropy, and all four methods were directionally consistent, regardless of significance.

The univariable IVW analysis estimated that each log-odds unit increase in genetic liability to depression was associated with a 1.51-fold higher odds of non-alcohol substance use disorders, consistent with a potential causal effect of depression liability on the risk of non-alcohol substance use disorders, under the assumptions of Mendelian randomization. Interestingly, this estimated effect is quite similar to the calculated 1.50-fold higher odds of AUD for each log-odds unit increase in depression liability. 

---

## 6. Multivariable Mendelian Randomization Analysis on the Effect of Liability to Risk-Taking and Depression on Risk for Non-Alcohol Substance Use Disorders

This analysis aims to understand whether liability for risk-taking and depression have independent effects on risk for substance use disorders (excluding alcohol), or if one of these exposures is capturing their shared genetic liability.

### Instrument Overlap

The same 67 unique SNPs from the multivariable analysis on AUD were used. Again, there were no overlapping SNPs between risk-taking and depression instrument sets, suggesting any observed attenuation reflects shared underlying genetic liability or correlated pathways.

### Conditional Instrument Strength

Again, both exposures maintained adequate conditional F-statistics (F > 10), with the conditional F-statistic for depression being noticeably larger than that of risk-taking. This indicates that once the analysis accounts for correlation between the two exposures, the risk-taking instruments retain comparatively less unique predictive power.

### Main MVMR Estimate

After mutual adjustment, the liability for depression retains a significant independent association with substance use disorders (excluding alcohol), while the estimate for risk-taking lost significance, with a much wider confidence interval (Figure 23). However, both point estimates remained positive, in alignment with the univariable results. More details on the specific estimates are outlined in the next subsection.

<br>

![MVMR Substance Use forest plot (OR)](Figures/mv_sub_forest_OR_plot.png)
<p align="center"><i>Figure 23: MVMR Odds Ratios for Risk-Taking + Depression → Substance Use Disorders (Excluding Alcohol)</i></p>

<br>

### Direct Comparison: Univariable vs. MVMR-Adjusted Estimates

As detailed in Table 6 and Figure 24, the point estimate for depression barely moves between the univariable and MVMR-adjusted models, while the estimate for risk-taking shifts substantially toward the null. 

<br>

<div>
  
| Analysis| β | β 95% CI | SE | OR | OR 95% CI | p-value |
|:-------:|:------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| Univariable Risk-Taking → Substance Use Disorders, Excluding Alcohol | 2.17 | (-0.0334, 4.38) | 1.13 | 8.77 | (0.967, 79.6) | 0.0536 |
| Multivariable Risk-Taking → Substance Use Disorders, Excluding Alcohol (adjusted for depression) | 0.543 | (-1.47, 2.55) | 1.03 | 1.72 | (0.230, 12.9) | 0.598 |
| Univariable Depression → Substance Use Disorders, Excluding Alcohol | 0.406 | (0.123, 0.689) | 0.144 | 1.50 | (1.13, 1.99) | 0.00491 |
| Multivariable Depression → Substance Use Disorders, Excluding Alcohol (adjusted for risk-taking) | 0.451 | (0.136, 0.766) | 0.161 | 1.57 | (1.15, 2.15) | 0.00660 |

</div>
<p><i>Table 6: Comparison of Univariable & Multivariable MR for Risk-Taking + Depression → Non-Alcohol Substance Use Disorders (all IVW)</i></p>

<br>
<br>

![Univariable vs MVMR — Substance Use (odds ratios)](Figures/univariable_vs_mv_sub_OR.png)
<p align="center"><i>Figure 24: Univariable vs. MVMR Odds Ratios for Risk-Taking + Depression → Non-Alcohol Substance Use Disorders</i></p>

<br>

### Heterogeneity, Pleiotropy, and Heterogeneity-Robust Estimates

As with the AUD MVMR analysis, this MVMR Q-statistic was significant, but the MVMR-Egger intercept was not significant. This suggests heterogeneity among the 67 joint instruments, but does not provide evidence of directional pleiotropy. To determine if the primary IVW MVMR estimate was compromised by the observed heterogeneity, a heterogeneity-robust `qhet_mvmr()` estimation was run, both under the default assumption of no correlation between SNP-exposure effects and with empirically estimated instrument correlation, included within a tested sensitivity range. This analysis produced results highly similar to the primary MVMR analysis. Genetic liability for depression remained significantly associated with increased odds of substance use disorders (excluding alcohol), while the independent effect of risk-taking remained non-significant. This suggests that the IVW MVMR results for substance use disorders (excluding alcohol) are not an artifact of the observed heterogeneity and are not sensitive to a plausible range of instrument correlations.

### Attenuation

The attenuation percent computed in Part 7 of [MVMR_Analysis.R](MVMR_Analysis.R) quantifies the change from each exposure's univariable to multivariable IVW estimate. Depression showed no attenuation, with the beta estimate actually increasing by 11.1% compared to risk-taking's extreme 75.0% decrease after adjustment. Although the univariable analysis suggested that liability for risk-taking has a causal relationship with the risk for non-alcohol substance use disorders, multivariable adjustment reveals that the observed effect of liability for risk-taking may be largely explained by the liability shared with depression.


### Interpretation

The MVMR analysis supports an independent causal relationship of depression on risk of non-alcohol substance use disorders, while the univariable association for risk-taking likely reflects genetic liability correlated with depression, rather than an independent pathway. This conclusion is substantiated by maintained conditional instrument strength, no evidence of MVMR-Egger pleiotropy, consistency across heterogeneity-robust estimates, a stable and significant depression effect, and an extremely attenuated and non-significant risk-taking effect. 

The multivariable IVW analysis estimated that each log-odds unit increase in genetic liability to depression was associated with a 1.57-fold higher odds of non-alcohol substance use disorders, consistent with the univariable analysis. The analysis did not find a significant effect from risk-taking liability.

Notably, the computed multivariable effect of depression liability on risk of non-alcohol substance use disorders (OR = 1.57) was strikingly similar to the univariable effect of depression liability on risk of non-alcohol substance use disorders (OR = 1.50), as well as the computed effects of depression liability on AUD risk (Univariable OR = 1.51, Multivariable OR = 1.48).

---

## 7. Cross-Outcome and Cross-Exposure Comparisons

### Instrument Overlap

Zero SNPs overlapped between the risk-taking and depression instrument SNPs, confirming that the two exposures are being assessed with distinct sets. This strengthens the interpretation of the MVMR attenuation results, as the exposure adjustments are separating distinct, correlated genetic sets, not simply adjusting for the same SNPs.

### Genome-Wide Genetic Correlation: Linkage Disequilibrium Score Regression (LDSC)

To understand whether the MVMR attenuation pattern is consistent with shared genetic liability between depression and risk-taking, bivariate linkage disequilibrium score regression (LDSC, via GenomicSEM) was run on the full depression and risk-taking GWAS summary statistics (1,197,272 overlapping SNPs after merging with the European LD reference panel). Both depression and risk-taking showed strong evidence of nonzero SNP heritability (h² Z = 26.1 for depression, h² Z = 21.8 for risk-taking). The LDSC intercept was approximately 1 for depression (0.996) and slightly higher for risk-taking (1.038), suggesting significant inflation from confounding is not likely.

LDSC estimated a significant positive genome-wide correlation between depression and risk-taking (rg = 0.150, SE = 0.0268, p = 2.04 * 10⁻⁸). Note that this correlation is a distinct value from `pcor ≈ 0.19` used in `qhet_mvmr()` for the previously discussed sensitivity analyses. The LDSC correlation (rg) describes the correlation essentially across the full genome, while the pcor value describes the correlation across the 67 evaluated SNPs.

The LDSC results provide direct and independent support for the MVMR attenuation findings detailed above, indicating that depression and risk-taking share a portion of their genetic makeup. The shared genetic makeup is consistent with the observation of risk-taking attenuation between the univariable and multivariable associations with both AUD and non-alcohol substance use disorders.

### SNP Effect Correlation between AUD and Substance Use Disorders (Excluding Alcohol)

For both depression and risk-taking, the SNP-level univariable effects on risk for AUD and substance use disorders (excluding alcohol) were moderately positively correlated (Figures 25-26), with most SNPs showing the same effect direction on both outcomes. For depression, there was slightly more discordance in the effect directions between the exposures, compared to that of risk-taking. Based on the univariable analyses, if an exposure's genetic instruments influence AUD risk, it would be reasonable to expect a similar influence on the risk for non-alcohol substance use disorders (and vice versa). This does not imply that the two outcomes experience the same effects, just that the exposure effects have similar patterns. This correlation is especially notable in the univariable and multivariable depression analyses, which produced very similar estimates for the effect of depression on both outcomes. The pattern for risk-taking is more challenging to articulate, given its attenuation when accounting for depression.

<br>

![Risk-taking SNP effects across outcomes](Figures/risk_comparison_plot.png)
<p align="center"><i>Figure 25: Risk-Taking Effects Across Outcomes (Univariable)</i></p>

<br>
<br>

![Depression SNP effects across outcomes](Figures/dep_comparison_plot.png)
<p align="center"><i>Figure 26: Depression Effects Across Outcomes (Univariable) </i></p>

<br>

### Overall Picture of Multivariable and Univariable Results

Across both outcomes, the univariable and MVMR-adjusted estimates for depression are nearly identical, with small CIs and minimal changes in the point estimates. In contrast, the univariable and MVMR-adjusted estimates across outcomes for risk-taking feature wide CIs with the point estimate shrinking toward the null after adjusting for depression. All eight point estimates are directionally positive, but only depression is consistently significant across these methods. See Figure 27.

<br>

![Combined forest plot: all four exposure-outcome pairs](Figures/combined_forest_plot.png)
<p align="center"><i>Figure 27: Univariable vs. Multivariable-Adjusted Odds Ratios</i></p>

<br>

---

## 8. Power Analysis

To contextualize the primary IVW findings of these analyses, the statistical power for detecting a binary MR effect was assessed for each of the four exposure-outcome combinations for both univariable and multivariable analyses, using the standard approach based on outcome sample size (N), outcome case proportion (K), and the variance in the exposure explained by its instruments (R²). A minimum power of 80% was used as the threshold. Note that the multivariable power values were computed using the univariable R² values as a proxy for conditional R² values, a key limitation. 

<br>

<div>

| Analysis | Observed OR | p-value | Observed Power | Minimum Detectable OR (80% power) |
|:---:|:---:|:---:|:---:|:---:|
| Risk-Taking → AUD | 4.81 | 0.00770 | 1.00 | 1.82 |
| Depression → AUD | 1.51 | 0.000300 | 0.637 | 1.65 |
| Risk-Taking → Substance Use | 8.77 | 0.0536 | 1.00 | 2.46 |
| Depression → Substance Use | 1.50 | 0.00491 | 0.326 | 2.12 |
| MVMR: Risk-Taking → AUD | 2.12 | 0.287 | 0.941 | 1.82 |
| MVMR: Depression → AUD | 1.48 | 0.000599 | 0.603 | 1.65 |
| MVMR: Risk-Taking → Substance Use | 1.72 | 0.598 | 0.393 | 2.46 |
| MVMR: Depression → Substance Use | 1.57 | 0.00660 | 0.388 | 2.12 |

</div>
<p><i>Table 7: Observed Power and Minimum Detectable OR (80% Power) for Univariable and Multivariable MR Results (all IVW)</i></p>

<br>

The minimum detectable OR at 80% power is computed independently of the observed effect size, using the instrument strength and outcome sample characteristics. All four depression analyses were significant despite being below the minimum detectable OR at 80% power. However, this is not contradictory, as the power reflects the probability of detecting a certain effect size, not a cutoff. Significance below 80% power is less probable but still valid. However, this does raise the concern of potential winner's curse, where significant findings in lower-powered analyses tend to overestimate the true effect size. These results should be interpreted with some caution and ideally replicated in larger samples.  

Conversely, three of the four risk-taking ORs were above the minimum detectable OR at 80% power, despite two of these results not having significant results. These two outcomes had large point estimates and wide confidence intervals, suggesting imprecision rather than a lack of effect. The high variability and attenuation across risk-taking instruments, as well as instances of heterogeneity and pleiotropy, likely contributed to inflated, imprecise estimates, leading to large observed power values.

Overall, the power estimates here should be treated cautiously and as a rough guide rather than a definitive measure of analysis credibility. Because power is probabilistic and not a rigid threshold, both univariable and multivariable depression results were able to reach significance despite falling below the minimum detectable OR. Risk-taking, however, had imprecise, non-significant results despite large point estimates that exceeded the minimum detectable OR. For the multivariable analyses specifically, there is an additional limitation of using the univariable R² values as proxies for conditional R² values. Independent replication with larger samples is the clear next step in validating these findings.

---

## Overall Discussion

Across both outcomes, this analysis supports a consistent narrative:

1. **Depression liability shows evidence of an independent causal effect** on risk of both AUD and non-alcohol substance use disorders. This effect is directionally correct (Steiger testing), estimated with strong instruments, significant across both univariable and multivariable methods, largely free of directional pleiotropy, and is virtually unaffected by the MVMR adjustment for risk-taking, even under a range of heterogeneity-robust re-estimation methods.

2. **Univariable associations of risk-taking liability with AUD and non-alcohol substance use disorder risk are less secure**, particularly for substance use disorders. Significant Cochran's Q and MR-PRESSO results, considerable variability among methods, and leave-one-out intervals overlapping the null all raise questions on the effect's precision (but not necessarily for the effect direction). More importantly, **both outcome associations with risk-taking largely attenuate after adjustments for depression**, dropping over 50% for AUD and nearly 75% for non-alcohol substance use disorders, leading to the loss of statistical significance in both MVMR models.

3. **Collectively, more evidence supports depression having a genuine independent causal pathway to alcohol and substance outcomes**, with the apparent univariable effects for risk-taking being significantly, but not completely, attributable to genetic liability shared with and/or mediated through depression.

4. **The potential for shared liability between depression and risk-taking is supported by genome-wide LDSC**, which found a significant positive genetic correlation between depression and risk-taking (rg = 0.150, p = 2.04 * 10⁻⁸) using the full GWAS summary statistics for each trait. Additionally, the SNP-level correlation used in the `qhet` sensitivity analysis (pcor ≈ 0.19) and the genome-wide LDSC estimate (rg = 0.150) have similar magnitudes, further supporting the rationale that MVMR risk-taking attenuation reflects shared genetic makeup between the two exposures, rather than being merely a feature of the specific 67 SNPs jointly evaluated.

## Limitations

- **Heterogeneity in risk-taking analyses**, especially for non-alcohol substance use disorders, was not resolved by identifying specific outliers. This heterogeneity may reflect genuine biological pleiotropy distributed across instruments or differences in underlying causal pathways. Consequently, the magnitude (and potentially the direction) of the estimated causal effects may be biased.
- **No multiple-testing correction** was used across the four primary IVW comparisons. However, all IVW p-values were sufficiently low that a Bonferroni correction for a 0.05/4 = 0.0125 significance level would not change any of the core results.
- **The `pcor` parameter in `qhet_mvmr()`** was approximated using instrument-level correlation, not phenotypic correlation between depression and risk-taking. A direct genome-wide LDSC estimate (rg = 0.150, p = 2.04 * 10⁻⁸) was obtained separately and was generally consistent in direction and magnitude with the instrument-level estimate. However, both are genetic, not phenotypic, and the true phenotypic correlation cannot be truly captured by either estimate.
- **Sample overlap between the UK Biobank–derived exposures and the FinnGen outcomes** is expected to be minimal, as these are independently recruited (UK vs. Finland). However, no formal quantification of the overlap was computed, so overlap bias cannot be ruled out.
- **Both exposure datasets reflect European populations, and both outcome datasets reflect Finnish populations**, which makes it challenging to generalize these results more globally. Ideally, these analyses would be performed on larger, global datasets to reduce the risk of geographic and demographic confounding (though this may violate MR assumptions).
- **Power calculations** summarize instrument strength as an overall R² value instead of accounting for the individual variation in SNP strength. Due to this limitation, the computed minimum detectable ORs at 80% power should be treated as rough benchmarks and not precise measurements.
- **Risk-taking is quantified as a single, binary, self-reported trait** (*"would you describe yourself as someone who takes risks?"*). This is a much broader exposure measure than a dataset with proven risk-taking behavior, which likely contributes to weaker instrument strength and the imprecision flagged in the risk-taking results.
- **Depression (log-odds) and risk-taking (standard deviation scale)** were derived using different approaches, so a single unit change in genetically predicted liability does not reflect an equivalent magnitude. Significance and attenuation are unaffected by this difference, but direct comparisons of effect magnitude between the exposures should account for scaling differences.
- **The outcome datasets have relatively small case counts**, compared to the exposure datasets. Small case counts can reduce precision and increase vulnerability to winner's curse, which was a potential concern.
- **Core MR assumptions (relevance, independence, exclusion restriction) cannot be fully verified.** Relevance was supported by strong F-statistics, but independence and exclusion restriction are not directly testable and only partially addressed (MR-Egger, MR-PRESSO, etc.).

---

## References

[1] ieu-b-102 - Major Depression. (n.d.). https://opengwas.io/datasets/ieu-b-102 

[2] ukb-b-14147 - Risk-Taking. (n.d.). https://opengwas.io/datasets/ukb-b-14147

[3] finn-b-AUD - Alcohol Use Disorder: ICD-Based. (n.d.). https://opengwas.io/datasets/finn-b-AUD

[4] finn-b-F5_SUBSNOALCO - Substance Use, Excluding Alcohol. (n.d.). https://opengwas.io/datasets/finn-b-F5_SUBSNOALCO

[5] European LD scores from 1000 genomes | Zenodo. (n.d.). https://zenodo.org/records/8182036 

[6] Howard et al. (2019). Genome-wide meta-analysis of depression identifies 102 independent variants and highlights the importance of the prefrontal brain regions. Nature neuroscience, 22(3), 343–352. https://doi.org/10.1038/s41593-018-0326-7

[7] Calarco, C. A., & Lobo, M. K. (2021). Depression and substance use disorders: Clinical comorbidity and shared neurobiology. International review of neurobiology, 157, 245–309. https://doi.org/10.1016/bs.irn.2020.09.004

[8] Strawbridge, R. J., Ward, J., Cullen, B., Tunbridge, E. M., Hartz, S., Bierut, L., Horton, A., Bailey, M. E. S., Graham, N., Ferguson, A., Lyall, D. M., Mackay, D., Pidgeon, L. M., Cavanagh, J., Pell, J. P., O'Donovan, M., Escott-Price, V., Harrison, P. J., & Smith, D. J. (2018). Genome-wide analysis of self-reported risk-taking behaviour and cross-disorder genetic correlations in the UK Biobank cohort. Translational psychiatry, 8(1), 39. https://doi.org/10.1038/s41398-017-0079-1

[9] Poore HE, Hatoum A, Mallard TT, et al. A multivariate approach to understanding the genetic overlap between externalizing phenotypes and substance use disorders. Addiction Biology. 2023; 28(9):e13319. doi:10.1111/adb.13319

[10] Alcohol use disorder: From risk to diagnosis to Recovery | National Institute on Alcohol Abuse and alcoholism (NIAAA). (n.d.). https://www.niaaa.nih.gov/health-professionals-communities/core-resource-on-alcohol/alcohol-use-disorder-risk-diagnosis-recovery 

[11] Alcohol. (n.d.). https://www.who.int/news-room/fact-sheets/detail/alcohol 

[12] Li, D., Sulovari, A., Cheng, C. et al. Association of Gamma-Aminobutyric Acid A Receptor α2 Gene (GABRA2) with Alcohol Use Disorder. Neuropsychopharmacol 39, 907–918 (2014). https://doi.org/10.1038/npp.2013.291

Also helpful: https://mr-dictionary.mrcieu.ac.uk/

---

## About

Inspired by my friend, Natalie, who is pursuing a career in clinical psychology focused on addiction and substance use, this project is my attempt to explore these topics from a more quantitative perspective. Genetic data was used to decipher whether depression and risk-taking actually play a causal role in alcohol and substance use disorders, or whether the connections often assumed are simply assumptions. Big shout-out to Natalie for her much-needed advice throughout this project 😎!!!
