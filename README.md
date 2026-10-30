A factor-augmented New Keynesian Phillips Curve proposed by Norkute and Westerlund (2024)[^1], based on the NSYC estimator by Norkutė et al. (2021)[^2], is used to estimate the hybrid NKPC method using Chinese month-on-month provincial inflation data. The data include two external variables—real estate investment and the output gap; the estimates below use the output gap. Because the Chinese government does not publish monthly data of compensation of employees and the labour cost index, real marginal cost and wage inflation at the monthly frequency are not available for the given period. The code provides IV estimates and their corresponding covariance matrices, as well as a function for the overidentifying restriction J-test. The details of the methodology closely follow Norkute and Westerlund (2024). The output-gap specification yields a negative forward-looking inflation coefficient and an output-gap coefficient that is not statistically significant at the 5% level.

Second-step IV estimates using four regressor factors and one residual factor. The estimation sample covers 31 provinces from 2013Q2 to 2023Q3 (1,302 observations).

| Variable | Coefficient | Standard error | 95% confidence interval |
| --- | ---: | ---: | ---: |
| Lagged inflation | -0.3153 | 0.2090 | [-0.7249, 0.0944] |
| Forward-looking inflation | -0.6885 | 0.3030 | [-1.2823, -0.0947] |
| Output gap | 1.1699 | 1.0002 | [-0.7903, 3.1302] |

[^1] Norkute, M. and Westerlund, J. (2024), A Factor-Augmented New Keynesian Phillips Curve for the European Union Countries. Oxf Bull Econ Stat, 86: 794-810. https://doi.org/10.1111/obes.12614

[^2] Norkutė, M., Sarafidis, V., Yamagata, T. and Cui, G. (2021). ‘Instrumental variable estimation of dynamic linear panel data models with defactored regressors and a multifactor error structure’, Journal of Econometrics, Vol. 220, pp. 416–446.
