A factor-augmented New Keynesian Phillips Curve proposed by Norkute and Westerlund (2024)[^1], based on the NSYC estimator by Norkutė et al. (2021)[^2], is used to estimate the hybrid NKPC method using Chinese month-on-month provincial inflation data. The data include two external variables—real estate investment and the output gap; the estimates below use the output gap. Because the Chinese government does not publish monthly data of compensation of employees and the labour cost index, real marginal cost and wage inflation at the monthly frequency are not available for the given period. The code provides IV estimates, their covariance matrices and an overidentification statistic. The details of the methodology closely follow Norkute and Westerlund (2024). The forward-looking inflation coefficient is significantly positive, while the lagged inflation and output-gap coefficients are not statistically significant at the 5% level.

| Variable | Coefficient | Standard error | 95% confidence interval |
| --- | ---: | ---: | ---: |
| Lagged inflation | 0.0477 | 0.1384 | [-0.2236, 0.3189] |
| Forward-looking inflation | 0.5352 | 0.2333 | [0.0781, 0.9924] |
| Output gap | -0.3384 | 0.3068 | [-0.9397, 0.2630] |

[^1] Norkute, M. and Westerlund, J. (2024), A Factor-Augmented New Keynesian Phillips Curve for the European Union Countries. Oxf Bull Econ Stat, 86: 794-810. https://doi.org/10.1111/obes.12614

[^2] Norkutė, M., Sarafidis, V., Yamagata, T. and Cui, G. (2021). ‘Instrumental variable estimation of dynamic linear panel data models with defactored regressors and a multifactor error structure’, Journal of Econometrics, Vol. 220, pp. 416–446.
