## BMI at birth

| Name | # Children | # Mothers | # Fathers | # Total |
| ---- | ---------- | --------- | --------- | ------- |
| bmi_birth | 78093 | 73475 | 51321 | 202889 |
| z_bmi_birth | 78070 | 73453 | 51310 | 202833 |

- Formula: `bmi_birth ~ fp(pregnancy_duration_1)`
- Sigma formula: ` ~ pregnancy_duration_1`
- Distribution: `LOGNO`
- Normalization: `centiles.pred` Z-scores
![](plots/z_bmi_birth_vs_bmi_birth_child.png)


