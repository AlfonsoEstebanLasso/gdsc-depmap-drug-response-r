# Regression: mean and 95 percent interval over the folds

Cells show the mean over the 15 held-out folds and a two-sided 95 percent t-based interval (mean plus or minus qt(0.975, n - 1) times sd over sqrt(n)). Folds of different repeats share cell lines, so the interval is descriptive, not an exact inferential statement.

## R squared

| Model | Erlotinib | Rapamycin | Sunitinib | Paclitaxel |
|---|---|---|---|---|
| rf | 0.087 (0.013 to 0.161) | 0.111 (0.071 to 0.152) | 0.059 (-0.093 to 0.211) | 0.053 (0.016 to 0.090) |
| ridge | 0.120 (0.045 to 0.196) | 0.087 (0.049 to 0.125) | 0.095 (-0.003 to 0.193) | 0.043 (0.013 to 0.073) |
| lasso | 0.095 (-0.009 to 0.199) | 0.093 (0.051 to 0.135) | 0.023 (-0.124 to 0.170) | 0.029 (0.005 to 0.053) |
| enet | 0.115 (0.024 to 0.206) | 0.088 (0.047 to 0.128) | 0.063 (-0.062 to 0.189) | 0.042 (0.014 to 0.070) |
| svm_linear | 0.137 (0.075 to 0.199) | 0.009 (-0.048 to 0.067) | -0.012 (-0.119 to 0.095) | -0.103 (-0.153 to -0.054) |

## RMSE

| Model | Erlotinib | Rapamycin | Sunitinib | Paclitaxel |
|---|---|---|---|---|
| rf | 0.074 (0.061 to 0.088) | 0.207 (0.195 to 0.219) | 0.143 (0.132 to 0.154) | 0.193 (0.189 to 0.198) |
| ridge | 0.073 (0.059 to 0.087) | 0.210 (0.198 to 0.222) | 0.142 (0.130 to 0.154) | 0.195 (0.190 to 0.199) |
| lasso | 0.074 (0.060 to 0.087) | 0.209 (0.196 to 0.223) | 0.147 (0.135 to 0.159) | 0.196 (0.191 to 0.201) |
| enet | 0.073 (0.060 to 0.087) | 0.210 (0.197 to 0.223) | 0.144 (0.132 to 0.156) | 0.195 (0.190 to 0.199) |
| svm_linear | 0.073 (0.059 to 0.087) | 0.218 (0.206 to 0.231) | 0.151 (0.137 to 0.164) | 0.209 (0.205 to 0.213) |

## MAE

| Model | Erlotinib | Rapamycin | Sunitinib | Paclitaxel |
|---|---|---|---|---|
| rf | 0.041 (0.038 to 0.045) | 0.152 (0.143 to 0.162) | 0.097 (0.091 to 0.102) | 0.165 (0.160 to 0.170) |
| ridge | 0.040 (0.037 to 0.044) | 0.152 (0.143 to 0.160) | 0.094 (0.090 to 0.099) | 0.167 (0.162 to 0.171) |
| lasso | 0.041 (0.038 to 0.044) | 0.152 (0.143 to 0.161) | 0.099 (0.093 to 0.104) | 0.168 (0.164 to 0.171) |
| enet | 0.041 (0.037 to 0.044) | 0.153 (0.144 to 0.161) | 0.096 (0.091 to 0.101) | 0.166 (0.162 to 0.170) |
| svm_linear | 0.038 (0.034 to 0.041) | 0.144 (0.135 to 0.153) | 0.097 (0.092 to 0.103) | 0.172 (0.168 to 0.176) |

