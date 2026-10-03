import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, RealData

# 線形モデル
def linear_model(params, x):
    a, b = params
    return a * x + b

# サンプルデータ生成（直線）
np.random.seed(42)
x_true = np.linspace(-10, 10, 20)
params_true = [3.0, 2.0]
y_true = linear_model(params_true, x_true)

# ノイズを加えた観測データ
x_obs = x_true + np.random.normal(0, 0.5, size=x_true.size)
y_obs = y_true + np.random.normal(0, 3.0, size=y_true.size)
x_err = np.full_like(x_obs, 0.5)  # x の誤差
y_err = np.full_like(y_obs, 3.0)  # y の誤差

# ODR 用データとモデルの設定
data = RealData(x_obs, y_obs, sx=x_err, sy=y_err)
model = Model(linear_model)

# ODR 実行
odr_instance = ODR(data, model, beta0=[10.0, 10.0])
output = odr_instance.run()

# フィッティング結果
popt = output.beta
perr = output.sd_beta
print(f"True param: y = {params_true[0]:.2f} x + {params_true[1]:.2f}")
print(f"Fitting result: y = {popt[0]:.2f} x + {popt[1]:.2f}")
print(f"Parameter uncertainty: a_err={perr[0]:.2f}, b_err={perr[1]:.2f}")

# 信頼区間の計算
nstd = 1.96  # 95% 信頼区間
x_fit = np.linspace(x_obs.min(), x_obs.max(), 100)
y_fit = linear_model(popt, x_fit)

# 答えの曲線をプロット
y_true_fit = linear_model(params_true, x_fit)

# 予測区間の計算
y_fit_err = np.sqrt((x_fit**2) * perr[0]**2 + perr[1]**2)
y_fit_up = y_fit + nstd * y_fit_err
y_fit_dw = y_fit - nstd * y_fit_err

# プロット
plt.errorbar(x_obs, y_obs, xerr=x_err, yerr=y_err, fmt='o', label="Data")
plt.plot(x_fit, y_fit, label=f"Fitted: y = {popt[0]:.2f} x + {popt[1]:.2f}", color='r')
plt.plot(x_fit, y_true_fit, label="True line", linestyle='--', color='k')
plt.fill_between(x_fit, y_fit_up, y_fit_dw, color='gray', alpha=0.3, label='95% CI')
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()
