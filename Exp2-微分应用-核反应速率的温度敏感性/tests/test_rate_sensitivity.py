import math
import matplotlib.pyplot as plt
import numpy as np


def reaction_rate(T):
    """
    编写一个函数，输入参数温度 T (单位 K)，返回 3-α 反应速率中与温度相关的部分 q_3a / (ρ²γ³)
    即计算: f(T) = 5.09×10¹¹T₈⁻³ e⁻⁴⁴.⁰²⁷/T₈，其中 T₈ = T / 10⁸。
    """
    T8 = T / 1e8
    return 5.09 * (10 ** 11) * (T8 ** (-3)) * math.exp(-44.027 / T8)


def temperature_sensitivity(T0, h=1e-8):
    """
    根据 v 的定义 v = (T/q) * (dq/dT)，使用数值微分（向前差分法）来近似导数 dq/dT
    选择一个合适的微小温度扰动 ΔT。为了平衡截断误差和舍入误差，一个常见的选择是 ΔT = h * T0，其中 h 是一个很小的数，例如 h = 10⁻⁸。
    将近似导数代入 v 的公式。
    """
    dT = h * T0
    q_T0 = reaction_rate(T0)
    q_T0_plus_dT = reaction_rate(T0 + dT)
    dq_dT = (q_T0_plus_dT - q_T0) / dT
    return (T0 / q_T0) * dq_dT


# 任务3
reference_temperatures = [1.0e8, 2.5e8, 5.0e8, 1.0e9, 2.5e9, 5.0e9]
for T0 in reference_temperatures:
    v = temperature_sensitivity(T0)
    print(f"温度: {T0} K, v值: {v}")

# 任务4
# 生成温度范围
T = np.logspace(8, 9, 100)  # 从10^8 到 10^9 取100个点
q_values = [reaction_rate(t) for t in T]

# 绘制对数-对数坐标图
plt.loglog(T, q_values)
plt.xlabel('温度 T (K)')
plt.ylabel('q_3a / (ρ²γ³)')
plt.title('反应速率与温度的关系')
plt.show()
