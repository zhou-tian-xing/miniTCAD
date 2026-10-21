"""基本物理常数。

搬运自 C2_v1/constants.py，改成带单位后缀的名字。
epsilon_0 用 F/cm（SI 的 F/m 数值不能直接代进 cm 公式）。
"""

from __future__ import annotations

# 元电荷 [C]
q = 1.602176634e-19

# Boltzmann 常数 [J/K]
k_b_j_per_k = 1.380649e-23

# Boltzmann 常数 [eV/K]
k_b_ev_per_k = 8.617333262e-5

# 真空介电常数 [F/cm]
epsilon_0_f_per_cm = 8.8541878128e-14

# 真空介电常数 [F/m]（只在适配层核对用）
epsilon_0_f_per_m = 8.8541878128e-12
