# 半导体器件建模、典型参数与归一化约定

## 1. 基本约定

- 若无特别说明，材料默认指 **晶体硅（Si）**。
- 温度记为 $T$，**不指定固定数值**。
- 半导体器件中的长度、浓度、迁移率、电流密度等，默认采用常见的 **cm 制工程单位**。
- 所有与温度有关的材料参数原则上写成 $T$ 的函数。
- 为建立数量级直觉，可同时给出常用的 **300 K 参考值**，但这不意味着后续计算默认固定 $T=300\ \mathrm{K}$。
- 归一化参考量由 $T$、$n_i(T)$、$D_0$ 等共同决定，因此原则上保留解析表达式。

---

## 2. 基本物理常数

### 2.1 元电荷

$$
\boxed{
q=1.602176634\times10^{-19}\ \mathrm{C}
}
$$

### 2.2 Boltzmann 常数

$$
\boxed{
k_B=1.380649\times10^{-23}\ \mathrm{J/K}
}
$$

也常写为

$$
\boxed{
k_B=8.617333262\times10^{-5}\ \mathrm{eV/K}
}
$$

### 2.3 真空介电常数

采用 cm 制时，

$$
\boxed{
\varepsilon_0
\approx
8.854\times10^{-14}\ \mathrm{F/cm}
}
$$

SI 制中，

$$
\varepsilon_0
\approx
8.854\times10^{-12}\ \mathrm{F/m}.
$$

---

## 3. 硅的典型材料参数

### 3.1 相对介电常数

晶体硅通常取

$$
\boxed{
\varepsilon_{r,\mathrm{Si}}\approx11.7
}
$$

因此

$$
\varepsilon_s
=
\varepsilon_{\mathrm{Si}}
=
\varepsilon_{r,\mathrm{Si}}\varepsilon_0,
$$

即

$$
\boxed{
\varepsilon_{\mathrm{Si}}
\approx
1.04\times10^{-12}\ \mathrm{F/cm}
}
$$

### 3.2 本征载流子浓度

本征载流子浓度写为

$$
\boxed{
n_i=n_i(T)
}
$$

它对温度非常敏感。

在室温附近，硅通常具有

$$
n_i\sim10^{10}\ \mathrm{cm^{-3}}.
$$

不同材料模型和参数集给出的 300 K 数值有所不同，常见参考值约为

$$
\boxed{
n_i(300\ \mathrm{K})
\sim
(1.0\text{--}1.5)\times10^{10}\ \mathrm{cm^{-3}}
}
$$

因此，在只讨论数量级时可记为

$$
\boxed{
n_i(300\ \mathrm{K})\sim10^{10}\ \mathrm{cm^{-3}}.
}
$$

---

## 4. 热电压

热电压定义为

$$
\boxed{
V_T(T)=\frac{k_B T}{q}
}
$$

单位为 V。

由于 $T$ 不固定，因此后续计算原则上保留 $V_T(T)$。

作为参考，在 300 K 时，

$$
\boxed{
V_T(300\ \mathrm{K})
\approx
25.85\ \mathrm{mV}
}
$$

即约为

$$
V_T\sim 26\ \mathrm{mV}.
$$

---

## 5. 变量、常用单位及典型数量级

下表中的典型值主要用于建立硅器件中的数量级直觉，并不是所有器件都必须落在这些范围内。

| 变量 | 符号 | 常用单位 | 典型值或数量级 |
|---|---|---|---|
| 空间坐标 | $x$ | cm、$\mu\mathrm{m}$、nm | nm 到数百 $\mu\mathrm{m}$，视器件而定 |
| 时间 | $t$ | s | ps 到 ms 甚至更长，视过程而定 |
| 静电势 | $\phi$ | V | 通常约 $10^{-2}$ 到数 V |
| 外加电压 | $V_A$ | V | 通常约 $10^{-1}$ 到数 V；功率器件可远高于此 |
| 电场强度 | $E$ | V/cm | 常见约 $10^3$--$10^5\ \mathrm{V/cm}$ |
| 电子浓度 | $n$ | $\mathrm{cm^{-3}}$ | 可从远低于 $10^{10}$ 到约 $10^{20}$ |
| 空穴浓度 | $p$ | $\mathrm{cm^{-3}}$ | 可从远低于 $10^{10}$ 到约 $10^{20}$ |
| 净掺杂浓度 | $N$ | $\mathrm{cm^{-3}}$ | 常见约 $10^{14}$--$10^{20}$ |
| 受主浓度 | $N_A$ | $\mathrm{cm^{-3}}$ | 常见约 $10^{14}$--$10^{20}$ |
| 施主浓度 | $N_D$ | $\mathrm{cm^{-3}}$ | 常见约 $10^{14}$--$10^{20}$ |
| 总电流密度 | $J$ | $\mathrm{A/cm^2}$ | 强烈依赖器件，常见可跨越多个数量级 |
| 电子电流密度 | $J_n$ | $\mathrm{A/cm^2}$ | 强烈依赖器件和偏压 |
| 空穴电流密度 | $J_p$ | $\mathrm{A/cm^2}$ | 强烈依赖器件和偏压 |
| 产生率 | $G$ | $\mathrm{cm^{-3}s^{-1}}$ | 强烈依赖产生机制 |
| 复合率 | $R$ | $\mathrm{cm^{-3}s^{-1}}$ | 强烈依赖载流子浓度及寿命 |
| 电子迁移率 | $\mu_n$ | $\mathrm{cm^2/(V\cdot s)}$ | 轻掺杂 Si、300 K 时约 $1350$--$1450$ |
| 空穴迁移率 | $\mu_p$ | $\mathrm{cm^2/(V\cdot s)}$ | 轻掺杂 Si、300 K 时约 $450$--$500$ |
| 电子扩散系数 | $D_n$ | $\mathrm{cm^2/s}$ | 300 K 时典型约 $35$--$37$ |
| 空穴扩散系数 | $D_p$ | $\mathrm{cm^2/s}$ | 300 K 时典型约 $11$--$13$ |

---

## 6. 典型掺杂浓度

硅器件中的掺杂浓度跨越很大的范围。

粗略地，可以使用以下数量级建立直觉：

| 区域或用途 | 典型掺杂浓度 |
|---|---|
| 轻掺杂漂移区 | $10^{14}$--$10^{16}\ \mathrm{cm^{-3}}$ |
| 普通体区、阱区 | $10^{15}$--$10^{18}\ \mathrm{cm^{-3}}$ |
| 较重掺杂区域 | $10^{18}$--$10^{19}\ \mathrm{cm^{-3}}$ |
| 源漏、欧姆接触附近 | $10^{19}$--$10^{20}\ \mathrm{cm^{-3}}$ |

因此，在普通教学型 PN 结问题中，一个很常见的代表性数量级是

$$
\boxed{
N_A,\ N_D\sim10^{16}\ \mathrm{cm^{-3}}
}
$$

但这只是方便建立数量级直觉的参考值，并不是统一默认值。

---

## 7. 迁移率

迁移率一般写成

$$
\mu_n=\mu_n(T,N,E,\ldots),
$$

$$
\mu_p=\mu_p(T,N,E,\ldots).
$$

因此迁移率不仅与温度有关，也与掺杂浓度、电场强度、散射机制等有关。

对于 300 K、低场、轻掺杂体硅，可采用典型值

$$
\boxed{
\mu_n\approx1400\ \mathrm{cm^2/(V\cdot s)}
}
$$

$$
\boxed{
\mu_p\approx450\ \mathrm{cm^2/(V\cdot s)}
}
$$

作为数量级参考。

---

## 8. Einstein 关系与扩散系数

在非简并、低场条件下，

$$
\boxed{
D_n=\mu_n V_T
}
$$

$$
\boxed{
D_p=\mu_p V_T
}
$$

更一般地，

$$
\boxed{
D=\mu V_T.
}
$$

由于

$$
V_T=\frac{k_B T}{q},
$$

因此扩散系数也是温度的函数。

采用 300 K 附近的典型迁移率，可得到数量级上

$$
\boxed{
D_n\approx35\text{--}37\ \mathrm{cm^2/s}
}
$$

$$
\boxed{
D_p\approx11\text{--}13\ \mathrm{cm^2/s}.
}
$$

---

## 9. Debye 长度

采用以下 Debye 长度作为长度归一化尺度：

$$
\boxed{
L_D=
\sqrt{
\frac{\varepsilon_s k_B T}
{q^2n_i}
}
}
$$

利用

$$
V_T=\frac{k_B T}{q},
$$

也可写成

$$
\boxed{
L_D=
\sqrt{
\frac{\varepsilon_s V_T}
{qn_i}
}.
}
$$

这里的 $L_D$ 是以本征载流子浓度 $n_i$ 为参考浓度定义的 Debye 长度。

由于

$$
n_i=n_i(T),
$$

因此

$$
\boxed{
L_D=L_D(T).
}
$$

若仅用于建立数量级直觉，在 300 K 且取

$$
n_i\sim10^{10}\ \mathrm{cm^{-3}}
$$

时，有

$$
\boxed{
L_D\sim 40\ \mu\mathrm{m}.
}
$$

这个长度明显大于典型掺杂区域的局部 Debye 长度，原因是这里使用了很低的本征浓度 $n_i$ 作为参考浓度。

---

## 10. 掺杂半导体中的局部 Debye 长度

如果使用某个实际载流子或掺杂浓度 $N$ 作为参考浓度，则可定义数量级上的局部 Debye 长度

$$
L_D(N)
=
\sqrt{
\frac{\varepsilon_s k_B T}
{q^2N}
}.
$$

在 300 K 附近，其典型数量级为：

| $N$ | $L_D(N)$ 的典型数量级 |
|---:|---:|
| $10^{15}\ \mathrm{cm^{-3}}$ | $\sim 130\ \mathrm{nm}$ |
| $10^{16}\ \mathrm{cm^{-3}}$ | $\sim 40\ \mathrm{nm}$ |
| $10^{17}\ \mathrm{cm^{-3}}$ | $\sim 13\ \mathrm{nm}$ |
| $10^{18}\ \mathrm{cm^{-3}}$ | $\sim 4\ \mathrm{nm}$ |

因此必须区分：

1. 以 $n_i$ 定义的归一化 Debye 长度；
2. 以实际载流子浓度或掺杂浓度定义的局部屏蔽长度。

---

## 11. 归一化参考量

定义特征扩散系数

$$
D_0.
$$

$D_0$ 是人为选取的参考扩散系数，不是固定的材料基本常数。

一种方便的选择是使特征迁移率

$$
\mu_0=\frac{D_0}{V_T}
$$

具有便于计算的数量级。

例如可以选择

$$
\boxed{
\mu_0=1000\ \mathrm{cm^2/(V\cdot s)}
}
$$

于是定义

$$
\boxed{
D_0=\mu_0V_T.
}
$$

因为 $T$ 不固定，所以 $D_0$ 仍然通过 $V_T(T)$ 保持温度依赖。

在 300 K 时，这一选择对应大约

$$
D_0\approx25.9\ \mathrm{cm^2/s}.
$$

归一化尺度规定如下：

| 变量 | 参数 | 归一化因子 |
|---|---|---|
| 坐标 | $x$ | $L_D$ |
| 时间 | $t$ | $L_D^2/D_0$ |
| 静电势、外加电压 | $\phi,\ V_A$ | $V_T$ |
| 电场强度 | $E$ | $V_T/L_D$ |
| 载流子浓度、杂质浓度 | $n,p,N,N_A,N_D$ | $n_i$ |
| 电流密度 | $J,J_n,J_p$ | $qD_0n_i/L_D$ |
| 产生率、复合率 | $G,R$ | $D_0n_i/L_D^2$ |
| 迁移率 | $\mu_n,\mu_p$ | $D_0/V_T$ |
| 扩散系数 | $D_n,D_p$ | $D_0$ |

---

## 12. 无量纲变量定义

### 12.1 空间

$$
\boxed{
\tilde{x}=\frac{x}{L_D}
}
$$

### 12.2 时间

$$
\boxed{
\tilde{t}
=
\frac{t}{L_D^2/D_0}
=
\frac{D_0t}{L_D^2}
}
$$

### 12.3 电势和电压

$$
\boxed{
\tilde{\phi}=\frac{\phi}{V_T}
}
$$

$$
\boxed{
\tilde{V}_A=\frac{V_A}{V_T}
}
$$

### 12.4 电场

$$
\boxed{
\tilde{E}
=
\frac{E}{V_T/L_D}
=
\frac{EL_D}{V_T}
}
$$

### 12.5 载流子浓度

$$
\boxed{
\tilde{n}=\frac{n}{n_i},
\qquad
\tilde{p}=\frac{p}{n_i}
}
$$

### 12.6 掺杂浓度

$$
\boxed{
\tilde{N}=\frac{N}{n_i},
\qquad
\tilde{N}_A=\frac{N_A}{n_i},
\qquad
\tilde{N}_D=\frac{N_D}{n_i}
}
$$

### 12.7 电流密度

$$
\boxed{
\tilde{J}
=
\frac{J}{qD_0n_i/L_D}
}
$$

$$
\boxed{
\tilde{J}_n
=
\frac{J_n}{qD_0n_i/L_D},
\qquad
\tilde{J}_p
=
\frac{J_p}{qD_0n_i/L_D}
}
$$

### 12.8 迁移率

$$
\boxed{
\tilde{\mu}_n
=
\frac{\mu_n}{D_0/V_T},
\qquad
\tilde{\mu}_p
=
\frac{\mu_p}{D_0/V_T}
}
$$

### 12.9 扩散系数

$$
\boxed{
\tilde{D}_n=\frac{D_n}{D_0},
\qquad
\tilde{D}_p=\frac{D_p}{D_0}
}
$$

---

## 13. 产生率与复合率

如果 $G$ 和 $R$ 表示标准连续性方程中的载流子数产生率和复合率，则

$$
[G]=[R]=\mathrm{cm^{-3}s^{-1}}.
$$

自然的归一化尺度为

$$
\boxed{
G_0=R_0=
\frac{D_0n_i}{L_D^2}
}
$$

因此

$$
\boxed{
\tilde{G}
=
\frac{G}{D_0n_i/L_D^2}
}
$$

$$
\boxed{
\tilde{R}
=
\frac{R}{D_0n_i/L_D^2}.
}
$$

产生率和复合率没有一个适用于所有器件的单一典型值。

例如，在简单寿命近似下，

$$
R\sim\frac{\Delta n}{\tau},
$$

因此如果

$$
\Delta n\sim10^{15}\ \mathrm{cm^{-3}},
\qquad
\tau\sim1\ \mu\mathrm{s},
$$

则数量级为

$$
R\sim10^{21}\ \mathrm{cm^{-3}s^{-1}}.
$$

这仅用于说明 $G$、$R$ 可以具有非常大的数值范围。

---

## 14. 电荷产生/复合源项

如果连续性方程中考虑的是电荷源项

$$
qG,\qquad qR,
$$

则单位为

$$
\mathrm{C\,cm^{-3}s^{-1}}
=
\mathrm{A/cm^3}.
$$

相应归一化尺度为

$$
\boxed{
\frac{qD_0n_i}{L_D^2}.
}
$$

因此：

- 对 $G,R$ 本身，应使用

$$
\frac{D_0n_i}{L_D^2};
$$

- 对 $qG,qR$，应使用

$$
\frac{qD_0n_i}{L_D^2}.
$$

二者不能混用。

---

## 15. 电流密度的典型数量级

电流密度高度依赖器件结构、偏置条件以及器件用途，因此不存在类似 $q$ 那样唯一的典型值。

器件分析中可能遇到从极低漏电流密度到

$$
10^0\text{--}10^2\ \mathrm{A/cm^2}
$$

甚至更高的工作电流密度。

因此对于

$$
J,\qquad J_n,\qquad J_p
$$

应优先根据具体器件和边界条件确定，而不是预设单一数值。

---

## 16. 常用参考参数汇总

在不固定 $T$ 的前提下，可以使用下面这组参数作为后续建模的参考约定：

| 参数 | 参考值或定义 |
|---|---|
| $q$ | $1.602176634\times10^{-19}\ \mathrm{C}$ |
| $k_B$ | $1.380649\times10^{-23}\ \mathrm{J/K}$ |
| $\varepsilon_0$ | $8.854\times10^{-14}\ \mathrm{F/cm}$ |
| $\varepsilon_{r,\mathrm{Si}}$ | $\approx11.7$ |
| $\varepsilon_{\mathrm{Si}}$ | $\approx1.04\times10^{-12}\ \mathrm{F/cm}$ |
| $T$ | **不指定** |
| $V_T$ | $k_BT/q$ |
| $n_i$ | $n_i(T)$；300 K 时约 $10^{10}\ \mathrm{cm^{-3}}$ |
| $\mu_n$ | 轻掺杂 Si、300 K 时约 $1400\ \mathrm{cm^2/(V\cdot s)}$ |
| $\mu_p$ | 轻掺杂 Si、300 K 时约 $450\ \mathrm{cm^2/(V\cdot s)}$ |
| $D_n$ | $\mu_nV_T$；300 K 时约 $35$--$37\ \mathrm{cm^2/s}$ |
| $D_p$ | $\mu_pV_T$；300 K 时约 $11$--$13\ \mathrm{cm^2/s}$ |
| $N_A,N_D$ | 常见约 $10^{14}$--$10^{20}\ \mathrm{cm^{-3}}$ |
| $E$ | 常见约 $10^3$--$10^5\ \mathrm{V/cm}$ |
| $L_D$ | $\sqrt{\varepsilon_sV_T/(qn_i)}$ |
| $\mu_0$ | 可取 $1000\ \mathrm{cm^2/(V\cdot s)}$ |
| $D_0$ | $\mu_0V_T$ |

---

## 17. 单位换算约定

常用 cm 制和 SI 制之间有

$$
\boxed{
1\ \mathrm{V/cm}
=
10^2\ \mathrm{V/m}
}
$$

$$
\boxed{
1\ \mathrm{A/cm^2}
=
10^4\ \mathrm{A/m^2}
}
$$

$$
\boxed{
1\ \mathrm{cm^2/(V\cdot s)}
=
10^{-4}\ \mathrm{m^2/(V\cdot s)}
}
$$

$$
\boxed{
1\ \mathrm{cm^{-3}}
=
10^6\ \mathrm{m^{-3}}
}
$$

以及

$$
\boxed{
1\ \mu\mathrm{m}=10^{-4}\ \mathrm{cm}
}
$$

$$
\boxed{
1\ \mathrm{nm}=10^{-7}\ \mathrm{cm}.
}
$$

在后续半导体器件推导中，若无特别说明，优先采用上述 cm 制工程单位。

---

## 18. 使用原则

后续建模中遵循以下原则：

1. $T$ 始终作为外部给定参数，不预设固定数值。
2. $V_T$ 由

   $$
   V_T=\frac{k_BT}{q}
   $$

   计算。

3. $n_i$ 使用与给定温度 $T$ 相匹配的材料模型，即

   $$
   n_i=n_i(T).
   $$

4. $\mu_n$、$\mu_p$ 根据所采用的迁移率模型确定；若只进行数量级估算，可参考轻掺杂 Si 在 300 K 下的典型值。

5. 扩散系数原则上由 Einstein 关系确定：

   $$
   D_n=\mu_nV_T,
   \qquad
   D_p=\mu_pV_T.
   $$

6. 归一化长度定义为

   $$
   L_D=
   \sqrt{
   \frac{\varepsilon_sV_T}
   {qn_i}
   }.
   $$

7. 归一化扩散系数取

   $$
   D_0=\mu_0V_T,
   $$

   其中可选

   $$
   \mu_0=1000\ \mathrm{cm^2/(V\cdot s)}.
   $$

8. 所有归一化参考量均由当前 $T$ 下的 $V_T$、$n_i(T)$ 和 $D_0$ 计算，而不是固定为 300 K 的数值。