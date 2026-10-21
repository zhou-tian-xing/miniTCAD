# W2 → W3 接口确认表（等待 A/B 双签）

**依据**：仓库根目录 `MiniTCAD_1D_PN_8周实施验收规范_cm工程单位版_VSCode.md`，尤其 §3、§4、§5、§6、§7、§10.2、§14.7。此表为修复后的实际 Python 接口清单，**不是双方已签署的冻结证明**。

| 归属 | 唯一实现路径 | 固定输入/输出 |
| --- | --- | --- |
| B：材料/网格/器件 | `minitcad.core.material.Material`、`minitcad.core.mesh.Mesh1D`、`minitcad.core.device.build_device` | cm 物理配置 → `NormalizedDevice`; `Mesh1D` 同时有 cm 和 nd 的 `M/M-1/M-2` 网格 |
| B：状态和接触 | `minitcad.core.state.State/Fluxes/SolveResult`、`minitcad.core.contact.contact_targets` | `State(u,n,p)` 无量纲且不可变；`contact_targets(device,bias_V)` 返回 `(uL,nL,pL,uR,nR,pR)`；`bias_V` 仍为 V |
| A：尺度 | `minitcad.physics.normalization.build_normalization_scales` | 当前 T 的 Si 参数；显式 `reference_diffusivity_cm2_s` 或 `reference_mobility_cm2_per_Vs` 二选一；`L0[cm],D0[cm²/s],VT[V]` 等 |
| A：Poisson | `minitcad.discretization.poisson.poisson_residual_nd` | `xi,u,n,p` 各 M；生产 `cbar_nd` 为内部 M-2；支持 W2 遗留测试 M；BC 两个无量纲数 → Fu[M] |
| A：SG | `minitcad.discretization.sg.edge_flux` | `State,NormalizedDevice → Fluxes(jn_edges,jp_edges)`；两条边数组均为 M-1 的无量纲向右常规电流 |
| A：SRH | `minitcad.physics.recombination.srh_recombination_nd` | `n_nd,p_nd,tau_n_nd,tau_p_nd → r_nd`，中能级、异寿命公式 |
| A：continuity | `minitcad.discretization.continuity.continuity_residual_nd` | `Fluxes[M-1], r[M], w[M-2] → (Fn[M-2],Fp[M-2])` |
| A：唯一 F | `minitcad.discretization.residual.residual` | `residual(State,NormalizedDevice,bias_V) → ndarray[3M]`，节点交错 `(Fu,Fn,Fp)`，六行 Dirichlet |
| 教学独立 oracle | `examples.poisson_si_oracle.poisson_residual_si_mms` | `x[m],phi[V],rho[C/m³],epsilon[F/m],BC[V]`；内部残差 C/m³，边界 V，**不得直接传入无量纲求解器** |
| B：单位适配底层 | `minitcad.io.units.convert_units` | 明示 `cm`/`si`，只做单位换算，未替代 B 将来负责的完整 YAML 适配入口 |

## 代码约束与风险

- Poisson 内部掺杂 `cbar_nd` 由 `DopingProfile.cv_averaged_net_doping_cm3` 对 CV 分段积分生成，**不得**将 `net_doping_nd` 在结位置的单侧采样传入。
- `State.pack()` 是 `[u0,n0,p0,u1,n1,p1,...]`；生产边界行 `Fu0/Fn0/Fp0/FuM-1/FnM-1/FpM-1` 均为无量纲差值。
- B 须复核并签字：`contact_targets` 正向极性；`cbar_nd` M−2 合同；`residual` 输出与未来 writer 字段映射；cm/SI 输入仅入口转换。**本次不存在 B 签字证据。**
- 不提供 `solve_gummel`、`solve_newton`、正式 `RunResult`/CLI/三文件 writer：这些属于 W4–W7，现阶段不能声称通过相应 Gate。

A 签字：**未签**  
B 签字：**未签**
