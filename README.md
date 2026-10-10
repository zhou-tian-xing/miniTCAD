# miniTCAD
课程作业


## ai_example：一维 PN 结漂移扩散模拟器

单文件示例，解一维稳态漂移扩散方程组（泊松方程 + 电子/空穴连续性方程），
用 Scharfetter-Gummel 离散、SRH 复合、阻尼 Newton 求解，再扫描偏压得到 I-V 曲线。

    python ai_example/pn_junction.py        # 求解 + 自检 + 偏压扫描
    python ai_example/test_pn_junction.py   # 16 个测试（不需要 pytest）

只依赖 numpy 和 scipy。测试对的是解析结果而非程序自己的输出：内建电势、
质量作用定律 np = ni^2、耗尽近似峰值电场、短二极管饱和电流、电流守恒、理想因子。

* `pn_junction.py`      —— 全部代码，按 13 节编号，注释多过代码
* `test_pn_junction.py` —— 测试，也可以直接用 pytest 收集
