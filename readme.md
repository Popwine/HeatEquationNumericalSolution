# 项目：求解一维热传导方程的加权残量法

本项目是为了解决 C.A.J. Fletcher 的书 *Computational Techniques for Fluid Dynamics 1* 第157页的问题5.2。该问题要求使用三种不同的加权残量法（伽辽金法、子区域法、配点法）求解一个无量纲化的一维热传导方程。

## 问题描述

通过适当的无量纲化，扩散或热传导方程可以写成：

$$
\frac{\partial \bar{\theta}}{\partial t} - \frac{\partial^2 \bar{\theta}}{\partial x^2} = 0 。
$$

给定初始条件 $\bar{\theta}(x,0) = \sin(\pi x) + x$，以及边界条件 $\bar{\theta}(0,t) = 0$，$\bar{\theta}(1,t) = 1.0$，在计算区域 $0 \leq x \leq 1.0$，$0 \leq t \leq 0.20$ 内求解。

应使用下述近似解形式：

$$
\theta_{近似} = \sin(\pi x) + x + \sum_{j=1}^{N} a_j(t) \left( x^j - x^{j+1} \right) 。
$$

目标是针对 $N = 1, 3, 5, 7$ 的情形，使用伽辽金法、子区域法和配点法分别求解，并比较其精度与收敛特性。

该问题的精确解为：

$$
\bar{\theta}_{精确} = \sin(\pi x) e^{-\pi^2 t} + x。
$$

## 项目结构

```
C:.
├───.vscode                 # VS Code 编辑器配置
├───build                   # CMake 构建输出目录 (被 .gitignore 忽略)
├───HENS_plot               # Python 绘图脚本及虚拟环境
│   ├───.venv               # Python 虚拟环境 (被 .gitignore 忽略)
│   ├───plot.py             # Python 脚本，用于读取 C++ 输出并绘图、计算误差
│   ├───pyproject.toml      # Python 项目配置文件 (如使用 Poetry, Hatch, etc.)
│   ├───README.md           # Python 绘图部分的说明
├───include                 # C++ 头文件目录
│   ├───HENS_math_core.h    # 求解器核心逻辑和数学工具的头文件
│   └───HENS_matrix.h       # 矩阵类定义头文件
├───src                     # C++ 源文件目录
│   ├───HENS_math_core.cpp  # 求解器核心逻辑和数学工具的实现
│   ├───HENS_matrix.cpp     # 矩阵类的实现
│   └───HENS.cpp            # 主程序 (main 函数)
├───.gitignore              # Git 忽略文件配置
├───CMakeLists.txt          # CMake 构建系统配置文件
├───output.txt              # C++ 程序输出的 a_j(t) 系数值 (被 .gitignore 忽略)
├───readme.md               # 本项目的主要说明文件
└───run.ps1                 # PowerShell 运行脚本
```

## 编译与运行

### C++ 部分

1.  确保已安装 CMake 和 C++ 编译器 (如 GCC, Clang, MSVC)。
2.  创建构建目录并进入：
    ```bash
    mkdir build
    cd build
    ```
3.  运行 CMake 配置项目：
    ```bash
    cmake ..
    ```
4.  编译项目：
    ```bash
    cmake --build .
    # 或者直接使用 make (Linux/macOS) 或 nmake/msbuild (Windows)
    ```
5.  运行生成的可执行文件 (通常在 `build` 目录下，名称可能为 `HENS` 或 `HENS.exe`)：
    ```bash
    ./HENS
    # 或者在 Windows 上
    .\HENS.exe
    ```
    这将生成 `output.txt` 文件在项目根目录。

### Python 绘图部分

1.  导航到 `HENS_plot` 目录：
    ```bash
    cd HENS_plot
    ```
2.  (推荐) 创建并激活 Python 虚拟环境：
    ```bash
    python -m venv .venv
    source .venv/bin/activate  # Linux/macOS
    # .venv\Scripts\activate   # Windows
    ```
3.  安装必要的 Python 包：
    ```bash
    pip install numpy matplotlib tabulate
    ```
4.  运行绘图脚本 `plot.py` (确保 `output.txt` 在其预期的相对路径，通常是 `../output.txt`)：
    ```bash
    python plot.py
    ```
    这将显示比较图，并将图像保存为 `solution_comparison_rmse.png`，同时在控制台输出 RMSE 误差表格。

## 结果

### 数值解与精确解比较图

下图展示了在 $t=0.20$ 时，不同方法和不同 $N$ 值下的数值解（红色虚线）与精确解（黑色实线）的比较，并在每个子图中显示了均方根误差 (RMSE)。

![数值解与精确解比较](HENS_plot/solution_comparison_rmse.png)

### 均方根误差 (RMSE) 总结

下表总结了在 $t=0.20$ 时，不同方法和不同 $N$ 值下的均方根误差：

| Method      |     N=1 |       N=3 |       N=5 |       N=7 |
| :---------- | ------: | --------: | --------: | --------: |
| Galerkin    | 2.366e-02 | 4.886e-04 | 5.376e-06 | 3.686e-08 |
| Subdomain   | 8.948e-02 | 3.166e-03 | 4.791e-05 | 5.268e-07 |
| Collocation | 1.125e-01 | 4.307e-03 | 1.010e-04 | 1.505e-06 |

### 结果分析

从上图和 RMSE 表格可以看出：

1.  **收敛性：** 对于所有三种方法，随着 $N$ 值的增加（即近似解中使用的基函数数量增加），数值解的精度都显著提高，RMSE 值减小。这表明这些方法对于此问题是收敛的。
2.  **精度比较：**
    *   **伽辽金法 (Galerkin)** 在相同的 $N$ 值下通常能提供最高的精度（最小的 RMSE）。它的收敛速度也最快，尤其是在 $N$ 较大时，误差迅速减小到非常低的水平。
    *   **子区域法 (Subdomain)** 的精度介于伽辽金法和配点法之间。
    *   **配点法 (Collocation)** 在相同的 $N$ 值下精度相对较低，但其实现最为简单，因为它避免了积分运算。
3.  **图像吻合度：** 当 $N$ 增加到 5 或 7 时，对于伽辽金法，数值解曲线与精确解曲线在视觉上几乎完全重合。其他方法在 $N=7$ 时也表现出很好的吻合度。

## 实现原理概述

本项目通过以下步骤实现数值求解：

1.  **近似解代入：** 将给定的近似解形式 $\theta_{近似}$ 代入原始的偏微分方程 $\frac{\partial \bar{\theta}}{\partial t} - \frac{\partial^2 \bar{\theta}}{\partial x^2} = 0$，得到一个包含未知系数 $a_j(t)$ 及其时间导数 $\dot{a}_j(t)$ 的**残差方程 $R(x,t) = 0$**。

2.  **加权残量法：** 为了确定 $a_j(t)$，我们应用加权残量法。其核心思想是选择一组权函数 $W_i(x)$，并要求残差在加权积分意义下为零：
    $$ \int_{0}^{1} R(x,t) W_i(x) \, dx = 0 \quad \text{for } i = 1, 2, \dots, N $$
    这会产生 $N$ 个关于 $a_j(t)$ 和 $\dot{a}_j(t)$ 的方程。

3.  **不同方法的权函数选择：**
    *   **伽辽金法：** 权函数 $W_i(x)$ 选择为近似解中的第 $i$ 个基函数，即 $W_i(x) = x^i - x^{i+1}$。积分在整个定义域 $[0,1]$ 上进行。本项目中，这些积分使用辛普森数值积分法计算。
    *   **子区域法：** 将求解域 $[0,1]$ 等分成 $N$ 个子区域。第 $i$ 个权函数 $W_i(x)$ 在第 $i$ 个子区域内为 1，其余为 0。因此，积分 $\int R(x,t) dx = 0$ 在每个子区域上分别执行。本项目中，这些积分是解析计算的。
    *   **配点法：** 在求解域内部选择 $N$ 个离散的配点 $x_i$ (本项目选择等间距点 $x_i = i/(N+1)$ for $i=1 \dots N$)。权函数可以看作是 $W_i(x) = \delta(x-x_i)$，这使得加权残量方程简化为 $R(x_i, t) = 0$，即要求残差在每个配点处精确为零，无需积分。

4.  **常微分方程组 (ODEs)：**
    上述加权残量过程将原始的偏微分方程转化为一个关于 $a_j(t)$ 的一阶常微分方程组，形式通常为：
    $$ \mathbf{M} \frac{d\mathbf{a}}{dt} + \mathbf{K} \mathbf{a} = \mathbf{F} $$
    其中 $\mathbf{a}(t) = [a_1(t), \dots, a_N(t)]^T$，$\mathbf{M}$ 是质量矩阵，$\mathbf{K}$ 是刚度矩阵，$\mathbf{F}$ 是力向量或常数项向量。它们的元素由第3步中的积分或求值过程确定。

5.  **时间积分：**
    将上述 ODE 组改写为标准形式 $\frac{d\mathbf{a}}{dt} = \mathbf{M}^{-1}(\mathbf{F} - \mathbf{K}\mathbf{a})$。
    初始条件由 $\bar{\theta}(x,0) = \sin(\pi x) + x$ 确定，这导致所有 $a_j(0) = 0$。
    然后，使用**四阶龙格-库塔法 (RK4)** 对这个 ODE 系统进行时间积分，从 $t=0$ 到 $t=0.20$，从而求解出 $a_j(t)$ 在每个时间步的值。本项目输出 $a_j(0.20)$。

6.  **精度评估：**
    得到 $a_j(0.20)$ 后，可以构建 $t=0.20$ 时的数值解 $\theta_{近似}(x, 0.20)$，并将其与精确解 $\bar{\theta}_{精确}(x, 0.20)$ 进行比较，通过计算均方根误差 (RMSE) 来量化精度。

## 未来工作与展望

*   实现自适应步长的龙格-库塔方法以提高时间积分效率和精度控制。
*   对配点法的配点选择策略进行研究（如使用高斯点）。
*   扩展到二维问题。
*   将 C++ 计算核心与 Python 绘图更紧密地集成（例如通过 pybind11）。

