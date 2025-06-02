本项目是为了解决 C.A.J. Fletcher 的书Computational Techniques for Fluid Dynamics 1第157页的问题5.2，描述如下：

**5.2** 通过适当的无量纲化，扩散或热传导方程可以写成：

$$
\frac{\partial \bar{\theta}}{\partial t} - \frac{\partial^2 \bar{\theta}}{\partial x^2} = 0 。
$$

给定初始条件 $\bar{\theta}(x,0) = \sin(\pi x) + x$，以及边界条件 $\bar{\theta}(0,t) = 0$，$\bar{\theta}(1,t) = 1.0$，在计算区域 $0 \leq x \leq 1.0$，$0 \leq t \leq 0.20$ 内，使用以下方法求解上述方程的数值解：

* (a) Galerkin 方法，
* (b) 子区域法（subdomain），
* (c) 配点法（collocation）。

应使用下述近似解形式：

$$
\theta = \sin(\pi x) + x + \sum_{j=1}^{N} a_j(t) \left( x^j - x^{j+1} \right) 。
$$

将上述加权残量法应用于无量纲热传导方程，可得到一组关于系数 $a_j(t)$ 的常微分方程。这些微分方程通过数值积分在时间上推进。

请编写程序，对 $N = 3, 5, 7$ 的情形，使用上述三种方法分别求解，并比较精度与收敛速度。

该问题的精确解为：

$$
\bar{\theta} = \sin(\pi x) e^{-\pi^2 t} + x。
$$


