# simulation

JSM2009 Fig1:

![pic](picture/JSM_Fig1.png)

$dv_i = \sum_{j\neq i}(J_{ij}v_j - J_{ji} v_i)dt + \sqrt{2}\sigma v_i dW_i$

如果令 $J_{ij} = \frac{1}{N-1}$,

$dv_i = (\tilde{v_i} - v_i)dt + \sqrt{2}\sigma v_i dW_i$

这里 $\tilde{v_i} = \frac{1}{N-1}\sum_{j\neq i}v_j$,并且 $\langle dW_i \rangle = 0,\langle dW_i dW_j \rangle = \delta_{ij} dt$,

显然，$\langle \sum_i dv_i \rangle = 0$

对这个方程进行数值模拟。

数值解方程：

$v_i(t+\Delta t) =v_i(t) + (\overline{v(t)} - \frac{N}{N-1} v_i)\Delta t + \sqrt{2}\sigma v_i(t)\eta \sqrt{\Delta t}, \bar{v}=\frac{1}{N-1}\sum_j v_j(t)$

$v_i(t+\Delta t) =v_i(t)*(1 - \frac{N}{N-1} \Delta t + \sqrt{2}\sigma \eta \sqrt{\Delta t}) + \overline{v(t)} \Delta t$

目标是计算 $\langle v_i^2(t) \rangle$, $\langle v_i(t) v_j(t) \rangle$

$var[v_i](t) = \langle v_i^2(t) \rangle -\langle v_i(t) \rangle^2$

$C_{ij}(t) =\frac{\langle v_i(t) v_j(t) \rangle - \langle v_i(t)\rangle \langle v_j(t) \rangle}{\sqrt{var[v_i](t)var[v_j](t)}} $

## J=-1

$v_i(t+\Delta t) =v_i(t)*(1 + \frac{N}{N-1} \Delta t + \sqrt{2}\sigma \eta \sqrt{\Delta t}) - \overline{v(t)} \Delta t$

## SDE simulation

一种随机微分方程数值模拟的另类方法 - Steven Li的文章 - 知乎
<https://zhuanlan.zhihu.com/p/28628912>

<https://frouah.com/finance%20notes/Euler%20and%20Milstein%20Discretization.pdf>

考虑一维 SDE, Ito 形式微分：

$dX_t = \mu(X_t,t) dt + \sigma(X_t,t) dW_t$

Euler 方法：

$X(t+\Delta t)=X(t)+\mu\left(X_{t}, t\right) \Delta t+\sigma\left(X_{t}, t\right) Z_{t} \sqrt{\Delta t}$

Milstein 方法：

$X(t+\Delta t)=X(t)+\mu\left(X_{t}, t\right) \Delta t+\sigma\left(X_{t}, t\right) Z_{t} \sqrt{\Delta t}+\frac{1}{2} \sigma\left(X_{t}, t\right) \frac{\partial}{\partial x} \sigma\left(X_{t}, t\right)\left[Z_{t}^{2} \Delta t-\Delta t\right]$
