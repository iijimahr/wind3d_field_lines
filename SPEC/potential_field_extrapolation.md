# 直交曲線座標系におけるポテンシャル磁場計算

2026-03-06: initial commit

## ターゲット問題

磁場はポテンシャルで評価(**磁場をポテンシャルから計算する限りrotB=0**):
$$
{\bf B}=-\nabla\Psi
$$
$\nabla\cdot{\bf B}=0$ 条件より:
$$
\nabla\cdot\left(\nabla\Psi\right)=0
$$
下部境界条件は表面磁場の動径成分 $B_0$ :
$$
B_n = B_0 \quad \mathrm{@}\ \xi_3=0
$$
上部境界条件は水平磁場ゼロ:
$$
\mathbf{B}_h = 0 \quad \mathrm{@}\ \xi_3 = L_3
$$

## 直交座標系の演算子

デカルト座標:
$$
\xi_j=(x,y,z),\ h_j=(1,1,1)
$$

局所球座標 (Iijima et al. 2023):
$$
\xi_j=(\theta,\phi,r)
h_j=(r\sqrt{f},r\sqrt{f},1)
$$

スカラーの勾配:
$$
\left(\nabla\Phi\right)_j
=\frac{1}{h_j}\frac{{\partial}{\Phi}}{{\partial}{\xi_j}}
$$

ベクトルの発散:
$$
\nabla\cdot{\bf A}
=J \sum_j \frac{\partial}{\partial \xi_j}
\left(\frac{A_j}{Jh_j}\right)
$$

ヤコビアン:
$$
J = \frac{1}{h_1 h_2 h_3}
$$

## 波数空間の定義と演算

水平波数と水平位置の積:
$$
k_j\xi_j=k_xx+k_yy=k_\xi\xi+k_\eta\eta
$$
FFTの定義:
$$
F(k_\xi,k_\eta)
=\mathcal{F}\left\{f(\xi,\eta)\right\}
=\iint{f(\xi,\eta)}\mathrm{e}^{-i(k_\xi\xi+k_\eta\eta)}d{\eta}d{\xi}
$$
FFTの逆変換:
$$
\mathcal{F^{-1}}\left\{F(k_\xi, k_\eta)\right\}=f(\xi, \eta)
$$
FFTを利用した変数の微分:
$$
\frac{{\partial}f(k_\xi,k_\eta)}{{\partial}{\xi}}
=\mathcal{F^{-1}}\left\{ik_\xi{\cdot}F(k_\xi,k_\eta)\right\}
$$

## 直交曲線座標系におけるポテンシャル磁場の定式化

- スケールファクタ $h_j$ ($j=1,2,3$) は全て $\xi_3$ のみの関数であると仮定する。
  - 局所球座標系はこの条件を満たす。
- 水平方向 $\xi_1$ と $\xi_2$ は周期的であると仮定する。
  - 水平FFTを利用して、水平微分を計算するため。
- 以下の「$f$ に対するポアソン方程式」を各波数 $(k_1, k_2)$ ごとに数値的に解けば、ポテンシャル磁場を計算可能。
  - 単なる2階線形常微分方程式の境界値問題であり、TDMAなどで素直に計算可能。
- 注: 上部境界で $f = 0$ を課すのはポテンシャルの基準点 (ゼロ波数の場合)、水平磁場ゼロ条件 (非ゼロ波数の場合) のため。

ポアソン方程式:

$$
\nabla^2\Psi = J \sum_j \frac{\partial}{\partial \xi_j}
\left( \frac{1}{Jh_j^2} \frac{\partial \Psi}{\partial \xi_j} \right) = 0
$$

境界磁場のフーリエ変換:

$$
A(k_1,k_2) = \mathcal{F}\{B_3(\xi_1,\xi_2, \xi_3=0)\}
$$

天下り式による解:

$$
\Psi(\xi_1,\xi_2,\xi_3) = \mathcal{F}^{-1}\{A(k_1,k_2) f(k_1,k_2,\xi_3)\}
$$

$$
\frac{\partial \Psi}{\partial \xi_1} = - h_1 B_1 = \mathcal{F}^{-1}\{ik_1 A f\}
$$

$$
\frac{\partial \Psi}{\partial \xi_3} = - h_3 B_3 = \mathcal{F}^{-1}\{A \frac{\partial f}{\partial \xi_3}\}
$$

$f$ に対するポアソン方程式:

$$
J \frac{\partial}{\partial \xi_3}
\left( \frac{1}{Jh_3^2} \frac{\partial f}{\partial \xi_3} \right) -
\left(\frac{k_1^2}{h_1^2} + \frac{k_2^2}{h_2^2} \right) f = 0
$$

$f$ の境界条件:

$$
\frac{\partial f}{\partial \xi_3} = - h_3 \quad \mathrm{@}\ \xi_3=0
$$

$$
f = 0 \quad \mathrm{@}\ \xi_3=L_3
$$

## 直交曲線座標系におけるポテンシャル磁場の離散化

- 方針: 行列形式に整理し、TDMAで解く。
- 境界条件では、ゴースト点を想定して外部に仮想の物理点を置き、内部点と同様の離散化を行う。
  - ただし、座標やスケールファクタはゴースト点においても定義されているものとする。
- 注意: 内部点は最小が $\xi_3 = 0$ 、最大が $\xi_3 = L_3 - \Delta \xi_3/2$ である。
  - セル内部点からポテンシャル磁場を計算する際の利便性のため。

座標系の取り方:

$$
\xi_{3, k} = k \Delta \xi_3 \quad (k=0,1,\cdots,N_3-1)
$$

格子刻み幅:

$$
\Delta \xi_3 = \frac{L_3}{N_3 - 1/2}
$$

ポアソン方程式の内部点:

$$
J_k \frac{1}{\Delta \xi_3}
\left(
  \frac{1}{J_{k+1/2} h_{3,k+1/2}^2} \frac{f_{k+1}-f_k}{\Delta \xi_3} -
  \frac{1}{J_{k-1/2} h_{3,k-1/2}^2} \frac{f_k-f_{k-1}}{\Delta \xi_3}
\right)
- \left(\frac{k_1^2}{h_{1,k}^2} + \frac{k_2^2}{h_{2,k}^2} \right) f_k = 0
$$

ポアソン方程式の下部境界:

$$
f_{k=-1} = f_{k=1} + 2 h_{3,k=0} \Delta \xi_3
$$

ポアソン方程式の上部境界:

$$
f_{k=N_3} = - f_{k=N_3-1}
$$

## 実装の検証

- デカルト座標系で解析解と整合する妥当な解が得られるか。
- 水平方向に積分した鉛直磁束量が、各高さで保存されるか、特に下部境界と一致するか。
- 電流 $\nabla\times\mathbf{B}$ の強度が内部点で十分小さいか。
- 磁場の発散 $\nabla\cdot\mathbf{B} = 0$ が内部点で十分小さいか。
  - ただし、上部境界で水平磁場ゼロ条件を課しているため、上部境界付近では誤差が発生することに注意。
