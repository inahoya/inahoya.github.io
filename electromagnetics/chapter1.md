---
layout: single
title: 電磁気学 Electromagnetics
mathjax: true
sidebar:
  nav: "main"
---



## Chapter 1 電磁気学ではなぜベクトル解析が重要なのか？

### 1.1 まずベクトルの演算を定義
#### 1.1.1
**ベクトル**とは、向きと大きさを持つ量である。対して、大きさを表す量を**スカラー**という。ベクトルの和について、交換則

$$\boldsymbol{A}+\boldsymbol{B}=\boldsymbol{B}+\boldsymbol{A}$$

と結合則

$$(\boldsymbol{A}+\boldsymbol{B})+\boldsymbol{C}=\boldsymbol{A}+(\boldsymbol{B}+\boldsymbol{C})$$

が成り立つ。また、スカラー $a$ 倍について分配則

$$a(\boldsymbol{A}+\boldsymbol{B})=a\boldsymbol{A}+a\boldsymbol{B}$$

も成り立つ。

内積

$$\boldsymbol{A}\cdot\boldsymbol{B}\equiv AB\cos\theta$$

を導入すると、交換則

$$\boldsymbol{A}\cdot\boldsymbol{B}=\boldsymbol{B}\cdot\boldsymbol{A}$$

と分配則

$$\boldsymbol{A}\cdot(\boldsymbol{B}+\boldsymbol{C})=\boldsymbol{A}\cdot\boldsymbol{B}+\boldsymbol{A}\cdot\boldsymbol{C}$$

が成り立つ。

外積

$$\boldsymbol{A}\times\boldsymbol{B}\equiv AB\sin\theta\boldsymbol{\hat{n}}$$

を導入する。内積はスカラーであったが、外積はベクトルである。分配則

$$\boldsymbol{A}\times(\boldsymbol{B}+\boldsymbol{C})=(\boldsymbol{A}\times\boldsymbol{B})+(\boldsymbol{A}\times\boldsymbol{C})$$

が成り立つが、交換則は成り立たない。

$$\boldsymbol{A}\times\boldsymbol{B}=-\boldsymbol{B}\times\boldsymbol{A}$$

#### 1.1.2
ベクトルの表現として、カルテシアン座標 $x,y,z$ を使うこともできる。

$$\boldsymbol{A}=A_x\boldsymbol{\hat{x}}+A_y\boldsymbol{\hat{y}}+A_z\boldsymbol{\hat{z}}$$

このとき、

$$\boldsymbol{A}\cdot\boldsymbol{B}=A_xB_x+A_yB_y+A_zB_z$$

$$\boldsymbol{A}\times\boldsymbol{B}=\begin{vmatrix}
\boldsymbol{\hat{x}} & \boldsymbol{\hat{y}} &\boldsymbol{\hat{z}} \\
A_x & A_y& A_z\\
B_x & B_y & B_z
\end{vmatrix}$$

#### 1.1.3
スカラー三重積 $\boldsymbol{A}\cdot(\boldsymbol{B}\times\boldsymbol{C})$ は、その大きさが３つのベクトルのなす平行六面体の体積であり、以下の性質を持つ。

$$\boldsymbol{A}\cdot(\boldsymbol{B}\times\boldsymbol{C})=\boldsymbol{B}\cdot(\boldsymbol{C}\times\boldsymbol{A})=\boldsymbol{C}\cdot(\boldsymbol{A}\times\boldsymbol{B})$$

$$\boldsymbol{A}\cdot(\boldsymbol{B}\times\boldsymbol{C})=\begin{vmatrix}
A_x & A_y& A_z\\
B_x & B_y & B_z\\
C_x & C_y & C_z
\end{vmatrix}$$


ベクトル三重積 $\boldsymbol{A}\times(\boldsymbol{B}\times\boldsymbol{C})$ は、以下の式を満たす。

$$\boldsymbol{A}\times(\boldsymbol{B}\times\boldsymbol{C})=\boldsymbol{B}(\boldsymbol{A}\cdot\boldsymbol{C})-\boldsymbol{C}(\boldsymbol{A}\cdot\boldsymbol{B})$$

#### 1.1.4
位置ベクトルを

$$\boldsymbol{r}\equiv x\boldsymbol{\hat{x}}+y\boldsymbol{\hat{y}}+z\boldsymbol{\hat{z}}$$

と定義する。大きさは

$$r=\sqrt{x^2+y^2+z^2}$$

単位ベクトルは

$$\boldsymbol{\hat{r}}=\frac{\boldsymbol{r}}{r}=\frac{x\boldsymbol{\hat{x}}+y\boldsymbol{\hat{y}}+z\boldsymbol{\hat{z}}}{\sqrt{x^2+y^2+z^2}}$$

となる。微小ベクトル要素を

$$d\boldsymbol{l}=dx\boldsymbol{\hat{x}}+dy\boldsymbol{\hat{y}}+dz\boldsymbol{\hat{z}}$$

と定義する。

電磁気学ではしばしば、観測点 $\boldsymbol{r}$ と源点 $\boldsymbol{r'}$ の距離が問題となる。そこで、以下のようなベクトルを定義する。

$$\mathscr{\boldsymbol{r}} \equiv \boldsymbol{r}-\boldsymbol{r'}$$

#### 1.1.5

ベクトルは、座標系の回転に対して以下の規則に従わなければならない。一般に3次元空間では、

$$\overline{A_i}=\sum_{j=1}^3R_{ij}A_j$$

さらに、9つの成分を持つ2階のテンソルでは、

$$\overline{T_{ij}}=\sum_{k=1}^3\sum_{l=1}^3R_{ik}R_{jl}T_{kl}$$

が成立する。

### 1.2 次にベクトルの微分を整理

#### 1.2.1

1変数の微分は、

$$df=\left(\frac{df}{dx}\right) dx$$

とあらわされ、微分係数はグラフの傾きをあらわす。

#### 1.2.2

ベクトルのような3変数の微分は、以下のようにあらわすことができる。

$$dT=\left(\frac{\partial T}{\partial x}\right)dx+\left(\frac{\partial T}{\partial y}\right)dy+\left(\frac{\partial T}{\partial z}\right)dz$$

微小ベクトル要素 $d\boldsymbol{l}$ を用いて、

$$dT=(\nabla T)\cdot(d\boldsymbol{l})$$

$$\nabla T\equiv\frac{\partial T}{\partial x}\boldsymbol {\hat{x}}+\frac{\partial T}{\partial y}\boldsymbol {\hat{y}}+\frac{\partial T}{\partial z}\boldsymbol {\hat{z}}$$

と表現すると、 $\nabla T$ は最も勾配が急な方向を、 $\lvert\nabla T\rvert$ はその傾きをあらわす。

#### 1.2.3

$\nabla$ は、以下のように

$$ \nabla = \boldsymbol {\hat{x}}\frac{\partial}{\partial x}+\boldsymbol {\hat{y}}\frac{\partial}{\partial y}+\boldsymbol {\hat{z}}\frac{\partial}{\partial z}$$

ベクトル演算子として見なすことができる。通常のベクトルと同じように、一般のベクトルとの内積・外積を考えることができ、それぞれダイバージェンス、カールと呼ぶ。

#### 1.2.4

$$ \nabla\cdot\boldsymbol{v} =\frac{\partial v_x}{\partial x}+\frac{\partial v_y}{\partial y}+\frac{\partial v_z}{\partial z}$$

ダイバージェンスは、ベクトル $\boldsymbol{v}$がどの程度発散しているかをあらわす指標である。

#### 1.2.5

$$ \nabla\times\boldsymbol{v}=\begin{vmatrix}
\boldsymbol {\hat{x}} & \boldsymbol {\hat{y}}& \boldsymbol {\hat{z}}\\
\frac{\partial}{\partial x} & \frac{\partial}{\partial y} & \frac{\partial}{\partial z}\\
v_x & v_y & v_z
\end{vmatrix}$$

カールは、ベクトル $\boldsymbol{v}$がある点を中心にどの程度回転しているかを表す指標である。

#### 1.2.6

非自明な微分公式を6本列挙する。

$$\nabla(fg)=f\nabla g+g\nabla f$$

$$\nabla(A\cdot B)=A\times(\nabla\times B)+B\times(\nabla\times A)+(A\cdot\nabla)B+(B\cdot\nabla)A$$

$$\nabla\cdot(fA)=f(\nabla\cdot A)+A\cdot(\nabla f)$$

$$\nabla\cdot(A\times B)=B\cdot(\nabla\times A)-A\cdot(\nabla\times B)$$

$$\nabla\times(fA)=f(\nabla\times A)-A\times(\nabla f)$$

$$\nabla\times(A\times B)=(B\cdot\nabla)A-(A\cdot\nabla)B+A(\nabla\cdot B)-B(\nabla\cdot A)$$

#### 1.2.7

ベクトルの二階微分も同じように考えることができる。全組み合わせは5通りある。

$$\nabla\cdot(\nabla T)=\frac{\partial^2 T}{\partial x^2}+\frac{\partial^2 T}{\partial y^2}+\frac{\partial^2 T}{\partial z^2}=\nabla^2 T$$

これをラプラシアンと呼ぶ。

$$\nabla\times(\nabla T)=\boldsymbol{0}$$

グラディアントのカールは常に０である。

$$\nabla(\nabla\cdot\boldsymbol{v})$$

はたまに見かけるものの、特別な名前はついていない。

$$\nabla\cdot(\nabla\times\boldsymbol{v})=0$$

カールのダイバージェンスは常に0である。

$$\nabla\times(\nabla\times\boldsymbol{v})=\nabla(\nabla\cdot\boldsymbol{v})-\nabla^2\boldsymbol{v}$$

ダイバージェンスのダイバージェンスは新しい項は出てこない。

### 1.3 ベクトルの積分は重要な定理を与える

#### 1.3.1

線積分

$$\int_a^b \boldsymbol{v}\cdot d\boldsymbol{l}$$

面積分

$$\int_S \boldsymbol{v}\cdot d\boldsymbol{a}$$

ただし、 $d\boldsymbol{a}$ の向きは面の外側を正に取るのが通例である。

体積分

$$\int_V Td\tau$$

ただし、 $d\tau=dxdydz$ と置く。

#### 1.3.2

微積分学の基本定理

$$\int_a^b \left(\frac{df}{dx}\right)dx=f(b)-f(a)$$

これは、線の傾きの総和は端点の差でのみ決まるという定理である。

#### 1.3.3

$$\int_a^b \left(\nabla T\right)\cdot d\boldsymbol{l}=T(b)-T(a)$$

微積分学の基本定理の線積分版である。
系1 $\int_a^b \left(\nabla T\right)\cdot d\boldsymbol{l}$ は経路によらない。
系2 $\oint \left(\nabla T\right)\cdot d\boldsymbol{l}=0$

#### 1.3.4

$$\int_V \left(\nabla \cdot\boldsymbol{v}\right) d\tau=\oint_S \boldsymbol{v}\cdot d\boldsymbol{a}$$

微積分学の基本定理の体積分版である。ガウスの発散定理などと呼ばれる。

#### 1.3.5

$$\int_S \left(\nabla \times\boldsymbol{v}\right) d\boldsymbol{a}=\oint_P \boldsymbol{v}\cdot d\boldsymbol{l}$$

微積分学の基本定理の面積分版である。ストークスの定理と呼ばれる。
系1 $\int \left(\nabla \times\boldsymbol{v}\right) d\boldsymbol{a}$ は面の取り方によらず、その境界線のみによって決まる。
系2 $\oint \left(\nabla \times\boldsymbol{v}\right) d\boldsymbol{a}=0$ ただし閉じた面についての積分をあらわす。

#### 1.3.6

部分積分

$$\int_a^b f\left(\frac{dg}{dx}\right)dx=-\int_a^b g\left(\frac{df}{dx}\right)dx+fg\rvert_a^b$$

をベクトルに拡張することができる。適切な微分公式を積分することによって、例えば

$$\nabla\cdot(fA)=f(\nabla\cdot A)+A\cdot(\nabla f)$$

を体積分して

$$\int_V f(\nabla\cdot A)d\tau = -\int_V A(\nabla f)d\tau +\oint_S fA\cdot d\boldsymbol{a}$$

類似の式は複数導かれる。

### 1.4 球面・円柱座標系におけるベクトルの取り扱い

#### 1.4.1

球面座標系において、

$$x=r\sin\theta\cos\phi,y=r\sin\theta\sin\phi, z=r\cos\theta$$

と変換される。ベクトル $A$ は、

$$ A=A_r\boldsymbol{\hat{r}}+A_\theta\boldsymbol{\hat{\theta}}+A_\phi\boldsymbol{\hat{\phi}}$$

ただし、

$$ \boldsymbol{\hat{r}}=\sin\theta\cos\phi\boldsymbol{\hat{x}}+\sin\theta\sin\phi\boldsymbol{\hat{y}}+\cos\theta\boldsymbol{\hat{z}}$$

$$ \boldsymbol{\hat{\theta}}=\cos\theta\cos\phi\boldsymbol{\hat{x}}+\cos\theta\sin\phi\boldsymbol{\hat{y}}-\sin\theta\boldsymbol{\hat{z}}$$

$$ \boldsymbol{\hat{\phi}}=-\sin\phi\boldsymbol{\hat{x}}+cos\phi\boldsymbol{\hat{y}}$$

である。

線積分の作用素は

$$ d\boldsymbol{l}=dr\boldsymbol{\hat{r}}+rd\theta\boldsymbol{\hat{\theta}}+r\sin\theta d\phi\boldsymbol{\hat{\phi}}$$

体積分の作用素は

$$d\tau=r^2\sin\theta drd\theta d\phi$$

面積分の作用素は位置によって異なる。

ベクトルの微分についても、係数が変化し、

$$\nabla T=\frac{\partial T}{\partial r}\boldsymbol{\hat{r}}+\frac{1}{r}\frac{\partial T}{\partial\theta}\boldsymbol{\hat{\theta}}+\frac{1}{r\sin\theta}\frac{\partial T}{\partial\phi}\boldsymbol{\hat{\phi}}$$

$$\nabla\cdot\boldsymbol{v}=\frac{1}{r^2}\frac{\partial}{\partial r}(r^2v_r)+\frac{1}{r\sin\theta}\frac{\partial}{\partial\theta}(\sin\theta v_\theta)+\frac{1}{r\sin\theta}\frac{\partial v_\phi}{\partial\phi}$$

$$\nabla\times\boldsymbol{v}=\frac{1}{r\sin\theta}\left[\frac{\partial}{\partial \theta}(\sin\theta v_\phi)-\frac{\partial v_\theta}{\partial\phi}\right]\boldsymbol{\hat{r}}+\frac{1}{r}\left[\frac{1}{\sin\theta}\frac{\partial v_r}{\partial \phi}-\frac{\partial}{\partial r}(rv_\phi)\right]\boldsymbol{\hat{\theta}}$$

$$+\frac{1}{r}\left[\frac{\partial}{\partial r}(r v_\theta)-\frac{\partial v_r}{\partial\theta}\right]\boldsymbol{\hat{\phi}}$$

$$\nabla^2T=\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial T}{\partial r}\right)+\frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta \frac{\partial T}{\partial \theta}\right)+\frac{1}{r^2\sin^2\theta}\frac{\partial^2 T}{\partial\phi^2}$$

#### 1.4.2

円柱座標系も同様に、

$$x=s\cos\phi,y=s\sin\phi, z=z$$

と変換される。ベクトル $A$ は、

$$ A=A_s\boldsymbol{\hat{s}}+A_\phi\boldsymbol{\hat{\phi}}+A_z\boldsymbol{\hat{z}}$$

ただし、

$$ \boldsymbol{\hat{s}}=\cos\phi\boldsymbol{\hat{x}}+\sin\phi\boldsymbol{\hat{y}}$$

$$ \boldsymbol{\hat{\phi}}=-\sin\phi\boldsymbol{\hat{x}}+\cos\phi\boldsymbol{\hat{y}}$$

$$ \boldsymbol{\hat{z}}=\boldsymbol{\hat{z}}$$

である。

線積分の作用素は

$$ d\boldsymbol{l}=ds\boldsymbol{\hat{s}}+sd\phi\boldsymbol{\hat{\phi}}+dz\boldsymbol{\hat{z}}$$

体積分の作用素は

$$d\tau=sdsd\phi dz$$

面積分の作用素は位置によって異なる。

ベクトルの微分についても、係数が変化し、

$$\nabla T=\frac{\partial T}{\partial s}\boldsymbol{\hat{s}}+\frac{1}{s}\frac{\partial T}{\partial\phi}\boldsymbol{\hat{\phi}}+\frac{\partial T}{\partial z}\boldsymbol{\hat{z}}$$

$$\nabla\cdot\boldsymbol{v}=\frac{1}{s}\frac{\partial}{\partial s}(sv_s)+\frac{1}{s}\frac{\partial v_\phi}{\partial\phi}+\frac{\partial v_z}{\partial z}$$

$$\nabla\times\boldsymbol{v}=\left(\frac{1}{s}\frac{\partial v_z}{\partial \phi}-\frac{\partial v_\phi}{\partial z}\right)\boldsymbol{\hat{s}}+\left(\frac{\partial v_s}{\partial z}-\frac{\partial v_z}{\partial s}\right)\boldsymbol{\hat{\phi}}+\frac{1}{s}\left[\frac{\partial}{\partial s}(s v_\phi)-\frac{\partial v_s}{\partial\phi}\right]\boldsymbol{\hat{z}}$$

$$\nabla^2T=\frac{1}{s}\frac{\partial}{\partial s}\left(s\frac{\partial T}{\partial s}\right)+\frac{1}{s^2}\frac{\partial^2 T}{\partial\phi^2}+\frac{\partial^2 T}{\partial z^2}$$
