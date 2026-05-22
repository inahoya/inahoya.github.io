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

微小ベクトル $d\boldsymbol{l}$ を用いて、

$$dT=(\nabla T)\cdot(d\boldsymbol{l})$$

$$\nabla T\equiv\frac{\partial T}{\partial x}\boldsymbol {\hat{x}}+\frac{\partial T}{\partial y}\boldsymbol {\hat{y}}+\frac{\partial T}{\partial z}\boldsymbol {\hat{z}}$$

と表現すると、 $\nabla T$ は最も勾配が急な方向を、 $|\nabla T|$ はその傾きをあらわす。

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

