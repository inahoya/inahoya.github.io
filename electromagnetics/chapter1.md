---
title: 電磁気学 Electromagnetics
layout: single
sidebar:
  nav: "main"
  
---

## Chapter 1 電磁気学ではなぜベクトル解析が重要なの？

### 1.1 ベクトルの計算
#### 1.1.1
**ベクトル**とは、向きと大きさを持つ量である。対して、大きさを表す量を**スカラー**という。ベクトルの和について、交換則
$$\boldsymbol{A}+\boldsymbol{B}=\boldsymbol{B}+\boldsymbol{A}$$
と結合則
$$(\boldsymbol{A}+\boldsymbol{B})+\boldsymbol{C}=\boldsymbol{A}+(\boldsymbol{B}+\boldsymbol{C})$$
が成り立つ。また、スカラー$a$倍について分配則
$$a(\boldsymbol{A}+\boldsymbol{B})=a\boldsymbol{A}+a\boldsymbol{B}$$も成り立つ。
内積
$$\boldsymbol{A}\cdot\boldsymbol{B}\equiv AB\cos\theta$$を導入すると、交換則
$$\boldsymbol{A}\cdot\boldsymbol{B}=\boldsymbol{B}\cdot\boldsymbol{A}$$と分配則
$$\boldsymbol{A}\cdot(\boldsymbol{B}+\boldsymbol{C})=\boldsymbol{A}\cdot\boldsymbol{B}+\boldsymbol{A}\cdot\boldsymbol{C}$$が成り立つ。
外積
$$\boldsymbol{A}\times\boldsymbol{B}\equiv AB\sin\theta\boldsymbol{\hat{n}}$$を導入する。内積はスカラーであったが、外積はベクトルである。計算する上で便利なのもそうだが、そもそも$\sin\theta$の符号が不定なので、向きを決めなければ定義できないという都合もある。分配則
$$\boldsymbol{A}\times(\boldsymbol{B}+\boldsymbol{C})=(\boldsymbol{A}\times\boldsymbol{B})+(\boldsymbol{A}\times\boldsymbol{C})$$が成り立つが、交換則は成り立たない。
$$\boldsymbol{A}\times\boldsymbol{B}=-\boldsymbol{B}\times\boldsymbol{A}$$
#### 1.1.2
ベクトルの表現として、カルテシアン座標$$x,y,z$$を使うこともできる。
$$\boldsymbol{A}=A_x\boldsymbol{\hat{x}}+A_y\boldsymbol{\hat{y}}+A_z\boldsymbol{\hat{z}}$$このとき、
$$\boldsymbol{A}\cdot\boldsymbol{B}=A_xB_x+A_yB_y+A_zB_z$$

$$\boldsymbol{A}\times\boldsymbol{B}=\begin{vmatrix}
\boldsymbol{\hat{x}} & \boldsymbol{\hat{y}} &\boldsymbol{\hat{z}} \\
A_x & A_y& A_z\\
B_x & B_y & B_z
\end{vmatrix}$$

#### 1.1.3
スカラー三重積$$\boldsymbol{A}\cdot(\boldsymbol{B}\times\boldsymbol{C})$$は、その大きさが３つのベクトルのなす平行六面体の体積であり、以下の性質を持つ。
$$\boldsymbol{A}\cdot(\boldsymbol{B}\times\boldsymbol{C})=\boldsymbol{B}\cdot(\boldsymbol{C}\times\boldsymbol{A})=\boldsymbol{C}\cdot(\boldsymbol{A}\times\boldsymbol{B})$$

$$\boldsymbol{A}\cdot(\boldsymbol{B}\times\boldsymbol{C})=\begin{vmatrix}
A_x & A_y& A_z\\
B_x & B_y & B_z\\
C_x & C_y & C_z
\end{vmatrix}$$


ベクトル三重積$$\boldsymbol{A}\times(\boldsymbol{B}\times\boldsymbol{C})$$は、以下の式を満たす。
$$\boldsymbol{A}\times(\boldsymbol{B}\times\boldsymbol{C})=\boldsymbol{B}(\boldsymbol{A}\cdot\boldsymbol{C})-\boldsymbol{C}(\boldsymbol{A}\cdot\boldsymbol{B})$$

#### 1.1.4
位置ベクトル$$\boldsymbol{r}\equiv x\boldsymbol{\hat{x}}+ y \boldsymbol{\hat{y}}+ z\boldsymbol{\hat{z}}$$
