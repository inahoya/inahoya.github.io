---
layout: single
title: 電磁気学 Electromagnetics
mathjax: true
toc: true
sidebar:
  nav: "main"
---

## Chapter 2 静止した電荷の理論

### 2.1 電場とは何か

#### 2.1.1

電荷 $q_1, q_2, q_3,\cdots$ が存在する空間に、電荷 $Q$ を置いたとする。 $Q$ に働く力 $F$ は、重ね合わせの原理により、 $F=F_1+F_2+F_3+\cdots$ とあらわすことができる。当たり前に思えるかもしれないが、これは自明なことではない。また、電荷 $q_1, q_2, q_3,\cdots$ を一つの電荷 $q$ で置き換えることもできる。厳密な議論はのちに行うとして、今はすべての電荷が静止している特別な場合を考える。

#### 2.1.2

実験的事実として、クーロンの法則がある。

$$F=\frac{1}{4\pi\varepsilon_0}\frac{qQ}{\mathscr{r}^2}\boldsymbol{\hat{\mathscr{r}}}$$

$$<span style="font-family: kaufmann; font-size: 20px;">normal</span>$$

ここで、 $\varepsilon_0$ は真空の誘電率

$$\varepsilon_0=8.85\times10^{-12}\frac{C^2}{N\cdot m^2}$$

である。クーロンの法則は重ね合わせの原理を保証する。

#### 2.1.3

クーロンの法則を書き換えると、

$$F=QE$$

$$E(\boldsymbol{r})\equiv\frac{1}{4\pi\varepsilon_0}\sum_{i=1}^n\frac{q_i}{\mathscr{r}_i^2}\boldsymbol{\hat{\mathscr{r}}}_i$$

この $E$ を電場と呼ぶ。

#### 2.1.4

ここまでは点電荷を考えたが、実際には線電荷、面電荷、体積電荷を考えることも多い。ここでは、体積電荷の与える電場の式を載せる。

$$E(\boldsymbol{r})\equiv\frac{1}{4\pi\varepsilon_0}\int\frac{\rho(\boldsymbol{r'})}{\mathscr{r}^2}\boldsymbol{\hat{\mathscr{r}}}d\tau'$$
