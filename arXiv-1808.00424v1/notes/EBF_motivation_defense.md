# 為 EBF 鏟子找泥土 — 對應 spatial-skew-t 專案的方法論辯護

> 對應論文:Morris, Reich & Thibaud (2018) *Exploration and inference in spatial extremes using empirical basis functions* (arXiv:1808.00424v1)
> 對照基準:Morris, Reich, Thibaud & Cooley (2017) *A space-time skew-t model for threshold exceedances*, Biometrics 73(3):749–758
> 配套筆記:[EBF_basis_functions.md](EBF_basis_functions.md)

---

## 0. 為什麼第一眼覺得「為做而做」

讀完 arXiv 版,直觀的不滿大致有三:

1. **模擬實驗自說自話**:資料用 Reich & Shaby (2012) 的正穩定模型(也就是 EBF 配適的模型)生成,然後配適同一個模型 — 沒測 misspecification。
2. **實證對 GKF 改善有限**:current 期 MAD 兩種方法幾乎打平,只有 future 期才見差距(第 405 行)。
3. **storm 物理比喻牽強**:把 $B_l(\mathbf{s})$ 解釋成「第 $l$ 種風暴的空間範圍」,接近事後敘事 — 從未檢驗各 EBF 的物理顯著性。

這些都成立。但是,如果只停在這層直覺,會錯過一條重要的線索 — **EBF 並非憑空的方法論演習,而是同一群作者在 Morris 2017 *自己的 future work* 段落中許下的承諾**。下面把這條線拉直,並且具體對應到你正在做的 ozone / spatial-skew-t 工作。

---

## 1. 鏟子的精準形狀 — EBF 真正獨特的三件事

剝掉 PR 式的包裝,EBF 只貢獻三件其他方法不能同時給的事:

### (a) 資料驅動的 *非平穩* 相依結構,不必指定 $\rho(\mathbf{s})$ 的函數型
傳統 Matérn / Brown–Resnick 要非平穩,得寫 $\rho(\mathbf{s})$、$\nu(\mathbf{s})$ 的參數族(常用 thin-plate spline、deformation),且必須先想清楚要如何隨地理變化。EBF 直接從成對 extremal coefficient 平滑後逆解 $B_l(\mathbf{s})$,**型態完全由資料決定**,沒有手寫的非平穩函數族。

### (b) 空間極值分析的第一個 EOF/PCA 類比
氣象/海洋科學的家常工具是 EOF。對 *相依* 結構而言(不是 mean field),極值版的 EOF 之前根本沒有 — 沒有可以畫的圖、沒有 "$v_l$ 排序" 之類能放進報告的東西。$v_l = \frac{1}{n_s}\sum_i \hat B_{il}$ 與 $\sum v_l = 1$ 提供了一個可解讀、可排序、可可視化的相依結構分解 — 這是真實的方法論缺口被補上。

### (c) 邊際 / 相依結構模組化解耦
Reich & Shaby (2012) 把所有參數一起塞進 MCMC,大資料下吃不消。EBF 把 $\hat{\mathbf{B}}, \hat\alpha$ **一次性 plug-in 固定**,MCMC 只跑 $A_{lt}$ 與邊際參數 — 10000 次 update ≈ $2 \times L \times n_t$ 秒。對於 $n_s$ 上千的資料這是質的差距,而非常數加速。

---

## 2. 泥土在哪裡 — 對應你 spatial-skew-t 專案的四個具體節點

這一節是辯護的主幹。

### (a) Morris 2017 *自己* 的 future work 段落就是 EBF

打開你專案的 [LaTeX/skewt_rev_2.tex](../../LaTeX/skewt_rev_2.tex),Discussion 第 581–583 行是這樣寫的:

> *"One possibility is the implementation of a different partition structure. We choose to define the random effects for a site by using an indicator function based on closeness to a knot. However, this indicator function could be replaced by **kernel function that would allow for multiple knots to impact each site, with the weight of each knot to be determined by some characteristic such as distance**."*

把這段話與 EBF 模型 (eq. 4 / [EBF_arxiv_v1.tex:273-275](../EBF_arxiv_v1.tex#L273-L275))

$$
\theta_t(\mathbf{s}) = \Big\{\textstyle\sum_{l=1}^L B_l(\mathbf{s})^{1/\alpha} A_{lt}\Big\}^{\alpha},\qquad B_l(\mathbf{s})\ge 0,\ \sum_l B_l(\mathbf{s})=1
$$

對齊一下 — $B_l(\mathbf{s})$ 就是 Morris 2017 想要的「site-to-knot 軟權重」,而且額外滿足非負、和為 1、由資料估出(而非預先指定 Gaussian 核)。**EBF 不是無泥土的鏟子,它是同作者群在 Morris 2017 留下的待辦事項上的後續實現**。從這個角度看,稱它「為做而做」其實是低估了它在這條研究脈絡上的位置。

### (b) Ozone 分析隱性的平穩假設可以被釋放

你的 skew-t 實作中 ([code/R/auxfunctions.R:160-186](../../code/R/auxfunctions.R#L160-L186)) `CorFx()` 是這樣的:

```r
CorFx <- function(d, gamma, rho, nu) {
  ...
  cor <- gamma * simple.cov.sp(D = d, sp.type = "matern",
                               sp.par = c(1, rho), ...)
  ...
}
```

`rho`、`nu`、`gamma` 是 **全域常數**,意味著協變異結構在整個美國本土被假設為平穩。但實際上 ozone 形成機制隨地理大幅變動 —
- 沿海邊界層動力 vs. 內陸自由對流
- 山區 vs. 平原(EBF 論文 fig. EBFs10current 中第一個 EBF 也對應到 Appalachian Mts.)
- 都會排放源(NYC、Atlanta)vs. 鄉村背景值

`CorFxDef` 雖然加了 geometric anisotropy,但仍是全域單一變換。**EBF 提供的價值是:給你一條「不必事先指定 $\rho(\mathbf{s})$ 函數族就能得到非平穩相依結構」的路徑**。對於你的 ozone 分析,這恰好是 Matérn 假設最薄弱的地方。

### (c) MRTS basis 目前只當協變量 — EBF 開啟「basis 也當相依結構」的可能

你在 [code/analysis/ozone/US-all-auto/0428temp.md](../../code/analysis/ozone/US-all-auto/0428temp.md) 對 MRTS basis $F$ 的處理是:

> *「決定 $F$ 是完全固定的... $F$ 固定,代表我們不把 basis construction 的不確定性塞進推論。」*

也就是說 basis function **只進到 mean function 當回歸協變量**,沒有進入相依結構。EBF 提供的視角是:同樣的「以少量空間 basis 控制空間複雜度」思維,可以**升一層用到相依結構本身**(即 $\theta_t(\mathbf{s})$ 的構造),而不只是 mean。對你來說,這是把現有 MRTS 經驗推到 covariance 端的自然路徑 — 但要注意它仍是 plug-in、固定的(同 Morris 2017 不傳基底不確定性的選擇,EBF 也不傳),所以方法論上的概念距離比想像中小。

### (d) 1089 站點的 scaling 痛點

你目前的設定是把 1089 個 AQS 站點 *次取樣到 800 站*(見 `US-ALL-RUN-GUIDE.md`)以求 MCMC 可行。原因之一是 skew-t 的相依結構吃 $n_s \times n_s$ Matérn 矩陣與 partition 計算。**EBF 模型的 Bayesian 階段計算量是 $O(L \times n_t)$**,$L$ 是 basis 個數(實證建議 10–15),完全與 $n_s$ 解耦。如果哪天你想處理整個美國本土上千甚至萬點的衛星 ozone retrievals 或 CMAQ 全格網,EBF 的低秩結構是少數可以開上去而不必子採樣的路徑之一。

---

## 3. 但鏟子也有瑕疵 — 誠實點名,特別對 ozone 場景

辯護不等於背書。下面是讀完論文後對 EBF 的具體保留:

### (i) 模擬實驗自說自話
真值 $B_l$ 用 Gaussian kernel (eq. 5) 在 $\sqrt L \times \sqrt L$ 格點生成,然後用同模型 EBF 配適。**沒測 model misspecification** — 如果真實過程是 Brown–Resnick 或 Smith,EBF 還能恢復多少?論文沒答。

### (ii) 未與真正競爭者比較
比較對象只有 GKF(同個模型族、不同 basis 選擇)。沒比 Brown–Resnick、非平穩 Matérn、Elliptical Pareto、或 Bernard (2013) 的 clustering 方法。所謂的「方法論優勢」只是同族內的內部比較。

### (iii) **對 ozone 場景的關鍵張力**:asymptotic independence
這點對你最重要。Morris 2017 之所以選 skew-t 而 *不選* max-stable,原因就在你論文 §3.1–3.3 反覆強調的:**ozone 在大距離展現 asymptotic independence**,而 max-stable 強制 asymptotic dependence 不會隨距離消失。EBF 完全建立在 max-stable (正穩定隨機效應) 架構上 — 把 EBF 直接套到 ozone 殘差,**就是把 Morris 2017 一開始就拒絕的東西放回來**。要在 ozone 上用 EBF 的精神,得先把它從 max-stable 框架移植到 skew-t 框架(等於是另一篇論文)。這個張力 EBF 論文自己沒處理,因為他們的應用是 *降水年最大值*,asymptotic dependence 在區域內合理。

### (iv) $T = 31$ 天的 basis 估計穩定性
EBF 第一步是 F-madogram 估 $\hat\vartheta_{ij}$,小樣本(你只有 31 天)雜訊很大。然後第二步用梯度下降在受約束流形上解,小樣本下估出來的 $\hat B_l$ 可能對哪幾天落入訓練集很敏感 — 等於把樣本誤差搬到 plug-in 的相依結構裡且不傳遞給 MCMC 後驗。論文用 $n_t \in \{50, 200\}$ 模擬,$n_t = 50$ 已經吃力,$n_t = 31$ 風險更高。

### (v) 正穩定密度仍需數值近似
即使 EBF 把 $\hat{\mathbf{B}}, \hat\alpha$ plug-in,$A_{lt}$ 的後驗仍需在每次 MCMC update 用中點法 + 50 個 Beta(0.5, 0.5) quantile 對 PS 密度作一維數值積分(附錄 A,[EBF_arxiv_v1.tex:464-478](../EBF_arxiv_v1.tex#L464-L478))。這不是乾淨的解析 MCMC,而且離散化誤差會影響後驗 calibration。

### (vi) current data 上 GKF ≈ EBF
論文 fig 1 自己顯示 current 期兩者 MAD 幾乎重疊。如果真實場景下的相依結構沒有特別非平穩 / 不規則,EBF 的資料驅動優勢不會顯現,卻多了一層 plug-in 不確定性。

---

## 4. 結論 — 何時動鏟,何時放下

**動鏟子的時機**
- 大空間維度($n_s \gtrsim 10^3$)、適中時間樣本($n_t \gtrsim 100$)
- 過程合理地接近 asymptotic dependence(年最大值、強對流降水、極值熱浪)
- 需要 *視覺化* 報告相依結構(政策、保險、氣候 attribution)
- 跨期 / 跨子集比較相依結構演化(climate-change-style comparison 是 EBF 的甜蜜點)

**放下鏟子的時機**
- 資料展現大距離 asymptotic independence(**ozone 中尾、PM2.5、許多環境化學變量**)
- $n_t$ 過短(< 50),F-madogram 不可靠
- 有強烈物理理由相信平穩 + 已知函數型 $\rho(\mathbf{s})$,且不需要視覺化
- 邊際與相依需要聯合不確定性傳遞(EBF 的 plug-in 結構在此先天劣勢)

**對你的專案的具體建議**

1. **不要直接搬 EBF 到 ozone 殘差** — asymptotic dependence 性質與 skew-t 動機矛盾。
2. **可以借的是 basis 思維** — 把 Morris 2017 的硬分配 partition 換成 site-to-knot 軟權重 $B_l(\mathbf{s})$(這就是 [skewt_rev_2.tex:581-583](../../LaTeX/skewt_rev_2.tex#L581-L583) 許下的承諾),但 $B_l$ 的估計與意義要重新放在 skew-t 框架下推導,不能直接抄 F-madogram 那條路。
3. **EBF 的 $v_l$ 排序視覺化** 倒是可以原樣借用 — 即使在 skew-t 框架下,若你有對應的 $B_l(\mathbf{s})$,把它畫出來搭配 $v_l$ 的排序,就是現有 skew-t 文獻沒有的探索工具。

EBF 對你不是現成的工具,而是一個方法論方向標 — 它示範了「資料驅動的軟空間基底」這條路可以走,缺的是把它從 max-stable 搬到 skew-t 的橋。
