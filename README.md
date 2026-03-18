# master-thesis-appendix

修士論文に関する附属のプログラムを公開するためのレポジトリです.

## BlockMaps.jl

subshift of finite type (SFT) の間のシフト可換連続写像 (sliding block code) を扱うライブラリです。
`BlockMap` の定義・評価・合成・逆写像計算に加えて、`extend`/`shorten` による表現の調整、`StrongShiftEquivalence` の計算、dimension group representation や orbit sign number などの関連量の計算機能を提供します。

## cascade-count.jl

フルシフト $\Sigma_2=\{1,2\}^{\mathbb{Z}}$ から与えられた SFT を除いた系における $(2,n)$-cascade の個数を計算するプログラムです。
forbidden words (ワイルドカード `"*"` を含む) を入力として受け取り、`count_cascade`, `gcd2_cascade`, `even1_cascade` の 3 つの方法で $n=1,\dots,m$ の個数列を計算できます。

## subshift-lattice.jl

フルシフト $\Sigma_2=\{0,1\}^{\mathbb{Z}}$ 上の複数の SFT について、包含関係を計算するプログラムです。
包含判定演算子 `⊏`, `⊐` に加えて、`shift_hasse_diagram` により与えた SFT 集合の Hasse diagram を可視化できます。

## Pruning.jl

保存系の実 H\'enon 写像に対する primary pruned region 推定アルゴリズム (Hagiwara-Shudo) の Julia 実装です。
双曲的不動点・安定/不安定多様体・横断的ホモクリニック点の計算、symbolic encoding、primary pruned region のブロック抽出、forbidden words と markerized 表現の生成までを一連のワークフローで実行できます。
