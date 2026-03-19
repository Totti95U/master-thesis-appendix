# Appendix of my master thesis

修士論文の Appendix に記載したプログラムのドキュメントページです。

## プログラム一覧

公開しているプログラムは以下の通りです。

- [BlockMaps.jl](@ref guide_of_BlockMaps.jl)
- [cascade-count.jl](@ref guide_of_cascade-count.jl)
- [subshift-lattice.jl](@ref guide_of_subshift-lattice.jl)
- [Pruning.jl](@ref guide_of_Pruning.jl)

## 導入方法

GitHub リポジトリ [https://github.com/Totti95U/master-thesis-appendix](https://github.com/Totti95U/master-thesis-appendix) の `src` ディレクトリに、各プログラムの `.jl` ファイルがあります。これらのファイルをダウンロードして、`include` 関数を使用して Julia で読み込むことができます。
`src/examples` ディレクトリに、各プログラムの使用例が記載された `.jl` ファイルがありますので、そちらも参考にしてください。

各プログラムの概要については、[ガイド](@ref)を参照してください。

各プログラムを使用するに必要なパッケージは以下の通りです。`Pkg.add` 関数を使用して、必要なパッケージをインストールしてください。

- BlockMaps.jl: `LinearAlgebra`, `Permutations`
- cascade-count.jl: なし
- subshift-lattice.jl: `CairoMakie`, `Graphs`, `GraphPlot`
- Pruning.jl: `CairoMakie`, `LinearAlgebra`, `SpatialIndexing`, `ProgressMeter`
