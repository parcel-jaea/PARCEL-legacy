[English](README.md) | 日本語

# PARCEL (旧版 1.0 〜 1.2)

このリポジトリは、Fortran (1.2 からは CUDA も) で書かれた並列 Krylov 部分空間法のライブラリ
PARCEL (PARallel Computing ELements) の過去の版を保存するものである。これらの版の開発は終了している。

**現在の PARCEL は SYCL 版 (PARCEL_SYCL) である:** <https://github.com/parcel-jaea/PARCEL-release>
(Release: <https://github.com/parcel-jaea/PARCEL-release/releases>、マニュアル: <https://ccse2.jaea.go.jp/software/parcel/index.html>)。

## 版

各版にタグと GitHub の Release があり、ソースの tarball と SHA-256 のファイルを置いている。
`main` ブランチは最後の版 (1.2) である。各版のファイルはそのタグで見られる。

| 版 | 公開日 | タグ・ダウンロード | マニュアル (HTML) |
|---|---|---|---|
| 1.2 | 2023-01-31 | [v1.2](https://github.com/parcel-jaea/PARCEL-legacy/releases/tag/v1.2) | [日本語](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.2/manual/PARCEL_manual_mathjax.html) / [英語](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.2/manual/PARCEL_manual_mathjax_eng.html) |
| 1.1 | 2022-01-13 | [v1.1](https://github.com/parcel-jaea/PARCEL-legacy/releases/tag/v1.1) | [日本語](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.1/manual/PARCEL_manual_mathjax.html) / [英語](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.1/manual/PARCEL_manual_mathjax_eng.html) |
| 1.01 | 2020-05-20 | [v1.01](https://github.com/parcel-jaea/PARCEL-legacy/releases/tag/v1.01) | [日本語](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.01/manual/PARCEL_manual_mathml.html) / [英語](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.01/manual/PARCEL_manual_mathml_eng.html) |
| 1.0 | 2019-03-29 | [v1.0](https://github.com/parcel-jaea/PARCEL-legacy/releases/tag/v1.0) | [日本語](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.0/manual/PARCEL_manual_mathml.html) / [英語](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.0/manual/PARCEL_manual_mathml_eng.html) |

公開日は配布サイトの過去の版のページ (廃止予定) の記載による。
ビルドと使い方は各版の `readme` を参照。

各版のマニュアルは tarball とタグの `manual/` に入っている (図と数式を埋め込んだ 1 ファイルの HTML で、インターネットから何も読み込まない)。
tarball を展開してウェブブラウザで開く。GitHub 上では HTML はソースとして表示される。

## 配布サイトの tarball との違い

ファイル名と展開先のディレクトリ名は配布サイトにあったものと同じで、違いは次だけである。

- ビルドの生成物と作業用のファイルを除いた (オブジェクト `*.o`、モジュール `*.mod`、実行ファイル、ジョブの出力ファイル、
  1.0 の `SPARSE/src/tmp.f90` と `tmp.sh`)。
- 開発者のアカウントの個人パスを含む行を除いた (いずれもコメントアウトされた行)。1.0 の `SPARSE/example5/Makefile` の
  ライブラリのパスは `<path-to>/libjaeri_ilesf77.a` に、出力例 `SPARSE/example/STDOUT_PARCEL-dev` の作業ディレクトリは `<path-to>/PARCEL-dev` に置き換えた。
- マニュアル (HTML、日本語・英語) を `manual/` に加えた。配布サイトでは tarball と別に配っていたもの。
  1.1・1.2 は、インターネットから読み込む MathJax で表示していた数式を MathML に変換した (本文は変えていない。ファイル名はそのまま)。
- ライセンス文 `COPYING` と `COPYING.LESSER` を加えた。1.0 には空のディレクトリを保つための `MPI/lib/.gitkeep` を加えた。
- ファイルの所有者は記録せず、時刻はすべてその版の日付にそろえた。

## ライセンス

配布サイトと同じ GNU Lesser General Public License v3.0 (LGPL-3.0)。`COPYING` と `COPYING.LESSER` を参照。

## 問い合わせ

質問は PARCEL_SYCL の GitHub の Issues へ: <https://github.com/parcel-jaea/PARCEL-release/issues>。
これらの版は保守していない。新しく使う場合は PARCEL_SYCL を使ってほしい。
