# 目的

現在の setuptools + custom build_ext + f2py 構成を廃止し、
**Meson + meson-python ベースのビルドへ移行する**。

---

## 背景

* 現在は setup.py の custom build_ext で `numpy.f2py -c` を呼び出している
* numpy.distutils は廃止方向
* Python 3.12 + NumPy 2.x 環境
* 将来を見据えて Meson ベースに統一したい

---

## ゴール

* `pip install .` でビルド・インストール可能
* wheel ビルド (`pip wheel .`) が成功する
* 拡張モジュール `wind3d_field_lines._bbtobln` が正しく import できる
* setup.py を完全削除する

---

## 現在の構成

* src layout

  * src/wind3d_field_lines/fortran/field_line_integrator.f90
* f2py を使って `_bbtobln` モジュールを生成
* Python package: `wind3d_field_lines`

---

## 要件

### 1. ビルドシステム

* `meson-python` を使用
* `pyproject.toml` の build-backend を変更する

```toml
[build-system]
requires = [
  "meson-python",
  "numpy>=2.0"
]
build-backend = "mesonpy"
```

---

### 2. Meson で Fortran + f2py を扱う

以下のいずれかで実装すること：

#### 推奨

* f2py を使ってラッパー C ファイルを生成
* Meson の `extension_module()` でビルド

#### 代替

* Meson の `custom_target()` で f2py を呼ぶ

---

### 3. モジュール名

最終的に Python から

```python
import wind3d_field_lines._bbtobln
```

で import できるようにすること

---

### 4. ディレクトリ構成

必要に応じて以下を作成：

* `meson.build`（ルート）
* `src/wind3d_field_lines/meson.build`

---

### 5. NumPy 連携

* include path は Meson の Python module を使って取得する
* numpy.get_include() を適切に利用

---

### 6. Python パッケージ連携

* src layout を維持
* pure Python 部分もインストールされるようにする

---

### 7. 削除対象

* setup.py
* custom build_ext

---

## 実装タスク

1. `pyproject.toml` を meson-python 用に書き換え
2. `meson.build` を新規作成
3. Fortran ファイルをビルド対象に追加
4. f2py ラッパー生成 or 直接ビルドの設計
5. Python 拡張モジュールとして `_bbtobln` を出力
6. `pip install .` で動作確認
7. `pip wheel .` で wheel 確認

---

## 注意点

* wheel に `.so` が含まれることを確認
* macOS / Linux 両方で動く設計にする
* Python ABI suffix を Meson に任せる（手動命名しない）
* 並列ビルドを有効化（MesonデフォルトでOK）

---

## 出力形式

以下を提示すること：

1. 完成した `pyproject.toml`
2. `meson.build`（ルート）
3. 必要ならサブディレクトリの `meson.build`
4. 変更されたディレクトリ構成
5. 動作確認手順

---

## 禁止事項

* setuptools を使い続ける
* setup.py を残す
* 手動コピーで .so を配置する
* numpy.distutils を使う

---

## 補足

最小構成でよい。過剰な abstraction は不要。
シンプルで保守しやすい Meson 構成を優先すること。
