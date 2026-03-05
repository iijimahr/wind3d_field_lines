# wind3d_field_lines

wind3dデータの磁力線トレーサの初期実装です。

## 公開API（現段階）

- `trace_field_lines`: 磁力線トレース
- `compute_open_field_fraction`: 開放磁場の面積充填率計算
- `map_field_lines_to_height`: 観測高さへのフットポイント写像

## 開発環境セットアップ

```bash
pip install -e .[dev]
```

## テスト

```bash
pytest
```
