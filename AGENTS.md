# AI Workflow Rules for Julia Lecture

このプロジェクトにおけるAIアシスタント（Agents）の行動指針とルールを以下に定めます。

## 1. コーディング規約 (Code Style)
- **Julia Style Guide準拠**: 公式の [Julia Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/) に従うこと。
- **命名規則**:
  - 変数・関数名: `snake_case` (例: `calculate_mean`)
  - 型・モジュール名: `UpperCamelCase` (例: `DataProcessor`)
- **型注釈**: 関数の引数には可能な限り具体的な型注釈を行うこと。

## 2. 品質保証 (Quality Assurance)
- **型安定性**: 型不安定なコードを避けること。
- **テスト**: `Test` パッケージを使用したテストを含めること。
- **ドキュメント**: 公開関数にはDocstringを記述すること。

## 3. AIの振る舞い
- コード変更の理由を説明すること。
- 可読性を最優先すること。
