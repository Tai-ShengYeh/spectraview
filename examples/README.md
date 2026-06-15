# 範例光譜庫：糖類與食品添加物（NIR）

`sugars_nir.speclib` 是一個現成的**光譜庫範例**，用來示範 SpectraView 的「光譜庫相似度比對」。

## 內容
- **9 種參考物質**：aspartame、benzoic acid、caffeine、fructose、glucose、lactose、
  maltose、sucralose、sucrose。
- 每條參考譜是該物質 **5 次重複量測的平均**（DLP Hadamard 近紅外，1600–2400 nm、200 點、吸光度）。
- 另附 3 個單次量測的「未知譜」查詢檔：`unknown_glucose.csv`、`unknown_sucrose.csv`、
  `unknown_caffeine.csv`。

## 怎麼用（在 SpectraView 裡）
1. `Library ▸ Load library…` 選 `examples/sugars_nir.speclib`
2. `File ▸ Open…` 載入一個 `unknown_*.csv`（或你自己的 NIR 譜）
3. 在左側清單選那條未知譜 → `Library ▸ Search selected against library…`
4. 命中清單依**相關係數**排序，最相近的物質排第一；可按「Overlay top hit」把參考譜疊上去比對。

## 預期結果與教學重點
- 每個 `unknown_*` 都會正確命中對應物質（top-1 相關係數 ≈ 1.00）。
- 相關係數會說一個**化學上合理的故事**：糖類彼此相似（蔗糖 vs 果糖 ≈ 0.99，因為蔗糖
  就是葡萄糖＋果糖），但和**添加物**（咖啡因、苯甲酸）明顯區隔（≈ 0.7–0.87）。
- 對比同一份命中清單裡的 cosine（全擠在 0.95–1.0）可以看出：**相關係數的鑑別力更好**，
  這也是為什麼預設用相關係數排序。

## 重建方式
此庫由 `*_a.csv`（每物質 5 重複的吸光度譜）以「每物質取平均」建成。資料來源為
使用者自有量測（DLP Hadamard NIR）。
