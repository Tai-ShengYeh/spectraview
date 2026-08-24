# 範例光譜庫：糖類與食品添加物（NIR）

`sugars_nir.speclib` 是一個現成的**光譜庫範例**，用來示範 SpectraView 的「光譜庫相似度比對」。

## 內容
- **9 種參考物質**：aspartame、benzoic acid、caffeine、fructose、glucose、lactose、
  maltose、sucralose、sucrose。
- 每條參考譜是該物質 **5 次重複量測的平均**（DLP Hadamard 近紅外，1600–2400 nm、200 點、吸光度）。
- 另附 3 個單次量測的「未知譜」查詢檔：`unknown_glucose.csv`、`unknown_sucrose.csv`、
  `unknown_caffeine.csv`。

## 怎麼用（在 Orange 裡）
orange-spectra **0.7.2 起這個庫已內建於套件**：Spectral Library widget 按
**Add built-in**、下拉選「Sugars & food additives (NIR, Hadamard)」即可一鍵載入，
不需要這個檔案。詳見[教學頁](https://tai-shengyeh.github.io/spectraview/orange.html#library)。

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

---

## CGL NIR — Mixture Analysis (NNLS) 真實範例

`cgl_mixture_nnls.py` 用 Eigenvector 的 **CGL 三成分 NIR 混合設計**
（casein／glucose／lactate，1104–2496 nm；來源：<https://eigenvector.com/resources/data-sets/>）
示範 **Mixture Analysis (NNLS)** 拆解混合譜：

1. 用校正集以古典最小平方（CLS）估出三個成分的**純光譜**（`S = pinv(C)·X`）；
2. 對一條**測試混合譜**做 NNLS（`orangespectra.core.mixture_nnls`）拆解；
3. 回推的比例與參考 wt% 比對。

**執行**（需要能連到 eigenvector.com，或自己下載 `CGL_nir.mat` 後傳路徑）：
```bash
python examples/cgl_mixture_nnls.py [path/to/CGL_nir.mat]
```
會產生 `docs/demo_mixture.png`，並輸出可直接在 Orange 裡試的兩個檔：
- `cgl_components.speclib` — 3 條純成分參考譜
- `cgl_mixture.csv` — 一條測試混合譜

**在 Orange 裡實跑**：用 **Load Spectra Files** 載入 `cgl_components.speclib`（接
Mixture Analysis 的 *References*），再載入 `cgl_mixture.csv`（接 *Mixture*），即得成分比例與 R²。

> `CGL_nir.mat` 只在 Eigenvector 網站，不隨 repo 附帶（已 gitignore）。
