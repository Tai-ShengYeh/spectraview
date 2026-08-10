# 拉曼光譜顏料分析入門

給大一零基礎學生的教學包。同一件事用三條路線各做一遍：
**Python**、**R**、**Orange Data Mining**（不寫程式）。

👉 **[線上閱讀教材](https://tai-shengyeh.github.io/spectraview/raman-tutorial/)**

## 內容

5 篇 26 章、8 張圖表、50 題互動測驗（含學習儀表板與成績單匯出）。
`index.html` 是單一自足檔案，下載後離線也能用。

| 篇 | 內容 |
|---|---|
| 一 | 觀念：光譜在說什麼（不寫程式）——扣基線、訊噪比、鑑定的三道門檻 |
| 二 | Python 路線：ALS 基線、find_peaks、譜庫比對、NNLS 解混 |
| 三 | R 路線：只用 base R + Matrix，自行實作 SG 濾波、prominence、Lawson–Hanson NNLS |
| 四 | Orange 路線：用 [orange-spectra](https://pypi.org/project/orange-spectra/) 的 widget 串流程 |
| 五 | 三個真實案例，以及程式一定會犯的三種錯 |

## 檔案

| 檔案 | 說明 |
|---|---|
| `index.html` | 主教材（單一自足 HTML，約 1 MB） |
| `raman_pipeline.py` | Python 版分析腳本 |
| `raman_pipeline.R` | R 版分析腳本（功能對等，判定結果逐一比對過皆一致） |
| `spectra/` | 10 個教學用光譜（CSV，表頭 `wavenumber,intensity`） |
| `spectra_matrix.csv` | 矩陣格式（10 列 × 801 欄，200–1800 cm⁻¹ 每 2 cm⁻¹），給 Orange 用 |
| `pigment_library.csv` | 內建 15 種顏料的譜庫表 |

## 快速開始

```bash
# Python（需要 numpy scipy matplotlib）
python raman_pipeline.py spectra/ -o results
python raman_pipeline.py spectra/UNK3_blue.csv \
       --refs spectra/PB15_phthalo_blue.csv spectra/rutile_titanium_white.csv
python raman_pipeline.py --list-library

# R（只需 base R + Matrix，Matrix 是官方隨附套件）
Rscript raman_pipeline.R spectra/ -o results
Rscript raman_pipeline.R --list-library

# Orange
pip install orange-spectra
python -m Orange.canvas
```

## 教學用樣品與正確答案

| 檔案 | 應得結果 | 教學重點 |
|---|---|---|
| `cinnabar_1` / `cinnabar_2` | 辰砂 α-HgS | 基本流程、量測重複性 |
| `rutile_titanium_white` | 金紅石 TiO₂ | 晶型辨別（vs 銳鈦礦） |
| `PB15_phthalo_blue` | 酞菁藍 PB15 | 有機顏料的多帶指紋 |
| `PB15_cobalt_blue_hue` | 也是 PB15 | **標籤不等於內容**（cobalt blue hue） |
| `carmine` | 譜庫裡沒有 | **庫外成分**該怎麼處理 |
| `PB15_failed_measurement` | 無法判定 | **資料品質把關** |
| `UNK1_blue` | PB15 + rutile | 年代下限 1938 年 |
| `UNK2_blue` | 普魯士藍 + 苯乙烯樹脂 | 斷代要看黏合劑，不看顏料 |
| `UNK3_blue` | PB15 + PG7 + rutile + 樹脂 | 用殘差分析找出第三種成分 |

三件未知樣為實際採樣資料，已匿名處理。

## 相關資源

- [orange-spectra 完整 widget 教學](https://tai-shengyeh.github.io/spectraview/orange.html) —— 本教材第四篇的延伸
- [spectraview 專案](https://github.com/Tai-ShengYeh/spectraview) —— 套件原始碼與問題回報
- [orange-spectra @ PyPI](https://pypi.org/project/orange-spectra/)

## 授權

雙授權，詳見 [`LICENSE`](LICENSE)。

| 對象 | 授權 | 白話 |
|---|---|---|
| 教材與資料（`index.html`、`spectra/`、`spectra_matrix.csv`、`pigment_library.csv`） | [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.zh-Hant) | 可自由使用與改作，須標示來源、不得商用、衍生作品同條款釋出 |
| 程式碼（`raman_pipeline.py`、`raman_pipeline.R`） | MIT | 想怎麼用就怎麼用，包含商業用途 |

課堂教學、學生自學、學術研究皆屬非商業使用，無須另外取得同意。

建議的標示方式：

> 「拉曼光譜顏料分析入門」，Tai-Sheng Yeh，
> https://tai-shengyeh.github.io/spectraview/raman-tutorial/ ，依 CC BY-NC-SA 4.0 授權。

## 文獻依據

Burgio, L. & Clark, R. J. H. (2001). Library of FT-Raman spectra of pigments, minerals,
pigment media and varnishes. *Spectrochimica Acta Part A*, 57, 1491–1521.
顏料年代取自 MFA CAMEO 材料資料庫。
