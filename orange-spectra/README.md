# orange-spectra — Orange Data Mining 光譜 widgets

五個給 [Orange Data Mining](https://orangedatamining.com/) 的光譜學 widgets，
與 [SpectraView](https://github.com/Tai-ShengYeh/spectraview) 共用演算法與
`.speclib` 光譜庫格式。Five spectroscopy widgets for Orange Data Mining,
sharing algorithms and the `.speclib` library format with SpectraView.

| Widget | 功能 |
|---|---|
| **Import Spectrum URL** | 貼上 **IRUG 編號／網址**或 **SOPRANO 網址**（也支援 JCAMP-DX/CSV 直接網址），下載並畫出光譜，輸出成 Orange Table（每列一條光譜、欄名＝波數） |
| **Spectra Similarity** | 兩組光譜間的相似度：correlation / cosine / 光譜角 SAM / Euclidean |
| **Spectral Library** | 建立參考光譜庫、存成 **`.speclib`（與 SpectraView 互通）**、對庫比對未知譜並輸出排名 |
| **Mixture Analysis** | 混合光譜的成分分析：非負最小平方（NNLS）解 mixture ≈ Σ cᵢ·refᵢ，回報係數、比例與 R² |
| **Aquagram** | Aquaphotomics：在水的 12 個特徵吸收帶（WAMACs）取正規化吸光度，畫 12 軸雷達圖（raw / SNV / aquagram 三種正規化） |

輸出的 Table 採 [Orange-Spectroscopy](https://orange-spectroscopy.readthedocs.io/)
慣例（欄名＝波長/波數、每列一條光譜），可直接接其 Spectra 檢視 widget 或
PCA / PLS 等化學計量學流程。

## 安裝 Install

> ⚠️ **先確認你的 Orange 是哪一種**——桌面版 App 與 pip 版是不同的 Python 環境，裝錯不會出現 widgets。

**A. 桌面版 App**（orangedatamining.com 下載的獨立程式）：
`Options ▸ Add-ons… ▸ Add more…` 輸入 **`orange-spectra`**（只吃 PyPI 名稱，**不吃** git 網址）→ 打勾 → OK → 重啟。
需先發佈到 PyPI（見 [`PUBLISHING.md`](PUBLISHING.md)）；發佈前可用 App 的 Python 執行下面 pip 指令裝進 App 環境。

**B. pip 版 Orange**（`python -m Orange.canvas` 啟動）：
```bash
pip install orange3 PyQt5 PyQtWebEngine   # Orange + Qt 綁定（缺 Qt 會開不了）
pip install "git+https://github.com/Tai-ShengYeh/spectraview.git#subdirectory=orange-spectra"
python -m Orange.canvas                    # 用這個開，別點 App 圖示（是另一個環境）
```

重新啟動 Orange，工具箱會出現 **Spectra** 分類（5 個 widgets）。
Restart Orange; a **Spectra** category (5 widgets) appears in the toolbox.

常見錯誤：`ImportError: PyQt5 … not available` = 少裝 Qt，補 `pip install PyQt5 PyQtWebEngine`。

## 快速上手 Quick start

1. 拖出 **Import Spectrum URL**，輸入 `4119`（IRUG 的 PB15 酞菁藍 Raman 譜）→ Fetch。
2. 多抓幾條參考譜 → 接 **Spectral Library** 的 *Spectra* 輸入 → *Add input spectra to library* → *Save…* 存成 `.speclib`。
3. 未知譜接 Library 的 *Query* 輸入 → *Hits* 輸出就是排名表。
4. 混合譜接 **Mixture Analysis** 的 *Mixture*、參考譜（或 Library 的 *Library* 輸出）接 *References* → 得成分比例與 R²。

## 測試 Tests

```bash
python tests/test_core.py                                   # 純邏輯（免 GUI）
QT_QPA_PLATFORM=offscreen python tests/test_widgets.py      # widget 煙霧測試
```

## 注意 Notes

- URL 匯入支援 IRUG 詳情頁（jqPlot 內嵌資料）、SOPRANO 頁（Dygraph 內嵌資料）、
  JCAMP-DX（AFFN 純數字格式）與兩欄 CSV/TSV。壓縮的 JCAMP（SQZ/DIF）請改用
  SpectraView 開啟後匯出。
- 教學網頁：<https://tai-shengyeh.github.io/spectraview/orange.html>
