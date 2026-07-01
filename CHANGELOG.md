# 更新說明（Changelog）

SpectraView 的更新紀錄，最新在上。每筆變更的完整技術細節見
[GitHub commit 紀錄](https://github.com/Tai-ShengYeh/spectraview/commits/main)。

格式參考 [Keep a Changelog](https://keepachangelog.com/)；目前以日期分節（尚未發行語意化版本號）。

## 未發行（Unreleased）

> 已在 `main` 分支、但尚未納入 [v0.2.0](https://github.com/Tai-ShengYeh/spectraview/releases/tag/v0.2.0) 發行標籤的變更。

### 新增（Added）
- **線上匯入（File ▸ Import from URL / IRUG…）**：直接從 [IRUG 光譜資料庫](http://www.irug.org)
  或任何網址抓光譜進來。輸入可以是 **IRUG 編號**（如 `4119`）、**IRUG 詳情頁網址**，
  或**直接的 JCAMP-DX／CSV 檔網址**；程式會下載並自動判斷來源：
  - **IRUG 詳情頁**：資料是寫在頁面互動圖（jqPlot）的 `<script>` 裡，為一串
    `"波數:強度"` 配對——直接解析還原成 x／y（波數軸 cm⁻¹），技術別（Raman／IR）
    用來判斷 y 單位。
  - 其他網址：JCAMP-DX 走 ASDF 解析、CSV/TXT 走 ASCII；詳情頁也會回退去找
    頁面內嵌的 JCAMP 或下載連結（含無副檔名的下載端點）。
  來源網址會記進光譜 meta。只用 Python 標準庫（urllib），免裝額外套件。
- 讀取**轉置／矩陣式 CSV**（波長軸在表頭、每一列是一條光譜）：自動辨識外部工具
  （如 Orange）與化學計量學 ML 匯出的矩陣檔，以及本程式自己的 `layout='rows'`
  合併匯出，逐列載入成多條光譜。支援多個前導中繼欄（如 Sample ID、Concentration），
  其值會帶進光譜名稱與 meta。

## 2026-06-24

### 新增（Added）
- **PerkinElmer `.sp` 讀取器**：解析 PerkinElmer「PEPE」標籤區塊二進位格式（FTIR / UV-Vis），
  自動對映座標與訊號單位。純 Python，免裝套件。
- **Shimadzu UV-Vis `.SPC` 讀取器**：`.spc` 副檔名會自動辨識 GRAMS 與 Shimadzu 兩種格式；
  Shimadzu（如 UV-1900）為純 Python 原生解析（120-byte 表頭 + float32 資料）。
- **PerkinElmer ASCII `.asc`（PEDS）讀取器**：直接讀檔案 `#GR` 區塊明寫的座標單位
  （如 CM-1 / A）與 `#DATA` 資料，不再從自由文字猜測。

### 修正（Fixed）
- 收緊 ASCII 單位判斷：`um`（微米）只在獨立單字或 `µm` / `micron` 時才成立，
  修掉標頭含 `Spectrum` / `aluminum` 時被誤判為微米的問題。

## 2026-06-23

### 新增（Added）
- **檢量線（calibration curve）**：用已知濃度的標準品建迴歸、回推未知樣品濃度。
  支援線性／過原點／二次模型，回報斜率、截距（含標準誤）、R²、殘差標準誤 s(y/x)、
  偵測極限 LOD、定量極限 LOQ，以及每個未知樣品的 95% 信賴區間（Miller & Miller）。
- **NeoSpectra `.Spectrum` 讀取器**（NIR、tab 分隔匯出）。
- 游標停在曲線上時浮現該光譜名稱。
- View 選單「顯示／隱藏左側光譜清單」（F9）。

## 2026-06-15

### 新增（Added）
- **螢光 EEM（2D / 3D）** 與 **2D 相關光譜（2D-COS / 2T2D）**。
- **EEM PARAFAC 分解** 與 **2D-COS 異相關（hetero-correlation）**。
- 範例光譜庫：9 種糖類／添加物（NIR）；教學頁新增互動式光譜庫比對 demo。

## 2026-06-13 ～ 06-14

### 文件（Docs）
- 新手安裝指南：Windows（git + Python 3.12，winget）、macOS（Homebrew）。

## 2026-06-12（初版）

### 新增（Added）
- **SpectraView 初版**：受 [SpectraGryph](https://www.effemm2.de/spectragryph/) 啟發、
  以 Python（PySide6 + pyqtgraph）寫的桌面光譜檢視與處理程式。
- 載入、疊圖檢視、座標軸換算（X / Y）；多格式 IO（ASCII / JCAMP-DX / JSON / MATLAB /
  GRAMS SPC / Bruker OPUS）。
- 前處理：平滑（Savitzky-Golay、移動平均）、基線（Rubberband／多項式／ALS／airPLS）、
  微分、正規化、散射校正（SNV／MSC／detrend）、光譜運算。
- 分析：尋峰（含 FWHM）、峰形擬合／去卷積（Gaussian／Lorentzian／Pseudo-Voigt）、
  區間積分、混合物 NNLS、XRF 元素鑑定。
- 光譜庫相似度搜尋（相關係數／cosine／SAM／歐氏距離）。
- MIT 授權 + CI；教學頁（`docs/index.html`）。
