# SpectraView

一個受 [SpectraGryph](https://www.effemm2.de/spectragryph/) 啟發、用 Python 寫的**桌面光譜檢視與處理程式**。
專為食品分析 / 儀器分析教學與研究設計，可載入、疊圖檢視、換算座標軸、做前處理的
FTIR、Raman、UV-Vis、NIR、XRF/LIBS 等光譜。

![demo](docs/demo.png)

> 上圖：同一條 FTIR 光譜的「原始 / 基線校正 / Savitzky-Golay 平滑」三條疊圖，X 軸依 IR 慣例由高波數往低波數。

---

## 安裝

需要 Python 3.10+。

```bash
cd spectra_viewer
pip install -r requirements.txt
```

核心套件：`numpy`、`scipy`、`PySide6`、`pyqtgraph`、`matplotlib`。
二進位格式（SPC、Bruker OPUS）為選用套件，未安裝時程式會給清楚的安裝提示而非崩潰。
（`spc-spectra` 的舊式建置會在隔離環境找不到 numpy，請用
`pip install --no-build-isolation spc-spectra`。）

## 執行

```bash
python run.py                 # 開啟程式
python run.py file1.dx a.csv  # 開啟並直接載入這些光譜
```

或在程式啟動後，直接把檔案**拖放**到視窗，或用 **File ▸ Load demo spectra** 載入內建示範光譜立即試玩。

---

## 功能（v1）

### 檢視核心
- **多格式載入**：ASCII（CSV/TXT/DAT/TSV/PRN）、JCAMP-DX（.dx/.jdx，含 ASDF 壓縮）、
  GRAMS SPC（.spc）、Bruker OPUS（.0/.1…）。一個檔可含多條光譜。
- **疊圖檢視**：多光譜同圖、各自顏色與圖例；左側清單可勾選顯示/隱藏、雙擊改色、改名。
- **互動**：滑鼠縮放/平移、自動縮放（Ctrl+0）、十字游標即時讀出座標（狀態列）。
- **座標軸換算**
  - X：波長 nm ↔ 波數 cm⁻¹ ↔ 波長 µm ↔ Raman 位移 cm⁻¹（需雷射波長）↔ 能量 eV ↔ 頻率 THz
  - Y：穿透率 ↔ 穿透率% ↔ 吸光度 ↔ 反射率 ↔ Kubelka-Munk ↔ log(1/R)
- **顯示選項**：堆疊位移（Scale & Shift）、翻轉 X 軸、格線、對數 Y、深色背景。
- **匯出**：圖檔 PNG（WYSIWYG）與 SVG / PDF / EPS 向量圖（publication 等級）；
  單一光譜存成 CSV 或 JCAMP-DX；**合併匯出**（File ▸ Export combined data）把多條
  光譜放到共用波長軸後輸出成一個 CSV——可選「波長欄 + 每條一欄」（可再載回）或
  「每條一列的 X-matrix」（直接餵 sklearn PLS/PCA）。

### 前處理（可作用於選取的光譜，未選取則套用全部）
- **平滑**：Savitzky-Golay、移動平均
- **基線校正**：Rubberband（凸包）、多項式（ModPoly 疊代）、非對稱最小平方 ALS（Eilers）
- **微分**：1–4 階（Savitzky-Golay）
- **正規化**：峰值=1、Min-Max(0..1)、面積=1、向量(L2)、指定 x₀ 處=1
- **散射校正**：SNV、MSC（對平均光譜）、detrend
- **光譜運算**：兩條相加/減/乘/除（自動對齊重疊區）、平均（mean/median）
- **轉換**：裁切 x 範圍、重採樣到均勻格點
- **復原/重做**：每步處理皆可 Undo/Redo（Ctrl+Z / Ctrl+Y）

### 分析（v2，Analyze 選單）
- **尋峰**：自動偵測峰（或谷）、量測 FWHM、近似面積；在圖上標峰位、並彈出
  可匯出 CSV 的峰表。可調最小高度/突出度/間距、選擇性預平滑。
- **峰形擬合 / 去卷積**：以**自動尋峰當初值**，擬合 Gaussian / Lorentzian /
  Pseudo-Voigt 的疊加（可選同時擬合線性基線）；圖上疊出各分量（虛線）與加總、
  附 R² 與每峰中心/高度/FWHM/面積的表；一鍵把**擬合分量存回光譜清單**。
- **區間積分**：指定 x 範圍、可減端點線性基線，算面積與重心，並在圖上標示積分區。

---

## 專案結構

```
spectra_viewer/
├── run.py                  啟動入口
├── requirements.txt
├── docs/demo.png
├── tests/
│   ├── smoke_test.py       核心邏輯自我測試（無需 GUI）
│   └── gui_smoke_test.py   headless（offscreen）GUI 測試
└── specview/
    ├── app.py              QApplication 設定 + 進入點
    ├── spectrum.py         Spectrum 資料模型 + SpectrumSet 文件
    ├── axes.py             座標軸換算（X、Y）
    ├── processing.py       所有前處理演算法
    ├── analysis.py         尋峰 / 峰形擬合 / 積分
    ├── demo.py             內建示範光譜
    ├── formats/            檔案 IO
    │   ├── __init__.py     依副檔名派發載入 + 存檔
    │   ├── ascii_io.py     ASCII 解析（偵測分隔符/表頭/多欄/歐規小數）
    │   ├── jcamp.py        自寫 JCAMP-DX 解析（AFFN/PAC/SQZ/DIF/DUP）
    │   └── binary_io.py    SPC / OPUS（選用套件，優雅降級）
    └── ui/                 PySide6 + pyqtgraph 介面
        ├── main_window.py  主視窗、選單、工具列、清單、復原
        ├── plotview.py     繪圖區、十字游標、堆疊、匯出
        └── dialogs.py      通用參數對話框
```

## 測試

```bash
python tests/smoke_test.py       # 模型 / 換算 / 處理 / 分析 / IO / JCAMP 解碼
python tests/gui_smoke_test.py   # 建視窗 / 載入 / 處理 / 換算 / 分析 / 匯出（offscreen）
```

---

## 設計說明

- **資料模型**：每條光譜 x 一律以遞增排序儲存；顯示端可獨立翻轉（IR/Raman 慣例由高往低）。
- **每個處理函式都回傳新的 Spectrum**，不改動輸入，配合文件層快照實作 Undo/Redo。
- **JCAMP-DX 自己解析**，不依賴第三方套件，支援真實 FTIR 檔常見的 DIF/SQZ/DUP 壓縮。
- **向量匯出走 matplotlib**（自帶字型、可輸出 PDF/EPS），比 pyqtgraph 的 SVG 寫出穩定。

## 尚未納入（下一步）

v1（檢視 + 前處理）與 v2（尋峰 / 峰形擬合 / 積分）已完成。後續可加：
- **檢量線濃度**：用峰高或峰面積對已知濃度建線、回推未知樣品
- **螢光 EEM**：激發-發射矩陣等高線 / 3D 圖、瑞利散射移除
- **資料庫比對**：光譜庫建立與相似度搜尋
- **批次自動化**：把一串處理步驟套用到整批檔案
- **化學計量建模**：PLS / PCA（接 scikit-learn，配合「合併匯出」的 X-matrix）
