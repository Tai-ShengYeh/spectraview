# 儀器連線指南 (SpectraView ▸ 儀器 ▸ 光譜儀擷取…)

## Ocean Optics / Ocean Insight(USB2000、USB4000、QE Pro 等)

```bash
pip install seabreeze
seabreeze_os_setup      # 一次性:Windows 裝 WinUSB 驅動說明 / Linux 裝 udev 規則
```

連線 → 積分時間/平均次數 → ▶ 即時預覽。「存為暗背景 D」(蓋住入光口)、
「存為參考 R」(空白/100%T)後可切 Raw / 扣暗 / %T / 吸收度。
「➕ 加入光譜清單」存進主視窗。注意:先關 OceanView,一次只能一個程式佔用。

## Vernier Go Direct SpectroVis Plus

```bash
pip install godirect
```

godirect 可連線識別裝置與一般感測通道,但 **Vernier 官方程式庫不支援光譜儀的
光譜串流**(官方文件明載)。可用流程:Spectral Analysis 取光譜 → 匯出 CSV 到
固定資料夾 → 在 Go Direct 頁籤「啟用監看」→ SpectraView 自動載入每個新檔。

## InnoSpectra NIRScan(NIR-S-G1 等)

官方 Python SDK(`ISC_NIRScan_PyQt`)的 `iscpy.pyd` 只能在 **32 位元
Python 3.11(Windows)** 執行,SpectraView 以橋接子行程呼叫它,主程式不受影響。

一次性設定:

1. 解壓 `ISC_NIRScan_PyQt-master.zip`,放在固定位置。
2. 安裝 32-bit Python 3.11:python.org ▸ Downloads ▸ Windows ▸
   「Windows installer (32-bit)」(橋接程式只用標準函式庫,不必再 pip 裝任何東西)。
3. InnoSpectra 頁籤:指定 SDK 資料夾與該 `python.exe` 路徑(會自動記住)。

之後:連線 → 平均次數/燈/參考模式 → 「📷 掃描並加入清單」。「掃描並存為新參考」
會把本次掃描寫回裝置作為之後的參考。輸出可勾吸收度/強度/反射率(各成一條光譜)。
掃描前請關閉 ISC NIRScan GUI。

## MATLAB .mat

- v7 以前:內建支援(具名變數 / N×2 / 矩陣多光譜)。
- **v7.3 (HDF5)**:`pip install h5py` 即可。
- **PLS_Toolbox / Eigenvector dataset 物件**(如 corn.mat):自動展開
  `data` 矩陣並從 `axisscale` 取波長軸,一列一條光譜。

## 測試

```bash
python tests/devices_mat_test.py   # 12 項,免硬體、免 PySide6
```
