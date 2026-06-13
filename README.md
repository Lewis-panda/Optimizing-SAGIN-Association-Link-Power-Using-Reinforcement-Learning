# SAGIN 關聯 + 功率最佳化 — Q-Learning（乾淨重寫版）

用表格式 Q-Learning，在 **太空–空中–地面整合網路 (SAGIN)** 中同時最佳化
**關聯 (association)** 與 **傳輸功率 (power)**，最大化系統總速率 (sum-rate)。

這是三年前一份 MATLAB 專案
([原始 repo](https://github.com/Lewis-panda/Optimizing-SAGIN-Association-Link-Power-Using-Reinforcement-Learning))
的乾淨、可重現、有 baseline 對照的重寫版本。

> ### ⚠️ 定位聲明（請先讀）
> 這個題目（satellite–HAPS–ground 的 joint 關聯 + 功率，用 RL）**已有大量 prior work**，
> 本專案**不主張任何原創研究貢獻**。它的定位是「正確的工程參考實作 / 作品集 / 別人論文的 baseline」。
> 代表性 prior work：
> - Alsharoa & Alouini, *Joint User Association and Beamforming in Integrated Satellite-HAPS-Ground Networks*, IEEE TWC, [arXiv:2204.13257](https://arxiv.org/abs/2204.13257)
> - *Deep Q-Learning-Based Transmission Power Control of a HAPS with Spectrum Sharing*, MDPI Sensors 2022, [link](https://www.mdpi.com/1424-8220/22/4/1630)
> - *Machine Learning-Based User Scheduling in Integrated Satellite-HAPS-Ground Networks*, [arXiv:2205.13958](https://arxiv.org/abs/2205.13958)
> - 綜述：*On the Interplay of AI and SAGIN: A Survey*, [arXiv:2402.00881](https://arxiv.org/abs/2402.00881)

---

## 問題定義

兩層獨立的下行鏈路（各自佔用一段頻譜）：

```
tier 1 (backhaul):  LEO  ──►  HAPS
tier 2 (access):    HAPS ──►  Ground User
```

每個接收端只關聯一個發射端（association）；每個發射端從離散功率碼本中選一個功率等級。
目標是最大化兩層的總速率 `tier1 + tier2`（bits/s/Hz）。

![SAGIN architecture](docs/architecture.png)

## 系統模型

| 元件 | 模型 |
|------|------|
| 幾何 | 100×100 km 服務區，LEO@300km、HAPS@20km、GU@0km，位置隨機 |
| 路徑損耗 | 自由空間 `FSPL(dB)=20log10(d)+20log10(f)+20log10(4π/c)`，距離以 **公尺** 計 |
| 通道增益 | **線性** `g = 10^(gain_dB/10) · |h|²`，`|h|²~Exp(1)`（Rayleigh 衰落）。距離越遠增益越小 |
| 干擾 | 全頻率重用：接收端 j 受 **所有其他在傳的發射端** 同頻干擾 |
| SINR | `SINR_j = P[a_j]·g[a_j,j] / ( Σ_{m≠a_j, active} P[m]·g[m,j] + N0 )` |
| 雜訊 | 熱雜訊 `N0 = kTB · 10^(NF/10)` |
| 小區內資源共享 | 每個發射端把頻寬平均分給它服務的接收端：`rate_j = (1/load)·log2(1+SINR_j)` |

**為什麼這個模型才「合理」**：因為干擾隨功率上升、且小區內要分享資源，所以
「大家都開最大功率」與「全部塞給同一個最佳發射端」都不是最佳解 —— 功率控制與負載平衡都變成
真正需要權衡的問題。原版缺了這兩點，導致最佳解 trivial（見下表）。

## 方法

- **表格式 Q-Learning**，每個決策者是獨立 learner（independent multi-agent Q-learning）。
- 狀態 = 自己上一次的動作；動作 = 新的關聯 / 新的功率等級（Bellman 更新）。
- **信用分配 (credit assignment)**：
  - 關聯 agent 用**自己鏈路的速率**（local reward）—— 接收端能直接判斷自己選的發射端好不好；
  - 功率 agent 用**整層的總速率**（global reward）—— 才能感受到自己加大功率對別人造成的干擾外部性。
- ε-greedy，ε 線性衰減；小尺度衰落每個 episode 重抽，學的是 **ergodic（期望）速率**。

![Q-learning methodology](docs/methodology.png)

## 相對原始 MATLAB 版修正了什麼

| 原版問題 | 後果 | 本版修正 |
|----------|------|----------|
| 把路徑損耗 (dB) 直接當通道增益 | 距離越遠訊號越強，物理顛倒 | 增益改線性 `10^(-PL/10)`，並補天線增益做合理 link budget |
| `20log10(4π/c)+147.55` 常數自我抵消、距離用 km 卻套用 m 的公式 | FSPL 數值錯 ~60–147 dB | 常數正確、單位統一在 `channel.py` 換算 |
| 干擾 = 增益和 − 收訊功率（不乘功率，可為負） | 功率越大 SINR 單調變大 → 最佳功率恆為 P_max，功率最佳化是假問題 | 干擾改成**其他發射端的** `P·g` 之和，產生真正的功率↔干擾權衡 |
| 雜訊 = 單一複數樣本的平方 | 物理無意義 | 改成熱雜訊功率 `kTB·NF` |
| `Power_Leos` vs `Power_Leo` 變數名打錯 | LEO 功率動作從未進入 reward（半套最佳化是 no-op） | 重寫消除此類接線錯誤 |
| 隨機探索 link 時沒建連線矩陣 | 探索其實沒套用到 reward | 動作直接決定關聯，無此 bug |
| 狀態 = 量化後的 reward | 退化 MDP（state≈reward） | 狀態 = 上一動作；reward 做 local/global 分流 |
| 無 baseline、單一拓撲、無多次平均 | 結論無法驗證 | 三種 baseline + 跨拓撲 Monte-Carlo + ±std |
| 無小區內資源共享 | 可把所有人塞到一個無干擾發射端，速率虛高 | 加入等量頻寬分享，association 變成真實負載平衡問題 |

## 專案結構

```
config.py      所有參數（dataclass，含 .quick() 快速版）
channel.py     幾何、FSPL、衰落、link budget（純物理，可單測）
env.py         sum_rate：SINR + 干擾 + 小區資源共享
agents.py      AgentGroup：向量化表格式 Q-learner + ε 衰減
trainer.py     單一拓撲的訓練與貪婪策略評估
baselines.py   Random / MaxPower+Greedy / BestUniform+Greedy
experiment.py  跨拓撲 Monte-Carlo 主程式，輸出圖表/CSV
tests.py       物理模型的健全性測試（守住已修的 bug）
docs/          架構 / 方法示意圖（.drawio 原始檔 + 匯出 .png，可用 draw.io 編輯）
```

## 如何執行

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python tests.py            # 5 個健全性測試
python experiment.py --quick   # 快速 smoke run（數秒）
python experiment.py           # 完整實驗（~15 秒，輸出到 results/）
```

常用旗標：`--seed`, `--episodes`, `--topologies`, `--eval-samples`。

## 結果

完整設定（5 LEO / 10 HAPS / 15 GU，20000 episodes，跨 10 個隨機拓撲）：

| 方法 | 總速率 (bits/s/Hz) |
|------|--------------------|
| Random | 5.10 ± 0.15 |
| MaxPower + Greedy | 7.08 ± 0.43 |
| BestUniform + Greedy | 7.12 ± 0.41 |
| **Q-Learning** | **13.93 ± 0.81** |

Q-Learning 在 **10 個拓撲全部勝出**，較最佳 baseline **+95.6%**。

![convergence](results/learning_curve.png)
![comparison](results/comparison.png)

**怎麼讀**：`MaxPower < BestUniform` 證明「最大功率不是最佳」——干擾權衡是真的。
Q-Learning 的主要增益來自 **tier-1 的負載平衡**（貪婪關聯把所有 HAPS 塞給單一最佳 LEO，
反而比隨機還差），在 tier-2 則與貪婪相當。

## 已知簡化與可延伸方向

- **獨立 multi-agent Q-learning**：收斂無理論保證；功率 agent 用全域 reward 仍有 credit-assignment 雜訊。可換 difference reward / VDN / QMIX。
- **兩層獨立、未做端到端耦合**：實務上每個 GU 的端到端速率應為 `min(backhaul, access)`。目前分開最佳化兩層的 sum-rate。
- **小區內等量共享、靜態拓撲、無 LEO 移動性 / handover**：可加入軌道動力學、duty cycle、QoS 約束。
- **方法升級**：狀態改連續特徵後可換 **DQN / SAC / multi-agent DRL**，對應上面 prior work 的現代做法。`agents.py` 的介面已預留替換空間。

## 授權與出處

教學 / 作品集用途。物理與評估方法見上表與 prior work 連結。
