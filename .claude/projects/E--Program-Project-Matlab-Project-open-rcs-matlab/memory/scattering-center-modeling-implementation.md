---
name: scattering-center-modeling-implementation
description: Implementation of PO-based scattering center modeling following the 6-step method document
metadata:
  type: project
---

# Scattering Center Modeling Implementation

Created `main_scattering_center_modeling.m` implementing the complete 6-step pipeline described in `基于PO结果的散射中心建模方法.md`:

1. **HRRP生成**: IFFT along frequency dimension with windowing and zero-padding
2. **CLEAN峰值检测**: Iterative frequency-domain subtraction, finding strongest HRRP peaks
3. **参数估计**: α via phase-demodulated log-log fit; A from peak magnitude; angle tracking for L and dR/dφ
4. **分类**: Decision tree → local (narrow angle, fixed R), distributed (broadened range/angle), sliding (R moves with φ)
5. **参数优化**: Levenberg-Marquardt nonlinear least squares (lsqnonlin) with bounds
6. **模型验证**: Correlation coefficient ρ, HRRP comparison, RCS comparison, frequency-domain phase/magnitude comparison

GTD model types:
- Local: S = A·(jf/fc)^α·exp(-j4πfR/c)
- Distributed: Same × sinc(2πfL·sin(φ-φ̄)/c)
- Sliding: Same with R(φ) = R₀ + dR/dφ·(φ-φ̄)

Input: `results/wideband_scattering_*.mat` (from main_wideband_scattering)
Usage: `main_scattering_center_modeling` or `main_scattering_center_modeling('file.mat', 'MaxCenters', 15)`

**Why:** User asked to implement the program following the method document after studying existing code.
**How to apply:** Run `main_wideband_scattering` first to generate wideband data, then run `main_scattering_center_modeling`.
