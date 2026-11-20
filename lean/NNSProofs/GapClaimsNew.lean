-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.GapSuperiority
public import NNSProofs.GapCertainty
public import NNSProofs.GapDominance

namespace NNS.GapClaimsNew

open MeasureTheory ProbabilityTheory Set Finset
open scoped ENNReal

@[expose] public noncomputable section

section SuperiorityClaims

open GapSuperiority

variable {Ω : Type*} [MeasurableSpace Ω]

lemma probMap {P : Measure Ω} [IsProbabilityMeasure P] {X : Ω → ℝ} (hX : Measurable X) :
    IsProbabilityMeasure (P.map X) :=
  (Measure.isProbabilityMeasure_map_iff hX.aemeasurable).mpr inferInstance

lemma indep_gt_add_gt_add_tie {P : Measure Ω} [IsProbabilityMeasure P] {X Y : Ω → ℝ}
    (hX : Measurable X) (hY : Measurable Y) (hXY : IndepFun X Y P) :
    P.real {ω | Y ω < X ω} = pGt (P.map X) (P.map Y) ∧
      P.real {ω | X ω < Y ω} = pGt (P.map Y) (P.map X) ∧
      P.real {ω | X ω = Y ω} = pTie (P.map X) (P.map Y) :=
  ⟨measureReal_gt_eq_pGt hX hY hXY, measureReal_gt_eq_pGt hY hX hXY.symm,
    measureReal_eq_eq_pTie hX hY hXY⟩

theorem claim_C284 :
    (∀ (μ ν : Measure ℝ) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν],
      NullSingletonClass μ ∨ NullSingletonClass ν → pTie μ ν = 0 ∧ pGt μ ν + pGt ν μ = 1) ∧
    (∀ (P : Measure Ω) [IsProbabilityMeasure P] (X Y : Ω → ℝ), Measurable X → Measurable Y →
      IndepFun X Y P → NullSingletonClass (P.map X) ∨ NullSingletonClass (P.map Y) →
      P.real {ω | X ω = Y ω} = 0 ∧ P.real {ω | Y ω < X ω} + P.real {ω | X ω < Y ω} = 1) := by
  refine ⟨fun μ ν _ _ h => ⟨?_, pGt_add_pGt_of_noAtoms h⟩, fun P _ X Y hX hY hXY h => ?_⟩
  · exact h.elim (fun _ => pTie_eq_zero_left) fun _ => pTie_eq_zero_right
  · have : IsProbabilityMeasure (P.map X) := probMap hX
    have : IsProbabilityMeasure (P.map Y) := probMap hY
    obtain ⟨h1, h2, h3⟩ := indep_gt_add_gt_add_tie hX hY hXY
    rw [h1, h2, h3]
    exact ⟨h.elim (fun _ => pTie_eq_zero_left) fun _ => pTie_eq_zero_right,
      pGt_add_pGt_of_noAtoms h⟩

theorem claim_C285 (μ ν : Measure ℝ) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    pStar μ ν = pGt μ ν + pTie μ ν / 2 ∧ pStar μ ν + pStar ν μ = 1 ∧
      (pStar μ ν = 1 / 2 ↔ pGt μ ν = pGt ν μ) ∧ pStar μ μ = 1 / 2 :=
  ⟨rfl, pStar_add_pStar, pStar_eq_half_iff, pStar_self⟩

theorem claim_C289 {n m : ℕ} (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ssGt x y = (∑ i, ∑ j, if y j < x i then (1 : ℝ) else 0) / ((n : ℝ) * m) ∧
      (0 < n → 0 < m → ssGt x y = pGt (empMeasure x) (empMeasure y)) :=
  ⟨ssGt_eq_sum x y, fun hn hm => ssGt_eq_pGt hn hm x y⟩

theorem claim_C319 (P : Measure Ω) [IsProbabilityMeasure P] {X Y : Ω → ℝ} (hX : Measurable X)
    (hY : Measurable Y) (hXY : IndepFun X Y P) :
    P.real {ω | Y ω < X ω} = pGt (P.map X) (P.map Y) ∧
      P.real {ω | X ω = Y ω} = pTie (P.map X) (P.map Y) ∧
      P.real {ω | Y ω < X ω} + P.real {ω | X ω < Y ω} + P.real {ω | X ω = Y ω} = 1 ∧
      (NullSingletonClass (P.map X) ∨ NullSingletonClass (P.map Y) →
        P.real {ω | X ω = Y ω} = 0 ∧ pGt (P.map X) (P.map Y) + pGt (P.map Y) (P.map X) = 1) := by
  have : IsProbabilityMeasure (P.map X) := probMap hX
  have : IsProbabilityMeasure (P.map Y) := probMap hY
  obtain ⟨h1, h2, h3⟩ := indep_gt_add_gt_add_tie hX hY hXY
  refine ⟨h1, h3, by rw [h1, h2, h3]; exact pGt_add_pGt_add_pTie, fun h => ?_⟩
  exact ⟨h3 ▸ h.elim (fun _ => pTie_eq_zero_left) fun _ => pTie_eq_zero_right,
    pGt_add_pGt_of_noAtoms h⟩

theorem claim_C320 (μ ν : Measure ℝ) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    pGt μ ν + pGt ν μ + pTie μ ν = 1 ∧ pStar μ ν = pGt μ ν + pTie μ ν / 2 ∧
      pStar μ ν + pStar ν μ = 1 ∧ (pStar μ ν = 1 / 2 ↔ pGt μ ν = pGt ν μ) :=
  ⟨pGt_add_pGt_add_pTie, rfl, pStar_add_pStar, pStar_eq_half_iff⟩

theorem claim_C321 (μ ν : Measure ℝ) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    cSS μ ν = |2 * pStar μ ν - 1| ∧ cSS μ ν = |pGt μ ν - pGt ν μ| ∧ pStar μ ν ∈ Icc (0 : ℝ) 1 ∧
      cSS μ ν ∈ Icc (0 : ℝ) 1 ∧ (cSS μ ν = 0 ↔ pGt μ ν = pGt ν μ) ∧
      (cSS μ ν = 1 ↔ pGt μ ν = 1 ∨ pGt ν μ = 1) :=
  ⟨rfl, cSS_eq, pStar_mem_Icc, cSS_mem_Icc, cSS_eq_zero_iff, cSS_eq_one_iff⟩

theorem claim_C323 {n m : ℕ} (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ssGt x y = (∑ i, ∑ j, if y j < x i then (1 : ℝ) else 0) / ((n : ℝ) * m) ∧
      ssStar x y = (∑ i, ∑ j, ((if y j < x i then (1 : ℝ) else 0) +
        1 / 2 * (if x i = y j then (1 : ℝ) else 0))) / ((n : ℝ) * m) ∧
      ssStar x y = ssGt x y + 1 / 2 * ssTie x y ∧
      (0 < n → 0 < m → ssGt x y = pGt (empMeasure x) (empMeasure y) ∧
        ssStar x y = pStar (empMeasure x) (empMeasure y)) :=
  ⟨ssGt_eq_sum x y, ssStar_eq_sum x y, rfl,
    fun hn hm => ⟨ssGt_eq_pGt hn hm x y, ssStar_eq_pStar hn hm x y⟩⟩

theorem claim_C438 :
    (∀ (P : Measure Ω) [IsProbabilityMeasure P] (X Y : Ω → ℝ), Measurable X → Measurable Y →
      IndepFun X Y P → P.real {ω | Y ω < X ω} = pGt (P.map X) (P.map Y)) ∧
    (∀ {n m : ℕ} (x : Fin n → ℝ) (y : Fin m → ℝ),
      ssGt x y = (∑ i, ∑ j, if y j < x i then (1 : ℝ) else 0) / ((n : ℝ) * m) ∧
      ssGt x y = (∑ j, UPM 0 (y j) x) / m ∧
      (0 < n → 0 < m → ssGt x y = pGt (empMeasure x) (empMeasure y))) :=
  ⟨fun _ _ _ _ hX hY hXY => measureReal_gt_eq_pGt hX hY hXY,
    fun x y => ⟨ssGt_eq_sum x y, ssGt_eq_UPM_zero x y, fun hn hm => ssGt_eq_pGt hn hm x y⟩⟩

theorem claim_C500 :
    (∀ (P : Measure Ω) [IsProbabilityMeasure P] (X Y : Ω → ℝ), Measurable X → Measurable Y →
      IndepFun X Y P →
      pGt (P.map X) (P.map Y) = P.real {ω | Y ω < X ω} ∧
      pStar (P.map X) (P.map Y) = P.real {ω | Y ω < X ω} + P.real {ω | X ω = Y ω} / 2) ∧
    (∀ {n m : ℕ} (x : Fin n → ℝ) (y : Fin m → ℝ), ssGt x y = (∑ j, UPM 0 (y j) x) / m) := by
  refine ⟨fun P _ X Y hX hY hXY => ?_, fun x y => ssGt_eq_UPM_zero x y⟩
  rw [measureReal_gt_eq_pGt hX hY hXY, measureReal_eq_eq_pTie hX hY hXY]
  exact ⟨rfl, rfl⟩

theorem claim_C501 (μ ν : Measure ℝ) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    pStar μ ν - 1 / 2 = (pGt μ ν - pGt ν μ) / 2 ∧ (pStar μ ν = 1 / 2 ↔ pGt μ ν = pGt ν μ) ∧
      (1 / 2 < pStar μ ν ↔ pGt ν μ < pGt μ ν) ∧ (pStar μ ν < 1 / 2 ↔ pGt μ ν < pGt ν μ) := by
  have h := pStar_sub_half (μ := μ) (ν := ν)
  refine ⟨h, pStar_eq_half_iff, half_lt_pStar_iff, ?_⟩
  constructor <;> intro h' <;> linarith

theorem claim_C271 :
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ] [NullSingletonClass μ] (c : ℝ) (n : ℕ),
      μ.map (fun v => 2 * c - v) = μ → 2 ≤ n →
      ∫ x, LPM 0 (mean x) x ∂(Measure.pi fun _ : Fin n => μ) = 1 / 2) ∧
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ],
      ∫ x, LPM 0 (mean x) x ∂(Measure.pi fun _ : Fin 1 => μ) = 1) :=
  ⟨fun μ _ _ _ _ hsym hn => integral_ecdf_mean_eq_half μ hsym hn,
    fun μ _ => integral_ecdf_mean_one μ⟩

end SuperiorityClaims

section CertaintyClaims

open GapCertainty

variable {n m : ℕ}

theorem claim_C301 (x : Fin n → ℝ) (y : Fin m → ℝ) :
    pooledMean x y = ((n : ℝ) * mean x + m * mean y) / (n + m) ∧
      pooledMean x y = mean (Fin.append x y) ∧ grandStat x y false = mean (Fin.append x y) :=
  ⟨rfl, pooledMean_eq_mean_append x y, grandStat_eq_mean_append x y⟩

theorem claim_C302 :
    (centralDist ![(3 : ℝ), 3] ![(3 : ℝ), 3] false
        (grandStat ![(3 : ℝ), 3] ![(3 : ℝ), 3] false) = 0 ∧
      min (1 / 2) (max |LPMratio 1 3 ![(3 : ℝ), 3] - 1 / 2| |LPMratio 1 3 ![(3 : ℝ), 3] - 1 / 2|)
        = 1 / 2) ∧
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ) (z : ℝ), 0 < n → 0 < m →
      centralDist x y false z =
          min (1 / 2) (max (if ∀ i, x i = z then 0 else |LPMratio 1 z x - 1 / 2|)
            (if ∀ j, y j = z then 0 else |LPMratio 1 z y - 1 / 2|)) ∧
        centralDist x y false z =
          max (if ∀ i, x i = z then 0 else |LPMratio 1 z x - 1 / 2|)
            (if ∀ j, y j = z then 0 else |LPMratio 1 z y - 1 / 2|)) ∧
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ) (z : ℝ), constVal x ≠ some z → constVal y ≠ some z →
      centralDist x y false z =
        min (1 / 2) (max |LPMratio 1 z x - 1 / 2| |LPMratio 1 z y - 1 / 2|)) :=
  ⟨centralDist_const_witness, fun x y z hn hm => centralDist_eq_piecewise hn hm x y z,
    fun x y _ hx hy => (centralDist_eq_book x y hx hy).1⟩

theorem claim_C303 :
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ),
      (defaultTargets x y).l25 = (LPMVaR 1 (1 / 4) x + LPMVaR 1 (1 / 4) y) / 2 ∧
      (defaultTargets x y).u25 = (UPMVaR 1 (1 / 4) x + UPMVaR 1 (1 / 4) y) / 2 ∧
      (defaultTargets x y).l125 = (LPMVaR 1 (1 / 8) x + LPMVaR 1 (1 / 8) y) / 2 ∧
      (defaultTargets x y).u125 = (UPMVaR 1 (1 / 8) x + UPMVaR 1 (1 / 8) y) / 2) ∧
    (∀ (x : Fin n → ℝ) (p : ℝ), 0 < n → p ∈ Icc (0 : ℝ) 1 →
      UPMVaR 1 p x = LPMVaR 1 (1 - p) x) ∧
    (∀ (x : Fin n → ℝ), (∃ i j, x i ≠ x j) → ∀ t p : ℝ,
      (LPMratio 1 t x = p ↔ (1 - p) * LPM 1 t x = p * UPM 1 t x) ∧
      (UPMratio 1 t x = p ↔ LPMratio 1 t x = 1 - p)) ∧
    (∀ (x : Fin n → ℝ), (∃ i j, x i ≠ x j) → ∀ p ∈ Icc (0 : ℝ) 1,
      |LPMratio 1 (LPMVaR 1 p x) x - p| ≤ rTol ∧ |UPMratio 1 (UPMVaR 1 p x) x - p| ≤ rTol) ∧
    (∀ (cx cy : Option ℝ) (r1 r2 t p : ℝ), 0 ≤ p →
      areaDist cx cy r1 r2 t p ∈ Icc 0 p ∧
      (cx ≠ some t → cy ≠ some t → areaDist cx cy r1 r2 t p = min p (max |r1 - p| |r2 - p|))) :=
  ⟨fun _ _ => ⟨rfl, rfl, rfl, rfl⟩, fun x _ hn hp => UPMVaR_one_eq hn x hp,
    fun x hx t p => ⟨LPMratio_one_eq_iff_balance hx t p,
      UPMratio_eq_iff_LPMratio (by obtain ⟨i, _⟩ := hx; exact Fin.pos i) x t p⟩,
    fun _ hx _ hp => ⟨lowerTarget_ratio hx hp, upperTarget_ratio hx hp⟩,
    fun _ _ _ _ _ _ hp => ⟨areaDist_mem_Icc _ _ _ _ _ hp,
      fun hx hy => areaDist_eq_book hx hy _ _ _⟩⟩

theorem claim_C304 :
    (∀ a b c d e : ℝ, rawScore false a b c d e =
      1 / (25 / 10) * ((1 / 2 - a) ^ 2 / (25 / 100) +
        5 / 10 * ((25 / 100 - b) ^ 2 / (625 / 10000)) +
        5 / 10 * ((25 / 100 - c) ^ 2 / (625 / 10000)) +
        25 / 100 * ((125 / 1000 - d) ^ 2 / (15625 / 1000000)) +
        25 / 100 * ((125 / 1000 - e) ^ 2 / (15625 / 1000000)))) ∧
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ) (medians : Bool) (z : ℝ) (T : Targets),
      certaintyAt x y false medians z T = min 1 (rawScore false (centralDist x y medians z)
        (dU25 x y T) (dL25 x y T) (dU125 x y T) (dL125 x y T) * popAdj n m)) :=
  ⟨rawScore_eq_book, fun _ _ _ _ _ => rfl⟩

theorem claim_C305 :
    (∀ meansOnly : Bool, rawScore meansOnly 0 0 0 0 0 = 1) ∧
      (1 : ℝ) + 1 / 2 + 1 / 2 + 1 / 4 + 1 / 4 = 5 / 2 ∧
      ∀ a b c d e : ℝ, rawScore true a b c d e = ((1 / 2 - a) ^ 2 / (1 / 4)) / 1 :=
  ⟨fun b => (rawScore_zero b).1, (rawScore_zero true).2, rawScore_meansOnly⟩

theorem claim_C306 :
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ] (t : ℝ),
      (1 / 2 ≤ μ.real (Iic t) ∧ 1 / 2 ≤ μ.real (Ici t)) ↔
        μ.real (Iio t) ≤ 1 / 2 ∧ 1 / 2 ≤ cdf μ t) ∧
    (∀ (x : Fin n → ℝ) (t : ℝ), 0 < n →
      (((n : ℝ) / 2 ≤ ∑ i, (if x i ≤ t then (1 : ℝ) else 0) ∧
          (n : ℝ) / 2 ≤ ∑ i, (if t ≤ x i then (1 : ℝ) else 0)) ↔ medianDist t x = 0) ∧
      (medianDist t x = 0 ↔ ltFrac t x ≤ 1 / 2 ∧ 1 / 2 ≤ LPM 0 t x)) := by
  refine ⟨isMedian_iff, fun x t hn => ⟨sample_isMedian_iff hn t x, ?_⟩⟩
  rw [medianDist_eq, medDist_eq_zero_iff]

theorem claim_C307 :
    (∀ (x : Fin n → ℝ) (t : ℝ),
      medianDist t x = max (max 0 (ltFrac t x - 1 / 2)) (1 / 2 - LPM 0 t x) ∧
      ltFrac t x ≤ LPM 0 t x ∧
      (∀ s ∈ Icc (ltFrac t x) (LPM 0 t x), medianDist t x ≤ |s - 1 / 2|) ∧
      (∃ s ∈ Icc (ltFrac t x) (LPM 0 t x), medianDist t x = |s - 1 / 2|) ∧
      (medianDist t x = 0 ↔ ltFrac t x ≤ 1 / 2 ∧ 1 / 2 ≤ LPM 0 t x) ∧
      ((∀ i, x i ≠ t) → medianDist t x = |LPM 0 t x - 1 / 2|)) ∧
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ] (t : ℝ), μ {t} = 0 →
      medDist (μ.real (Iio t)) (cdf μ t) = |cdf μ t - 1 / 2|) ∧
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ) (z : ℝ),
      centralDist x y true z = max (medianDist z x) (medianDist z y) ∧
      controlCDF x true z = LPM 0 z x ∧ controlCDF y true z = LPM 0 z y) := by
  refine ⟨fun x t => ?_, fun μ _ t ht => medDist_eq_abs_of_no_atom μ ht,
    fun x y z => ⟨rfl, controlCDF_medians x z, controlCDF_medians y z⟩⟩
  have hab := ltFrac_le_LPM_zero t x
  rw [medianDist_eq]
  exact ⟨rfl, hab, fun s hs => medDist_le s hs, exists_medDist_eq hab,
    medDist_eq_zero_iff _ _,
    fun ht => by rw [← medianDist_eq]; exact medianDist_eq_abs_of_no_atom ht⟩

theorem claim_C309 :
    (∀ (x : Fin n → ℝ) (c : ℝ), (∀ i, x i = c) →
      LPM 1 c x = 0 ∧ UPM 1 c x = 0 ∧ ∀ p, (1 - p) * LPM 1 c x = p * UPM 1 c x) ∧
    (∀ (t : ℝ) (cy : Option ℝ) (r1 r2 p : ℝ),
      areaDist (some t) cy r1 r2 t p = min p (if cy = some t then 0 else |r2 - p|)) ∧
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ) (c : ℝ), 0 < n → 0 < m → (∀ i, x i = c) →
      (∀ j, y j = c) → ∀ meansOnly medians, certainty x y meansOnly medians = popAdj n m) ∧
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ) (c : ℝ), 0 < n → 0 < m → (∀ i, x i = c) →
      (∃ i j, y i ≠ y j) → mean y = c →
      certainty x y true false = popAdj n m ∧ certainty x y false false = 2 / 5 * popAdj n m) ∧
    (certainty (fun _ : Fin 10 => (0 : ℝ)) (fun j : Fin 10 => (-1 : ℝ) ^ (j : ℕ)) true false
        = 81 / 100 ∧
      certainty (fun _ : Fin 10 => (0 : ℝ)) (fun j : Fin 10 => (-1 : ℝ) ^ (j : ℕ)) false false
        = 81 / 250) ∧
    (mean c01ctrl = c01 ∧ mean c01trt = c01 ∧ certainty c01ctrl c01trt true false = 4 / 9) := by
  refine ⟨fun x c hx => ⟨(LPM_one_const hx).1, (LPM_one_const hx).2, balance_const hx⟩,
    areaDist_const_target,
    fun x y c hn hm hx hy mo me => certainty_const_const hn hm hx hy mo me,
    fun x y c hn hm hx hy hmean => ⟨certainty_const_var_meansOnly hn hm hx hy hmean,
      certainty_const_var_full hn hm hx hy hmean⟩,
    certainty_chapter_example,
    ⟨certainty_c01.1, certainty_c01.2.1, certainty_c01.2.2.2⟩⟩

theorem claim_C310 :
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ) (meansOnly medians : Bool) (z : ℝ) (T : Targets),
      certaintyAt x y meansOnly medians z T = min 1 (rawScore meansOnly (centralDist x y medians z)
        (dU25 x y T) (dL25 x y T) (dU125 x y T) (dL125 x y T) * popAdj n m) ∧
      (0 < n → 0 < m → certaintyAt x y meansOnly medians z T =
        rawScore meansOnly (centralDist x y medians z) (dU25 x y T) (dL25 x y T) (dU125 x y T)
          (dL125 x y T) * popAdj n m)) ∧
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ) (meansOnly medians : Bool),
      certainty x y meansOnly medians =
        certaintyAt x y meansOnly medians (grandStat x y medians) (defaultTargets x y)) ∧
    popAdj n m = (((n : ℝ) + m - 2) / (n + m)) ^ 2 :=
  ⟨fun x y mo me z T => ⟨rfl, fun hn hm => certaintyAt_eq_mul hn hm x y mo me z T⟩,
    fun _ _ _ _ => rfl, rfl⟩

theorem claim_C312 :
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ) (meansOnly medians : Bool) (z : ℝ) (T : Targets),
      rawScore meansOnly (centralDist x y medians z) (dU25 x y T) (dL25 x y T) (dU125 x y T)
          (dL125 x y T) ∈ Icc (0 : ℝ) 1 ∧
      (rawScore meansOnly (centralDist x y medians z) (dU25 x y T) (dL25 x y T) (dU125 x y T)
          (dL125 x y T) = 1 ↔ centralDist x y medians z = 0 ∧ (meansOnly = true ∨
            (dU25 x y T = 0 ∧ dL25 x y T = 0 ∧ dU125 x y T = 0 ∧ dL125 x y T = 0)))) ∧
    (0 < n → 0 < m →
      (∀ (x : Fin n → ℝ) (y : Fin m → ℝ) (meansOnly medians : Bool) (z : ℝ) (T : Targets),
        certaintyAt x y meansOnly medians z T ≤ popAdj n m) ∧ popAdj n m < 1 ∧
      ∀ meansOnly medians : Bool,
        certainty (fun _ : Fin n => (0 : ℝ)) (fun _ : Fin m => (0 : ℝ)) meansOnly medians =
          popAdj n m) := by
  refine ⟨fun x y mo me z T => ⟨rawScore_at_mem x y mo me z T, ?_⟩,
    fun hn hm => ⟨fun x y mo me z T => (certaintyAt_le_popAdj hn hm x y mo me z T).2.1,
      (popAdj_mem hn hm).2,
    fun mo me => certainty_const_const hn hm (fun _ => rfl) (fun _ => rfl) mo me⟩⟩
  obtain ⟨h1, h2, h3, h4⟩ := tails_mem x y T
  have := rawScore_eq_one_iff mo (centralDist_mem_Icc x y me z) h1 h2 h3 h4
  rw [this]

theorem claim_C517 (hn : 0 < n) (hm : 0 < m) :
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ) (meansOnly medians : Bool) (z : ℝ) (T : Targets),
      certaintyAt x y meansOnly medians z T ≤ popAdj n m) ∧ popAdj n m < 1 ∧
    ∀ meansOnly medians : Bool,
      certainty (fun _ : Fin n => (0 : ℝ)) (fun _ : Fin m => (0 : ℝ)) meansOnly medians =
        popAdj n m :=
  ⟨fun x y mo me z T => (certaintyAt_le_popAdj hn hm x y mo me z T).2.1, (popAdj_mem hn hm).2,
    fun mo me => certainty_const_const hn hm (fun _ => rfl) (fun _ => rfl) mo me⟩

theorem claim_C612 :
    (certainty GapCertainty.sepCtrl sepTrt true false = 4 / 9 ∧
      certainty sepCtrl sepTrt true true = 16 / 81 ∧
      certainty sepCtrl sepTrt false false = 8 / 45) ∧
    (∀ (x : Fin n → ℝ) (t : ℝ), medianDist t x = medDist (ltFrac t x) (LPM 0 t x) ∧
      (∃ s ∈ Icc (ltFrac t x) (LPM 0 t x), medianDist t x = |s - 1 / 2|) ∧
      (∀ s ∈ Icc (ltFrac t x) (LPM 0 t x), medianDist t x ≤ |s - 1 / 2|) ∧
      (medianDist t x = 0 ↔ ltFrac t x ≤ 1 / 2 ∧ 1 / 2 ≤ LPM 0 t x)) ∧
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ] (t : ℝ),
      (1 / 2 ≤ μ.real (Iic t) ∧ 1 / 2 ≤ μ.real (Ici t)) ↔
        μ.real (Iio t) ≤ 1 / 2 ∧ 1 / 2 ≤ cdf μ t) ∧
    (∀ (x : Fin n → ℝ) (c : ℝ), (∀ i, x i = c) → ∀ p, (1 - p) * LPM 1 c x = p * UPM 1 c x) ∧
    (∀ (t : ℝ) (cy : Option ℝ) (r1 r2 p : ℝ),
      areaDist (some t) cy r1 r2 t p = min p (if cy = some t then 0 else |r2 - p|)) := by
  refine ⟨certainty_modes_distinct, fun x t => ?_, isMedian_iff, fun x c hx => balance_const hx,
    areaDist_const_target⟩
  have hab := ltFrac_le_LPM_zero t x
  rw [medianDist_eq]
  exact ⟨rfl, exists_medDist_eq hab, fun s hs => medDist_le s hs, medDist_eq_zero_iff _ _⟩

theorem claim_C621 :
    (∀ (x : Fin n → ℝ) (t : ℝ), medianDist t x = medDist (ltFrac t x) (LPM 0 t x) ∧
      (∃ s ∈ Icc (ltFrac t x) (LPM 0 t x), medianDist t x = |s - 1 / 2|) ∧
      (medianDist t x = 0 ↔ ltFrac t x ≤ 1 / 2 ∧ 1 / 2 ≤ LPM 0 t x)) ∧
    (∀ (x : Fin n → ℝ), 0 < n → centralDist x x true (grandStat x x true) = 0) ∧
    (∀ (t : ℝ) (cy : Option ℝ) (r1 r2 p : ℝ),
      areaDist (some t) cy r1 r2 t p = min p (if cy = some t then 0 else |r2 - p|)) ∧
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ) (c : ℝ), 0 < n → 0 < m → (∀ i, x i = c) →
      (∀ j, y j = c) → ∀ meansOnly medians, certainty x y meansOnly medians = popAdj n m) := by
  refine ⟨fun x t => ?_, fun x hn => centralDist_median_self hn x,
    areaDist_const_target,
    fun x y c hn hm hx hy mo me => certainty_const_const hn hm hx hy mo me⟩
  have hab := ltFrac_le_LPM_zero t x
  rw [medianDist_eq]
  exact ⟨rfl, exists_medDist_eq hab, medDist_eq_zero_iff _ _⟩

end CertaintyClaims

section DominanceClaims

open GapDominance

variable {n m : ℕ}

theorem claim_C288 (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ksStat x y = sSup (Set.range fun t => |LPM 0 t x - LPM 0 t y|) ∧
      (∀ t, LPM 0 t x = (#{i | x i ≤ t} : ℝ) / n) ∧
      (∀ t, |LPM 0 t x - LPM 0 t y| ≤ ksStat x y) ∧
      (0 < n → ∃ u, IsKnot x y u ∧ ksStat x y = |LPM 0 u x - LPM 0 u y|) :=
  ⟨rfl, fun t => LPM_zero_eq_ecdf t x, le_ksStat x y, fun hn => ksStat_attained hn x y⟩

theorem claim_C328 :
    (∀ μ ν : Measure ℝ, LawFSD μ ν ↔ (∀ t, cdf μ t ≤ cdf ν t) ∧ ∃ t, cdf μ t < cdf ν t) ∧
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ),
      LPMDominates 0 x y ↔ (∀ t, LPM 0 t x ≤ LPM 0 t y) ∧ ∃ t, LPM 0 t x < LPM 0 t y) ∧
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ), 0 < n → 0 < m →
      (sdDomPairCentred x y 1 true = 1 ↔
        (∀ t, (#{i | x i ≤ t} : ℝ) / n ≤ (#{j | y j ≤ t} : ℝ) / m) ∧
          ∃ t, (#{i | x i ≤ t} : ℝ) / n < (#{j | y j ≤ t} : ℝ) / m)) := by
  refine ⟨fun μ ν => ⟨fun ⟨h, t, ht⟩ => ⟨h, t, lt_of_le_of_ne (h t) ht⟩,
    fun ⟨h, t, ht⟩ => ⟨h, t, ht.ne⟩⟩, fsd_iff_lt, fun x y hn hm => fsdUni_eq_one_iff hn hm x y⟩

theorem claim_C343 :
    (∀ μ : Measure ℝ, ¬ LawFSD μ μ) ∧
    (∀ μ ν κ : Measure ℝ, LawFSD μ ν → LawFSD ν κ → LawFSD μ κ) ∧
    (∀ {k : ℕ} (M : Fin k → Measure ℝ) (i : Fin k),
      (i ∈ effLaws M ↔ ∀ j, ¬ LawFSD (M j) (M i)) ∧
      (i ∉ effLaws M → ∃ e ∈ effLaws M, LawFSD (M e) (M i))) ∧
    (∀ (μ ν : Measure ℝ) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν],
      LawFSD μ ν → ∀ u : ℝ → ℝ, Monotone u → Integrable u μ → Integrable u ν →
      ∫ x, u x ∂ν ≤ ∫ x, u x ∂μ) ∧
    (∀ {k r : ℕ} (X : Fin k → Fin r → ℝ), 0 < r → ∀ (R : Finset (Fin k)) (v : Fin k),
      (v ∈ effSet (domMatrix X 1 true) R ↔
        v ∈ R ∧ ∀ j ∈ R, sdDomPairCentred (X j) (X v) 1 true ≠ 1) ∧
      (v ∈ R → v ∉ effSet (domMatrix X 1 true) R →
        ∃ e ∈ effSet (domMatrix X 1 true) R, sdDomPairCentred (X e) (X v) 1 true = 1 ∧
          ∀ u : ℝ → ℝ, Monotone u → mean (fun i => u (X v i)) ≤ mean (fun i => u (X e i)))) := by
  refine ⟨fun μ => curveDom_irrefl _, fun _ _ _ h1 h2 => curveDom_trans h1 h2,
    fun M i => ⟨by simp [effLaws], exists_effLaws_dominates M⟩,
    fun μ ν _ _ h u hu hiμ hiν => integral_le_of_cdf_le h.1 hu hiμ hiν,
    fun X hr R v => ⟨effSet_eq hr X (Or.inl rfl) true R v, fun hv hnot => ?_⟩⟩
  obtain ⟨e, he, hev⟩ := exists_eff_dominates hr X (Or.inl rfl) true hv hnot
  exact ⟨e, he, hev, fun u hu => mean_comp_le_of_fsdUni hr hr hev hu⟩

theorem claim_C344 :
    (∀ {k r : ℕ} (X : Fin k → Fin r → ℝ), 0 < r → ∀ d, (d = 1 ∨ d = 2 ∨ d = 3) → ∀ b,
      let D := domMatrix X d b
      Acyclic D ∧
      (∀ (R : Finset (Fin k)) (i j : Fin k), i ∈ effSet D R → j ∈ effSet D R → D i j ≠ 1) ∧
      (∀ minC : ℕ, 0 < minC → ∀ S ∈ (sdCluster D minC).dropLast, ∀ i ∈ S, ∀ j ∈ S, D i j ≠ 1)) ∧
    (sdCluster (domMatrix chainFamily 1 true) 2 = [{0}, {1, 2}] ∧
      sdDomPairCentred (chainFamily 1) (chainFamily 2) 1 true = 1) ∧
    ((0 : Fin 3) ∈ peelLevel (domMatrix peelFamily 1 true) 0 ∧
      (2 : Fin 3) ∈ peelLevel (domMatrix peelFamily 1 true) 0 ∧
      (1 : Fin 3) ∈ peelLevel (domMatrix peelFamily 1 true) 1 ∧
      sdDomPairCentred (peelFamily 0) (peelFamily 1) 1 true = 1 ∧
      sdDomPairCentred (peelFamily 2) (peelFamily 1) 1 true ≠ 1) := by
  refine ⟨?_, cluster_witness, peel_witness⟩
  intro k r X hr d hd b
  have hD := acyclic_domMatrix hr X hd b
  exact ⟨hD, fun R _ _ hi hj => effSet_antichain_of_acyclic hD R hi hj,
    fun minC _ _ hS _ hi _ hj => sdCluster_antichain_of_acyclic hD minC hS hi hj⟩

theorem claim_C519 :
    (∀ σ τ : Equiv.Perm (Fin 1000),
      nnsFSD (col c519xN ∘ σ) (col c519yN ∘ τ) true = .yDom ∧
        nnsFSD (col c519yN ∘ τ) (col c519xN ∘ σ) true = .xDom) ∧
    (∀ (x : Fin n → ℝ) (y : Fin m → ℝ), 0 < n → 0 < m →
      (nnsFSD x y true = .yDom ↔ LPMDominates 0 y x)) := by
  refine ⟨c519_verdicts, fun x y hn hm => ?_⟩
  rw [nnsFSD_eq_yDom_iff hn hm]
  simp only [domCurve, ↓reduceIte]
  rfl

theorem claim_C520 :
    ((sdDomPairCentred (c520Cols 1) (c520Cols 0) 1 true = 1 ∧
      sdDomPairCentred (c520Cols 3) (c520Cols 2) 1 true = 1 ∧
      sdDomPairCentred (c520Cols 5) (c520Cols 4) 1 true = 1 ∧
      sdDomPairCentred (c520Cols 7) (c520Cols 6) 1 true = 1) ∧
    (∀ i ∈ ({1, 3, 5, 7} : Finset (Fin 8)), ∀ j ∈ ({1, 3, 5, 7} : Finset (Fin 8)),
      sdDomPairCentred (c520Cols i) (c520Cols j) 1 true ≠ 1) ∧
    effSet (domMatrix c520Cols 1 true) univ = {1, 3, 5, 7}) ∧
    (∀ (x : Fin n → ℝ) (c : ℝ), 0 < n → 0 < c →
      sdDomPairCentred (fun i => x i + c) x 1 true = 1) :=
  ⟨c520_effSet, fun x _ hn hc => fsdUni_shift hn x hc⟩

theorem claim_C591 :
    sdDomPairCentred c591X c591Y 1 false = 1 ∧ sdDomPairCentred c591X c591Y 1 true = 0 ∧
      Monotone (fun v : ℝ => if 3 < v then (1 : ℝ) else 0) ∧
      mean (fun i => if 3 < c591X i then (1 : ℝ) else 0) <
        mean (fun j => if 3 < c591Y j then (1 : ℝ) else 0) :=
  ⟨c591_omega, c591_not_fsd, c591_utility.1, c591_utility.2⟩

end DominanceClaims

end

end NNS.GapClaimsNew
