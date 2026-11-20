module

public import NNSProofs.GapMoments

open MeasureTheory ProbabilityTheory Finset NNS.GapMoments

namespace NNS.ClaimsA

@[expose] public noncomputable section

variable {n : ℕ}

theorem claim_C025 (P : Measure ℝ) [IsProbabilityMeasure P] (h : MemLp id 3 P) :
    (∫ x, (x - ∫ y, y ∂P) ^ 3 ∂P) = mUPM 3 (∫ y, y ∂P) P - mLPM 3 (∫ y, y ∂P) P :=
  thirdCentralMoment_eq_mUPM_three_sub_mLPM_three h

theorem claim_C106 (P : Measure ℝ) [IsProbabilityMeasure P] (_h : MemLp id 2 P) :
    (∫ x, (x - ∫ y, y ∂P) ^ 2 ∂P) = variance id P :=
  (variance_eq_integral aemeasurable_id').symm

theorem claim_C113 (P : Measure ℝ) [IsProbabilityMeasure P] (h : Integrable id P) :
    let m := ∫ y, y ∂P
    (∫ x, x - m ∂P) = mUPM 1 m P - mLPM 1 m P ∧
      mUPM 1 m P - mLPM 1 m P = 0 ∧ mUPM 1 m P = mLPM 1 m P := by
  dsimp only
  have heq := (mUPM_one_eq_mLPM_one_iff h _).2 rfl
  have hi : Integrable (fun x : ℝ => x) P := h
  simp [heq, integral_sub hi (integrable_const _)]

theorem claim_C116 (x : Fin n → ℝ) (_hn : 0 < n) :
    let v := UPM 2 (mean x) x + LPM 2 (mean x) x
    v = popVariance x ∧ (1 < n → v * (n / ((n : ℝ) - 1)) =
      (∑ i, (x i - mean x) ^ 2) / ((n : ℝ) - 1)) := by
  dsimp only
  rw [add_comm (UPM ..), ← popVariance_eq_PM]
  exact ⟨rfl, fun _ => sampleVariance_eq x⟩

theorem claim_C119 (P : Measure ℝ) [IsProbabilityMeasure P] (h : MemLp id 3 P) :
    (∫ x, (x - ∫ y, y ∂P) ^ 3 ∂P) = mUPM 3 (∫ y, y ∂P) P - mLPM 3 (∫ y, y ∂P) P :=
  claim_C025 P h

theorem claim_C122 :
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ], MemLp (fun x : ℝ => x) 3 μ →
      (0 < mSkewness μ ↔ mLPM 3 (∫ y, y ∂μ) μ < mUPM 3 (∫ y, y ∂μ) μ)) ∧
    ∃ (x : Fin 5 → ℝ) (lo hi : ℝ), (∀ i, lo ≤ x i ∧ x i ≤ hi) ∧
      (∃ i, x i = lo) ∧ (∃ i, x i = hi) ∧
      hi - mean x > mean x - lo ∧ 0 < variance id (empMeasure x) ∧
      mUPM 3 (mean x) (empMeasure x) < mLPM 3 (mean x) (empMeasure x) ∧
      mSkewness (empMeasure x) < 0 := by
  refine ⟨?_, ?_⟩
  · intro μ _ h
    exact mSkewness_pos_iff h
  let x : Fin 5 → ℝ := ![-5, -5, 1, 3, 6]
  have hm : mean x = 0 := by norm_num [mean, x, Fin.sum_univ_five]
  have hv : 0 < popVariance x := by norm_num [popVariance, centralMoment, hm, x, Fin.sum_univ_five]
  have hd : UPM 3 (mean x) x < LPM 3 (mean x) x := by
    norm_num [hm, UPM, LPM, upperDev, lowerDev, x, Fin.sum_univ_five]
  refine ⟨x, -5, 6, ?_, ⟨0, rfl⟩, ⟨4, rfl⟩, by norm_num [hm], ?_, ?_, ?_⟩
  · intro i
    fin_cases i <;> norm_num [x]
  · exact (variance_empMeasure x) ▸ hv
  · simpa only [mUPM_empMeasure, mLPM_empMeasure] using hd
  · rw [mSkewness_empMeasure]
    exact (popSkewness_neg_iff x hv).2 hd

theorem claim_C150 (P : Measure ℝ) [IsProbabilityMeasure P] (h : MemLp id 3 P) :
    (∫ x, (x - ∫ y, y ∂P) ^ 3 ∂P) = mUPM 3 (∫ y, y ∂P) P - mLPM 3 (∫ y, y ∂P) P :=
  claim_C025 P h

theorem claim_C161 (P : Measure ℝ) [IsProbabilityMeasure P] (h : Integrable id P) (t : ℝ) :
    mUPM 1 t P = mLPM 1 t P ↔ t = ∫ y, y ∂P :=
  mUPM_one_eq_mLPM_one_iff h t

theorem claim_C162 (P : Measure ℝ) [IsProbabilityMeasure P] (h : MemLp id 2 P) :
    variance id P = mUPM 2 (∫ y, y ∂P) P + mLPM 2 (∫ y, y ∂P) P :=
  variance_eq_mUPM_two_add_mLPM_two h

theorem claim_C178 (P : Measure ℝ) [IsProbabilityMeasure P] (h : MemLp id 3 P) :
    (∫ x, (x - ∫ y, y ∂P) ^ 3 ∂P) = mUPM 3 (∫ y, y ∂P) P - mLPM 3 (∫ y, y ∂P) P :=
  claim_C025 P h

theorem claim_C183 :
    (∀ (P : Measure ℝ) [IsProbabilityMeasure P], MemLp id 4 P → 0 < variance id P →
      mKurtosis P = (∫ x, (x - ∫ y, y ∂P) ^ 4 ∂P) / (variance id P) ^ 2) ∧
    (∀ {k : ℕ} (x : Fin k → ℝ), 0 < popVariance x →
      mKurtosis (empMeasure x) = (moments x).kurtosis + 3) :=
  ⟨fun _ _ _ _ => rfl, fun x _ => mKurtosis_empMeasure x⟩

theorem claim_C435 (x : Fin n → ℝ) (population : Bool) :
    moments x population = momentsPM x population :=
  moments_eq_momentsPM x population

theorem claim_C460 (x : Fin n → ℝ) (hn : 1 < n) :
    (∑ i, (x i - mean x) ^ 2) / ((n : ℝ) - 1) =
      (UPM 2 (mean x) x + LPM 2 (mean x) x) * (n / ((n : ℝ) - 1)) :=
  ((claim_C116 x (by omega)).2 hn).symm

theorem claim_C463 (x : Fin n → ℝ) (hn : 1 < n) :
    Real.sqrt ((∑ i, (x i - mean x) ^ 2) / ((n : ℝ) - 1)) =
      ((UPM 2 (mean x) x + LPM 2 (mean x) x) * (n / ((n : ℝ) - 1))) ^ (1 / 2 : ℝ) := by
  rw [claim_C460 x hn, Real.sqrt_eq_rpow]

theorem claim_C464 (x : Fin n → ℝ) :
    moments x true = ⟨mean x, popVariance x, popSkewness x, popKurtosis x⟩ ∧
      moments x false = ⟨mean x, sampleVariance x, sampleSkewness x, sampleKurtosis x⟩ :=
  ⟨rfl, rfl⟩

theorem claim_C486 (x : Fin n → ℝ) (hn : 1 < n) :
    (LPM 2 (mean x) x + UPM 2 (mean x) x) * (n / ((n : ℝ) - 1)) =
      (∑ i, (x i - mean x) ^ 2) / ((n : ℝ) - 1) := by
  simpa only [add_comm (UPM ..)] using (claim_C460 x hn).symm

theorem claim_C490 :
    ∃ x y : Fin 7 → ℝ, 0 < popVariance x ∧ 0 < popVariance y ∧ moments x = moments y ∧
      UPM 3 (mean x) x ≠ UPM 3 (mean y) y ∧ LPM 3 (mean x) x ≠ LPM 3 (mean y) y ∧
      UPM 4 (mean x) x ≠ UPM 4 (mean y) y := by
  let x : Fin 7 → ℝ := ![0, 1, 1, 6, 6, 6, 10]
  let y : Fin 7 → ℝ := ![0, 0, 4, 4, 4, 9, 9]
  have hx : mean x = 30 / 7 := by norm_num [mean, x, Fin.sum_univ_succ]
  have hy : mean y = 30 / 7 := by norm_num [mean, y, Fin.sum_univ_succ]
  refine ⟨x, y, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · norm_num [popVariance, centralMoment, hx, x, Fin.sum_univ_succ]
  · norm_num [popVariance, centralMoment, hy, y, Fin.sum_univ_succ]
  · rw [moments_population, moments_population]
    norm_num [centralMoment, hx, hy, x, y, Fin.sum_univ_succ]
  · norm_num [UPM, upperDev, hx, hy, x, y, Fin.sum_univ_succ]
  · norm_num [LPM, lowerDev, hx, hy, x, y, Fin.sum_univ_succ]
  · norm_num [UPM, upperDev, hx, hy, x, y, Fin.sum_univ_succ]

theorem claim_C525 (x : Fin n → ℝ) (_hn : 0 < n) :
    LPM 1 (mean x) x = UPM 1 (mean x) x :=
  (UPM_one_mean_eq_LPM_one_mean x).symm

theorem claim_C526 (x : Fin n → ℝ) (_hn : 0 < n) :
    UPM 3 (mean x) x - LPM 3 (mean x) x = centralMoment 3 (mean x) x ∧
      UPM 4 (mean x) x + LPM 4 (mean x) x = centralMoment 4 (mean x) x :=
  ⟨(centralMoment_three_eq _ x).symm, (centralMoment_four_eq _ x).symm⟩

theorem claim_C585 (x : Fin n → ℝ) (_hn : 0 < n) :
    (moments x).mean = (∑ i, x i) / n ∧ (moments x).variance = centralMoment 2 (mean x) x ∧
    (0 < popVariance x →
      (moments x).skewness =
        centralMoment 3 (mean x) x / (centralMoment 2 (mean x) x) ^ (3 / 2 : ℝ) ∧
      (moments x).kurtosis = centralMoment 4 (mean x) x / (centralMoment 2 (mean x) x) ^ 2 - 3) :=
  ⟨rfl, rfl, fun _ => ⟨rfl, rfl⟩⟩

theorem claim_C586 (x : Fin n → ℝ) (_hn : 0 < n) :
    (moments x false).mean = mean x ∧
    (1 < n → (moments x false).variance = (∑ i, (x i - mean x) ^ 2) / ((n : ℝ) - 1)) ∧
    (2 < n → 0 < popVariance x → (moments x false).skewness =
      Real.sqrt ((n : ℝ) * ((n : ℝ) - 1)) / ((n : ℝ) - 2) * popSkewness x) ∧
    (3 < n → 0 < popVariance x → (moments x false).kurtosis =
      ((n : ℝ) - 1) / (((n : ℝ) - 2) * ((n : ℝ) - 3)) * (((n : ℝ) + 1) * popKurtosis x + 6)) :=
  ⟨rfl, fun _ => sampleVariance_eq x, fun h _ => sampleSkewness_eq_G1 h x,
    fun h _ => sampleKurtosis_eq_G2 h x⟩

end

end NNS.ClaimsA
