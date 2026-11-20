-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.GapCoMoments
public import NNSProofs.MeasurePM
public import NNSProofs.Moments
public import NNSProofs.Ratios

namespace NNS.ClaimsA

open Finset Matrix MeasureTheory ProbabilityTheory GapCoMoments

@[expose] public noncomputable section

variable {n k m : ℕ}

theorem claim_C024 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : MemLp (fun v : ℝ => v) 2 μ)
    (x : Fin n → ℝ) (_hn : 1 < n) :
    variance id μ = mUPM 2 (∫ v, v ∂μ) μ + mLPM 2 (∫ v, v ∂μ) μ ∧
    sampleVariance x = (n / ((n : ℝ) - 1)) * (UPM 2 (mean x) x + LPM 2 (mean x) x) := by
  refine ⟨variance_eq_mUPM_two_add_mLPM_two h, ?_⟩
  rw [sampleVariance, popVariance_eq_PM]
  ring

theorem claim_C163 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : MemLp (fun v : ℝ => v) 2 μ)
    (x : Fin n → ℝ) (hn : 1 < n) :
    variance id μ = mUPM 2 (∫ v, v ∂μ) μ + mLPM 2 (∫ v, v ∂μ) μ ∧
    sampleVariance x = (n / ((n : ℝ) - 1)) * (UPM 2 (mean x) x + LPM 2 (mean x) x) :=
  claim_C024 μ h x hn

theorem claim_C391 {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω) [IsProbabilityMeasure P]
    (X : Fin k → Ω → ℝ) (h : ∀ i, MemLp (X i) 2 P) (i j : Fin k) :
    let μ := P.map fun ω => (X i ω, X j ω)
    let a := ∫ ω, X i ω ∂P
    let b := ∫ ω, X j ω ∂P
    cov[X i, X j; P] = mCoLPM 1 1 a b μ + mCoUPM 1 1 a b μ - mDLPM 1 1 a b μ
      - mDUPM 1 1 a b μ := by
  let f := fun ω => (X i ω, X j ω)
  have hf : AEMeasurable f P := (h i).aemeasurable.prodMk (h j).aemeasurable
  have hi := (memLp_map_measure_iff (g := Prod.fst) (by fun_prop) hf).2 (h i)
  have hj := (memLp_map_measure_iff (g := Prod.snd) (by fun_prop) hf).2 (h j)
  have hh := covariance_decomposition_measure
    (integrable_coMoment_of_memLp hi hj (∫ q, q.1 ∂P.map f) (∫ q, q.2 ∂P.map f))
  rw [covariance_map_fun (by fun_prop) (by fun_prop) hf,
    integral_map hf (by fun_prop), integral_map hf (by fun_prop)] at hh
  dsimp only
  exact hh.symm.trans (by ac_rfl)

theorem claim_C392 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ) :
    let t := ∫ x, x ∂μ
    cov[id, id; μ] = variance id μ ∧
    (∫ x, upperDev 1 t x * lowerDev 1 t x ∂μ) = 0 ∧
    (∫ x, lowerDev 1 t x * upperDev 1 t x ∂μ) = 0 ∧
    variance id μ = mUPM 2 t μ + mLPM 2 t μ ∧
    (∀ (l u a v : ℝ), rLowerDev l a v * rUpperDev u a v = 0) := by
  intro t
  have hz (x : ℝ) : upperDev 1 t x * lowerDev 1 t x = 0 := by
    rcases le_or_gt x t with hh | hh
    · simp [upperDev, hh.not_gt]
    · simp [lowerDev, hh.not_ge]
  exact ⟨covariance_self measurable_id.aemeasurable, by simp [hz],
    by simp [mul_comm (lowerDev 1 t _), hz], variance_eq_mUPM_two_add_mLPM_two h, fun l u a v => by
      rcases le_or_gt v a with hv | hv
      · simp [rUpperDev, hv.not_gt]
      · simp [rLowerDev, hv.not_ge]⟩

theorem claim_C394 (hn : 1 < n) (X : Matrix (Fin n) (Fin k) ℝ) :
    let Q := PMmatrix 1 1 (colMeans X) X true false
    Q.clpm + Q.cupm - Q.dlpm - Q.dupm = sampleCovMatrix X ∧ Q.cov = sampleCovMatrix X := by
  intro Q
  refine ⟨?_, PMmatrix_cov_eq_sampleCov hn X⟩
  rw [← PMmatrix_cov_eq_sampleCov hn, PMmatrix_cov_eq]
  dsimp only [Q]
  abel

namespace A9_CoMoments

def devMatrix (g : ℝ → ℝ → ℝ → ℝ) (d : ℝ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) :
    Matrix (Fin n) (Fin k) ℝ := Matrix.of fun r j => g d (t j) (X r j)

def coMatrix (f : ℝ → ℝ → ℝ → ℝ → (Fin n → ℝ) → (Fin n → ℝ) → ℝ) (d : ℝ)
    (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) : Matrix (Fin k) (Fin k) ℝ :=
  Matrix.of fun i j => f d d (t i) (t j) (colOf X i) (colOf X j)

end A9_CoMoments

open A9_CoMoments

theorem claim_C396 (d : ℝ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) :
    let L := devMatrix rLowerDev d t X
    let U := devMatrix rUpperDev d t X
    coMatrix rCoLPM d t X = (1 / (n : ℝ)) • (Lᵀ * L) ∧
    coMatrix rCoUPM d t X = (1 / (n : ℝ)) • (Uᵀ * U) ∧
    (∀ (v : ℕ) (adj : Bool),
      let factor := if adj = true ∧ 1 < n ∧ 0 < v then 1 / ((n : ℝ) - 1) else 1 / (n : ℝ)
      (PMmatrix v v t X adj false).clpm =
        factor • ((lowerMat v t X)ᵀ * lowerMat v t X) ∧
      (PMmatrix v v t X adj false).cupm =
        factor • ((upperMat v t X)ᵀ * upperMat v t X)) ∧
    (let Z : Matrix (Fin 2) (Fin 1) ℝ := Matrix.of fun i _ => if i = 0 then 0 else 2
     (PMmatrix 1 1 (colMeans Z) Z false false).clpm 0 0 = 1 / 2 ∧
     (PMmatrix 1 1 (colMeans Z) Z true false).clpm 0 0 = 1) := by
  dsimp only
  refine ⟨?_, ?_, ?_, ?_⟩
  · ext i j
    simp [devMatrix, coMatrix, rCoLPM, Matrix.mul_apply, colOf, div_eq_inv_mul]
  · ext i j
    simp [devMatrix, coMatrix, rCoUPM, Matrix.mul_apply, colOf, div_eq_inv_mul]
  · intro v adj
    have hs : pmScale adj n v v / n =
        if adj = true ∧ 1 < n ∧ 0 < v then 1 / ((n : ℝ) - 1) else 1 / (n : ℝ) := by
      simp only [pmScale, and_self]
      split_ifs with h
      · have hn : (n : ℝ) ≠ 0 := by exact_mod_cast (by omega : n ≠ 0)
        field_simp
      · rfl
    rw [← hs]
    exact ⟨clpmM_eq_gram v v t X adj, cupmM_eq_gram v v t X adj⟩
  · norm_num [PMmatrix, clpmM, CoLPM, colOf, colMeans, mean, pmScale, lowerDev, Fin.sum_univ_two]

theorem claim_C397 (d : ℝ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (w : Fin k → ℝ) :
    let L := coMatrix rCoLPM d t X
    let U := coMatrix rCoUPM d t X
    L.PosSemidef ∧ U.PosSemidef ∧
    w ⬝ᵥ (L *ᵥ w) = (1 / (n : ℝ)) * ∑ r, (∑ j, w j * rLowerDev d (t j) (X r j)) ^ 2 ∧
    w ⬝ᵥ (U *ᵥ w) = (1 / (n : ℝ)) * ∑ r, (∑ j, w j * rUpperDev d (t j) (X r j)) ^ 2 := by
  dsimp only
  obtain ⟨hl, hu, _⟩ := claim_C396 d t X
  rw [hl, hu]
  exact ⟨posSemidef_gram (by positivity) _, posSemidef_gram (by positivity) _,
    quadForm_gram .., quadForm_gram ..⟩

theorem claim_C433 (_hn : 0 < n) (X : Matrix (Fin n) (Fin k) ℝ) :
    let P := PMmatrix 1 1 (colMeans X) X false false
    let S := PMmatrix 1 1 (colMeans X) X true false
    (P.clpm + P.cupm - P.dlpm - P.dupm = popCovMatrix X ∧ P.cov = popCovMatrix X) ∧
    (1 < n → S.clpm + S.cupm - S.dlpm - S.dupm = sampleCovMatrix X ∧
      S.cov = sampleCovMatrix X) := by
  intro P S
  refine ⟨⟨?_, PMmatrix_cov_eq_popCov X⟩, fun hn => claim_C394 hn X⟩
  rw [← PMmatrix_cov_eq_popCov, PMmatrix_cov_eq]
  dsimp only [P]
  abel

theorem claim_C467 (hn : 1 < n) (X : Matrix (Fin n) (Fin k) ℝ) :
    let Q := PMmatrix 1 1 (colMeans X) X true false
    sampleCovMatrix X = Q.clpm + Q.cupm - Q.dlpm - Q.dupm :=
  (claim_C394 hn X).1.symm

theorem claim_C471 (x y : Fin n → ℝ) (a b : Fin m → ℝ) :
    (fun j => CoLPM 0 0 (a j) (b j) x y) =
      fun j => (#{i | x i ≤ a j ∧ y i ≤ b j} : ℝ) / n :=
  funext fun j => CoLPM_zero_eq_card (a j) (b j) x y

theorem claim_C511 (d : ℝ) (e : Option ℝ) (a b : ℝ) (x y : Fin n → ℝ) :
    rCoLPM d (e.getD d) a b x y =
      (∑ i, lowKernel d a (x i) * lowKernel (e.getD d) b (y i)) / n ∧
    rCoUPM d (e.getD d) a b x y =
      (∑ i, upKernel d a (x i) * upKernel (e.getD d) b (y i)) / n := by
  obtain ⟨hl, hu, _⟩ := rCo_eq_kernels d (e.getD d) d d a b x y
  exact ⟨hl, hu⟩

private lemma strictIndicator (x t : Fin k → ℝ) :
    (if ∀ j, t j < x j then (1 : ℝ) else 0) =
      ∑ A ∈ (univ : Finset (Fin k)).powerset,
        (-1 : ℝ) ^ A.card * (if ∀ j ∈ A, x j ≤ t j then 1 else 0) := by
  have h (j : Fin k) : (1 : ℝ) - (if x j ≤ t j then 1 else 0) =
      if t j < x j then 1 else 0 := by
    by_cases hj : x j ≤ t j <;> simp [hj, not_lt.mpr, lt_of_not_ge]
  simpa [h, Finset.prod_boole] using
    Finset.prod_sub (fun _ : Fin k => (1 : ℝ)) (fun j => if x j ≤ t j then 1 else 0) univ

open scoped Classical in
theorem claim_C544 (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) :
    coUPMnD_R X t 0 true = (#{i | ∀ j, t j < X i j} : ℝ) / n ∧
    coUPMnD_R X t 0 true = ∑ A ∈ (univ : Finset (Fin k)).powerset,
      (-1 : ℝ) ^ A.card * ((#{i | ∀ j ∈ A, X i j ≤ t j} : ℝ) / n) := by
  refine ⟨coUPMnD_R_zero .., ?_⟩
  rw [coUPMnD_R_zero, ← Finset.sum_boole]
  simp_rw [strictIndicator]
  rw [Finset.sum_comm, Finset.sum_div]
  simp_rw [← Finset.mul_sum, Finset.sum_boole, mul_div_assoc]

theorem claim_C545 {d : ℝ} (hd : 0 < d) (X : Matrix (Fin n) (Fin k) ℝ)
    (T : Matrix (Fin m) (Fin k) ℝ) (r : Fin m) :
    coLPMnDBatch (-X) (-T) d true true r = coUPMnD_R X (T r) d true :=
  coLPMnDBatch_neg_eq_coUPMnD_R hd.le X T true r

namespace A9_CoMoments

def adjustment (adj : Bool) (n : ℕ) (l u : ℝ) : ℝ :=
  if adj = true ∧ 1 < n ∧ 0 < l ∧ 0 < u then (n : ℝ) / (n - 1) else 1

def quadrants (l u a b : ℝ) (x y : Fin n → ℝ) : Fin 4 → ℝ :=
  ![rCoLPM l l a b x y, rCoUPM u u a b x y, rDLPM l u a b x y, rDUPM l u a b x y]

def adjustedQuadrants (l u a b : ℝ) (x y : Fin n → ℝ) (adj : Bool) : Fin 4 → ℝ :=
  adjustment adj n l u • quadrants l u a b x y

def normalizedQuadrants (l u a b : ℝ) (x y : Fin n → ℝ) (adj : Bool) : Fin 4 → ℝ :=
  let q := adjustedQuadrants l u a b x y adj
  fun c => normQ0 (∑ j, q j) (q c) (quadrants 0 0 a b x y c)

def normalizedMatrix (l u : ℝ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ)
    (adj : Bool) : PMResult k :=
  let q := fun c => Matrix.of fun i j => normalizedQuadrants l u (t i) (t j) (colOf X i) (colOf X j)
    adj c
  ⟨q 1, q 3, q 2, q 0, q 1 + q 0 - q 3 - q 2⟩

private lemma quadrant_nonneg (l u a b : ℝ) (x y : Fin n → ℝ) (c : Fin 4) :
    0 ≤ quadrants l u a b x y c := by
  fin_cases c <;> refine div_nonneg (sum_nonneg fun i _ => mul_nonneg ?_ ?_) n.cast_nonneg <;>
    first | apply rLowerDev_nonneg | apply rUpperDev_nonneg

private lemma adjustment_pos (adj : Bool) (n : ℕ) (l u : ℝ) : 0 < adjustment adj n l u := by
  unfold adjustment
  split_ifs with h
  · have hn : (1 : ℝ) < n := by exact_mod_cast h.2.1
    exact div_pos (by linarith) (by linarith)
  · exact one_pos

private lemma zero_quadrants (hn : 0 < n) (a b : ℝ) (x y : Fin n → ℝ) :
    ∑ c, quadrants 0 0 a b x y c = 1 := by
  obtain ⟨hl, hu, hdl, hdu⟩ := rCo_natCast 0 0 a b x y
  simpa [quadrants, Fin.sum_univ_succ, add_assoc, ← hl, ← hu, ← hdl, ← hdu] using
    quadrants_sum_eq_one hn a b x y

private lemma normalized_sum (hn : 0 < n) (l u a b : ℝ) (x y : Fin n → ℝ) (adj : Bool) :
    ∑ c, normalizedQuadrants l u a b x y adj c = 1 := by
  let q := fun c => adjustment adj n l u * quadrants l u a b x y c
  have h : 0 ≤ ∑ c, q c := sum_nonneg fun c _ =>
    mul_nonneg (adjustment_pos ..).le (quadrant_nonneg ..)
  change (∑ c, normQ0 (∑ j, q j) (q c) (quadrants 0 0 a b x y c)) = 1
  rcases h.eq_or_lt with he | hp
  · simp [← he, normQ0, zero_quadrants hn]
  · simp [normQ0, hp, ← sum_div, hp.ne']

private lemma normalized_nonneg (l u a b : ℝ) (x y : Fin n → ℝ) (adj : Bool) (c : Fin 4) :
    0 ≤ normalizedQuadrants l u a b x y adj c := by
  unfold normalizedQuadrants adjustedQuadrants normQ0
  split_ifs with h
  · exact div_nonneg (mul_nonneg (adjustment_pos ..).le (quadrant_nonneg ..)) h.le
  · exact quadrant_nonneg ..
  · exact le_rfl

def rawMatrix (l u : ℝ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ)
    (adj : Bool) : PMResult k :=
  let q := fun c => Matrix.of fun i j =>
    adjustedQuadrants l u (t i) (t j) (colOf X i) (colOf X j) adj c
  ⟨q 1, q 3, q 2, q 0, q 1 + q 0 - q 3 - q 2⟩

end A9_CoMoments

theorem claim_C571 (X : Matrix (Fin n) (Fin k) ℝ) (l u : ℝ) :
    (∀ (t : Fin k → ℝ) (adj : Bool) (i j : Fin k),
      let Q := rawMatrix l u t X adj
      ![Q.clpm i j, Q.cupm i j, Q.dlpm i j, Q.dupm i j] =
        (if adj = true ∧ 1 < n ∧ 0 < l ∧ 0 < u then (n : ℝ) / (n - 1) else 1) •
          quadrants l u (t i) (t j) (colOf X i) (colOf X j)) ∧
    (∀ (v w : ℕ) (t : Fin k → ℝ) (adj : Bool),
      rawMatrix v w t X adj = PMmatrix v w t X adj false) ∧
    (1 < n → (rawMatrix 1 1 (colMeans X) X true).cov = sampleCovMatrix X) := by
  have hnat (v w : ℕ) (t : Fin k → ℝ) (adj : Bool) :
      rawMatrix v w t X adj = PMmatrix v w t X adj false := by
    have hs : adjustment adj n v w = pmScale adj n v w := by simp [adjustment, pmScale]
    simp only [rawMatrix, PMmatrix, Bool.false_eq_true, ite_false, PMResult.mk.injEq]
    refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> ext i j <;>
      simp [adjustedQuadrants, quadrants, hs, covM, cupmM, dupmM, dlpmM, clpmM,
        (rCo_natCast v v ..).1, (rCo_natCast w w ..).2.1,
        (rCo_natCast v w ..).2.2.1, (rCo_natCast v w ..).2.2.2]
  refine ⟨?_, hnat, fun hn => ?_⟩
  · intro t adj i j Q
    funext c
    fin_cases c <;> rfl
  · rw [← Nat.cast_one (R := ℝ), hnat]
    exact PMmatrix_cov_eq_sampleCov hn X

theorem claim_C573 (hn : 0 < n) (l u : ℝ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ)
    (adj : Bool) (i j : Fin k) :
    let Q := normalizedMatrix l u t X adj
    let x := colOf X i
    let y := colOf X j
    Q.clpm i j + Q.cupm i j + Q.dlpm i j + Q.dupm i j = 1 ∧
    ((∑ c, adjustment adj n l u * quadrants l u (t i) (t j) x y c) = 0 →
      ![Q.clpm i j, Q.cupm i j, Q.dlpm i j, Q.dupm i j] = quadrants 0 0 (t i) (t j) x y) := by
  refine ⟨?_, ?_⟩
  · simpa [normalizedMatrix, Fin.sum_univ_succ, add_assoc] using
      normalized_sum hn l u (t i) (t j) (colOf X i) (colOf X j) adj
  · intro h
    funext c
    fin_cases c <;> simp [normalizedMatrix, normalizedQuadrants, adjustedQuadrants, h, normQ0]

theorem claim_C575 (hn : 0 < n) (l u : ℝ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ)
    (adj : Bool) (i j : Fin k) :
    let Q := normalizedMatrix l u t X adj
    Q.cov = Q.cupm + Q.clpm - Q.dupm - Q.dlpm ∧ Q.cov i j ∈ Set.Icc (-1) 1 := by
  refine ⟨rfl, ?_⟩
  let q := normalizedQuadrants l u (t i) (t j) (colOf X i) (colOf X j) adj
  have hq (c : Fin 4) : 0 ≤ q c := normalized_nonneg ..
  have hs : ∑ c, q c = 1 := normalized_sum hn ..
  simp only [Fin.sum_univ_four] at hs
  change -1 ≤ q 1 + q 0 - q 3 - q 2 ∧ q 1 + q 0 - q 3 - q 2 ≤ 1
  constructor <;> linarith [hq 0, hq 1, hq 2, hq 3]

namespace A9_CoMoments

def ndValues (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (norm : Bool) :
    Fin 3 → ℝ := ![coLPMnD_R X t d norm, coUPMnD_R X t d norm, dpmND_R X t d norm]

end A9_CoMoments

theorem claim_C578 (hn : 0 < n) (hk : 0 < k) {d : ℝ} (hd : 0 ≤ d)
    (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (c : Fin 3) :
    let raw := ndValues X t d false
    let total := ∑ j, raw j
    ndValues X t d true c =
      (if d = 0 then raw c else if total = 0 then ndValues X t 0 false c else raw c / total) ∧
    ndValues X t d true c ∈ Set.Icc 0 1 := by
  have hraw (e : ℝ) (he : 0 ≤ e) (j : Fin 3) : ndValues X t e false j = cellND X t e j := by
    fin_cases j <;> simp [ndValues, cellND, coLPMnD_R_eq he, coUPMnD_R_eq he, dpmND_R_eq he]
  have hs : ∑ j, ndValues X t d false j = totalND X t d := by
    simp only [hraw d hd, Fin.sum_univ_three]
    exact cellND_sum ..
  dsimp only
  rw [hs, hraw d hd, hraw 0 le_rfl]
  have hz : totalND X t 0 = 1 := by rw [totalND_eq_absND hk, absND_zero hn]
  have ht : 0 ≤ totalND X t d :=
    add_nonneg (add_nonneg (clpmND_nonneg ..) (cupmND_nonneg ..)) (dpmND_nonneg ..)
  rcases hd.eq_or_lt with he | hp
  · subst d
    have h := normND_mem_Icc X t 0 c
    simp [normND, hz] at h
    fin_cases c <;> simpa [ndValues, cellND, coLPMnD_R_eq le_rfl, coUPMnD_R_eq le_rfl,
      dpmND_R_eq le_rfl] using h
  · have hval : ndValues X t d true c = normND X t d c := by
      fin_cases c <;> simp [ndValues, coLPMnD_R_eq hd, coUPMnD_R_eq hd, dpmND_R_eq hd,
        hp.ne', normND, cellND]
    rw [hval]
    refine ⟨?_, normND_mem_Icc ..⟩
    rcases ht.eq_or_lt with he | h
    · simp [normND, ← he, hp.ne']
    · simp [normND, h, h.ne', hp.ne']

open scoped Classical in
theorem claim_C579 (hn : 0 < n) (hk : 0 < k) (X : Matrix (Fin n) (Fin k) ℝ)
    (t : Fin k → ℝ) (norm : Bool) {d : ℝ} (hd : 0 < d) :
    coLPMnD_R X t 0 norm = (#{i | ∀ j, X i j ≤ t j} : ℝ) / n ∧
    coUPMnD_R X t 0 norm = (#{i | ∀ j, t j < X i j} : ℝ) / n ∧
    coLPMnD_R X t 0 norm + coUPMnD_R X t 0 norm + dpmND_R X t 0 norm = 1 ∧
    coLPMnD_R X t d false + coUPMnD_R X t d false + dpmND_R X t d false =
      (∑ i, ∏ j, |X i j - t j| ^ d) / n := by
  refine ⟨coLPMnD_R_zero .., coUPMnD_R_zero .., degree_zero_partition hn hk X t norm, ?_⟩
  simpa [coLPMnD_R_eq hd.le, coUPMnD_R_eq hd.le, dpmND_R_eq hd.le, totalND, absND, absRow] using
    totalND_eq_absND hk X t d

theorem claim_C580 {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin 2) ℝ) (t : Fin 2 → ℝ) :
    let x := colOf X 0
    let y := colOf X 1
    (∀ (k : ℕ) (Y : Matrix (Fin n) (Fin k) ℝ) (a : Fin k → ℝ),
      coLPMnD_R Y a d false = (∑ i, ∏ j, rLowerDev d (a j) (Y i j)) / n) ∧
    coLPMnD_R X t d false = rCoLPM d d (t 0) (t 1) x y ∧
    coUPMnD_R X t d false = rCoUPM d d (t 0) (t 1) x y ∧
    dpmND_R X t d false = rDLPM d d (t 0) (t 1) x y + rDUPM d d (t 0) (t 1) x y := by
  refine ⟨fun k Y a => ?_, ?_⟩
  · simp [coLPMnD_R_eq hd, clpmND, lowRow]
  simp only [coLPMnD_R_eq hd, coUPMnD_R_eq hd, dpmND_R_eq hd, Bool.false_eq_true, and_false,
    ite_false, dpmND_eq_abs_sub (by decide : 0 < 2)]
  simp only [clpmND, cupmND, absND, absRow, lowRow, upRow, Fin.prod_univ_two, rCoLPM, rCoUPM,
    rDLPM, rDUPM, colOf, true_and]
  rw [div_sub_div_same, div_sub_div_same, ← add_div, ← sum_sub_distrib, ← sum_sub_distrib,
    ← sum_add_distrib]
  congr 1
  apply sum_congr rfl
  intro i _
  rw [abs_sub_comm (X i 0), abs_sub_comm (X i 1), ← rLowerDev_add_rUpperDev d (t 0),
    ← rLowerDev_add_rUpperDev d (t 1)]
  ring

theorem claim_C583 {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) :
    dpmND_R X t d false =
      (∑ i, if (∀ j, X i j ≤ t j) ∨ (∀ j, t j < X i j) then 0 else ∏ j, |X i j - t j| ^ d) /
        n := by
  simp [dpmND_R_eq hd, dpmND, dpmRow, absRow]

theorem claim_C601 (dx dy l u a b : ℝ) (x y : Fin n → ℝ) :
    rCoLPM dx dy a b x y = (∑ i, lowKernel dx a (x i) * lowKernel dy b (y i)) / n ∧
    rCoUPM dx dy a b x y = (∑ i, upKernel dx a (x i) * upKernel dy b (y i)) / n ∧
    rDLPM l u a b x y = (∑ i, upKernel u a (x i) * lowKernel l b (y i)) / n ∧
    rDUPM l u a b x y = (∑ i, lowKernel l a (x i) * upKernel u b (y i)) / n :=
  rCo_eq_kernels ..

namespace A9_CoMoments

def rootStop (eps lo hi b c fb : ℝ) : Prop :=
  |(c - b) / 2| ≤ 2 * eps * |b| + eps ^ 2 * (hi - lo) / 2 ∨ fb = 0

def rootValue (d : ℕ) (p b : ℝ) (x : Fin n → ℝ) : ℝ :=
  (n : ℝ) * ((1 - p) * LPM d b x - p * UPM d b x)

private lemma scaledRootValue (d : ℕ) (p t : ℝ) (x : Fin n → ℝ) {a : ℝ} (ha : 0 < a) :
    rootValue d p (a * t) (fun i => a * x i) = a ^ d * rootValue d p t x := by
  have hl (v : ℝ) : lowerDev d (a * t) (a * v) = a ^ d * lowerDev d t v := by
    simp [lowerDev, mul_le_mul_iff_right₀ ha, ← mul_sub, mul_pow, mul_ite]
  have hu (v : ℝ) : upperDev d (a * t) (a * v) = a ^ d * upperDev d t v := by
    simp [upperDev, mul_lt_mul_iff_right₀ ha, ← mul_sub, mul_pow, mul_ite]
  simp only [rootValue, LPM, UPM, hl, hu, ← mul_sum]
  ring

end A9_CoMoments

theorem claim_C615 (dx dy l u tx ty : ℝ) (x y : Fin n → ℝ)
    (d : ℕ) (p eps gapTol lo hi b c : ℝ) {a : ℝ} (ha : 0 < a) :
    (rCoLPM dx dy tx ty x y = (∑ i, lowKernel dx tx (x i) * lowKernel dy ty (y i)) / n ∧
     rCoUPM dx dy tx ty x y = (∑ i, upKernel dx tx (x i) * upKernel dy ty (y i)) / n ∧
     rDLPM l u tx ty x y = (∑ i, upKernel u tx (x i) * lowKernel l ty (y i)) / n ∧
     rDUPM l u tx ty x y = (∑ i, lowKernel l tx (x i) * upKernel u ty (y i)) / n) ∧
    (|LPMratio d (a * b) (fun i => a * x i) - p| ≤ gapTol ↔ |LPMratio d b x - p| ≤ gapTol) ∧
    (rootStop eps (a * lo) (a * hi) (a * b) (a * c) (rootValue d p (a * b) (fun i => a * x i)) ↔
      rootStop eps lo hi b c (rootValue d p b x)) := by
  refine ⟨claim_C601 .., ?_, ?_⟩
  · simpa using Iff.of_eq (congrArg (fun v : ℝ => |v - p| ≤ gapTol) (LPMratio_affine d ha 0 b x))
  · rw [scaledRootValue d p b x ha]
    have hstep : |(a * c - a * b) / 2| = a * |(c - b) / 2| := by
      rw [← mul_sub, mul_div_assoc, abs_mul, abs_of_pos ha]
    have htol : 2 * eps * |a * b| + eps ^ 2 * (a * hi - a * lo) / 2 =
        a * (2 * eps * |b| + eps ^ 2 * (hi - lo) / 2) := by
      rw [abs_mul, abs_of_pos ha]
      ring
    simp only [rootStop, hstep, htol, mul_le_mul_iff_right₀ ha, mul_eq_zero,
      pow_ne_zero d ha.ne', false_or]

end

end NNS.ClaimsA
