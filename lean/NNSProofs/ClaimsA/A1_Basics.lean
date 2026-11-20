module

public import NNSProofs.GapCoMoments
public import NNSProofs.DominanceAlgo
public import NNSProofs.GapProjection
public import NNSProofs.GapCertainty
public import Mathlib.Analysis.SpecialFunctions.Sqrt

namespace NNS.ClaimsA

open MeasureTheory ProbabilityTheory Finset Filter Set Topology
open scoped ENNReal

@[expose] public noncomputable section

variable {n m : ℕ}

theorem claim_C002 (x t : ℝ) : 0 ≤ max (x - t) 0 ∧ 0 ≤ max (t - x) 0 :=
  ⟨le_max_right _ _, le_max_right _ _⟩

theorem claim_C013 (μ : Measure ℝ) [IsProbabilityMeasure μ] (t : ℝ) :
    mLPM 0 t μ = μ.real (Iic t) ∧ mUPM 0 t μ = μ.real (Ioi t) ∧
      mLPM 0 t μ + mUPM 0 t μ = 1 :=
  ⟨(mLPM_zero_eq_cdf t).trans (cdf_eq_real _ _), mUPM_zero_eq t, mLPM_zero_add_mUPM_zero t⟩

theorem claim_C039 (μ : Measure ℝ) [IsProbabilityMeasure μ] (t : ℝ) :
    mLPM 0 t μ + mUPM 0 t μ = 1 ∧ cdf μ t + mUPM 0 t μ = 1 := by
  rw [← mLPM_zero_eq_cdf]
  exact ⟨mLPM_zero_add_mUPM_zero t, mLPM_zero_add_mUPM_zero t⟩

theorem claim_C077 {Ω α : Type*} [MeasurableSpace Ω] (P : Measure Ω) [IsProbabilityMeasure P]
    (X : Ω → α) (f : α → ℝ) (h : Integrable (f ∘ X) P) :
    ∫ ω, f (X ω) ∂P = mUPM 1 0 (P.map (f ∘ X)) - mLPM 1 0 (P.map (f ∘ X)) := by
  have hm := h.aemeasurable
  have := (Measure.isProbabilityMeasure_map_iff hm).mpr inferInstance
  have hi : Integrable (fun y : ℝ => y) (P.map (f ∘ X)) :=
    (integrable_map_measure aestronglyMeasurable_id hm).mpr h
  rw [← mean_eq_mUPM_one_sub_mLPM_one hi]
  exact (integral_map hm aestronglyMeasurable_id).symm

theorem claim_C112 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : Integrable id μ) :
    ∫ x, x ∂μ = mUPM 1 0 μ - mLPM 1 0 μ := mean_eq_mUPM_one_sub_mLPM_one h

theorem claim_C159 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : Integrable id μ) :
    ∫ x, x ∂μ = mUPM 1 0 μ - mLPM 1 0 μ := claim_C112 μ h

theorem claim_C160 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : Integrable id μ) (t : ℝ) :
    ∫ x, x - t ∂μ = mUPM 1 t μ - mLPM 1 t μ := by
  change Integrable (fun x : ℝ => x) μ at h
  rw [mUPM_one_sub_mLPM_one h, integral_sub h (integrable_const t)]
  simp

theorem claim_C193 (hn : 0 < n) (x : Fin n → ℝ) (t : ℝ) :
    cdf (empMeasure x) t = (∑ i, if x i ≤ t then (1 : ℝ) else 0) / n ∧
      cdf (empMeasure x) t = (#{i | x i ≤ t} : ℝ) / n := by
  rw [cdf_empMeasure hn]
  exact ⟨by simp [LPM, lowerDev], LPM_zero_eq_ecdf t x⟩

theorem claim_C195 (hn : 0 < n) (x : Fin n → ℝ) (t : ℝ) :
    cdf (empMeasure x) t = LPM 0 t x ∧
      cdf (empMeasure x) t = (∑ i, lowerDev 0 t (x i)) / n := by
  exact ⟨cdf_empMeasure hn x t, cdf_empMeasure hn x t⟩

theorem claim_C198 (hn : 0 < n) (x : Fin n → ℝ) :
    (fun t => LPM 0 t x) = (fun t => cdf (empMeasure x) t) :=
  funext fun t => (cdf_empMeasure hn x t).symm

theorem claim_C202 (hn : 0 < n) (x : Fin n → ℝ) (t : ℝ) :
    (t ∉ Set.range x → HasDerivAt (fun s => cdf (empMeasure x) s) 0 t) ∧
    (t ∈ Set.range x → ¬ ContinuousAt (fun s => cdf (empMeasure x) s) t) := by
  constructor
  · intro ht
    have each (i : Fin n) : ∀ᶠ s in 𝓝 t, x i ≤ s ↔ x i ≤ t := by
      rcases lt_or_gt_of_ne (fun h => ht ⟨i, h⟩) with h | h
      · filter_upwards [lt_mem_nhds h] with s hs using iff_of_true hs.le h.le
      · filter_upwards [gt_mem_nhds h] with s hs using iff_of_false (not_le.mpr hs) (not_le.mpr h)
    apply (hasDerivAt_const t (cdf (empMeasure x) t)).congr_of_eventuallyEq
    filter_upwards [Filter.eventually_all.mpr each] with s hs
    simp only [(claim_C193 hn x _).2, hs]
  · rintro ⟨i, rfl⟩ hc
    have hl := cdfCurve_zero_tendsto_left (x i) x
    have he : (fun s => cdf (empMeasure x) s) = (fun s => cdfCurve 0 s x) :=
      funext (cdf_empMeasure hn x)
    rw [he] at hc
    have heq := tendsto_nhds_unique hl hc.continuousWithinAt.tendsto
    rw [cdfCurve_zero_eq_lt_add_eq] at heq
    have hp : 0 < (#{j | x j = x i} : ℝ) / n := by
      apply div_pos _ (by exact_mod_cast hn)
      exact_mod_cast Finset.card_pos.mpr ⟨i, by simp⟩
    linarith

theorem claim_C215 (hn : 0 < n) (x : Fin n → ℝ) (j : Fin n) (v t : ℝ) :
    let y := Function.update x j v
    cdf (empMeasure y) t - cdf (empMeasure x) t =
      ((if v ≤ t then 1 else 0) - (if x j ≤ t then 1 else 0)) / (n : ℝ) ∧
    |cdf (empMeasure y) t - cdf (empMeasure x) t| ≤ 1 / (n : ℝ) := by
  dsimp only
  have he : cdf (empMeasure (Function.update x j v)) t - cdf (empMeasure x) t =
      ((if v ≤ t then 1 else 0) - (if x j ≤ t then 1 else 0)) / (n : ℝ) := by
    rw [(claim_C193 hn _ _).1, (claim_C193 hn _ _).1, ← sub_div, ← sum_sub_distrib]
    congr 1
    rw [sum_eq_single j]
    · simp
    · intro i _ hij
      simp [Function.update_of_ne hij]
    · simp
  refine ⟨he, ?_⟩
  rw [he, abs_div, abs_of_nonneg (n.cast_nonneg (α := ℝ))]
  apply div_le_div_of_nonneg_right _ n.cast_nonneg
  split_ifs <;> norm_num

theorem claim_C217 (hn : 0 < n) (x : Fin n → ℝ) (t : ℝ) :
    cdf (empMeasure x) t = (∑ i, if x i ≤ t then (1 : ℝ) else 0) / n :=
  (claim_C193 hn x t).1

namespace A1_Basics

def tail (upper : Bool) (t : ℝ) : Set ℝ := if upper then Ioi t else Iic t

def cell (μ : Measure (ℝ × ℝ)) (a b : ℝ) : Bool → Bool → ℝ
  | false, false => mCoLPM 0 0 a b μ
  | false, true => mDUPM 0 0 a b μ
  | true, false => mDLPM 0 0 a b μ
  | true, true => mCoUPM 0 0 a b μ

private theorem partition (a b : ℝ) (p : ℝ × ℝ) :
    ∃! q : Bool × Bool, p ∈ tail q.1 a ×ˢ tail q.2 b := by
  refine ⟨(decide (a < p.1), decide (b < p.2)), ?_, ?_⟩
  · by_cases h : a < p.1 <;> by_cases k : b < p.2 <;> simp [tail, h, k, le_of_not_gt, Prod.le_def]
  · rintro ⟨i, j⟩ h
    cases i <;> cases j <;> simp_all [tail, Set.mem_prod, not_lt, Prod.le_def]

end A1_Basics

open A1_Basics

theorem claim_C254 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] (a b : ℝ) :
    (∀ p : ℝ × ℝ, ∃! q : Bool × Bool, p ∈ tail q.1 a ×ˢ tail q.2 b) ∧
    mCoLPM 0 0 a b μ + mCoUPM 0 0 a b μ + mDUPM 0 0 a b μ + mDLPM 0 0 a b μ = 1 ∧
    ∀ i j : Bool,
      let X := if i then xGt a else xLe a
      let Y := if j then yGt b else yLe b
      let px := if i then mUPMx 0 a μ else mLPMx 0 a μ
      let py := if j then mUPMy 0 b μ else mLPMy 0 b μ
      cell μ a b i j = μ.real (tail i a ×ˢ tail j b) ∧
      px = μ.real X ∧ py = μ.real Y ∧
      (0 < px → (μ[|X]).real Y = cell μ a b i j / px) ∧
      (0 < py → (μ[|Y]).real X = cell μ a b i j / py) ∧
      (μ[|X]).real Y * px = cell μ a b i j ∧
      (μ[|Y]).real X * py = cell μ a b i j := by
  refine ⟨partition a b, quadrants_sum_eq_one_measure μ a b, ?_⟩
  intro i j
  cases i <;> cases j
  · exact ⟨mCoLPM_zero .., mLPMx_zero .., mLPMy_zero ..,
      fun _ => cond_yLe_xLe .., fun _ => cond_xLe_yLe .., bayes_lower μ a b⟩
  · exact ⟨mDUPM_zero .., mLPMx_zero .., mUPMy_zero ..,
      fun _ => cond_yGt_xLe .., fun _ => cond_xLe_yGt .., bayes_DUPM μ a b⟩
  · exact ⟨mDLPM_zero .., mUPMx_zero .., mLPMy_zero ..,
      fun _ => cond_yLe_xGt .., fun _ => cond_xGt_yLe .., bayes_DLPM μ a b⟩
  · exact ⟨mCoUPM_zero .., mUPMx_zero .., mUPMy_zero ..,
      fun _ => cond_yGt_xGt .., fun _ => cond_xGt_yGt .., bayes_upper μ a b⟩

theorem claim_C264 (x t : ℝ) : max (t - x) 0 - max (x - t) 0 = t - x := by
  rw [← lowerDev_one, ← upperDev_one]
  linarith [upperDev_sub_lowerDev_one t x]

theorem claim_C265 (hn : 0 < n) (x : Fin n → ℝ) (t : ℝ) :
    LPM 1 t x - UPM 1 t x = t - mean x ∧
    ((∃ i, x i ≠ t) → LPMratio 1 t x - 1 / 2 =
      (t - mean x) / (2 * (LPM 1 t x + UPM 1 t x))) := by
  have hd := UPM_one_sub_LPM_one hn t x
  refine ⟨by linarith, fun h => ?_⟩
  have hpos := (ratio_denom_pos (d := 1) h).ne'
  have hnt : ¬ (0 < n ∧ ∀ i, x i = t) := fun H => h.elim fun i hi => hi (H.2 i)
  simp only [LPMratio, one_ne_zero, hnt, ite_false]
  field_simp
  linarith

private theorem scaled_distance (x : Fin n → ℝ) (t : ℝ) (h : ∃ i, x i ≠ t) :
    LPMratio 1 t x - 1 / 2 = (t - mean x) / (2 * mean (fun i => |x i - t|)) := by
  have hsum : LPM 1 t x + UPM 1 t x = mean (fun i => |x i - t|) := by
    rw [add_comm, UPM, LPM, mean, ← add_div, ← sum_add_distrib]
    simp only [upperDev_one_add_lowerDev_one]
  exact hsum ▸ (claim_C265 (Fin.pos h.choose) x t).2 h

theorem claim_C279 (x : Fin n → ℝ) (y : Fin m → ℝ) (t : ℝ) :
    ((∃ i, x i ≠ t) → LPMratio 1 t x - 1 / 2 =
      (t - mean x) / (2 * mean (fun i => |x i - t|))) ∧
    ((∃ j, y j ≠ t) → LPMratio 1 t y - 1 / 2 =
      (t - mean y) / (2 * mean (fun j => |y j - t|))) :=
  ⟨scaled_distance x t, scaled_distance y t⟩

theorem claim_C296 (x : Fin n → ℝ) (y : Fin m → ℝ) :
    let z := GapCertainty.pooledMean x y
    ((∃ i, x i ≠ z) → LPMratio 1 z x - 1 / 2 =
      (z - mean x) / (2 * mean (fun i => |x i - z|))) ∧
    ((∃ j, y j ≠ z) → LPMratio 1 z y - 1 / 2 =
      (z - mean y) / (2 * mean (fun j => |y j - z|))) :=
  claim_C279 x y _

private theorem mean_order {μ ν : Measure ℝ} [IsProbabilityMeasure μ] [IsProbabilityMeasure ν]
    (hμ : MemLp id 2 μ) (hν : MemLp id 2 ν) (h : ∀ t, mLPM 2 t μ ≤ mLPM 2 t ν) :
    (∫ x, x ∂ν) ≤ ∫ x, x ∂μ := by
  let gap := (∫ x, x ∂ν) - ∫ x, x ∂μ
  let c := variance (fun x : ℝ => x) ν + (∫ x, x ∂ν) ^ 2 -
    variance (fun x : ℝ => x) μ - (∫ x, x ∂μ) ^ 2 + mUPM 2 0 μ
  have bound (t : ℝ) (ht : 0 ≤ t) : t * (2 * gap) ≤ c := by
    have upper := mUPM_antitone (r := 2) (μ := μ) hμ ht
    dsimp only [gap, c]
    change mUPM 2 t μ ≤ mUPM 2 0 μ at upper
    nlinarith [mUPM_two_add_mLPM_two hμ t, mUPM_two_add_mLPM_two hν t, h t, mUPM_nonneg 2 t ν]
  by_contra! hbad
  have hg : 0 < 2 * gap := mul_pos (by norm_num) (sub_pos.mpr hbad)
  have habs := bound ((|c| + 1) / (2 * gap)) (div_nonneg (by positivity) hg.le)
  rw [div_mul_cancel₀ _ hg.ne'] at habs
  linarith [le_abs_self c]

private theorem prudent_not_dara :
    ∃ u : ℝ → ℝ, (∀ x ∈ Icc (0 : ℝ) 1,
      0 < deriv u x ∧ deriv (deriv u) x < 0 ∧ 0 < deriv (deriv (deriv u)) x) ∧
      ¬ AntitoneOn (fun x => -deriv (deriv u) x / deriv u x) (Icc (0 : ℝ) 1) := by
  let u (x : ℝ) := 3 * x - x ^ 2 + x ^ 3 / 12
  have h₁ : deriv u = fun x => 3 - 2 * x + x ^ 2 / 4 := by
    funext x
    apply HasDerivAt.deriv
    convert (((hasDerivAt_id x).const_mul 3).sub ((hasDerivAt_id x).pow 2)).add
      (((hasDerivAt_id x).pow 3).div_const 12) using 1 <;> first | rfl | (simp; ring)
  have h₂ : deriv (deriv u) = fun x => -2 + x / 2 := by
    rw [h₁]
    funext x
    apply HasDerivAt.deriv
    convert ((hasDerivAt_const x 3).sub ((hasDerivAt_id x).const_mul 2)).add
      (((hasDerivAt_id x).pow 2).div_const 4) using 1 <;> first | rfl | (simp; ring)
  have h₃ : deriv (deriv (deriv u)) = fun _ => 1 / 2 := by
    rw [h₂]
    funext x
    apply HasDerivAt.deriv
    convert ((hasDerivAt_id x).div_const 2).const_add (-2) using 1
    rfl
  refine ⟨u, ?_, ?_⟩
  · intro x hx
    rw [h₃, h₂, h₁]
    obtain ⟨_, hx⟩ := hx
    constructor
    · nlinarith [sq_nonneg x]
    · constructor <;> linarith
  · intro h
    have bad := h (a := 0) (b := 1) (by simp) (by simp) (by norm_num)
    rw [h₂, h₁] at bad
    norm_num at bad

namespace A1_Basics

structure PrudentUtility (u : ℝ → ℝ) : Prop where
  diff : Differentiable ℝ u
  diff₁ : Differentiable ℝ (deriv u)
  diff₂ : Differentiable ℝ (deriv (deriv u))
  first : ∀ x, 0 ≤ deriv u x
  second : ∀ x, deriv (deriv u) x ≤ 0
  third : ∀ x, 0 ≤ deriv (deriv (deriv u)) x

private def smoothHinge (t e x : ℝ) : ℝ :=
  -((Real.sqrt ((x - t)^2 + e) - (x - t)) / 2)^2

private theorem prudent_smoothHinge (t : ℝ) {e : ℝ} (he : 0 < e) :
    PrudentUtility (smoothHinge t e) := by
  let y (x : ℝ) := x - t
  let r (x : ℝ) := Real.sqrt ((y x)^2 + e)
  let q (x : ℝ) := 1 - y x / r x
  have hpos (x) : 0 < (y x)^2 + e := by positivity
  have hr (x) : 0 < r x := Real.sqrt_pos.2 (hpos x)
  have hsq (x) : (r x)^2 = (y x)^2 + e := Real.sq_sqrt (hpos x).le
  have hy (x) : HasDerivAt y 1 x := (hasDerivAt_id x).sub_const t
  have hdr (x) : HasDerivAt r (y x / r x) x := by
    convert ((((hy x).pow 2).add_const e).sqrt (hpos x).ne') using 1
    dsimp [r]
    ring
  have hq (x) : 0 ≤ q x ∧ q x ≤ 2 := by
    have hyabs : |y x| ≤ r x := by
      rw [← sq_le_sq₀ (abs_nonneg _) (hr x).le, sq_abs, hsq]
      linarith
    obtain ⟨hl, hu⟩ := abs_le.mp hyabs
    dsimp [q]
    constructor
    · linarith [(div_le_one (hr x)).2 hu]
    · linarith [(le_div_iff₀ (hr x)).2 (by linarith : -1 * r x ≤ y x)]
  have hdq (x) : HasDerivAt q (-(q x * (2 - q x)) / r x) x := by
    convert ((hy x).div (hdr x) (hr x).ne').const_sub 1 using 1
    dsimp [q]
    field_simp [(hr x).ne']; ring
  let v x := (r x - y x) * q x / 2
  let w x := -(q x)^2 * (3 - q x) / 2
  have h₁ (x) : HasDerivAt (smoothHinge t e) (v x) x := by
    convert (((((hdr x).sub (hy x)).div_const 2).pow 2).neg) using 1
    · rfl
    · simp only [v, q, Pi.sub_apply]; ring
  have h₂ (x) : HasDerivAt v (w x) x := by
    convert ((((hdr x).sub (hy x)).mul (hdq x)).div_const 2) using 1
    dsimp [w, q]
    field_simp [(hr x).ne']; ring
  have h₃ (x) : HasDerivAt w (3 * (q x * (2 - q x))^2 / (2 * r x)) x := by
    convert (((((hdq x).pow 2).neg).mul ((hdq x).const_sub 3)).div_const 2) using 1
    simp only [Pi.neg_apply, Pi.pow_apply]
    ring
  have d₁ : deriv (smoothHinge t e) = v := funext fun x => (h₁ x).deriv
  have d₂ : deriv v = w := funext fun x => (h₂ x).deriv
  constructor
  · exact fun x => (h₁ x).differentiableAt
  · rw [d₁]; exact fun x => (h₂ x).differentiableAt
  · rw [d₁, d₂]; exact fun x => (h₃ x).differentiableAt
  · intro x
    rw [d₁]
    have h : r x - y x = r x * q x := by dsimp [q]; field_simp [(hr x).ne']
    dsimp only [v]
    rw [h]
    exact div_nonneg (mul_nonneg (mul_nonneg (hr x).le (hq x).1) (hq x).1) (by norm_num)
  · intro x
    rw [d₁, d₂]
    exact div_nonpos_of_nonpos_of_nonneg
      (mul_nonpos_of_nonpos_of_nonneg (neg_nonpos.mpr (sq_nonneg _)) (by linarith [(hq x).2]))
      (by norm_num)
  · intro x
    rw [d₁, d₂, (h₃ x).deriv]
    positivity

private theorem hinge_bound (t x : ℝ) {e : ℝ} (he : 0 < e) (he' : e ≤ 1) :
    ‖smoothHinge t e x‖ ≤ (|x| + |t| + 1)^2 := by
  let y := x - t
  let r := Real.sqrt (y^2 + e)
  have hr : 0 ≤ r := Real.sqrt_nonneg _
  have hs : r^2 = y^2 + e := Real.sq_sqrt (by positivity)
  have h₁ : |y| ≤ r := by
    rw [← sq_le_sq₀ (abs_nonneg _) hr, sq_abs, hs]; linarith
  have h₂ : r ≤ |y| + 1 := by
    rw [← sq_le_sq₀ hr (by positivity), hs]
    nlinarith [sq_abs y, abs_nonneg y]
  rw [Real.norm_eq_abs, smoothHinge, abs_neg, abs_of_nonneg (sq_nonneg _)]
  change ((r - y) / 2)^2 ≤ _
  apply pow_le_pow_left₀ (by linarith [(abs_le.mp h₁).2])
  linarith [neg_le_abs y, abs_sub x t]

private theorem moment_of_utility_order {μ ν : Measure ℝ} [IsProbabilityMeasure μ]
    [IsProbabilityMeasure ν] (hμ : MemLp id 2 μ) (hν : MemLp id 2 ν)
    (pref : ∀ u, PrudentUtility u → Integrable u μ → Integrable u ν →
      (∫ x, u x ∂ν) ≤ ∫ x, u x ∂μ) (t : ℝ) : mLPM 2 t μ ≤ mLPM 2 t ν := by
  let e (n : ℕ) : ℝ := ((n : ℝ) + 1)⁻¹
  have he (n) : 0 < e n := by dsimp [e]; positivity
  have he' (n) : e n ≤ 1 := inv_le_one_of_one_le₀ (by linarith [Nat.cast_nonneg (α := ℝ) n])
  have ht : Tendsto e atTop (𝓝 0) := tendsto_one_div_add_atTop_nhds_zero_nat.congr
    (fun n => by simp [e, one_div])
  have cont (n) := (prudent_smoothHinge t (he n)).diff.continuous
  have hi {ρ : Measure ℝ} [IsProbabilityMeasure ρ] (hρ : MemLp id 2 ρ) (n) :
      Integrable (smoothHinge t (e n)) ρ :=
    (integrable_abs_add_pow hρ t).mono' (cont n).aestronglyMeasurable
      (ae_of_all _ fun x => hinge_bound t x (he n) (he' n))
  have lim {ρ : Measure ℝ} [IsProbabilityMeasure ρ] (hρ : MemLp id 2 ρ) :
      Tendsto (fun n => ∫ x, smoothHinge t (e n) x ∂ρ) atTop (𝓝 (-mLPM 2 t ρ)) := by
    rw [mLPM, ← integral_neg]
    apply tendsto_integral_of_dominated_convergence (fun x => (|x| + |t| + 1)^2)
      (fun n => (cont n).aestronglyMeasurable) (integrable_abs_add_pow hρ t)
      (fun n => ae_of_all _ fun x => hinge_bound t x (he n) (he' n))
    apply ae_of_all
    intro x
    have hc : Continuous (fun z : ℝ => smoothHinge t z x) := by unfold smoothHinge; fun_prop
    have hz : smoothHinge t 0 x = -lowerDev 2 t x := by
      rw [smoothHinge, add_zero, Real.sqrt_sq_eq_abs, lowerDev_eq_max_pow (by decide)]
      rcases le_total x t with h | h
      · rw [abs_of_nonpos (sub_nonpos.mpr h), max_eq_left (sub_nonneg.mpr h)]; ring
      · rw [abs_of_nonneg (sub_nonneg.mpr h), max_eq_right (sub_nonpos.mpr h)]; ring
    exact hz ▸ (hc.continuousAt.tendsto.comp ht)
  have h := le_of_tendsto_of_tendsto' (lim hν) (lim hμ)
    (fun n => pref _ (prudent_smoothHinge t (he n)) (hi hμ n) (hi hν n))
  linarith
private def tangent (u : ℝ → ℝ) (b x : ℝ) := u b + deriv u b * (x - b)

private def extension (u : ℝ → ℝ) (b x : ℝ) := if x < b then u x else tangent u b x

private theorem tangent_bounds {u : ℝ → ℝ} (hu : PrudentUtility u) (b x : ℝ) :
    u x ≤ tangent u b x := by
  have anti := antitone_of_deriv_nonpos hu.diff₁ hu.second
  rcases lt_trichotomy b x with h | rfl | h
  · obtain ⟨c, hc, eq⟩ := exists_deriv_eq_slope u h hu.diff.continuous.continuousOn
      hu.diff.differentiableOn
    have := mul_le_mul_of_nonneg_right (anti hc.1.le) (sub_nonneg.mpr h.le)
    rw [eq_div_iff (sub_ne_zero.mpr h.ne')] at eq
    unfold tangent
    linarith
  · simp [tangent]
  · obtain ⟨c, hc, eq⟩ := exists_deriv_eq_slope u h hu.diff.continuous.continuousOn
      hu.diff.differentiableOn
    have := mul_le_mul_of_nonneg_right (anti hc.2.le) (sub_nonneg.mpr h.le)
    rw [eq_div_iff (sub_ne_zero.mpr h.ne')] at eq
    unfold tangent
    linarith

private theorem extension_bounds {u : ℝ → ℝ} (hu : PrudentUtility u) {b : ℝ} (hb : 0 ≤ b)
    (x : ℝ) : u x ≤ extension u b x ∧ extension u b x ≤ tangent u 0 x := by
  unfold extension
  split_ifs with hx
  · exact ⟨le_rfl, tangent_bounds hu 0 x⟩
  refine ⟨tangent_bounds hu b x, ?_⟩
  have h := tangent_bounds hu 0 b
  have d := mul_le_mul_of_nonneg_right (antitone_of_deriv_nonpos hu.diff₁ hu.second hb)
    (sub_nonneg.mpr (le_of_not_gt hx))
  unfold tangent at *
  linarith

private theorem int_tangent {ρ : Measure ℝ} [IsProbabilityMeasure ρ] (hρ : MemLp id 2 ρ)
    (u : ℝ → ℝ) (b : ℝ) : Integrable (tangent u b) ρ := by
  have hi : Integrable (fun x : ℝ => x) ρ := hρ.integrable (by norm_num)
  exact (integrable_const _).add ((hi.sub (integrable_const b)).const_mul _)

private theorem int_extension {ρ : Measure ℝ} [IsProbabilityMeasure ρ] (hρ : MemLp id 2 ρ)
    {u : ℝ → ℝ} (hi : Integrable u ρ) (b : ℝ) :
    Integrable (extension u b) ρ :=
  Integrable.piecewise (s := Iio b) measurableSet_Iio.nullMeasurableSet hi.integrableOn
    (int_tangent hρ u b).integrableOn

private theorem tail_primitive {κ : Measure ℝ} [SFinite κ] (d : ℕ) (b : ℝ)
    (hi : ∀ x, Integrable (upperDev (d + 1) x) κ) {F g : ℝ → ℝ}
    (hF : ∀ x, HasDerivAt F (g x) x) (hg : Continuous g)
    (h : ∀ x, mUPM d x κ = (Iio b).indicator g x) (x : ℝ) :
    mUPM (d + 1) x κ = ((d : ℝ) + 1) * (if x < b then F b - F x else 0) := by
  rw [mUPM_succ d x (hi x)]
  simp_rw [h]
  rw [setIntegral_indicator measurableSet_Iio, Ici_inter_Iio]
  congr 1
  split_ifs with hx
  · rw [integral_Ico_eq_integral_Ioo, ← integral_Ioc_eq_integral_Ioo,
      ← intervalIntegral.integral_of_le hx.le]
    exact intervalIntegral.integral_eq_sub_of_hasDerivAt (fun z _ => hF z)
      (hg.intervalIntegrable ..)
  · simp [Ico_eq_empty hx]

private theorem curvature_representation {u : ℝ → ℝ} (hu : PrudentUtility u) (b : ℝ) :
    ∃ κ : Measure ℝ, SFinite κ ∧ (∀ x, Integrable (upperDev 2 x) κ) ∧
      ∀ x, mUPM 2 x κ = 2 * (tangent u b x - extension u b x) := by
  -- The jump to zero at b supplies the endpoint curvature term.
  let g := deriv (deriv u)
  let k : StieltjesFunction ℝ := {
    toFun := fun x => if x < b then g x else 0
    mono' := by
      intro x y hxy
      dsimp only
      split_ifs with hx hy hy
      · exact monotone_of_deriv_nonneg hu.diff₂ hu.third hxy
      · exact hu.second x
      · exact False.elim (hx (hxy.trans_lt hy))
      · exact le_rfl
    right_continuous' := by
      intro x
      by_cases hx : x < b
      · apply (hu.diff₂.continuous.continuousWithinAt (x := x)).congr_of_eventuallyEq_of_mem
            _ (self_mem_Ici)
        filter_upwards [mem_nhdsWithin_of_mem_nhds (gt_mem_nhds hx)] with z hz
        simp [hz, g]
      · apply (continuousWithinAt_const (b := (0 : ℝ))).congr_of_eventuallyEq_of_mem _ self_mem_Ici
        filter_upwards [self_mem_nhdsWithin] with z hz
        simp [not_lt.mpr ((le_of_not_gt hx).trans hz)] }
  let κ := k.measure.restrict (Iic b)
  have hi (d : ℕ) (x : ℝ) : Integrable (upperDev d x) κ := by
    rw [upperDev_eq_indicator, integrable_indicator_iff measurableSet_Ioi, IntegrableOn]
    change Integrable _ ((k.measure.restrict _).restrict _)
    rw [Measure.restrict_restrict measurableSet_Ioi, Ioi_inter_Iic]
    exact (by fun_prop : Continuous fun z : ℝ => (z - x)^d).integrableOn_Icc.mono_set
      Ioc_subset_Icc_self
  have h₀ (x : ℝ) : mUPM 0 x κ = (Iio b).indicator (fun z => -g z) x := by
    simp only [mUPM_eq_setIntegral, pow_zero, setIntegral_const, smul_eq_mul, mul_one]
    change (k.measure.restrict (Iic b)).real (Ioi x) = _
    rw [Measure.real, Measure.restrict_apply measurableSet_Ioi, Ioi_inter_Iic, k.measure_Ioc]
    by_cases hx : x < b
    · simp [k, hx, indicator_of_mem, ENNReal.toReal_ofReal (neg_nonneg.mpr (hu.second x)), g]
    · simp [k, hx]
  have h₁ (x : ℝ) : mUPM 1 x κ =
      (Iio b).indicator (fun z => deriv u z - deriv u b) x := by
    have h := tail_primitive 0 b (hi 1) (fun z => (hu.diff₁ z).hasDerivAt.neg)
      hu.diff₂.continuous.neg h₀ x
    simpa [Set.indicator, sub_eq_add_neg, add_comm] using h
  refine ⟨κ, inferInstance, hi 2, fun x => ?_⟩
  have hd (z : ℝ) : HasDerivAt (fun s => u s - deriv u b * s) (deriv u z - deriv u b) z := by
    convert (hu.diff z).hasDerivAt.sub ((hasDerivAt_id z).const_mul (deriv u b)) using 1
    · rfl
    · ring
  have h := tail_primitive 1 b (hi 2) hd (hu.diff₁.continuous.sub continuous_const) h₁ x
  rw [h]
  unfold extension tangent
  split_ifs <;> push_cast <;> ring

private theorem integral_extension_order {μ ν : Measure ℝ} [IsProbabilityMeasure μ]
    [IsProbabilityMeasure ν] (hμ : MemLp id 2 μ) (hν : MemLp id 2 ν)
    (ord : ∀ t, mLPM 2 t μ ≤ mLPM 2 t ν) (hm : (∫ x, x ∂ν) ≤ ∫ x, x ∂μ)
    {u : ℝ → ℝ} (hu : PrudentUtility u) (iμ : Integrable u μ) (iν : Integrable u ν) (b : ℝ) :
    (∫ x, extension u b x ∂ν) ≤ ∫ x, extension u b x ∂μ := by
  obtain ⟨κ, hκ, hi, rep⟩ := curvature_representation hu b
  have swap {ρ : Measure ℝ} [IsProbabilityMeasure ρ] (hρ : MemLp id 2 ρ) (iu : Integrable u ρ) :
      (∫ t, mLPM 2 t ρ ∂κ) = 2 * ((∫ x, tangent u b x ∂ρ) - ∫ x, extension u b x ∂ρ) ∧
      Integrable (fun t => mLPM 2 t ρ) κ := by
    have eq (x t : ℝ) : upperDev 2 x t = lowerDev 2 t x := by
      rw [upperDev_eq_max_pow (by decide), lowerDev_eq_max_pow (by decide)]
    have prod : Integrable (Function.uncurry (upperDev 2)) (ρ.prod κ) := by
      apply (integrable_prod_iff (measurable_upperDev_uncurry 2).aestronglyMeasurable).2
      refine ⟨ae_of_all _ hi, ?_⟩
      simp_rw [Real.norm_eq_abs, abs_of_nonneg (upperDev_nonneg 2 _ _)]
      change Integrable (fun x => mUPM 2 x κ) ρ
      simp_rw [rep]
      exact ((int_tangent hρ u b).sub (int_extension hρ iu b)).const_mul 2
    have F := integral_integral_swap prod
    change (∫ x, mUPM 2 x κ ∂ρ) = ∫ t, ∫ x, upperDev 2 x t ∂ρ ∂κ at F
    simp_rw [eq, rep] at F
    rw [integral_const_mul, integral_sub (int_tangent hρ u b) (int_extension hρ iu b)] at F
    refine ⟨F.symm, ?_⟩
    have hj := prod.integral_prod_right
    change Integrable (fun t => ∫ x, upperDev 2 x t ∂ρ) κ at hj
    simpa only [eq, mLPM] using hj
  obtain ⟨eμ, jμ⟩ := swap hμ iμ
  obtain ⟨eν, jν⟩ := swap hν iν
  have hgap := integral_mono jμ jν ord
  rw [eμ, eν] at hgap
  have affine {ρ : Measure ℝ} [IsProbabilityMeasure ρ] (hρ : MemLp id 2 ρ) :
      (∫ x, tangent u b x ∂ρ) = u b + deriv u b * ((∫ x, x ∂ρ) - b) := by
    have hi' : Integrable (fun x : ℝ => x) ρ := hρ.integrable (by norm_num)
    have hlin : Integrable (fun x : ℝ => deriv u b * (x - b)) ρ :=
      (hi'.sub (integrable_const b)).const_mul _
    unfold tangent
    rw [integral_add (integrable_const (u b)) hlin,
      integral_const_mul, integral_sub hi' (integrable_const _)]
    simp
  rw [affine hμ, affine hν] at hgap
  linarith [mul_le_mul_of_nonneg_left hm (hu.first b)]

private theorem utility_of_moment_order {μ ν : Measure ℝ} [IsProbabilityMeasure μ]
    [IsProbabilityMeasure ν] (hμ : MemLp id 2 μ) (hν : MemLp id 2 ν)
    (ord : ∀ t, mLPM 2 t μ ≤ mLPM 2 t ν) (hm : (∫ x, x ∂ν) ≤ ∫ x, x ∂μ)
    {u : ℝ → ℝ} (hu : PrudentUtility u) (iμ : Integrable u μ) (iν : Integrable u ν) :
    (∫ x, u x ∂ν) ≤ ∫ x, u x ∂μ := by
  have lim {ρ : Measure ℝ} [IsProbabilityMeasure ρ] (hρ : MemLp id 2 ρ) (iu : Integrable u ρ) :
      Tendsto (fun n : ℕ => ∫ x, extension u n x ∂ρ) atTop (𝓝 (∫ x, u x ∂ρ)) := by
    apply tendsto_integral_of_dominated_convergence (fun x => ‖u x‖ + ‖tangent u 0 x‖)
      (fun n => (int_extension hρ iu n).aestronglyMeasurable)
      (iu.norm.add (int_tangent hρ u 0).norm)
    · intro n
      apply ae_of_all
      intro x
      obtain ⟨hl, hr⟩ := extension_bounds hu (Nat.cast_nonneg n) x
      simp only [Real.norm_eq_abs]
      exact abs_le.mpr ⟨by linarith [neg_abs_le (u x), abs_nonneg (tangent u 0 x)],
        by linarith [le_abs_self (tangent u 0 x), abs_nonneg (u x)]⟩
    · apply ae_of_all
      intro x
      apply tendsto_nhds_of_eventually_eq
      filter_upwards [tendsto_natCast_atTop_atTop.eventually (eventually_gt_atTop x)] with n hn
      simp [extension, hn]
  exact le_of_tendsto_of_tendsto' (lim hν iν) (lim hμ iμ)
    (fun n => integral_extension_order hμ hν ord hm hu iμ iν n)

private theorem finite_crossing :
    ∃ (μ ν : Measure ℝ) (_ : IsProbabilityMeasure μ) (_ : IsProbabilityMeasure ν),
      MemLp id 2 μ ∧ MemLp id 2 ν ∧ (∀ t, mLPM 2 t μ ≤ mLPM 2 t ν) ∧
      (∃ t, mLPM 2 t μ < mLPM 2 t ν) ∧ (∃ t, cdf μ t < cdf ν t) ∧
      ∃ t, cdf ν t < cdf μ t := by
  let x : Fin 1 → ℝ := ![0]
  let y : Fin 2 → ℝ := ![-1, 1]
  have hx := isProbabilityMeasure_empMeasure (by decide) x
  have hy := isProbabilityMeasure_empMeasure (by decide) y
  have h₂x := memLp_empMeasure (by decide) x 2
  have h₂y := memLp_empMeasure (by decide) y 2
  refine ⟨empMeasure x, empMeasure y, hx, hy, h₂x, h₂y, moment_of_utility_order h₂x h₂y ?_,
    ⟨0, ?_⟩, ⟨-1, ?_⟩, ⟨0, ?_⟩⟩
  · intro u hu _ _
    rw [integral_empMeasure, integral_empMeasure]
    simp [x, y, Fin.sum_univ_two]
    have hl := tangent_bounds hu 0 (-1)
    have hr := tangent_bounds hu 0 1
    unfold tangent at hl hr
    linarith
  · rw [mLPM_empMeasure, mLPM_empMeasure]
    norm_num [x, y, LPM, lowerDev, Fin.sum_univ_two, Fin.sum_univ_one]
  all_goals
    rw [cdf_empMeasure (by decide), cdf_empMeasure (by decide)]
    simp only [x, y, LPM, lowerDev, Fin.sum_univ_two, Fin.sum_univ_one]
    norm_num

end A1_Basics

theorem claim_C332 (μ ν : Measure ℝ) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν]
    (hμ : MemLp id 2 μ) (hν : MemLp id 2 ν) :
    ((∀ t : ℝ, (∫ s in Iic t, ∫ u in Iic s, cdf μ u) ≤
      ∫ s in Iic t, ∫ u in Iic s, cdf ν u) ↔ ∀ t, mLPM 2 t μ ≤ mLPM 2 t ν) ∧
    ((∀ t, mLPM 2 t μ ≤ mLPM 2 t ν) → (∫ x, x ∂ν) ≤ ∫ x, x ∂μ) ∧
    ((∀ t, mLPM 2 t μ ≤ mLPM 2 t ν) ↔ ∀ u : ℝ → ℝ, PrudentUtility u →
      Integrable u μ → Integrable u ν → (∫ x, u x ∂ν) ≤ ∫ x, u x ∂μ) ∧
    (∃ (ρ σ : Measure ℝ) (_ : IsProbabilityMeasure ρ) (_ : IsProbabilityMeasure σ),
      MemLp id 2 ρ ∧ MemLp id 2 σ ∧ (∀ t, mLPM 2 t ρ ≤ mLPM 2 t σ) ∧
      (∃ t, mLPM 2 t ρ < mLPM 2 t σ) ∧ (∃ t, cdf ρ t < cdf σ t) ∧ ∃ t, cdf σ t < cdf ρ t) ∧
    ∃ u : ℝ → ℝ, (∀ x ∈ Icc (0 : ℝ) 1,
      0 < deriv u x ∧ deriv (deriv u) x < 0 ∧ 0 < deriv (deriv (deriv u)) x) ∧
      ¬ AntitoneOn (fun x => -deriv (deriv u) x / deriv u x) (Icc (0 : ℝ) 1) := by
  refine ⟨?_, mean_order hμ hν, ⟨fun h _ hu => utility_of_moment_order hμ hν h
    (mean_order hμ hν h) hu, moment_of_utility_order hμ hν⟩, finite_crossing, prudent_not_dara⟩
  apply forall_congr'
  intro t
  rw [mLPM_two_eq_iterated_integral_cdf t (integrable_lowerDev hμ t),
    mLPM_two_eq_iterated_integral_cdf t (integrable_lowerDev hν t)]
  exact (mul_le_mul_iff_right₀ (by norm_num : (0 : ℝ) < 2)).symm

theorem claim_C333 (μ ν : Measure ℝ) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    ((∀ t, mLPM 0 t μ ≤ mLPM 0 t μ) ∧ ¬ ∃ t, cdf μ t < cdf μ t) ∧
    (((∀ t, cdf μ t ≤ cdf ν t) ∧ ∃ t, cdf μ t < cdf ν t) ↔
      (∀ t, mLPM 0 t μ ≤ mLPM 0 t ν) ∧ ∃ t, mLPM 0 t μ < mLPM 0 t ν) := by
  simp [mLPM_zero_eq_cdf]

theorem claim_C388 (x y a b : ℝ) :
    let A := max (x - a) 0
    let B := max (a - x) 0
    let C := max (y - b) 0
    let D := max (b - y) 0
    (x - a) * (y - b) = (A - B) * (C - D) ∧
      (A - B) * (C - D) = A * C + B * D - A * D - B * C := by
  dsimp only
  constructor
  · rw [← upperDev_one, ← lowerDev_one, ← upperDev_one, ← lowerDev_one,
      upperDev_sub_lowerDev_one, upperDev_sub_lowerDev_one]
  · ring

theorem claim_C390 :
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ], MemLp id 2 μ →
      variance id μ = mUPM 2 (∫ x, x ∂μ) μ + mLPM 2 (∫ x, x ∂μ) μ) ∧
    (∀ (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ],
      Integrable Prod.fst μ → Integrable Prod.snd μ →
      Integrable (fun p : ℝ × ℝ => (p.1 - ∫ q, q.1 ∂μ) * (p.2 - ∫ q, q.2 ∂μ)) μ →
      let a := ∫ p, p.1 ∂μ
      let b := ∫ p, p.2 ∂μ
      cov[Prod.fst, Prod.snd; μ] =
        mCoUPM 1 1 a b μ + mCoLPM 1 1 a b μ - mDLPM 1 1 a b μ - mDUPM 1 1 a b μ) :=
  ⟨fun _ _ h => variance_eq_mUPM_two_add_mLPM_two h,
   fun _ _ _ _ h => (covariance_decomposition_measure h).symm⟩

theorem claim_C407 (r s a b x y : ℝ) (hr : 0 < r) (hs : 0 < s) (h : x = a ∨ y = b) :
    rLowerDev r a x * rLowerDev s b y = 0 ∧ rUpperDev r a x * rUpperDev s b y = 0 ∧
    rUpperDev r a x * rLowerDev s b y = 0 ∧ rLowerDev r a x * rUpperDev s b y = 0 := by
  rcases h with rfl | rfl <;> simp [rLowerDev, rUpperDev, hr.ne', hs.ne']

namespace A1_Basics

def signedCoMoment (x y : Fin n → ℝ) : ℝ :=
  CoUPM 1 1 (mean x) (mean y) x y + CoLPM 1 1 (mean x) (mean y) x y -
    DLPM 1 1 (mean x) (mean y) x y - DUPM 1 1 (mean x) (mean y) x y

end A1_Basics

theorem claim_C455 (x y : Fin n → ℝ) :
    popVariance x = LPM 2 (mean x) x + UPM 2 (mean x) x ∧
      covariance x y = signedCoMoment x y :=
  ⟨popVariance_eq_PM x, (covariance_decomposition x y).symm⟩

theorem claim_C459 (hn : 0 < n) (x : Fin n → ℝ) : mean x = UPM 1 0 x - LPM 1 0 x := by
  simpa using (UPM_one_sub_LPM_one hn 0 x).symm

theorem claim_C461 (hn : 1 < n) (x : Fin n → ℝ) :
    popVariance x = UPM 2 (mean x) x + LPM 2 (mean x) x ∧
    popVariance x = sampleVariance x * (((n : ℝ) - 1) / n) := by
  refine ⟨by simpa [add_comm] using popVariance_eq_PM x, ?_⟩
  have h : (1 : ℝ) < n := by exact_mod_cast hn
  rw [sampleVariance]
  field_simp

theorem claim_C462 (x : Fin n → ℝ) :
    popVariance x = covariance x x ∧ covariance x x = signedCoMoment x x :=
  ⟨(GapProjection.sample_covariance_self x).symm, (covariance_decomposition x x).symm⟩

theorem claim_C466 (hn : 1 < n) (x y : Fin n → ℝ) :
    sampleCov x y = signedCoMoment x y * ((n : ℝ) / (n - 1)) := by
  rw [signedCoMoment, covariance_decomposition]
  unfold sampleCov covariance coMoment
  have h : (1 : ℝ) < n := by exact_mod_cast hn
  field_simp

theorem claim_C469 (hn : 1 < n) (x y : Fin n → ℝ) :
    let k := (n : ℝ) / (n - 1)
    (0 < popVariance x ∧ 0 < popVariance y) →
    GapProjection.sampleCorr x y = (signedCoMoment x y * k) /
      (Real.sqrt ((UPM 2 (mean x) x + LPM 2 (mean x) x) * k) *
       Real.sqrt ((UPM 2 (mean y) y + LPM 2 (mean y) y) * k)) := by
  dsimp only
  intro _
  have hn' : (1 : ℝ) < n := by exact_mod_cast hn
  have hk : 0 < (n : ℝ) / (n - 1) := div_pos (by linarith) (by linarith)
  rw [signedCoMoment, covariance_decomposition, ← (claim_C461 hn x).1, ← (claim_C461 hn y).1,
    Real.sqrt_mul (popVariance_nonneg x), Real.sqrt_mul (popVariance_nonneg y),
    mul_mul_mul_comm, Real.mul_self_sqrt hk.le, mul_div_mul_right _ _ hk.ne',
    GapProjection.sampleCorr]

theorem claim_C470 (hn : 0 < n) (x : Fin n → ℝ) :
    (∀ t, LPM 0 t x = cdf (empMeasure x) t) ∧
    (∀ {k : ℕ} (ts : Fin k → ℝ), (fun j => rLPMVecExact 0 (ts j) x) =
      fun j => cdf (empMeasure x) (ts j)) ∧
    ∀ j, (cdfTable 0 x j).2 = cdf (empMeasure x) ((cdfTable 0 x j).1) :=
  ⟨fun t => (cdf_empMeasure hn x t).symm,
   fun ts => funext fun j => by
     rw [rLPMVecExact_eq_rLPM, rLPM_zero, cdf_empMeasure hn],
   fun j => (cdf_empMeasure hn x (x (Tuple.sort x j))).symm⟩

theorem claim_C484 (x : Fin n → ℝ) :
    popVariance x = UPM 2 (mean x) x + LPM 2 (mean x) x := by
  simpa [add_comm] using popVariance_eq_PM x

theorem claim_C521 (x : Fin n → ℝ) (t : ℝ) :
    LPM 0 t x = (#{i | x i ≤ t} : ℝ) / n ∧ LPM 1 t x = (∑ i, max (t - x i) 0) / n :=
  ⟨LPM_zero_eq_ecdf t x, by simp only [LPM, lowerDev_one]⟩

theorem claim_C523 (hn : 0 < n) (x : Fin n → ℝ) (t : ℝ) :
    LPM 0 t x + UPM 0 t x = 1 := LPM_zero_add_UPM_zero hn t x

theorem claim_C524 (x : Fin n → ℝ) :
    LPM 2 (mean x) x + UPM 2 (mean x) x = popVariance x := (popVariance_eq_PM x).symm

theorem claim_C569 (hn : 0 < n) (x y : Fin n → ℝ) (a b : ℝ) :
    (∀ i, ∃! q : Bool × Bool, (x i, y i) ∈ tail q.1 a ×ˢ tail q.2 b) ∧
    CoLPM 0 0 a b x y + CoUPM 0 0 a b x y + DLPM 0 0 a b x y + DUPM 0 0 a b x y = 1 ∧
    CoLPM 0 0 a b x y = (#{i | x i ≤ a ∧ y i ≤ b} : ℝ) / n ∧
    CoUPM 0 0 a b x y = (#{i | a < x i ∧ b < y i} : ℝ) / n ∧
    DLPM 0 0 a b x y = (#{i | a < x i ∧ y i ≤ b} : ℝ) / n ∧
    DUPM 0 0 a b x y = (#{i | x i ≤ a ∧ b < y i} : ℝ) / n :=
  ⟨fun i => partition a b (x i, y i), quadrants_sum_eq_one hn a b x y,
    GapCoMoments.CoLPM_zero_eq_card .., GapCoMoments.CoUPM_zero_eq_card ..,
    GapCoMoments.DLPM_zero_eq_card .., GapCoMoments.DUPM_zero_eq_card ..⟩

namespace A1_Basics

def ssdChecks (x : Fin n → ℝ) (y : Fin m → ℝ) : Prop :=
  mean y ≤ mean x ∧ (∀ t, IsKnot x y t → LPM 1 t x ≤ LPM 1 t y) ∧
    ∃ t, IsKnot x y t ∧ LPM 1 t x < LPM 1 t y

def tsdChecks (x : Fin n → ℝ) (y : Fin m → ℝ) : Prop :=
  mean y ≤ mean x ∧ (∀ t, IsKnot x y t → LPM 2 t x ≤ LPM 2 t y) ∧
  (∀ u w, IsKnot x y u → IsKnot x y w → (∀ i, x i ∉ Ioo u w) →
    (∀ j, y j ∉ Ioo u w) → tsdCoefA u x y < 0 →
    let v := tsdCoefB u x y / tsdCoefA u x y
    v ∈ Ioo u w → LPM 2 v x ≤ LPM 2 v y) ∧
  ∃ t, IsKnot x y t ∧ LPM 2 t x < LPM 2 t y

private theorem strict_knot {d : ℕ} {x : Fin n → ℝ} {y : Fin m → ℝ}
    (he : (∀ t, IsKnot x y t → LPM d t x = LPM d t y) → ∀ t, LPM d t x = LPM d t y)
    (h : LPMDominates d x y) : ∃ t, IsKnot x y t ∧ LPM d t x < LPM d t y := by
  by_contra! hnot
  obtain ⟨t, ht⟩ := h.2
  exact ht (he (fun u hu => le_antisymm (h.1 u) (hnot u hu)) t)

private theorem ssd_iff (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    sdDomPairCentred x y 2 false = 1 ↔ ssdChecks x y := by
  rw [sdDomPairCentred_eq hn hm x y (by decide) false (by decide), sdDomPair_ssd_iff hn hm]
  constructor
  · intro h
    exact ⟨mean_le_of_LPM_one_le hn hm h.1, fun t _ => h.1 t,
      strict_knot (LPM_one_eq_of_eq_at_knots hn hm) h⟩
  · rintro ⟨_, hk, t, _, ht⟩
    exact ⟨LPM_one_le_of_le_at_knots hm hk, t, ht.ne⟩

private theorem tsd_iff (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    sdDomPairCentred x y 3 false = 1 ↔ tsdChecks x y := by
  rw [sdDomPairCentred_eq hn hm x y (by decide) false (by decide), sdDomPair_tsd_iff hn hm]
  constructor
  · intro h
    exact ⟨mean_le_of_LPM_two_le hn hm h.1, fun t _ => h.1 t,
      fun _ _ _ _ _ _ _ _ => h.1 _, strict_knot (LPM_two_eq_of_eq_at_knots hn hm) h⟩
  · rintro ⟨hme, hk, hv, t, _, ht⟩
    exact ⟨LPM_two_le_of_le_at_knots hm hme hk hv, t, ht.ne⟩

end A1_Basics

theorem claim_C595 (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    (sdDomPairCentred x y 2 false = 1 ↔ ssdChecks x y) ∧
    (ssdChecks x y ↔ LPMDominates 1 x y) ∧
    (sdDomPairCentred y x 2 false = 1 ↔ ssdChecks y x) ∧
    (ssdChecks y x ↔ LPMDominates 1 y x) := by
  have exactTest {k l : ℕ} (hk : 0 < k) (hl : 0 < l) (v : Fin k → ℝ) (w : Fin l → ℝ) :
      ssdChecks v w ↔ LPMDominates 1 v w := by
    rw [← ssd_iff hk hl, sdDomPairCentred_eq hk hl v w (by decide) false (by decide)]
    exact sdDomPair_ssd_iff hk hl v w false
  exact ⟨ssd_iff hn hm x y, exactTest hn hm x y, ssd_iff hm hn y x, exactTest hm hn y x⟩

theorem claim_C596 (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    (sdDomPairCentred x y 3 false = 1 ↔ tsdChecks x y) ∧
    (tsdChecks x y ↔ LPMDominates 2 x y) ∧
    (sdDomPairCentred y x 3 false = 1 ↔ tsdChecks y x) ∧
    (tsdChecks y x ↔ LPMDominates 2 y x) := by
  have exactTest {k l : ℕ} (hk : 0 < k) (hl : 0 < l) (v : Fin k → ℝ) (w : Fin l → ℝ) :
      tsdChecks v w ↔ LPMDominates 2 v w := by
    rw [← tsd_iff hk hl, sdDomPairCentred_eq hk hl v w (by decide) false (by decide)]
    exact sdDomPair_tsd_iff hk hl v w false
  exact ⟨tsd_iff hn hm x y, exactTest hn hm x y, tsd_iff hm hn y x, exactTest hm hn y x⟩

theorem claim_C623 (hn : 0 < n) (x y : Fin n → ℝ) (a b : ℝ) :
    covariance x y = signedCoMoment x y ∧
    CoLPM 0 0 a b x y + CoUPM 0 0 a b x y + DLPM 0 0 a b x y + DUPM 0 0 a b x y = 1 :=
  ⟨(covariance_decomposition x y).symm, quadrants_sum_eq_one hn a b x y⟩

end

end NNS.ClaimsA
