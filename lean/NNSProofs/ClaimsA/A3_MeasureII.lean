-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import Mathlib.Analysis.Convex.Approximation
public import Mathlib.Analysis.Convex.Continuous
public import NNSProofs.GapClaims

namespace NNS.ClaimsA

open MeasureTheory ProbabilityTheory Set Filter Topology
open scoped ENNReal

@[expose] public noncomputable section

variable {μ ν : Measure ℝ} {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} {X : Ω → ℝ}
variable {n k : ℕ}

theorem claim_C142 (hX : Measurable X) (r : ℕ) (t : ℝ) :
    mUPM r t (P.map X) = ∫ w, (if r = 0 then (if t < X w then 1 else 0)
      else max (X w - t) 0 ^ r) ∂P ∧
    mLPM r t (P.map X) = ∫ w, (if r = 0 then (if X w ≤ t then 1 else 0)
      else max (t - X w) 0 ^ r) ∂P := by
  rw [mUPM, mLPM, integral_map hX.aemeasurable (upperDev_measurable r t).aestronglyMeasurable,
    integral_map hX.aemeasurable (lowerDev_measurable r t).aestronglyMeasurable]
  by_cases hr : r = 0
  · simp [hr, upperDev, lowerDev]
  · simp [hr, upperDev_eq_max_pow hr, lowerDev_eq_max_pow hr]

theorem claim_C143 [IsProbabilityMeasure μ] (r : ℕ) (t : ℝ) (h : MemLp id r μ) :
    Measurable (lowerDev r t) ∧ Measurable (upperDev r t) ∧
    (∀ x, 0 ≤ lowerDev r t x ∧ 0 ≤ upperDev r t x) ∧
    Integrable (lowerDev r t) μ ∧ Integrable (upperDev r t) μ ∧
    mLPM r t μ = ∫ x, lowerDev r t x ∂μ ∧ mUPM r t μ = ∫ x, upperDev r t x ∂μ :=
  ⟨lowerDev_measurable .., upperDev_measurable ..,
    fun x => ⟨lowerDev_nonneg r t x, upperDev_nonneg r t x⟩,
    integrable_lowerDev h t, integrable_upperDev h t, rfl, rfl⟩

theorem claim_C144 (hX : Measurable X) (r : ℕ) (t : ℝ) :
    mLPM r t (P.map X) = ∫ w in {w | X w ≤ t}, (t - X w) ^ r ∂P ∧
    mUPM r t (P.map X) = ∫ w in {w | t < X w}, (X w - t) ^ r ∂P := by
  rw [mLPM_eq_setIntegral, mUPM_eq_setIntegral]
  exact ⟨setIntegral_map measurableSet_Iic (by fun_prop) hX.aemeasurable,
    setIntegral_map measurableSet_Ioi (by fun_prop) hX.aemeasurable⟩

theorem claim_C146 [IsProbabilityMeasure μ] (r : ℕ) (h : MemLp id r μ) (t : ℝ) :
    ∫ x, (x - t) ^ r ∂μ = mUPM r t μ + (-1) ^ r * mLPM r t μ :=
  (mUPM_add_neg_one_pow_mul_mLPM h t).symm

theorem claim_C147 : ∃ μ : Measure ℝ, IsProbabilityMeasure μ ∧ Integrable id μ ∧
    eUPM 3 0 μ = ∞ ∧ eLPM 3 0 μ = ∞ ∧ ¬ Integrable (fun x : ℝ => x ^ 3) μ := by
  obtain ⟨μ, hp, hi, h⟩ := exists_integrable_eUPM_eLPM_eq_top
  obtain ⟨hu, hl, hs, -⟩ := h 3 (by decide)
  exact ⟨μ, hp, hi, hu, hl, hs⟩

theorem claim_C148 [IsProbabilityMeasure μ] (h : Integrable id μ) :
    ∫ x, x ∂μ = mUPM 1 0 μ - mLPM 1 0 μ := mean_eq_mUPM_one_sub_mLPM_one h

theorem claim_C149 [IsProbabilityMeasure μ] (h : MemLp id 2 μ) :
    variance id μ = mUPM 2 (∫ x, x ∂μ) μ + mLPM 2 (∫ x, x ∂μ) μ :=
  variance_eq_mUPM_two_add_mLPM_two h

theorem claim_C155 [IsProbabilityMeasure μ] (t : ℝ) : 1 - cdf μ t = μ.real (Ioi t) := by
  have h := mLPM_zero_add_mUPM_zero (μ := μ) t
  rw [mLPM_zero_eq_cdf, mUPM_zero_eq] at h
  linarith

theorem claim_C156 [IsProbabilityMeasure μ] (t : ℝ) :
    mLPM 0 t μ = μ.real (Iic t) ∧ mUPM 0 t μ = μ.real (Ioi t) :=
  ⟨(mLPM_zero_eq_cdf t).trans (cdf_eq_real μ t), mUPM_zero_eq t⟩

private theorem signedMomentMap [IsProbabilityMeasure P] (hX : Measurable X) (r : ℕ)
    (h : MemLp X r P) (t : ℝ) :
    ∫ w, (X w - t) ^ r ∂P = mUPM r t (P.map X) + (-1) ^ r * mLPM r t (P.map X) := by
  have := (Measure.isProbabilityMeasure_map_iff hX.aemeasurable).mpr ‹IsProbabilityMeasure P›
  have hm : MemLp id r (P.map X) :=
    (memLp_map_measure_iff aestronglyMeasurable_id hX.aemeasurable).mpr h
  rw [← claim_C146 r hm, integral_map hX.aemeasurable (by fun_prop)]

theorem claim_C157 [IsProbabilityMeasure P] (hX : Measurable X) (t : ℝ) :
    MeasurableSet {w | X w ≤ t} ∧ MeasurableSet {w | t < X w} ∧
    {w | X w ≤ t}ᶜ = {w | t < X w} ∧
    (∀ g : Ω → ℝ, Integrable g P → ∫ w, g w ∂P =
      (∫ w in {w | X w ≤ t}, g w ∂P) + ∫ w in {w | t < X w}, g w ∂P) ∧
    (∀ r : ℕ, mLPM r t (P.map X) = ∫ w in {w | X w ≤ t}, (t - X w) ^ r ∂P ∧
      mUPM r t (P.map X) = ∫ w in {w | t < X w}, (X w - t) ^ r ∂P) ∧
    ∀ r : ℕ, MemLp X r P → ∫ w, (X w - t) ^ r ∂P =
      mUPM r t (P.map X) + (-1) ^ r * mLPM r t (P.map X) := by
  have hc : {w | X w ≤ t}ᶜ = {w | t < X w} := by ext; simp
  refine ⟨hX measurableSet_Iic, hX measurableSet_Ioi, hc, fun g hg => ?_,
    fun r => claim_C144 hX r t, fun r h => signedMomentMap hX r h t⟩
  rw [← hc]
  exact (integral_add_compl (s := {w | X w ≤ t}) (hX measurableSet_Iic) hg).symm

theorem claim_C158 [IsProbabilityMeasure P] (hX : Measurable X) (r : ℕ) (t : ℝ) :
    (∀ w, X w - t = max (X w - t) 0 - max (t - X w) 0) ∧
    mLPM r t (P.map X) = ∫ w in {w | X w ≤ t}, (t - X w) ^ r ∂P ∧
    mUPM r t (P.map X) = ∫ w in {w | t < X w}, (X w - t) ^ r ∂P ∧
    (MemLp X r P → ∫ w, (X w - t) ^ r ∂P =
      mUPM r t (P.map X) + (-1) ^ r * mLPM r t (P.map X)) := by
  obtain ⟨hl, hu⟩ := claim_C144 (P := P) hX r t
  exact ⟨fun w => by simpa only [neg_sub] using (max_sub_max_neg (X w - t)).symm,
    hl, hu, fun h => signedMomentMap hX r h t⟩

theorem claim_C167 [IsProbabilityMeasure μ] (h : MemLp id 2 μ) :
    variance id μ = mUPM 2 (∫ x, x ∂μ) μ + mLPM 2 (∫ x, x ∂μ) μ := claim_C149 h

theorem claim_C184 [IsProbabilityMeasure μ] (h : MemLp id 4 μ) :
    ∫ x, (x - ∫ y, y ∂μ) ^ 4 ∂μ = mUPM 4 (∫ y, y ∂μ) μ + mLPM 4 (∫ y, y ∂μ) μ := by
  have hs := claim_C146 4 h (∫ y, y ∂μ)
  norm_num at hs
  exact hs

private theorem upperIntegrableOfLE [IsProbabilityMeasure μ] {r j : ℕ} (hj : j ≤ r) (t : ℝ)
    (h : Integrable (upperDev r t) μ) : Integrable (upperDev j t) μ := by
  apply (h.add (integrable_const 1)).mono' (upperDev_measurable j t).aestronglyMeasurable
  apply Eventually.of_forall
  intro x
  rw [Real.norm_of_nonneg (upperDev_nonneg ..)]
  change upperDev j t x ≤ upperDev r t x + 1
  by_cases hx : t < x
  · simp only [upperDev, ite_eq_left hx]
    rcases le_total (x - t) 1 with hle | hge
    · linarith [pow_le_one₀ (sub_nonneg.mpr hx.le) hle (n := j),
        pow_nonneg (sub_nonneg.mpr hx.le) r]
    · linarith [pow_le_pow_right₀ hge hj]
  · simp [upperDev, hx]

theorem claim_C188 [IsProbabilityMeasure μ] (r : ℕ) (t : ℝ) :
    (∀ j ≤ r, Integrable (lowerDev j t) μ ∧ Integrable (upperDev j t) μ) ↔
      eLPM r t μ < ∞ ∧ eUPM r t μ < ∞ := by
  rw [eLPM_lt_top_iff, eUPM_lt_top_iff]
  exact ⟨fun h => h r le_rfl, fun ⟨hl, hu⟩ j hj =>
    ⟨integrable_lowerDev_of_le hj le_rfl hl, upperIntegrableOfLE hj t hu⟩⟩

theorem claim_C192 [IsProbabilityMeasure μ] (t : ℝ) : mLPM 0 t μ = cdf μ t :=
  mLPM_zero_eq_cdf t

theorem claim_C194 [IsProbabilityMeasure μ] (t : ℝ) : cdf μ t = mLPM 0 t μ :=
  (claim_C192 t).symm

theorem claim_C199 [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] (hn : 0 < n) (x : Fin n → ℝ) :
    ((∀ t, mLPM 0 t μ = mLPM 0 t ν) ↔ μ = ν) ∧
    ((∀ t, mLPM 0 t μ = LPM 0 t x) ↔ μ = empMeasure x) := by
  have := isProbabilityMeasure_empMeasure hn x
  have extLaw (ρ : Measure ℝ) [IsProbabilityMeasure ρ] :
      (∀ t, mLPM 0 t μ = mLPM 0 t ρ) ↔ μ = ρ := by
    constructor
    · intro h
      apply Measure.eq_of_cdf μ ρ
      exact StieltjesFunction.ext fun t => by simpa [mLPM_zero_eq_cdf] using h t
    · rintro rfl t
      rfl
  exact ⟨extLaw ν, by simpa only [mLPM_empMeasure] using extLaw (empMeasure x)⟩

theorem claim_C225 [IsProbabilityMeasure μ] (t : ℝ) :
    cdf μ t = mLPM 0 t μ ∧ mLPM 0 t μ = μ.real (Iic t) := ⟨claim_C194 t, (claim_C156 t).1⟩

theorem claim_C291 [IsProbabilityMeasure μ] (t : ℝ) :
    mLPM 0 t μ = μ.real (Iic t) ∧
    mLPM 1 t μ = ∫ x, max (t - x) 0 ∂μ ∧
    mLPM 2 t μ = ∫ x, max (t - x) 0 ^ 2 ∂μ := by
  exact ⟨(claim_C156 t).1, by simp [mLPM, lowerDev_one],
    by simp [mLPM, lowerDev_eq_max_pow (by decide : 2 ≠ 0)]⟩

open scoped Finset in
theorem claim_C299 {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j) :
    Continuous (fun t => cdfCurve 1 t x) ∧
    ∀ t, Tendsto (fun s => cdfCurve 0 s x) (𝓝[<] t)
      (𝓝 ((Finset.univ.filter fun i => x i < t).card / (n : ℝ))) ∧
      cdfCurve 0 t x = ((Finset.univ.filter fun i => x i < t).card : ℝ) / n +
        ((Finset.univ.filter fun i => x i = t).card : ℝ) / n :=
  ⟨cdfCurve_continuous (by decide) hx, fun t =>
    ⟨cdfCurve_zero_tendsto_left t x, cdfCurve_zero_eq_lt_add_eq t x⟩⟩

section Dominance

variable [IsProbabilityMeasure μ] [IsProbabilityMeasure ν]

private def cappedOrder (μ ν : Measure ℝ) (f : ℝ → ℝ) : Prop :=
  ∀ t : ℝ, Integrable (fun x => f (min x t)) μ ∧ Integrable (fun x => f (min x t)) ν ∧
    (∫ x, f (min x t) ∂ν) ≤ ∫ x, f (min x t) ∂μ

private theorem cappedConst (b : ℝ) : cappedOrder μ ν (fun _ => b) :=
  fun _ => ⟨integrable_const _, integrable_const _, by simp⟩

omit [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] in
private theorem cappedAddLinear {f : ℝ → ℝ} (hid : cappedOrder μ ν id) (hf : cappedOrder μ ν f)
    {a : ℝ} (ha : 0 ≤ a) : cappedOrder μ ν (fun x => a * x + f x) := by
  intro t
  obtain ⟨hm, hn, h⟩ := hid t
  obtain ⟨fm, fn_, hf'⟩ := hf t
  refine ⟨(hm.const_mul a).add fm, (hn.const_mul a).add fn_, ?_⟩
  simp only [id_eq] at hm hn h
  rw [integral_add (hn.const_mul a) fn_, integral_add (hm.const_mul a) fm,
    integral_const_mul, integral_const_mul]
  exact add_le_add (mul_le_mul_of_nonneg_left h ha) hf'

private theorem cappedMin {f : ℝ → ℝ} (hf : cappedOrder μ ν f) (hm : Monotone f)
    (hc : Continuous f) (b : ℝ) : cappedOrder μ ν (fun x => min (f x) b) := by
  by_cases h₁ : ∀ x, f x ≤ b
  · simpa only [min_eq_left (h₁ _)] using hf
  by_cases h₂ : ∀ x, b ≤ f x
  · simpa only [min_eq_right (h₂ _)] using (cappedConst (μ := μ) (ν := ν) b)
  obtain ⟨x, hx⟩ := not_forall.mp h₂
  obtain ⟨y, hy⟩ := not_forall.mp h₁
  obtain ⟨z, rfl⟩ := intermediate_value_univ x y hc ⟨(lt_of_not_ge hx).le, (lt_of_not_ge hy).le⟩
  intro t
  simpa only [← hm.map_min, min_assoc] using hf (min t z)

private def envelope (s : Finset ℕ) (hs : s.Nonempty) (a b : ℕ → ℝ) (x : ℝ) : ℝ :=
  s.inf' hs fun i => a i * x + b i

private theorem envelopeOrder (hid : cappedOrder μ ν id) (s : Finset ℕ) (hs : s.Nonempty)
    (a b : ℕ → ℝ) (ha : ∀ i ∈ s, 0 ≤ a i) : cappedOrder μ ν (envelope s hs a b) := by
  classical
  -- Removing the smallest slope leaves a nondecreasing envelope capped at its intercept.
  obtain ⟨i, hi, hmin⟩ := s.exists_min_image a hs
  by_cases he : (s.erase i).Nonempty
  · let g := envelope (s.erase i) he (fun j => a j - a i) b
    have hsmall : cappedOrder μ ν g := envelopeOrder hid _ he _ b
      fun j hj => sub_nonneg.mpr (hmin j (Finset.mem_of_mem_erase hj))
    have hg : Monotone g := by
      intro x y hxy
      apply Finset.le_inf'
      intro j hj
      apply (Finset.inf'_le _ hj).trans
      exact add_le_add_left (mul_le_mul_of_nonneg_left hxy
        (sub_nonneg.mpr (hmin j (Finset.mem_of_mem_erase hj)))) _
    have hc : Continuous g := Continuous.finset_inf'_apply he fun _ _ => by fun_prop
    have heq : envelope s hs a b = fun x => a i * x + min (g x) (b i) := by
      funext x
      obtain ⟨j, hj, hjx⟩ := Finset.exists_mem_eq_inf' he (fun j => (a j - a i) * x + b j)
      have hj' : g x = (a j - a i) * x + b j := hjx
      apply le_antisymm
      · rw [add_min]
        apply le_min
        · rw [hj']
          unfold envelope
          linarith [Finset.inf'_le (fun k => a k * x + b k) (Finset.mem_of_mem_erase hj)]
        · exact Finset.inf'_le _ hi
      · apply Finset.le_inf'
        intro j hj
        by_cases hij : j = i
        · subst j
          linarith [min_le_right (g x) (b i)]
        · have hle : g x ≤ (a j - a i) * x + b j :=
            Finset.inf'_le _ (Finset.mem_erase.mpr ⟨hij, hj⟩)
          linarith [min_le_left (g x) (b i)]
    rw [heq]
    exact cappedAddLinear hid (cappedMin hsmall hg hc _) (ha i hi)
  · have heq : s = {i} := by
      have := Finset.insert_erase hi
      simpa [Finset.not_nonempty_iff_eq_empty.mp he] using this.symm
    subst s
    have hfun : envelope {i} hs a b = fun x => a i * x + b i := by funext; simp [envelope]
    rw [hfun]
    exact cappedAddLinear hid (cappedConst (b i)) (ha i hi)
termination_by s.card
decreasing_by exact Finset.card_erase_lt_of_mem hi

private theorem integralUtility (hid : cappedOrder μ ν id) (hμ : Integrable id μ)
    (hν : Integrable id ν) {u : ℝ → ℝ} (hm : Monotone u) (hc : ConcaveOn ℝ univ u)
    (huμ : Integrable u μ) (huν : Integrable u ν) : ∫ x, u x ∂ν ≤ ∫ x, u x ∂μ := by
  have hcont : Continuous u := continuousOn_univ.mp (hc.continuousOn isOpen_univ)
  obtain ⟨l, c, hle, hsup⟩ := hc.neg.real_univ_sSup_of_nat_affine_eq hcont.neg.lowerSemicontinuous
  let a (i : ℕ) := -(l i 1)
  let b (i : ℕ) := -(c i)
  have hl (i : ℕ) (x : ℝ) : l i x = l i 1 * x := by
    calc l i x = l i (x • (1 : ℝ)) := by simp
         _ = l i 1 * x := by rw [map_smul, smul_eq_mul, mul_comm]
  have hab (i : ℕ) (x : ℝ) : u x ≤ a i * x + b i := by
    have h := hle i x
    change l i x + c i ≤ -u x at h
    rw [hl] at h
    dsimp [a, b]
    linarith
  have ha (i : ℕ) : 0 ≤ a i := by
    by_contra! h
    let y := (u 0 - b i - 1) / a i
    have hxy := mul_le_mul_of_nonpos_left (le_max_right (0 : ℝ) y) h.le
    have hdiv : a i * y = u 0 - b i - 1 := by dsimp [y]; field_simp [h.ne]
    linarith [hm (le_max_left 0 y), hab i (max 0 y)]
  let s (m : ℕ) := Finset.range (m + 1)
  have hmem (m : ℕ) : 0 ∈ s m := Finset.mem_range.mpr (by omega)
  let F (m : ℕ) := envelope (s m) ⟨0, hmem m⟩ a b
  have hlow (m : ℕ) (x : ℝ) : u x ≤ F m x := Finset.le_inf' _ _ fun i _ => hab i x
  have hupp (m i : ℕ) (hi : i ≤ m) (x : ℝ) : F m x ≤ a i * x + b i :=
    Finset.inf'_le _ (Finset.mem_range.mpr (by omega : i < m + 1))
  have hF (m : ℕ) : Continuous (F m) := Continuous.finset_inf'_apply _ fun _ _ => by fun_prop
  have hlim (x : ℝ) : Tendsto (fun m : ℕ => F m (min x (m : ℝ))) atTop (𝓝 (u x)) := by
    obtain ⟨N, hN⟩ := exists_nat_ge x
    have he : ∀ᶠ m : ℕ in atTop, min x (m : ℝ) = x := by
      filter_upwards [eventually_ge_atTop N] with m hm'
      exact min_eq_left (hN.trans (Nat.cast_le.mpr hm'))
    apply (tendsto_congr' (he.mono fun m h => congrArg (F m) h)).mpr
    apply tendsto_order.mpr
    refine ⟨fun z hz => Eventually.of_forall fun m => hz.trans_le (hlow m x), fun z hz => ?_⟩
    have hex : ∃ i, a i * x + b i < z := by
      by_contra! h
      have hsx := congrFun hsup x
      simp only [iSup_apply, Pi.add_apply, Function.const_apply, Pi.neg_apply] at hsx
      have hb : (⨆ i, l i x + c i) ≤ -z := by
        apply ciSup_le
        intro i
        have hh := h i
        rw [hl]
        dsimp [a, b] at hh
        linarith
      linarith
    obtain ⟨i, hi⟩ := hex
    filter_upwards [eventually_ge_atTop i] with m hm'
    exact (hupp m i hm' x).trans_lt hi
  -- The first line bounds every envelope above, and the utility bounds it below.
  let G (m : ℕ) (x : ℝ) := F m (min x (m : ℝ))
  have lim (ρ : Measure ℝ) [IsProbabilityMeasure ρ] (hi : Integrable id ρ) (hu : Integrable u ρ) :
      Tendsto (fun m => ∫ x, G m x ∂ρ) atTop (𝓝 (∫ x, u x ∂ρ)) := by
    apply tendsto_integral_of_dominated_convergence
      (fun x => |u x| + |u 0| + (a 0 * |x| + |b 0|))
    · intro m
      exact ((hF m).comp (continuous_id.min continuous_const)).aestronglyMeasurable
    · exact (hu.abs.add (integrable_const _)).add ((hi.abs.const_mul _).add (integrable_const _))
    · intro m
      apply Eventually.of_forall
      intro x
      rw [Real.norm_eq_abs, abs_le]
      have h₀ := hlow m (min x (m : ℝ))
      have h₁ := hupp m 0 (Nat.zero_le m) (min x (m : ℝ))
      have h₂ := mul_le_mul_of_nonneg_left ((min_le_left x (m : ℝ)).trans (le_abs_self x)) (ha 0)
      have h₃ := mul_nonneg (ha 0) (abs_nonneg x)
      change _ ≤ F m (min x (m : ℝ)) ∧ F m (min x (m : ℝ)) ≤ _
      rw [hm.map_min] at h₀
      have hm₀ := hm (Nat.cast_nonneg m : (0 : ℝ) ≤ m)
      constructor
      · rcases min_choice (u x) (u m) with hh | hh <;> rw [hh] at h₀ <;>
          linarith [neg_abs_le (u x), neg_abs_le (u 0), abs_nonneg (u x),
            abs_nonneg (u 0), abs_nonneg (b 0)]
      · linarith [le_abs_self (b 0), abs_nonneg (u x), abs_nonneg (u 0)]
    · exact Eventually.of_forall hlim
  apply le_of_tendsto_of_tendsto' (lim ν hν huν) (lim μ hμ huμ)
  intro m
  exact (envelopeOrder hid (s m) ⟨0, hmem m⟩ a b (fun i _ => ha i) m).2.2

private lemma capUtilityLimit {u : ℝ → ℝ} (hm : Monotone u) (hc : Continuous u)
    (ρ : Measure ℝ) [IsProbabilityMeasure ρ] (hu : Integrable u ρ) :
    Tendsto (fun m : ℕ => ∫ x, u (min x m) ∂ρ) atTop (𝓝 (∫ x, u x ∂ρ)) := by
  apply tendsto_integral_of_dominated_convergence (fun x => |u x| + |u 0|)
  · intro m
    exact (hc.measurable.comp (measurable_id.min measurable_const)).aestronglyMeasurable
  · exact hu.abs.add (integrable_const _)
  · intro m
    apply Eventually.of_forall
    intro x
    simp only [hm.map_min, Real.norm_eq_abs, abs_le]
    have h0 := hm (Nat.cast_nonneg m : (0 : ℝ) ≤ m)
    constructor
    · exact le_min (by linarith [neg_abs_le (u x), abs_nonneg (u 0)])
        (by linarith [neg_abs_le (u 0), abs_nonneg (u x)])
    · exact (min_le_left ..).trans (by linarith [le_abs_self (u x), abs_nonneg (u 0)])
  · apply Eventually.of_forall
    intro x
    obtain ⟨N, hN⟩ := exists_nat_ge x
    apply tendsto_const_nhds.congr'
    filter_upwards [eventually_ge_atTop N] with m hm
    simp [min_eq_left (hN.trans (Nat.cast_le.mpr hm))]

private theorem utilityOrderOfCappedOrder (hid : cappedOrder μ ν id) {u : ℝ → ℝ}
    (hm : Monotone u) (hc : ConcaveOn ℝ univ u) (huμ : Integrable u μ) (huν : Integrable u ν) :
    ∫ x, u x ∂ν ≤ ∫ x, u x ∂μ := by
  have cont := continuousOn_univ.mp (hc.continuousOn isOpen_univ)
  apply le_of_tendsto_of_tendsto' (capUtilityLimit hm cont ν huν) (capUtilityLimit hm cont μ huμ)
  intro m
  let f (x : ℝ) := min x (m : ℝ)
  have hf : Measurable f := measurable_id.min measurable_const
  have prob (ρ : Measure ℝ) [IsProbabilityMeasure ρ] : IsProbabilityMeasure (ρ.map f) :=
    (Measure.isProbabilityMeasure_map_iff hf.aemeasurable).mpr inferInstance
  have intId (ρ : Measure ℝ) (hi : Integrable f ρ) : Integrable id (ρ.map f) :=
    (integrable_map_measure aestronglyMeasurable_id hf.aemeasurable).mpr hi
  have intU (ρ : Measure ℝ) [IsProbabilityMeasure ρ] (hu : Integrable u ρ) :
      Integrable u (ρ.map f) := by
    apply (integrable_map_measure cont.aestronglyMeasurable hf.aemeasurable).mpr
    have hi : Integrable (fun x => min (u x) (u m)) ρ := hu.inf (integrable_const (u m))
    simpa only [Function.comp_def, f, hm.map_min] using hi
  have capped : cappedOrder (μ.map f) (ν.map f) id := by
    intro t
    have ht : Measurable (fun x : ℝ => id (min x t)) := measurable_id.min measurable_const
    simpa only [integrable_map_measure ht.aestronglyMeasurable hf.aemeasurable,
      integral_map hf.aemeasurable ht.aestronglyMeasurable, Function.comp_def, f, min_assoc]
      using hid (min (m : ℝ) t)
  have ord := integralUtility capped (intId μ (hid m).1) (intId ν (hid m).2.1) hm hc
    (intU μ huμ) (intU ν huν)
  simpa only [integral_map hf.aemeasurable cont.aestronglyMeasurable] using ord

end Dominance

theorem claim_C331 [IsProbabilityMeasure μ] [IsProbabilityMeasure ν]
    (hμ : ∀ t, Integrable (lowerDev 1 t) μ) (hν : ∀ t, Integrable (lowerDev 1 t) ν) :
    ((∀ u : ℝ → ℝ, Monotone u → ConcaveOn ℝ univ u → Integrable u μ → Integrable u ν →
      ∫ x, u x ∂ν ≤ ∫ x, u x ∂μ) ↔
      ∀ t, ∫ s in Iic t, cdf μ s ≤ ∫ s in Iic t, cdf ν s) ∧
    (∃ ρ σ : Measure ℝ, IsProbabilityMeasure ρ ∧ IsProbabilityMeasure σ ∧
      Integrable id ρ ∧ Integrable id σ ∧
      (∀ t, ∫ s in Iic t, cdf ρ s ≤ ∫ s in Iic t, cdf σ s) ∧
      (∃ t, cdf σ t < cdf ρ t) ∧ ∃ t, cdf ρ t < cdf σ t) := by
  have hmin (t : ℝ) : (fun x : ℝ => min x t) = fun x => t - lowerDev 1 t x := by
    funext x
    rw [lowerDev_one]
    rcases le_total x t with hx | hx
    · rw [min_eq_left hx, max_eq_left (sub_nonneg.mpr hx)]; ring
    · rw [min_eq_right hx, max_eq_right (sub_nonpos.mpr hx)]; ring
  have hi (t : ℝ) (ρ : Measure ℝ) [IsProbabilityMeasure ρ] (h : Integrable (lowerDev 1 t) ρ) :
      Integrable (fun x : ℝ => min x t) ρ := hmin t ▸ (integrable_const t).sub h
  have integralMin (t : ℝ) (ρ : Measure ℝ) [IsProbabilityMeasure ρ]
      (h : Integrable (lowerDev 1 t) ρ) :
      (∫ x, min x t ∂ρ) = t - ∫ s in Iic t, cdf ρ s := by
    rw [hmin t, integral_sub (integrable_const _) h, ← mLPM_one_eq_integral_cdf t h]
    simp [mLPM]
  constructor
  · constructor
    · intro h t
      have hv := h (fun x => min x t) (monotone_id.min monotone_const)
        ((concaveOn_id convex_univ).inf (concaveOn_const t convex_univ))
        (hi t μ (hμ t)) (hi t ν (hν t))
      simpa only [integralMin t ν (hν t), integralMin t μ (hμ t), sub_le_sub_iff_left] using hv
    · intro h u hm hc huμ huν
      apply utilityOrderOfCappedOrder _ hm hc huμ huν
      intro t
      refine ⟨hi t μ (hμ t), hi t ν (hν t), ?_⟩
      change (∫ x, min x t ∂ν) ≤ ∫ x, min x t ∂μ
      rw [integralMin t ν (hν t), integralMin t μ (hμ t)]
      exact sub_le_sub_left (h t) t
  · let ρ := empMeasure (![0] : Fin 1 → ℝ)
    let σ := empMeasure (![-1, 1] : Fin 2 → ℝ)
    have hρ : IsProbabilityMeasure ρ := isProbabilityMeasure_empMeasure (by decide) _
    have hσ : IsProbabilityMeasure σ := isProbabilityMeasure_empMeasure (by decide) _
    refine ⟨ρ, σ, hρ, hσ, GapClaims.integrable_empMeasure (by decide : 0 < 1) _ id,
      GapClaims.integrable_empMeasure (by decide : 0 < 2) _ id, fun t => ?_, ⟨0, ?_⟩, ⟨-1, ?_⟩⟩
    · rw [← mLPM_one_eq_integral_cdf t (GapClaims.integrable_empMeasure (by decide : 0 < 1) _ _),
        ← mLPM_one_eq_integral_cdf t (GapClaims.integrable_empMeasure (by decide : 0 < 2) _ _)]
      simp [mLPM_empMeasure, LPM, Fin.sum_univ_two, lowerDev_one]
      grind
    · simp only [ρ, σ, cdf_empMeasure (by decide : 0 < 1), cdf_empMeasure (by decide : 0 < 2),
        LPM, Fin.sum_univ_one, Fin.sum_univ_two]
      norm_num [lowerDev]
    · simp only [ρ, σ, cdf_empMeasure (by decide : 0 < 1), cdf_empMeasure (by decide : 0 < 2),
        LPM, Fin.sum_univ_one, Fin.sum_univ_two]
      norm_num [lowerDev]

theorem claim_C334 [IsProbabilityMeasure μ] [IsProbabilityMeasure ν]
    (hμ : ∀ t, Integrable (lowerDev 1 t) μ) (hν : ∀ t, Integrable (lowerDev 1 t) ν) :
    (∀ t, mLPM 1 t μ = ∫ s in Iic t, cdf μ s) ∧
    ((∀ t, mLPM 1 t μ ≤ mLPM 1 t ν) ↔ ∀ t, ∫ s in Iic t, cdf μ s ≤ ∫ s in Iic t, cdf ν s) :=
  ⟨fun t => mLPM_one_eq_integral_cdf t (hμ t), forall_mLPM_one_le_iff hμ hν⟩

theorem claim_C335 [IsProbabilityMeasure μ] [IsProbabilityMeasure ν]
    (hμ : ∀ t, Integrable (lowerDev 2 t) μ) (hν : ∀ t, Integrable (lowerDev 2 t) ν) :
    (∀ t, mLPM 2 t μ ≤ mLPM 2 t ν) ↔
      ∀ t, ∫ s in Iic t, ∫ u in Iic s, cdf μ u ≤ ∫ s in Iic t, ∫ u in Iic s, cdf ν u :=
  forall_mLPM_two_le_iff hμ hν

theorem claim_C355 [IsProbabilityMeasure μ] {p : ℝ} (hp : p ∈ Ioo 0 1) :
    sInf {t | p ≤ mLPM 0 t μ} = GapQuantile.quantile μ p ∧
    {t | p ≤ mLPM 0 t μ} = Ici (GapQuantile.quantile μ p) := by
  simp only [mLPM_zero_eq_cdf]
  exact ⟨(GapQuantile.quantile_eq_sInf μ p).symm, GapQuantile.setOf_le_cdf_eq_Ici μ hp.1 hp.2⟩

theorem claim_C378 [IsProbabilityMeasure μ] (h : MemLp id 2 μ) :
    variance id μ = (∫ x, max (x - ∫ y, y ∂μ) 0 ^ 2 ∂μ) +
      ∫ x, max ((∫ y, y ∂μ) - x) 0 ^ 2 ∂μ := by
  simpa [mUPM, mLPM, upperDev_eq_max_pow (by decide : 2 ≠ 0),
    lowerDev_eq_max_pow (by decide : 2 ≠ 0)] using claim_C149 h

theorem claim_C383 [IsProbabilityMeasure μ] :
    (∀ t, cdf μ t = mLPM 0 t μ) ∧ ∀ (r : ℕ) (t : ℝ), MemLp id r μ →
      ∫ x, (x - t) ^ r ∂μ = mUPM r t μ + (-1) ^ r * mLPM r t μ :=
  ⟨claim_C194, fun r t h => claim_C146 r h t⟩

theorem claim_C387 (x : ℝ) : x = max x 0 - max (-x) 0 := (max_sub_max_neg x).symm

open GapOperators

theorem claim_C425 (μ : Measure ℝ) (x : Fin n → ℝ) (t : ℝ) {d : ℝ} (hd : 0 < d) :
    rlpm d t μ = ∫ v, max (t - v) 0 ^ d ∂μ ∧
    rLPM d t x = (∑ i, max (t - x i) 0 ^ d) / n ∧
    rLPM 0 t x = ((Finset.univ.filter fun i => x i ≤ t).card : ℝ) / n ∧
    rlpm 0 t μ = μ.real (Iic t) ∧
    rlpm 0 0 (Measure.dirac 1) ≠ ∫ v, max (0 - v) 0 ^ (0 : ℝ) ∂Measure.dirac (1 : ℝ) := by
  obtain ⟨h0, h1⟩ := rlpm_zero_ne_literal
  exact ⟨rlpm_eq_integral_max_rpow hd .., rLPM_eq_max hd.ne' .., rLPM_zero_eq_ecdf .., rlpm_zero ..,
    by rw [h0, h1]; norm_num⟩

theorem claim_C426 (μ : Measure ℝ) (x : Fin n → ℝ) (t : ℝ) {d : ℝ} (hd : 0 < d) :
    rupm d t μ = ∫ v, max (v - t) 0 ^ d ∂μ ∧
    rUPM d t x = (∑ i, max (x i - t) 0 ^ d) / n ∧
    rUPM 0 t x = ((Finset.univ.filter fun i => t < x i).card : ℝ) / n ∧
    rupm 0 t μ = μ.real (Ioi t) ∧
    rupm 0 0 (Measure.dirac 0) ≠ ∫ v, max (v - 0) 0 ^ (0 : ℝ) ∂Measure.dirac (0 : ℝ) :=
  ⟨rupm_eq_integral_max_rpow hd.ne' .., rUPM_eq_max hd.ne' ..,
    by simp [rUPM, rUpperDev_zero, Finset.sum_boole],
    rupm_zero ..,
    by simp [rupm, rUpperDev]⟩

private theorem regionReal (d t : ℝ) (μ : Measure ℝ) :
    rlpm d t μ = ∫ v in Iic t, (t - v) ^ d ∂μ ∧
    rupm d t μ = ∫ v in Ioi t, (v - t) ^ d ∂μ := by
  constructor
  · change ∫ v, (Iic t).indicator (fun v => (t - v) ^ d) v ∂μ = _
    exact integral_indicator measurableSet_Iic
  · change ∫ v, (Ioi t).indicator (fun v => (v - t) ^ d) v ∂μ = _
    exact integral_indicator measurableSet_Ioi

theorem claim_C481 [IsProbabilityMeasure μ] (d t : ℝ) :
    rlpm d t μ = ∫ v in Iic t, (t - v) ^ d ∂(cdf μ).measure ∧
    rupm d t μ = ∫ v in Ioi t, (v - t) ^ d ∂(cdf μ).measure ∧
    (∫⁻ v, ENNReal.ofReal (rLowerDev d t v) ∂μ) =
      ∫⁻ v in Iic t, ENNReal.ofReal ((t - v) ^ d) ∂(cdf μ).measure ∧
    (∫⁻ v, ENNReal.ofReal (rUpperDev d t v) ∂μ) =
      ∫⁻ v in Ioi t, ENNReal.ofReal ((v - t) ^ d) ∂(cdf μ).measure := by
  rw [measure_cdf]
  obtain ⟨hl, hu⟩ := regionReal d t μ
  refine ⟨hl, hu, ?_, ?_⟩
  · rw [← lintegral_indicator measurableSet_Iic]
    exact lintegral_congr fun v => by by_cases hv : v ≤ t <;> simp [rLowerDev, hv]
  · rw [← lintegral_indicator measurableSet_Ioi]
    exact lintegral_congr fun v => by by_cases hv : t < v <;> simp [rUpperDev, hv]

private theorem empiricalRegions (x : Fin n → ℝ) (d t : ℝ) :
    rLPM d t x = ∫ v in Iic t, (t - v) ^ d ∂empMeasure x ∧
    rUPM d t x = ∫ v in Ioi t, (v - t) ^ d ∂empMeasure x := by
  have hl : rlpm d t (empMeasure x) = rLPM d t x := integral_empMeasure x _
  have hu : rupm d t (empMeasure x) = rUPM d t x := integral_empMeasure x _
  simpa only [hl, hu] using regionReal d t (empMeasure x)

theorem claim_C483 (hX : Measurable X) (x : Fin n → ℝ) (d t : ℝ) (hd : 0 ≤ d) :
    MeasurableSet {w | X w ≤ t} ∧ MeasurableSet {w | t < X w} ∧
    rlpm d t (P.map X) = ∫ w in {w | X w ≤ t}, (t - X w) ^ d ∂P ∧
    rupm d t (P.map X) = ∫ w in {w | t < X w}, (X w - t) ^ d ∂P ∧
    rLPM d t x = ∫ v in Iic t, (t - v) ^ d ∂empMeasure x ∧
    rUPM d t x = ∫ v in Ioi t, (v - t) ^ d ∂empMeasure x := by
  obtain ⟨hl, hu⟩ := regionReal d t (P.map X)
  exact ⟨hX measurableSet_Iic, hX measurableSet_Ioi,
    hl.trans (setIntegral_map measurableSet_Iic (by fun_prop) hX.aemeasurable),
    hu.trans (setIntegral_map measurableSet_Ioi (by fun_prop) hX.aemeasurable), empiricalRegions ..⟩

theorem claim_C508 [IsProbabilityMeasure μ] (d t : ℝ) :
    rlpm d t μ = ∫ v in Iic t, (t - v) ^ d ∂(cdf μ).measure ∧
    rupm d t μ = ∫ v in Ioi t, (v - t) ^ d ∂(cdf μ).measure ∧
    (∫⁻ v, ENNReal.ofReal (rLowerDev d t v) ∂μ) =
      ∫⁻ v in Iic t, ENNReal.ofReal ((t - v) ^ d) ∂(cdf μ).measure ∧
    (∫⁻ v, ENNReal.ofReal (rUpperDev d t v) ∂μ) =
      ∫⁻ v in Ioi t, ENNReal.ofReal ((v - t) ^ d) ∂(cdf μ).measure := claim_C481 d t

theorem claim_C509 (x : Fin n → ℝ) (d t : ℝ) :
    rLPM d t x = ∫ v in Iic t, (t - v) ^ d ∂empMeasure x ∧
    rUPM d t x = ∫ v in Ioi t, (v - t) ^ d ∂empMeasure x := empiricalRegions ..

theorem claim_C510 [IsProbabilityMeasure μ] (h : MemLp id 2 μ) (x : Fin n → ℝ) :
    variance id μ = mUPM 2 (∫ v, v ∂μ) μ + mLPM 2 (∫ v, v ∂μ) μ ∧
    popVariance x = LPM 2 (mean x) x + UPM 2 (mean x) x := ⟨claim_C149 h, popVariance_eq_PM x⟩

namespace A3_MeasureII

def sampleRatios (d t : ℝ) (x : Fin n → ℝ) : ℝ × ℝ :=
  let total := rLPM d t x + rUPM d t x
  if total = 0 then (rLPM 0 t x, rUPM 0 t x) else (rLPM d t x / total, rUPM d t x / total)

private theorem sampleRatios_zero (hn : 0 < n) (x : Fin n → ℝ) (t : ℝ) :
    sampleRatios 0 t x = ((empMeasure x).real (Iic t), (empMeasure x).real (Ioi t)) := by
  have := isProbabilityMeasure_empMeasure hn x
  rw [sampleRatios, rLPM_zero_add_rUPM_zero hn]
  simp only [one_ne_zero, ↓reduceIte, div_one, rLPM_zero, rUPM_zero]
  rw [← mLPM_empMeasure, ← mUPM_empMeasure, (claim_C156 t).1, (claim_C156 t).2]

private theorem sampleRatios_const (hn : 0 < n) {d : ℝ} (hd : 0 < d) (c t : ℝ) :
    sampleRatios d t (fun _ : Fin n => c) =
      ((Measure.dirac c).real (Iic t), (Measure.dirac c).real (Ioi t)) := by
  simp only [sampleRatios, rLPM_const hn, rUPM_const hn, measureReal_def,
    Measure.dirac_apply' _ measurableSet_Iic, Measure.dirac_apply' _ measurableSet_Ioi,
    indicator_apply, mem_Iic, mem_Ioi]
  rcases lt_trichotomy c t with h | rfl | h
  · have hp := Real.rpow_pos_of_pos (sub_pos.mpr h) d
    simp [rLowerDev, rUpperDev, h.le, h.not_gt, hp.ne']
  · simp [rLowerDev, rUpperDev, Real.zero_rpow hd.ne']
  · have hp := Real.rpow_pos_of_pos (sub_pos.mpr h) d
    simp [rLowerDev, rUpperDev, h.not_ge, h, hp.ne']

private theorem sampleRatioLaw (hn : 0 < n) (x : Fin n → ℝ) {d : ℝ} (hd : 0 ≤ d) :
    ∃ ρ : Measure ℝ, IsProbabilityMeasure ρ ∧ ∀ t,
      sampleRatios d t x = (ρ.real (Iic t), ρ.real (Ioi t)) := by
  rcases hd.eq_or_lt with rfl | hd
  · exact ⟨empMeasure x, isProbabilityMeasure_empMeasure hn x, sampleRatios_zero hn x⟩
  by_cases hx : ∃ i j, x i ≠ x j
  · have := isProbabilityMeasure_empMeasure hn x
    have hi := GapClaims.integrable_empMeasure hn x fun v => |v| ^ d
    have ha := GapQuantile.empMeasure_compl_singleton_ne_zero hx
    let F : StieltjesFunction ℝ := ⟨fun t => GapClaims.rPMRatio d t (empMeasure x),
      GapClaims.rPMRatio_monotone hd hi,
      fun _ => (GapClaims.continuous_rPMRatio hd hi ha).continuousWithinAt⟩
    have h0 : Tendsto F atBot (𝓝 0) := GapClaims.tendsto_rPMRatio_atBot hd hi
    have h1 : Tendsto F atTop (𝓝 1) := GapClaims.tendsto_rPMRatio_atTop hd hi
    have hp : IsProbabilityMeasure F.measure :=
      ⟨by rw [F.measure_univ h0 h1, sub_zero, ENNReal.ofReal_one]⟩
    refine ⟨F.measure, hp, fun t => ?_⟩
    have hL : rlpm d t (empMeasure x) = rLPM d t x := integral_empMeasure x _
    have hU : rupm d t (empMeasure x) = rUPM d t x := integral_empMeasure x _
    have ht := GapClaims.rlpm_add_rupm_pos hd hi (ha t)
    rw [hL, hU] at ht
    have hc : F.measure.real (Iic t) = rLPM d t x / (rLPM d t x + rUPM d t x) := by
      rw [← cdf_eq_real, cdf_measure_stieltjesFunction F h0 h1]
      change GapClaims.rPMRatio d t (empMeasure x) = _
      rw [GapClaims.rPMRatio, hL, hU]
    rw [sampleRatios, ite_eq_right ht.ne', ← claim_C155, cdf_eq_real, hc]
    congr 1
    field_simp
    ring
  · push Not at hx
    let c := x ⟨0, hn⟩
    have heq : x = fun _ => c := funext fun i => hx i _
    exact ⟨Measure.dirac c, inferInstance, fun t => heq ▸ sampleRatios_const hn hd c t⟩

private theorem matrixLaw (hn : 0 < n) (x : Matrix (Fin n) (Fin k) ℝ) :
    ∃ ρ : Measure (Fin k → ℝ), IsProbabilityMeasure ρ ∧ ∀ t,
      coLPMnD_R x t 0 true = ρ.real {v | ∀ j, v j ≤ t j} ∧
      coUPMnD_R x t 0 true = ρ.real {v | ∀ j, t j < v j} := by
  classical
  let : Nonempty (Fin n) := ⟨⟨0, hn⟩⟩
  let P := uniformOn (univ : Set (Fin n))
  have hP : IsProbabilityMeasure P := inferInstanceAs (IsProbabilityMeasure (uniformOn univ))
  let rows : Fin n → (Fin k → ℝ) := fun i j => x i j
  have hx : Measurable rows := measurable_of_finite rows
  refine ⟨P.map rows, (Measure.isProbabilityMeasure_map_iff hx.aemeasurable).mpr hP, fun t => ?_⟩
  have mass (s : Set (Fin k → ℝ)) (hs : MeasurableSet s) :
      (P.map rows).real s = ((Finset.univ.filter fun i => rows i ∈ s).card : ℝ) / n := by
    rw [map_measureReal_apply hx hs, measureReal_def, uniformOn_univ]
    have heq : rows ⁻¹' s = (↑(Finset.univ.filter fun i => rows i ∈ s) : Set (Fin n)) := by
      ext; simp [mem_preimage]
    rw [heq, Measure.count_apply_finset, ENNReal.toReal_div]
    simp
  have lower : MeasurableSet {v : Fin k → ℝ | ∀ j, v j ≤ t j} := by
    rw [ofPred_forall]
    exact MeasurableSet.iInter fun j => (measurable_pi_apply j) measurableSet_Iic
  have upper : MeasurableSet {v : Fin k → ℝ | ∀ j, t j < v j} := by
    rw [ofPred_forall]
    exact MeasurableSet.iInter fun j => (measurable_pi_apply j) measurableSet_Ioi
  rw [coLPMnD_R_zero, coUPMnD_R_zero, mass _ lower, mass _ upper]
  simp only [mem_ofPred_eq, rows, and_self]

end A3_MeasureII

open A3_MeasureII

theorem claim_C538 (hn : 0 < n) (x : Fin n → ℝ) {d : ℝ} (hd : 0 ≤ d) :
    Monotone (fun t => (sampleRatios d t x).1) ∧ ∀ t, (sampleRatios d t x).1 ∈ Icc 0 1 := by
  obtain ⟨ρ, hp, h⟩ := sampleRatioLaw hn x hd
  simp_rw [h]
  exact ⟨fun _ _ hst => measureReal_mono (Iic_subset_Iic.mpr hst),
    fun _ => ⟨measureReal_nonneg, by simp⟩⟩

theorem claim_C548 (hn : 0 < n) (X : Matrix (Fin n) (Fin k) ℝ) (x : Fin n → ℝ) :
    (∃ ρ : Measure (Fin k → ℝ), IsProbabilityMeasure ρ ∧ ∀ t,
      coLPMnD_R X t 0 true = ρ.real {v | ∀ j, v j ≤ t j} ∧
      coUPMnD_R X t 0 true = ρ.real {v | ∀ j, t j < v j}) ∧
    ∀ d : ℝ, 0 ≤ d → ∃ σ : Measure ℝ, IsProbabilityMeasure σ ∧ ∀ t,
      sampleRatios d t x = (σ.real (Iic t), σ.real (Ioi t)) :=
  ⟨matrixLaw hn X, fun _ hd => sampleRatioLaw hn x hd⟩

end

end NNS.ClaimsA
