-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.GapOperators
public import NNSProofs.GapProjection
public import NNSProofs.GapMoments
public import NNSProofs.Integration
public import NNSProofs.GapDensity
public import NNSProofs.GapCoMoments

namespace NNS.GapClaims

open MeasureTheory ProbabilityTheory Set Filter Topology
open scoped ENNReal

@[expose] public noncomputable section

section DensityClaims

open GapDensity

theorem claim_C045 {f : ℝ → ℝ} [IsProbabilityMeasure (densityLaw f)] (hf : Measurable f)
    (hf0 : ∀ x, 0 ≤ f x) :
    (∀ᵐ t, 0 < mUPM 0 t (densityLaw f) →
      HasDerivAt (fun s => -Real.log (mUPM 0 s (densityLaw f)))
        (f t / mUPM 0 t (densityLaw f)) t) ∧
    ∀ t, ContinuousAt f t → 0 < mUPM 0 t (densityLaw f) →
      HasDerivAt (fun s => -Real.log (mUPM 0 s (densityLaw f)))
        (f t / mUPM 0 t (densityLaw f)) t :=
  ⟨ae_hasDerivAt_cumHazard hf hf0, fun _ hc hS => hasDerivAt_cumHazard hf hf0 hc hS⟩

theorem claim_C047 {f : ℝ → ℝ} [IsProbabilityMeasure (densityLaw f)] (hf : Measurable f)
    (hf0 : ∀ x, 0 ≤ f x) :
    (∀ t, 0 < mUPM 0 t (densityLaw f) →
      -Real.log (mUPM 0 t (densityLaw f)) = ∫ s in Set.Iic t, f s / mUPM 0 s (densityLaw f)) ∧
    (densityLaw f (Set.Iio 0) = 0 → ∀ t, 0 ≤ t → 0 < mUPM 0 t (densityLaw f) →
      -Real.log (mUPM 0 t (densityLaw f)) = ∫ s in (0)..t, f s / mUPM 0 s (densityLaw f)) :=
  ⟨fun _ hS => cumHazard_eq_integral_Iic_of_density hf hf0 hS,
    fun h0 _ ht hS => cumHazard_eq_intervalIntegral_of_density hf hf0 h0 ht hS⟩

theorem claim_C063 {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} {X : Ω → ℝ} {f g : ℝ → ℝ}
    (hX : AEMeasurable X P) (hf : Measurable f) (hf0 : ∀ x, 0 ≤ f x)
    (hlaw : P.map X = volume.withDensity fun x => ENNReal.ofReal (f x))
    (hg : Integrable g (P.map X)) : ∫ ω, g (X ω) ∂P = ∫ x, g x * f x :=
  expectation_eq_integral_mul_density hX hf hf0 hlaw hg.aestronglyMeasurable

theorem claim_C075 {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P]
    {X : Ω → ℝ} {g : ℝ → ℝ} (hX : AEMeasurable X P) (hg : Integrable g (P.map X)) :
    ∫ ω, g (X ω) ∂P = ∫ x, g x ∂(cdf (P.map X)).measure :=
  expectation_eq_integral_cdf_measure hX hg.aestronglyMeasurable

theorem claim_C076 {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} {X : Ω → ℝ} {a b : ℝ}
    {f : ℝ → ℝ} (hab : a < b) (hX : pdf.IsUniform X (Set.Icc a b) P)
    (hf : IntegrableOn f (Set.Icc a b)) :
    ∫ ω, f (X ω) ∂P = (b - a)⁻¹ * ∫ x in a..b, f x ∧
      ∫ x in a..b, f x = (b - a) * ∫ ω, f (X ω) ∂P :=
  ⟨IsUniform.expectation_eq hab.le hX hf.aestronglyMeasurable,
    IsUniform.intervalIntegral_eq hab hX hf.aestronglyMeasurable⟩

theorem claim_C080 {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} {a b : ℝ} {f : ℝ → ℝ}
    (hab : a < b) (hf : IntervalIntegrable f volume a b) :
    (∀ X : Ω → ℝ, pdf.IsUniform X (Set.Icc a b) P →
      ∫ x in a..b, f x = (b - a) *
        (mUPM 1 0 (P.map fun ω => f (X ω)) - mLPM 1 0 (P.map fun ω => f (X ω))) ∧
      ∫ x in a..b, |f x| = (b - a) *
        (mUPM 1 0 (P.map fun ω => f (X ω)) + mLPM 1 0 (P.map fun ω => f (X ω)))) ∧
    (∀ X : ℕ → Ω → ℝ,
      (∀ i, pdf.IsUniform (X i) (Set.Icc a b) P) → (Pairwise fun i j => X i ⟂ᵢ[P] X j) →
      ∀ᵐ ω ∂P,
        Tendsto (fun n : ℕ => (b - a) * (UPM 1 0 (fun i : Fin n => f (X i ω)) -
          LPM 1 0 (fun i : Fin n => f (X i ω)))) atTop (𝓝 (∫ x in a..b, f x)) ∧
        Tendsto (fun n : ℕ => (b - a) * (UPM 1 0 (fun i : Fin n => f (X i ω)) +
          LPM 1 0 (fun i : Fin n => f (X i ω)))) atTop (𝓝 (∫ x in a..b, |f x|))) :=
  ⟨fun _ hX => ⟨IsUniform.intervalIntegral_eq_partialMoments hab hX hf,
      IsUniform.intervalIntegral_abs_eq hab hX hf⟩,
    fun _ hX hind => ae_tendsto_uniform_estimate_of_intervalIntegrable hab hX hind hf⟩

theorem claim_C206 {n : ℕ} (K : ℝ → ℝ) (h : ℝ) (x : Fin n → ℝ) :
    (∀ t, kde K h x t = ((n : ℝ) * h)⁻¹ * ∑ i, K ((t - x i) / h)) ∧
    (0 < n → 0 < h → (∀ u, 0 ≤ K u) → ∫ u, K u = 1 →
      (∀ t, 0 ≤ kde K h x t) ∧ ∫ t, kde K h x t = 1) :=
  ⟨fun _ => rfl, fun hn hh hK h1 => ⟨kde_nonneg hh.le hK x, integral_kde hn hh h1 x⟩⟩

theorem claim_C272 :
    (∀ r : ℝ, 0 < r → cdf (expMeasure r) (∫ x, x ∂expMeasure r) = 1 - Real.exp (-1)) ∧
    0.632 < 1 - Real.exp (-1) ∧ 1 - Real.exp (-1) < 0.633 :=
  ⟨fun _ hr => cdf_expMeasure_mean hr, one_sub_exp_neg_one_bounds⟩

theorem claim_C418 :
    negExpProfile 1 = 1 ∧ 7.75 < negExpProfile 2 ∧ negExpProfile 2 < 7.76 ∧
      113.3 < negExpProfile 3 ∧ negExpProfile 3 < 113.5 :=
  negExpProfile_bounds

end DensityClaims

private lemma expectile_error {n : ℕ} (x : Fin n → ℝ) (hx : ∃ i j, x i ≠ x j) (t : ℝ) :
    |GapQuantile.threshold 1 (empMeasure x) (1 / 20) - t| ≤
      |19 * LPM 1 t x - UPM 1 t x| := by
  let q := GapQuantile.threshold 1 (empMeasure x) (1 / 20)
  obtain ⟨i, _, _⟩ := id hx
  have hn := Fin.pos i
  have := isProbabilityMeasure_empMeasure hn x
  have h := GapQuantile.threshold_balance (d := 1) one_ne_zero (memLp_empMeasure hn x _)
    (GapQuantile.empMeasure_compl_singleton_ne_zero hx) (a := 1 / 20) (by norm_num) (by norm_num)
  rw [mLPM_empMeasure, mUPM_empMeasure] at h
  change (1 - 1 / 20) * LPM 1 q x = 1 / 20 * UPM 1 q x at h
  change |q - t| ≤ _
  have hq := UPM_one_sub_LPM_one hn q x
  have ht := UPM_one_sub_LPM_one hn t x
  rcases le_total q t with hqt | htq
  · rw [abs_of_nonpos (sub_nonpos.mpr hqt)]
    exact le_trans (by linarith [LPM_monotone_target 1 x hqt]) (le_abs_self _)
  · rw [abs_of_nonneg (sub_nonneg.mpr htq)]
    exact le_trans (by linarith [LPM_monotone_target 1 x htq]) (neg_le_abs _)

private lemma nonconstant_of_ecdf {n : ℕ} {x : Fin n → ℝ} {t : ℝ}
    (hl : 0 < LPM 0 t x) (hu : LPM 0 t x < 1) : ∃ i j, x i ≠ x j := by
  classical
  by_contra! hall
  cases n with
  | zero => simp [LPM] at hl
  | succ n =>
    by_cases h : x 0 ≤ t
    · have hle (i) : x i ≤ t := (hall i 0).le.trans h
      simp [LPM, lowerDev, hle, div_self (by positivity : (n : ℝ) + 1 ≠ 0)] at hu
    · have hgt (i) : ¬ x i ≤ t := by rwa [hall i 0]
      simp [LPM, lowerDev, hgt] at hl

theorem claim_C274 (r : ℝ) (hr : 0 < r) :
    let t := GapQuantile.threshold 1 (expMeasure r) (1 / 20)
    0.25455 < cdf (expMeasure r) t ∧ cdf (expMeasure r) t < 0.25465 ∧
    ∀ {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P]
      (X : ℕ → Ω → ℝ), Pairwise (fun i j => X i ⟂ᵢ[P] X j) →
      (∀ i, IdentDistrib (X i) (X 0) P P) → P.map (X 0) = expMeasure r →
      ∀ᵐ ω ∂P,
        let x := fun n (i : Fin n) => X i ω
        let q := fun n => GapQuantile.threshold 1 (empMeasure (x n)) (1 / 20)
        Tendsto (fun n => LPM 0 t (x n)) atTop (𝓝 (cdf (expMeasure r) t)) ∧
        Tendsto q atTop (𝓝 t) ∧
        Tendsto (fun n => LPM 0 (q n) (x n)) atTop (𝓝 (cdf (expMeasure r) t)) ∧
        Tendsto (fun n => GapSampling.ecdfLeft (q n) (x n)) atTop (𝓝 (cdf (expMeasure r) t)) := by
  let t := GapQuantile.threshold 1 (expMeasure r) (1 / 20)
  have := isProbabilityMeasure_expMeasure hr
  obtain ⟨ht, _, hl, hu⟩ := GapDensity.threshold_one_expMeasure hr
  refine ⟨hl, hu, fun X hi he hlaw => ?_⟩
  have hm := GapDensity.memLp_one_expMeasure hr
  have hX := (he 0).aemeasurable_fst
  have hmX := (memLp_map_measure_iff aestronglyMeasurable_id hX).mp (hlaw ▸ hm)
  have hL := GapSampling.ae_tendsto_LPM X 1 t (GapSampling.integrable_lowerDev_comp hX hmX t) hi he
  have hU := GapSampling.ae_tendsto_UPM X 1 t (GapSampling.integrable_upperDev_comp hX hmX t) hi he
  have hF := GapSampling.ae_tendsto_ecdf X hi he t
  have hGC := GapSampling.glivenko_cantelli_eventually X hi he
  rw [hlaw] at hL hU hF hGC
  have ht0 : 0 < t := pos_of_mul_pos_right (lt_trans (by norm_num) ht) hr.le
  have balance := GapQuantile.threshold_balance one_ne_zero hm
    (GapDensity.expMeasure_compl_singleton_ne_zero hr) (a := 1 / 20) (by norm_num) (by norm_num)
  have hb : 19 * mLPM 1 t (expMeasure r) - mUPM 1 t (expMeasure r) = 0 := by
    dsimp only [t]; linarith
  filter_upwards [hL, hU, hF, hGC] with ω hLω hUω hFω hGCω
  let x := fun n (i : Fin n) => X i ω
  let q := fun n => GapQuantile.threshold 1 (empMeasure (x n)) (1 / 20)
  have hscore := ((hLω.const_mul 19).sub hUω).abs
  rw [hb, abs_zero] at hscore
  have hq : Tendsto q atTop (𝓝 t) := by
    rw [tendsto_iff_dist_tendsto_zero]
    apply squeeze_zero' (Eventually.of_forall fun _ => dist_nonneg) _ hscore
    filter_upwards [hFω.eventually_const_lt (lt_trans (by norm_num : (0 : ℝ) < 0.25455) hl),
      hFω.eventually_lt_const (lt_trans hu (by norm_num : (0.25465 : ℝ) < 1))] with n h1 h2
    exact expectile_error (x n) (nonconstant_of_ecdf h1 h2) t
  have inclusive (u : ℕ → ℝ) (hu' : Tendsto u atTop (𝓝 t)) :
      Tendsto (fun n => LPM 0 (u n) (x n)) atTop (𝓝 (cdf (expMeasure r) t)) := by
    have hc : Tendsto (fun n => 1 - Real.exp (-(r * u n))) atTop
        (𝓝 (1 - Real.exp (-(r * t)))) :=
      tendsto_const_nhds.sub (Real.continuous_exp.continuousAt.tendsto.comp (hu'.const_mul r).neg)
    have hcdf : Tendsto (fun n => cdf (expMeasure r) (u n)) atTop
        (𝓝 (cdf (expMeasure r) t)) := by
      rw [cdf_expMeasure_eq hr t, ite_eq_left ht0.le]
      refine hc.congr' ?_
      filter_upwards [hu'.eventually_const_lt ht0] with n hn
      rw [cdf_expMeasure_eq hr, ite_eq_left hn.le]
    apply hcdf.congr_dist
    rw [Metric.tendsto_nhds]
    intro ε hε
    filter_upwards [hGCω (ε / 2) (half_pos hε)] with n hn
    simpa [Real.dist_eq, abs_sub_comm, abs_abs] using (hn (u n)).trans_lt (half_lt_self hε)
  have hinc := inclusive q hq
  refine ⟨hFω, hq, hinc, ?_⟩
  let v := fun n : ℕ => q n - 1 / ((n : ℝ) + 1)
  have hv : Tendsto v atTop (𝓝 t) := by
    simpa [v] using hq.sub (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ))
  apply tendsto_of_tendsto_of_tendsto_of_le_of_le (inclusive v hv) hinc
  · intro n
    dsimp only
    rw [GapSampling.ecdfLeft_eq_card, LPM_zero_eq_ecdf]
    apply div_le_div_of_nonneg_right _ n.cast_nonneg
    exact_mod_cast Finset.card_le_card (Finset.monotone_filter_right Finset.univ
      fun i _ (hi' : x n i ≤ v n) => hi'.trans_lt (sub_lt_self _ (by positivity)))
  · intro n
    dsimp only
    rw [GapSampling.ecdfLeft_eq_card, LPM_zero_eq_ecdf]
    apply div_le_div_of_nonneg_right _ n.cast_nonneg
    exact_mod_cast Finset.card_le_card (Finset.monotone_filter_right Finset.univ
      fun i _ (hi' : x n i < q n) => hi'.le)

theorem claim_C555 (r : ℝ) (hr : 0 < r) :
    let t := GapQuantile.threshold 1 (expMeasure r) (1 / 20)
    0.25455 < cdf (expMeasure r) t ∧ cdf (expMeasure r) t < 0.25465 ∧
    ∀ {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P]
      (X : ℕ → Ω → ℝ), Pairwise (fun i j => X i ⟂ᵢ[P] X j) →
      (∀ i, IdentDistrib (X i) (X 0) P P) → P.map (X 0) = expMeasure r →
      ∀ᵐ ω ∂P,
        let x := fun n (i : Fin n) => X i ω
        let q := fun n => GapQuantile.threshold 1 (empMeasure (x n)) (1 / 20)
        Tendsto (fun n => LPM 0 t (x n)) atTop (𝓝 (cdf (expMeasure r) t)) ∧
        Tendsto q atTop (𝓝 t) ∧
        Tendsto (fun n => LPM 0 (q n) (x n)) atTop (𝓝 (cdf (expMeasure r) t)) ∧
        Tendsto (fun n => GapSampling.ecdfLeft (q n) (x n)) atTop (𝓝 (cdf (expMeasure r) t)) :=
  claim_C274 r hr

section SamplingClaims

open GapSampling

variable {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}

theorem claim_C133 [IsFiniteMeasure P] {X T : Ω → ℝ} (hX : Measurable X) :
    (∀ t : ℝ, MeasurableSet {ω | X ω ≤ t} ∧ MeasurableSet {ω | t < X ω}) ∧
    (∀ (r : ℕ) (t : ℝ), MemLp X r P →
      Integrable (fun ω => lowerDev r t (X ω)) P ∧ Integrable (fun ω => upperDev r t (X ω)) P) ∧
    (Measurable T →
      MeasurableSet {ω | X ω ≤ T ω} ∧ MeasurableSet {ω | T ω < X ω} ∧
      (∀ r : ℕ, Measurable (fun ω => lowerDev r (T ω) (X ω)) ∧
        Measurable (fun ω => upperDev r (T ω) (X ω))) ∧
      ∀ r : ℕ, MemLp X r P → MemLp T r P →
        Integrable (fun ω => lowerDev r (T ω) (X ω)) P ∧
          Integrable (fun ω => upperDev r (T ω) (X ω)) P) :=
  ⟨fun t => ⟨measurableSet_le_event hX t, measurableSet_gt_event hX t⟩,
    fun _ t h => ⟨integrable_lowerDev_comp hX.aemeasurable h t,
      integrable_upperDev_comp hX.aemeasurable h t⟩,
    fun hT => ⟨measurableSet_le_benchmark hX hT, measurableSet_gt_benchmark hX hT,
      fun r => ⟨measurable_lowerDev_benchmark hX hT r, measurable_upperDev_benchmark hX hT r⟩,
      fun _ hXr hTr => ⟨integrable_lowerDev_benchmark hX hT hXr hTr,
        integrable_upperDev_benchmark hX hT hXr hTr⟩⟩⟩

theorem claim_C134 [IsProbabilityMeasure P] {X : Ω → ℝ} (hX : AEMeasurable X P) (t : ℝ) :
    cdf (P.map X) t = P.real {ω | X ω ≤ t} :=
  cdf_map_eq_prob hX t

theorem claim_C135 {X : Ω → ℝ} (hX : AEMeasurable X P) :
    (∀ g : ℝ → ℝ, AEStronglyMeasurable g (P.map X) → ∫ x, g x ∂P.map X = ∫ ω, g (X ω) ∂P) ∧
    ∀ g : ℝ → ℝ≥0∞, AEMeasurable g (P.map X) → ∫⁻ x, g x ∂P.map X = ∫⁻ ω, g (X ω) ∂P :=
  ⟨fun _ hg => integral_law_eq_integral hX hg, fun _ hg => lintegral_law_eq_lintegral hX hg⟩

theorem claim_C154 [IsProbabilityMeasure P] {A : Set Ω} (hA : MeasurableSet A) :
    P.real A + P.real Aᶜ = 1 ∧ P A + P Aᶜ = 1 :=
  ⟨prob_event_add_prob_compl hA, measure_add_measure_compl_eq_one hA⟩

theorem claim_C507 {X : Ω → ℝ} (hX : Measurable X) (t : ℝ) :
    {ω | X ω ≤ t} = X ⁻¹' Set.Iic t ∧ {ω | t < X ω} = X ⁻¹' Set.Ioi t ∧
      MeasurableSet {ω | X ω ≤ t} ∧ MeasurableSet {ω | t < X ω} :=
  ⟨le_event_eq_preimage X t, gt_event_eq_preimage X t, measurableSet_le_event hX t,
    measurableSet_gt_event hX t⟩

theorem claim_C251 {Θ : Type*} [MeasurableSpace Θ] {π : Measure Θ} {κ : Θ → Measure ℝ}
    (hκ : AEMeasurable κ π) (r : ℕ) (t : ℝ) :
    ([IsProbabilityMeasure π] → (∀ θ, IsProbabilityMeasure (κ θ)) →
      IsProbabilityMeasure (π.bind κ)) ∧
    eLPM r t (π.bind κ) = ∫⁻ θ, eLPM r t (κ θ) ∂π ∧
    eUPM r t (π.bind κ) = ∫⁻ θ, eUPM r t (κ θ) ∂π ∧
    (Integrable (lowerDev r t) (π.bind κ) → mLPM r t (π.bind κ) = ∫ θ, mLPM r t (κ θ) ∂π) ∧
    (Integrable (upperDev r t) (π.bind κ) → mUPM r t (π.bind κ) = ∫ θ, mUPM r t (κ θ) ∂π) ∧
    ∀ {n : ℕ} (x : Fin n → ℝ),
      mLPM r t (empMeasure x) = LPM r t x ∧ mUPM r t (empMeasure x) = UPM r t x :=
  ⟨fun h => isProbabilityMeasure_bind hκ h, eLPM_bind hκ r t, eUPM_bind hκ r t,
    mLPM_bind hκ, mUPM_bind hκ,
    fun x => ⟨mLPM_frequentist r t x, mUPM_frequentist r t x⟩⟩

theorem claim_C416 [IsProbabilityMeasure P] (X Y : ℕ → Ω → ℝ) (r s tx ty : ℝ)
    (hindep : Pairwise fun i j => (fun ω => (X i ω, Y i ω)) ⟂ᵢ[P] fun ω => (X j ω, Y j ω))
    (hident : ∀ i, IdentDistrib (fun ω => (X i ω, Y i ω)) (fun ω => (X 0 ω, Y 0 ω)) P P)
    (h1 : Integrable (fun ω => rLowerDev r tx (X 0 ω) * rLowerDev s ty (Y 0 ω)) P)
    (h2 : Integrable (fun ω => rUpperDev r tx (X 0 ω) * rUpperDev s ty (Y 0 ω)) P)
    (h3 : Integrable (fun ω => rUpperDev r tx (X 0 ω) * rLowerDev s ty (Y 0 ω)) P)
    (h4 : Integrable (fun ω => rLowerDev r tx (X 0 ω) * rUpperDev s ty (Y 0 ω)) P) :
    ∀ᵐ ω ∂P,
      Tendsto (fun n : ℕ => GapCoMoments.rCoLPM r s tx ty (fun i : Fin n => X i ω)
          fun i : Fin n => Y i ω) atTop
        (𝓝 (∫ p, rLowerDev r tx p.1 * rLowerDev s ty p.2 ∂P.map fun ω => (X 0 ω, Y 0 ω))) ∧
      Tendsto (fun n : ℕ => GapCoMoments.rCoUPM r s tx ty (fun i : Fin n => X i ω)
          fun i : Fin n => Y i ω) atTop
        (𝓝 (∫ p, rUpperDev r tx p.1 * rUpperDev s ty p.2 ∂P.map fun ω => (X 0 ω, Y 0 ω))) ∧
      Tendsto (fun n : ℕ => GapCoMoments.rDLPM s r tx ty (fun i : Fin n => X i ω)
          fun i : Fin n => Y i ω) atTop
        (𝓝 (∫ p, rUpperDev r tx p.1 * rLowerDev s ty p.2 ∂P.map fun ω => (X 0 ω, Y 0 ω))) ∧
      Tendsto (fun n : ℕ => GapCoMoments.rDUPM r s tx ty (fun i : Fin n => X i ω)
          fun i : Fin n => Y i ω) atTop
        (𝓝 (∫ p, rLowerDev r tx p.1 * rUpperDev s ty p.2 ∂P.map fun ω => (X 0 ω, Y 0 ω))) := by
  have mL := GapCoMoments.measurable_rLowerDev'
  have mU := GapCoMoments.measurable_rUpperDev'
  filter_upwards [ae_tendsto_average (fun i ω => (X i ω, Y i ω))
      (((mL r tx).comp measurable_fst).mul ((mL s ty).comp measurable_snd)) h1 hindep hident,
    ae_tendsto_average (fun i ω => (X i ω, Y i ω))
      (((mU r tx).comp measurable_fst).mul ((mU s ty).comp measurable_snd)) h2 hindep hident,
    ae_tendsto_average (fun i ω => (X i ω, Y i ω))
      (((mU r tx).comp measurable_fst).mul ((mL s ty).comp measurable_snd)) h3 hindep hident,
    ae_tendsto_average (fun i ω => (X i ω, Y i ω))
      (((mL r tx).comp measurable_fst).mul ((mU s ty).comp measurable_snd)) h4 hindep hident]
    with ω a1 a2 a3 a4
  exact ⟨a1, a2, a3, a4⟩

theorem claim_C200 [IsProbabilityMeasure P] (X : ℕ → Ω → ℝ)
    (hindep : Pairwise fun i j => X i ⟂ᵢ[P] X j) (hident : ∀ i, IdentDistrib (X i) (X 0) P P) :
    ∀ᵐ ω ∂P, Tendsto (fun n : ℕ => ⨆ t : ℝ, |LPM 0 t (fun i : Fin n => X i ω) -
      cdf (P.map (X 0)) t|) atTop (𝓝 0) :=
  glivenko_cantelli X hindep hident

theorem claim_C218 [IsProbabilityMeasure P] (X : ℕ → Ω → ℝ)
    (hindep : Pairwise fun i j => X i ⟂ᵢ[P] X j) (hident : ∀ i, IdentDistrib (X i) (X 0) P P) :
    (∀ᵐ ω ∂P, Tendsto (fun n : ℕ => ⨆ t : ℝ, |LPM 0 t (fun i : Fin n => X i ω) -
      cdf (P.map (X 0)) t|) atTop (𝓝 0)) ∧
    ∀ᵐ ω ∂P, Tendsto (fun n : ℕ => ⨆ t : ℝ, |cdf (empMeasure fun i : Fin (n + 1) => X i ω) t -
      cdf (P.map (X 0)) t|) atTop (𝓝 0) :=
  ⟨glivenko_cantelli X hindep hident, glivenko_cantelli_empMeasure X hindep hident⟩

end SamplingClaims

section Moments

open GapMoments

theorem claim_C177 :
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ], MemLp (fun x : ℝ => x) 3 μ →
      0 < variance (fun x : ℝ => x) μ →
      mSkewness μ = (∫ x, (x - ∫ y, y ∂μ) ^ 3 ∂μ) / variance (fun x : ℝ => x) μ ^ ((3 : ℝ) / 2) ∧
      mSkewness μ = (mUPM 3 (∫ y, y ∂μ) μ - mLPM 3 (∫ y, y ∂μ) μ) /
        (mUPM 2 (∫ y, y ∂μ) μ + mLPM 2 (∫ y, y ∂μ) μ) ^ ((3 : ℝ) / 2)) ∧
    (∀ (n : ℕ) (x : Fin n → ℝ), 0 < popVariance x →
      popSkewness x = mSkewness (empMeasure x) ∧
      popSkewness x = centralMoment 3 (mean x) x / popVariance x ^ ((3 : ℝ) / 2)) :=
  ⟨fun _ _ h _ => ⟨rfl, mSkewness_eq_PM h⟩,
    fun _ x _ => ⟨(mSkewness_empMeasure x).symm, rfl⟩⟩

theorem claim_C185 :
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ], MemLp (fun x : ℝ => x) 4 μ →
      0 < variance (fun x : ℝ => x) μ →
      mKurtosis μ = (∫ x, (x - ∫ y, y ∂μ) ^ 4 ∂μ) / variance (fun x : ℝ => x) μ ^ 2 ∧
      mKurtosis μ = (mUPM 4 (∫ y, y ∂μ) μ + mLPM 4 (∫ y, y ∂μ) μ) /
        (mUPM 2 (∫ y, y ∂μ) μ + mLPM 2 (∫ y, y ∂μ) μ) ^ 2) ∧
    (∀ (n : ℕ) (x : Fin n → ℝ), 0 < popVariance x →
      popKurtosis x + 3 = mKurtosis (empMeasure x) ∧
      popKurtosis x + 3 = (UPM 4 (mean x) x + LPM 4 (mean x) x) /
        (LPM 2 (mean x) x + UPM 2 (mean x) x) ^ 2) :=
  ⟨fun _ _ h _ => ⟨rfl, mKurtosis_eq_PM h⟩,
    fun _ x _ => ⟨(mKurtosis_empMeasure x).symm, by rw [popKurtosis_eq_PM]; ring⟩⟩

theorem claim_C191 :
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ],
      (Integrable (fun x : ℝ => x) μ →
        (∀ t : ℝ, mUPM 1 t μ = mLPM 1 t μ ↔ t = ∫ y, y ∂μ) ∧
        ∫ y, y ∂μ = mUPM 1 0 μ - mLPM 1 0 μ) ∧
      (MemLp (fun x : ℝ => x) 2 μ →
        variance (fun x : ℝ => x) μ = mUPM 2 (∫ y, y ∂μ) μ + mLPM 2 (∫ y, y ∂μ) μ) ∧
      (MemLp (fun x : ℝ => x) 3 μ → 0 < variance (fun x : ℝ => x) μ →
        mSkewness μ = (mUPM 3 (∫ y, y ∂μ) μ - mLPM 3 (∫ y, y ∂μ) μ) /
          (mUPM 2 (∫ y, y ∂μ) μ + mLPM 2 (∫ y, y ∂μ) μ) ^ ((3 : ℝ) / 2)) ∧
      (MemLp (fun x : ℝ => x) 4 μ → 0 < variance (fun x : ℝ => x) μ →
        mKurtosis μ = (mUPM 4 (∫ y, y ∂μ) μ + mLPM 4 (∫ y, y ∂μ) μ) /
          (mUPM 2 (∫ y, y ∂μ) μ + mLPM 2 (∫ y, y ∂μ) μ) ^ 2)) ∧
    (∀ (n : ℕ) (x : Fin n → ℝ),
      (0 < n →
        (∀ t : ℝ, UPM 1 t x = LPM 1 t x ↔ t = mean x) ∧
        mean x = UPM 1 0 x - LPM 1 0 x ∧
        popVariance x = LPM 2 (mean x) x + UPM 2 (mean x) x) ∧
      (0 < popVariance x →
        popSkewness x = (UPM 3 (mean x) x - LPM 3 (mean x) x) /
          (LPM 2 (mean x) x + UPM 2 (mean x) x) ^ ((3 : ℝ) / 2) ∧
        popKurtosis x + 3 = (UPM 4 (mean x) x + LPM 4 (mean x) x) /
          (LPM 2 (mean x) x + UPM 2 (mean x) x) ^ 2)) := by
  refine ⟨fun _ _ => ⟨fun h1 => ⟨mUPM_one_eq_mLPM_one_iff h1, mean_eq_mUPM_one_sub_mLPM_one h1⟩,
    variance_eq_mUPM_two_add_mLPM_two, fun h3 _ => mSkewness_eq_PM h3,
    fun h4 _ => mKurtosis_eq_PM h4⟩, fun n x => ⟨fun hn => ?_, fun hv => ?_⟩⟩
  · obtain ⟨h1, h2, h3, -⟩ := classical_moments_from_pairs_sample hn x
    exact ⟨h1, h2, h3⟩
  · rcases Nat.eq_zero_or_pos n with rfl | hn
    · simp [popVariance, centralMoment] at hv
    · exact (classical_moments_from_pairs_sample hn x).2.2.2

theorem claim_C114 (μ : Measure ℝ) [IsProbabilityMeasure μ]
    (h : Integrable (fun x : ℝ => x) μ) (t : ℝ) :
    (mUPM 1 t μ = mLPM 1 t μ ↔ t = ∫ x, x ∂μ) ∧
    (mLPM 1 t μ < mUPM 1 t μ ↔ t < ∫ x, x ∂μ) ∧ (mUPM 1 t μ < mLPM 1 t μ ↔ ∫ x, x ∂μ < t) :=
  ⟨mUPM_one_eq_mLPM_one_iff h t, mLPM_one_lt_mUPM_one_iff h t, mUPM_one_lt_mLPM_one_iff h t⟩

theorem claim_C117 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ) :
    variance (fun x : ℝ => x) μ = mUPM 2 (∫ y, y ∂μ) μ + mLPM 2 (∫ y, y ∂μ) μ ∧
    (0 < variance (fun x : ℝ => x) μ →
      mLPM 2 (∫ y, y ∂μ) μ = μ.real (Set.Iic (∫ y, y ∂μ)) *
        (condVar μ (Set.Iic (∫ y, y ∂μ)) + (condMean μ (Set.Iic (∫ y, y ∂μ)) - ∫ y, y ∂μ) ^ 2) ∧
      mUPM 2 (∫ y, y ∂μ) μ = μ.real (Set.Ioi (∫ y, y ∂μ)) *
        (condVar μ (Set.Ioi (∫ y, y ∂μ)) + (condMean μ (Set.Ioi (∫ y, y ∂μ)) - ∫ y, y ∂μ) ^ 2)) ∧
    ∀ A : Set ℝ, MeasurableSet A →
      (μ A ≠ 0 → μ Aᶜ ≠ 0 → variance (fun x : ℝ => x) μ =
        (μ.real A * condVar μ A + μ.real Aᶜ * condVar μ Aᶜ) +
        (μ.real A * (condMean μ A - ∫ y, y ∂μ) ^ 2 +
          μ.real Aᶜ * (condMean μ Aᶜ - ∫ y, y ∂μ) ^ 2)) ∧
      0 ≤ μ.real A * (condMean μ A - ∫ y, y ∂μ) ^ 2 +
        μ.real Aᶜ * (condMean μ Aᶜ - ∫ y, y ∂μ) ^ 2 ∧
      μ.real A * condVar μ A + μ.real Aᶜ * condVar μ Aᶜ ≤ variance (fun x : ℝ => x) μ :=
  ⟨variance_eq_mUPM_two_add_mLPM_two h,
    fun hv => ⟨(variance_split_at_mean h hv).1, (variance_split_at_mean h hv).2.1⟩,
    fun _ hA => ⟨variance_eq_within_add_between h hA, between_nonneg μ _ _,
      within_le_variance h hA⟩⟩

theorem claim_C121 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 3 μ)
    (_hv : 0 < variance (fun x : ℝ => x) μ) :
    mSkewness μ = (∫ x, (x - ∫ y, y ∂μ) ^ 3 ∂μ) / variance (fun x : ℝ => x) μ ^ ((3 : ℝ) / 2) ∧
    mSkewness μ = (mUPM 3 (∫ y, y ∂μ) μ - mLPM 3 (∫ y, y ∂μ) μ) /
      (mUPM 2 (∫ y, y ∂μ) μ + mLPM 2 (∫ y, y ∂μ) μ) ^ ((3 : ℝ) / 2) :=
  ⟨rfl, mSkewness_eq_PM h⟩

theorem claim_C179 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 3 μ)
    (_hv : 0 < variance (fun x : ℝ => x) μ) :
    mSkewness μ = (mUPM 3 (∫ y, y ∂μ) μ - mLPM 3 (∫ y, y ∂μ) μ) /
      (mUPM 2 (∫ y, y ∂μ) μ + mLPM 2 (∫ y, y ∂μ) μ) ^ ((3 : ℝ) / 2) :=
  mSkewness_eq_PM h

theorem claim_C123 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 4 μ)
    (_hv : 0 < variance (fun x : ℝ => x) μ) :
    mKurtosis μ = (∫ x, (x - ∫ y, y ∂μ) ^ 4 ∂μ) / variance (fun x : ℝ => x) μ ^ 2 ∧
    mKurtosis μ = (mUPM 4 (∫ y, y ∂μ) μ + mLPM 4 (∫ y, y ∂μ) μ) /
      (mUPM 2 (∫ y, y ∂μ) μ + mLPM 2 (∫ y, y ∂μ) μ) ^ 2 :=
  ⟨rfl, mKurtosis_eq_PM h⟩

theorem claim_C180 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 3 μ) :
    (0 < mSkewness μ ↔ mLPM 3 (∫ y, y ∂μ) μ < mUPM 3 (∫ y, y ∂μ) μ) ∧
    (mSkewness μ < 0 ↔ mUPM 3 (∫ y, y ∂μ) μ < mLPM 3 (∫ y, y ∂μ) μ) :=
  ⟨mSkewness_pos_iff h, mSkewness_neg_iff h⟩

theorem claim_C127 (μ : Measure ℝ) [IsProbabilityMeasure μ] :
    (Integrable (fun x : ℝ => x) μ →
      (eLPM 2 (∫ y, y ∂μ) μ = 0 → (∀ᵐ x ∂μ, ∫ y, y ∂μ ≤ x) ∧ ∀ᵐ x ∂μ, x = ∫ y, y ∂μ) ∧
      (eUPM 2 (∫ y, y ∂μ) μ = 0 → (∀ᵐ x ∂μ, x ≤ ∫ y, y ∂μ) ∧ ∀ᵐ x ∂μ, x = ∫ y, y ∂μ)) ∧
    (MemLp (fun x : ℝ => x) 2 μ →
      (mLPM 2 (∫ y, y ∂μ) μ = 0 ↔ ∀ᵐ x ∂μ, x = ∫ y, y ∂μ) ∧
      (mUPM 2 (∫ y, y ∂μ) μ = 0 ↔ ∀ᵐ x ∂μ, x = ∫ y, y ∂μ)) :=
  ⟨fun h => ⟨ae_eq_mean_of_eLPM_two_eq_zero h, ae_eq_mean_of_eUPM_two_eq_zero h⟩,
    fun h => ⟨mLPM_two_mean_eq_zero_iff h, mUPM_two_mean_eq_zero_iff h⟩⟩

theorem claim_C128 :
    (∀ c : ℝ, ∫ x, x ∂Measure.dirac c = c ∧ mUPM 2 c (Measure.dirac c) = 0 ∧
      mLPM 2 c (Measure.dirac c) = 0 ∧ variance (fun x : ℝ => x) (Measure.dirac c) = 0) ∧
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ], MemLp (fun x : ℝ => x) 2 μ →
      0 < variance (fun x : ℝ => x) μ →
      0 < mLPM 2 (∫ y, y ∂μ) μ ∧ 0 < mUPM 2 (∫ y, y ∂μ) μ) :=
  ⟨dirac_degenerate, fun _ _ h hv => mLPM_two_pos_and_mUPM_two_pos h hv⟩

theorem claim_C166 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ) :
    μ.real (Set.Ici (∫ y, y ∂μ)) * condVar μ (Set.Ici (∫ y, y ∂μ)) +
        (1 - μ.real (Set.Ici (∫ y, y ∂μ))) * condVar μ (Set.Iio (∫ y, y ∂μ)) ≤
      variance (fun x : ℝ => x) μ ∧
    (0 < variance (fun x : ℝ => x) μ → variance (fun x : ℝ => x) μ =
      μ.real (Set.Ici (∫ y, y ∂μ)) * condVar μ (Set.Ici (∫ y, y ∂μ)) +
        (1 - μ.real (Set.Ici (∫ y, y ∂μ))) * condVar μ (Set.Iio (∫ y, y ∂μ)) +
        μ.real (Set.Ici (∫ y, y ∂μ)) * (1 - μ.real (Set.Ici (∫ y, y ∂μ))) *
          (condMean μ (Set.Ici (∫ y, y ∂μ)) - condMean μ (Set.Iio (∫ y, y ∂μ))) ^ 2) :=
  ⟨within_le_variance_Ici h _, variance_eq_total_variance_mean h⟩

theorem claim_C485 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ)
    (hv : 0 < variance (fun x : ℝ => x) μ) :
    mLPM 2 (∫ y, y ∂μ) μ = μ.real (Set.Iic (∫ y, y ∂μ)) *
        (condVar μ (Set.Iic (∫ y, y ∂μ)) + (condMean μ (Set.Iic (∫ y, y ∂μ)) - ∫ y, y ∂μ) ^ 2) ∧
      mUPM 2 (∫ y, y ∂μ) μ = μ.real (Set.Ioi (∫ y, y ∂μ)) *
        (condVar μ (Set.Ioi (∫ y, y ∂μ)) + (condMean μ (Set.Ioi (∫ y, y ∂μ)) - ∫ y, y ∂μ) ^ 2) ∧
      variance (fun x : ℝ => x) μ =
        (μ.real (Set.Iic (∫ y, y ∂μ)) * condVar μ (Set.Iic (∫ y, y ∂μ)) +
          μ.real (Set.Ioi (∫ y, y ∂μ)) * condVar μ (Set.Ioi (∫ y, y ∂μ))) +
        (μ.real (Set.Iic (∫ y, y ∂μ)) * (condMean μ (Set.Iic (∫ y, y ∂μ)) - ∫ y, y ∂μ) ^ 2 +
          μ.real (Set.Ioi (∫ y, y ∂μ)) * (condMean μ (Set.Ioi (∫ y, y ∂μ)) - ∫ y, y ∂μ) ^ 2) ∧
      μ.real (Set.Iic (∫ y, y ∂μ)) * condVar μ (Set.Iic (∫ y, y ∂μ)) +
          μ.real (Set.Iic (∫ y, y ∂μ))ᶜ * condVar μ (Set.Iic (∫ y, y ∂μ))ᶜ ≤
        variance (fun x : ℝ => x) μ :=
  ⟨(variance_split_at_mean h hv).1, (variance_split_at_mean h hv).2.1,
    (variance_split_at_mean h hv).2.2, within_le_variance h measurableSet_Iic⟩

theorem claim_C266 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : Integrable (fun x : ℝ => x) μ)
    (hμ : ¬ ∃ c : ℝ, ∀ᵐ x ∂μ, x = c) : mPMRatio 1 (∫ y, y ∂μ) μ = 1 / 2 :=
  mPMRatio_one_mean_of_not_ae_const h hμ

theorem claim_C187 :
    0 < popVariance tailSample ∧ popVariance (-tailSample) = popVariance tailSample ∧
      popKurtosis (-tailSample) = popKurtosis tailSample ∧
      UPM 4 (mean tailSample) tailSample = 16 / 3 ∧
      LPM 4 (mean tailSample) tailSample = 2 / 3 ∧
      UPM 4 (mean (-tailSample)) (-tailSample) = 2 / 3 ∧
      LPM 4 (mean (-tailSample)) (-tailSample) = 16 / 3 :=
  tailSample_witness

theorem claim_C186 :
    (∀ a : ℝ, a ≠ 0 → (∀ i, |threePoint a i| ≤ |a|) ∧
      mKurtosis (empMeasure (threePoint a)) = 9) ∧
    (mKurtosis scaledT5 = 9 ∧ ∀ M : ℝ, scaledT5 {y | M < |y|} ≠ 0) ∧
    (popKurtosis (-tailSample) = popKurtosis tailSample ∧
      UPM 4 (mean tailSample) tailSample = 16 / 3 ∧
      LPM 4 (mean tailSample) tailSample = 2 / 3 ∧
      UPM 4 (mean (-tailSample)) (-tailSample) = 2 / 3 ∧
      LPM 4 (mean (-tailSample)) (-tailSample) = 16 / 3 ∧
      LPM 4 (mean tailSample) tailSample < UPM 4 (mean tailSample) tailSample ∧
      UPM 4 (mean (-tailSample)) (-tailSample) < LPM 4 (mean (-tailSample)) (-tailSample)) := by
  obtain ⟨-, -, hk, hU, hL, hU', hL'⟩ := tailSample_witness
  refine ⟨fun a ha => ⟨threePoint_bounded a, (threePoint_kurtosis ha).2.2⟩,
    ⟨mKurtosis_scaledT5, scaledT5_unbounded⟩, hk, hU, hL, hU', hL', ?_, ?_⟩
  · rw [hU, hL]; norm_num
  · rw [hU', hL']; norm_num

end Moments

section Projection

open GapProjection

lemma sampleCorr_self {n : ℕ} (x : Fin n → ℝ) (hv : 0 < popVariance x) : sampleCorr x x = 1 := by
  rw [sampleCorr, sample_covariance_self, Real.mul_self_sqrt hv.le, div_self hv.ne']

abbrev StrictRefinement (r : ℕ) (t : ℝ) : Prop :=
  (∀ μ ν : Measure ℝ, IsProbabilityMeasure μ → IsProbabilityMeasure ν →
      MemLp (fun x : ℝ => x) r μ → MemLp (fun x : ℝ => x) r ν →
      (mUPM r t μ, mLPM r t μ) = (mUPM r t ν, mLPM r t ν) →
      ∫ x, (x - t) ^ r ∂μ = ∫ x, (x - t) ^ r ∂ν) ∧
    ∃ μ ν : Measure ℝ, IsProbabilityMeasure μ ∧ IsProbabilityMeasure ν ∧
      MemLp (fun x : ℝ => x) r μ ∧ MemLp (fun x : ℝ => x) r ν ∧
      ∫ x, (x - t) ^ r ∂μ = ∫ x, (x - t) ^ r ∂ν ∧
      (mUPM r t μ, mLPM r t μ) ≠ (mUPM r t ν, mLPM r t ν)

lemma strictRefinement (r : ℕ) (t : ℝ) : StrictRefinement r t :=
  ⟨(pair_partition_strictly_refines r t).1, exists_laws_moment_eq_pair_ne r t⟩

theorem claim_C026 : ∀ (r : ℕ) (t : ℝ),
    (∀ (μ : Measure ℝ) [IsFiniteMeasure μ], MemLp (fun x : ℝ => x) r μ →
      ∫ x, (x - t) ^ r ∂μ = symProj r (mUPM r t μ, mLPM r t μ)) ∧
    StrictRefinement r t ∧
    ¬ ∃ g : ℝ → ℝ × ℝ, ∀ μ : Measure ℝ, IsProbabilityMeasure μ →
      MemLp (fun x : ℝ => x) r μ → g (∫ x, (x - t) ^ r ∂μ) = (mUPM r t μ, mLPM r t μ) :=
  fun r t => ⟨fun _ _ h => integral_sub_pow_eq_symProj h t, strictRefinement r t,
    not_exists_recover_pair r t⟩

theorem claim_C125 : ∀ (r : ℕ) (t : ℝ),
    (∀ (μ : Measure ℝ) [IsFiniteMeasure μ], MemLp (fun x : ℝ => x) r μ →
      ∫ x, (x - t) ^ r ∂μ = symProj r (mUPM r t μ, mLPM r t μ)) ∧
    StrictRefinement r t ∧
    ¬ ∃ g : ℝ → ℝ × ℝ, ∀ μ : Measure ℝ, IsProbabilityMeasure μ →
      MemLp (fun x : ℝ => x) r μ → g (∫ x, (x - t) ^ r ∂μ) = (mUPM r t μ, mLPM r t μ) :=
  claim_C026

theorem claim_C151 : ∀ (r : ℕ) (t : ℝ),
    (∀ (μ : Measure ℝ) [IsFiniteMeasure μ], MemLp (fun x : ℝ => x) r μ →
      ∫ x, (x - t) ^ r ∂μ = symProj r (mUPM r t μ, mLPM r t μ)) ∧
    StrictRefinement r t ∧
    ¬ ∃ g : ℝ → ℝ × ℝ, ∀ μ : Measure ℝ, IsProbabilityMeasure μ →
      MemLp (fun x : ℝ => x) r μ → g (∫ x, (x - t) ^ r ∂μ) = (mUPM r t μ, mLPM r t μ) :=
  claim_C026

theorem claim_C152 : ∀ (r : ℕ) (t : ℝ), StrictRefinement r t := strictRefinement

theorem claim_C028 : ∀ (r : ℕ) (t : ℝ),
    ¬ ∃ g : ℝ → ℝ × ℝ, ∀ μ : Measure ℝ, IsProbabilityMeasure μ →
      MemLp (fun x : ℝ => x) r μ → g (∫ x, (x - t) ^ r ∂μ) = (mUPM r t μ, mLPM r t μ) :=
  not_exists_recover_pair

theorem claim_C129 : ∀ (r : ℕ) (t : ℝ),
    (∀ (μ : Measure ℝ) [IsFiniteMeasure μ], MemLp (fun x : ℝ => x) r μ →
      ∫ x, (x - t) ^ r ∂μ = symProj r (mUPM r t μ, mLPM r t μ)) ∧
    ¬ ∃ g : ℝ → ℝ × ℝ, ∀ μ : Measure ℝ, IsProbabilityMeasure μ →
      MemLp (fun x : ℝ => x) r μ → g (∫ x, (x - t) ^ r ∂μ) = (mUPM r t μ, mLPM r t μ) :=
  fun r t => ⟨(claim_C026 r t).1, not_exists_recover_pair r t⟩

theorem claim_C380 :
    (∀ (r : ℕ) (t : ℝ), StrictRefinement r t ∧
      ¬ ∃ g : ℝ → ℝ × ℝ, ∀ μ : Measure ℝ, IsProbabilityMeasure μ →
        MemLp (fun x : ℝ => x) r μ → g (∫ x, (x - t) ^ r ∂μ) = (mUPM r t μ, mLPM r t μ)) ∧
    (∀ μ ν : Measure ℝ, IsProbabilityMeasure μ → IsProbabilityMeasure ν →
      MemLp (fun x : ℝ => x) 2 μ → MemLp (fun x : ℝ => x) 2 ν →
      (mUPM 2 (∫ x, x ∂μ) μ, mLPM 2 (∫ x, x ∂μ) μ)
        = (mUPM 2 (∫ x, x ∂ν) ν, mLPM 2 (∫ x, x ∂ν) ν) →
      Var[fun x : ℝ => x; μ] = Var[fun x : ℝ => x; ν]) ∧
    ¬ ∃ g : ℝ → ℝ × ℝ, ∀ μ : Measure ℝ, IsProbabilityMeasure μ → MemLp (fun x : ℝ => x) 2 μ →
      g Var[fun x : ℝ => x; μ] = (mUPM 2 (∫ x, x ∂μ) μ, mLPM 2 (∫ x, x ∂μ) μ) :=
  ⟨fun r t => ⟨strictRefinement r t, not_exists_recover_pair r t⟩,
    pair_partition_strictly_refines_variance.1, not_exists_recover_pair_of_variance⟩

theorem claim_C027 :
    ∃ μ ν : Measure ℝ, IsProbabilityMeasure μ ∧ IsProbabilityMeasure ν ∧
      MemLp (fun x : ℝ => x) 2 μ ∧ MemLp (fun x : ℝ => x) 2 ν ∧
      Var[fun x : ℝ => x; μ] = 10 ∧ Var[fun x : ℝ => x; ν] = 10 ∧
      mUPM 2 (∫ x, x ∂μ) μ = 9 ∧ mLPM 2 (∫ x, x ∂μ) μ = 1 ∧
      mUPM 2 (∫ x, x ∂ν) ν = 1 ∧ mLPM 2 (∫ x, x ∂ν) ν = 9 :=
  exists_laws_variance_eq_pair_swapped

abbrev VarianceBlind : Prop :=
  (∃ μ ν : Measure ℝ, IsProbabilityMeasure μ ∧ IsProbabilityMeasure ν ∧
      MemLp (fun x : ℝ => x) 2 μ ∧ MemLp (fun x : ℝ => x) 2 ν ∧
      Var[fun x : ℝ => x; μ] = Var[fun x : ℝ => x; ν] ∧
      (mUPM 2 (∫ x, x ∂μ) μ, mLPM 2 (∫ x, x ∂μ) μ) ≠ (mUPM 2 (∫ x, x ∂ν) ν, mLPM 2 (∫ x, x ∂ν) ν)) ∧
    ¬ ∃ g : ℝ → ℝ × ℝ, ∀ μ : Measure ℝ, IsProbabilityMeasure μ → MemLp (fun x : ℝ => x) 2 μ →
      g Var[fun x : ℝ => x; μ] = (mUPM 2 (∫ x, x ∂μ) μ, mLPM 2 (∫ x, x ∂μ) μ)

theorem claim_C118 : VarianceBlind :=
  ⟨exists_laws_variance_eq_pair_ne, not_exists_recover_pair_of_variance⟩

theorem claim_C379 : VarianceBlind := claim_C118

abbrev TableAB : Prop :=
  ((mUPM 2 (∫ x, x ∂empMeasure tableA) (empMeasure tableA) = 9 ∧
      mLPM 2 (∫ x, x ∂empMeasure tableA) (empMeasure tableA) = 1 ∧
      Var[fun x : ℝ => x; empMeasure tableA] = 10) ∧
    (mUPM 2 (∫ x, x ∂empMeasure tableB) (empMeasure tableB) = 5 ∧
      mLPM 2 (∫ x, x ∂empMeasure tableB) (empMeasure tableB) = 5 ∧
      Var[fun x : ℝ => x; empMeasure tableB] = 10)) ∧
    (UPM 2 (mean tableA) tableA = 9 ∧ LPM 2 (mean tableA) tableA = 1 ∧ popVariance tableA = 10) ∧
    (UPM 2 (mean tableB) tableB = 5 ∧ LPM 2 (mean tableB) tableB = 5 ∧ popVariance tableB = 10)

theorem claim_C126 : TableAB := ⟨table_laws, table_A, table_B⟩

theorem claim_C169 : TableAB := claim_C126

theorem claim_C170 :
    TableAB ∧ ¬ ∃ g : ℝ → ℝ × ℝ, ∀ μ : Measure ℝ, IsProbabilityMeasure μ →
      MemLp (fun x : ℝ => x) 2 μ →
      g Var[fun x : ℝ => x; μ] = (mUPM 2 (∫ x, x ∂μ) μ, mLPM 2 (∫ x, x ∂μ) μ) :=
  ⟨claim_C126, not_exists_recover_pair_of_variance⟩

theorem claim_C399 :
    (∀ {Ω : Type*} {_ : MeasurableSpace Ω} (P : Measure Ω) [IsProbabilityMeasure P]
      (X : Ω → ℝ), P.map (fun ω => -X ω) = P.map X → MemLp X 4 P →
      0 < Var[fun ω => X ω ^ 2; P] →
      0 < Var[X; P] ∧ cov[X, fun ω => X ω ^ 2; P] = 0 ∧ corr X (fun ω => X ω ^ 2) P = 0) ∧
    corr (fun x : ℝ => x) (fun x => x ^ 2) (gaussianReal 0 1) = 0 :=
  ⟨fun _ _ X hsym hX hv => ⟨variance_pos_of_variance_sq_pos (hX.mono_exponent (by norm_num)) hv,
    covariance_sq_eq_zero_of_symm hsym hX, corr_sq_eq_zero_of_symm hsym hX⟩,
    corr_sq_gaussianReal⟩

theorem claim_C402 :
    (∃ X X' : Matrix (Fin 3) (Fin 2) ℝ, popCovMatrix X = popCovMatrix X' ∧
      (PMmatrix_v8898017 1 1 (colMeans X) X false false).clpm
        ≠ (PMmatrix_v8898017 1 1 (colMeans X') X' false false).clpm) ∧
    (∀ {n k : ℕ} (X : Matrix (Fin n) (Fin k) ℝ) (d : Fin k → ℝ), (∀ i, 0 < d i) →
      popCovMatrix (Matrix.of fun r j => d j * X r j)
        = (Matrix.of fun i j => d i * d j * popCovMatrix X i j) ∧
      corrMatrix (popCovMatrix (Matrix.of fun r j => d j * X r j))
        = corrMatrix (popCovMatrix X)) ∧
    (∀ i, 0 < popCovMatrix (!![1, 1; -1, -1] : Matrix (Fin 2) (Fin 2) ℝ) i i) ∧
    popCovMatrix (Matrix.of fun r j => 2 * (!![1, 1; -1, -1] : Matrix (Fin 2) (Fin 2) ℝ) r j)
      ≠ popCovMatrix (!![1, 1; -1, -1] : Matrix (Fin 2) (Fin 2) ℝ) := by
  have hdiag : ∀ i, popCovMatrix (!![1, 1; -1, -1] : Matrix (Fin 2) (Fin 2) ℝ) i i = 1 := by
    intro i
    fin_cases i <;> norm_num [popCovMatrix, covariance, coMoment, mean, colOf, Fin.sum_univ_two]
  have hraw : ∀ {m k : ℕ} (t : Fin k → ℝ) (X : Matrix (Fin m) (Fin k) ℝ),
      PMmatrix_v8898017 1 1 t X false false = PMmatrix 1 1 t X false false := fun _ _ => rfl
  simp only [hraw]
  refine ⟨exists_popCovMatrix_eq_clpm_ne, fun X d hd => ⟨popCovMatrix_scale X d, ?_⟩,
    fun i => by rw [hdiag]; norm_num, fun h => ?_⟩
  · rw [popCovMatrix_scale, corrMatrix_scale _ hd]
  · have h00 := congrFun (congrFun h 0) 0
    rw [popCovMatrix_scale _ fun _ => 2] at h00
    simp only [Matrix.of_apply, hdiag] at h00
    norm_num at h00

theorem claim_C381 :
    (∀ {Ω : Type*} {_ : MeasurableSpace Ω} (P : Measure Ω) [IsProbabilityMeasure P]
      (X Y : Ω → ℝ), MemLp X 2 P → MemLp Y 2 P → 0 < Var[X; P] → 0 < Var[Y; P] →
      corr X Y P = (∫ ω, (X ω - ∫ ω', X ω' ∂P) * (Y ω - ∫ ω', Y ω' ∂P) ∂P) /
        (Real.sqrt Var[X; P] * Real.sqrt Var[Y; P])) ∧
    (∀ (n : ℕ) (x y : Fin n → ℝ), 0 < popVariance x → 0 < popVariance y →
      sampleCorr x y = covariance x y / (Real.sqrt (popVariance x) * Real.sqrt (popVariance y)) ∧
      sampleCorr x y = corr x y (unifFin n)) :=
  ⟨fun _ _ _ _ _ _ _ _ => rfl, fun _ x y _ _ => ⟨rfl, sampleCorr_eq_corr x y⟩⟩

theorem claim_C384 :
    (∀ {Ω : Type*} {_ : MeasurableSpace Ω} (P : Measure Ω) [IsProbabilityMeasure P]
      (X Y : Ω → ℝ), MemLp X 2 P → MemLp Y 2 P → 0 < Var[X; P] → 0 < Var[Y; P] →
      cov[X, Y; P] = ∫ ω, (X ω - ∫ ω', X ω' ∂P) * (Y ω - ∫ ω', Y ω' ∂P) ∂P ∧
      corr X Y P = cov[X, Y; P] / (Real.sqrt Var[X; P] * Real.sqrt Var[Y; P]) ∧
      -1 ≤ corr X Y P ∧ corr X Y P ≤ 1) ∧
    (∀ (n : ℕ) (x y : Fin n → ℝ), 0 < popVariance x → 0 < popVariance y →
      -1 ≤ sampleCorr x y ∧ sampleCorr x y ≤ 1) :=
  ⟨fun _ _ _ _ hX hY _ _ => ⟨rfl, rfl, corr_mem_Icc hX hY⟩,
    fun _ x y _ _ => abs_le.mp (abs_sampleCorr_le_one x y)⟩

theorem claim_C398 {n k : ℕ} (X : Matrix (Fin n) (Fin k) ℝ) :
    (PMmatrix_v8898017 1 1 (colMeans X) X false false).cupm
        + (PMmatrix_v8898017 1 1 (colMeans X) X false false).clpm
        - (PMmatrix_v8898017 1 1 (colMeans X) X false false).dupm
        - (PMmatrix_v8898017 1 1 (colMeans X) X false false).dlpm = popCovMatrix X ∧
    (∀ i j, corrMatrix (popCovMatrix X) i j =
      popCovMatrix X i j / Real.sqrt (popCovMatrix X i i * popCovMatrix X j j)) ∧
    (∀ i j, corrMatrix (popCovMatrix X) i j = sampleCorr (colOf X i) (colOf X j)) ∧
    ((∀ i, 0 < popCovMatrix X i i) →
      (∀ i, corrMatrix (popCovMatrix X) i i = 1) ∧
      ∀ i j, -1 ≤ corrMatrix (popCovMatrix X) i j ∧ corrMatrix (popCovMatrix X) i j ≤ 1) := by
  have hraw : ∀ {m k : ℕ} (t : Fin k → ℝ) (X : Matrix (Fin m) (Fin k) ℝ),
      PMmatrix_v8898017 1 1 t X false false = PMmatrix 1 1 t X false false := fun _ _ => rfl
  simp only [hraw]
  refine ⟨by rw [← PMmatrix_cov_eq, PMmatrix_cov_eq_popCov], fun _ _ => rfl,
    corrMatrix_popCovMatrix X, fun hpos => ⟨fun i => ?_, fun i j => ?_⟩⟩
  · rw [corrMatrix_popCovMatrix]
    refine sampleCorr_self _ ?_
    have := hpos i
    simpa [popCovMatrix, sample_covariance_self] using this
  · rw [corrMatrix_popCovMatrix]
    exact abs_le.mp (abs_sampleCorr_le_one _ _)

end Projection

section Quantile

open GapQuantile

section RealDegreeThreshold

open GapOperators

variable {μ : Measure ℝ} {d : ℝ}

def rPMRatio (d t : ℝ) (μ : Measure ℝ) : ℝ := rlpm d t μ / (rlpm d t μ + rupm d t μ)

def rThreshold (d : ℝ) (μ : Measure ℝ) (a : ℝ) : ℝ := ginv (fun t => rPMRatio d t μ) a

lemma rThreshold_natCast (k : ℕ) (μ : Measure ℝ) (a : ℝ) :
    rThreshold k μ a = threshold k μ a := by
  simp only [rThreshold, threshold, rPMRatio, mPMRatio, rlpm_eq_mLPM, rupm_eq_mUPM]

lemma add_rpow_le_two_rpow_mul {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (hd : 0 ≤ d) :
    (a + b) ^ d ≤ 2 ^ d * (a ^ d + b ^ d) := by
  have hm : a + b ≤ 2 * max a b := by linarith [le_max_left a b, le_max_right a b]
  calc (a + b) ^ d ≤ (2 * max a b) ^ d := Real.rpow_le_rpow (by positivity) hm hd
    _ = 2 ^ d * max a b ^ d := Real.mul_rpow (by norm_num) (le_max_of_le_left ha)
    _ ≤ 2 ^ d * (a ^ d + b ^ d) := by
      gcongr
      rcases max_choice a b with h | h <;> rw [h]
      · linarith [Real.rpow_nonneg hb d]
      · linarith [Real.rpow_nonneg ha d]

lemma rDev_le_bound (hd : 0 ≤ d) {t c : ℝ} (htc : |t| ≤ c) (x : ℝ) :
    rLowerDev d t x ≤ 2 ^ d * (|x| ^ d + c ^ d) ∧
      rUpperDev d t x ≤ 2 ^ d * (|x| ^ d + c ^ d) := by
  have hc : 0 ≤ c := (abs_nonneg t).trans htc
  have h : |x - t| ^ d ≤ 2 ^ d * (|x| ^ d + c ^ d) :=
    (Real.rpow_le_rpow (abs_nonneg _) ((abs_sub x t).trans (by linarith)) hd).trans
      (add_rpow_le_two_rpow_mul (abs_nonneg x) hc hd)
  exact ⟨(GapCoMoments.rLowerDev_le_abs_rpow d t x).trans h,
    (GapCoMoments.rUpperDev_le_abs_rpow d t x).trans h⟩

lemma integrable_rLowerDev_of_rpow [IsFiniteMeasure μ] (hd : 0 ≤ d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) (t : ℝ) : Integrable (rLowerDev d t) μ :=
  ((hI.add (integrable_const _)).const_mul _).mono'
    (GapCoMoments.measurable_rLowerDev' d t).aestronglyMeasurable (ae_of_all _ fun x => by
      rw [Real.norm_of_nonneg (rLowerDev_nonneg d t x)]; exact (rDev_le_bound hd le_rfl x).1)

lemma integrable_rUpperDev_of_rpow [IsFiniteMeasure μ] (hd : 0 ≤ d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) (t : ℝ) : Integrable (rUpperDev d t) μ :=
  ((hI.add (integrable_const _)).const_mul _).mono'
    (GapCoMoments.measurable_rUpperDev' d t).aestronglyMeasurable (ae_of_all _ fun x => by
      rw [Real.norm_of_nonneg (rUpperDev_nonneg d t x)]; exact (rDev_le_bound hd le_rfl x).2)

lemma rLowerDev_mono_target (hd : 0 < d) (x : ℝ) : Monotone fun t => rLowerDev d t x := by
  intro s t hst
  simp only [rLowerDev_eq_max hd.ne']
  exact Real.rpow_le_rpow (le_max_right _ _) (max_le_max (by linarith) le_rfl) hd.le

lemma rUpperDev_anti_target (hd : 0 < d) (x : ℝ) : Antitone fun t => rUpperDev d t x := by
  intro s t hst
  simp only [rUpperDev_eq_max hd.ne']
  exact Real.rpow_le_rpow (le_max_right _ _) (max_le_max (by linarith) le_rfl) hd.le

lemma rlpm_monotone [IsFiniteMeasure μ] (hd : 0 < d) (hI : Integrable (fun x : ℝ => |x| ^ d) μ) :
    Monotone fun t => rlpm d t μ := fun s t hst =>
  integral_mono (integrable_rLowerDev_of_rpow hd.le hI s) (integrable_rLowerDev_of_rpow hd.le hI t)
    fun x => rLowerDev_mono_target hd x hst

lemma rupm_antitone [IsFiniteMeasure μ] (hd : 0 < d) (hI : Integrable (fun x : ℝ => |x| ^ d) μ) :
    Antitone fun t => rupm d t μ := fun s t hst =>
  integral_mono (integrable_rUpperDev_of_rpow hd.le hI t) (integrable_rUpperDev_of_rpow hd.le hI s)
    fun x => rUpperDev_anti_target hd x hst

lemma abs_le_of_mem_ball {t t₀ : ℝ} (ht : t ∈ Metric.ball t₀ 1) : |t| ≤ |t₀| + 1 := by
  rw [Metric.mem_ball, Real.dist_eq] at ht
  linarith [abs_sub_abs_le_abs_sub t t₀]

lemma continuous_rlpm [IsFiniteMeasure μ] (hd : 0 < d) (hI : Integrable (fun x : ℝ => |x| ^ d) μ) :
    Continuous fun t => rlpm d t μ := by
  refine continuous_iff_continuousAt.mpr fun t₀ => ?_
  change ContinuousAt (fun t => ∫ v, rLowerDev d t v ∂μ) t₀
  refine continuousAt_of_dominated (bound := fun x => 2 ^ d * (|x| ^ d + (|t₀| + 1) ^ d))
    (Eventually.of_forall fun t =>
      (GapCoMoments.measurable_rLowerDev' d t).aestronglyMeasurable) ?_
    ((hI.add (integrable_const _)).const_mul _) (ae_of_all _ fun x => ?_)
  · filter_upwards [Metric.ball_mem_nhds t₀ one_pos] with t ht
    refine ae_of_all _ fun x => ?_
    rw [Real.norm_of_nonneg (rLowerDev_nonneg d t x)]
    exact (rDev_le_bound hd.le (abs_le_of_mem_ball ht) x).1
  · simp only [rLowerDev_eq_max hd.ne']
    exact (((continuous_id.sub continuous_const).max continuous_const).rpow_const
      fun _ => Or.inr hd.le).continuousAt

lemma continuous_rupm [IsFiniteMeasure μ] (hd : 0 < d) (hI : Integrable (fun x : ℝ => |x| ^ d) μ) :
    Continuous fun t => rupm d t μ := by
  refine continuous_iff_continuousAt.mpr fun t₀ => ?_
  change ContinuousAt (fun t => ∫ v, rUpperDev d t v ∂μ) t₀
  refine continuousAt_of_dominated (bound := fun x => 2 ^ d * (|x| ^ d + (|t₀| + 1) ^ d))
    (Eventually.of_forall fun t =>
      (GapCoMoments.measurable_rUpperDev' d t).aestronglyMeasurable) ?_
    ((hI.add (integrable_const _)).const_mul _) (ae_of_all _ fun x => ?_)
  · filter_upwards [Metric.ball_mem_nhds t₀ one_pos] with t ht
    refine ae_of_all _ fun x => ?_
    rw [Real.norm_of_nonneg (rUpperDev_nonneg d t x)]
    exact (rDev_le_bound hd.le (abs_le_of_mem_ball ht) x).2
  · simp only [rUpperDev_eq_max hd.ne']
    exact (((continuous_const.sub continuous_id).max continuous_const).rpow_const
      fun _ => Or.inr hd.le).continuousAt

lemma tendsto_rlpm_atBot [IsFiniteMeasure μ] (hd : 0 < d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) : Tendsto (fun t => rlpm d t μ) atBot (𝓝 0) := by
  have := tendsto_integral_filter_of_dominated_convergence (μ := μ) (l := atBot)
    (F := fun t x => rLowerDev d t x) (f := fun _ => (0 : ℝ)) (rLowerDev d 0)
    (Eventually.of_forall fun t =>
      (GapCoMoments.measurable_rLowerDev' d t).aestronglyMeasurable) ?_
    (integrable_rLowerDev_of_rpow hd.le hI 0) (Eventually.of_forall fun x => ?_)
  · rw [integral_zero] at this; exact this
  · filter_upwards [eventually_le_atBot 0] with t ht
    refine Eventually.of_forall fun x => ?_
    rw [Real.norm_of_nonneg (rLowerDev_nonneg d t x)]
    exact rLowerDev_mono_target hd x ht
  · refine tendsto_const_nhds.congr' ?_
    filter_upwards [eventually_lt_atBot x] with t ht
    simp [rLowerDev, not_le.mpr ht]

lemma tendsto_rupm_atTop [IsFiniteMeasure μ] (hd : 0 < d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) : Tendsto (fun t => rupm d t μ) atTop (𝓝 0) := by
  have := tendsto_integral_filter_of_dominated_convergence (μ := μ) (l := atTop)
    (F := fun t x => rUpperDev d t x) (f := fun _ => (0 : ℝ)) (rUpperDev d 0)
    (Eventually.of_forall fun t =>
      (GapCoMoments.measurable_rUpperDev' d t).aestronglyMeasurable) ?_
    (integrable_rUpperDev_of_rpow hd.le hI 0) (Eventually.of_forall fun x => ?_)
  · rw [integral_zero] at this; exact this
  · filter_upwards [eventually_ge_atTop 0] with t ht
    refine Eventually.of_forall fun x => ?_
    rw [Real.norm_of_nonneg (rUpperDev_nonneg d t x)]
    exact rUpperDev_anti_target hd x ht
  · refine tendsto_const_nhds.congr' ?_
    filter_upwards [eventually_gt_atTop x] with t ht
    simp [rUpperDev, not_lt.mpr ht.le]

lemma exists_rupm_pos [IsProbabilityMeasure μ] (hd : 0 < d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) : ∃ t, 0 < rupm d t μ := by
  obtain ⟨t, ht⟩ := ((tendsto_cdf_atBot μ).eventually (gt_mem_nhds one_pos)).exists
  refine ⟨t, ?_⟩
  have hI' : μ (Ioi t) ≠ 0 := by
    intro h0
    rw [← compl_Iic, prob_compl_eq_zero_iff measurableSet_Iic] at h0
    rw [cdf_eq_real, measureReal_def, h0, ENNReal.toReal_one] at ht
    exact lt_irrefl _ ht
  change 0 < ∫ x, rUpperDev d t x ∂μ
  rw [integral_pos_iff_support_of_nonneg (rUpperDev_nonneg d t)
    (integrable_rUpperDev_of_rpow hd.le hI t)]
  refine (pos_iff_ne_zero.mpr hI').trans_le (measure_mono fun x hx => ?_)
  have hx' : t < x := hx
  simp only [Function.mem_support, rUpperDev, hx', ↓reduceIte]
  exact (Real.rpow_pos_of_pos (sub_pos.mpr hx') d).ne'

lemma exists_rlpm_pos [IsProbabilityMeasure μ] (hd : 0 < d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) : ∃ t, 0 < rlpm d t μ := by
  obtain ⟨t, ht⟩ := ((tendsto_cdf_atTop μ).eventually (lt_mem_nhds one_pos)).exists
  refine ⟨t + 1, ?_⟩
  have hI' : μ (Iic t) ≠ 0 := by
    intro h0
    rw [cdf_eq_real, measureReal_def, h0, ENNReal.toReal_zero] at ht
    exact lt_irrefl _ ht
  change 0 < ∫ x, rLowerDev d (t + 1) x ∂μ
  rw [integral_pos_iff_support_of_nonneg (rLowerDev_nonneg d _)
    (integrable_rLowerDev_of_rpow hd.le hI _)]
  refine (pos_iff_ne_zero.mpr hI').trans_le (measure_mono fun x hx => ?_)
  have hx' : x ≤ t := hx
  have hx1 : x ≤ t + 1 := by linarith
  simp only [Function.mem_support, rLowerDev, hx1, ↓reduceIte]
  exact (Real.rpow_pos_of_pos (by linarith) d).ne'

lemma rlpm_add_rupm_pos [IsFiniteMeasure μ] (hd : 0 < d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) {t : ℝ} (ht : μ {t}ᶜ ≠ 0) :
    0 < rlpm d t μ + rupm d t μ := by
  have h1 := integrable_rLowerDev_of_rpow hd.le hI t
  have h2 := integrable_rUpperDev_of_rpow hd.le hI t
  change 0 < (∫ x, rLowerDev d t x ∂μ) + ∫ x, rUpperDev d t x ∂μ
  rw [← integral_add h1 h2]
  simp_rw [rLowerDev_add_rUpperDev]
  rw [integral_pos_iff_support_of_nonneg (fun x => Real.rpow_nonneg (abs_nonneg _) d)
    ((h1.add h2).congr (ae_of_all _ fun x => rLowerDev_add_rUpperDev d t x))]
  refine (pos_iff_ne_zero.mpr ht).trans_le (measure_mono fun x hx => ?_)
  exact (Real.rpow_pos_of_pos (abs_pos.mpr (sub_ne_zero.mpr (Ne.symm hx))) d).ne'

lemma rPMRatio_nonneg (d t : ℝ) (μ : Measure ℝ) : 0 ≤ rPMRatio d t μ :=
  div_nonneg (rlpm_nonneg d t μ) (add_nonneg (rlpm_nonneg d t μ) (integral_nonneg
    (rUpperDev_nonneg d t)))

lemma rPMRatio_le_one (d t : ℝ) (μ : Measure ℝ) : rPMRatio d t μ ≤ 1 :=
  div_le_one_of_le₀ (le_add_of_nonneg_right (integral_nonneg (rUpperDev_nonneg d t)))
    (add_nonneg (rlpm_nonneg d t μ) (integral_nonneg (rUpperDev_nonneg d t)))

lemma rPMRatio_monotone [IsFiniteMeasure μ] (hd : 0 < d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) : Monotone fun t => rPMRatio d t μ := by
  intro s t hst
  have hL : rlpm d s μ ≤ rlpm d t μ := rlpm_monotone hd hI hst
  have hU : rupm d t μ ≤ rupm d s μ := rupm_antitone hd hI hst
  have hLs := rlpm_nonneg d s μ
  have hLt := rlpm_nonneg d t μ
  have hUs : 0 ≤ rupm d s μ := integral_nonneg (rUpperDev_nonneg d s)
  have hUt : 0 ≤ rupm d t μ := integral_nonneg (rUpperDev_nonneg d t)
  change rPMRatio d s μ ≤ rPMRatio d t μ
  rcases (add_nonneg hLs hUs).eq_or_lt with h0 | h0
  · rw [rPMRatio, ← h0, div_zero]; exact rPMRatio_nonneg d t μ
  rcases (add_nonneg hLt hUt).eq_or_lt with h1 | h1
  · have hs0 : rlpm d s μ = 0 := by linarith
    rw [rPMRatio, hs0, zero_div]; exact rPMRatio_nonneg d t μ
  rw [rPMRatio, rPMRatio, div_le_div_iff₀ h0 h1]
  nlinarith [mul_le_mul hL hU hUt hLt]

lemma continuous_rPMRatio [IsFiniteMeasure μ] (hd : 0 < d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) (hμ : ∀ t, μ {t}ᶜ ≠ 0) :
    Continuous fun t => rPMRatio d t μ :=
  (continuous_rlpm hd hI).div ((continuous_rlpm hd hI).add (continuous_rupm hd hI))
    fun t => (rlpm_add_rupm_pos hd hI (hμ t)).ne'

lemma tendsto_rPMRatio_atBot [IsProbabilityMeasure μ] (hd : 0 < d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) :
    Tendsto (fun t => rPMRatio d t μ) atBot (𝓝 0) := by
  obtain ⟨t₀, ht₀⟩ := exists_rupm_pos hd hI
  have hlim : Tendsto (fun t => rlpm d t μ / rupm d t₀ μ) atBot (𝓝 0) := by
    simpa using (tendsto_rlpm_atBot hd hI).div_const (rupm d t₀ μ)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hlim
    (Eventually.of_forall fun t => rPMRatio_nonneg d t μ) ?_
  filter_upwards [eventually_le_atBot t₀] with t ht
  exact div_le_div_of_nonneg_left (rlpm_nonneg d t μ) ht₀
    ((rupm_antitone hd hI ht).trans (le_add_of_nonneg_left (rlpm_nonneg d t μ)))

lemma tendsto_rPMRatio_atTop [IsProbabilityMeasure μ] (hd : 0 < d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) :
    Tendsto (fun t => rPMRatio d t μ) atTop (𝓝 1) := by
  obtain ⟨t₀, ht₀⟩ := exists_rlpm_pos hd hI
  have hlim : Tendsto (fun t => 1 - rupm d t μ / rlpm d t₀ μ) atTop (𝓝 1) := by
    simpa using ((tendsto_rupm_atTop hd hI).div_const (rlpm d t₀ μ)).const_sub 1
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' hlim tendsto_const_nhds ?_
    (Eventually.of_forall fun t => rPMRatio_le_one d t μ)
  filter_upwards [eventually_ge_atTop t₀] with t ht
  have hL : rlpm d t₀ μ ≤ rlpm d t μ := rlpm_monotone hd hI ht
  have hU : 0 ≤ rupm d t μ := integral_nonneg (rUpperDev_nonneg d t)
  have hS : 0 < rlpm d t μ + rupm d t μ := by linarith
  have h1 : rupm d t μ / (rlpm d t μ + rupm d t μ) ≤ rupm d t μ / rlpm d t₀ μ :=
    div_le_div_of_nonneg_left hU ht₀ (by linarith)
  have h2 : rPMRatio d t μ = 1 - rupm d t μ / (rlpm d t μ + rupm d t μ) := by
    rw [rPMRatio, eq_sub_iff_add_eq, ← add_div, div_self hS.ne']
  linarith

lemma rlpm_lt_rlpm [IsFiniteMeasure μ] (hd : 0 < d) (hI : Integrable (fun x : ℝ => |x| ^ d) μ)
    {s t : ℝ} (hst : s < t) (hs : 0 < rlpm d s μ) : rlpm d s μ < rlpm d t μ := by
  have hIs := integrable_rLowerDev_of_rpow hd.le hI s
  have hIt := integrable_rLowerDev_of_rpow hd.le hI t
  have hpos : 0 < μ (Iic s) := by
    refine pos_iff_ne_zero.mpr fun h0 => hs.ne' ?_
    refine integral_eq_zero_of_ae ?_
    filter_upwards [measure_eq_zero_iff_ae_notMem.mp h0] with x hx
    simp [rLowerDev, show ¬ x ≤ s from hx]
  have hsub : 0 < ∫ x, (rLowerDev d t x - rLowerDev d s x) ∂μ := by
    rw [integral_pos_iff_support_of_nonneg
      (fun x => sub_nonneg.mpr (rLowerDev_mono_target hd x hst.le)) (hIt.sub hIs)]
    refine hpos.trans_le (measure_mono fun x hx => ?_)
    have hx' : x ≤ s := hx
    have hlt : rLowerDev d s x < rLowerDev d t x := by
      simp only [rLowerDev, hx', hx'.trans hst.le, ↓reduceIte]
      exact Real.rpow_lt_rpow (sub_nonneg.mpr hx') (by linarith) hd
    exact (sub_pos.mpr hlt).ne'
  rw [integral_sub hIt hIs] at hsub
  exact sub_pos.mp hsub

lemma rPMRatio_lt_rPMRatio [IsFiniteMeasure μ] (hd : 0 < d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) {s t : ℝ} (hst : s < t)
    (hs0 : 0 < rPMRatio d s μ) (hs1 : rPMRatio d s μ < 1) :
    rPMRatio d s μ < rPMRatio d t μ := by
  have hLs0 := rlpm_nonneg d s μ
  have hUs0 : 0 ≤ rupm d s μ := integral_nonneg (rUpperDev_nonneg d s)
  have hLs : 0 < rlpm d s μ := by
    refine hLs0.lt_of_ne fun h0 => hs0.ne' ?_
    rw [rPMRatio, ← h0, zero_div]
  have hUs : 0 < rupm d s μ := by
    refine hUs0.lt_of_ne fun h0 => hs1.ne ?_
    rw [rPMRatio, ← h0, add_zero, div_self hLs.ne']
  have hL := rlpm_lt_rlpm hd hI hst hLs
  have hU : rupm d t μ ≤ rupm d s μ := rupm_antitone hd hI hst.le
  have hUt0 : 0 ≤ rupm d t μ := integral_nonneg (rUpperDev_nonneg d t)
  rw [rPMRatio, rPMRatio, div_lt_div_iff₀ (by linarith) (by linarith)]
  nlinarith [mul_le_mul_of_nonneg_left hU hLs0, mul_lt_mul_of_pos_right hL hUs]

theorem rThreshold_root [IsProbabilityMeasure μ] (hd : 0 < d)
    (hI : Integrable (fun x : ℝ => |x| ^ d) μ) (hμ : ∀ t, μ {t}ᶜ ≠ 0) {a : ℝ}
    (ha0 : 0 < a) (ha1 : a < 1) :
    rPMRatio d (rThreshold d μ a) μ = a ∧
      (∀ t, rPMRatio d t μ = a → t = rThreshold d μ a) ∧
      ∀ t, rThreshold d μ a ≤ t ↔ a ≤ rPMRatio d t μ := by
  have hmono := rPMRatio_monotone hd hI
  have hc := continuous_rPMRatio hd hI hμ
  have h0 := tendsto_rPMRatio_atBot hd hI
  have h1 := tendsto_rPMRatio_atTop hd hI
  have hr : rPMRatio d (rThreshold d μ a) μ = a := apply_ginv_eq hmono hc h0 h1 ha0 ha1
  refine ⟨hr, fun s hs => ?_, fun t => ginv_le_iff hmono (fun _ => hc.continuousWithinAt) h0 h1
    ha0 ha1⟩
  rcases lt_trichotomy s (rThreshold d μ a) with hlt | heq | hgt
  · have := rPMRatio_lt_rPMRatio hd hI hlt (hs.symm ▸ ha0) (hs.symm ▸ ha1)
    rw [hs, hr] at this
    exact absurd this (lt_irrefl _)
  · exact heq
  · have := rPMRatio_lt_rPMRatio hd hI hgt (hr.symm ▸ ha0) (hr.symm ▸ ha1)
    rw [hs, hr] at this
    exact absurd this (lt_irrefl _)

lemma integrable_empMeasure {n : ℕ} (hn : 0 < n) (x : Fin n → ℝ) (f : ℝ → ℝ) :
    Integrable f (empMeasure x) := by
  unfold empMeasure
  refine Integrable.smul_measure ?_ (ENNReal.inv_ne_top.mpr (by exact_mod_cast hn.ne'))
  exact integrable_finsetSum_measure.mpr fun i _ => integrable_dirac enorm_lt_top

theorem rThreshold_empMeasure {n : ℕ} (hd : 0 < d) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    {a : ℝ} (ha0 : 0 < a) (ha1 : a < 1) :
    rThreshold d (empMeasure x) a ∈ Icc (varMin x) (varMax x) ∧
      rLPM d (rThreshold d (empMeasure x) a) x /
        (rLPM d (rThreshold d (empMeasure x) a) x + rUPM d (rThreshold d (empMeasure x) a) x)
        = a := by
  obtain ⟨i, _, _⟩ := id hx
  have hn : 0 < n := Fin.pos i
  have := isProbabilityMeasure_empMeasure hn x
  have hI := integrable_empMeasure hn x fun v => |v| ^ d
  have hμ := empMeasure_compl_singleton_ne_zero hx
  have hL : ∀ t, rlpm d t (empMeasure x) = rLPM d t x := fun t => integral_empMeasure x _
  have hU : ∀ t, rupm d t (empMeasure x) = rUPM d t x := fun t => integral_empMeasure x _
  have hF : ∀ t, rPMRatio d t (empMeasure x) = rLPM d t x / (rLPM d t x + rUPM d t x) :=
    fun t => by rw [rPMRatio, hL, hU]
  obtain ⟨hr, -, hle⟩ := rThreshold_root hd hI hμ ha0 ha1
  refine ⟨⟨?_, ?_⟩, (hF _).symm.trans hr⟩
  · by_contra hlt
    rw [not_le] at hlt
    have h0 : rLPM d (rThreshold d (empMeasure x) a) x = 0 := by
      refine div_eq_zero_iff.mpr (Or.inl (Finset.sum_eq_zero fun j _ => ?_))
      simp [rLowerDev, not_le.mpr (hlt.trans_le (varMin_le x j))]
    rw [hF, h0, zero_div] at hr
    exact ha0.ne' hr.symm
  · refine (hle _).mpr ?_
    have hU0 : rUPM d (varMax x) x = 0 := by
      refine div_eq_zero_iff.mpr (Or.inl (Finset.sum_eq_zero fun j _ => ?_))
      simp [rUpperDev, not_lt.mpr (le_varMax x j)]
    have hpos := rlpm_add_rupm_pos hd hI (hμ (varMax x))
    rw [hL, hU, hU0, add_zero] at hpos
    rw [hF, hU0, add_zero, div_self hpos.ne']
    exact ha1.le

end RealDegreeThreshold

theorem claim_C363 :
    (∀ (d : ℕ) (μ : Measure ℝ) (α : ℝ), threshold d μ α = sInf {t | α ≤ mPMRatio d t μ}) ∧
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ] (α : ℝ), threshold 0 μ α = quantile μ α) ∧
    (∀ {n : ℕ} (x : Fin n → ℝ) (α : ℝ), 0 < n → 0 < α → α ≤ 1 →
      threshold 0 (empMeasure x) α = ecdfInv α x) ∧
    (LPMVaR0 (1 / 10) ![(0 : ℝ), 1, 2] = 1 / 5 ∧
      threshold 0 (empMeasure ![(0 : ℝ), 1, 2]) (1 / 10) = 0) ∧
    (∀ (d : ℕ) (μ : Measure ℝ) [IsProbabilityMeasure μ] (α : ℝ), d ≠ 0 →
      MemLp (fun x : ℝ => x) d μ → (∀ t, μ {t}ᶜ ≠ 0) → 0 < α → α < 1 →
      mPMRatio d (threshold d μ α) μ = α ∧
      (∀ t, mPMRatio d t μ = α → t = threshold d μ α) ∧
      ∀ t, threshold d μ α ≤ t ↔ α ≤ mPMRatio d t μ) ∧
    (∀ (d : ℕ) {n : ℕ} (x : Fin n → ℝ) (α : ℝ), d ≠ 0 → (∃ i j, x i ≠ x j) → 0 < α → α < 1 →
      threshold d (empMeasure x) α ∈ Icc (varMin x) (varMax x) ∧
      LPMratio d (threshold d (empMeasure x) α) x = α) ∧
    (∀ (d : ℕ) {n : ℕ} (x : Fin n → ℝ) (α : ℝ), 1 ≤ d → d ≤ 4 → (∃ i j, x i ≠ x j) →
      0 < α → α < 1 →
      LPMVaR d α x ∈ Icc (varMin x) (varMax x) ∧ |LPMratio d (LPMVaR d α x) x - α| ≤ rTol ∧
      |LPMVaR d α x - threshold d (empMeasure x) α| ≤
        uniqAt (shiftedArr x) n (locate d (shiftedArr x) n α) -
          uniqAt (shiftedArr x) n (locate d (shiftedArr x) n α - 1)) ∧
    (LPMVaR 4 (1 / 10 ^ 12) ![(0 : ℝ), 1] = 0 ∧
      LPMVaR 4 (1 / 10 ^ 12) ![(0 : ℝ), 1] ≠
        threshold 4 (empMeasure ![(0 : ℝ), 1]) (1 / 10 ^ 12)) ∧
    (∀ (d : ℕ) {n : ℕ} (x : Fin n → ℝ) (α : ℝ), 5 ≤ d → (∃ i j, x i ≠ x j) → 0 < α → α < 1 →
      LPMVaR d α x = threshold d (empMeasure x) α) ∧
    (∀ (d : ℝ) (μ : Measure ℝ) (α : ℝ), rThreshold d μ α = sInf {t | α ≤ rPMRatio d t μ}) ∧
    (∀ (k : ℕ) (μ : Measure ℝ) (α : ℝ), rThreshold k μ α = threshold k μ α) ∧
    (∀ (d : ℝ) (μ : Measure ℝ) [IsProbabilityMeasure μ] (α : ℝ), 0 < d →
      Integrable (fun x : ℝ => |x| ^ d) μ → (∀ t, μ {t}ᶜ ≠ 0) → 0 < α → α < 1 →
      rPMRatio d (rThreshold d μ α) μ = α ∧
      (∀ t, rPMRatio d t μ = α → t = rThreshold d μ α) ∧
      ∀ t, rThreshold d μ α ≤ t ↔ α ≤ rPMRatio d t μ) ∧
    (∀ (d : ℝ) {n : ℕ} (x : Fin n → ℝ) (α : ℝ), 0 < d → (∃ i j, x i ≠ x j) → 0 < α → α < 1 →
      rThreshold d (empMeasure x) α ∈ Icc (varMin x) (varMax x) ∧
      rLPM d (rThreshold d (empMeasure x) α) x /
        (rLPM d (rThreshold d (empMeasure x) α) x + rUPM d (rThreshold d (empMeasure x) α) x)
        = α) := by
  refine ⟨threshold_eq_sInf, fun _ _ α => threshold_zero α, fun x α hn h0 h1 => ?_, ?_,
    fun d μ _ α hd h hμ h0 h1 => ⟨mPMRatio_threshold hd h hμ h0 h1,
      fun _ ht => eq_threshold_of_mPMRatio_eq hd h hμ h0 h1 ht,
      fun _ => threshold_le_iff hd h hμ h0 h1⟩,
    fun d _ x α hd hx h0 h1 => threshold_empMeasure hd hx h0 h1,
    fun d _ x α hd1 hd4 hx h0 h1 => LPMVaR_near_threshold_of_le_four hd1 hd4 hx h0 h1,
    LPMVaR_ne_threshold_snap,
    fun d _ x α hd hx h0 h1 => LPMVaR_eq_threshold_of_five_le hd hx h0 h1,
    fun _ _ _ => rfl, rThreshold_natCast,
    fun d μ _ α hd hI hμ h0 h1 => rThreshold_root hd hI hμ h0 h1,
    fun d _ x α hd hx h0 h1 => rThreshold_empMeasure hd hx h0 h1⟩
  · have := isProbabilityMeasure_empMeasure hn x
    rw [threshold_zero, quantile_empMeasure hn x h0 h1]
  · refine ⟨LPMVaR0_not_minimal.1, ?_⟩
    have := isProbabilityMeasure_empMeasure (n := 3) (by norm_num) ![(0 : ℝ), 1, 2]
    have hc : ⌈((3 : ℕ) : ℝ) * (1 / 10)⌉₊ = 1 := by
      rw [Nat.ceil_eq_iff (by norm_num)]; norm_num
    rw [threshold_zero, quantile_empMeasure (by norm_num) _ (by norm_num) (by norm_num), ecdfInv,
      hc]
    change varMin ![(0 : ℝ), 1, 2] = 0
    obtain ⟨i, hi⟩ := exists_eq_varMin (by norm_num) ![(0 : ℝ), 1, 2]
    have h0 := varMin_le ![(0 : ℝ), 1, 2] 0
    rw [← hi] at h0 ⊢
    fin_cases i <;> simp at h0 ⊢ <;> norm_num at h0

theorem claim_C049 :
    (∀ (μ : Measure ℝ) (p : ℝ), quantile μ p = sInf {t | p ≤ cdf μ t}) ∧
    (∀ (μ : Measure ℝ) (p : ℝ), 0 < p → p < 1 →
      {t | p ≤ cdf μ t} = Ici (quantile μ p) ∧ ∀ t, quantile μ p ≤ t ↔ p ≤ cdf μ t) ∧
    (∀ {n : ℕ} (x : Fin n → ℝ) (p : ℝ), 0 < n → 0 < p → p ≤ 1 →
      quantile (empMeasure x) p = ecdfInv p x) :=
  ⟨quantile_eq_sInf, fun μ _ h0 h1 => ⟨setOf_le_cdf_eq_Ici μ h0 h1,
      fun _ => quantile_le_iff μ h0 h1⟩,
    fun x _ hn h0 h1 => quantile_empMeasure hn x h0 h1⟩

theorem claim_C352 :
    (∀ (μ : Measure ℝ) (p : ℝ), quantile μ p = sInf {t | p ≤ cdf μ t}) ∧
    (∀ (μ : Measure ℝ) (p : ℝ), 0 < p → p < 1 →
      {t | p ≤ cdf μ t} = Ici (quantile μ p) ∧ ∀ t, quantile μ p ≤ t ↔ p ≤ cdf μ t) ∧
    (∀ {n : ℕ} (x : Fin n → ℝ) (p : ℝ), 0 < n → 0 < p → p ≤ 1 →
      quantile (empMeasure x) p = ecdfInv p x) :=
  claim_C049

theorem claim_C353 :
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ] (α : ℝ), 0 < α → α < 1 →
      μ.real (Iio (quantile μ α)) ≤ α ∧ α ≤ μ.real (Iic (quantile μ α))) ∧
    (quantile (empMeasure ![(0 : ℝ), 0, 1]) (1 / 2) = 0 ∧
      (empMeasure ![(0 : ℝ), 0, 1]).real (Iio 0) = 0 ∧
      (empMeasure ![(0 : ℝ), 0, 1]).real (Iic 0) = 2 / 3) :=
  ⟨fun _ _ _ h0 h1 => measureReal_Iio_quantile_le_and_le_Iic h0 h1, tied_sample⟩

theorem claim_C356 :
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ] (α : ℝ), threshold 0 μ α = quantile μ α) ∧
    (∀ (μ : Measure ℝ) (α : ℝ), Continuous (cdf μ) → 0 < α → α < 1 →
      cdf μ (quantile μ α) = α) ∧
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ] (α : ℝ), 0 < α → α < 1 →
      μ.real (Iio (quantile μ α)) ≤ α ∧ α ≤ μ.real (Iic (quantile μ α))) ∧
    cdf (empMeasure ![(0 : ℝ), 0, 1]) (quantile (empMeasure ![(0 : ℝ), 0, 1]) (1 / 2)) = 2 / 3 :=
  ⟨fun _ _ α => threshold_zero α, fun μ _ hc h0 h1 => cdf_quantile_eq μ hc h0 h1,
    fun _ _ _ h0 h1 => measureReal_Iio_quantile_le_and_le_Iic h0 h1, cdf_quantile_tied_sample⟩

theorem claim_C330 :
    (∀ μ ν : Measure ℝ, ((∀ t, cdf μ t ≤ cdf ν t) ∧ ∃ t, cdf μ t < cdf ν t) ↔
      (∀ p ∈ Ioo (0 : ℝ) 1, quantile ν p ≤ quantile μ p) ∧
        ∃ p ∈ Ioo (0 : ℝ) 1, quantile ν p < quantile μ p) ∧
    (∀ {n m : ℕ} (x : Fin n → ℝ) (y : Fin m → ℝ), 0 < n → 0 < m →
      (((∀ t, LPM 0 t x ≤ LPM 0 t y) ∧ ∃ t, LPM 0 t x < LPM 0 t y) ↔
        (∀ p ∈ Ioo (0 : ℝ) 1, ecdfInv p y ≤ ecdfInv p x) ∧
          ∃ p ∈ Ioo (0 : ℝ) 1, ecdfInv p y < ecdfInv p x)) :=
  ⟨fsd_iff_quantile, fun x y hn hm => and_congr (forall_LPM_zero_le_iff_forall_ecdfInv_le hn hm x y)
    (exists_LPM_zero_lt_iff_exists_ecdfInv_lt hn hm x y)⟩

theorem claim_C057 :
    (∀ {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω) [IsProbabilityMeasure P] (X : Ω → ℝ),
      AEMeasurable X P → Continuous (cdf (P.map X)) →
      P.map (fun ω => cdf (P.map X) (X ω)) = volume.restrict (Icc 0 1)) ∧
    (bernoulliHalf.map (cdf bernoulliHalf) =
        (2 : ℝ≥0∞)⁻¹ • Measure.dirac (1 / 2) + (2 : ℝ≥0∞)⁻¹ • Measure.dirac 1 ∧
      bernoulliHalf.map (cdf bernoulliHalf) ≠ volume.restrict (Icc 0 1)) ∧
    (∀ {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω) [IsProbabilityMeasure P] (X V : Ω → ℝ),
      AEMeasurable X P → AEMeasurable V P → IndepFun X V P →
      P.map V = volume.restrict (Icc 0 1) →
      P.map (fun ω => randTransform (P.map X) (X ω, V ω)) = volume.restrict (Icc 0 1)) :=
  ⟨fun _ _ _ hX hc => map_cdf_comp_eq_uniform hX hc,
    ⟨map_cdf_bernoulliHalf, map_cdf_bernoulliHalf_ne_uniform⟩,
    fun _ _ _ _ hX hV hXV hVu => map_randTransform_comp_eq_uniform hX hV hXV hVu⟩

theorem claim_C204 :
    (∀ (μ : Measure (ℝ × ℝ)), HingeFinite μ → ∀ a b : ℝ,
      derivWithin (fun b' => derivWithin (fun a' => mCoLPM 1 1 a' b' μ) (Ici a) a) (Ici b) b
        = μ.real (Iic a ×ˢ Iic b)) ∧
    (∀ (μ ν : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν],
      HingeFinite μ → HingeFinite ν → (∀ a b, mCoLPM 1 1 a b μ = mCoLPM 1 1 a b ν) → μ = ν) ∧
    (∀ (μ ν : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν],
      (∀ a b, mCoLPM 0 0 a b μ = mCoLPM 0 0 a b ν) → μ = ν) ∧
    (∀ (μ ν : Measure ℝ) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν],
      (∀ t, Integrable (lowerDev 1 t) μ) → (∀ t, Integrable (lowerDev 1 t) ν) →
      (∀ t, mLPM 1 t μ = mLPM 1 t ν) → μ = ν) :=
  ⟨fun _ hH a b => hinge_recovers_joint_cdf hH a b, fun _ _ _ _ hμ hν h => ext_of_hinge hμ hν h,
    fun _ _ _ _ h => ext_of_mCoLPM_zero h, fun _ _ _ _ hμ hν h => eq_of_mLPM_one_eq hμ hν h⟩

end Quantile

section GaussianDensity

abbrev stdNormal : ℝ → ℝ := gaussianPDFReal 0 1

lemma stdNormal_eq (x : ℝ) : stdNormal x = stdNormal 0 * Real.exp (-(1 / 2) * x ^ 2) := by
  simp only [stdNormal, gaussianPDFReal, NNReal.coe_one, sub_zero, mul_one]
  rw [show (0 : ℝ) ^ 2 = 0 by norm_num, neg_zero, zero_div, Real.exp_zero, mul_one]
  congr 2
  ring

lemma stdNormal_le (x : ℝ) : stdNormal x ≤ stdNormal 0 := by
  rw [stdNormal_eq x]
  exact mul_le_of_le_one_right (gaussianPDFReal_nonneg 0 1 0)
    (Real.exp_le_one_iff.mpr (by nlinarith [sq_nonneg x]))

lemma continuous_stdNormal : Continuous stdNormal := by
  unfold stdNormal gaussianPDFReal
  fun_prop

lemma integrable_abs_mul_stdNormal : Integrable fun x => |x| * stdNormal x := by
  refine (((integrable_mul_exp_neg_mul_sq (b := 1 / 2) (by norm_num)).abs).const_mul
    (stdNormal 0)).congr (ae_of_all _ fun x => ?_)
  simp only
  rw [stdNormal_eq x, abs_mul, abs_of_pos (Real.exp_pos _)]
  ring

lemma integrable_hinge_mul_stdNormal (c : ℝ) :
    Integrable fun x => max (c - x) 0 * stdNormal x := by
  refine (((integrable_gaussianPDFReal 0 1).const_mul |c|).add
    integrable_abs_mul_stdNormal).mono'
    (((continuous_const.sub continuous_id).max continuous_const).mul
      continuous_stdNormal).aestronglyMeasurable (ae_of_all _ fun x => ?_)
  have h0 := gaussianPDFReal_nonneg 0 1 x
  rw [Real.norm_of_nonneg (mul_nonneg (le_max_right _ _) h0)]
  have : max (c - x) 0 ≤ |c| + |x| :=
    max_le (by linarith [le_abs_self c, neg_abs_le x]) (by positivity)
  change max (c - x) 0 * stdNormal x ≤ |c| * stdNormal x + |x| * stdNormal x
  nlinarith

theorem stdNormal_prod_hypotheses :
    Continuous (fun p : ℝ × ℝ => stdNormal p.1 * stdNormal p.2) ∧
    (∀ p : ℝ × ℝ, 0 ≤ stdNormal p.1 * stdNormal p.2) ∧
    ¬ HasCompactSupport (fun p : ℝ × ℝ => stdNormal p.1 * stdNormal p.2) ∧
    (∀ a b : ℝ, Integrable fun p : ℝ × ℝ =>
      stdNormal p.1 * stdNormal p.2 * (max (a - p.1) 0 * max (b - p.2) 0)) ∧
    (∀ b x₀ : ℝ, ∃ φ : ℝ → ℝ, Integrable φ ∧
      ∀ᶠ x in 𝓝 x₀, ∀ y, max (b - y) 0 * (stdNormal (x, y).1 * stdNormal (x, y).2) ≤ φ y) := by
  have hpos : ∀ x, 0 < stdNormal x := fun x => gaussianPDFReal_pos 0 1 x one_ne_zero
  refine ⟨(continuous_stdNormal.comp continuous_fst).mul (continuous_stdNormal.comp
      continuous_snd), fun p => (mul_pos (hpos _) (hpos _)).le, fun h => ?_, fun a b => ?_,
    fun b x₀ => ⟨fun y => stdNormal 0 * (max (b - y) 0 * stdNormal y),
      (integrable_hinge_mul_stdNormal b).const_mul _, Eventually.of_forall fun x y => ?_⟩⟩
  · have hs : Function.support (fun p : ℝ × ℝ => stdNormal p.1 * stdNormal p.2) = univ := by
      ext p
      simp [Function.mem_support, (hpos p.1).ne', (hpos p.2).ne']
    have hc : IsCompact (univ : Set (ℝ × ℝ)) := by
      simpa [HasCompactSupport, tsupport, hs] using h
    obtain ⟨M, hM⟩ := (hc.image continuous_fst).bddAbove
    have := hM (a := M + 1) ⟨(M + 1, 0), mem_univ _, rfl⟩
    linarith
  · rw [Measure.volume_eq_prod]
    refine ((integrable_hinge_mul_stdNormal a).mul_prod (integrable_hinge_mul_stdNormal b)).congr
      (ae_of_all _ fun p => ?_)
    simp only
    ring
  · change max (b - y) 0 * (stdNormal x * stdNormal y) ≤ stdNormal 0 * (max (b - y) 0 * stdNormal y)
    have := mul_le_mul_of_nonneg_right (stdNormal_le x)
      (mul_nonneg (le_max_right (b - y) 0) (hpos y).le)
    linarith [show max (b - y) 0 * (stdNormal x * stdNormal y) =
      stdNormal x * (max (b - y) 0 * stdNormal y) by ring]

end GaussianDensity

section CoMoments

open GapCoMoments

theorem claim_C240 (tx ty : ℝ) :
    (∀ μ : Measure (ℝ × ℝ),
      mCoLPM 1 1 tx ty μ = ∫ p, max (tx - p.1) 0 * max (ty - p.2) 0 ∂μ ∧
      mCoUPM 1 1 tx ty μ = ∫ p, max (p.1 - tx) 0 * max (p.2 - ty) 0 ∂μ ∧
      mDLPM 1 1 tx ty μ = ∫ p, max (p.1 - tx) 0 * max (ty - p.2) 0 ∂μ ∧
      mDUPM 1 1 tx ty μ = ∫ p, max (tx - p.1) 0 * max (p.2 - ty) 0 ∂μ) ∧
    ∀ {n : ℕ} (x y : Fin n → ℝ),
      CoLPM 1 1 tx ty x y = (∑ i, max (tx - x i) 0 * max (ty - y i) 0) / n ∧
      CoUPM 1 1 tx ty x y = (∑ i, max (x i - tx) 0 * max (y i - ty) 0) / n ∧
      DLPM 1 1 tx ty x y = (∑ i, max (x i - tx) 0 * max (ty - y i) 0) / n ∧
      DUPM 1 1 tx ty x y = (∑ i, max (tx - x i) 0 * max (y i - ty) 0) / n ∧
      mCoLPM 1 1 tx ty (pairEmp x y) = (∑ i, max (tx - x i) 0 * max (ty - y i) 0) / n :=
  ⟨fun μ => ⟨mCoLPM_one_one_eq_hinge μ tx ty, mCoUPM_one_one_eq_hinge μ tx ty,
      mDLPM_one_one_eq_hinge μ tx ty, mDUPM_one_one_eq_hinge μ tx ty⟩,
    fun x y => ⟨CoLPM_one_one_eq_hinge tx ty x y, CoUPM_one_one_eq_hinge tx ty x y,
      DLPM_one_one_eq_hinge tx ty x y, DUPM_one_one_eq_hinge tx ty x y, hinge_pairEmp tx ty x y⟩⟩

theorem claim_C405 :
    CoLPM 1 1 0 0 ch10X ch10Y = 7 / 5 ∧ CoUPM 1 1 0 0 ch10X ch10Y = 23 / 5 ∧
      DLPM 1 1 0 0 ch10X ch10Y = 0 ∧ DUPM 1 1 0 0 ch10X ch10Y = 0 :=
  ch10_worked_example

abbrev CoKernels : Prop :=
  ∀ {n : ℕ} (x y : Fin n → ℝ) (dx dy dl du tx ty : ℝ),
    (rCoLPM dx dy tx ty x y = (∑ i, lowKernel dx tx (x i) * lowKernel dy ty (y i)) / n ∧
      rCoUPM dx dy tx ty x y = (∑ i, upKernel dx tx (x i) * upKernel dy ty (y i)) / n ∧
      rDLPM dl du tx ty x y = (∑ i, upKernel du tx (x i) * lowKernel dl ty (y i)) / n ∧
      rDUPM dl du tx ty x y = (∑ i, lowKernel dl tx (x i) * upKernel du ty (y i)) / n) ∧
    ∀ a b : ℕ, rCoLPM a b tx ty x y = CoLPM a b tx ty x y ∧
      rCoUPM a b tx ty x y = CoUPM a b tx ty x y ∧ rDLPM a b tx ty x y = DLPM a b tx ty x y ∧
      rDUPM a b tx ty x y = DUPM a b tx ty x y

lemma coKernels : CoKernels := fun x y dx dy dl du tx ty =>
  ⟨rCo_eq_kernels dx dy dl du tx ty x y, fun a b => rCo_natCast a b tx ty x y⟩

theorem claim_C431 : CoKernels := coKernels

theorem claim_C432 : CoKernels := coKernels

theorem claim_C622 :
    (∀ {n : ℕ} (x y : Fin n → ℝ) (dl du tx ty : ℝ),
      rDLPM dl du tx ty x y = (∑ i, upKernel du tx (x i) * lowKernel dl ty (y i)) / n ∧
      ∀ a b : ℕ, rDLPM a b tx ty x y = DLPM a b tx ty x y) ∧
    DLPM 1 2 (-2) 2 ch10X ch10Y = 7 / 5 ∧ DLPM 2 1 (-2) 2 ch10X ch10Y = 11 / 5 :=
  ⟨fun x y dl du tx ty => ⟨(rCo_eq_kernels 0 0 dl du tx ty x y).2.2.1,
      fun a b => (rCo_natCast a b tx ty x y).2.2.1⟩, DLPM_degree_order⟩

theorem claim_C249 :
    (∀ (μ ν : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν],
      (∀ a b, mCoLPM 0 0 a b μ = mCoLPM 0 0 a b ν) → μ = ν) ∧
    (∀ (μ ν : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν],
      (∀ a b, mCoLPM 0 0 a b μ = mCoLPM 0 0 a b ν) → ∀ (r s : ℕ) (tx ty : ℝ),
      mCoLPM r s tx ty μ = mCoLPM r s tx ty ν ∧ mCoUPM r s tx ty μ = mCoUPM r s tx ty ν ∧
        mDLPM r s tx ty μ = mDLPM r s tx ty ν ∧ mDUPM r s tx ty μ = mDUPM r s tx ty ν) ∧
    (∀ (μ ν : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν],
      HingeFinite μ → HingeFinite ν → (∀ a b, mCoLPM 1 1 a b μ = mCoLPM 1 1 a b ν) → μ = ν) :=
  ⟨fun _ _ _ _ h => ext_of_mCoLPM_zero h, fun _ _ _ _ h => coPartialMoments_eq_of_degree_zero_eq h,
    fun _ _ _ _ hμ hν h => ext_of_hinge hμ hν h⟩

open Filter Topology Finset in
theorem claim_C247 {n : ℕ} (x y : Fin n → ℝ) (a b : ℝ) :
    mCoLPM 1 1 a b (pairEmp x y) = (∑ i, max (a - x i) 0 * max (b - y i) 0) / n ∧
    ∀ᶠ h in 𝓝[>] (0 : ℝ), (CoLPM 1 1 (a + h) (b + h) x y - CoLPM 1 1 (a + h) b x y -
      CoLPM 1 1 a (b + h) x y + CoLPM 1 1 a b x y) / h ^ 2 = CoLPM 0 0 a b x y ∧
      CoLPM 0 0 a b x y = (#{i | x i ≤ a ∧ y i ≤ b} : ℝ) / n := by
  refine ⟨hinge_pairEmp a b x y, ?_⟩
  filter_upwards [eventually_hinge_finite_difference a b x y] with h hh
  exact ⟨hh, CoLPM_zero_eq_card a b x y⟩

open Filter Topology in
theorem claim_C243 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] (a b : ℝ) :
    (∀ q : ℝ, Integrable (fun p : ℝ × ℝ => (|p.1 - a| * |p.2 - b|) ^ q) μ ↔
      Integrable (fun p : ℝ × ℝ => rLowerDev q a p.1 * rLowerDev q b p.2) μ ∧
      Integrable (fun p : ℝ × ℝ => rUpperDev q a p.1 * rUpperDev q b p.2) μ ∧
      Integrable (fun p : ℝ × ℝ => rUpperDev q a p.1 * rLowerDev q b p.2) μ ∧
      Integrable (fun p : ℝ × ℝ => rLowerDev q a p.1 * rUpperDev q b p.2) μ) ∧
    ((∃ ε : ℝ, 0 < ε ∧ Integrable (fun p : ℝ × ℝ => (|p.1 - a| * |p.2 - b|) ^ ε) μ) →
      Tendsto (fun q => quadLL q a b μ) (𝓝[>] 0) (𝓝 (μ.real {p | p.1 < a ∧ p.2 < b})) ∧
      Tendsto (fun q => quadUU q a b μ) (𝓝[>] 0) (𝓝 (μ.real {p | a < p.1 ∧ b < p.2})) ∧
      Tendsto (fun q => quadUL q a b μ) (𝓝[>] 0) (𝓝 (μ.real {p | a < p.1 ∧ p.2 < b})) ∧
      Tendsto (fun q => quadLU q a b μ) (𝓝[>] 0) (𝓝 (μ.real {p | p.1 < a ∧ b < p.2})) ∧
      Tendsto (fun q => quadTotal q a b μ) (𝓝[>] 0) (𝓝 (μ.real {p | p.1 ≠ a ∧ p.2 ≠ b}))) ∧
    μ.real {p : ℝ × ℝ | p.1 ≠ a ∧ p.2 ≠ b} = 1 - μ.real {p | p.1 = a ∨ p.2 = b} ∧
    ((∃ ε : ℝ, 0 < ε ∧ Integrable (fun p : ℝ × ℝ => (|p.1 - a| * |p.2 - b|) ^ ε) μ) →
      μ {p | p.1 = a} = 0 → μ {p | p.2 = b} = 0 →
      Tendsto (fun q => quadLL q a b μ / quadTotal q a b μ) (𝓝[>] 0)
        (𝓝 (μ.real {p | p.1 ≤ a ∧ p.2 ≤ b})) ∧
      Tendsto (fun q => quadUU q a b μ / quadTotal q a b μ) (𝓝[>] 0)
        (𝓝 (μ.real {p | a < p.1 ∧ b < p.2})) ∧
      Tendsto (fun q => quadUL q a b μ / quadTotal q a b μ) (𝓝[>] 0)
        (𝓝 (μ.real {p | a < p.1 ∧ p.2 ≤ b})) ∧
      Tendsto (fun q => quadLU q a b μ / quadTotal q a b μ) (𝓝[>] 0)
        (𝓝 (μ.real {p | p.1 ≤ a ∧ b < p.2}))) :=
  ⟨fun q => integrable_rpow_iff_quad q a b,
    fun ⟨_, hε, hμ⟩ => tendsto_quad_degree_zero_of_rpow hε hμ, off_axis_eq_one_sub a b,
    fun ⟨_, hε, hμ⟩ hax hay => tendsto_quad_weights_degree_zero_of_rpow hε hμ hax hay⟩

theorem claim_C246 :
    (∀ f : ℝ × ℝ → ℝ, Continuous f → (∀ p, 0 ≤ f p) →
      (∀ a b : ℝ, Integrable fun p : ℝ × ℝ => f p * (max (a - p.1) 0 * max (b - p.2) 0)) →
      (∀ b x₀ : ℝ, ∃ φ : ℝ → ℝ, Integrable φ ∧
        ∀ᶠ x in 𝓝 x₀, ∀ y, max (b - y) 0 * f (x, y) ≤ φ y) →
      (∃ H₁ H₂ H₃ : ℝ → ℝ → ℝ,
        (∀ a b, HasDerivAt (fun a' => mCoLPM 1 1 a' b
          (volume.withDensity fun p => ENNReal.ofReal (f p))) (H₁ a b) a) ∧
        (∀ a b, HasDerivAt (fun a' => H₁ a' b) (H₂ a b) a) ∧
        (∀ a b, HasDerivAt (fun b' => H₂ a b') (H₃ a b) b) ∧
        (∀ a b, HasDerivAt (fun b' => H₃ a b') (f (a, b)) b)) ∧
      ∀ a b : ℝ,
        deriv (fun b' => deriv (fun b'' => deriv (fun a' => deriv (fun a'' =>
          mCoLPM 1 1 a'' b'' (volume.withDensity fun p => ENNReal.ofReal (f p))) a') a) b') b
          = f (a, b)) ∧
    (∀ f : ℝ × ℝ → ℝ, Continuous f → HasCompactSupport f → (∀ p, 0 ≤ f p) →
      (∃ H₁ H₂ H₃ : ℝ → ℝ → ℝ,
        (∀ a b, HasDerivAt (fun a' => mCoLPM 1 1 a' b
          (volume.withDensity fun p => ENNReal.ofReal (f p))) (H₁ a b) a) ∧
        (∀ a b, HasDerivAt (fun a' => H₁ a' b) (H₂ a b) a) ∧
        (∀ a b, HasDerivAt (fun b' => H₂ a b') (H₃ a b) b) ∧
        (∀ a b, HasDerivAt (fun b' => H₃ a b') (f (a, b)) b)) ∧
      ∀ a b : ℝ,
        deriv (fun b' => deriv (fun b'' => deriv (fun a' => deriv (fun a'' =>
          mCoLPM 1 1 a'' b'' (volume.withDensity fun p => ENNReal.ofReal (f p))) a') a) b') b
          = f (a, b)) ∧
    (Continuous (fun p : ℝ × ℝ => stdNormal p.1 * stdNormal p.2) ∧
      (∀ p : ℝ × ℝ, 0 ≤ stdNormal p.1 * stdNormal p.2) ∧
      ¬ HasCompactSupport (fun p : ℝ × ℝ => stdNormal p.1 * stdNormal p.2) ∧
      (∀ a b : ℝ, Integrable fun p : ℝ × ℝ =>
        stdNormal p.1 * stdNormal p.2 * (max (a - p.1) 0 * max (b - p.2) 0)) ∧
      (∀ b x₀ : ℝ, ∃ φ : ℝ → ℝ, Integrable φ ∧
        ∀ᶠ x in 𝓝 x₀, ∀ y, max (b - y) 0 * (stdNormal (x, y).1 * stdNormal (x, y).2) ≤ φ y)) ∧
    ∀ a b : ℝ,
      deriv (fun b' => deriv (fun b'' => deriv (fun a' => deriv (fun a'' =>
        mCoLPM 1 1 a'' b'' (volume.withDensity fun p =>
          ENNReal.ofReal (stdNormal p.1 * stdNormal p.2))) a') a) b') b
        = stdNormal a * stdNormal b := by
  obtain ⟨hc, h0, -, hH, hdom⟩ := stdNormal_prod_hypotheses
  exact ⟨fun _ hf hf0 hH hdom => ⟨hinge_fourth_deriv_eq_density_of_dominated hf hf0 hH hdom,
      deriv_hinge_fourth_eq_density_of_dominated hf hf0 hH hdom⟩,
    fun _ hf hfc hf0 => ⟨hinge_fourth_deriv_eq_density hf hfc hf0,
      deriv_hinge_fourth_eq_density hf hfc hf0⟩, stdNormal_prod_hypotheses,
    deriv_hinge_fourth_eq_density_of_dominated hc h0 hH hdom⟩

open Finset in
theorem claim_C498 {n : ℕ} (x y : Fin n → ℝ) :
    (∀ u v : ℝ, CoLPM 0 0 u v (fun i => LPM 0 (x i) x) (fun i => LPM 0 (y i) y)
      = (#{i | (#{j | x j ≤ x i} : ℝ) / n ≤ u ∧ (#{j | y j ≤ y i} : ℝ) / n ≤ v} : ℝ) / n) ∧
    (∀ u v : ℝ, u < 1 / (n : ℝ) →
      CoLPM 0 0 u v (fun i => LPM 0 (x i) x) (fun i => LPM 0 (y i) y) = 0) ∧
    (∀ u v : ℝ, 1 ≤ u →
      CoLPM 0 0 u v (fun i => LPM 0 (x i) x) (fun i => LPM 0 (y i) y)
        = LPM 0 v fun i => LPM 0 (y i) y) ∧
    (0 < n → CoLPM 0 0 1 1 (fun i => LPM 0 (x i) x) (fun i => LPM 0 (y i) y) = 1) :=
  ⟨fun u v => CoLPM_zero_pseudo_eq_empirical_copula u v x y,
    fun _ v hu => empirical_copula_of_lt hu v x y,
    fun _ v hu => empirical_copula_of_one_le hu v x y,
    fun hn => empirical_copula_one_one hn x y⟩

theorem claim_C547 :
    (∃ (X : Matrix (Fin 3) (Fin 2) ℝ) (a b : Fin 2 → ℝ), (∀ j, a j < b j) ∧
      coLPMnDBatch_v8898017 X (Matrix.of fun _ : Fin 1 => b) 1 true 0
        - coLPMnDBatch_v8898017 X (Matrix.of fun _ : Fin 1 => ![a 0, b 1]) 1 true 0
        - coLPMnDBatch_v8898017 X (Matrix.of fun _ : Fin 1 => ![b 0, a 1]) 1 true 0
        + coLPMnDBatch_v8898017 X (Matrix.of fun _ : Fin 1 => a) 1 true 0 < 0) ∧
    (∃ (X : Matrix (Fin 3) (Fin 2) ℝ) (a b : Fin 2 → ℝ), (∀ j, a j < b j) ∧
      coLPMnDBatch_v8898017 (-X) (-(Matrix.of fun _ : Fin 1 => a)) 1 true 0
        - coLPMnDBatch_v8898017 (-X) (-(Matrix.of fun _ : Fin 1 => ![b 0, a 1])) 1 true 0
        - coLPMnDBatch_v8898017 (-X) (-(Matrix.of fun _ : Fin 1 => ![a 0, b 1])) 1 true 0
        + coLPMnDBatch_v8898017 (-X) (-(Matrix.of fun _ : Fin 1 => b)) 1 true 0 < 0) :=
  ⟨exists_negative_rectangle_batch_v8898017, exists_negative_survival_rectangle_v8898017⟩

end CoMoments

section Operators

open GapOperators Finset

variable {n : ℕ}

theorem claim_C000 : ∀ x : ℝ, x⁺ = max x 0 := posPart_eq_max

theorem claim_C001 : ∀ t x : ℝ,
    upperDev 1 t x = (x - t)⁺ ∧ (x - t)⁺ = max (x - t) 0 ∧
      lowerDev 1 t x = (t - x)⁺ ∧ (t - x)⁺ = max (t - x) 0 := fun t x =>
  ⟨(upperDev_one_eq_posPart t x).1, (upperDev_one_eq_posPart t x).2,
    (lowerDev_one_eq_posPart t x).1, (lowerDev_one_eq_posPart t x).2⟩

theorem claim_C107 : ∀ t x : ℝ,
    upperDev 1 t x = (x - t)⁺ ∧ (x - t)⁺ = max (x - t) 0 ∧
      lowerDev 1 t x = (t - x)⁺ ∧ (t - x)⁺ = max (t - x) 0 := claim_C001

theorem claim_C423 :
    (∀ x : ℝ, x⁺ = max x 0) ∧ ∀ (x : Fin n → ℝ) (i : Fin n), x⁺ i = max (x i) 0 :=
  ⟨posPart_eq_max, pmax_zero_apply⟩

theorem claim_C424 :
    (∀ (r t v : ℝ), 0 < r →
      rUpperDev r t v = (v - t)⁺ ^ r ∧ rLowerDev r t v = (t - v)⁺ ^ r) ∧
    (∀ (r t : ℝ) (x : Fin n → ℝ), 0 < r →
      rUPM r t x = (∑ i, (x i - t)⁺ ^ r) / n ∧ rLPM r t x = (∑ i, (t - x i)⁺ ^ r) / n) ∧
    (∀ t v : ℝ, rUpperDev 0 t v = (if t < v then 1 else 0) ∧
      rLowerDev 0 t v = (if v ≤ t then 1 else 0)) ∧
    (∀ t : ℝ, rUpperDev 0 t t = 0 ∧ rLowerDev 0 t t = 1) ∧
    (∀ (d : ℕ) (t : ℝ) (x : Fin n → ℝ), rUPM d t x = UPM d t x ∧ rLPM d t x = LPM d t x) ∧
    ∀ (t : ℝ) (x : Fin n → ℝ),
      UPM 1 t x = (∑ i, (x i - t)⁺) / n ∧ LPM 1 t x = (∑ i, (t - x i)⁺) / n := by
  refine ⟨fun r t v hr => ?_, fun r t x hr => ?_,
    fun t v => ⟨rUpperDev_zero t v, rLowerDev_zero t v⟩,
    fun t => ⟨by simp [rUpperDev_zero], by simp [rLowerDev_zero]⟩,
    fun d t x => ⟨rUPM_natCast d t x, rLPM_natCast d t x⟩, fun t x => ?_⟩
  · simp only [rUpperDev_eq_max hr.ne', rLowerDev_eq_max hr.ne', posPart_eq_max, and_self]
  · simp only [rUPM_eq_max hr.ne', rLPM_eq_max hr.ne', posPart_eq_max, and_self]
  · simp only [UPM, LPM, fun i => (upperDev_one_eq_posPart t (x i)).1,
      fun i => (lowerDev_one_eq_posPart t (x i)).1, and_self]

theorem claim_C014 :
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ] (t : ℝ),
      ∫ x, max (t - x) 0 ^ 0 ∂μ = 1 ∧ ∫ x, max (x - t) 0 ^ 0 ∂μ = 1 ∧
        mLPM 0 t μ + mUPM 0 t μ = 1) ∧
    (∀ (t : ℝ) (x : Fin n → ℝ), 0 < n →
      (∑ i, max (t - x i) 0 ^ 0) / n = 1 ∧ (∑ i, max (x i - t) 0 ^ 0) / n = 1 ∧
        LPM 0 t x + UPM 0 t x = 1) :=
  ⟨fun _ _ t => ⟨(literal_degree_zero_integral t).1, (literal_degree_zero_integral t).2,
      (literal_degree_zero_sum_integral t).2⟩,
    fun t x hn => ⟨(literal_degree_zero_sample hn t x).1, (literal_degree_zero_sample hn t x).2,
      (literal_degree_zero_sum_sample hn t x).2⟩⟩

theorem claim_C087 :
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ] (t : ℝ),
      mLPM 1 t μ = ∫ x, max (t - x) 0 ∂μ ∧
      mLPM 1 t μ = μ.real (Set.Iic t) * ∫ x, (t - x) ∂μ[|Set.Iic t] ∧
      mLPM 1 t μ ≤ ∫ x, (t - x) ∂μ[|Set.Iic t] ∧
      mLPM 0 t μ = μ.real (Set.Iic t) ∧
      (μ (Set.Iic t) ≠ 0 → mLPM 1 t μ / mLPM 0 t μ = ∫ x, (t - x) ∂μ[|Set.Iic t])) ∧
    (mLPM 1 1 (empMeasure ch07Sample) = 6 / 5 ∧
      ∫ v, (1 - v) ∂(empMeasure ch07Sample)[|Set.Iic 1] = 2) ∧
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ] (k : ℕ), MemLp (fun x : ℝ => x) k μ →
      ProbabilityTheory.centralMoment id k μ =
        mUPM k (∫ x, x ∂μ) μ + (-1) ^ k * mLPM k (∫ x, x ∂μ) μ) :=
  ⟨fun _ _ t => ⟨mLPM_one_eq_integral_max t _, mLPM_one_eq_real_mul_cond t,
      mLPM_one_le_cond t, mLPM_zero_eq_real t _, mLPM_one_div_mLPM_zero t⟩,
    ch07_L1_ne_conditional, fun _ _ _ h => centralMoment_eq_mUPM_add_mLPM h⟩

theorem claim_C088 :
    (∀ (r τ : ℝ) (μ : Measure ℝ), 0 < r → rlpm r τ μ = ∫ x, max (τ - x) 0 ^ r ∂μ) ∧
    (∀ (r τ : ℝ) (x : Fin n → ℝ), 0 < r → rLPM r τ x = (∑ i, max (τ - x i) 0 ^ r) / n) ∧
    (∀ (d : ℕ) (τ : ℝ) (μ : Measure ℝ), rlpm d τ μ = mLPM d τ μ) ∧
    (∀ (τ : ℝ) (μ : Measure ℝ), rlpm 0 τ μ = μ.real (Set.Iic τ)) ∧
    (rlpm 0 0 (Measure.dirac (1 : ℝ)) = 0 ∧
      ∫ x, max (0 - x) 0 ^ (0 : ℝ) ∂(Measure.dirac (1 : ℝ)) = 1) :=
  ⟨fun _ τ μ hr => rlpm_eq_integral_max_rpow hr τ μ, fun _ τ x hr => rLPM_eq_sum_max_rpow hr τ x,
    rlpm_eq_mLPM, rlpm_zero, rlpm_zero_ne_literal⟩

theorem claim_C104 :
    (∀ (r : ℕ) (t : ℝ) (μ : Measure ℝ),
      ProbabilityTheory.moment (fun x : ℝ => x - t) r μ = ∫ x, (x - t) ^ r ∂μ) ∧
    (∀ (r : ℕ) (t : ℝ) (μ : Measure ℝ) [IsFiniteMeasure μ], MemLp (fun x : ℝ => x) r μ →
      Integrable (fun x : ℝ => (x - t) ^ r) μ) ∧
    (∀ (r : ℕ) (μ : Measure ℝ),
      ∫ x, (x - ∫ y, y ∂μ) ^ r ∂μ = ProbabilityTheory.centralMoment id r μ) ∧
    (∀ (r : ℕ) (t : ℝ) (x : Fin n → ℝ),
      ∫ v, (v - t) ^ r ∂empMeasure x = (∑ i, (x i - t) ^ r) / n) := by
  refine ⟨fun _ _ _ => rfl, fun r t μ _ h => ?_, integral_sub_pow_mean_eq_centralMoment,
    fun r t x => by rw [integral_sub_pow_empMeasure, centralMoment_eq_sum]⟩
  refine (integrable_abs_sub_pow h t).mono'
    ((measurable_id.sub_const t).pow_const r).aestronglyMeasurable (ae_of_all _ fun x => ?_)
  rw [Real.norm_eq_abs, abs_pow]

theorem claim_C172 :
    (∀ (t : ℝ) (μ : Measure ℝ), mUPM 1 t μ = ∫ x, max (x - t) 0 ∂μ ∧
      mLPM 1 t μ = ∫ x, max (t - x) 0 ∂μ) ∧
    (∀ (t : ℝ) (x : Fin n → ℝ), UPM 1 t x = 1 / n * ∑ i, max (x i - t) 0 ∧
      LPM 1 t x = 1 / n * ∑ i, max (t - x i) 0) :=
  ⟨fun t μ => ⟨mUPM_one_eq_integral_max t μ, mLPM_one_eq_integral_max t μ⟩,
    fun t x => ⟨UPM_one_eq_sum_max t x, LPM_one_eq_sum_max t x⟩⟩

theorem claim_C173 :
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ] (t : ℝ),
      mUPM 1 t μ = μ.real (Set.Ioi t) * ∫ x, (x - t) ∂μ[|Set.Ioi t] ∧
      mLPM 1 t μ = μ.real (Set.Iic t) * ∫ x, (t - x) ∂μ[|Set.Iic t] ∧
      mUPM 1 t μ ≤ ∫ x, (x - t) ∂μ[|Set.Ioi t] ∧ mLPM 1 t μ ≤ ∫ x, (t - x) ∂μ[|Set.Iic t]) ∧
    (∀ (t : ℝ) (x : Fin n → ℝ), UPM 1 t x = 1 / n * ∑ i, max (x i - t) 0 ∧
      LPM 1 t x = 1 / n * ∑ i, max (t - x i) 0) ∧
    (mLPM 1 1 (empMeasure ch07Sample) = 6 / 5 ∧
      ∫ v, (1 - v) ∂(empMeasure ch07Sample)[|Set.Iic 1] = 2) :=
  ⟨fun _ _ t => ⟨mUPM_one_eq_real_mul_cond t, mLPM_one_eq_real_mul_cond t, mUPM_one_le_cond t,
      mLPM_one_le_cond t⟩, (claim_C172 (n := n)).2, ch07_L1_ne_conditional⟩

theorem claim_C214 :
    ∀ (t : ℝ) (x : Fin n → ℝ),
      LPM 1 t x = 1 / n * ∑ i, max (t - x i) 0 ∧ UPM 1 t x = 1 / n * ∑ i, max (x i - t) 0 ∧
      (0 < #{i | x i ≤ t} →
        LPM 1 t x = LPM 0 t x * ((∑ i with x i ≤ t, (t - x i)) / #{i | x i ≤ t})) ∧
      (0 < #{i | t < x i} →
        UPM 1 t x = UPM 0 t x * ((∑ i with t < x i, (x i - t)) / #{i | t < x i})) :=
  fun t x => ⟨LPM_one_eq_sum_max t x, UPM_one_eq_sum_max t x,
    LPM_one_eq_LPM_zero_mul_condMean t x, UPM_one_eq_UPM_zero_mul_condMean t x⟩

theorem claim_C196 :
    (∀ (r t : ℝ) (x : Fin n → ℝ), 0 < r →
      rLPM r t x = (∑ i, max (t - x i) 0 ^ r) / n ∧ rUPM r t x = (∑ i, max (x i - t) 0 ^ r) / n) ∧
    (∀ (t : ℝ) (x : Fin n → ℝ),
      rLPM 0 t x = (#{i | x i ≤ t} : ℝ) / n ∧ rUPM 0 t x = (#{i | t < x i} : ℝ) / n) ∧
    (∀ (d : ℕ) (t : ℝ) (x : Fin n → ℝ), rLPM d t x = LPM d t x ∧ rUPM d t x = UPM d t x) :=
  ⟨fun _ t x hr => ⟨rLPM_eq_max hr.ne' t x, rUPM_eq_max hr.ne' t x⟩,
    fun t x => ⟨rLPM_zero_eq_ecdf t x, rUPM_zero_eq_card t x⟩,
    fun d t x => ⟨rLPM_natCast d t x, rUPM_natCast d t x⟩⟩

theorem claim_C359 :
    (∀ (d t : ℝ) (μ : Measure ℝ), 1 ≤ d →
      rlpm d t μ = ∫ x, max (t - x) 0 ^ d ∂μ ∧ rupm d t μ = ∫ x, max (x - t) 0 ^ d ∂μ) ∧
    (∀ (d : ℕ) (t : ℝ) (μ : Measure ℝ), rlpm d t μ = mLPM d t μ ∧ rupm d t μ = mUPM d t μ) :=
  ⟨fun _ t μ hd => ⟨rlpm_eq_integral_max_rpow (by linarith) t μ,
      rupm_eq_integral_max_rpow (by positivity) t μ⟩,
    fun d t μ => ⟨rlpm_eq_mLPM d t μ, rupm_eq_mUPM d t μ⟩⟩

theorem claim_C216 :
    ∀ (d : Fin n → ℝ) (k : Fin n), (∀ j, 0 ≤ d j) → (∀ j, d j ≤ d k) → 0 < d k →
      MonotoneOn (fun r : ℝ => d k ^ r / ∑ j, d j ^ r) (Set.Ioi 0) ∧
      ((∃ j, 0 < d j ∧ d j < d k) →
        StrictMonoOn (fun r : ℝ => d k ^ r / ∑ j, d j ^ r) (Set.Ioi 0)) :=
  fun d k h0 hmax hk => ⟨share_monotoneOn d k h0 hmax hk, share_strictMonoOn d k h0 hmax⟩

theorem claim_C377 {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} {X : Ω → ℝ}
    [IsProbabilityMeasure P] (hX : MemLp X 2 P) :
    variance X P = ∫ ω, (X ω - ∫ ω', X ω' ∂P) ^ 2 ∂P :=
  variance_eq_integral_sq hX.aestronglyMeasurable.aemeasurable

theorem claim_C072 :
    (fun i => lowerDev 1 0 (ch04Sample i)) = ![4, 2, 1, 0, 0, 0] ∧ LPM 1 0 ch04Sample = 7 / 6 :=
  ⟨ch04_deviations, ch04_LPM_one⟩

theorem claim_C175 : LPM 1 1 ch07Sample = 6 / 5 ∧ UPM 1 1 ch07Sample = 6 / 5 :=
  ⟨ch07_LPM_one, ch07_UPM_one⟩

theorem claim_C176 :
    (∀ i, ch07Sample i ≠ 1) ∧ LPM 0 1 ch07Sample = 3 / 5 ∧ UPM 0 1 ch07Sample = 2 / 5 ∧
      LPM 1 1 ch07Sample = UPM 1 1 ch07Sample :=
  ⟨ch07_no_tie, ch07_LPM_zero, ch07_UPM_zero, ch07_LPM_one.trans ch07_UPM_one.symm⟩

theorem claim_C208 : LPM 0 1 ch08Sample = 3 / 5 := ch08_at_one.1

theorem claim_C209 : LPM 1 1 ch08Sample = 7 / 5 := ch08_at_one.2.1

theorem claim_C210 : UPM 1 1 ch08Sample = 4 / 5 := ch08_at_one.2.2

theorem claim_C211 :
    (∀ i, ch08Sample i ≠ 1) ∧ LPM 0 1 ch08Sample = 3 / 5 ∧ LPM 1 1 ch08Sample = 7 / 5 ∧
      UPM 1 1 ch08Sample = 4 / 5 :=
  ⟨ch08_no_tie_at_one, ch08_at_one⟩

theorem claim_C212 :
    LPM 0 0 ch08Sample = 3 / 5 ∧ LPM 1 0 ch08Sample = 4 / 5 ∧ UPM 1 0 ch08Sample = 6 / 5 :=
  ch08_at_zero

theorem claim_C213 :
    LPM 1 1 ch08Sample = 7 / 5 ∧ UPM 1 1 ch08Sample = 4 / 5 ∧ LPM 1 0 ch08Sample = 4 / 5 ∧
      UPM 1 0 ch08Sample = 6 / 5 ∧ UPM 1 1 ch08Sample < LPM 1 1 ch08Sample ∧
      LPM 1 0 ch08Sample < UPM 1 0 ch08Sample ∧
      ∀ t : ℝ, UPM 1 t ch08Sample < LPM 1 t ch08Sample ↔ 2 / 5 < t :=
  ⟨ch08_at_one.2.1, ch08_at_one.2.2, ch08_at_zero.2.1, ch08_at_zero.2.2, ch08_comparisons.1,
    ch08_comparisons.2, ch08_shortfall_dominates_iff⟩

end Operators

end

end NNS.GapClaims
