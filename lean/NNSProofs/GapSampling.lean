-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.MeasurePM
public import Mathlib.Probability.StrongLaw

namespace NNS

namespace GapSampling

open MeasureTheory ProbabilityTheory Set Filter Topology
open scoped ENNReal Function

@[expose] public noncomputable section

theorem le_event_eq_preimage {Ω : Type*} (X : Ω → ℝ) (t : ℝ) : {ω | X ω ≤ t} = X ⁻¹' Iic t :=
  rfl

theorem gt_event_eq_preimage {Ω : Type*} (X : Ω → ℝ) (t : ℝ) : {ω | t < X ω} = X ⁻¹' Ioi t :=
  rfl

variable {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} {X : Ω → ℝ}

theorem measurableSet_le_event (hX : Measurable X) (t : ℝ) : MeasurableSet {ω | X ω ≤ t} :=
  measurableSet_le hX measurable_const

theorem measurableSet_gt_event (hX : Measurable X) (t : ℝ) : MeasurableSet {ω | t < X ω} :=
  measurableSet_lt measurable_const hX

theorem prob_event_add_prob_compl [IsProbabilityMeasure P] {A : Set Ω} (hA : MeasurableSet A) :
    P.real A + P.real Aᶜ = 1 :=
  probReal_add_probReal_compl hA

theorem measure_add_measure_compl_eq_one [IsProbabilityMeasure P] {A : Set Ω}
    (hA : MeasurableSet A) : P A + P Aᶜ = 1 :=
  MeasureTheory.prob_add_prob_compl hA

theorem cdf_map_eq_prob [IsProbabilityMeasure P] (hX : AEMeasurable X P) (t : ℝ) :
    cdf (P.map X) t = P.real {ω | X ω ≤ t} := by
  have := (Measure.isProbabilityMeasure_map_iff hX).mpr ‹_›
  rw [cdf_eq_real, map_measureReal_apply_of_aemeasurable hX measurableSet_Iic]; rfl

theorem integral_law_eq_integral (hX : AEMeasurable X P) {g : ℝ → ℝ}
    (hg : AEStronglyMeasurable g (P.map X)) : ∫ x, g x ∂P.map X = ∫ ω, g (X ω) ∂P :=
  integral_map hX hg

theorem lintegral_law_eq_lintegral (hX : AEMeasurable X P) {g : ℝ → ℝ≥0∞}
    (hg : AEMeasurable g (P.map X)) : ∫⁻ x, g x ∂P.map X = ∫⁻ ω, g (X ω) ∂P :=
  lintegral_map' hg hX

theorem integrable_lowerDev_comp [IsFiniteMeasure P] {r : ℕ} (hX : AEMeasurable X P)
    (h : MemLp X r P) (t : ℝ) : Integrable (fun ω => lowerDev r t (X ω)) P :=
  ((integrable_map_measure (lowerDev_measurable r t).aestronglyMeasurable hX).mp
    (integrable_lowerDev ((memLp_map_measure_iff aestronglyMeasurable_id hX).mpr h) t))

theorem integrable_upperDev_comp [IsFiniteMeasure P] {r : ℕ} (hX : AEMeasurable X P)
    (h : MemLp X r P) (t : ℝ) : Integrable (fun ω => upperDev r t (X ω)) P :=
  ((integrable_map_measure (upperDev_measurable r t).aestronglyMeasurable hX).mp
    (integrable_upperDev ((memLp_map_measure_iff aestronglyMeasurable_id hX).mpr h) t))

theorem measurableSet_le_benchmark {T : Ω → ℝ} (hX : Measurable X) (hT : Measurable T) :
    MeasurableSet {ω | X ω ≤ T ω} :=
  measurableSet_le hX hT

theorem measurableSet_gt_benchmark {T : Ω → ℝ} (hX : Measurable X) (hT : Measurable T) :
    MeasurableSet {ω | T ω < X ω} :=
  measurableSet_lt hT hX

theorem measurable_lowerDev_benchmark {T : Ω → ℝ} (hX : Measurable X) (hT : Measurable T)
    (r : ℕ) : Measurable fun ω => lowerDev r (T ω) (X ω) :=
  (measurable_lowerDev_uncurry r).comp (hT.prodMk hX)

theorem measurable_upperDev_benchmark {T : Ω → ℝ} (hX : Measurable X) (hT : Measurable T)
    (r : ℕ) : Measurable fun ω => upperDev r (T ω) (X ω) :=
  (measurable_upperDev_uncurry r).comp (hT.prodMk hX)

theorem integrable_lowerDev_benchmark [IsFiniteMeasure P] {T : Ω → ℝ} {r : ℕ}
    (hX : Measurable X) (hT : Measurable T) (hXr : MemLp X r P) (hTr : MemLp T r P) :
    Integrable (fun ω => lowerDev r (T ω) (X ω)) P := by
  refine (integrable_lowerDev_comp (hX.sub hT).aemeasurable (hXr.sub hTr) 0).congr
    (Eventually.of_forall fun ω => ?_)
  simp [lowerDev]

theorem integrable_upperDev_benchmark [IsFiniteMeasure P] {T : Ω → ℝ} {r : ℕ}
    (hX : Measurable X) (hT : Measurable T) (hXr : MemLp X r P) (hTr : MemLp T r P) :
    Integrable (fun ω => upperDev r (T ω) (X ω)) P := by
  refine (integrable_upperDev_comp (hX.sub hT).aemeasurable (hXr.sub hTr) 0).congr
    (Eventually.of_forall fun ω => ?_)
  simp [upperDev]

theorem mLPM_frequentist {n : ℕ} (r : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    mLPM r t (empMeasure x) = LPM r t x :=
  mLPM_empMeasure r t x

theorem mUPM_frequentist {n : ℕ} (r : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    mUPM r t (empMeasure x) = UPM r t x :=
  mUPM_empMeasure r t x

section Bayes

variable {Θ : Type*} [MeasurableSpace Θ] {π : Measure Θ} {κ : Θ → Measure ℝ}

theorem isProbabilityMeasure_bind [IsProbabilityMeasure π] (hκ : AEMeasurable κ π)
    (h : ∀ θ, IsProbabilityMeasure (κ θ)) : IsProbabilityMeasure (π.bind κ) := by
  constructor
  rw [Measure.bind_apply .univ hκ]
  simp

theorem eLPM_bind (hκ : AEMeasurable κ π) (r : ℕ) (t : ℝ) :
    eLPM r t (π.bind κ) = ∫⁻ θ, eLPM r t (κ θ) ∂π :=
  Measure.lintegral_bind hκ (lowerDev_measurable r t).ennreal_ofReal.aemeasurable

theorem eUPM_bind (hκ : AEMeasurable κ π) (r : ℕ) (t : ℝ) :
    eUPM r t (π.bind κ) = ∫⁻ θ, eUPM r t (κ θ) ∂π :=
  Measure.lintegral_bind hκ (upperDev_measurable r t).ennreal_ofReal.aemeasurable

theorem mLPM_bind (hκ : AEMeasurable κ π) {r : ℕ} {t : ℝ}
    (h : Integrable (lowerDev r t) (π.bind κ)) :
    mLPM r t (π.bind κ) = ∫ θ, mLPM r t (κ θ) ∂π := by
  have hm : AEMeasurable (fun θ => eLPM r t (κ θ)) π :=
    (Measure.measurable_lintegral (lowerDev_measurable r t).ennreal_ofReal).comp_aemeasurable hκ
  have hfin : ∫⁻ θ, eLPM r t (κ θ) ∂π ≠ ∞ := by
    rw [← eLPM_bind hκ]; exact ((eLPM_lt_top_iff r t).mpr h).ne
  simp_rw [mLPM_eq_toReal]
  rw [integral_toReal hm (ae_lt_top' hm hfin), eLPM_bind hκ]

theorem mUPM_bind (hκ : AEMeasurable κ π) {r : ℕ} {t : ℝ}
    (h : Integrable (upperDev r t) (π.bind κ)) :
    mUPM r t (π.bind κ) = ∫ θ, mUPM r t (κ θ) ∂π := by
  have hm : AEMeasurable (fun θ => eUPM r t (κ θ)) π :=
    (Measure.measurable_lintegral (upperDev_measurable r t).ennreal_ofReal).comp_aemeasurable hκ
  have hfin : ∫⁻ θ, eUPM r t (κ θ) ∂π ≠ ∞ := by
    rw [← eUPM_bind hκ]; exact ((eUPM_lt_top_iff r t).mpr h).ne
  simp_rw [mUPM_eq_toReal]
  rw [integral_toReal hm (ae_lt_top' hm hfin), eUPM_bind hκ]

end Bayes

section LLN

variable {E : Type*} [MeasurableSpace E]

theorem ae_tendsto_average (W : ℕ → Ω → E) {g : E → ℝ} (hg : Measurable g)
    (hint : Integrable (fun ω => g (W 0 ω)) P) (hindep : Pairwise ((· ⟂ᵢ[P] ·) on W))
    (hident : ∀ i, IdentDistrib (W i) (W 0) P P) :
    ∀ᵐ ω ∂P, Tendsto (fun n : ℕ => (∑ i : Fin n, g (W i ω)) / n) atTop
      (𝓝 (∫ x, g x ∂P.map (W 0))) := by
  have h := strong_law_ae_real (fun i ω => g (W i ω)) hint
    (fun i j hij => (hindep hij).comp hg hg) (fun i => (hident i).comp hg)
  rw [integral_map (hident 0).aemeasurable_fst hg.aestronglyMeasurable]
  filter_upwards [h] with ω hω
  simpa only [Fin.sum_univ_eq_sum_range (fun i => g (W i ω))] using hω

theorem ae_tendsto_LPM (X : ℕ → Ω → ℝ) (r : ℕ) (t : ℝ)
    (hint : Integrable (fun ω => lowerDev r t (X 0 ω)) P)
    (hindep : Pairwise ((· ⟂ᵢ[P] ·) on X)) (hident : ∀ i, IdentDistrib (X i) (X 0) P P) :
    ∀ᵐ ω ∂P, Tendsto (fun n : ℕ => LPM r t fun i : Fin n => X i ω) atTop
      (𝓝 (mLPM r t (P.map (X 0)))) :=
  ae_tendsto_average X (lowerDev_measurable r t) hint hindep hident

variable (X Y : ℕ → Ω → ℝ) (r s : ℕ) (tx ty : ℝ)

end LLN

section GlivenkoCantelli

def ecdfLeft {n : ℕ} (t : ℝ) (x : Fin n → ℝ) : ℝ :=
  (∑ i, (Iio t).indicator (1 : ℝ → ℝ) (x i)) / n

private lemma LPM_zero_le_ecdfLeft {n : ℕ} {t u : ℝ} (htu : t < u) (x : Fin n → ℝ) :
    LPM 0 t x ≤ ecdfLeft u x := by
  unfold LPM ecdfLeft
  gcongr with i
  by_cases h : x i ≤ t
  · simp [lowerDev, h, lt_of_le_of_lt h htu]
  · simp only [lowerDev, h, ite_false]
    exact indicator_nonneg (fun a _ => by simp) _

private lemma LPM_zero_le_one {n : ℕ} (t : ℝ) (x : Fin n → ℝ) : LPM 0 t x ≤ 1 := by
  rw [LPM_zero_eq_ecdf]
  refine div_le_one_of_le₀ ?_ n.cast_nonneg
  exact_mod_cast (Finset.card_le_univ _).trans_eq (Fintype.card_fin n)

variable {μ : Measure ℝ} [IsProbabilityMeasure μ]

private lemma measureReal_Iio_le {q p : ℝ} (h : ∀ s < q, cdf μ s ≤ p) : μ.real (Iio q) ≤ p := by
  have hp : 0 ≤ p := (cdf_nonneg μ (q - 1)).trans (h _ (sub_one_lt q))
  have hlim : Function.leftLim (cdf μ) q ≤ p :=
    le_of_tendsto ((monotone_cdf μ).tendsto_leftLim q)
      (eventually_nhdsWithin_of_forall fun s hs => h s hs)
  rw [measureReal_def, ← measure_cdf μ, StieltjesFunction.measure_Iio _ (tendsto_cdf_atBot μ),
    sub_zero, ENNReal.toReal_ofReal']
  exact max_le hlim hp

private lemma exists_quantile (μ : Measure ℝ) {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) :
    ∃ q : ℝ, p ≤ cdf μ q ∧ ∀ s < q, cdf μ s < p := by
  set S : Set ℝ := {s | p ≤ cdf μ s} with hS
  have hne : S.Nonempty := ((tendsto_cdf_atTop μ).eventually (eventually_ge_nhds hp1)).exists
  have hbdd : BddBelow S := by
    obtain ⟨a, ha⟩ := eventually_atBot.mp ((tendsto_cdf_atBot μ).eventually (Iio_mem_nhds hp0))
    refine ⟨a, fun s hs => ?_⟩
    by_contra hlt
    exact absurd hs (not_le.mpr (ha s (not_le.mp hlt).le))
  refine ⟨sInf S, ?_, fun s hs => not_le.mp fun h => notMem_of_lt_csInf hs hbdd h⟩
  by_contra hlt
  have hev : ∀ᶠ s in 𝓝[≥] sInf S, cdf μ s < p :=
    ((cdf μ).right_continuous (sInf S)).eventually (Iio_mem_nhds (not_le.mp hlt))
  obtain ⟨u, hu, hIco⟩ := mem_nhdsGE_iff_exists_Ico_subset.mp hev
  obtain ⟨s, hsS, hsu⟩ := exists_lt_of_csInf_lt hne hu
  have hps : p ≤ cdf μ s := hsS
  have hsp : cdf μ s < p := hIco ⟨csInf_le hbdd hsS, hsu⟩
  exact absurd hps (not_le.mpr hsp)

theorem abs_ecdf_sub_cdf_le_of_grid {n : ℕ} (x : Fin n → ℝ) {k : ℕ} (hk : 0 < k) {δ : ℝ}
    (hδ : 0 ≤ δ) (q : ℕ → ℝ)
    (hq : ∀ j, 0 < j → j < k → (j : ℝ) / k ≤ cdf μ (q j) ∧ ∀ s < q j, cdf μ s < (j : ℝ) / k)
    (h : ∀ j, 0 < j → j < k → |LPM 0 (q j) x - cdf μ (q j)| ≤ δ)
    (hl : ∀ j, 0 < j → j < k → |ecdfLeft (q j) x - μ.real (Iio (q j))| ≤ δ) (t : ℝ) :
    |LPM 0 t x - cdf μ t| ≤ 1 / k + δ := by
  have hk0 : (0 : ℝ) < k := Nat.cast_pos.mpr hk
  have hF0 := cdf_nonneg μ t
  have hF1 := cdf_le_one μ t
  have hG0 := LPM_nonneg 0 t x
  have hG1 := LPM_zero_le_one t x
  have hup : ∀ j : ℕ, 0 < j → j < k → cdf μ t < (j : ℝ) / k → LPM 0 t x ≤ (j : ℝ) / k + δ := by
    intro j hj0 hjk hlt
    have htq : t < q j := not_le.mp fun hle =>
      absurd ((hq j hj0 hjk).1.trans (monotone_cdf μ hle)) (not_le.mpr hlt)
    have h1 := LPM_zero_le_ecdfLeft htq x
    have h2 := (abs_le.mp (hl j hj0 hjk)).2
    have h3 : μ.real (Iio (q j)) ≤ (j : ℝ) / k :=
      measureReal_Iio_le fun s hs => ((hq j hj0 hjk).2 s hs).le
    linarith
  have hlo : ∀ j : ℕ, 0 < j → j < k → (j : ℝ) / k ≤ cdf μ t → (j : ℝ) / k - δ ≤ LPM 0 t x := by
    intro j hj0 hjk hle
    have hqt : q j ≤ t := not_lt.mp fun hlt =>
      absurd hle (not_le.mpr ((hq j hj0 hjk).2 t hlt))
    have h1 := LPM_zero_monotone x hqt
    have h2 := (abs_le.mp (h j hj0 hjk)).1
    have h3 := (hq j hj0 hjk).1
    dsimp only at h1
    linarith
  have hkF : 0 ≤ k * cdf μ t := by positivity
  have hfl1 : (⌊k * cdf μ t⌋₊ : ℝ) ≤ k * cdf μ t := Nat.floor_le hkF
  have hfl2 : k * cdf μ t < ⌊k * cdf μ t⌋₊ + 1 := Nat.lt_floor_add_one _
  set m := ⌊k * cdf μ t⌋₊ with hm
  have hmk : (m : ℝ) / k ≤ cdf μ t := by rw [div_le_iff₀ hk0]; linarith
  have hmk' : cdf μ t < ((m + 1 : ℕ) : ℝ) / k := by
    rw [lt_div_iff₀ hk0]; push_cast; linarith
  have hsucc : ((m + 1 : ℕ) : ℝ) / k = (m : ℝ) / k + 1 / k := by push_cast; ring
  rw [abs_le]
  constructor
  · rcases Nat.eq_zero_or_pos m with hm0 | hm0
    · rw [hm0] at hmk' hsucc
      simp only [Nat.cast_zero, zero_div, zero_add] at hsucc
      linarith
    · rcases lt_or_ge m k with hmlt | hmge
      · have := hlo m hm0 hmlt hmk
        linarith
      · rcases Nat.lt_or_ge 1 k with hk2 | hk1
        · have hcast : ((k - 1 : ℕ) : ℝ) / k = 1 - 1 / k := by
            rw [Nat.cast_sub hk, Nat.cast_one, sub_div, div_self hk0.ne']
          have hle : ((k - 1 : ℕ) : ℝ) / k ≤ cdf μ t := by
            refine le_trans ?_ hmk
            gcongr
            exact_mod_cast (Nat.sub_le k 1).trans hmge
          have := hlo (k - 1) (Nat.sub_pos_of_lt hk2) (Nat.sub_lt hk one_pos) hle
          linarith
        · have hk1' : k = 1 := le_antisymm hk1 hk
          rw [hk1']
          simp only [Nat.cast_one, div_one]
          linarith
  · rcases lt_or_ge (m + 1) k with hmlt | hmge
    · have := hup (m + 1) (Nat.succ_pos m) hmlt hmk'
      linarith
    · have : (1 : ℝ) ≤ ((m + 1 : ℕ) : ℝ) / k := by
        rw [le_div_iff₀ hk0, one_mul]; exact_mod_cast hmge
      linarith

end GlivenkoCantelli

section GlivenkoCantelliAE

variable [IsProbabilityMeasure P] (X : ℕ → Ω → ℝ)

theorem ae_tendsto_ecdf (hindep : Pairwise ((· ⟂ᵢ[P] ·) on X))
    (hident : ∀ i, IdentDistrib (X i) (X 0) P P) (t : ℝ) :
    ∀ᵐ ω ∂P, Tendsto (fun n : ℕ => LPM 0 t fun i : Fin n => X i ω) atTop
      (𝓝 (cdf (P.map (X 0)) t)) := by
  have hX := (hident 0).aemeasurable_fst
  have := (Measure.isProbabilityMeasure_map_iff hX).mpr ‹_›
  rw [← mLPM_zero_eq_cdf]
  refine ae_tendsto_LPM X 0 t (integrable_lowerDev_comp hX ?_ t) hindep hident
  simpa using hX.aestronglyMeasurable

theorem ae_tendsto_ecdf_left (hindep : Pairwise ((· ⟂ᵢ[P] ·) on X))
    (hident : ∀ i, IdentDistrib (X i) (X 0) P P) (t : ℝ) :
    ∀ᵐ ω ∂P, Tendsto (fun n : ℕ => ecdfLeft t fun i : Fin n => X i ω) atTop
      (𝓝 ((P.map (X 0)).real (Iio t))) := by
  have hX := (hident 0).aemeasurable_fst
  have hg : Measurable ((Iio t).indicator (1 : ℝ → ℝ)) := measurable_one.indicator measurableSet_Iio
  rw [← integral_indicator_one measurableSet_Iio]
  refine ae_tendsto_average X hg ?_ hindep hident
  exact (integrable_map_measure hg.aestronglyMeasurable hX).mp
    ((integrable_const (1 : ℝ)).indicator measurableSet_Iio)

theorem glivenko_cantelli_eventually (hindep : Pairwise ((· ⟂ᵢ[P] ·) on X))
    (hident : ∀ i, IdentDistrib (X i) (X 0) P P) :
    ∀ᵐ ω ∂P, ∀ ε > 0, ∀ᶠ n : ℕ in atTop, ∀ t : ℝ,
      |LPM 0 t (fun i : Fin n => X i ω) - cdf (P.map (X 0)) t| ≤ ε := by
  have := (Measure.isProbabilityMeasure_map_iff (hident 0).aemeasurable_fst).mpr ‹_›
  have hex : ∀ k j : ℕ, ∃ q : ℝ, 0 < j → j < k →
      (j : ℝ) / k ≤ cdf (P.map (X 0)) q ∧ ∀ s < q, cdf (P.map (X 0)) s < (j : ℝ) / k := by
    intro k j
    by_cases hjk : 0 < j ∧ j < k
    · have hk0 : (0 : ℝ) < k := Nat.cast_pos.mpr (hjk.1.trans hjk.2)
      obtain ⟨q, hq⟩ := exists_quantile (P.map (X 0)) (p := (j : ℝ) / k)
        (div_pos (Nat.cast_pos.mpr hjk.1) hk0)
        ((div_lt_one hk0).mpr (Nat.cast_lt.mpr hjk.2))
      exact ⟨q, fun _ _ => hq⟩
    · exact ⟨0, fun h1 h2 => absurd ⟨h1, h2⟩ hjk⟩
  choose q hq using hex
  have h1 : ∀ᵐ ω ∂P, ∀ k j : ℕ, Tendsto (fun n : ℕ => LPM 0 (q k j) fun i : Fin n => X i ω)
      atTop (𝓝 (cdf (P.map (X 0)) (q k j))) :=
    ae_all_iff.mpr fun k => ae_all_iff.mpr fun j => ae_tendsto_ecdf X hindep hident (q k j)
  have h2 : ∀ᵐ ω ∂P, ∀ k j : ℕ, Tendsto (fun n : ℕ => ecdfLeft (q k j) fun i : Fin n => X i ω)
      atTop (𝓝 ((P.map (X 0)).real (Iio (q k j)))) :=
    ae_all_iff.mpr fun k => ae_all_iff.mpr fun j => ae_tendsto_ecdf_left X hindep hident (q k j)
  filter_upwards [h1, h2] with ω hω1 hω2 ε hε
  obtain ⟨k, hk⟩ := exists_nat_one_div_lt (half_pos hε)
  have e1 : ∀ᶠ n : ℕ in atTop, ∀ j ∈ Finset.range (k + 1),
      |LPM 0 (q (k + 1) j) (fun i : Fin n => X i ω) - cdf (P.map (X 0)) (q (k + 1) j)|
        ≤ ε / 2 :=
    (eventually_all_finset _).mpr fun j _ =>
      (Metric.tendsto_nhds.mp (hω1 (k + 1) j) (ε / 2) (half_pos hε)).mono fun n hn =>
        (Real.dist_eq _ _ ▸ hn).le
  have e2 : ∀ᶠ n : ℕ in atTop, ∀ j ∈ Finset.range (k + 1),
      |ecdfLeft (q (k + 1) j) (fun i : Fin n => X i ω) - (P.map (X 0)).real (Iio (q (k + 1) j))|
        ≤ ε / 2 :=
    (eventually_all_finset _).mpr fun j _ =>
      (Metric.tendsto_nhds.mp (hω2 (k + 1) j) (ε / 2) (half_pos hε)).mono fun n hn =>
        (Real.dist_eq _ _ ▸ hn).le
  filter_upwards [e1, e2] with n hn1 hn2 t
  have hgrid := abs_ecdf_sub_cdf_le_of_grid (μ := P.map (X 0)) (fun i : Fin n => X i ω)
    (Nat.succ_pos k) (half_pos hε).le (q (k + 1)) (hq (k + 1))
    (fun j _ hjk => hn1 j (Finset.mem_range.mpr hjk))
    (fun j _ hjk => hn2 j (Finset.mem_range.mpr hjk)) t
  push_cast at hgrid
  linarith

theorem glivenko_cantelli (hindep : Pairwise ((· ⟂ᵢ[P] ·) on X))
    (hident : ∀ i, IdentDistrib (X i) (X 0) P P) :
    ∀ᵐ ω ∂P, Tendsto
      (fun n : ℕ => ⨆ t : ℝ, |LPM 0 t (fun i : Fin n => X i ω) - cdf (P.map (X 0)) t|) atTop
      (𝓝 0) := by
  filter_upwards [glivenko_cantelli_eventually X hindep hident] with ω hω
  refine tendsto_order.2 ⟨fun a ha => Eventually.of_forall fun n =>
    lt_of_lt_of_le ha (Real.iSup_nonneg fun t => abs_nonneg _), fun a ha => ?_⟩
  filter_upwards [hω (a / 2) (half_pos ha)] with n hn
  exact lt_of_le_of_lt (ciSup_le hn) (half_lt_self ha)

theorem glivenko_cantelli_empMeasure (hindep : Pairwise ((· ⟂ᵢ[P] ·) on X))
    (hident : ∀ i, IdentDistrib (X i) (X 0) P P) :
    ∀ᵐ ω ∂P, Tendsto
      (fun n : ℕ => ⨆ t : ℝ,
        |cdf (empMeasure fun i : Fin (n + 1) => X i ω) t - cdf (P.map (X 0)) t|) atTop
      (𝓝 0) := by
  filter_upwards [glivenko_cantelli X hindep hident] with ω hω
  refine (hω.comp (tendsto_add_atTop_nat 1)).congr fun n => ?_
  exact iSup_congr fun t => by rw [cdf_empMeasure (Nat.succ_pos n)]

end GlivenkoCantelliAE

end

end GapSampling

end NNS
