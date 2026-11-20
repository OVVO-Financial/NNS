-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.MeasurePM
public import NNSProofs.VaRRoot

namespace NNS.GapQuantile

open MeasureTheory ProbabilityTheory Set Filter Topology
open scoped ENNReal Finset

@[expose] public noncomputable section

def ginv (G : ℝ → ℝ) (p : ℝ) : ℝ := sInf {t | p ≤ G t}

section ginv

variable {G : ℝ → ℝ} {p t : ℝ}

lemma ginv_set_nonempty (h1 : Tendsto G atTop (𝓝 1)) (hp1 : p < 1) : {t | p ≤ G t}.Nonempty :=
  (h1.eventually (eventually_ge_nhds hp1)).exists

lemma ginv_set_bddBelow (hG : Monotone G) (h0 : Tendsto G atBot (𝓝 0)) (hp0 : 0 < p) :
    BddBelow {t | p ≤ G t} := by
  obtain ⟨s, hs⟩ := (h0.eventually (eventually_lt_nhds hp0)).exists
  refine ⟨s, fun t ht => le_of_lt ?_⟩
  by_contra hts
  exact absurd ((hG (not_lt.mp hts)).trans_lt hs) (not_lt.mpr ht)

theorem le_apply_ginv (hG : Monotone G) (hr : ∀ x, ContinuousWithinAt G (Ici x) x)
    (h1 : Tendsto G atTop (𝓝 1)) (hp1 : p < 1) : p ≤ G (ginv G p) := by
  have hlim : Tendsto G (𝓝[>] (ginv G p)) (𝓝 (G (ginv G p))) :=
    ((hr _).mono Ioi_subset_Ici_self).tendsto
  refine ge_of_tendsto hlim ?_
  filter_upwards [self_mem_nhdsWithin] with t ht
  obtain ⟨s, hs, hst⟩ := exists_lt_of_csInf_lt (ginv_set_nonempty h1 hp1) ht
  exact le_trans hs (hG hst.le)

theorem apply_lt_of_lt_ginv (hG : Monotone G) (h0 : Tendsto G atBot (𝓝 0)) (hp0 : 0 < p)
    (ht : t < ginv G p) : G t < p := by
  by_contra h
  exact absurd (csInf_le (ginv_set_bddBelow hG h0 hp0) (not_lt.mp h)) (not_le.mpr ht)

theorem ginv_le_iff (hG : Monotone G) (hr : ∀ x, ContinuousWithinAt G (Ici x) x)
    (h0 : Tendsto G atBot (𝓝 0)) (h1 : Tendsto G atTop (𝓝 1)) (hp0 : 0 < p) (hp1 : p < 1) :
    ginv G p ≤ t ↔ p ≤ G t :=
  ⟨fun h => (le_apply_ginv hG hr h1 hp1).trans (hG h),
    fun h => csInf_le (ginv_set_bddBelow hG h0 hp0) h⟩

theorem apply_ginv_eq (hG : Monotone G) (hc : Continuous G) (h0 : Tendsto G atBot (𝓝 0))
    (h1 : Tendsto G atTop (𝓝 1)) (hp0 : 0 < p) (hp1 : p < 1) : G (ginv G p) = p := by
  refine le_antisymm ?_ (le_apply_ginv hG (fun x => hc.continuousWithinAt) h1 hp1)
  have hlim : Tendsto G (𝓝[<] (ginv G p)) (𝓝 (G (ginv G p))) := hc.continuousWithinAt.tendsto
  refine le_of_tendsto hlim ?_
  filter_upwards [self_mem_nhdsWithin] with t ht
  exact (apply_lt_of_lt_ginv hG h0 hp0 ht).le

end ginv

def quantile (μ : Measure ℝ) (p : ℝ) : ℝ := ginv (cdf μ) p

variable {μ ν : Measure ℝ} {p t : ℝ}

theorem quantile_eq_sInf (μ : Measure ℝ) (p : ℝ) : quantile μ p = sInf {t | p ≤ cdf μ t} := rfl

theorem quantile_set_nonempty (μ : Measure ℝ) (hp1 : p < 1) : {t | p ≤ cdf μ t}.Nonempty :=
  ginv_set_nonempty (tendsto_cdf_atTop μ) hp1

theorem quantile_set_bddBelow (μ : Measure ℝ) (hp0 : 0 < p) : BddBelow {t | p ≤ cdf μ t} :=
  ginv_set_bddBelow (monotone_cdf μ) (tendsto_cdf_atBot μ) hp0

theorem le_cdf_quantile (μ : Measure ℝ) (hp1 : p < 1) : p ≤ cdf μ (quantile μ p) :=
  le_apply_ginv (monotone_cdf μ) (cdf μ).right_continuous (tendsto_cdf_atTop μ) hp1

theorem quantile_le_iff (μ : Measure ℝ) (hp0 : 0 < p) (hp1 : p < 1) :
    quantile μ p ≤ t ↔ p ≤ cdf μ t :=
  ginv_le_iff (monotone_cdf μ) (cdf μ).right_continuous (tendsto_cdf_atBot μ)
    (tendsto_cdf_atTop μ) hp0 hp1

theorem cdf_lt_of_lt_quantile (μ : Measure ℝ) (hp0 : 0 < p) (ht : t < quantile μ p) :
    cdf μ t < p :=
  apply_lt_of_lt_ginv (monotone_cdf μ) (tendsto_cdf_atBot μ) hp0 ht

theorem setOf_le_cdf_eq_Ici (μ : Measure ℝ) (hp0 : 0 < p) (hp1 : p < 1) :
    {t | p ≤ cdf μ t} = Ici (quantile μ p) :=
  Set.ext fun _ => (quantile_le_iff μ hp0 hp1).symm

theorem quantile_empMeasure {n : ℕ} (hn : 0 < n) (x : Fin n → ℝ) (hp0 : 0 < p) (hp1 : p ≤ 1) :
    quantile (empMeasure x) p = ecdfInv p x := by
  have h : {t | p ≤ cdf (empMeasure x) t} = Ici (ecdfInv p x) := by
    ext t
    rw [mem_ofPred_eq, cdf_empMeasure hn, mem_Ici, ← le_LPMratio_zero_iff hp0 hp1 hn]
    simp [LPMratio]
  rw [quantile, ginv, h, csInf_Ici]

lemma measureReal_Iio_le [IsProbabilityMeasure μ] {q a : ℝ} (h : ∀ t < q, cdf μ t ≤ a) :
    μ.real (Iio q) ≤ a := by
  have ha : 0 ≤ a := (cdf_nonneg μ (q - 1)).trans (h _ (by linarith))
  have hU : Iio q = ⋃ n : ℕ, Iic (q - 1 / ((n : ℝ) + 1)) := by
    ext t
    simp only [mem_Iio, mem_iUnion, mem_Iic]
    constructor
    · intro ht
      obtain ⟨n, hn⟩ := exists_nat_one_div_lt (sub_pos.mpr ht)
      exact ⟨n, by linarith⟩
    · rintro ⟨n, hn⟩
      have : (0 : ℝ) < 1 / ((n : ℝ) + 1) := by positivity
      linarith
  have hmono : Monotone fun n : ℕ => Iic (q - 1 / ((n : ℝ) + 1)) := by
    intro m n hmn
    refine Iic_subset_Iic.mpr ?_
    have : 1 / ((n : ℝ) + 1) ≤ 1 / ((m : ℝ) + 1) :=
      one_div_le_one_div_of_le (by positivity) (by exact_mod_cast Nat.succ_le_succ hmn)
    linarith
  have hle : μ (Iio q) ≤ ENNReal.ofReal a := by
    rw [hU, hmono.measure_iUnion]
    refine iSup_le fun n => ?_
    rw [← ofReal_cdf]
    exact ENNReal.ofReal_le_ofReal (h _ (by
      have : (0 : ℝ) < 1 / ((n : ℝ) + 1) := by positivity
      linarith))
  exact ENNReal.toReal_le_of_le_ofReal ha hle

theorem measureReal_Iio_quantile_le [IsProbabilityMeasure μ] (hp0 : 0 < p) :
    μ.real (Iio (quantile μ p)) ≤ p :=
  measureReal_Iio_le fun _ ht => (cdf_lt_of_lt_quantile μ hp0 ht).le

theorem le_measureReal_Iic_quantile [IsProbabilityMeasure μ] (hp1 : p < 1) :
    p ≤ μ.real (Iic (quantile μ p)) := by
  rw [← cdf_eq_real]
  exact le_cdf_quantile μ hp1

theorem measureReal_Iio_quantile_le_and_le_Iic [IsProbabilityMeasure μ] (hp0 : 0 < p)
    (hp1 : p < 1) :
    μ.real (Iio (quantile μ p)) ≤ p ∧ p ≤ μ.real (Iic (quantile μ p)) :=
  ⟨measureReal_Iio_quantile_le hp0, le_measureReal_Iic_quantile hp1⟩

theorem tied_sample :
    quantile (empMeasure ![(0 : ℝ), 0, 1]) (1 / 2) = 0 ∧
      (empMeasure ![(0 : ℝ), 0, 1]).real (Iio 0) = 0 ∧
      (empMeasure ![(0 : ℝ), 0, 1]).real (Iic 0) = 2 / 3 := by
  have hn : 0 < 3 := by norm_num
  have := isProbabilityMeasure_empMeasure hn ![(0 : ℝ), 0, 1]
  have hcdf : ∀ t, cdf (empMeasure ![(0 : ℝ), 0, 1]) t =
      ((if 0 ≤ t then 1 else 0) + (if 0 ≤ t then 1 else 0) + (if 1 ≤ t then 1 else 0)) / 3 := by
    intro t
    rw [cdf_empMeasure hn, LPM, Fin.sum_univ_three]
    simp [lowerDev]
  refine ⟨le_antisymm ?_ ?_, ?_, ?_⟩
  · refine (quantile_le_iff _ (by norm_num) (by norm_num)).mpr ?_
    rw [hcdf]; norm_num
  · by_contra h
    have h1 := le_cdf_quantile (empMeasure ![(0 : ℝ), 0, 1]) (p := 1 / 2) (by norm_num)
    have h2 : ¬ (1 : ℝ) ≤ quantile (empMeasure ![(0 : ℝ), 0, 1]) (1 / 2) := fun h' =>
      h (zero_le_one.trans h')
    rw [hcdf] at h1
    simp only [h, h2, ↓reduceIte] at h1
    norm_num at h1
  · simp [empMeasure, measureReal_def, Fin.sum_univ_three]
  · rw [← cdf_eq_real, hcdf]; norm_num

theorem cdf_quantile_eq (μ : Measure ℝ) (hc : Continuous (cdf μ)) (hp0 : 0 < p) (hp1 : p < 1) :
    cdf μ (quantile μ p) = p :=
  apply_ginv_eq (monotone_cdf μ) hc (tendsto_cdf_atBot μ) (tendsto_cdf_atTop μ) hp0 hp1

theorem cdf_quantile_tied_sample :
    cdf (empMeasure ![(0 : ℝ), 0, 1]) (quantile (empMeasure ![(0 : ℝ), 0, 1]) (1 / 2)) = 2 / 3 := by
  have := isProbabilityMeasure_empMeasure (by norm_num : 0 < 3) ![(0 : ℝ), 0, 1]
  rw [cdf_eq_real, tied_sample.1, tied_sample.2.2]

theorem forall_cdf_le_iff_forall_quantile_le (μ ν : Measure ℝ) :
    (∀ t, cdf μ t ≤ cdf ν t) ↔ ∀ p ∈ Ioo (0 : ℝ) 1, quantile ν p ≤ quantile μ p := by
  constructor
  · intro h p hp
    exact (quantile_le_iff ν hp.1 hp.2).mpr ((le_cdf_quantile μ hp.2).trans (h _))
  · intro h t
    by_contra hlt
    have hlt' : cdf ν t < cdf μ t := not_le.mp hlt
    have hp0 : 0 < (cdf ν t + cdf μ t) / 2 := by have := cdf_nonneg ν t; linarith
    have hp1 : (cdf ν t + cdf μ t) / 2 < 1 := by have := cdf_le_one μ t; linarith
    have h1 : quantile μ ((cdf ν t + cdf μ t) / 2) ≤ t :=
      (quantile_le_iff μ hp0 hp1).mpr (by linarith)
    have h2 := (quantile_le_iff ν hp0 hp1).mp ((h _ ⟨hp0, hp1⟩).trans h1)
    linarith

theorem exists_cdf_lt_iff_exists_quantile_lt (μ ν : Measure ℝ) :
    (∃ t, cdf μ t < cdf ν t) ↔ ∃ p ∈ Ioo (0 : ℝ) 1, quantile ν p < quantile μ p := by
  constructor
  · rintro ⟨t, ht⟩
    have hp0 : 0 < (cdf μ t + cdf ν t) / 2 := by have := cdf_nonneg μ t; linarith
    have hp1 : (cdf μ t + cdf ν t) / 2 < 1 := by have := cdf_le_one ν t; linarith
    have h1 : quantile ν ((cdf μ t + cdf ν t) / 2) ≤ t :=
      (quantile_le_iff ν hp0 hp1).mpr (by linarith)
    refine ⟨_, ⟨hp0, hp1⟩, lt_of_le_of_lt h1 ?_⟩
    by_contra h
    have := (quantile_le_iff μ hp0 hp1).mp (not_lt.mp h)
    linarith
  · rintro ⟨p, hp, h⟩
    exact ⟨quantile ν p, (cdf_lt_of_lt_quantile μ hp.1 h).trans_le (le_cdf_quantile ν hp.2)⟩

theorem fsd_iff_quantile (μ ν : Measure ℝ) :
    ((∀ t, cdf μ t ≤ cdf ν t) ∧ ∃ t, cdf μ t < cdf ν t) ↔
      (∀ p ∈ Ioo (0 : ℝ) 1, quantile ν p ≤ quantile μ p) ∧
        ∃ p ∈ Ioo (0 : ℝ) 1, quantile ν p < quantile μ p :=
  and_congr (forall_cdf_le_iff_forall_quantile_le μ ν) (exists_cdf_lt_iff_exists_quantile_lt μ ν)

theorem forall_LPM_zero_le_iff_forall_ecdfInv_le {n m : ℕ} (hn : 0 < n) (hm : 0 < m)
    (x : Fin n → ℝ) (y : Fin m → ℝ) :
    (∀ t, LPM 0 t x ≤ LPM 0 t y) ↔ ∀ p ∈ Ioo (0 : ℝ) 1, ecdfInv p y ≤ ecdfInv p x := by
  have h := forall_cdf_le_iff_forall_quantile_le (empMeasure x) (empMeasure y)
  simp only [cdf_empMeasure hn, cdf_empMeasure hm] at h
  rw [h]
  refine forall₂_congr fun p hp => ?_
  rw [quantile_empMeasure hn x hp.1 hp.2.le, quantile_empMeasure hm y hp.1 hp.2.le]

theorem exists_LPM_zero_lt_iff_exists_ecdfInv_lt {n m : ℕ} (hn : 0 < n) (hm : 0 < m)
    (x : Fin n → ℝ) (y : Fin m → ℝ) :
    (∃ t, LPM 0 t x < LPM 0 t y) ↔ ∃ p ∈ Ioo (0 : ℝ) 1, ecdfInv p y < ecdfInv p x := by
  have h := exists_cdf_lt_iff_exists_quantile_lt (empMeasure x) (empMeasure y)
  simp only [cdf_empMeasure hn, cdf_empMeasure hm] at h
  rw [h]
  refine exists_congr fun p => and_congr_right fun hp => ?_
  rw [quantile_empMeasure hn x hp.1 hp.2.le, quantile_empMeasure hm y hp.1 hp.2.le]

def threshold (d : ℕ) (μ : Measure ℝ) (a : ℝ) : ℝ := ginv (fun t => mPMRatio d t μ) a

variable {d : ℕ} {a : ℝ}

theorem threshold_eq_sInf (d : ℕ) (μ : Measure ℝ) (a : ℝ) :
    threshold d μ a = sInf {t | a ≤ mPMRatio d t μ} := rfl

theorem threshold_zero [IsProbabilityMeasure μ] (a : ℝ) : threshold 0 μ a = quantile μ a := by
  simp only [threshold, quantile, mPMRatio_zero_eq_cdf]

theorem mPMRatio_threshold [IsProbabilityMeasure μ] (hd : d ≠ 0)
    (h : MemLp (fun x : ℝ => x) d μ) (hμ : ∀ t, μ {t}ᶜ ≠ 0) (ha0 : 0 < a) (ha1 : a < 1) :
    mPMRatio d (threshold d μ a) μ = a :=
  apply_ginv_eq (G := fun t => mPMRatio d t μ) (mPMRatio_monotone h)
    (continuous_mPMRatio hd h hμ) (tendsto_mPMRatio_atBot h) (tendsto_mPMRatio_atTop h) ha0 ha1

theorem threshold_le_iff [IsProbabilityMeasure μ] (hd : d ≠ 0)
    (h : MemLp (fun x : ℝ => x) d μ) (hμ : ∀ t, μ {t}ᶜ ≠ 0) (ha0 : 0 < a) (ha1 : a < 1) :
    threshold d μ a ≤ t ↔ a ≤ mPMRatio d t μ :=
  ginv_le_iff (G := fun t => mPMRatio d t μ) (mPMRatio_monotone h)
    (fun _ => (continuous_mPMRatio hd h hμ).continuousWithinAt) (tendsto_mPMRatio_atBot h)
    (tendsto_mPMRatio_atTop h) ha0 ha1

theorem threshold_balance [IsProbabilityMeasure μ] (hd : d ≠ 0)
    (h : MemLp (fun x : ℝ => x) d μ) (hμ : ∀ t, μ {t}ᶜ ≠ 0) (ha0 : 0 < a) (ha1 : a < 1) :
    (1 - a) * mLPM d (threshold d μ a) μ = a * mUPM d (threshold d μ a) μ := by
  have hpos := mLPM_add_mUPM_pos h (hμ (threshold d μ a))
  have hr := mPMRatio_threshold hd h hμ ha0 ha1
  rw [mPMRatio, div_eq_iff hpos.ne'] at hr
  linarith

lemma mLPM_lt_mLPM [IsFiniteMeasure μ] (hd : d ≠ 0) (h : MemLp (fun x : ℝ => x) d μ)
    {s t : ℝ} (hst : s < t) (hs : 0 < mLPM d s μ) : mLPM d s μ < mLPM d t μ := by
  have hI := integrable_lowerDev h
  have hpos : 0 < μ (Iic s) := by
    refine pos_iff_ne_zero.mpr fun h0 => hs.ne' ?_
    rw [mLPM_eq_setIntegral, Measure.restrict_eq_zero.mpr h0, integral_zero_measure]
  have hsub : 0 < ∫ x, (lowerDev d t x - lowerDev d s x) ∂μ := by
    rw [integral_pos_iff_support_of_nonneg
      (fun x => sub_nonneg.mpr (lowerDev_mono d x hst.le)) ((hI t).sub (hI s))]
    refine hpos.trans_le (measure_mono fun x hx => ?_)
    have hx' : x ≤ s := hx
    have hlt : lowerDev d s x < lowerDev d t x := by
      simp only [lowerDev, hx', hx'.trans hst.le, ↓reduceIte]
      exact pow_lt_pow_left₀ (by linarith) (sub_nonneg.mpr hx') hd
    exact (sub_pos.mpr hlt).ne'
  rw [integral_sub (hI t) (hI s)] at hsub
  exact sub_pos.mp hsub

theorem mPMRatio_lt_mPMRatio [IsFiniteMeasure μ] (hd : d ≠ 0) (h : MemLp (fun x : ℝ => x) d μ)
    {s t : ℝ} (hst : s < t) (hs0 : 0 < mPMRatio d s μ) (hs1 : mPMRatio d s μ < 1) :
    mPMRatio d s μ < mPMRatio d t μ := by
  have hLs0 := mLPM_nonneg d s μ
  have hUs0 := mUPM_nonneg d s μ
  have hUt0 := mUPM_nonneg d t μ
  have hLs : 0 < mLPM d s μ := by
    refine hLs0.lt_of_ne fun h0 => hs0.ne' ?_
    rw [mPMRatio, ← h0, zero_div]
  have hUs : 0 < mUPM d s μ := by
    refine hUs0.lt_of_ne fun h0 => hs1.ne ?_
    rw [mPMRatio, ← h0, add_zero, div_self hLs.ne']
  have hL := mLPM_lt_mLPM hd h hst hLs
  have hU : mUPM d t μ ≤ mUPM d s μ := mUPM_antitone h hst.le
  rw [mPMRatio, mPMRatio, div_lt_div_iff₀ (by linarith) (by linarith)]
  nlinarith [mul_le_mul_of_nonneg_left hU hLs0, mul_lt_mul_of_pos_right hL hUs]

theorem existsUnique_mPMRatio_eq [IsProbabilityMeasure μ] (hd : d ≠ 0)
    (h : MemLp (fun x : ℝ => x) d μ) (hμ : ∀ t, μ {t}ᶜ ≠ 0) (ha0 : 0 < a) (ha1 : a < 1) :
    ∃! t, mPMRatio d t μ = a := by
  have hr := mPMRatio_threshold hd h hμ ha0 ha1
  refine ⟨threshold d μ a, hr, fun s hs => ?_⟩
  rcases lt_trichotomy s (threshold d μ a) with hlt | heq | hgt
  · have := mPMRatio_lt_mPMRatio hd h hlt (hs.symm ▸ ha0) (hs.symm ▸ ha1)
    rw [hs, hr] at this
    exact absurd this (lt_irrefl _)
  · exact heq
  · have := mPMRatio_lt_mPMRatio hd h hgt (hr.symm ▸ ha0) (hr.symm ▸ ha1)
    rw [hs, hr] at this
    exact absurd this (lt_irrefl _)

theorem eq_threshold_of_mPMRatio_eq [IsProbabilityMeasure μ] (hd : d ≠ 0)
    (h : MemLp (fun x : ℝ => x) d μ) (hμ : ∀ t, μ {t}ᶜ ≠ 0) (ha0 : 0 < a) (ha1 : a < 1)
    (ht : mPMRatio d t μ = a) : t = threshold d μ a :=
  (existsUnique_mPMRatio_eq hd h hμ ha0 ha1).unique ht (mPMRatio_threshold hd h hμ ha0 ha1)

lemma empMeasure_compl_singleton_ne_zero {n : ℕ} {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    (t : ℝ) : empMeasure x {t}ᶜ ≠ 0 := by
  obtain ⟨i, j, hij⟩ := hx
  have hne : ∃ k, x k ≠ t := by
    by_cases h : x i = t
    · exact ⟨j, fun e => hij (h.trans e.symm)⟩
    · exact ⟨i, h⟩
  obtain ⟨k, hk⟩ := hne
  have hle : (n : ℝ≥0∞)⁻¹ ≤ empMeasure x {t}ᶜ := by
    rw [empMeasure, Measure.smul_apply, Measure.coe_finsetSum, Finset.sum_apply, smul_eq_mul]
    refine le_mul_of_one_le_right' ((Finset.single_le_sum (f := fun i => Measure.dirac (x i) {t}ᶜ)
      (fun _ _ => zero_le) (Finset.mem_univ k)).trans' ?_)
    simp [hk]
  exact (lt_of_lt_of_le (ENNReal.inv_pos.mpr (ENNReal.natCast_ne_top n)) hle).ne'

theorem mPMRatio_empMeasure {n : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    (t : ℝ) : mPMRatio d t (empMeasure x) = LPMratio d t x := by
  have hne : ¬ (0 < n ∧ ∀ i, x i = t) := by
    rintro ⟨-, hall⟩
    obtain ⟨i, j, hij⟩ := hx
    exact hij ((hall i).trans (hall j).symm)
  simp only [mPMRatio, mLPM_empMeasure, mUPM_empMeasure, LPMratio, hd, hne, ↓reduceIte]

theorem threshold_empMeasure {n : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    (ha0 : 0 < a) (ha1 : a < 1) :
    threshold d (empMeasure x) a ∈ Icc (varMin x) (varMax x) ∧
      LPMratio d (threshold d (empMeasure x) a) x = a := by
  obtain ⟨i, _, _⟩ := id hx
  have hn : 0 < n := Fin.pos i
  have := isProbabilityMeasure_empMeasure hn x
  have hL := memLp_empMeasure hn x (d : ℝ≥0∞)
  have hμ := empMeasure_compl_singleton_ne_zero hx
  obtain ⟨t₀, ⟨ht₀, hr⟩, -⟩ := existsUnique_LPMratio_eq hd hx (p := a) ⟨ha0.le, ha1.le⟩
  rw [← eq_threshold_of_mPMRatio_eq hd hL hμ ha0 ha1 ((mPMRatio_empMeasure hd hx t₀).trans hr)]
  exact ⟨ht₀, hr⟩

theorem LPMVaR_near_threshold_of_le_four {n : ℕ} (hd1 : 1 ≤ d) (hd4 : d ≤ 4) {x : Fin n → ℝ}
    (hx : ∃ i j, x i ≠ x j) (ha0 : 0 < a) (ha1 : a < 1) :
    LPMVaR d a x ∈ Icc (varMin x) (varMax x) ∧ |LPMratio d (LPMVaR d a x) x - a| ≤ rTol ∧
      |LPMVaR d a x - threshold d (empMeasure x) a| ≤
        uniqAt (shiftedArr x) n (locate d (shiftedArr x) n a) -
          uniqAt (shiftedArr x) n (locate d (shiftedArr x) n a - 1) := by
  have hs := LPMVaR_spec_of_le_four hd1 hd4 hx a
  rw [clampP_of_mem ⟨ha0.le, ha1.le⟩] at hs
  have ht := threshold_empMeasure (d := d) (by omega) hx ha0 ha1
  exact ⟨hs.1, hs.2, LPMVaR_dist_root_le_of_le_four hd1 hd4 hx ha0 ha1 ht.1 ht.2⟩

theorem LPMVaR_eq_threshold_of_five_le {n : ℕ} (hd : 5 ≤ d) {x : Fin n → ℝ}
    (hx : ∃ i j, x i ≠ x j) (ha0 : 0 < a) (ha1 : a < 1) :
    LPMVaR d a x = threshold d (empMeasure x) a :=
  (existsUnique_LPMratio_eq (by omega) hx ⟨ha0.le, ha1.le⟩).unique
    (LPMVaR_spec_of_five_le hd hx ⟨ha0.le, ha1.le⟩)
    (threshold_empMeasure (d := d) (by omega) hx ha0 ha1)

theorem LPMVaR_ne_threshold_snap :
    LPMVaR 4 (1 / 10 ^ 12) ![(0 : ℝ), 1] = 0 ∧
      LPMVaR 4 (1 / 10 ^ 12) ![(0 : ℝ), 1] ≠
        threshold 4 (empMeasure ![(0 : ℝ), 1]) (1 / 10 ^ 12) := by
  have hv : LPMVaR 4 (1 / 10 ^ 12) ![(0 : ℝ), 1] =
      intLPMVaR rTol 4 (1 / 10 ^ 12) ![(0 : ℝ), 1] := by
    rw [LPMVaR, ite_eq_right (by norm_num), ite_eq_left (by norm_num),
      clampP_of_mem (by norm_num)]
  have hc := intLPMVaR_snap_counterexample
  refine ⟨hv.trans hc.1, fun he => ?_⟩
  have ht := (threshold_empMeasure (d := 4) (by norm_num) (x := ![(0 : ℝ), 1]) ⟨0, 1, by simp⟩
    (a := 1 / 10 ^ 12) (by norm_num) (by norm_num)).2
  rw [← he, hv, hc.2.1] at ht
  norm_num at ht

lemma lowerDev_one_sub_bounds {s t : ℝ} (hst : t ≤ s) (x : ℝ) :
    (s - t) * lowerDev 0 t x ≤ lowerDev 1 s x - lowerDev 1 t x ∧
      lowerDev 1 s x - lowerDev 1 t x ≤ (s - t) * lowerDev 0 s x := by
  simp only [lowerDev, pow_zero, pow_one]
  split_ifs with h1 h2 h2 <;> constructor <;> linarith

theorem mLPM_one_sub_bounds [IsProbabilityMeasure μ] (hμ : ∀ t, Integrable (lowerDev 1 t) μ)
    {s t : ℝ} (hst : t ≤ s) :
    (s - t) * cdf μ t ≤ mLPM 1 s μ - mLPM 1 t μ ∧
      mLPM 1 s μ - mLPM 1 t μ ≤ (s - t) * cdf μ s := by
  have h0 : ∀ u, Integrable (lowerDev 0 u) μ := fun u =>
    integrable_lowerDev_of_le (Nat.zero_le 1) le_rfl (hμ u)
  rw [← mLPM_zero_eq_cdf, ← mLPM_zero_eq_cdf, mLPM, mLPM, mLPM, mLPM,
    ← integral_sub (hμ s) (hμ t), ← integral_const_mul, ← integral_const_mul]
  exact ⟨integral_mono ((h0 t).const_mul _) ((hμ s).sub (hμ t))
      fun x => (lowerDev_one_sub_bounds hst x).1,
    integral_mono ((hμ s).sub (hμ t)) ((h0 s).const_mul _)
      fun x => (lowerDev_one_sub_bounds hst x).2⟩

theorem hasDerivWithinAt_mLPM_one [IsProbabilityMeasure μ]
    (hμ : ∀ t, Integrable (lowerDev 1 t) μ) (t : ℝ) :
    HasDerivWithinAt (fun s => mLPM 1 s μ) (cdf μ t) (Ici t) t := by
  rw [← hasDerivWithinAt_Ioi_iff_Ici, hasDerivWithinAt_iff_tendsto_slope' self_notMem_Ioi]
  have hlim : Tendsto (cdf μ) (𝓝[>] t) (𝓝 (cdf μ t)) :=
    (((cdf μ).right_continuous t).mono Ioi_subset_Ici_self).tendsto
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hlim ?_ ?_
  · filter_upwards [self_mem_nhdsWithin] with s hs
    have hs' : t < s := hs
    rw [slope_def_field, le_div_iff₀ (sub_pos.mpr hs')]
    linarith [(mLPM_one_sub_bounds hμ hs'.le).1]
  · filter_upwards [self_mem_nhdsWithin] with s hs
    have hs' : t < s := hs
    rw [slope_def_field, div_le_iff₀ (sub_pos.mpr hs')]
    linarith [(mLPM_one_sub_bounds hμ hs'.le).2]

theorem eq_of_mLPM_one_eq [IsProbabilityMeasure μ] [IsProbabilityMeasure ν]
    (hμ : ∀ t, Integrable (lowerDev 1 t) μ) (hν : ∀ t, Integrable (lowerDev 1 t) ν)
    (h : ∀ t, mLPM 1 t μ = mLPM 1 t ν) : μ = ν := by
  refine Measure.eq_of_cdf μ ν (StieltjesFunction.ext fun t => ?_)
  have h1 := hasDerivWithinAt_mLPM_one hμ t
  rw [funext h] at h1
  exact (uniqueDiffWithinAt_Ici t).eq_deriv _ h1 (hasDerivWithinAt_mLPM_one hν t)

theorem isProbabilityMeasure_unitInterval :
    IsProbabilityMeasure (volume.restrict (Icc (0 : ℝ) 1)) := ⟨by simp⟩

lemma Iic_inter_Icc_zero_one {u : ℝ} (hu : u ≤ 1) : Iic u ∩ Icc (0 : ℝ) 1 = Icc 0 u := by
  ext x
  simp only [mem_inter_iff, mem_Iic, mem_Icc]
  exact ⟨fun h => ⟨h.2.1, h.1⟩, fun h => ⟨h.2, h.1, h.2.trans hu⟩⟩

lemma eq_uniform_of_measureReal_Iic {P : Measure ℝ} [IsProbabilityMeasure P]
    (h : ∀ u ∈ Ioo (0 : ℝ) 1, P.real (Iic u) = u) : P = volume.restrict (Icc 0 1) := by
  have := isProbabilityMeasure_unitInterval
  refine Measure.ext_of_Iic _ _ fun u => ?_
  rw [← measureReal_eq_measureReal_iff, measureReal_restrict_apply measurableSet_Iic]
  rcases le_or_gt u 0 with hu | hu
  · have hR : volume.real (Iic u ∩ Icc (0 : ℝ) 1) = 0 := by
      refine le_antisymm ((measureReal_mono (s₂ := Icc 0 u)
        (fun x hx => ⟨hx.2.1, hx.1⟩) measure_Icc_lt_top.ne).trans ?_) measureReal_nonneg
      rw [Real.volume_real_Icc, sub_zero, max_eq_right hu]
    rw [hR]
    refine le_antisymm (not_lt.mp fun hm => ?_) measureReal_nonneg
    have hε : min (P.real (Iic u)) 1 / 2 ∈ Ioo (0 : ℝ) 1 :=
      ⟨half_pos (lt_min hm one_pos), by have := min_le_right (P.real (Iic u)) 1; linarith⟩
    have hle := measureReal_mono (μ := P) (Iic_subset_Iic.mpr (hu.trans hε.1.le))
    rw [h _ hε] at hle
    have := min_le_left (P.real (Iic u)) 1
    linarith
  · rcases lt_or_ge u 1 with hu1 | hu1
    · rw [h u ⟨hu, hu1⟩, Iic_inter_Icc_zero_one hu1.le, Real.volume_real_Icc_of_le hu.le,
        sub_zero]
    · have hI : Iic u ∩ Icc (0 : ℝ) 1 = Icc 0 1 := inter_eq_right.mpr fun x hx => hx.2.trans hu1
      rw [hI, Real.volume_real_Icc_of_le zero_le_one, sub_zero]
      refine le_antisymm measureReal_le_one (not_lt.mp fun hm => ?_)
      have h0 := measureReal_nonneg (μ := P) (s := Iic u)
      have hε : (P.real (Iic u) + 1) / 2 ∈ Ioo (0 : ℝ) 1 := ⟨by linarith, by linarith⟩
      have hle := measureReal_mono (μ := P) (Iic_subset_Iic.mpr (hε.2.le.trans hu1))
      rw [h _ hε] at hle
      linarith

lemma map_eq_uniform_of_sandwich {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    [IsProbabilityMeasure P] {T : Ω → ℝ} (hT : AEMeasurable T P)
    (h : ∀ u ∈ Ioo (0 : ℝ) 1, P.real {ω | T ω < u} ≤ u ∧ u ≤ P.real {ω | T ω ≤ u}) :
    P.map T = volume.restrict (Icc 0 1) := by
  have := (Measure.isProbabilityMeasure_map_iff hT).mpr inferInstance
  refine eq_uniform_of_measureReal_Iic fun u hu => ?_
  rw [map_measureReal_apply_of_aemeasurable hT measurableSet_Iic]
  change P.real {ω | T ω ≤ u} = u
  refine le_antisymm (not_lt.mp fun hm => ?_) (h u hu).2
  have h1 := measureReal_le_one (μ := P) (s := {ω | T ω ≤ u})
  have hu' : (u + P.real {ω | T ω ≤ u}) / 2 ∈ Ioo (0 : ℝ) 1 :=
    ⟨by linarith [hu.1], by linarith [hu.2]⟩
  have h2 : P.real {ω | T ω ≤ u} ≤ P.real {ω | T ω < (u + P.real {ω | T ω ≤ u}) / 2} :=
    measureReal_mono fun ω hω => by
      have hω' : T ω ≤ u := hω
      change T ω < _
      linarith
  linarith [(h _ hu').1]

theorem map_cdf_eq_uniform [IsProbabilityMeasure μ] (hc : Continuous (cdf μ)) :
    μ.map (cdf μ) = volume.restrict (Icc 0 1) := by
  refine map_eq_uniform_of_sandwich (monotone_cdf μ).measurable.aemeasurable fun u hu => ⟨?_, ?_⟩
  · refine (measureReal_mono (s₂ := Iio (quantile μ u)) fun x hx => ?_).trans
      (measureReal_Iio_quantile_le hu.1)
    have hx' : cdf μ x < u := hx
    exact lt_of_not_ge fun hq => absurd ((quantile_le_iff μ hu.1 hu.2).mp hq) (not_le.mpr hx')
  · have hq := cdf_quantile_eq μ hc hu.1 hu.2
    calc u = μ.real (Iic (quantile μ u)) := by rw [← cdf_eq_real, hq]
      _ ≤ μ.real {x | cdf μ x ≤ u} :=
        measureReal_mono fun x hx => (monotone_cdf μ hx).trans hq.le

theorem map_cdf_comp_eq_uniform {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    [IsProbabilityMeasure P] {X : Ω → ℝ} (hX : AEMeasurable X P)
    (hc : Continuous (cdf (P.map X))) :
    P.map (fun ω => cdf (P.map X) (X ω)) = volume.restrict (Icc 0 1) := by
  have := (Measure.isProbabilityMeasure_map_iff hX).mpr inferInstance
  rw [← map_cdf_eq_uniform hc,
    AEMeasurable.map_map_of_aemeasurable (monotone_cdf _).measurable.aemeasurable hX]
  rfl

def bernoulliHalf : Measure ℝ :=
  (2 : ℝ≥0∞)⁻¹ • Measure.dirac 0 + (2 : ℝ≥0∞)⁻¹ • Measure.dirac 1

lemma bernoulliHalf_apply (s : Set ℝ) :
    bernoulliHalf s = 2⁻¹ * s.indicator 1 0 + 2⁻¹ * s.indicator 1 1 := by
  simp [bernoulliHalf]

instance : IsProbabilityMeasure bernoulliHalf :=
  ⟨by rw [bernoulliHalf_apply]; simp [ENNReal.inv_two_add_inv_two]⟩

theorem cdf_bernoulliHalf_zero : cdf bernoulliHalf 0 = 1 / 2 := by
  rw [cdf_eq_real, measureReal_def, bernoulliHalf_apply]
  simp

theorem cdf_bernoulliHalf_one : cdf bernoulliHalf 1 = 1 := by
  rw [cdf_eq_real, measureReal_def, bernoulliHalf_apply]
  simp [ENNReal.inv_two_add_inv_two]

theorem map_cdf_bernoulliHalf :
    bernoulliHalf.map (cdf bernoulliHalf) =
      (2 : ℝ≥0∞)⁻¹ • Measure.dirac (1 / 2) + (2 : ℝ≥0∞)⁻¹ • Measure.dirac 1 := by
  have hm : Measurable (cdf bernoulliHalf) := (monotone_cdf _).measurable
  have h : bernoulliHalf.map (cdf bernoulliHalf) =
      ((2 : ℝ≥0∞)⁻¹ • Measure.dirac (0 : ℝ) + (2 : ℝ≥0∞)⁻¹ • Measure.dirac (1 : ℝ)).map
        (cdf bernoulliHalf) := rfl
  rw [h, Measure.map_add _ _ hm, Measure.map_smul, Measure.map_smul, Measure.map_dirac' hm,
    Measure.map_dirac' hm, cdf_bernoulliHalf_zero, cdf_bernoulliHalf_one]
  all_goals exact hm.aemeasurable

theorem map_cdf_bernoulliHalf_ne_uniform :
    bernoulliHalf.map (cdf bernoulliHalf) ≠ volume.restrict (Icc 0 1) := by
  intro h
  have h1 : bernoulliHalf.map (cdf bernoulliHalf) {1} = 2⁻¹ := by
    rw [map_cdf_bernoulliHalf]
    simp
  have h2 : (volume.restrict (Icc (0 : ℝ) 1)) {1} = 0 :=
    measure_mono_null (fun _ hx => hx) (by simp)
  rw [h, h2] at h1
  exact absurd h1.symm (by simp)

def cdfLeft (μ : Measure ℝ) (x : ℝ) : ℝ := μ.real (Iio x)

def randTransform (μ : Measure ℝ) (z : ℝ × ℝ) : ℝ :=
  cdfLeft μ z.1 + z.2 * (cdf μ z.1 - cdfLeft μ z.1)

lemma cdfLeft_mono (μ : Measure ℝ) [IsFiniteMeasure μ] : Monotone (cdfLeft μ) :=
  fun _ _ h => measureReal_mono (Iio_subset_Iio h)

lemma cdfLeft_le_cdf [IsProbabilityMeasure μ] (x : ℝ) : cdfLeft μ x ≤ cdf μ x := by
  rw [cdf_eq_real]
  exact measureReal_mono Iio_subset_Iic_self

lemma cdf_le_cdfLeft [IsProbabilityMeasure μ] {x y : ℝ} (h : x < y) : cdf μ x ≤ cdfLeft μ y := by
  rw [cdf_eq_real]
  exact measureReal_mono (Iic_subset_Iio.mpr h)

lemma measureReal_singleton_eq [IsProbabilityMeasure μ] (x : ℝ) :
    μ.real {x} = cdf μ x - cdfLeft μ x := by
  rw [cdf_eq_real, ← Iio_union_right, measureReal_union (by simp) (measurableSet_singleton x),
    cdfLeft]
  ring

lemma measurable_randTransform (μ : Measure ℝ) [IsFiniteMeasure μ] :
    Measurable (randTransform μ) :=
  ((cdfLeft_mono μ).measurable.comp measurable_fst).add (measurable_snd.mul
    (((monotone_cdf μ).measurable.comp measurable_fst).sub
      ((cdfLeft_mono μ).measurable.comp measurable_fst)))

lemma randTransform_sandwich (μ I : Measure ℝ) [IsProbabilityMeasure μ] [IsProbabilityMeasure I]
    (hI : ∀ c ∈ Icc (0 : ℝ) 1, I.real (Icc 0 c) = c) {u : ℝ} (hu : u ∈ Ioo (0 : ℝ) 1) :
    (μ.prod I).real {z | randTransform μ z < u} ≤ u ∧
      u ≤ (μ.prod I).real {z | randTransform μ z ≤ u} := by
  have ha : cdfLeft μ (quantile μ u) ≤ u := measureReal_Iio_quantile_le hu.1
  have hb : u ≤ cdf μ (quantile μ u) := le_cdf_quantile μ hu.2
  have hm := measureReal_singleton_eq (μ := μ) (quantile μ u)
  have hc : (u - cdfLeft μ (quantile μ u)) / (cdf μ (quantile μ u) - cdfLeft μ (quantile μ u)) ∈
      Icc (0 : ℝ) 1 :=
    ⟨div_nonneg (by linarith) (by linarith), div_le_one_of_le₀ (by linarith) (by linarith)⟩
  have hmc : (cdf μ (quantile μ u) - cdfLeft μ (quantile μ u)) *
      ((u - cdfLeft μ (quantile μ u)) / (cdf μ (quantile μ u) - cdfLeft μ (quantile μ u))) =
      u - cdfLeft μ (quantile μ u) := by
    by_cases h0 : cdf μ (quantile μ u) - cdfLeft μ (quantile μ u) = 0
    · rw [h0, zero_mul]; linarith
    · rw [mul_div_cancel₀ _ h0]
  have hI1 : I.real (Icc 0 1) = 1 := hI 1 ⟨zero_le_one, le_rfl⟩
  have hIc : I.real (Icc 0 1)ᶜ = 0 := by
    rw [probReal_compl_eq_one_sub measurableSet_Icc, hI1, sub_self]
  generalize (u - cdfLeft μ (quantile μ u)) / (cdf μ (quantile μ u) - cdfLeft μ (quantile μ u))
    = c at hc hmc
  constructor
  · have hsub : {z : ℝ × ℝ | randTransform μ z < u} ⊆
        (Iio (quantile μ u) ×ˢ univ ∪ {quantile μ u} ×ˢ Icc 0 c) ∪ univ ×ˢ (Icc 0 1)ᶜ := by
      rintro ⟨x, v⟩ hz
      have hz' : cdfLeft μ x + v * (cdf μ x - cdfLeft μ x) < u := hz
      by_cases hv : v ∈ Icc (0 : ℝ) 1
      · refine Or.inl ?_
        rcases lt_trichotomy x (quantile μ u) with hx | hx | hx
        · exact Or.inl ⟨hx, mem_univ _⟩
        · refine Or.inr ⟨hx, hv.1, le_of_not_gt fun hcv => ?_⟩
          rw [hx] at hz'
          nlinarith [mul_le_mul_of_nonneg_right hcv.le
            (sub_nonneg.mpr (cdfLeft_le_cdf (μ := μ) (quantile μ u)))]
        · have h1 := cdf_le_cdfLeft (μ := μ) hx
          have h2 := mul_nonneg hv.1 (sub_nonneg.mpr (cdfLeft_le_cdf (μ := μ) x))
          linarith
      · exact Or.inr ⟨mem_univ _, hv⟩
    refine (measureReal_mono hsub).trans (((measureReal_union_le _ _).trans
      (add_le_add_left (measureReal_union_le _ _) _)).trans ?_)
    rw [measureReal_prod_prod, measureReal_prod_prod, measureReal_prod_prod, hI c hc, hIc, hm,
      hmc, probReal_univ, probReal_univ]
    change cdfLeft μ (quantile μ u) * 1 + _ + _ ≤ u
    linarith
  · have hsub : Iio (quantile μ u) ×ˢ Icc (0 : ℝ) 1 ∪ {quantile μ u} ×ˢ Icc 0 c ⊆
        {z : ℝ × ℝ | randTransform μ z ≤ u} := by
      rintro ⟨x, v⟩ (⟨hx, hv⟩ | ⟨hx, hv⟩)
      · have hx' : x < quantile μ u := hx
        have h1 := cdf_lt_of_lt_quantile μ hu.1 hx'
        have h2 := mul_le_mul_of_nonneg_right hv.2 (sub_nonneg.mpr (cdfLeft_le_cdf (μ := μ) x))
        change cdfLeft μ x + v * (cdf μ x - cdfLeft μ x) ≤ u
        linarith
      · have hx' : x = quantile μ u := hx
        have h2 := mul_le_mul_of_nonneg_right hv.2
          (sub_nonneg.mpr (cdfLeft_le_cdf (μ := μ) (quantile μ u)))
        change cdfLeft μ x + v * (cdf μ x - cdfLeft μ x) ≤ u
        rw [hx']
        linarith
    refine le_trans ?_ (measureReal_mono hsub)
    have hd : Disjoint (Iio (quantile μ u) ×ˢ Icc (0 : ℝ) 1) ({quantile μ u} ×ˢ Icc 0 c) :=
      Set.disjoint_left.mpr fun z hz hz' => (ne_of_lt hz.1) hz'.1
    rw [measureReal_union hd (by measurability),
      measureReal_prod_prod, measureReal_prod_prod, hI c hc, hI1, hm, hmc]
    change u ≤ cdfLeft μ (quantile μ u) * 1 + _
    linarith

theorem map_randTransform_eq_uniform (μ : Measure ℝ) [IsProbabilityMeasure μ] :
    (μ.prod (volume.restrict (Icc (0 : ℝ) 1))).map (randTransform μ) =
      volume.restrict (Icc 0 1) := by
  have := isProbabilityMeasure_unitInterval
  refine map_eq_uniform_of_sandwich (measurable_randTransform μ).aemeasurable fun u hu =>
    randTransform_sandwich μ _ (fun c hc => ?_) hu
  rw [measureReal_restrict_apply measurableSet_Icc,
    inter_eq_left.mpr (Icc_subset_Icc_right hc.2), Real.volume_real_Icc_of_le hc.1, sub_zero]

theorem map_randTransform_comp_eq_uniform {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    [IsProbabilityMeasure P] {X V : Ω → ℝ} (hX : AEMeasurable X P) (hV : AEMeasurable V P)
    (hXV : IndepFun X V P) (hVu : P.map V = volume.restrict (Icc 0 1)) :
    P.map (fun ω => randTransform (P.map X) (X ω, V ω)) = volume.restrict (Icc 0 1) := by
  have := (Measure.isProbabilityMeasure_map_iff hX).mpr inferInstance
  have hprod := (indepFun_iff_map_prod_eq_prod_map_map hX hV).mp hXV
  rw [hVu] at hprod
  rw [← map_randTransform_eq_uniform (P.map X), ← hprod,
    AEMeasurable.map_map_of_aemeasurable (measurable_randTransform _).aemeasurable
      (hX.prodMk hV)]
  rfl

end

end NNS.GapQuantile
