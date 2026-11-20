-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.PartialMoments
public import Mathlib.Analysis.SpecialFunctions.Log.Basic
public import Mathlib.Data.Fin.Tuple.Sort

namespace NNS

open Finset Filter Topology

@[expose] public noncomputable section

variable {n : ℕ}

def cdfCurve (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : ℝ :=
  if d = 0 then LPM 0 t x
  else if 0 < n ∧ ∀ i, x i = t then 1
  else LPM d t x / (LPM d t x + UPM d t x)

def survivalCurve (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : ℝ := 1 - cdfCurve d t x

def cdfTable (d : ℕ) (x : Fin n → ℝ) (k : Fin n) : ℝ × ℝ :=
  (x (Tuple.sort x k), cdfCurve d (x (Tuple.sort x k)) x)

def cdfTarget (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : ℝ := cdfCurve d t x

private lemma lowerDev_eq_max_pow {d : ℕ} (hd : d ≠ 0) (t v : ℝ) :
    lowerDev d t v = (max (t - v) 0) ^ d := by
  unfold lowerDev
  split_ifs with h
  · rw [max_eq_left (sub_nonneg.mpr h)]
  · rw [max_eq_right (by linarith [not_le.mp h]), zero_pow hd]

private lemma upperDev_eq_max_pow {d : ℕ} (hd : d ≠ 0) (t v : ℝ) :
    upperDev d t v = (max (v - t) 0) ^ d := by
  unfold upperDev
  split_ifs with h
  · rw [max_eq_left (sub_nonneg.mpr h.le)]
  · rw [max_eq_right (by linarith [not_lt.mp h]), zero_pow hd]

private lemma LPM_cont {d : ℕ} (hd : d ≠ 0) (x : Fin n → ℝ) :
    Continuous fun t => LPM d t x := by
  unfold LPM
  simp_rw [lowerDev_eq_max_pow hd]
  fun_prop

private lemma UPM_cont {d : ℕ} (hd : d ≠ 0) (x : Fin n → ℝ) :
    Continuous fun t => UPM d t x := by
  unfold UPM
  simp_rw [upperDev_eq_max_pow hd]
  fun_prop

private lemma denom_eq_zero_iff {d : ℕ} (hd : d ≠ 0) (t : ℝ) (x : Fin n → ℝ) :
    LPM d t x + UPM d t x = 0 ↔ ∀ i, x i = t := by
  constructor
  · intro h i
    have hn : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Fin.pos i).ne'
    have hL := (add_eq_zero_iff_of_nonneg (LPM_nonneg d t x) (UPM_nonneg d t x)).mp h
    have hl : ∑ j, lowerDev d t (x j) = 0 := by
      have := hL.1; unfold LPM at this; simpa [hn] using this
    have hu : ∑ j, upperDev d t (x j) = 0 := by
      have := hL.2; unfold UPM at this; simpa [hn] using this
    have hl' := (sum_eq_zero_iff_of_nonneg fun j _ => lowerDev_nonneg d t (x j)).mp hl i
      (mem_univ i)
    have hu' := (sum_eq_zero_iff_of_nonneg fun j _ => upperDev_nonneg d t (x j)).mp hu i
      (mem_univ i)
    rw [lowerDev_eq_max_pow hd, pow_eq_zero_iff hd] at hl'
    rw [upperDev_eq_max_pow hd, pow_eq_zero_iff hd] at hu'
    have h1 : t - x i ≤ 0 := by rw [← hl']; exact le_max_left _ _
    have h2 : x i - t ≤ 0 := by rw [← hu']; exact le_max_left _ _
    linarith
  · intro h
    unfold LPM UPM
    simp [h, lowerDev, upperDev, zero_pow hd]

private lemma denom_pos {d : ℕ} (hd : d ≠ 0) {t : ℝ} {x : Fin n → ℝ}
    (h : ¬ ∀ i, x i = t) : 0 < LPM d t x + UPM d t x :=
  lt_of_le_of_ne (add_nonneg (LPM_nonneg d t x) (UPM_nonneg d t x))
    (fun h0 => h ((denom_eq_zero_iff hd t x).mp h0.symm))

private lemma exists_pos_le_of_pos (g : Fin n → ℝ) :
    ∃ ε > 0, ∀ i, 0 < g i → ε ≤ g i := by
  by_cases hS : ({i | 0 < g i} : Finset (Fin n)).Nonempty
  · obtain ⟨i₀, hi₀, hmin⟩ := exists_min_image _ g hS
    refine ⟨g i₀, (mem_filter.mp hi₀).2, fun i hi => hmin i ?_⟩
    exact mem_filter.mpr ⟨mem_univ i, hi⟩
  · refine ⟨1, one_pos, fun i hi => absurd ⟨i, ?_⟩ hS⟩
    exact mem_filter.mpr ⟨mem_univ i, hi⟩

theorem cdfCurve_eq_ratio {d : ℕ} (hd : d ≠ 0) {t : ℝ} {x : Fin n → ℝ}
    (h : ∃ i, x i ≠ t) :
    cdfCurve d t x = LPM d t x / (LPM d t x + UPM d t x) := by
  obtain ⟨i, hi⟩ := h
  have : ¬ (0 < n ∧ ∀ i, x i = t) := fun h' => hi (h'.2 i)
  simp [cdfCurve, hd, this]

theorem cdfCurve_eq_one_of_forall_le (hn : 0 < n) (d : ℕ) {t : ℝ} {x : Fin n → ℝ}
    (h : ∀ i, x i ≤ t) : cdfCurve d t x = 1 := by
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  rcases eq_or_ne d 0 with rfl | hd
  · simp only [cdfCurve, ↓reduceIte, LPM_zero_eq_ecdf]
    rw [filter_true_of_mem fun i _ => h i, card_univ, Fintype.card_fin, div_self hn']
  · by_cases hc : ∀ i, x i = t
    · simp [cdfCurve, hd, hn, hc]
    · have hc' : ∃ i, x i ≠ t := by simpa using hc
      rw [cdfCurve_eq_ratio hd hc']
      have hU : UPM d t x = 0 := by
        unfold UPM upperDev
        simp [fun i => not_lt.mpr (h i)]
      have hpos := denom_pos hd (t := t) (x := x) hc
      rw [hU, add_zero] at hpos ⊢
      exact div_self hpos.ne'

theorem cdfCurve_eq_zero_of_forall_lt (d : ℕ) {t : ℝ} {x : Fin n → ℝ}
    (h : ∀ i, t < x i) : cdfCurve d t x = 0 := by
  have hL : ∀ d, LPM d t x = 0 := by
    intro d
    unfold LPM lowerDev
    simp [fun i => not_le.mpr (h i)]
  rcases eq_or_ne d 0 with rfl | hd
  · simp [cdfCurve, hL]
  · have : ¬ (0 < n ∧ ∀ i, x i = t) := fun h' =>
      (h ⟨0, h'.1⟩).ne' (h'.2 ⟨0, h'.1⟩)
    simp [cdfCurve, hd, this, hL]

theorem cdfCurve_eq_zero_of_le_min {d : ℕ} (hd : d ≠ 0) {t : ℝ} {x : Fin n → ℝ}
    (hmin : ∀ i, t ≤ x i) (h : ∃ i, x i ≠ t) : cdfCurve d t x = 0 := by
  rw [cdfCurve_eq_ratio hd h]
  have hL : LPM d t x = 0 := by
    unfold LPM
    rw [sum_eq_zero fun i _ => ?_, zero_div]
    rw [lowerDev_eq_max_pow hd, max_eq_right (by linarith [hmin i]), zero_pow hd]
  rw [hL, zero_div]

theorem cdfCurve_of_const (hn : 0 < n) (d : ℕ) {c : ℝ} {x : Fin n → ℝ}
    (hx : ∀ i, x i = c) (t : ℝ) : cdfCurve d t x = if c ≤ t then 1 else 0 := by
  split_ifs with h
  · exact cdfCurve_eq_one_of_forall_le hn d fun i => (hx i).le.trans h
  · exact cdfCurve_eq_zero_of_forall_lt d fun i => (hx i).symm ▸ not_le.mp h

theorem not_continuousAt_cdfCurve_of_const (hn : 0 < n) (d : ℕ) {c : ℝ}
    {x : Fin n → ℝ} (hx : ∀ i, x i = c) : ¬ ContinuousAt (fun t => cdfCurve d t x) c := by
  intro hc
  obtain ⟨δ, hδ, hball⟩ := Metric.continuousAt_iff.mp hc (1 / 2) (by norm_num)
  have hdist : dist (c - δ / 2) c < δ := by
    rw [Real.dist_eq, show c - δ / 2 - c = -(δ / 2) by ring, abs_neg,
      abs_of_pos (by linarith)]
    linarith
  have := hball hdist
  rw [cdfCurve_of_const hn d hx, cdfCurve_of_const hn d hx, Real.dist_eq] at this
  simp only [show ¬ c ≤ c - δ / 2 by linarith, le_refl, ↓reduceIte] at this
  norm_num at this

theorem cdfCurve_continuous {d : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ}
    (hx : ∃ i j, x i ≠ x j) : Continuous fun t => cdfCurve d t x := by
  obtain ⟨i, j, hij⟩ := hx
  have hne : ∀ t, ∃ k, x k ≠ t := fun t => by
    by_cases h : x i = t
    · exact ⟨j, fun h' => hij (h.trans h'.symm)⟩
    · exact ⟨i, h⟩
  have hfun : (fun t => cdfCurve d t x) =
      fun t => LPM d t x / (LPM d t x + UPM d t x) :=
    funext fun t => cdfCurve_eq_ratio hd (hne t)
  rw [hfun]
  exact (LPM_cont hd x).div ((LPM_cont hd x).add (UPM_cont hd x))
    fun t => (denom_pos hd (by simpa using hne t)).ne'

theorem cdfCurve_zero_eq_ecdf (t : ℝ) (x : Fin n → ℝ) :
    cdfCurve 0 t x = (#{i | x i ≤ t} : ℝ) / n := by
  simp only [cdfCurve, ↓reduceIte, LPM_zero_eq_ecdf]

theorem cdfCurve_zero_eq_lt_add_eq (t : ℝ) (x : Fin n → ℝ) :
    cdfCurve 0 t x = (#{i | x i < t} : ℝ) / n + (#{i | x i = t} : ℝ) / n := by
  rw [cdfCurve_zero_eq_ecdf, ← add_div]
  congr 1
  rw [← Nat.cast_add, ← card_union_of_disjoint]
  · congr 2
    ext i
    simp [le_iff_lt_or_eq]
  · rw [disjoint_filter]
    intro i _ h1 h2
    exact h1.ne h2

theorem cdfCurve_zero_eq_of_no_obs (x : Fin n → ℝ) {t s : ℝ}
    (h : ∀ i, x i ≤ s → x i ≤ t) (hts : t ≤ s) : cdfCurve 0 s x = cdfCurve 0 t x := by
  rw [cdfCurve_zero_eq_ecdf, cdfCurve_zero_eq_ecdf]
  congr 3
  ext i
  simp only [mem_filter, mem_univ, true_and]
  exact ⟨h i, fun hi => hi.trans hts⟩

theorem cdfCurve_zero_eventually_const_left (t : ℝ) (x : Fin n → ℝ) :
    ∃ ε > 0, ∀ s, t - ε < s → s < t → cdfCurve 0 s x = (#{i | x i < t} : ℝ) / n := by
  obtain ⟨ε, hε, hgap⟩ := exists_pos_le_of_pos fun i => t - x i
  refine ⟨ε, hε, fun s hs hst => ?_⟩
  rw [cdfCurve_zero_eq_ecdf]
  congr 3
  ext i
  simp only [mem_filter, mem_univ, true_and]
  constructor
  · intro hi; linarith
  · intro hi
    have := hgap i (by linarith)
    linarith

theorem cdfCurve_zero_tendsto_left (t : ℝ) (x : Fin n → ℝ) :
    Tendsto (fun s => cdfCurve 0 s x) (𝓝[<] t) (𝓝 ((#{i | x i < t} : ℝ) / n)) := by
  obtain ⟨ε, hε, h⟩ := cdfCurve_zero_eventually_const_left t x
  have hev : (fun s => cdfCurve 0 s x) =ᶠ[𝓝[<] t] fun _ => (#{i | x i < t} : ℝ) / n :=
    Filter.mem_of_superset (Ioo_mem_nhdsLT (by linarith : t - ε < t))
      fun s hs => h s hs.1 hs.2
  exact tendsto_const_nhds.congr' hev.symm

end

end NNS
