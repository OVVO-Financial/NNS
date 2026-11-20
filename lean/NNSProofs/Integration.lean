-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
public import NNSProofs.PartialMoments

namespace NNS

open Finset Set

@[expose] public noncomputable section

variable {n m : ℕ}

lemma integral_congr_of_Ioo_aux {f g : ℝ → ℝ} {a b : ℝ} (hab : a ≤ b)
    (h : ∀ s ∈ Ioo a b, f s = g s) : ∫ s in a..b, f s = ∫ s in a..b, g s := by
  rw [intervalIntegral.integral_of_le hab, intervalIntegral.integral_of_le hab,
    MeasureTheory.integral_Ioc_eq_integral_Ioo, MeasureTheory.integral_Ioc_eq_integral_Ioo]
  exact MeasureTheory.setIntegral_congr_fun measurableSet_Ioo h

lemma lowerDev_monotone (d : ℕ) (v : ℝ) : Monotone fun s => lowerDev d s v := by
  intro s t hst
  dsimp only
  unfold lowerDev
  split_ifs with h1 h2 h2
  · exact pow_le_pow_left₀ (sub_nonneg.mpr h1) (by linarith) d
  · exact absurd (le_trans h1 hst) h2
  · exact pow_nonneg (sub_nonneg.mpr h2) d
  · exact le_refl 0

lemma integral_shifted_pow_aux (d : ℕ) (u t : ℝ) :
    ∫ s in u..t, (s - u) ^ d = (t - u) ^ (d + 1) / (d + 1) := by
  have h := intervalIntegral.integral_comp_sub_right (a := u) (b := t) (fun x : ℝ => x ^ d) u
  rw [h, sub_self, integral_pow, zero_pow (Nat.succ_ne_zero d), sub_zero]

theorem lowerDev_succ_eq_integral (d : ℕ) {a v : ℝ} (hav : a ≤ v) (t : ℝ) :
    lowerDev (d + 1) t v = (d + 1) * ∫ s in a..t, lowerDev d s v := by
  have hd : ((d : ℝ) + 1) ≠ 0 := Nat.cast_add_one_ne_zero d
  by_cases hvt : v ≤ t
  · have hI : IntervalIntegrable (fun s => lowerDev d s v) MeasureTheory.volume a v :=
      (lowerDev_monotone d v).intervalIntegrable
    have hJ : IntervalIntegrable (fun s => lowerDev d s v) MeasureTheory.volume v t :=
      (lowerDev_monotone d v).intervalIntegrable
    rw [← intervalIntegral.integral_add_adjacent_intervals hI hJ]
    have h1 : ∫ s in a..v, lowerDev d s v = ∫ s in a..v, (0 : ℝ) :=
      integral_congr_of_Ioo_aux hav fun s hs => by
        simp [lowerDev, not_le.mpr hs.2]
    have h2 : ∫ s in v..t, lowerDev d s v = ∫ s in v..t, (s - v) ^ d := by
      refine intervalIntegral.integral_congr fun s hs => ?_
      rw [uIcc_of_le hvt] at hs
      simp [lowerDev, hs.1]
    rw [h1, h2, intervalIntegral.integral_zero, zero_add, integral_shifted_pow_aux,
      mul_div_cancel₀ _ hd]
    simp [lowerDev, hvt]
  · have hvt' : t < v := not_le.mp hvt
    have h1 : ∫ s in a..t, lowerDev d s v = 0 := by
      rcases le_total a t with h | h
      · rw [integral_congr_of_Ioo_aux h (g := fun _ => 0) fun s hs => by
          simp [lowerDev, not_le.mpr (hs.2.trans hvt')], intervalIntegral.integral_zero]
      · rw [intervalIntegral.integral_symm, integral_congr_of_Ioo_aux h (g := fun _ => 0)
          fun s hs => by simp [lowerDev, not_le.mpr (hs.2.trans_le hav)],
          intervalIntegral.integral_zero, neg_zero]
    rw [h1, mul_zero]
    simp [lowerDev, hvt]

theorem LPM_monotone_target (d : ℕ) (x : Fin n → ℝ) : Monotone fun t => LPM d t x := by
  intro s t hst
  dsimp only
  unfold LPM
  gcongr with i
  exact lowerDev_monotone d (x i) hst

lemma LPM_intervalIntegrable (d : ℕ) (x : Fin n → ℝ) (a b : ℝ) :
    IntervalIntegrable (fun s => LPM d s x) MeasureTheory.volume a b :=
  (LPM_monotone_target d x).intervalIntegrable

theorem LPM_eq_zero_of_lt {d : ℕ} {t : ℝ} {x : Fin n → ℝ} (h : ∀ i, t < x i) :
    LPM d t x = 0 := by
  unfold LPM
  rw [sum_eq_zero fun i _ => by simp [lowerDev, not_le.mpr (h i)], zero_div]

theorem UPM_eq_zero_of_le {d : ℕ} {t : ℝ} {x : Fin n → ℝ} (h : ∀ i, x i ≤ t) :
    UPM d t x = 0 := by
  unfold UPM
  rw [sum_eq_zero fun i _ => by simp [upperDev, not_lt.mpr (h i)], zero_div]

lemma exists_sample_lowerBound (x : Fin n → ℝ) : ∃ a, ∀ i, a ≤ x i := by
  obtain ⟨a, ha⟩ := (Set.finite_range x).bddBelow
  exact ⟨a, fun i => ha ⟨i, rfl⟩⟩

lemma exists_sample_upperBound (x : Fin n → ℝ) : ∃ b, ∀ i, x i ≤ b := by
  obtain ⟨b, hb⟩ := (Set.finite_range x).bddAbove
  exact ⟨b, fun i => hb ⟨i, rfl⟩⟩

theorem LPM_succ_eq_integral (d : ℕ) {a : ℝ} {x : Fin n → ℝ} (ha : ∀ i, a ≤ x i) (t : ℝ) :
    LPM (d + 1) t x = (d + 1) * ∫ s in a..t, LPM d s x := by
  unfold LPM
  rw [intervalIntegral.integral_div, intervalIntegral.integral_finsetSum
    fun i _ => (lowerDev_monotone d (x i)).intervalIntegrable, mul_div_assoc', mul_sum]
  congr 1
  exact sum_congr rfl fun i _ => lowerDev_succ_eq_integral d (ha i) t

theorem LPM_succ_sub_eq_integral (d : ℕ) (x : Fin n → ℝ) (u t : ℝ) :
    LPM (d + 1) t x - LPM (d + 1) u x = (d + 1) * ∫ s in u..t, LPM d s x := by
  obtain ⟨a, ha⟩ := exists_sample_lowerBound x
  rw [LPM_succ_eq_integral d ha t, LPM_succ_eq_integral d ha u, ← mul_sub,
    intervalIntegral.integral_interval_sub_left (LPM_intervalIntegrable d x a t)
      (LPM_intervalIntegrable d x a u)]

theorem LPM_one_of_forall_le (hn : 0 < n) {t : ℝ} {x : Fin n → ℝ} (h : ∀ i, x i ≤ t) :
    LPM 1 t x = t - mean x := by
  have h1 := UPM_one_sub_LPM_one hn t x
  rw [UPM_eq_zero_of_le h] at h1
  linarith

theorem LPM_one_le_iff_mean_le (hn : 0 < n) (hm : 0 < m) {t : ℝ} {x : Fin n → ℝ}
    {y : Fin m → ℝ} (hx : ∀ i, x i ≤ t) (hy : ∀ j, y j ≤ t) :
    LPM 1 t x ≤ LPM 1 t y ↔ mean y ≤ mean x := by
  rw [LPM_one_of_forall_le hn hx, LPM_one_of_forall_le hm hy]
  constructor <;> intro h <;> linarith

theorem LPM_pos_of_lt {d : ℕ} {t : ℝ} {x : Fin n → ℝ} {i : Fin n} (hi : x i < t) :
    0 < LPM d t x := by
  have hn : (0 : ℝ) < n := Nat.cast_pos.mpr (Fin.pos i)
  unfold LPM
  refine div_pos (sum_pos' (fun j _ => lowerDev_nonneg d t (x j)) ⟨i, mem_univ i, ?_⟩) hn
  simp only [lowerDev, hi.le, ↓reduceIte]
  exact pow_pos (sub_pos.mpr hi) d

theorem exists_le_of_forall_LPM_le {d : ℕ} {x : Fin n → ℝ} {y : Fin m → ℝ}
    (h : ∀ s, LPM d s x ≤ LPM d s y) (i : Fin n) : ∃ j, y j ≤ x i := by
  by_contra hcon
  push Not at hcon
  have h1 : ∀ᶠ t in nhdsWithin (x i) (Ioi (x i)), t ∈ Ioi (x i) := eventually_mem_nhdsWithin
  have h2 : ∀ᶠ t in nhdsWithin (x i) (Ioi (x i)), ∀ j, t < y j :=
    Filter.eventually_all.mpr fun j => nhdsWithin_le_nhds (Iio_mem_nhds (hcon j))
  obtain ⟨t, hxt, hty⟩ := (h1.and h2).exists
  have hpos := LPM_pos_of_lt (d := d) (show x i < t from hxt)
  rw [← LPM_eq_zero_of_lt (d := d) hty] at hpos
  exact absurd (h t) (not_le.mpr hpos)

theorem mean_le_of_LPM_one_le (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ} {y : Fin m → ℝ}
    (h : ∀ s, LPM 1 s x ≤ LPM 1 s y) : mean y ≤ mean x := by
  obtain ⟨bx, hbx⟩ := exists_sample_upperBound x
  obtain ⟨by', hby⟩ := exists_sample_upperBound y
  exact (LPM_one_le_iff_mean_le hn hm (t := max bx by')
    (fun i => (hbx i).trans (le_max_left _ _))
    (fun j => (hby j).trans (le_max_right _ _))).mp (h _)

lemma LPM_of_fin_zero (d : ℕ) (t : ℝ) (x : Fin 0 → ℝ) : LPM d t x = 0 := by
  simp [LPM]

theorem LPM_le_on_Iic_of_le_on_Icc {d : ℕ} {a b : ℝ} {x : Fin n → ℝ} {y : Fin m → ℝ}
    (ha : ∀ i, a ≤ x i) (h : ∀ s ∈ Icc a b, LPM d s x ≤ LPM d s y) :
    ∀ s ≤ b, LPM d s x ≤ LPM d s y := by
  intro s hsb
  by_cases has : a ≤ s
  · exact h s ⟨has, hsb⟩
  · rw [LPM_eq_zero_of_lt fun i => lt_of_lt_of_le (not_le.mp has) (ha i)]
    exact LPM_nonneg _ _ _

theorem LPM_one_le_of_le_on_Iic (hm : 0 < m) {b : ℝ} {x : Fin n → ℝ} {y : Fin m → ℝ}
    (hx : ∀ i, x i ≤ b) (hy : ∀ j, y j ≤ b) (h : ∀ s ≤ b, LPM 1 s x ≤ LPM 1 s y) (t : ℝ) :
    LPM 1 t x ≤ LPM 1 t y := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · rw [LPM_of_fin_zero]; exact LPM_nonneg _ _ _
  by_cases htb : t ≤ b
  · exact h t htb
  · have htb' : b ≤ t := (not_le.mp htb).le
    have hb := (LPM_one_le_iff_mean_le hn hm hx hy).mp (h b le_rfl)
    exact (LPM_one_le_iff_mean_le hn hm (fun i => (hx i).trans htb')
      (fun j => (hy j).trans htb')).mpr hb

theorem LPM_two_le_of_le_on_Iic (hm : 0 < m) {b : ℝ} {x : Fin n → ℝ} {y : Fin m → ℝ}
    (hx : ∀ i, x i ≤ b) (hy : ∀ j, y j ≤ b) (hmean : mean y ≤ mean x)
    (h : ∀ s ≤ b, LPM 2 s x ≤ LPM 2 s y) (t : ℝ) : LPM 2 t x ≤ LPM 2 t y := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · rw [LPM_of_fin_zero]; exact LPM_nonneg _ _ _
  by_cases htb : t ≤ b
  · exact h t htb
  · have htb' : b ≤ t := (not_le.mp htb).le
    have ex := LPM_succ_sub_eq_integral 1 x b t
    have ey := LPM_succ_sub_eq_integral 1 y b t
    have hmono : ∫ s in b..t, LPM 1 s x ≤ ∫ s in b..t, LPM 1 s y :=
      intervalIntegral.integral_mono_on htb' (LPM_intervalIntegrable 1 x b t)
        (LPM_intervalIntegrable 1 y b t) fun s hs => by
          rw [LPM_one_of_forall_le hn fun i => (hx i).trans hs.1,
            LPM_one_of_forall_le hm fun j => (hy j).trans hs.1]
          linarith
    have hb := h b le_rfl
    norm_num at ex ey
    linarith

theorem LPM_two_le_of_le_on_Icc (hm : 0 < m) {a b : ℝ} {x : Fin n → ℝ} {y : Fin m → ℝ}
    (ha : ∀ i, a ≤ x i) (hx : ∀ i, x i ≤ b) (hy : ∀ j, y j ≤ b) (hmean : mean y ≤ mean x)
    (h : ∀ s ∈ Icc a b, LPM 2 s x ≤ LPM 2 s y) (t : ℝ) : LPM 2 t x ≤ LPM 2 t y :=
  LPM_two_le_of_le_on_Iic hm hx hy hmean (LPM_le_on_Iic_of_le_on_Icc ha h) t

def lowerCount (t : ℝ) (x : Fin n → ℝ) : ℕ := #(univ.filter fun i => x i ≤ t)

def lowerSum1 (t : ℝ) (x : Fin n → ℝ) : ℝ := ∑ i ∈ univ.filter (fun i => x i ≤ t), x i

def lowerSum2 (t : ℝ) (x : Fin n → ℝ) : ℝ := ∑ i ∈ univ.filter (fun i => x i ≤ t), x i ^ 2

def tsdCoefA (t : ℝ) (x : Fin n → ℝ) (y : Fin m → ℝ) : ℝ :=
  (lowerCount t x : ℝ) / n - (lowerCount t y : ℝ) / m

def tsdCoefB (t : ℝ) (x : Fin n → ℝ) (y : Fin m → ℝ) : ℝ :=
  lowerSum1 t x / n - lowerSum1 t y / m

def IsKnot (x : Fin n → ℝ) (y : Fin m → ℝ) (t : ℝ) : Prop := (∃ i, x i = t) ∨ ∃ j, y j = t

lemma lowerDev_succ_eq_on_Icc (d : ℕ) {u w s v : ℝ} (hv : v ∉ Ioo u w) (hs : s ∈ Icc u w) :
    lowerDev (d + 1) s v = if v ≤ u then (s - v) ^ (d + 1) else 0 := by
  unfold lowerDev
  by_cases hvu : v ≤ u
  · simp [hvu, hvu.trans hs.1]
  · have hwv : w ≤ v := by
      by_contra h
      exact hv ⟨not_le.mp hvu, not_le.mp h⟩
    simp only [hvu, ↓reduceIte]
    split_ifs with hvs
    · have : s = v := le_antisymm (hs.2.trans hwv) hvs
      simp [this]
    · rfl

theorem LPM_succ_eq_sum_filter (d : ℕ) {u w s : ℝ} {x : Fin n → ℝ}
    (hx : ∀ i, x i ∉ Ioo u w) (hs : s ∈ Icc u w) :
    LPM (d + 1) s x = (∑ i ∈ univ.filter (fun i => x i ≤ u), (s - x i) ^ (d + 1)) / n := by
  unfold LPM
  rw [sum_filter]
  congr 1
  exact sum_congr rfl fun i _ => lowerDev_succ_eq_on_Icc d (hx i) hs

theorem LPM_one_eq_on_Icc {u w s : ℝ} {x : Fin n → ℝ} (hx : ∀ i, x i ∉ Ioo u w)
    (hs : s ∈ Icc u w) : LPM 1 s x = ((lowerCount u x : ℝ) * s - lowerSum1 u x) / n := by
  rw [LPM_succ_eq_sum_filter 0 hx hs, lowerCount, lowerSum1]
  simp [sum_sub_distrib]

theorem LPM_two_eq_on_Icc {u w s : ℝ} {x : Fin n → ℝ} (hx : ∀ i, x i ∉ Ioo u w)
    (hs : s ∈ Icc u w) :
    LPM 2 s x = ((lowerCount u x : ℝ) * s ^ 2 - 2 * s * lowerSum1 u x + lowerSum2 u x) / n := by
  rw [LPM_succ_eq_sum_filter 1 hx hs, lowerCount, lowerSum1, lowerSum2]
  congr 1
  rw [sum_congr rfl fun i _ => (show (s - x i) ^ (1 + 1) = s ^ 2 - 2 * s * x i + x i ^ 2 by
    ring), sum_add_distrib, sum_sub_distrib, sum_const, nsmul_eq_mul, ← mul_sum]

lemma affine_nonpos_on_Icc {α β u w s : ℝ} (hs : s ∈ Icc u w) (hu : α * u + β ≤ 0)
    (hw : α * w + β ≤ 0) : α * s + β ≤ 0 := by
  rcases le_total 0 α with h | h
  · nlinarith [mul_le_mul_of_nonneg_left hs.2 h]
  · nlinarith [mul_le_mul_of_nonpos_left hs.1 h]

lemma quad_nonpos_on_Icc {a b c u w : ℝ} (hu : a * u ^ 2 - 2 * b * u + c ≤ 0)
    (hw : a * w ^ 2 - 2 * b * w + c ≤ 0)
    (hv : a < 0 → b / a ∈ Ioo u w → a * (b / a) ^ 2 - 2 * b * (b / a) + c ≤ 0) :
    ∀ s ∈ Icc u w, a * s ^ 2 - 2 * b * s + c ≤ 0 := by
  intro s hs
  have e1 : a * s ^ 2 - 2 * b * s + c - (a * u ^ 2 - 2 * b * u + c)
      = (s - u) * (a * (s + u) - 2 * b) := by ring
  have e2 : a * s ^ 2 - 2 * b * s + c - (a * w ^ 2 - 2 * b * w + c)
      = (s - w) * (a * (s + w) - 2 * b) := by ring
  have hsu : 0 ≤ s - u := sub_nonneg.mpr hs.1
  have hsw : 0 ≤ w - s := sub_nonneg.mpr hs.2
  rcases le_or_gt 0 a with ha | ha
  · rcases le_total (a * (s + u) - 2 * b) 0 with h | h
    · nlinarith [mul_nonneg hsu (neg_nonneg.mpr h)]
    · have h' : 0 ≤ a * (s + w) - 2 * b := by nlinarith [mul_nonneg ha (add_nonneg hsu hsw)]
      nlinarith [mul_nonneg hsw h']
  · have hb : b = a * (b / a) := (mul_div_cancel₀ b ha.ne).symm
    set v := b / a with hvdef
    by_cases hin : v ∈ Ioo u w
    · have hv' := hv ha hin
      have e3 : a * s ^ 2 - 2 * b * s + c - (a * v ^ 2 - 2 * b * v + c) = a * (s - v) ^ 2 := by
        rw [hb]; ring
      nlinarith [mul_nonneg (neg_nonneg.mpr ha.le) (sq_nonneg (s - v))]
    · rw [Set.mem_Ioo, not_and_or, not_lt, not_lt] at hin
      rcases hin with h | h
      · have h' : a * (s + u) - 2 * b ≤ 0 := by
          rw [hb]
          nlinarith [mul_nonneg (neg_nonneg.mpr ha.le) (by linarith : 0 ≤ s + u - 2 * v)]
        nlinarith [mul_nonneg hsu (neg_nonneg.mpr h')]
      · have h' : 0 ≤ a * (s + w) - 2 * b := by
          rw [hb]
          nlinarith [mul_nonneg (neg_nonneg.mpr ha.le) (by linarith : 0 ≤ 2 * v - s - w)]
        nlinarith [mul_nonneg hsw h']

theorem LPM_one_le_on_Icc_of_knots {u w : ℝ} {x : Fin n → ℝ} {y : Fin m → ℝ}
    (hx : ∀ i, x i ∉ Ioo u w) (hy : ∀ j, y j ∉ Ioo u w)
    (hu : LPM 1 u x ≤ LPM 1 u y) (hw : LPM 1 w x ≤ LPM 1 w y) :
    ∀ s ∈ Icc u w, LPM 1 s x ≤ LPM 1 s y := by
  intro s hs
  have huw : u ≤ w := hs.1.trans hs.2
  have key : ∀ r ∈ Icc u w,
      LPM 1 r x - LPM 1 r y = tsdCoefA u x y * r + -tsdCoefB u x y := fun r hr => by
    rw [LPM_one_eq_on_Icc hx hr, LPM_one_eq_on_Icc hy hr, tsdCoefA, tsdCoefB]
    ring
  have h := affine_nonpos_on_Icc hs
    (by rw [← key u ⟨le_rfl, huw⟩]; linarith) (by rw [← key w ⟨huw, le_rfl⟩]; linarith)
  rw [← key s hs] at h
  linarith

theorem LPM_two_le_on_Icc_of_knots {u w : ℝ} {x : Fin n → ℝ} {y : Fin m → ℝ}
    (hx : ∀ i, x i ∉ Ioo u w) (hy : ∀ j, y j ∉ Ioo u w)
    (hu : LPM 2 u x ≤ LPM 2 u y) (hw : LPM 2 w x ≤ LPM 2 w y)
    (hv : tsdCoefA u x y < 0 → tsdCoefB u x y / tsdCoefA u x y ∈ Ioo u w →
      LPM 2 (tsdCoefB u x y / tsdCoefA u x y) x ≤ LPM 2 (tsdCoefB u x y / tsdCoefA u x y) y) :
    ∀ s ∈ Icc u w, LPM 2 s x ≤ LPM 2 s y := by
  intro s hs
  have huw : u ≤ w := hs.1.trans hs.2
  set c := lowerSum2 u x / n - lowerSum2 u y / m
  have key : ∀ r ∈ Icc u w, LPM 2 r x - LPM 2 r y
      = tsdCoefA u x y * r ^ 2 - 2 * tsdCoefB u x y * r + c := fun r hr => by
    rw [LPM_two_eq_on_Icc hx hr, LPM_two_eq_on_Icc hy hr, tsdCoefA, tsdCoefB]
    ring
  have h := quad_nonpos_on_Icc (c := c)
    (by rw [← key u ⟨le_rfl, huw⟩]; linarith) (by rw [← key w ⟨huw, le_rfl⟩]; linarith)
    (fun ha hin => by
      rw [← key _ (Ioo_subset_Icc_self hin)]
      linarith [hv ha hin]) s hs
  rw [← key s hs] at h
  linarith

lemma exists_bracket (K : Finset ℝ) {s : ℝ} (h1 : ∃ z ∈ K, z ≤ s) (h2 : ∃ z ∈ K, s ≤ z) :
    ∃ u ∈ K, ∃ w ∈ K, s ∈ Icc u w ∧ ∀ z ∈ K, z ∉ Ioo u w := by
  have hn1 : (K.filter (· ≤ s)).Nonempty := by
    obtain ⟨z, hz, hzs⟩ := h1
    exact ⟨z, mem_filter.mpr ⟨hz, hzs⟩⟩
  have hn2 : (K.filter (s ≤ ·)).Nonempty := by
    obtain ⟨z, hz, hzs⟩ := h2
    exact ⟨z, mem_filter.mpr ⟨hz, hzs⟩⟩
  have hu := mem_filter.mp ((K.filter (· ≤ s)).max'_mem hn1)
  have hw := mem_filter.mp ((K.filter (s ≤ ·)).min'_mem hn2)
  refine ⟨_, hu.1, _, hw.1, ⟨hu.2, hw.2⟩, fun z hz hzi => ?_⟩
  rcases le_total z s with h | h
  · exact absurd (le_max' _ z (mem_filter.mpr ⟨hz, h⟩)) (not_le.mpr hzi.1)
  · exact absurd (min'_le _ z (mem_filter.mpr ⟨hz, h⟩)) (not_le.mpr hzi.2)

lemma isKnot_iff_mem {x : Fin n → ℝ} {y : Fin m → ℝ} {t : ℝ} :
    IsKnot x y t ↔ t ∈ Finset.univ.image x ∪ Finset.univ.image y := by
  simp [IsKnot, Finset.mem_union, Finset.mem_image]

lemma LPM_le_on_Iic_of_knots {d : ℕ} {x : Fin n → ℝ} {y : Fin m → ℝ} {b : ℝ}
    (hb : IsKnot x y b)
    (hloc : ∀ u w, IsKnot x y u → IsKnot x y w → (∀ i, x i ∉ Ioo u w) →
      (∀ j, y j ∉ Ioo u w) → ∀ s ∈ Icc u w, LPM d s x ≤ LPM d s y) :
    ∀ s ≤ b, LPM d s x ≤ LPM d s y := by
  intro s hsb
  by_cases hlow : ∀ i, s < x i
  · rw [LPM_eq_zero_of_lt hlow]
    exact LPM_nonneg _ _ _
  push Not at hlow
  obtain ⟨i, hi⟩ := hlow
  obtain ⟨u, hu, w, hw, hs, hK⟩ := exists_bracket (Finset.univ.image x ∪ Finset.univ.image y)
    ⟨x i, isKnot_iff_mem.mp (Or.inl ⟨i, rfl⟩), hi⟩ ⟨b, isKnot_iff_mem.mp hb, hsb⟩
  exact hloc u w (isKnot_iff_mem.mpr hu) (isKnot_iff_mem.mpr hw)
    (fun i => hK _ (isKnot_iff_mem.mp (Or.inl ⟨i, rfl⟩)))
    (fun j => hK _ (isKnot_iff_mem.mp (Or.inr ⟨j, rfl⟩))) s hs

theorem LPM_one_le_of_le_at_knots (hm : 0 < m) {x : Fin n → ℝ} {y : Fin m → ℝ}
    (h : ∀ t, IsKnot x y t → LPM 1 t x ≤ LPM 1 t y) (t : ℝ) : LPM 1 t x ≤ LPM 1 t y := by
  set K := Finset.univ.image x ∪ Finset.univ.image y
  have hK : K.Nonempty := ⟨y ⟨0, hm⟩, mem_union_right _ (Finset.mem_image_of_mem _ (mem_univ _))⟩
  have hb : IsKnot x y (K.max' hK) := isKnot_iff_mem.mpr (K.max'_mem hK)
  refine LPM_one_le_of_le_on_Iic hm
    (fun i => le_max' K _ (mem_union_left _ (Finset.mem_image_of_mem _ (mem_univ i))))
    (fun j => le_max' K _ (mem_union_right _ (Finset.mem_image_of_mem _ (mem_univ j))))
    (LPM_le_on_Iic_of_knots hb fun u w hu hw hx hy =>
      LPM_one_le_on_Icc_of_knots hx hy (h u hu) (h w hw)) t

theorem LPM_two_le_of_le_at_knots (hm : 0 < m) {x : Fin n → ℝ} {y : Fin m → ℝ}
    (hmean : mean y ≤ mean x) (h : ∀ t, IsKnot x y t → LPM 2 t x ≤ LPM 2 t y)
    (hv : ∀ u w, IsKnot x y u → IsKnot x y w → (∀ i, x i ∉ Ioo u w) →
      (∀ j, y j ∉ Ioo u w) → tsdCoefA u x y < 0 →
      tsdCoefB u x y / tsdCoefA u x y ∈ Ioo u w →
      LPM 2 (tsdCoefB u x y / tsdCoefA u x y) x ≤ LPM 2 (tsdCoefB u x y / tsdCoefA u x y) y)
    (t : ℝ) : LPM 2 t x ≤ LPM 2 t y := by
  set K := Finset.univ.image x ∪ Finset.univ.image y
  have hK : K.Nonempty := ⟨y ⟨0, hm⟩, mem_union_right _ (Finset.mem_image_of_mem _ (mem_univ _))⟩
  have hb : IsKnot x y (K.max' hK) := isKnot_iff_mem.mpr (K.max'_mem hK)
  refine LPM_two_le_of_le_on_Iic hm
    (fun i => le_max' K _ (mem_union_left _ (Finset.mem_image_of_mem _ (mem_univ i))))
    (fun j => le_max' K _ (mem_union_right _ (Finset.mem_image_of_mem _ (mem_univ j)))) hmean
    (LPM_le_on_Iic_of_knots hb fun u w hu hw hx hy =>
      LPM_two_le_on_Icc_of_knots hx hy (h u hu) (h w hw) (hv u w hu hw hx hy)) t

def tsdVertexCexX : Fin 5 → ℝ := ![0, 4, 5, 1, 4]

def tsdVertexCexY : Fin 4 → ℝ := ![4, 2, 4, 0]

theorem tsdVertexCex_le_at_knots :
    ∀ t, IsKnot tsdVertexCexX tsdVertexCexY t →
      LPM 2 t tsdVertexCexX ≤ LPM 2 t tsdVertexCexY := by
  rintro t (⟨i, rfl⟩ | ⟨j, rfl⟩)
  · fin_cases i <;>
      norm_num [LPM, lowerDev, Fin.sum_univ_five, Fin.sum_univ_four, tsdVertexCexX,
        tsdVertexCexY]
  · fin_cases j <;>
      norm_num [LPM, lowerDev, Fin.sum_univ_five, Fin.sum_univ_four, tsdVertexCexX,
        tsdVertexCexY]

theorem tsdVertexCex_mean_le : mean tsdVertexCexY ≤ mean tsdVertexCexX := by
  norm_num [mean, Fin.sum_univ_five, Fin.sum_univ_four, tsdVertexCexX, tsdVertexCexY]

theorem tsdVertexCex_not_le_at_three :
    LPM 2 3 tsdVertexCexY < LPM 2 3 tsdVertexCexX := by
  norm_num [LPM, lowerDev, Fin.sum_univ_five, Fin.sum_univ_four, tsdVertexCexX,
    tsdVertexCexY]

end

end NNS
