-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.PartialMoments
public import Mathlib.Analysis.SpecialFunctions.Pow.Continuity

namespace NNS

open Finset Filter Topology

@[expose] public noncomputable section

variable {n : ℕ}

def rLowerDev (d t v : ℝ) : ℝ := if v ≤ t then (t - v) ^ d else 0

def rUpperDev (d t v : ℝ) : ℝ := if t < v then (v - t) ^ d else 0

def rLPM (d t : ℝ) (x : Fin n → ℝ) : ℝ := (∑ i, rLowerDev d t (x i)) / n

def rUPM (d t : ℝ) (x : Fin n → ℝ) : ℝ := (∑ i, rUpperDev d t (x i)) / n

lemma rLowerDev_nonneg (d t v : ℝ) : 0 ≤ rLowerDev d t v := by
  unfold rLowerDev
  split_ifs with h
  · exact Real.rpow_nonneg (sub_nonneg.mpr h) d
  · exact le_rfl

lemma rUpperDev_nonneg (d t v : ℝ) : 0 ≤ rUpperDev d t v := by
  unfold rUpperDev
  split_ifs with h
  · exact Real.rpow_nonneg (sub_nonneg.mpr h.le) d
  · exact le_rfl

lemma rLowerDev_zero (t v : ℝ) : rLowerDev 0 t v = if v ≤ t then 1 else 0 := by
  simp [rLowerDev]

lemma rUpperDev_zero (t v : ℝ) : rUpperDev 0 t v = if t < v then 1 else 0 := by
  simp [rUpperDev]

lemma rLowerDev_natCast (d : ℕ) (t v : ℝ) : rLowerDev d t v = lowerDev d t v := by
  simp [rLowerDev, lowerDev, Real.rpow_natCast]

lemma rUpperDev_natCast (d : ℕ) (t v : ℝ) : rUpperDev d t v = upperDev d t v := by
  simp [rUpperDev, upperDev, Real.rpow_natCast]

lemma rLowerDev_eq_max {d : ℝ} (hd : d ≠ 0) (t v : ℝ) :
    rLowerDev d t v = max (t - v) 0 ^ d := by
  unfold rLowerDev
  split_ifs with h
  · rw [max_eq_left (sub_nonneg.mpr h)]
  · rw [max_eq_right (by linarith [not_le.mp h]), Real.zero_rpow hd]

lemma rUpperDev_eq_max {d : ℝ} (hd : d ≠ 0) (t v : ℝ) :
    rUpperDev d t v = max (v - t) 0 ^ d := by
  unfold rUpperDev
  split_ifs with h
  · rw [max_eq_left (sub_nonneg.mpr h.le)]
  · rw [max_eq_right (by linarith [not_lt.mp h]), Real.zero_rpow hd]

theorem rLPM_natCast (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : rLPM d t x = LPM d t x := by
  unfold rLPM LPM
  simp only [rLowerDev_natCast]

theorem rUPM_natCast (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : rUPM d t x = UPM d t x := by
  unfold rUPM UPM
  simp only [rUpperDev_natCast]

theorem rLPM_zero (t : ℝ) (x : Fin n → ℝ) : rLPM 0 t x = LPM 0 t x := by
  exact_mod_cast rLPM_natCast 0 t x

theorem rUPM_zero (t : ℝ) (x : Fin n → ℝ) : rUPM 0 t x = UPM 0 t x := by
  exact_mod_cast rUPM_natCast 0 t x

theorem rLPM_zero_add_rUPM_zero (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    rLPM 0 t x + rUPM 0 t x = 1 := by
  rw [rLPM_zero, rUPM_zero, LPM_zero_add_UPM_zero hn]

theorem rLPM_nonneg (d t : ℝ) (x : Fin n → ℝ) : 0 ≤ rLPM d t x :=
  div_nonneg (sum_nonneg fun i _ => rLowerDev_nonneg d t (x i)) n.cast_nonneg

theorem rUPM_nonneg (d t : ℝ) (x : Fin n → ℝ) : 0 ≤ rUPM d t x :=
  div_nonneg (sum_nonneg fun i _ => rUpperDev_nonneg d t (x i)) n.cast_nonneg

lemma rLowerDev_add_rUpperDev (d t v : ℝ) : rLowerDev d t v + rUpperDev d t v = |t - v| ^ d := by
  unfold rLowerDev rUpperDev
  rcases lt_trichotomy v t with h | h | h
  · rw [ite_eq_left h.le, ite_eq_right (not_lt.mpr h.le), add_zero,
      abs_of_pos (sub_pos.mpr h)]
  · subst h
    simp
  · rw [ite_eq_right (not_le.mpr h), ite_eq_left h, zero_add, abs_of_neg (sub_neg.mpr h),
      neg_sub]

theorem rLPM_add_rUPM (d t : ℝ) (x : Fin n → ℝ) :
    rLPM d t x + rUPM d t x = (∑ i, |x i - t| ^ d) / n := by
  unfold rLPM rUPM
  rw [← add_div, ← sum_add_distrib]
  congr 1
  exact sum_congr rfl fun i _ => by rw [rLowerDev_add_rUpperDev, abs_sub_comm]

theorem rLPM_eq_max {d : ℝ} (hd : d ≠ 0) (t : ℝ) (x : Fin n → ℝ) :
    rLPM d t x = (∑ i, max (t - x i) 0 ^ d) / n := by
  simp only [rLPM, rLowerDev_eq_max hd]

theorem rUPM_eq_max {d : ℝ} (hd : d ≠ 0) (t : ℝ) (x : Fin n → ℝ) :
    rUPM d t x = (∑ i, max (x i - t) 0 ^ d) / n := by
  simp only [rUPM, rUpperDev_eq_max hd]

theorem rUPM_eq_rLPM_neg {d : ℝ} (hd : d ≠ 0) (t : ℝ) (x : Fin n → ℝ) :
    rUPM d t x = rLPM d (-t) (fun i => -x i) := by
  rw [rUPM_eq_max hd, rLPM_eq_max hd]
  congr 1
  exact sum_congr rfl fun i _ => by ring_nf

theorem rLPM_mul_const {a : ℝ} (ha : 0 < a) (d t : ℝ) (x : Fin n → ℝ) :
    rLPM d (a * t) (fun i => a * x i) = a ^ d * rLPM d t x := by
  unfold rLPM rLowerDev
  rw [mul_div_assoc', mul_sum]
  congr 1
  refine sum_congr rfl fun i _ => ?_
  by_cases h : x i ≤ t
  · have h' : a * x i ≤ a * t := mul_le_mul_of_nonneg_left h ha.le
    rw [ite_eq_left h', ite_eq_left h, ← mul_sub, Real.mul_rpow ha.le (sub_nonneg.mpr h)]
  · have h' : ¬ a * x i ≤ a * t := fun h' => h (le_of_mul_le_mul_left h' ha)
    rw [ite_eq_right h', ite_eq_right h, mul_zero]

theorem rUPM_mul_const {a : ℝ} (ha : 0 < a) (d t : ℝ) (x : Fin n → ℝ) :
    rUPM d (a * t) (fun i => a * x i) = a ^ d * rUPM d t x := by
  unfold rUPM rUpperDev
  rw [mul_div_assoc', mul_sum]
  congr 1
  refine sum_congr rfl fun i _ => ?_
  by_cases h : t < x i
  · have h' : a * t < a * x i := mul_lt_mul_of_pos_left h ha
    rw [ite_eq_left h', ite_eq_left h, ← mul_sub, Real.mul_rpow ha.le (sub_nonneg.mpr h.le)]
  · have h' : ¬ a * t < a * x i := fun h' => h (lt_of_mul_lt_mul_left h' ha.le)
    rw [ite_eq_right h', ite_eq_right h, mul_zero]

theorem rLPM_const (hn : 0 < n) (d t c : ℝ) :
    rLPM d t (fun _ : Fin n => c) = rLowerDev d t c := by
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  simp only [rLPM, sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul]
  field_simp

theorem rUPM_const (hn : 0 < n) (d t c : ℝ) :
    rUPM d t (fun _ : Fin n => c) = rUpperDev d t c := by
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  simp only [rUPM, sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul]
  field_simp

theorem rLPM_continuous {d : ℝ} (hd : 0 < d) (x : Fin n → ℝ) :
    Continuous fun t => rLPM d t x := by
  simp only [rLPM_eq_max hd.ne']
  refine Continuous.div_const (continuous_finsetSum _ fun i _ => ?_) _
  exact ((continuous_id.sub continuous_const).max continuous_const).rpow_const
    fun _ => Or.inr hd.le

theorem rUPM_continuous {d : ℝ} (hd : 0 < d) (x : Fin n → ℝ) :
    Continuous fun t => rUPM d t x := by
  simp only [rUPM_eq_max hd.ne']
  refine Continuous.div_const (continuous_finsetSum _ fun i _ => ?_) _
  exact ((continuous_const.sub continuous_id).max continuous_const).rpow_const
    fun _ => Or.inr hd.le

def prefixDegreeOKExact (d : ℝ) : Prop :=
  0 ≤ d ∧ (round d : ℝ) = d ∧ d ≤ 32

open Classical in
def rLPMVecExact (d t : ℝ) (x : Fin n → ℝ) : ℝ :=
  if prefixDegreeOKExact d then LPM (round d).toNat t x else rLPM d t x

private lemma toNat_round_cast {d : ℝ} (h : prefixDegreeOKExact d) :
    (((round d).toNat : ℕ) : ℝ) = d := by
  have h0 : 0 ≤ round d := by exact_mod_cast h.2.1 ▸ h.1
  calc (((round d).toNat : ℕ) : ℝ) = ((round d : ℤ) : ℝ) := by
        exact_mod_cast Int.toNat_of_nonneg h0
    _ = d := h.2.1

theorem rLPMVecExact_eq_rLPM (d t : ℝ) (x : Fin n → ℝ) : rLPMVecExact d t x = rLPM d t x := by
  by_cases h : prefixDegreeOKExact d
  · rw [rLPMVecExact, ite_eq_left h, ← rLPM_natCast, toNat_round_cast h]
  · rw [rLPMVecExact, ite_eq_right h]

end

end NNS
