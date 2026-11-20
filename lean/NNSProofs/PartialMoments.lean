-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import Mathlib.Algebra.BigOperators.Ring.Finset
public import Mathlib.Algebra.Order.BigOperators.Group.Finset
public import Mathlib.Basic.Real.Basic
public import Mathlib.Tactic.Linarith
public import Mathlib.Tactic.Ring

namespace NNS

open Finset

@[expose] public noncomputable section

variable {n : ℕ}

-- Count ties on the lower side to recover the empirical CDF (man/LPM.Rd:10 at 98bc42ea).
def lowerDev (d : ℕ) (t v : ℝ) : ℝ := if v ≤ t then (t - v) ^ d else 0

def upperDev (d : ℕ) (t v : ℝ) : ℝ := if t < v then (v - t) ^ d else 0

def mean (x : Fin n → ℝ) : ℝ := (∑ i, x i) / n

def LPM (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : ℝ :=
  (∑ i, lowerDev d t (x i)) / n

def UPM (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : ℝ :=
  (∑ i, upperDev d t (x i)) / n

def centralMoment (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : ℝ :=
  (∑ i, (x i - t) ^ d) / n

lemma lowerDev_nonneg (d : ℕ) (t v : ℝ) : 0 ≤ lowerDev d t v := by
  unfold lowerDev
  split_ifs with h
  · exact pow_nonneg (sub_nonneg.mpr h) d
  · exact le_refl 0

lemma upperDev_nonneg (d : ℕ) (t v : ℝ) : 0 ≤ upperDev d t v := by
  unfold upperDev
  split_ifs with h
  · exact pow_nonneg (sub_nonneg.mpr h.le) d
  · exact le_refl 0

lemma upperDev_add_neg_one_pow_mul_lowerDev (d : ℕ) (t v : ℝ) :
    upperDev d t v + (-1) ^ d * lowerDev d t v = (v - t) ^ d := by
  unfold upperDev lowerDev
  by_cases h : t < v
  · have h' : ¬ v ≤ t := not_le.mpr h
    simp [h, h']
  · have h' : v ≤ t := not_lt.mp h
    simp only [h, h', ↓reduceIte, zero_add]
    rw [← mul_pow]
    ring

theorem UPM_add_neg_one_pow_mul_LPM (d : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    UPM d t x + (-1) ^ d * LPM d t x = centralMoment d t x := by
  unfold UPM LPM centralMoment
  rw [mul_div_assoc', ← add_div, mul_sum, ← sum_add_distrib]
  congr 1
  exact sum_congr rfl fun i _ => upperDev_add_neg_one_pow_mul_lowerDev d t (x i)

theorem LPM_zero_add_UPM_zero (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    LPM 0 t x + UPM 0 t x = 1 := by
  have h := UPM_add_neg_one_pow_mul_LPM 0 t x
  simp only [pow_zero, one_mul, centralMoment, sum_const, card_univ,
    Fintype.card_fin, nsmul_eq_mul, mul_one] at h
  rw [add_comm, h, div_self (Nat.cast_ne_zero.mpr hn.ne')]

theorem UPM_one_sub_LPM_one (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    UPM 1 t x - LPM 1 t x = mean x - t := by
  have h := UPM_add_neg_one_pow_mul_LPM 1 t x
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  simp only [pow_one, neg_one_mul, ← sub_eq_add_neg, centralMoment, sum_sub_distrib,
    sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul] at h
  rw [h, mean, sub_div, mul_div_cancel_left₀ t hn']

theorem LPM_nonneg (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : 0 ≤ LPM d t x :=
  div_nonneg (sum_nonneg fun i _ => lowerDev_nonneg d t (x i)) n.cast_nonneg

theorem UPM_nonneg (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : 0 ≤ UPM d t x :=
  div_nonneg (sum_nonneg fun i _ => upperDev_nonneg d t (x i)) n.cast_nonneg

theorem LPM_zero_eq_ecdf (t : ℝ) (x : Fin n → ℝ) :
    LPM 0 t x = (#{i | x i ≤ t} : ℝ) / n := by
  unfold LPM lowerDev
  simp only [pow_zero]
  rw [sum_boole]

theorem LPM_zero_monotone (x : Fin n → ℝ) : Monotone fun t => LPM 0 t x := by
  intro s t hst
  simp only [LPM_zero_eq_ecdf]
  gcongr

end

end NNS
