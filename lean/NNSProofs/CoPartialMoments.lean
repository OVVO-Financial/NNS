-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.PartialMoments

namespace NNS

open Finset

variable {n : ℕ}

@[expose] public noncomputable section

def CoLPM (dx dy : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, lowerDev dx tx (x i) * lowerDev dy ty (y i)) / n

def CoUPM (dx dy : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, upperDev dx tx (x i) * upperDev dy ty (y i)) / n

def DLPM (dLower dUpper : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, upperDev dUpper tx (x i) * lowerDev dLower ty (y i)) / n

def DUPM (dLower dUpper : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, lowerDev dLower tx (x i) * upperDev dUpper ty (y i)) / n

def coMoment (tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, (x i - tx) * (y i - ty)) / n

def covariance (x y : Fin n → ℝ) : ℝ := coMoment (mean x) (mean y) x y

lemma upperDev_sub_lowerDev_one (t v : ℝ) : upperDev 1 t v - lowerDev 1 t v = v - t := by
  have h := upperDev_add_neg_one_pow_mul_lowerDev 1 t v
  rw [pow_one, pow_one, neg_one_mul, ← sub_eq_add_neg] at h
  exact h

lemma lowerDev_zero_add_upperDev_zero (t v : ℝ) : lowerDev 0 t v + upperDev 0 t v = 1 := by
  have h := upperDev_add_neg_one_pow_mul_lowerDev 0 t v
  rw [pow_zero, one_mul, pow_zero, add_comm] at h
  exact h

theorem coMoment_decomposition (tx ty : ℝ) (x y : Fin n → ℝ) :
    CoUPM 1 1 tx ty x y + CoLPM 1 1 tx ty x y - DLPM 1 1 tx ty x y - DUPM 1 1 tx ty x y
      = coMoment tx ty x y := by
  unfold CoUPM CoLPM DLPM DUPM coMoment
  rw [← add_div, ← sub_div, ← sub_div, ← sum_add_distrib, ← sum_sub_distrib,
    ← sum_sub_distrib]
  congr 1
  refine sum_congr rfl fun i _ => ?_
  rw [← upperDev_sub_lowerDev_one tx (x i), ← upperDev_sub_lowerDev_one ty (y i)]
  ring

theorem covariance_decomposition (x y : Fin n → ℝ) :
    CoUPM 1 1 (mean x) (mean y) x y + CoLPM 1 1 (mean x) (mean y) x y
      - DLPM 1 1 (mean x) (mean y) x y - DUPM 1 1 (mean x) (mean y) x y
      = covariance x y :=
  coMoment_decomposition (mean x) (mean y) x y

theorem quadrants_sum_eq_one (hn : 0 < n) (tx ty : ℝ) (x y : Fin n → ℝ) :
    CoLPM 0 0 tx ty x y + CoUPM 0 0 tx ty x y + DLPM 0 0 tx ty x y + DUPM 0 0 tx ty x y
      = 1 := by
  unfold CoLPM CoUPM DLPM DUPM
  rw [← add_div, ← add_div, ← add_div, ← sum_add_distrib, ← sum_add_distrib,
    ← sum_add_distrib]
  have hterm : ∀ i : Fin n,
      lowerDev 0 tx (x i) * lowerDev 0 ty (y i) + upperDev 0 tx (x i) * upperDev 0 ty (y i)
        + upperDev 0 tx (x i) * lowerDev 0 ty (y i)
        + lowerDev 0 tx (x i) * upperDev 0 ty (y i) = 1 := fun i => by
    have hx := lowerDev_zero_add_upperDev_zero tx (x i)
    have hy := lowerDev_zero_add_upperDev_zero ty (y i)
    calc _ = (lowerDev 0 tx (x i) + upperDev 0 tx (x i))
              * (lowerDev 0 ty (y i) + upperDev 0 ty (y i)) := by ring
      _ = 1 := by rw [hx, hy, one_mul]
  rw [sum_congr rfl fun i _ => hterm i, sum_const, card_univ, Fintype.card_fin,
    nsmul_eq_mul, mul_one, div_self (Nat.cast_ne_zero.mpr hn.ne')]

end

end NNS
