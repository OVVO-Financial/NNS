-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.PartialMoments
public import Mathlib.Analysis.SpecialFunctions.Pow.Real

namespace NNS

open Finset

@[expose] public noncomputable section

variable {n : ℕ}

def popVariance (x : Fin n → ℝ) : ℝ := centralMoment 2 (mean x) x

def popSkewness (x : Fin n → ℝ) : ℝ :=
  centralMoment 3 (mean x) x / popVariance x ^ ((3 : ℝ) / 2)

def popKurtosis (x : Fin n → ℝ) : ℝ :=
  centralMoment 4 (mean x) x / popVariance x ^ 2 - 3

def sampleVariance (x : Fin n → ℝ) : ℝ := popVariance x * (n / ((n : ℝ) - 1))

def sampleSkewness (x : Fin n → ℝ) : ℝ :=
  (n / (((n : ℝ) - 1) * ((n : ℝ) - 2))) *
    ((n * centralMoment 3 (mean x) x) /
      (popVariance x * (n / ((n : ℝ) - 1))) ^ ((3 : ℝ) / 2))

def sampleKurtosis (x : Fin n → ℝ) : ℝ :=
  ((n * ((n : ℝ) + 1)) / (((n : ℝ) - 1) * ((n : ℝ) - 2) * ((n : ℝ) - 3))) *
      ((n * centralMoment 4 (mean x) x) / (popVariance x * (n / ((n : ℝ) - 1))) ^ 2) -
    (3 * ((n : ℝ) - 1) ^ 2) / (((n : ℝ) - 2) * ((n : ℝ) - 3))

structure Moments where
  mean : ℝ
  variance : ℝ
  skewness : ℝ
  kurtosis : ℝ

def moments (x : Fin n → ℝ) (population : Bool := true) : Moments :=
  if population then ⟨mean x, popVariance x, popSkewness x, popKurtosis x⟩
  else ⟨mean x, sampleVariance x, sampleSkewness x, sampleKurtosis x⟩

def momentsPM (x : Fin n → ℝ) (population : Bool := true) : Moments :=
  let μ := mean x
  let V := LPM 2 μ x + UPM 2 μ x
  let S := UPM 3 μ x - LPM 3 μ x
  let K := UPM 4 μ x + LPM 4 μ x
  let c : ℝ := n / ((n : ℝ) - 1)
  if population then ⟨μ, V, S / V ^ ((3 : ℝ) / 2), K / V ^ 2 - 3⟩
  else ⟨μ, V * c,
    (n / (((n : ℝ) - 1) * ((n : ℝ) - 2))) * ((n * S) / (V * c) ^ ((3 : ℝ) / 2)),
    ((n * ((n : ℝ) + 1)) / (((n : ℝ) - 1) * ((n : ℝ) - 2) * ((n : ℝ) - 3))) *
        ((n * K) / (V * c) ^ 2) -
      (3 * ((n : ℝ) - 1) ^ 2) / (((n : ℝ) - 2) * ((n : ℝ) - 3))⟩

theorem centralMoment_three_eq (t : ℝ) (x : Fin n → ℝ) :
    centralMoment 3 t x = UPM 3 t x - LPM 3 t x := by
  rw [← UPM_add_neg_one_pow_mul_LPM]
  ring

theorem centralMoment_four_eq (t : ℝ) (x : Fin n → ℝ) :
    centralMoment 4 t x = UPM 4 t x + LPM 4 t x := by
  rw [← UPM_add_neg_one_pow_mul_LPM]
  ring

theorem centralMoment_two_eq (t : ℝ) (x : Fin n → ℝ) :
    centralMoment 2 t x = LPM 2 t x + UPM 2 t x := by
  rw [← UPM_add_neg_one_pow_mul_LPM]
  ring

theorem UPM_one_mean_eq_LPM_one_mean (x : Fin n → ℝ) :
    UPM 1 (mean x) x = LPM 1 (mean x) x := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [UPM, LPM]
  · linarith [UPM_one_sub_LPM_one hn (mean x) x]

theorem popVariance_eq_PM (x : Fin n → ℝ) :
    popVariance x = LPM 2 (mean x) x + UPM 2 (mean x) x :=
  centralMoment_two_eq _ x

theorem popSkewness_eq_PM (x : Fin n → ℝ) :
    popSkewness x = (UPM 3 (mean x) x - LPM 3 (mean x) x) /
      (LPM 2 (mean x) x + UPM 2 (mean x) x) ^ ((3 : ℝ) / 2) := by
  rw [popSkewness, centralMoment_three_eq, popVariance_eq_PM]

theorem popKurtosis_eq_PM (x : Fin n → ℝ) :
    popKurtosis x = (UPM 4 (mean x) x + LPM 4 (mean x) x) /
      (LPM 2 (mean x) x + UPM 2 (mean x) x) ^ 2 - 3 := by
  rw [popKurtosis, centralMoment_four_eq, popVariance_eq_PM]

theorem moments_eq_momentsPM (x : Fin n → ℝ) (population : Bool) :
    moments x population = momentsPM x population := by
  cases population <;>
    simp only [moments, momentsPM, Bool.false_eq_true, ↓reduceIte, sampleVariance,
      sampleSkewness, sampleKurtosis, popSkewness, popKurtosis, popVariance_eq_PM,
      centralMoment_three_eq, centralMoment_four_eq]

theorem popVariance_nonneg (x : Fin n → ℝ) : 0 ≤ popVariance x :=
  div_nonneg (sum_nonneg fun _ _ => sq_nonneg _) n.cast_nonneg

theorem sampleVariance_eq (x : Fin n → ℝ) :
    sampleVariance x = (∑ i, (x i - mean x) ^ 2) / ((n : ℝ) - 1) := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [sampleVariance, popVariance, centralMoment]
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  rw [sampleVariance, popVariance, centralMoment]
  rw [div_mul_div_comm, mul_comm (∑ i, (x i - mean x) ^ 2), mul_div_mul_left _ _ hn']

theorem rpow_three_halves {v : ℝ} (hv : 0 ≤ v) : v ^ ((3 : ℝ) / 2) = Real.sqrt v ^ 3 := by
  rw [show (3 : ℝ) / 2 = 1 / 2 * (3 : ℕ) by norm_num, Real.rpow_mul hv, Real.rpow_natCast,
    Real.sqrt_eq_rpow]

theorem sampleSkewness_eq_G1 (hn : 2 < n) (x : Fin n → ℝ) :
    sampleSkewness x =
      Real.sqrt ((n : ℝ) * ((n : ℝ) - 1)) / ((n : ℝ) - 2) * popSkewness x := by
  have hn1 : (1 : ℝ) < n := by exact_mod_cast (by omega : 1 < n)
  have hn2 : (2 : ℝ) < n := by exact_mod_cast hn
  have hv := popVariance_nonneg x
  have hc : (0 : ℝ) ≤ n / ((n : ℝ) - 1) := by positivity
  rw [sampleSkewness, popSkewness, rpow_three_halves (mul_nonneg hv hc),
    rpow_three_halves hv, Real.sqrt_mul hv]
  set m3 := centralMoment 3 (mean x) x
  set sv := Real.sqrt (popVariance x)
  set sc := Real.sqrt ((n : ℝ) / ((n : ℝ) - 1))
  have hsc : 0 < sc := Real.sqrt_pos.mpr (by positivity)
  have hsc2 : sc ^ 2 = (n : ℝ) / ((n : ℝ) - 1) := Real.sq_sqrt hc
  have hs1 : Real.sqrt ((n : ℝ) * ((n : ℝ) - 1)) = n / sc := by
    rw [eq_div_iff hsc.ne', ← Real.sqrt_mul (by positivity)]
    rw [show (n : ℝ) * ((n : ℝ) - 1) * ((n : ℝ) / ((n : ℝ) - 1)) = (n : ℝ) ^ 2 by
      field_simp]
    exact Real.sqrt_sq (by positivity)
  rw [hs1]
  have hsv0 : 0 ≤ sv := Real.sqrt_nonneg _
  rcases hsv0.eq_or_lt with h0 | hsv
  · rw [← h0]; simp
  have h1 : (n : ℝ) - 1 ≠ 0 := by linarith
  have h2 : (n : ℝ) - 2 ≠ 0 := by linarith
  field_simp
  rw [hsc2]
  field_simp

theorem sampleKurtosis_eq_G2 (hn : 3 < n) (x : Fin n → ℝ) :
    sampleKurtosis x =
      ((n : ℝ) - 1) / (((n : ℝ) - 2) * ((n : ℝ) - 3)) *
        (((n : ℝ) + 1) * popKurtosis x + 6) := by
  have hn3 : (3 : ℝ) < n := by exact_mod_cast hn
  have h0 : (n : ℝ) ≠ 0 := by linarith
  have h1 : (n : ℝ) - 1 ≠ 0 := by linarith
  have h2 : (n : ℝ) - 2 ≠ 0 := by linarith
  have h3 : (n : ℝ) - 3 ≠ 0 := by linarith
  rw [sampleKurtosis, popKurtosis]
  rcases eq_or_ne (popVariance x) 0 with hv | hv
  · rw [hv]
    simp only [zero_mul, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, div_zero,
      mul_zero, zero_sub]
    field_simp
    ring
  field_simp
  ring

theorem popSkewness_neg_iff (x : Fin n → ℝ) (hv : 0 < popVariance x) :
    popSkewness x < 0 ↔ UPM 3 (mean x) x < LPM 3 (mean x) x := by
  have hd : 0 < popVariance x ^ ((3 : ℝ) / 2) := by positivity
  rw [popSkewness, centralMoment_three_eq, div_neg_iff, sub_neg]
  constructor
  · rintro (⟨_, h⟩ | ⟨h, _⟩)
    · exact absurd h (not_lt.mpr hd.le)
    · linarith
  · intro h
    exact Or.inr ⟨by linarith, hd⟩

section Affine

variable (a c : ℝ) (x : Fin n → ℝ)

theorem mean_affine (hn : 0 < n) : mean (fun i => a * x i + c) = a * mean x + c := by
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  simp only [mean, sum_add_distrib, ← mul_sum, sum_const, card_univ, Fintype.card_fin,
    nsmul_eq_mul]
  field_simp

theorem centralMoment_mean_affine (k : ℕ) :
    centralMoment k (mean fun i => a * x i + c) (fun i => a * x i + c) =
      a ^ k * centralMoment k (mean x) x := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [centralMoment]
  simp only [centralMoment, mean_affine a c x hn, mul_div_assoc', mul_sum, ← mul_pow]
  congr 1
  exact sum_congr rfl fun i _ => by ring

theorem popVariance_affine :
    popVariance (fun i => a * x i + c) = a ^ 2 * popVariance x :=
  centralMoment_mean_affine a c x 2

theorem popKurtosis_affine (ha : a ≠ 0) :
    popKurtosis (fun i => a * x i + c) = popKurtosis x := by
  rw [popKurtosis, popKurtosis, centralMoment_mean_affine, popVariance_affine, mul_pow,
    ← pow_mul, mul_div_mul_left _ _ (pow_ne_zero _ ha)]

end Affine

theorem moments_population (x : Fin n → ℝ) :
    moments x true = ⟨mean x, centralMoment 2 (mean x) x,
      centralMoment 3 (mean x) x / centralMoment 2 (mean x) x ^ ((3 : ℝ) / 2),
      centralMoment 4 (mean x) x / centralMoment 2 (mean x) x ^ 2 - 3⟩ := rfl

end

end NNS
