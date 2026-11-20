-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.CoPartialMoments
public import Mathlib.LinearAlgebra.Matrix.PosDef
public import Mathlib.Algebra.Order.Star.Real

namespace NNS

open Finset Matrix

variable {n k : ℕ}

@[expose] public noncomputable section

lemma sq_sum_div_le (f g : Fin n → ℝ) :
    ((∑ i, f i * g i) / n) ^ 2 ≤ ((∑ i, f i ^ 2) / n) * ((∑ i, g i ^ 2) / n) := by
  rw [div_pow, div_mul_div_comm, ← sq]
  exact div_le_div_of_nonneg_right (sum_mul_sq_le_sq_mul_sq _ f g) (sq_nonneg _)

def colOf (X : Matrix (Fin n) (Fin k) ℝ) (j : Fin k) : Fin n → ℝ := fun r => X r j

def colMeans (X : Matrix (Fin n) (Fin k) ℝ) : Fin k → ℝ := fun j => mean (colOf X j)

def pmScale (popAdj : Bool) (n dL dU : ℕ) : ℝ :=
  if popAdj = true ∧ 1 < n ∧ 0 < dL ∧ 0 < dU then (n : ℝ) / ((n : ℝ) - 1) else 1

def lowerMat (dL : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) :
    Matrix (Fin n) (Fin k) ℝ :=
  Matrix.of fun r j => lowerDev dL (t j) (X r j)

def upperMat (dU : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) :
    Matrix (Fin n) (Fin k) ℝ :=
  Matrix.of fun r j => upperDev dU (t j) (X r j)

def clpmM (dL dU : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (popAdj : Bool) :
    Matrix (Fin k) (Fin k) ℝ :=
  Matrix.of fun i j => pmScale popAdj n dL dU * CoLPM dL dL (t i) (t j) (colOf X i) (colOf X j)

def cupmM (dL dU : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (popAdj : Bool) :
    Matrix (Fin k) (Fin k) ℝ :=
  Matrix.of fun i j => pmScale popAdj n dL dU * CoUPM dU dU (t i) (t j) (colOf X i) (colOf X j)

def dupmM (dL dU : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (popAdj : Bool) :
    Matrix (Fin k) (Fin k) ℝ :=
  Matrix.of fun i j => pmScale popAdj n dL dU * DUPM dL dU (t i) (t j) (colOf X i) (colOf X j)

def dlpmM (dL dU : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (popAdj : Bool) :
    Matrix (Fin k) (Fin k) ℝ :=
  Matrix.of fun i j => pmScale popAdj n dL dU * DLPM dL dU (t i) (t j) (colOf X i) (colOf X j)

def covM (dL dU : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (popAdj : Bool) :
    Matrix (Fin k) (Fin k) ℝ :=
  cupmM dL dU t X popAdj + clpmM dL dU t X popAdj - dupmM dL dU t X popAdj
    - dlpmM dL dU t X popAdj

def pmTotal (dL dU : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (popAdj : Bool) :
    Matrix (Fin k) (Fin k) ℝ :=
  cupmM dL dU t X popAdj + dupmM dL dU t X popAdj + dlpmM dL dU t X popAdj
    + clpmM dL dU t X popAdj

def normQ (total q : ℝ) : ℝ := if 0 < total then q / total else 0

def normMat (T M : Matrix (Fin k) (Fin k) ℝ) : Matrix (Fin k) (Fin k) ℝ :=
  Matrix.of fun i j => normQ (T i j) (M i j)

-- Use degree-zero frequencies when all four moments vanish (man/PM.matrix.Rd:20 at 98bc42ea).
def normQ0 (total q q0 : ℝ) : ℝ := if 0 < total then q / total else if total = 0 then q0 else 0

def normMat0 (T M M0 : Matrix (Fin k) (Fin k) ℝ) : Matrix (Fin k) (Fin k) ℝ :=
  Matrix.of fun i j => normQ0 (T i j) (M i j) (M0 i j)

structure PMResult (k : ℕ) where
  cupm : Matrix (Fin k) (Fin k) ℝ
  dupm : Matrix (Fin k) (Fin k) ℝ
  dlpm : Matrix (Fin k) (Fin k) ℝ
  clpm : Matrix (Fin k) (Fin k) ℝ
  cov : Matrix (Fin k) (Fin k) ℝ

def PMmatrix_v8898017 (dL dU : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ)
    (popAdj norm : Bool) : PMResult k :=
  if norm then
    let T := pmTotal dL dU t X popAdj
    { cupm := normMat T (cupmM dL dU t X popAdj)
      dupm := normMat T (dupmM dL dU t X popAdj)
      dlpm := normMat T (dlpmM dL dU t X popAdj)
      clpm := normMat T (clpmM dL dU t X popAdj)
      cov := normMat T (cupmM dL dU t X popAdj) + normMat T (clpmM dL dU t X popAdj)
        - normMat T (dupmM dL dU t X popAdj) - normMat T (dlpmM dL dU t X popAdj) }
  else
    { cupm := cupmM dL dU t X popAdj
      dupm := dupmM dL dU t X popAdj
      dlpm := dlpmM dL dU t X popAdj
      clpm := clpmM dL dU t X popAdj
      cov := covM dL dU t X popAdj }

def PMmatrix (dL dU : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ)
    (popAdj norm : Bool) : PMResult k :=
  if norm then
    let T := pmTotal dL dU t X popAdj
    { cupm := normMat0 T (cupmM dL dU t X popAdj) (cupmM 0 0 t X false)
      dupm := normMat0 T (dupmM dL dU t X popAdj) (dupmM 0 0 t X false)
      dlpm := normMat0 T (dlpmM dL dU t X popAdj) (dlpmM 0 0 t X false)
      clpm := normMat0 T (clpmM dL dU t X popAdj) (clpmM 0 0 t X false)
      cov := normMat0 T (cupmM dL dU t X popAdj) (cupmM 0 0 t X false)
        + normMat0 T (clpmM dL dU t X popAdj) (clpmM 0 0 t X false)
        - normMat0 T (dupmM dL dU t X popAdj) (dupmM 0 0 t X false)
        - normMat0 T (dlpmM dL dU t X popAdj) (dlpmM 0 0 t X false) }
  else
    { cupm := cupmM dL dU t X popAdj
      dupm := dupmM dL dU t X popAdj
      dlpm := dlpmM dL dU t X popAdj
      clpm := clpmM dL dU t X popAdj
      cov := covM dL dU t X popAdj }

theorem pmScale_false (n dL dU : ℕ) : pmScale false n dL dU = 1 := by
  simp [pmScale]

theorem pmScale_true (h : 1 < n) {dL dU : ℕ} (hL : 0 < dL) (hU : 0 < dU) :
    pmScale true n dL dU = (n : ℝ) / ((n : ℝ) - 1) := by
  simp [pmScale, h, hL, hU]

section entries

variable (dL dU : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (popAdj : Bool)

theorem covM_apply (i j : Fin k) :
    covM dL dU t X popAdj i j = pmScale popAdj n dL dU *
      (CoUPM dU dU (t i) (t j) (colOf X i) (colOf X j)
        + CoLPM dL dL (t i) (t j) (colOf X i) (colOf X j)
        - DUPM dL dU (t i) (t j) (colOf X i) (colOf X j)
        - DLPM dL dU (t i) (t j) (colOf X i) (colOf X j)) := by
  simp only [covM, cupmM, clpmM, dupmM, dlpmM, Matrix.sub_apply, Matrix.add_apply, of_apply]
  ring

end entries

theorem covM_one_apply (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (popAdj : Bool)
    (i j : Fin k) :
    covM 1 1 t X popAdj i j = pmScale popAdj n 1 1 * coMoment (t i) (t j) (colOf X i) (colOf X j)
    := by
  rw [covM_apply, ← coMoment_decomposition]
  ring

def popCovMatrix (X : Matrix (Fin n) (Fin k) ℝ) : Matrix (Fin k) (Fin k) ℝ :=
  Matrix.of fun i j => covariance (colOf X i) (colOf X j)

def sampleCov (x y : Fin n → ℝ) : ℝ :=
  (∑ i, (x i - mean x) * (y i - mean y)) / ((n : ℝ) - 1)

def sampleCovMatrix (X : Matrix (Fin n) (Fin k) ℝ) : Matrix (Fin k) (Fin k) ℝ :=
  Matrix.of fun i j => sampleCov (colOf X i) (colOf X j)

theorem PMmatrix_cov_eq_popCov (X : Matrix (Fin n) (Fin k) ℝ) :
    (PMmatrix 1 1 (colMeans X) X false false).cov = popCovMatrix X := by
  ext i j
  simp only [PMmatrix, Bool.false_eq_true, ↓reduceIte, covM_one_apply, pmScale_false,
    one_mul, popCovMatrix, of_apply, covariance, colMeans]

theorem PMmatrix_cov_eq_sampleCov (hn : 1 < n) (X : Matrix (Fin n) (Fin k) ℝ) :
    (PMmatrix 1 1 (colMeans X) X true false).cov = sampleCovMatrix X := by
  ext i j
  simp only [PMmatrix, Bool.false_eq_true, ↓reduceIte, covM_one_apply,
    pmScale_true hn one_pos one_pos, sampleCovMatrix, of_apply, sampleCov, coMoment,
    colMeans]
  have h1 : (1 : ℝ) < n := by exact_mod_cast hn
  have h2 : (n : ℝ) - 1 ≠ 0 := by linarith
  have h3 : (n : ℝ) ≠ 0 := by linarith
  field_simp

section gram

variable (dL dU : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (popAdj : Bool)

theorem clpmM_eq_gram :
    clpmM dL dU t X popAdj
      = (pmScale popAdj n dL dU / n) • ((lowerMat dL t X)ᵀ * lowerMat dL t X) := by
  ext i j
  simp only [clpmM, CoLPM, lowerMat, colOf, Matrix.smul_apply, mul_apply, transpose_apply, of_apply,
    smul_eq_mul]
  ring

theorem cupmM_eq_gram :
    cupmM dL dU t X popAdj
      = (pmScale popAdj n dL dU / n) • ((upperMat dU t X)ᵀ * upperMat dU t X) := by
  ext i j
  simp only [cupmM, CoUPM, upperMat, colOf, Matrix.smul_apply, mul_apply, transpose_apply, of_apply,
    smul_eq_mul]
  ring

lemma quadForm_gram (c : ℝ) (A : Matrix (Fin n) (Fin k) ℝ) (v : Fin k → ℝ) :
    v ⬝ᵥ ((c • (Aᵀ * A)) *ᵥ v) = c * ∑ r, (∑ j, v j * A r j) ^ 2 := by
  rw [smul_mulVec, ← mulVec_mulVec, dotProduct_smul, dotProduct_mulVec, vecMul_transpose,
    smul_eq_mul]
  congr 1
  simp only [dotProduct, mulVec, sq]
  refine sum_congr rfl fun r _ => ?_
  simp only [mul_comm (v _)]

lemma posSemidef_gram {c : ℝ} (hc : 0 ≤ c) (A : Matrix (Fin n) (Fin k) ℝ) :
    (c • (Aᵀ * A)).PosSemidef := by
  have h := posSemidef_conjTranspose_mul_self A
  rw [conjTranspose_eq_transpose_of_trivial] at h
  exact h.smul hc

end gram

section output

variable (dL dU : ℕ) (t : Fin k → ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (popAdj : Bool)

theorem PMmatrix_cov_eq (norm : Bool) :
    (PMmatrix dL dU t X popAdj norm).cov = (PMmatrix dL dU t X popAdj norm).cupm
      + (PMmatrix dL dU t X popAdj norm).clpm - (PMmatrix dL dU t X popAdj norm).dupm
      - (PMmatrix dL dU t X popAdj norm).dlpm := by
  cases norm <;> rfl

end output

end

end NNS
