-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.RealDegree
public import NNSProofs.CoMomentsMatrix

namespace NNS

open Finset

@[expose] public noncomputable section

variable {n k m : ℕ}

def cppIntCast (d : ℝ) : ℤ := if 0 ≤ d then ⌊d⌋ else ⌈d⌉

def repeatMul (v : ℝ) (p : ℤ) : ℝ := (fun r => r * v)^[p.toNat] 1

def cppPow (v d : ℝ) : ℝ := if d = (cppIntCast d : ℝ) then repeatMul v (cppIntCast d) else v ^ d

lemma repeatMul_natCast (v : ℝ) (p : ℕ) : repeatMul v p = v ^ p := by
  unfold repeatMul
  rw [Int.toNat_natCast]
  induction p with
  | zero => simp
  | succ p ih => rw [Function.iterate_succ_apply', ih, pow_succ]

lemma cppPow_eq_rpow {d : ℝ} (hd : 0 ≤ d) (v : ℝ) : cppPow v d = v ^ d := by
  unfold cppPow
  have hc : cppIntCast d = ⌊d⌋ := ite_eq_left hd
  split_ifs with h
  · rw [hc] at h ⊢
    obtain ⟨p, hp⟩ := Int.eq_ofNat_of_zero_le (Int.floor_nonneg.mpr hd)
    rw [hp] at h ⊢
    rw [repeatMul_natCast, h, Int.cast_natCast, Real.rpow_natCast]
  · rfl

def lowerComp (diff d : ℝ) : ℝ :=
  if d = 0 then (if 0 ≤ diff then 1 else 0) else if diff < 0 then 0 else cppPow diff d

def upperComp (diff d : ℝ) : ℝ :=
  if d = 0 then (if 0 < diff then 1 else 0) else if diff < 0 then 0 else cppPow diff d

lemma lowerComp_eq {d : ℝ} (hd : 0 ≤ d) (t v : ℝ) : lowerComp (t - v) d = rLowerDev d t v := by
  unfold lowerComp rLowerDev
  rw [cppPow_eq_rpow hd]
  by_cases h : v ≤ t
  · have h1 : 0 ≤ t - v := sub_nonneg.mpr h
    by_cases hd0 : d = 0
    · simp [hd0, h, h1]
    · simp [hd0, h, not_lt.mpr h1]
  · have h1 : t - v < 0 := by linarith [not_le.mp h]
    simp [h, h1, not_le.mpr h1]

lemma upperComp_eq {d : ℝ} (hd : 0 ≤ d) (t v : ℝ) : upperComp (v - t) d = rUpperDev d t v := by
  unfold upperComp rUpperDev
  rw [cppPow_eq_rpow hd]
  rcases lt_trichotomy v t with h | h | h
  · have h1 : v - t < 0 := by linarith
    simp [h1, not_lt.mpr h.le, not_lt.mpr h1.le]
  · subst h
    by_cases hd0 : d = 0
    · simp [hd0]
    · simp [hd0, Real.zero_rpow hd0]
  · have h1 : 0 < v - t := by linarith
    by_cases hd0 : d = 0
    · simp [hd0, h, h1]
    · simp [hd0, h, not_lt.mpr h1.le]

def allLoop {α : Type*} (ok : α → Prop) [DecidablePred ok] : List α → Bool
  | [] => true
  | a :: l => if ok a then allLoop ok l else false

lemma allLoop_eq_true {α : Type*} (ok : α → Prop) [DecidablePred ok] (l : List α) :
    allLoop ok l = true ↔ ∀ a ∈ l, ok a := by
  induction l with
  | nil => simp [allLoop]
  | cons a l ih => by_cases h : ok a <;> simp [allLoop, h, ih]

def prodLoop {α : Type*} (kill : α → Prop) [DecidablePred kill] (g : α → ℝ) :
    List α → ℝ → ℝ
  | [], p => p
  | a :: l, p => if kill a then 0 else prodLoop kill g l (p * g a)

lemma prodLoop_eq {α : Type*} (kill : α → Prop) [DecidablePred kill] (g : α → ℝ)
    (l : List α) (p : ℝ) :
    prodLoop kill g l p = p * (l.map fun a => if kill a then 0 else g a).prod := by
  induction l generalizing p with
  | nil => simp [prodLoop]
  | cons a l ih => by_cases h : kill a <;> simp [prodLoop, h, ih, mul_assoc]

def flagLoop {α : Type*} (nb na : α → Prop) [DecidablePred nb] [DecidablePred na] :
    List α → Bool × Bool → Bool × Bool
  | [], s => s
  | a :: l, s =>
    let s' := (s.1 && !decide (nb a), s.2 && !decide (na a))
    if !s'.1 && !s'.2 then s' else flagLoop nb na l s'

lemma flagLoop_eq {α : Type*} (nb na : α → Prop) [DecidablePred nb] [DecidablePred na]
    (l : List α) (s : Bool × Bool) :
    flagLoop nb na l s = (s.1 && decide (∀ a ∈ l, ¬ nb a), s.2 && decide (∀ a ∈ l, ¬ na a)) := by
  induction l generalizing s with
  | nil => simp [flagLoop]
  | cons a l ih =>
    obtain ⟨b, u⟩ := s
    by_cases h1 : nb a <;> by_cases h2 : na a <;> cases b <;> cases u <;>
      simp [flagLoop, h1, h2, ih]

lemma foldl_mul_eq {α : Type*} (g : α → ℝ) (l : List α) (p : ℝ) :
    l.foldl (fun q a => q * g a) p = p * (l.map g).prod := by
  induction l generalizing p with
  | nil => simp
  | cons a l ih => simp [ih, mul_assoc]

lemma foldl_add_eq {α : Type*} (g : α → ℝ) (l : List α) (s : ℝ) :
    l.foldl (fun q a => q + g a) s = s + (l.map g).sum := by
  induction l generalizing s with
  | nil => simp
  | cons a l ih => simp [ih, add_assoc]

def sugarSum (v : Fin n → ℝ) : ℝ := (List.finRange n).foldl (fun s i => s + v i) 0

lemma sugarSum_eq (v : Fin n → ℝ) : sugarSum v = ∑ i, v i := by
  rw [sugarSum, foldl_add_eq, zero_add, Fin.sum_univ_def]

lemma finRange_prod_ite {k : ℕ} (P : Fin k → Prop) [DecidablePred P] (g : Fin k → ℝ) :
    ((List.finRange k).map fun j => if P j then (0 : ℝ) else g j).prod
      = ∏ j, if P j then 0 else g j := (Fin.prod_univ_def _).symm

def clpmCountRow (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (i : Fin n) : ℝ :=
  if allLoop (fun j => ¬ X i j > t j) (List.finRange k) then 1 else 0

def clpmSumRow (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (i : Fin n) : ℝ :=
  prodLoop (fun j => t j - X i j < 0) (fun j => cppPow (t j - X i j) d) (List.finRange k) 1

def cupmCountRow (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (i : Fin n) : ℝ :=
  if allLoop (fun j => ¬ X i j ≤ t j) (List.finRange k) then 1 else 0

def cupmSumRow (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (i : Fin n) : ℝ :=
  prodLoop (fun j => X i j - t j < 0) (fun j => cppPow (X i j - t j) d) (List.finRange k) 1

def dpmCountRow (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (i : Fin n) : ℝ :=
  let f := flagLoop (fun j => X i j - t j > 0) (fun j => X i j - t j ≤ 0) (List.finRange k)
    (true, true)
  if !f.1 && !f.2 then 1 else 0

def dpmSumRow (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (i : Fin n) : ℝ :=
  let f := flagLoop (fun j => X i j - t j ≥ 0) (fun j => X i j - t j ≤ 0) (List.finRange k)
    (true, true)
  if f.1 || f.2 then 0 else (List.finRange k).foldl (fun p j => p * cppPow |X i j - t j| d) 1

def clpmRaw (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) : ℝ :=
  if d = 0 then sugarSum (clpmCountRow X t) / n else sugarSum (clpmSumRow X t d) / n

def cupmRaw (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) : ℝ :=
  if d = 0 then sugarSum (cupmCountRow X t) / n else sugarSum (cupmSumRow X t d) / n

def dpmRaw (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) : ℝ :=
  if d = 0 then sugarSum (dpmCountRow X t) / n else sugarSum (dpmSumRow X t d) / n

def normConstRaw (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) : ℝ :=
  clpmRaw X t d + cupmRaw X t d + dpmRaw X t d

def coLPMnD_R_v8898017 (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (norm : Bool) :
    ℝ :=
  if d = 0 then clpmRaw X t d
  else if norm then normQ (normConstRaw X t d) (clpmRaw X t d) else clpmRaw X t d

def coUPMnD_R_v8898017 (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (norm : Bool) :
    ℝ :=
  if d = 0 then cupmRaw X t d
  else if norm then normQ (normConstRaw X t d) (cupmRaw X t d) else cupmRaw X t d

def dpmND_R_v8898017 (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (norm : Bool) :
    ℝ :=
  if d = 0 then dpmRaw X t d
  else if norm then normQ (normConstRaw X t d) (dpmRaw X t d) else dpmRaw X t d

def lowRow (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (i : Fin n) : ℝ :=
  ∏ j, rLowerDev d (t j) (X i j)

def upRow (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (i : Fin n) : ℝ :=
  ∏ j, rUpperDev d (t j) (X i j)

def absRow (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (i : Fin n) : ℝ :=
  ∏ j, |X i j - t j| ^ d

def dpmRow (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (i : Fin n) : ℝ :=
  if (∀ j, X i j ≤ t j) ∨ (∀ j, t j < X i j) then 0 else absRow X t d i

def clpmND (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) : ℝ :=
  (∑ i, lowRow X t d i) / n

def cupmND (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) : ℝ :=
  (∑ i, upRow X t d i) / n

def dpmND (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) : ℝ :=
  (∑ i, dpmRow X t d i) / n

def absND (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) : ℝ :=
  (∑ i, absRow X t d i) / n

def totalND (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) : ℝ :=
  clpmND X t d + cupmND X t d + dpmND X t d

lemma clpmCountRow_eq (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (i : Fin n) :
    clpmCountRow X t i = lowRow X t 0 i := by
  unfold clpmCountRow lowRow
  simp only [rLowerDev_zero]
  rw [prod_ite_zero, prod_const_one]
  congr 1
  apply propext
  rw [allLoop_eq_true]
  simp

lemma clpmSumRow_eq {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ)
    (i : Fin n) : clpmSumRow X t d i = lowRow X t d i := by
  unfold clpmSumRow lowRow
  rw [prodLoop_eq, one_mul, finRange_prod_ite]
  refine prod_congr rfl fun j _ => ?_
  rw [cppPow_eq_rpow hd, rLowerDev]
  by_cases h : X i j ≤ t j
  · simp [h, not_lt.mpr (sub_nonneg.mpr h)]
  · simp [h, sub_neg.mpr (not_le.mp h)]

lemma cupmCountRow_eq (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (i : Fin n) :
    cupmCountRow X t i = upRow X t 0 i := by
  unfold cupmCountRow upRow
  simp only [rUpperDev_zero]
  rw [prod_ite_zero, prod_const_one]
  congr 1
  apply propext
  rw [allLoop_eq_true]
  simp

lemma cupmSumRow_eq {d : ℝ} (hd : 0 ≤ d) (hd0 : d ≠ 0) (X : Matrix (Fin n) (Fin k) ℝ)
    (t : Fin k → ℝ) (i : Fin n) : cupmSumRow X t d i = upRow X t d i := by
  unfold cupmSumRow upRow
  rw [prodLoop_eq, one_mul, finRange_prod_ite]
  refine prod_congr rfl fun j _ => ?_
  rw [cppPow_eq_rpow hd, rUpperDev]
  rcases lt_trichotomy (X i j) (t j) with h | h | h
  · simp [sub_neg.mpr h, not_lt.mpr h.le]
  · simp [h, Real.zero_rpow hd0]
  · simp [h, not_lt.mpr (sub_nonneg.mpr h.le)]

lemma dpmCountRow_eq (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (i : Fin n) :
    dpmCountRow X t i = dpmRow X t 0 i := by
  unfold dpmCountRow dpmRow absRow
  rw [flagLoop_eq]
  simp only [Real.rpow_zero, prod_const_one, Bool.true_and, gt_iff_lt, sub_pos, not_lt,
    sub_nonpos, not_le, List.mem_finRange, forall_const]
  by_cases hL : ∀ j, X i j ≤ t j <;> by_cases hU : ∀ j, t j < X i j <;> simp [hL, hU]

lemma dpmSumRow_eq {d : ℝ} (hd : 0 ≤ d) (hd0 : d ≠ 0) (X : Matrix (Fin n) (Fin k) ℝ)
    (t : Fin k → ℝ) (i : Fin n) : dpmSumRow X t d i = dpmRow X t d i := by
  unfold dpmSumRow dpmRow
  rw [flagLoop_eq, foldl_mul_eq, one_mul, ← Fin.prod_univ_def]
  simp only [Bool.true_and, ge_iff_le, sub_nonneg, not_le, sub_nonpos, not_le,
    List.mem_finRange, forall_const, cppPow_eq_rpow hd]
  change (if decide (∀ j, X i j < t j) || decide (∀ j, t j < X i j) = true then 0
    else absRow X t d i) = _
  by_cases hA : ∀ j, X i j < t j
  · have hL : ∀ j, X i j ≤ t j := fun j => (hA j).le
    simp [hA, hL]
  · by_cases hU : ∀ j, t j < X i j
    · simp [hU]
    · by_cases hL : ∀ j, X i j ≤ t j
      · obtain ⟨j, hj⟩ := not_forall.mp hA
        have hjt : X i j = t j := le_antisymm (hL j) (not_lt.mp hj)
        have h0 : absRow X t d i = 0 :=
          prod_eq_zero (mem_univ j) (by rw [hjt, sub_self, abs_zero, Real.zero_rpow hd0])
        simp [hA, hU, hL, h0]
      · simp [hA, hU, hL]

theorem clpmRaw_eq {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) :
    clpmRaw X t d = clpmND X t d := by
  unfold clpmRaw clpmND
  split_ifs with hd0
  · subst hd0
    rw [sugarSum_eq]
    simp only [clpmCountRow_eq]
  · rw [sugarSum_eq]
    simp only [clpmSumRow_eq hd]

theorem cupmRaw_eq {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) :
    cupmRaw X t d = cupmND X t d := by
  unfold cupmRaw cupmND
  split_ifs with hd0
  · subst hd0
    rw [sugarSum_eq]
    simp only [cupmCountRow_eq]
  · rw [sugarSum_eq]
    simp only [cupmSumRow_eq hd hd0]

theorem dpmRaw_eq {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) :
    dpmRaw X t d = dpmND X t d := by
  unfold dpmRaw dpmND
  split_ifs with hd0
  · subst hd0
    rw [sugarSum_eq]
    simp only [dpmCountRow_eq]
  · rw [sugarSum_eq]
    simp only [dpmSumRow_eq hd hd0]

theorem normConstRaw_eq {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) :
    normConstRaw X t d = totalND X t d := by
  rw [normConstRaw, totalND, clpmRaw_eq hd, cupmRaw_eq hd, dpmRaw_eq hd]

theorem coLPMnD_R_eq_v8898017 {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ)
    (t : Fin k → ℝ) (norm : Bool) : coLPMnD_R_v8898017 X t d norm =
      if d ≠ 0 ∧ norm = true then normQ (totalND X t d) (clpmND X t d) else clpmND X t d := by
  unfold coLPMnD_R_v8898017
  rw [normConstRaw_eq hd, clpmRaw_eq hd]
  by_cases hd0 : d = 0 <;> cases norm <;> simp [hd0]

theorem coUPMnD_R_eq_v8898017 {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ)
    (t : Fin k → ℝ) (norm : Bool) : coUPMnD_R_v8898017 X t d norm =
      if d ≠ 0 ∧ norm = true then normQ (totalND X t d) (cupmND X t d) else cupmND X t d := by
  unfold coUPMnD_R_v8898017
  rw [normConstRaw_eq hd, cupmRaw_eq hd]
  by_cases hd0 : d = 0 <;> cases norm <;> simp [hd0]

theorem dpmND_R_eq_v8898017 {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ)
    (t : Fin k → ℝ) (norm : Bool) : dpmND_R_v8898017 X t d norm =
      if d ≠ 0 ∧ norm = true then normQ (totalND X t d) (dpmND X t d) else dpmND X t d := by
  unfold dpmND_R_v8898017
  rw [normConstRaw_eq hd, dpmRaw_eq hd]
  by_cases hd0 : d = 0 <;> cases norm <;> simp [hd0]

lemma lowRow_nonneg (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (i : Fin n) :
    0 ≤ lowRow X t d i := prod_nonneg fun _ _ => rLowerDev_nonneg _ _ _

lemma upRow_nonneg (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (i : Fin n) :
    0 ≤ upRow X t d i := prod_nonneg fun _ _ => rUpperDev_nonneg _ _ _

lemma absRow_nonneg (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (i : Fin n) :
    0 ≤ absRow X t d i := prod_nonneg fun _ _ => Real.rpow_nonneg (abs_nonneg _) _

lemma dpmRow_nonneg (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (i : Fin n) :
    0 ≤ dpmRow X t d i := by
  unfold dpmRow
  split_ifs
  · exact le_rfl
  · exact absRow_nonneg X t d i

theorem clpmND_nonneg (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) :
    0 ≤ clpmND X t d := div_nonneg (sum_nonneg fun i _ => lowRow_nonneg X t d i) n.cast_nonneg

theorem cupmND_nonneg (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) :
    0 ≤ cupmND X t d := div_nonneg (sum_nonneg fun i _ => upRow_nonneg X t d i) n.cast_nonneg

theorem dpmND_nonneg (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) :
    0 ≤ dpmND X t d := div_nonneg (sum_nonneg fun i _ => dpmRow_nonneg X t d i) n.cast_nonneg

lemma lowRow_add_upRow_add_dpmRow (hk : 0 < k) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ)
    (d : ℝ) (i : Fin n) : lowRow X t d i + upRow X t d i + dpmRow X t d i = absRow X t d i := by
  unfold lowRow upRow dpmRow absRow
  by_cases hL : ∀ j, X i j ≤ t j
  · have h1 : ∏ j, rLowerDev d (t j) (X i j) = ∏ j, |X i j - t j| ^ d :=
      prod_congr rfl fun j _ => by
        rw [rLowerDev, ite_eq_left (hL j), abs_of_nonpos (sub_nonpos.mpr (hL j)), neg_sub]
    have h2 : ∏ j, rUpperDev d (t j) (X i j) = 0 :=
      prod_eq_zero (mem_univ ⟨0, hk⟩) (by rw [rUpperDev, ite_eq_right (not_lt.mpr (hL _))])
    rw [h1, h2, ite_eq_left (Or.inl hL), add_zero, add_zero]
  · obtain ⟨j0, hj0⟩ := not_forall.mp hL
    have h1 : ∏ j, rLowerDev d (t j) (X i j) = 0 :=
      prod_eq_zero (mem_univ j0) (by rw [rLowerDev, ite_eq_right hj0])
    by_cases hU : ∀ j, t j < X i j
    · have h2 : ∏ j, rUpperDev d (t j) (X i j) = ∏ j, |X i j - t j| ^ d :=
        prod_congr rfl fun j _ => by
          rw [rUpperDev, ite_eq_left (hU j), abs_of_pos (sub_pos.mpr (hU j))]
      rw [h1, h2, ite_eq_left (Or.inr hU), zero_add, add_zero]
    · obtain ⟨j1, hj1⟩ := not_forall.mp hU
      have h2 : ∏ j, rUpperDev d (t j) (X i j) = 0 :=
        prod_eq_zero (mem_univ j1) (by rw [rUpperDev, ite_eq_right hj1])
      rw [h1, h2, ite_eq_right (not_or.mpr ⟨hL, hU⟩), zero_add, zero_add]

theorem totalND_eq_absND (hk : 0 < k) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ)
    (d : ℝ) : totalND X t d = absND X t d := by
  unfold totalND clpmND cupmND dpmND absND
  rw [← add_div, ← add_div, ← sum_add_distrib, ← sum_add_distrib]
  congr 1
  exact sum_congr rfl fun i _ => lowRow_add_upRow_add_dpmRow hk X t d i

theorem dpmND_eq_abs_sub (hk : 0 < k) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) :
    dpmND X t d = absND X t d - clpmND X t d - cupmND X t d := by
  rw [← totalND_eq_absND hk, totalND]
  ring

theorem absND_zero (hn : 0 < n) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) :
    absND X t 0 = 1 := by
  unfold absND absRow
  simp only [Real.rpow_zero, prod_const_one, sum_const, card_univ, Fintype.card_fin,
    nsmul_eq_mul, mul_one]
  exact div_self (Nat.cast_ne_zero.mpr hn.ne')

theorem degree_zero_partition_v8898017 (hn : 0 < n) (hk : 0 < k) (X : Matrix (Fin n) (Fin k) ℝ)
    (t : Fin k → ℝ) (norm : Bool) :
    coLPMnD_R_v8898017 X t 0 norm + coUPMnD_R_v8898017 X t 0 norm
      + dpmND_R_v8898017 X t 0 norm = 1 := by
  simp only [coLPMnD_R_eq_v8898017 le_rfl, coUPMnD_R_eq_v8898017 le_rfl,
    dpmND_R_eq_v8898017 le_rfl, ne_eq, not_true_eq_false, false_and, ite_false]
  rw [← totalND, totalND_eq_absND hk, absND_zero hn]

open scoped Classical in
theorem coLPMnD_R_zero_v8898017 (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (norm : Bool) :
    coLPMnD_R_v8898017 X t 0 norm = (#{i | ∀ j, X i j ≤ t j} : ℝ) / n := by
  rw [coLPMnD_R_eq_v8898017 le_rfl]
  simp only [ne_eq, not_true_eq_false, false_and, ite_false, clpmND, lowRow, rLowerDev_zero]
  simp only [prod_ite_zero, prod_const_one, sum_boole, mem_univ, forall_const]

open scoped Classical in
theorem coUPMnD_R_zero_v8898017 (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (norm : Bool) :
    coUPMnD_R_v8898017 X t 0 norm = (#{i | ∀ j, t j < X i j} : ℝ) / n := by
  rw [coUPMnD_R_eq_v8898017 le_rfl]
  simp only [ne_eq, not_true_eq_false, false_and, ite_false, cupmND, upRow, rUpperDev_zero]
  simp only [prod_ite_zero, prod_const_one, sum_boole, mem_univ, forall_const]

theorem zero_cols (hn : 0 < n) (X : Matrix (Fin n) (Fin 0) ℝ) (t : Fin 0 → ℝ) (d : ℝ) :
    clpmND X t d = 1 ∧ cupmND X t d = 1 ∧ dpmND X t d = 0 ∧ absND X t d = 1 := by
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  simp [clpmND, cupmND, dpmND, absND, lowRow, upRow, dpmRow, absRow, hn']

structure FusedAcc where
  lower : ℝ
  upper : ℝ
  dpm : ℝ
  belowStrict : Bool
  aboveStrict : Bool

def fusedStep (d diff : ℝ) (s : FusedAcc) : FusedAcc where
  lower := s.lower * lowerComp (-diff) d
  upper := s.upper * upperComp diff d
  dpm := s.dpm * cppPow |diff| d
  belowStrict := s.belowStrict && !decide (diff ≥ 0)
  aboveStrict := s.aboveStrict && !decide (diff ≤ 0)

def fusedRow (d : ℝ) (x tr : Fin k → ℝ) : FusedAcc :=
  (List.finRange k).foldl (fun s j => fusedStep d (x j - tr j) s) ⟨1, 1, 1, true, true⟩

structure BatchSums where
  clpm : ℝ
  cupm : ℝ
  dpm : ℝ

def batchStep (d : ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (tr : Fin k → ℝ) (s : BatchSums)
    (i : Fin n) : BatchSums :=
  let f := fusedRow d (X i) tr
  { clpm := s.clpm + f.lower, cupm := s.cupm + f.upper,
    dpm := if !(f.belowStrict || f.aboveStrict) then s.dpm + f.dpm else s.dpm }

def batchCount_v8898017 (X : Matrix (Fin n) (Fin k) ℝ) (tr : Fin k → ℝ) : ℝ :=
  (List.finRange n).foldl
    (fun c i => if allLoop (fun j => ¬ X i j > tr j) (List.finRange k) then c + 1 else c) 0

def batchAt_v8898017 (X : Matrix (Fin n) (Fin k) ℝ) (tr : Fin k → ℝ) (d : ℝ) (norm : Bool) :
    ℝ :=
  if d = 0 then batchCount_v8898017 X tr / n
  else
    let s := (List.finRange n).foldl (batchStep d X tr) ⟨0, 0, 0⟩
    let invN := 1 / (n : ℝ)
    let clpmUn := s.clpm * invN
    if !norm then clpmUn else normQ (clpmUn + s.cupm * invN + s.dpm * invN) clpmUn

def coLPMnDBatch_v8898017 (X : Matrix (Fin n) (Fin k) ℝ) (T : Matrix (Fin m) (Fin k) ℝ)
    (d : ℝ) (norm : Bool) : Fin m → ℝ :=
  fun r => batchAt_v8898017 X (T r) d norm

lemma not_decide_nonneg (x : ℝ) : (!decide (x ≥ 0)) = decide (x < 0) := by
  by_cases h : 0 ≤ x
  · simp [h, not_lt.mpr h]
  · simp [h, not_le.mp h]

lemma not_decide_nonpos (x : ℝ) : (!decide (x ≤ 0)) = decide (0 < x) := by
  by_cases h : x ≤ 0
  · simp [h, not_lt.mpr h]
  · simp [h, not_le.mp h]

lemma foldl_fusedStep {α : Type*} (d : ℝ) (g : α → ℝ) (l : List α) (s : FusedAcc) :
    l.foldl (fun s a => fusedStep d (g a) s) s =
      { lower := s.lower * (l.map fun a => lowerComp (-g a) d).prod
        upper := s.upper * (l.map fun a => upperComp (g a) d).prod
        dpm := s.dpm * (l.map fun a => cppPow |g a| d).prod
        belowStrict := s.belowStrict && decide (∀ a ∈ l, ¬ g a ≥ 0)
        aboveStrict := s.aboveStrict && decide (∀ a ∈ l, ¬ g a ≤ 0) } := by
  induction l generalizing s with
  | nil => simp
  | cons a l ih =>
    rw [List.foldl_cons, ih]
    simp only [fusedStep, not_decide_nonneg, not_decide_nonpos]
    simp [mul_assoc, Bool.and_assoc]

lemma fusedRow_lower {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ) (tr : Fin k → ℝ)
    (i : Fin n) : (fusedRow d (X i) tr).lower = lowRow X tr d i := by
  rw [fusedRow, foldl_fusedStep]
  dsimp only
  rw [one_mul, ← Fin.prod_univ_def, lowRow]
  exact prod_congr rfl fun j _ => by rw [neg_sub, lowerComp_eq hd]

lemma fusedRow_upper {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ) (tr : Fin k → ℝ)
    (i : Fin n) : (fusedRow d (X i) tr).upper = upRow X tr d i := by
  rw [fusedRow, foldl_fusedStep]
  dsimp only
  rw [one_mul, ← Fin.prod_univ_def, upRow]
  exact prod_congr rfl fun j _ => upperComp_eq hd _ _

lemma fusedRow_dpm (d : ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (tr : Fin k → ℝ) (i : Fin n) :
    (if ((fusedRow d (X i) tr).belowStrict || (fusedRow d (X i) tr).aboveStrict) then 0
      else (fusedRow d (X i) tr).dpm) = dpmSumRow X tr d i := by
  rw [fusedRow, foldl_fusedStep, dpmSumRow, flagLoop_eq, foldl_mul_eq]

lemma foldl_batchStep (d : ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (tr : Fin k → ℝ)
    (l : List (Fin n)) (s : BatchSums) :
    l.foldl (batchStep d X tr) s =
      { clpm := s.clpm + (l.map fun i => (fusedRow d (X i) tr).lower).sum
        cupm := s.cupm + (l.map fun i => (fusedRow d (X i) tr).upper).sum
        dpm := s.dpm + (l.map fun i =>
          if ((fusedRow d (X i) tr).belowStrict || (fusedRow d (X i) tr).aboveStrict) then 0
          else (fusedRow d (X i) tr).dpm).sum } := by
  induction l generalizing s with
  | nil => simp
  | cons a l ih =>
    rw [List.foldl_cons, ih]
    simp only [batchStep, List.map_cons, List.sum_cons, BatchSums.mk.injEq]
    refine ⟨by ring, by ring, ?_⟩
    split_ifs <;> simp_all; ring

lemma batchCount_eq_v8898017 (X : Matrix (Fin n) (Fin k) ℝ) (tr : Fin k → ℝ) :
    batchCount_v8898017 X tr = sugarSum (clpmCountRow X tr) := by
  unfold batchCount_v8898017 sugarSum clpmCountRow
  congr 1
  funext c i
  split_ifs <;> simp

theorem coLPMnDBatch_eq_v8898017 {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ)
    (T : Matrix (Fin m) (Fin k) ℝ) (norm : Bool) (r : Fin m) :
    coLPMnDBatch_v8898017 X T d norm r = coLPMnD_R_v8898017 X (T r) d norm := by
  unfold coLPMnDBatch_v8898017 batchAt_v8898017
  by_cases hd0 : d = 0
  · rw [ite_eq_left hd0, coLPMnD_R_v8898017, ite_eq_left hd0, clpmRaw, ite_eq_left hd0,
      batchCount_eq_v8898017]
  · rw [ite_eq_right hd0, coLPMnD_R_eq_v8898017 hd]
    simp only [foldl_batchStep, zero_add, ← Fin.sum_univ_def, fusedRow_lower hd,
      fusedRow_upper hd, fusedRow_dpm, dpmSumRow_eq hd hd0, mul_one_div]
    cases norm <;> simp [hd0, clpmND, cupmND, dpmND, totalND]

lemma rUpperDev_eq_rLowerDev_neg {d : ℝ} (hd : d ≠ 0) (t v : ℝ) :
    rUpperDev d t v = rLowerDev d (-t) (-v) := by
  unfold rUpperDev rLowerDev
  rcases lt_trichotomy v t with h | h | h
  · rw [ite_eq_right (not_lt.mpr h.le), ite_eq_right (by linarith)]
  · subst h
    simp [Real.zero_rpow hd]
  · rw [ite_eq_left h, ite_eq_left (by linarith), neg_sub_neg]

theorem clpmND_neg {d : ℝ} (hd : d ≠ 0) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) :
    clpmND (-X) (-t) d = cupmND X t d := by
  unfold clpmND cupmND lowRow upRow
  simp only [Matrix.neg_apply, Pi.neg_apply, rUpperDev_eq_rLowerDev_neg hd]

theorem cupmND_neg {d : ℝ} (hd : d ≠ 0) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) :
    cupmND (-X) (-t) d = clpmND X t d := by
  unfold clpmND cupmND lowRow upRow
  simp only [Matrix.neg_apply, Pi.neg_apply, rUpperDev_eq_rLowerDev_neg hd, neg_neg]

theorem absND_neg (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) :
    absND (-X) (-t) d = absND X t d := by
  unfold absND absRow
  simp only [Matrix.neg_apply, Pi.neg_apply, neg_sub_neg, abs_sub_comm]

theorem dpmND_neg {d : ℝ} (hd : d ≠ 0) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) :
    dpmND (-X) (-t) d = dpmND X t d := by
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · simp [dpmND, dpmRow]
  · rw [dpmND_eq_abs_sub hk, dpmND_eq_abs_sub hk, absND_neg, clpmND_neg hd, cupmND_neg hd]
    ring

theorem totalND_neg {d : ℝ} (hd : d ≠ 0) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) :
    totalND (-X) (-t) d = totalND X t d := by
  rw [totalND, totalND, clpmND_neg hd, cupmND_neg hd, dpmND_neg hd]
  ring

theorem coLPMnDBatch_neg_eq_coUPMnD_R_v8898017 {d : ℝ} (hd : 0 < d)
    (X : Matrix (Fin n) (Fin k) ℝ) (T : Matrix (Fin m) (Fin k) ℝ) (norm : Bool) (r : Fin m) :
    coLPMnDBatch_v8898017 (-X) (-T) d norm r = coUPMnD_R_v8898017 X (T r) d norm := by
  rw [coLPMnDBatch_eq_v8898017 hd.le, coLPMnD_R_eq_v8898017 hd.le, coUPMnD_R_eq_v8898017 hd.le]
  have hT : (-T) r = -(T r) := rfl
  rw [hT, totalND_neg hd.ne', clpmND_neg hd.ne']

def coLPMnD_R (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (norm : Bool) : ℝ :=
  if d ≠ 0 ∧ norm = true then normQ0 (normConstRaw X t d) (clpmRaw X t d) (clpmRaw X t 0)
  else clpmRaw X t d

def coUPMnD_R (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (norm : Bool) : ℝ :=
  if d ≠ 0 ∧ norm = true then normQ0 (normConstRaw X t d) (cupmRaw X t d) (cupmRaw X t 0)
  else cupmRaw X t d

def dpmND_R (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (norm : Bool) : ℝ :=
  if d ≠ 0 ∧ norm = true then normQ0 (normConstRaw X t d) (dpmRaw X t d) (dpmRaw X t 0)
  else dpmRaw X t d

def cellND (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) : Fin 3 → ℝ :=
  ![clpmND X t d, cupmND X t d, dpmND X t d]

def normND (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (c : Fin 3) : ℝ :=
  if 0 < totalND X t d then cellND X t d c / totalND X t d else cellND X t 0 c

lemma cellND_sum (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) :
    cellND X t d 0 + cellND X t d 1 + cellND X t d 2 = totalND X t d := rfl

lemma cellND_nonneg (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (c : Fin 3) :
    0 ≤ cellND X t d c := by
  fin_cases c
  · exact clpmND_nonneg X t d
  · exact cupmND_nonneg X t d
  · exact dpmND_nonneg X t d

private theorem normQ0_total (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (c : Fin 3) :
    normQ0 (totalND X t d) (cellND X t d c) (cellND X t 0 c) = normND X t d c := by
  have h : 0 ≤ totalND X t d := add_nonneg (add_nonneg (clpmND_nonneg ..) (cupmND_nonneg ..))
    (dpmND_nonneg ..)
  rcases h.eq_or_lt with he | hp
  · simp [normQ0, normND, ← he]
  · simp [normQ0, normND, hp]

theorem coLPMnD_R_eq {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ)
    (norm : Bool) : coLPMnD_R X t d norm =
      if d ≠ 0 ∧ norm = true then normND X t d 0 else clpmND X t d := by
  simp only [coLPMnD_R, normConstRaw_eq hd, clpmRaw_eq hd, clpmRaw_eq (le_refl 0)]
  split_ifs
  · exact normQ0_total X t d 0
  · rfl

theorem coUPMnD_R_eq {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ)
    (norm : Bool) : coUPMnD_R X t d norm =
      if d ≠ 0 ∧ norm = true then normND X t d 1 else cupmND X t d := by
  simp only [coUPMnD_R, normConstRaw_eq hd, cupmRaw_eq hd, cupmRaw_eq (le_refl 0)]
  split_ifs
  · exact normQ0_total X t d 1
  · rfl

private theorem dpmND_R_closed {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ)
    (norm : Bool) : dpmND_R X t d norm =
      if d ≠ 0 ∧ norm = true then normND X t d 2 else dpmND X t d := by
  simp only [dpmND_R, normConstRaw_eq hd, dpmRaw_eq hd, dpmRaw_eq (le_refl 0)]
  split_ifs
  · exact normQ0_total X t d 2
  · rfl

theorem dpmND_R_eq {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ)
    (norm : Bool) : dpmND_R X t d norm =
      if d ≠ 0 ∧ norm = true then
        if 0 < totalND X t d then dpmND X t d / totalND X t d else dpmND X t 0
      else dpmND X t d :=
  dpmND_R_closed hd X t norm

theorem coLPMnD_R_eq_v8898017_of {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ)
    (t : Fin k → ℝ) (norm : Bool) (h : d = 0 ∨ norm = false ∨ 0 < totalND X t d) :
    coLPMnD_R X t d norm = coLPMnD_R_v8898017 X t d norm := by
  rw [coLPMnD_R_eq hd, coLPMnD_R_eq_v8898017 hd]
  rcases h with h | h | h
  · simp [h]
  · simp [h]
  · simp [normND, normQ, h, cellND]

theorem coUPMnD_R_eq_v8898017_of {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ)
    (t : Fin k → ℝ) (norm : Bool) (h : d = 0 ∨ norm = false ∨ 0 < totalND X t d) :
    coUPMnD_R X t d norm = coUPMnD_R_v8898017 X t d norm := by
  rw [coUPMnD_R_eq hd, coUPMnD_R_eq_v8898017 hd]
  rcases h with h | h | h
  · simp [h]
  · simp [h]
  · simp [normND, normQ, h, cellND]

theorem dpmND_R_eq_v8898017_of {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ)
    (t : Fin k → ℝ) (norm : Bool) (h : d = 0 ∨ norm = false ∨ 0 < totalND X t d) :
    dpmND_R X t d norm = dpmND_R_v8898017 X t d norm := by
  rw [dpmND_R_closed hd, dpmND_R_eq_v8898017 hd]
  rcases h with h | h | h
  · simp [h]
  · simp [h]
  · simp [normND, normQ, h, cellND]

lemma normND_nonneg (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (c : Fin 3) :
    0 ≤ normND X t d c := by
  unfold normND
  split_ifs with h
  · exact div_nonneg (cellND_nonneg X t d c) h.le
  · exact cellND_nonneg X t 0 c

theorem normND_sum (hn : 0 < n) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) :
    normND X t d 0 + normND X t d 1 + normND X t d 2 = 1 := by
  unfold normND
  by_cases htot : 0 < totalND X t d
  · simp only [htot, ite_true]
    rw [← add_div, ← add_div, cellND_sum, div_self htot.ne']
  · simp only [htot, ite_false]
    rw [cellND_sum]
    rcases Nat.eq_zero_or_pos k with rfl | hk
    · obtain ⟨h1, h2, h3, -⟩ := zero_cols hn X t d
      exact absurd (by rw [totalND, h1, h2, h3]; norm_num) htot
    · rw [totalND_eq_absND hk, absND_zero hn]

theorem degree_zero_partition (hn : 0 < n) (hk : 0 < k) (X : Matrix (Fin n) (Fin k) ℝ)
    (t : Fin k → ℝ) (norm : Bool) :
    coLPMnD_R X t 0 norm + coUPMnD_R X t 0 norm + dpmND_R X t 0 norm = 1 := by
  rw [coLPMnD_R_eq_v8898017_of le_rfl X t norm (Or.inl rfl),
    coUPMnD_R_eq_v8898017_of le_rfl X t norm (Or.inl rfl),
    dpmND_R_eq_v8898017_of le_rfl X t norm (Or.inl rfl)]
  exact degree_zero_partition_v8898017 hn hk X t norm

lemma normND_mem_Icc (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (d : ℝ) (c : Fin 3) :
    normND X t d c ∈ Set.Icc 0 1 := by
  refine ⟨normND_nonneg X t d c, ?_⟩
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · fin_cases c <;> simp [normND, cellND, clpmND, cupmND, dpmND, totalND]
  · have hs := normND_sum hn X t d
    have h0 := normND_nonneg X t d 0
    have h1 := normND_nonneg X t d 1
    have h2 := normND_nonneg X t d 2
    fin_cases c
    · change normND X t d 0 ≤ 1
      linarith
    · change normND X t d 1 ≤ 1
      linarith
    · change normND X t d 2 ≤ 1
      linarith

open scoped Classical in
theorem coLPMnD_R_zero (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (norm : Bool) :
    coLPMnD_R X t 0 norm = (#{i | ∀ j, X i j ≤ t j} : ℝ) / n := by
  rw [coLPMnD_R_eq_v8898017_of le_rfl X t norm (Or.inl rfl), coLPMnD_R_zero_v8898017]

open scoped Classical in
theorem coUPMnD_R_zero (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) (norm : Bool) :
    coUPMnD_R X t 0 norm = (#{i | ∀ j, t j < X i j} : ℝ) / n := by
  rw [coUPMnD_R_eq_v8898017_of le_rfl X t norm (Or.inl rfl), coUPMnD_R_zero_v8898017]

def batchCount (strict : Bool) (X : Matrix (Fin n) (Fin k) ℝ) (tr : Fin k → ℝ) : ℝ :=
  sugarSum fun i => if allLoop (fun j => if strict then X i j < tr j else X i j ≤ tr j)
    (List.finRange k) then 1 else 0

def batchAt (X : Matrix (Fin n) (Fin k) ℝ) (tr : Fin k → ℝ) (d : ℝ) (norm strict : Bool) : ℝ :=
  if d = 0 then batchCount strict X tr / n
  else
    let s := (List.finRange n).foldl (batchStep d X tr) ⟨0, 0, 0⟩
    let invN := 1 / (n : ℝ)
    let lower := s.clpm * invN
    if !norm then lower
    else normQ0 (lower + s.cupm * invN + s.dpm * invN) lower (batchCount strict X tr / n)

def coLPMnDBatch (X : Matrix (Fin n) (Fin k) ℝ) (T : Matrix (Fin m) (Fin k) ℝ) (d : ℝ)
    (norm strict : Bool) : Fin m → ℝ :=
  fun r => batchAt X (T r) d norm strict

private theorem batchCount_neg (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) :
    batchCount true (-X) (-t) / n = cupmRaw X t 0 := by
  simp only [cupmRaw, batchCount]
  congr 2
  funext i
  simp [cupmCountRow, not_le]

private theorem batchAt_pos {d : ℝ} (hd : 0 < d) (X : Matrix (Fin n) (Fin k) ℝ)
    (t : Fin k → ℝ) (norm strict : Bool) :
    batchAt X t d norm strict = if norm then
      normQ0 (totalND X t d) (clpmND X t d) (batchCount strict X t / n)
      else clpmND X t d := by
  simp only [batchAt, hd.ne', ite_false, foldl_batchStep, zero_add, ← Fin.sum_univ_def,
    fusedRow_lower hd.le, fusedRow_upper hd.le, fusedRow_dpm, dpmSumRow_eq hd.le hd.ne',
    mul_one_div]
  cases norm <;> rfl

theorem coLPMnDBatch_neg_eq_coUPMnD_R {d : ℝ} (hd : 0 ≤ d) (X : Matrix (Fin n) (Fin k) ℝ)
    (T : Matrix (Fin m) (Fin k) ℝ) (norm : Bool) (r : Fin m) :
    coLPMnDBatch (-X) (-T) d norm true r = coUPMnD_R X (T r) d norm := by
  change batchAt (-X) (-(T r)) d norm true = _
  rcases hd.eq_or_lt with h | hp
  · subst d
    simp [batchAt, coUPMnD_R, batchCount_neg]
  · simp [batchAt_pos hp, coUPMnD_R, hp.ne', batchCount_neg, normConstRaw_eq hp.le,
      cupmRaw_eq hp.le, clpmND_neg hp.ne', totalND_neg hp.ne']

end

end NNS
