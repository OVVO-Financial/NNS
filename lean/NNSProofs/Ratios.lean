-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.PartialMoments
public import Mathlib.Topology.MetricSpace.Pseudo.Defs

namespace NNS

open Finset Filter Topology

@[expose] public noncomputable section

variable {n : ℕ}

def LPMratio (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : ℝ :=
  if d = 0 then LPM 0 t x
  else if 0 < n ∧ ∀ i, x i = t then 1
  else LPM d t x / (LPM d t x + UPM d t x)

def UPMratio (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : ℝ :=
  if d = 0 then UPM 0 t x
  else if 0 < n ∧ ∀ i, x i = t then 0
  else UPM d t x / (LPM d t x + UPM d t x)

private lemma lowerDev_mono (d : ℕ) {s t : ℝ} (hst : s ≤ t) (v : ℝ) :
    lowerDev d s v ≤ lowerDev d t v := by
  unfold lowerDev
  split_ifs with h1 h2
  · exact pow_le_pow_left₀ (sub_nonneg.mpr h1) (by linarith) d
  · exact absurd (le_trans h1 hst) h2
  · exact pow_nonneg (sub_nonneg.mpr ‹v ≤ t›) d
  · exact le_refl 0

private lemma upperDev_anti (d : ℕ) {s t : ℝ} (hst : s ≤ t) (v : ℝ) :
    upperDev d t v ≤ upperDev d s v := by
  unfold upperDev
  split_ifs with h1 h2
  · exact pow_le_pow_left₀ (sub_nonneg.mpr h1.le) (by linarith) d
  · exact absurd (lt_of_le_of_lt hst h1) h2
  · exact pow_nonneg (sub_nonneg.mpr ‹s < v›.le) d
  · exact le_refl 0

private lemma LPM_mono (d : ℕ) (x : Fin n → ℝ) {s t : ℝ} (hst : s ≤ t) :
    LPM d s x ≤ LPM d t x := by
  unfold LPM
  gcongr with i
  exact lowerDev_mono d hst (x i)

private lemma UPM_anti (d : ℕ) (x : Fin n → ℝ) {s t : ℝ} (hst : s ≤ t) :
    UPM d t x ≤ UPM d s x := by
  unfold UPM
  gcongr with i
  exact upperDev_anti d hst (x i)

private lemma lowerDev_add_upperDev_pos {d : ℕ} {t v : ℝ} (hv : v ≠ t) :
    0 < lowerDev d t v + upperDev d t v := by
  unfold lowerDev upperDev
  rcases lt_or_gt_of_ne hv with h | h
  · have h' : ¬ t < v := not_lt.mpr h.le
    simp only [h.le, h', ↓reduceIte, add_zero]
    exact pow_pos (sub_pos.mpr h) d
  · have h' : ¬ v ≤ t := not_le.mpr h
    simp only [h', h, ↓reduceIte, zero_add]
    exact pow_pos (sub_pos.mpr h) d

private lemma lowerDev_self {d : ℕ} (hd : d ≠ 0) (t : ℝ) : lowerDev d t t = 0 := by
  simp [lowerDev, hd]

private lemma upperDev_self (d : ℕ) (t : ℝ) : upperDev d t t = 0 := by
  simp [upperDev]

theorem ratio_denom_pos {d : ℕ} {t : ℝ} {x : Fin n → ℝ} (h : ∃ i, x i ≠ t) :
    0 < LPM d t x + UPM d t x := by
  obtain ⟨i, hi⟩ := h
  have hn : (0 : ℝ) < n := Nat.cast_pos.mpr (Fin.pos i)
  unfold LPM UPM
  rw [← add_div, ← sum_add_distrib]
  refine div_pos ?_ hn
  refine sum_pos' (fun j _ => add_nonneg (lowerDev_nonneg _ _ _) (upperDev_nonneg _ _ _))
    ⟨i, mem_univ _, lowerDev_add_upperDev_pos hi⟩

theorem ratio_denom_eq_zero_iff {d : ℕ} (hd : d ≠ 0) (t : ℝ) (x : Fin n → ℝ) :
    LPM d t x + UPM d t x = 0 ↔ ∀ i, x i = t := by
  constructor
  · intro h
    by_contra hc
    push Not at hc
    exact (ratio_denom_pos hc).ne' h
  · intro h
    unfold LPM UPM
    rw [← add_div, ← sum_add_distrib]
    simp [h, lowerDev_self hd, upperDev_self]

theorem LPMratio_of_n_eq_zero (d : ℕ) (t : ℝ) (x : Fin 0 → ℝ) : LPMratio d t x = 0 := by
  simp [LPMratio, LPM, UPM]

theorem UPMratio_of_n_eq_zero (d : ℕ) (t : ℝ) (x : Fin 0 → ℝ) : UPMratio d t x = 0 := by
  simp [UPMratio, LPM, UPM]

theorem LPMratio_eq_sum_div {d : ℕ} (hd : d ≠ 0) {t : ℝ} {x : Fin n → ℝ}
    (h : ∃ i, x i ≠ t) :
    LPMratio d t x =
      (∑ i, lowerDev d t (x i)) / ((∑ i, lowerDev d t (x i)) + ∑ i, upperDev d t (x i)) := by
  obtain ⟨i, hi⟩ := h
  have hn : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Fin.pos i).ne'
  have hdeg : ¬ (0 < n ∧ ∀ j, x j = t) := fun h' => hi (h'.2 i)
  simp only [LPMratio, hd, hdeg, ↓reduceIte, LPM, UPM]
  rw [← add_div, div_div_div_cancel_right₀ hn]

theorem LPMratio_nonneg (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : 0 ≤ LPMratio d t x := by
  unfold LPMratio
  split_ifs
  · exact LPM_nonneg 0 t x
  · exact zero_le_one
  · exact div_nonneg (LPM_nonneg d t x) (add_nonneg (LPM_nonneg d t x) (UPM_nonneg d t x))

private lemma LPM_zero_le_one (t : ℝ) (x : Fin n → ℝ) : LPM 0 t x ≤ 1 := by
  rw [LPM_zero_eq_ecdf]
  refine div_le_one_of_le₀ ?_ n.cast_nonneg
  have := card_filter_le (univ : Finset (Fin n)) (fun i => x i ≤ t)
  simp only [card_univ, Fintype.card_fin] at this
  exact_mod_cast this

theorem LPMratio_le_one (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : LPMratio d t x ≤ 1 := by
  unfold LPMratio
  split_ifs
  · exact LPM_zero_le_one t x
  · exact le_refl 1
  · exact div_le_one_of_le₀ (le_add_of_nonneg_right (UPM_nonneg d t x))
      (add_nonneg (LPM_nonneg d t x) (UPM_nonneg d t x))

theorem LPMratio_mem_Icc (d : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    LPMratio d t x ∈ Set.Icc (0 : ℝ) 1 :=
  ⟨LPMratio_nonneg d t x, LPMratio_le_one d t x⟩

theorem LPMratio_add_UPMratio (hn : 0 < n) (d : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    LPMratio d t x + UPMratio d t x = 1 := by
  unfold LPMratio UPMratio
  split_ifs with hd hdeg
  · exact LPM_zero_add_UPM_zero hn t x
  · exact add_zero 1
  · have hex : ∃ i, x i ≠ t := by
      by_contra hc
      push Not at hc
      exact hdeg ⟨hn, hc⟩
    rw [← add_div, div_self (ratio_denom_pos hex).ne']

theorem UPMratio_eq_one_sub_LPMratio (hn : 0 < n) (d : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    UPMratio d t x = 1 - LPMratio d t x := by
  rw [← LPMratio_add_UPMratio hn d t x]
  ring

theorem LPMratio_zero_eq_ecdf (t : ℝ) (x : Fin n → ℝ) :
    LPMratio 0 t x = (#{i | x i ≤ t} : ℝ) / n := by
  simp only [LPMratio, ↓reduceIte]
  exact LPM_zero_eq_ecdf t x

theorem UPMratio_zero_eq_survival (t : ℝ) (x : Fin n → ℝ) :
    UPMratio 0 t x = (#{i | t < x i} : ℝ) / n := by
  simp only [UPMratio, ↓reduceIte, UPM, upperDev, pow_zero]
  rw [sum_boole]

theorem LPMratio_of_forall_eq {d : ℕ} (hd : d ≠ 0) (hn : 0 < n) {t : ℝ} {x : Fin n → ℝ}
    (h : ∀ i, x i = t) : LPMratio d t x = 1 := by
  simp [LPMratio, hd, hn, h]

theorem LPMratio_eq_one_of_forall_le (hn : 0 < n) (d : ℕ) {t : ℝ} {x : Fin n → ℝ}
    (h : ∀ i, x i ≤ t) : LPMratio d t x = 1 := by
  by_cases hd : d = 0
  · subst hd
    rw [LPMratio_zero_eq_ecdf]
    have hfilt : ({i | x i ≤ t} : Finset (Fin n)) = univ := by
      ext i; simp [h i]
    rw [hfilt, card_univ, Fintype.card_fin]
    exact div_self (Nat.cast_ne_zero.mpr hn.ne')
  by_cases hall : ∀ i, x i = t
  · exact LPMratio_of_forall_eq hd hn hall
  push Not at hall
  have hU : UPM d t x = 0 := by
    unfold UPM upperDev
    simp [fun i => not_lt.mpr (h i)]
  have hpos := ratio_denom_pos (d := d) hall
  have hdeg : ¬ (0 < n ∧ ∀ j, x j = t) := fun h' => hall.elim fun i hi => hi (h'.2 i)
  rw [hU, add_zero] at hpos
  simp only [LPMratio, hd, hdeg, ↓reduceIte, hU, add_zero]
  exact div_self hpos.ne'

theorem LPMratio_eq_zero_of_forall_lt (d : ℕ) {t : ℝ} {x : Fin n → ℝ}
    (h : ∀ i, t < x i) : LPMratio d t x = 0 := by
  have hL : ∀ e, LPM e t x = 0 := by
    intro e
    unfold LPM lowerDev
    simp [fun i => not_le.mpr (h i)]
  unfold LPMratio
  split_ifs with hd hdeg
  · exact hL 0
  · obtain ⟨hn, hall⟩ := hdeg
    have := h ⟨0, hn⟩
    rw [hall ⟨0, hn⟩] at this
    exact absurd this (lt_irrefl t)
  · rw [hL d, zero_div]

theorem LPMratio_eq_zero_of_forall_ge {d : ℕ} (hd : d ≠ 0) {t : ℝ} {x : Fin n → ℝ}
    (h : ∀ i, t ≤ x i) (hlt : ∃ i, t < x i) : LPMratio d t x = 0 := by
  have hex : ∃ i, x i ≠ t := hlt.imp fun _ hi => hi.ne'
  rw [LPMratio_eq_sum_div hd hex]
  have : ∀ i, lowerDev d t (x i) = 0 := by
    intro i
    unfold lowerDev
    split_ifs with hi
    · rw [le_antisymm hi (h i), sub_self, zero_pow hd]
    · rfl
  simp [this]

theorem UPMratio_eq_zero_of_forall_le (d : ℕ) {t : ℝ} {x : Fin n → ℝ}
    (h : ∀ i, x i ≤ t) : UPMratio d t x = 0 := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · exact UPMratio_of_n_eq_zero d t x
  rw [UPMratio_eq_one_sub_LPMratio hn, LPMratio_eq_one_of_forall_le hn d h, sub_self]

theorem UPMratio_eq_one_of_forall_lt (hn : 0 < n) (d : ℕ) {t : ℝ} {x : Fin n → ℝ}
    (h : ∀ i, t < x i) : UPMratio d t x = 1 := by
  rw [UPMratio_eq_one_sub_LPMratio hn, LPMratio_eq_zero_of_forall_lt d h, sub_zero]

theorem LPMratio_monotone (d : ℕ) (x : Fin n → ℝ) : Monotone fun t => LPMratio d t x := by
  intro s t hst
  dsimp only
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [LPMratio_of_n_eq_zero]
  by_cases hd : d = 0
  · subst hd
    simp only [LPMratio, ↓reduceIte]
    exact LPM_zero_monotone x hst
  by_cases hs : ∀ i, x i = s
  · rw [LPMratio_of_forall_eq hd hn hs,
      LPMratio_eq_one_of_forall_le hn d fun i => (hs i).le.trans hst]
  by_cases ht : ∀ i, x i = t
  · rw [LPMratio_of_forall_eq hd hn ht]
    exact LPMratio_le_one d s x
  push Not at hs ht
  have hps := ratio_denom_pos (d := d) hs
  have hpt := ratio_denom_pos (d := d) ht
  have hdegs : ¬ (0 < n ∧ ∀ j, x j = s) := fun h' => hs.elim fun i hi => hi (h'.2 i)
  have hdegt : ¬ (0 < n ∧ ∀ j, x j = t) := fun h' => ht.elim fun i hi => hi (h'.2 i)
  simp only [LPMratio, hd, hdegs, hdegt, ↓reduceIte]
  rw [div_le_div_iff₀ hps hpt]
  have hL := LPM_mono d x hst
  have hU := UPM_anti d x hst
  have hLs := LPM_nonneg d s x
  have hUt := UPM_nonneg d t x
  nlinarith [mul_le_mul hL hU hUt (LPM_nonneg d t x)]

theorem LPMratio_lt_LPMratio {d : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ} {s t : ℝ}
    (hs : ∃ i, x i ≤ s) (ht : ∃ j, t ≤ x j) (hst : s < t) :
    LPMratio d s x < LPMratio d t x := by
  obtain ⟨i, hi⟩ := hs
  obtain ⟨j, hj⟩ := ht
  have hn : (0 : ℝ) < n := Nat.cast_pos.mpr (Fin.pos i)
  have hexs : ∃ k, x k ≠ s := ⟨j, (lt_of_lt_of_le hst hj).ne'⟩
  have hext : ∃ k, x k ≠ t := ⟨i, (lt_of_le_of_lt hi hst).ne⟩
  have hps := ratio_denom_pos (d := d) hexs
  have hpt := ratio_denom_pos (d := d) hext
  have hdegs : ¬ (0 < n ∧ ∀ k, x k = s) := fun h' => hexs.elim fun k hk => hk (h'.2 k)
  have hdegt : ¬ (0 < n ∧ ∀ k, x k = t) := fun h' => hext.elim fun k hk => hk (h'.2 k)
  have hL : LPM d s x < LPM d t x := by
    unfold LPM
    refine div_lt_div_of_pos_right ?_ hn
    refine sum_lt_sum (fun k _ => lowerDev_mono d hst.le (x k)) ⟨i, mem_univ _, ?_⟩
    simp only [lowerDev, hi, (hi.trans hst.le), ↓reduceIte]
    exact pow_lt_pow_left₀ (by linarith) (sub_nonneg.mpr hi) hd
  have hUs : 0 < UPM d s x := by
    unfold UPM
    refine div_pos (sum_pos' (fun k _ => upperDev_nonneg _ _ _) ⟨j, mem_univ _, ?_⟩) hn
    have hsj : s < x j := lt_of_lt_of_le hst hj
    simp only [upperDev, hsj, ↓reduceIte]
    exact pow_pos (sub_pos.mpr hsj) d
  have hU := UPM_anti d x hst.le
  have hLs := LPM_nonneg d s x
  simp only [LPMratio, hd, hdegs, hdegt, ↓reduceIte]
  rw [div_lt_div_iff₀ hps hpt]
  nlinarith [mul_le_mul_of_nonneg_left hU hLs, mul_lt_mul_of_pos_right hL hUs]

theorem LPMratio_strictMonoOn {d : ℕ} (hd : d ≠ 0) (x : Fin n → ℝ) :
    StrictMonoOn (fun t => LPMratio d t x) {t | (∃ i, x i ≤ t) ∧ ∃ j, t ≤ x j} :=
  fun _ hs _ ht hst => LPMratio_lt_LPMratio hd hs.1 ht.2 hst

theorem LPMratio_one_eq_half_iff {x : Fin n → ℝ} (h : ∃ i j, x i ≠ x j) (t : ℝ) :
    LPMratio 1 t x = 1 / 2 ↔ t = mean x := by
  obtain ⟨i, j, hij⟩ := h
  have hn : 0 < n := Fin.pos i
  have hex : ∃ k, x k ≠ t := by
    by_contra hc
    push Not at hc
    exact hij ((hc i).trans (hc j).symm)
  have hdeg : ¬ (0 < n ∧ ∀ k, x k = t) := fun h' => hex.elim fun k hk => hk (h'.2 k)
  have hpos := ratio_denom_pos (d := 1) hex
  have hdiff := UPM_one_sub_LPM_one hn t x
  simp only [LPMratio, one_ne_zero, hdeg, ↓reduceIte]
  rw [div_eq_div_iff hpos.ne' two_ne_zero]
  constructor
  · intro h; linarith
  · intro h; subst h; linarith

theorem LPMratio_one_mean {x : Fin n → ℝ} (h : ∃ i j, x i ≠ x j) :
    LPMratio 1 (mean x) x = 1 / 2 :=
  (LPMratio_one_eq_half_iff h (mean x)).mpr rfl

theorem UPMratio_one_mean {x : Fin n → ℝ} (h : ∃ i j, x i ≠ x j) :
    UPMratio 1 (mean x) x = 1 / 2 := by
  obtain ⟨i, _, _⟩ := h
  rw [UPMratio_eq_one_sub_LPMratio (Fin.pos i), LPMratio_one_mean ⟨i, _, ‹_›⟩]
  norm_num

theorem LPMratio_mean_of_const (hn : 0 < n) (d : ℕ) {x : Fin n → ℝ}
    (h : ∀ i j, x i = x j) : LPMratio d (mean x) x = 1 := by
  have hm : mean x = x ⟨0, hn⟩ := by
    unfold mean
    have : ∀ i, x i = x ⟨0, hn⟩ := fun i => h i _
    simp only [this, sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul]
    field_simp
  exact LPMratio_eq_one_of_forall_le hn d fun i => (h i _).le.trans hm.symm.le

private lemma lowerDev_affine (d : ℕ) {a : ℝ} (ha : 0 < a) (b t v : ℝ) :
    lowerDev d (a * t + b) (a * v + b) = a ^ d * lowerDev d t v := by
  unfold lowerDev
  have hiff : a * v + b ≤ a * t + b ↔ v ≤ t := by
    rw [add_le_add_iff_right, mul_le_mul_iff_right₀ ha]
  by_cases h : v ≤ t
  · simp only [hiff, h, ↓reduceIte]
    rw [← mul_pow]; ring
  · simp [hiff, h]

private lemma upperDev_affine (d : ℕ) {a : ℝ} (ha : 0 < a) (b t v : ℝ) :
    upperDev d (a * t + b) (a * v + b) = a ^ d * upperDev d t v := by
  unfold upperDev
  have hiff : a * t + b < a * v + b ↔ t < v := by
    rw [add_lt_add_iff_right, mul_lt_mul_iff_right₀ ha]
  by_cases h : t < v
  · simp only [hiff, h, ↓reduceIte]
    rw [← mul_pow]; ring
  · simp [hiff, h]

private lemma LPM_affine (d : ℕ) {a : ℝ} (ha : 0 < a) (b t : ℝ) (x : Fin n → ℝ) :
    LPM d (a * t + b) (fun i => a * x i + b) = a ^ d * LPM d t x := by
  simp only [LPM, lowerDev_affine d ha, ← mul_sum, mul_div_assoc]

private lemma UPM_affine (d : ℕ) {a : ℝ} (ha : 0 < a) (b t : ℝ) (x : Fin n → ℝ) :
    UPM d (a * t + b) (fun i => a * x i + b) = a ^ d * UPM d t x := by
  simp only [UPM, upperDev_affine d ha, ← mul_sum, mul_div_assoc]

theorem LPMratio_affine (d : ℕ) {a : ℝ} (ha : 0 < a) (b t : ℝ) (x : Fin n → ℝ) :
    LPMratio d (a * t + b) (fun i => a * x i + b) = LPMratio d t x := by
  have hiff : (∀ i, a * x i + b = a * t + b) ↔ ∀ i, x i = t := by
    refine forall_congr' fun i => ?_
    rw [add_left_inj, mul_right_inj' ha.ne']
  unfold LPMratio
  simp only [hiff, LPM_affine _ ha, UPM_affine _ ha, pow_zero, one_mul]
  split_ifs
  · rfl
  · rfl
  · rw [← mul_add, mul_div_mul_left _ _ (pow_ne_zero d ha.ne')]

end

end NNS
