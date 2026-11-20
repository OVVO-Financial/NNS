-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.Ratios
public import NNSProofs.CDF
public import Mathlib.Tactic.IntervalCases
public import Mathlib.Basic.Real.Sign

namespace NNS

open Finset

@[expose] public noncomputable section

variable {n : ℕ}

def sortedAt (x : Fin n → ℝ) (k : ℕ) : ℝ :=
  if h : k < n then x (Tuple.sort x ⟨k, h⟩) else 0

def VaRSorted (a : ℕ → ℝ) (m : ℕ) : Prop := ∀ i j, i ≤ j → j < m → a i ≤ a j

def countLe (a : ℕ → ℝ) (m : ℕ) (t : ℝ) : ℕ := #{r ∈ range m | a r ≤ t}

lemma sortedAt_sorted (x : Fin n → ℝ) : VaRSorted (sortedAt x) n := by
  intro i j hij hj
  have hi : i < n := lt_of_le_of_lt hij hj
  simp only [sortedAt, dite_eq_left hi, dite_eq_left hj]
  exact Tuple.monotone_sort x (show (⟨i, hi⟩ : Fin n) ≤ ⟨j, hj⟩ from hij)

lemma sum_sortedAt (x : Fin n → ℝ) (f : ℝ → ℝ) :
    ∑ k ∈ range n, f (sortedAt x k) = ∑ i, f (x i) := by
  rw [← Fin.sum_univ_eq_sum_range (fun k => f (sortedAt x k)),
    ← Equiv.sum_comp (Tuple.sort x) (fun i => f (x i))]
  refine sum_congr rfl fun i _ => ?_
  simp [sortedAt]

lemma exists_sortedAt_eq (x : Fin n → ℝ) (i : Fin n) : ∃ k < n, sortedAt x k = x i := by
  refine ⟨((Tuple.sort x).symm i).val, ((Tuple.sort x).symm i).isLt, ?_⟩
  simp [sortedAt]

lemma exists_eq_sortedAt (x : Fin n → ℝ) {k : ℕ} (hk : k < n) : ∃ i, x i = sortedAt x k :=
  ⟨Tuple.sort x ⟨k, hk⟩, by simp [sortedAt, hk]⟩

lemma lt_countLe_iff {a : ℕ → ℝ} {m : ℕ} (hs : VaRSorted a m) (t : ℝ) {r : ℕ} (hr : r < m) :
    r < countLe a m t ↔ a r ≤ t := by
  constructor
  · intro h
    by_contra hc
    push Not at hc
    have hsub : ({r' ∈ range m | a r' ≤ t} : Finset ℕ) ⊆ range r := by
      intro r' hr'
      rw [mem_filter, mem_range] at hr'
      rw [mem_range]
      by_contra hge
      push Not at hge
      exact absurd (hc.trans_le (hs _ _ hge hr'.1)) (not_lt.mpr hr'.2)
    have := card_le_card hsub
    rw [card_range] at this
    exact absurd h (not_lt.mpr this)
  · intro h
    have hsub : range (r + 1) ⊆ ({r' ∈ range m | a r' ≤ t} : Finset ℕ) := by
      intro r' hr'
      rw [mem_range] at hr'
      rw [mem_filter, mem_range]
      exact ⟨by omega, (hs _ _ (by omega) hr).trans h⟩
    have := card_le_card hsub
    rw [card_range] at this
    exact this

lemma countLe_le (a : ℕ → ℝ) (m : ℕ) (t : ℝ) : countLe a m t ≤ m := by
  unfold countLe
  exact (card_filter_le _ _).trans (card_range m).le

lemma countLe_sortedAt (x : Fin n → ℝ) (t : ℝ) :
    countLe (sortedAt x) n t = #{i | x i ≤ t} := by
  have h := sum_sortedAt x (fun v => if v ≤ t then (1 : ℝ) else 0)
  rw [sum_boole, sum_boole] at h
  unfold countLe
  exact_mod_cast h

def varMin (x : Fin n → ℝ) : ℝ := sortedAt x 0

def varMax (x : Fin n → ℝ) : ℝ := sortedAt x (n - 1)

lemma varMin_le (x : Fin n → ℝ) (i : Fin n) : varMin x ≤ x i := by
  obtain ⟨k, hk, hki⟩ := exists_sortedAt_eq x i
  rw [← hki]
  exact sortedAt_sorted x 0 k (Nat.zero_le k) hk

lemma le_varMax (x : Fin n → ℝ) (i : Fin n) : x i ≤ varMax x := by
  obtain ⟨k, hk, hki⟩ := exists_sortedAt_eq x i
  rw [← hki]
  exact sortedAt_sorted x k (n - 1) (by omega) (by omega)

lemma exists_eq_varMin (hn : 0 < n) (x : Fin n → ℝ) : ∃ i, x i = varMin x :=
  exists_eq_sortedAt x hn

lemma exists_eq_varMax (hn : 0 < n) (x : Fin n → ℝ) : ∃ i, x i = varMax x :=
  exists_eq_sortedAt x (by omega)

lemma varMin_lt_varMax_iff (hn : 0 < n) (x : Fin n → ℝ) :
    varMin x < varMax x ↔ ∃ i j, x i ≠ x j := by
  constructor
  · intro h
    obtain ⟨i, hi⟩ := exists_eq_varMin hn x
    obtain ⟨j, hj⟩ := exists_eq_varMax hn x
    exact ⟨i, j, by rw [hi, hj]; exact h.ne⟩
  · rintro ⟨i, j, hij⟩
    by_contra hc
    push Not at hc
    have h1 := le_antisymm ((le_varMax x i).trans hc) (varMin_le x i)
    have h2 := le_antisymm ((le_varMax x j).trans hc) (varMin_le x j)
    exact hij (h1.trans h2.symm)

theorem LPMratio_continuous {d : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j) :
    Continuous fun t => LPMratio d t x :=
  cdfCurve_continuous hd hx

theorem LPMratio_varMin {d : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j) :
    LPMratio d (varMin x) x = 0 := by
  obtain ⟨i, j, hij⟩ := hx
  have hn : 0 < n := Fin.pos i
  refine LPMratio_eq_zero_of_forall_ge hd (varMin_le x) ?_
  by_cases h : varMin x < x i
  · exact ⟨i, h⟩
  · refine ⟨j, lt_of_le_of_ne (varMin_le x j) fun e => hij ?_⟩
    rw [← e]
    exact le_antisymm (not_lt.mp h) (varMin_le x i)

theorem LPMratio_varMax (hn : 0 < n) (d : ℕ) (x : Fin n → ℝ) :
    LPMratio d (varMax x) x = 1 :=
  LPMratio_eq_one_of_forall_le hn d (le_varMax x)

theorem LPMratio_strictMonoOn_Icc {d : ℕ} (hd : d ≠ 0) (hn : 0 < n) (x : Fin n → ℝ) :
    StrictMonoOn (fun t => LPMratio d t x) (Set.Icc (varMin x) (varMax x)) := by
  refine (LPMratio_strictMonoOn hd x).mono fun t ht => ⟨?_, ?_⟩
  · obtain ⟨i, hi⟩ := exists_eq_varMin hn x
    exact ⟨i, hi ▸ ht.1⟩
  · obtain ⟨j, hj⟩ := exists_eq_varMax hn x
    exact ⟨j, hj ▸ ht.2⟩

theorem existsUnique_LPMratio_eq {d : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ}
    (hx : ∃ i j, x i ≠ x j) {p : ℝ} (hp : p ∈ Set.Icc (0 : ℝ) 1) :
    ∃! t, t ∈ Set.Icc (varMin x) (varMax x) ∧ LPMratio d t x = p := by
  obtain ⟨i, _, _⟩ := id hx
  have hn : 0 < n := Fin.pos i
  have hle : varMin x ≤ varMax x := (varMin_le x i).trans (le_varMax x i)
  have hivt := intermediate_value_Icc hle (LPMratio_continuous hd hx).continuousOn
  rw [LPMratio_varMin hd hx, LPMratio_varMax hn] at hivt
  obtain ⟨t, ht, htp⟩ := hivt hp
  refine ⟨t, ⟨ht, htp⟩, fun s hs => ?_⟩
  exact (LPMratio_strictMonoOn_Icc hd hn x).injOn hs.1 ht (hs.2.trans htp.symm)

theorem eq_of_LPMratio_eq {d : ℕ} (hd : d ≠ 0) (hn : 0 < n) {x : Fin n → ℝ} {s t : ℝ}
    (hs : s ∈ Set.Icc (varMin x) (varMax x))
    (ht : t ∈ Set.Icc (varMin x) (varMax x))
    (h : LPMratio d s x = LPMratio d t x) : s = t :=
  (LPMratio_strictMonoOn_Icc hd hn x).injOn hs ht h

theorem LPMratio_eq_iff_balance {d : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    (t p : ℝ) :
    LPMratio d t x = p ↔
      (1 - p) * ∑ i, lowerDev d t (x i) = p * ∑ i, upperDev d t (x i) := by
  obtain ⟨i0, j0, hij⟩ := hx
  have hne : ∃ i, x i ≠ t := by
    by_cases h : x i0 = t
    · exact ⟨j0, fun e => hij (h.trans e.symm)⟩
    · exact ⟨i0, h⟩
  have hpos : 0 < (∑ i, lowerDev d t (x i)) + ∑ i, upperDev d t (x i) := by
    have h := ratio_denom_pos (d := d) hne
    unfold LPM UPM at h
    rw [← add_div] at h
    exact (div_pos_iff_of_pos_right (Nat.cast_pos.mpr (Fin.pos i0))).mp h
  rw [LPMratio_eq_sum_div hd hne, div_eq_iff hpos.ne']
  constructor <;> intro h <;> linarith

theorem LPMratio_one_eq_iff_expectile {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j) (t p : ℝ) :
    LPMratio 1 t x = p ↔
      (1 - p) * ∑ i, lowerDev 1 t (x i) = p * ∑ i, upperDev 1 t (x i) :=
  LPMratio_eq_iff_balance one_ne_zero hx t p

def clampP (p : ℝ) : ℝ := min (max p 0) 1

def quantile7 (p : ℝ) (x : Fin n → ℝ) : ℝ :=
  let q := max 0 (min 1 p)
  let index : ℝ := 1 + ((n - 1 : ℕ) : ℝ) * q
  let lo := ⌊index⌋₊
  let hi := ⌈index⌉₊
  let qs := sortedAt x (lo - 1)
  if (lo : ℝ) < index ∧ sortedAt x (hi - 1) ≠ qs then
    (1 - (index - lo)) * qs + (index - lo) * sortedAt x (hi - 1)
  else qs

def LPMVaR0 (p : ℝ) (x : Fin n → ℝ) : ℝ := quantile7 (clampP p) x

def UPMVaR0 (p : ℝ) (x : Fin n → ℝ) : ℝ := quantile7 (1 - clampP p) x

lemma clampP_mem (p : ℝ) : clampP p ∈ Set.Icc (0 : ℝ) 1 :=
  ⟨le_min (le_max_right _ _) zero_le_one, min_le_right _ _⟩

lemma clampP_of_mem {p : ℝ} (hp : p ∈ Set.Icc (0 : ℝ) 1) : clampP p = p := by
  simp [clampP, hp.1, hp.2]

theorem UPMVaR0_eq_LPMVaR0 (p : ℝ) (x : Fin n → ℝ) : UPMVaR0 p x = LPMVaR0 (1 - clampP p) x := by
  have h := clampP_mem p
  rw [LPMVaR0, clampP_of_mem ⟨by linarith [h.2], by linarith [h.1]⟩, UPMVaR0]

theorem quantile7_eq {p : ℝ} (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    quantile7 p x = sortedAt x ⌊((n - 1 : ℕ) : ℝ) * p⌋₊ +
      (((n - 1 : ℕ) : ℝ) * p - ⌊((n - 1 : ℕ) : ℝ) * p⌋₊) *
        (sortedAt x ⌈((n - 1 : ℕ) : ℝ) * p⌉₊ - sortedAt x ⌊((n - 1 : ℕ) : ℝ) * p⌋₊) := by
  set h : ℝ := ((n - 1 : ℕ) : ℝ) * p with hh
  have h0 : 0 ≤ h := mul_nonneg (Nat.cast_nonneg _) hp.1
  have hq : max 0 (min 1 p) = p := by simp [hp.1, hp.2]
  have hfl : ⌊1 + h⌋₊ = ⌊h⌋₊ + 1 := by rw [add_comm]; exact Nat.floor_add_one h0
  have hce : ⌈1 + h⌉₊ = ⌈h⌉₊ + 1 := by rw [add_comm]; exact Nat.ceil_add_one h0
  simp only [quantile7, hq, ← hh, hfl, hce, Nat.add_sub_cancel, Nat.cast_add, Nat.cast_one]
  split_ifs with hc
  · ring
  · push Not at hc
    by_cases hlt : ((⌊h⌋₊ : ℝ) + 1) < 1 + h
    · rw [hc hlt]; ring
    · have : h = ⌊h⌋₊ := le_antisymm (by linarith [not_lt.mp hlt]) (Nat.floor_le h0)
      rw [← this]; ring

private lemma sortedAt_fin3 (k : ℕ) (hk : k < 3) : sortedAt ![(0 : ℝ), 1, 2] k = k := by
  have hmono : Monotone ![(0 : ℝ), 1, 2] := by
    refine Fin.monotone_iff_le_succ.mpr fun i => ?_
    fin_cases i <;> simp
  have hs := Tuple.sort_eq_refl_iff_monotone.mpr hmono
  simp only [sortedAt, dite_eq_left hk, hs, Equiv.refl_apply]
  interval_cases k <;> simp

theorem LPMVaR0_type7_example :
    LPMVaR0 (9 / 10) ![(0 : ℝ), 1, 2] = 9 / 5 ∧
      LPMratio 0 (9 / 5) ![(0 : ℝ), 1, 2] = 2 / 3 ∧
      LPMratio 0 (LPMVaR0 (9 / 10) ![(0 : ℝ), 1, 2]) ![(0 : ℝ), 1, 2] < 9 / 10 := by
  have hv : LPMVaR0 (9 / 10) ![(0 : ℝ), 1, 2] = 9 / 5 := by
    rw [LPMVaR0, clampP_of_mem ⟨by norm_num, by norm_num⟩, quantile7_eq ⟨by norm_num, by norm_num⟩]
    have hf : ⌊((3 - 1 : ℕ) : ℝ) * (9 / 10)⌋₊ = 1 := by
      rw [Nat.floor_eq_iff (by norm_num)]; norm_num
    have hc : ⌈((3 - 1 : ℕ) : ℝ) * (9 / 10)⌉₊ = 2 := by
      rw [Nat.ceil_eq_iff (by norm_num)]; norm_num
    rw [hf, hc, sortedAt_fin3 1 (by norm_num), sortedAt_fin3 2 (by norm_num)]
    norm_num
  have hr : LPMratio 0 (9 / 5) ![(0 : ℝ), 1, 2] = 2 / 3 := by
    simp only [LPMratio, ↓reduceIte, LPM, lowerDev, Fin.sum_univ_three, pow_zero]
    norm_num
  refine ⟨hv, hr, ?_⟩
  rw [hv, hr]
  norm_num

theorem LPMVaR0_not_minimal :
    LPMVaR0 (1 / 10) ![(0 : ℝ), 1, 2] = 1 / 5 ∧
      1 / 10 ≤ LPMratio 0 0 ![(0 : ℝ), 1, 2] ∧ (0 : ℝ) < 1 / 5 := by
  refine ⟨?_, ?_, by norm_num⟩
  · rw [LPMVaR0, clampP_of_mem ⟨by norm_num, by norm_num⟩,
      quantile7_eq ⟨by norm_num, by norm_num⟩]
    have hf : ⌊((3 - 1 : ℕ) : ℝ) * (1 / 10)⌋₊ = 0 := by
      rw [Nat.floor_eq_iff (by norm_num)]; norm_num
    have hc : ⌈((3 - 1 : ℕ) : ℝ) * (1 / 10)⌉₊ = 1 := by
      rw [Nat.ceil_eq_iff (by norm_num)]; norm_num
    rw [hf, hc, sortedAt_fin3 0 (by norm_num), sortedAt_fin3 1 (by norm_num)]
    norm_num
  · simp only [LPMratio, ↓reduceIte, LPM, lowerDev, Fin.sum_univ_three, pow_zero]
    norm_num

def ecdfInv (p : ℝ) (x : Fin n → ℝ) : ℝ := sortedAt x (⌈(n : ℝ) * p⌉₊ - 1)

theorem le_LPMratio_zero_iff {p : ℝ} (hp0 : 0 < p) (hp1 : p ≤ 1) (hn : 0 < n)
    (x : Fin n → ℝ) (t : ℝ) : p ≤ LPMratio 0 t x ↔ ecdfInv p x ≤ t := by
  have hnr : (0 : ℝ) < n := Nat.cast_pos.mpr hn
  have hc1 : 1 ≤ ⌈(n : ℝ) * p⌉₊ := Nat.one_le_iff_ne_zero.mpr
    (Nat.ceil_eq_zero.not.mpr (not_le.mpr (mul_pos hnr hp0)))
  have hcn : ⌈(n : ℝ) * p⌉₊ ≤ n :=
    Nat.ceil_le.mpr (by simpa using mul_le_of_le_one_right hnr.le hp1)
  rw [LPMratio_zero_eq_ecdf, ← countLe_sortedAt, le_div_iff₀ hnr, ecdfInv,
    ← lt_countLe_iff (sortedAt_sorted x) t (by omega), mul_comm, ← Nat.ceil_le]
  omega

def rTol : ℝ := 1 / 2 ^ 39

def uniqSet (a : ℕ → ℝ) (m : ℕ) : Finset ℝ := (range m).image a

def nUniq (a : ℕ → ℝ) (m : ℕ) : ℕ := (uniqSet a m).card

def uniqAt (a : ℕ → ℝ) (m : ℕ) (k : ℕ) : ℝ :=
  if h : k < nUniq a m then (uniqSet a m).orderEmbOfFin rfl ⟨k, h⟩ else 0

def kBreak (a : ℕ → ℝ) (m : ℕ) (k : ℕ) : ℕ := countLe a m (uniqAt a m k)

def prefixPower (a : ℕ → ℝ) (K : ℕ) (j : ℕ) : ℝ := ∑ r ∈ range K, a r ^ j

def totalPower (a : ℕ → ℝ) (m : ℕ) (j : ℕ) : ℝ := ∑ r ∈ range m, a r ^ j

def lpmRaw (d : ℕ) (t : ℝ) (P : ℕ → ℝ) : ℝ :=
  ∑ j ∈ range (d + 1), (d.choose j : ℝ) * (-1) ^ j * t ^ (d - j) * P j

def upmRaw (d : ℕ) (t : ℝ) (P T : ℕ → ℝ) : ℝ :=
  ∑ j ∈ range (d + 1), (d.choose j : ℝ) * (-1) ^ (d - j) * t ^ (d - j) * (T j - P j)

def ratioBreakRaw (d : ℕ) (a : ℕ → ℝ) (m : ℕ) (k : ℕ) : ℝ :=
  lpmRaw d (uniqAt a m k) (prefixPower a (kBreak a m k)) /
    (lpmRaw d (uniqAt a m k) (prefixPower a (kBreak a m k)) +
      upmRaw d (uniqAt a m k) (prefixPower a (kBreak a m k)) (totalPower a m))

def cummaxSeq (r : ℕ → ℝ) : ℕ → ℝ
  | 0 => r 0
  | k + 1 => max (cummaxSeq r k) (r (k + 1))

def ratioBreak (d : ℕ) (a : ℕ → ℝ) (m : ℕ) (k : ℕ) : ℝ :=
  if k = nUniq a m - 1 then 1
  else if k = 0 then 0
  else cummaxSeq (fun l => clampP (ratioBreakRaw d a m l)) k

def locate (d : ℕ) (a : ℕ → ℝ) (m : ℕ) (p : ℝ) : ℕ :=
  let N := nUniq a m
  let fi := countLe (ratioBreak d a m) N p
  let fi' := if fi = N ∧ p = ratioBreak d a m (N - 1) then N - 1 else fi
  min (max fi' 1) (N - 1)

open Classical in
def unirootIdeal (f : ℝ → ℝ) (lo hi : ℝ) : ℝ :=
  if h : ∃ t ∈ Set.Icc lo hi, f t = 0 then h.choose else lo

def solveInterval (tol : ℝ) (d : ℕ) (a : ℕ → ℝ) (m : ℕ) (p : ℝ) (i : ℕ) : ℝ :=
  let lower := uniqAt a m (i - 1)
  let upper := uniqAt a m i
  let P := prefixPower a (kBreak a m (i - 1))
  let T := totalPower a m
  let f : ℝ → ℝ := fun t => (1 - p) * lpmRaw d t P - p * upmRaw d t P T
  let gap : ℝ → ℝ := fun t => lpmRaw d t P / (lpmRaw d t P + upmRaw d t P T) - p
  if lower = upper then lower
  else if |gap lower| ≤ tol then lower
  else if |gap upper| ≤ tol then upper
  else if f lower * f upper ≤ 0 then min (max (unirootIdeal f lower upper) lower) upper
  else if |gap lower| ≤ |gap upper| then lower else upper

def intVaRCore (tol : ℝ) (d : ℕ) (a : ℕ → ℝ) (m : ℕ) (p : ℝ) : ℝ :=
  if p ≤ 0 then a 0
  else if 1 ≤ p then a (m - 1)
  else solveInterval tol d a m p (locate d a m p)

def vshift (x : Fin n → ℝ) : ℝ := sortedAt x ((n + 1) / 2 - 1)

def shiftedArr (x : Fin n → ℝ) : ℕ → ℝ := fun r => sortedAt x r - vshift x

def intLPMVaR (tol : ℝ) (d : ℕ) (p : ℝ) (x : Fin n → ℝ) : ℝ :=
  if n = 1 ∨ shiftedArr x 0 = shiftedArr x (n - 1) then shiftedArr x 0 + vshift x
  else intVaRCore tol d (shiftedArr x) n (clampP p) + vshift x

def intUPMVaR (tol : ℝ) (d : ℕ) (p : ℝ) (x : Fin n → ℝ) : ℝ :=
  intLPMVaR tol d (1 - clampP p) x

section Uniq

variable {a : ℕ → ℝ} {m : ℕ}

lemma uniqAt_strictMono {k l : ℕ} (hkl : k < l) (hl : l < nUniq a m) :
    uniqAt a m k < uniqAt a m l := by
  have hk : k < nUniq a m := hkl.trans hl
  simp only [uniqAt, dite_eq_left hk, dite_eq_left hl]
  exact ((uniqSet a m).orderEmbOfFin rfl).strictMono (show (⟨k, hk⟩ : Fin _) < ⟨l, hl⟩ from hkl)

lemma uniqAt_mono {k l : ℕ} (hkl : k ≤ l) (hl : l < nUniq a m) :
    uniqAt a m k ≤ uniqAt a m l := by
  rcases hkl.lt_or_eq with h | h
  · exact (uniqAt_strictMono h hl).le
  · rw [h]

lemma exists_eq_uniqAt {k : ℕ} (hk : k < nUniq a m) : ∃ r < m, a r = uniqAt a m k := by
  have hmem := (uniqSet a m).orderEmbOfFin_mem rfl ⟨k, hk⟩
  simp only [uniqSet, mem_image, mem_range] at hmem
  obtain ⟨r, hr, hra⟩ := hmem
  exact ⟨r, hr, by rw [uniqAt, dite_eq_left hk]; exact hra⟩

lemma exists_uniqAt_eq {r : ℕ} (hr : r < m) : ∃ k < nUniq a m, uniqAt a m k = a r := by
  have hmem : a r ∈ Set.range ((uniqSet a m).orderEmbOfFin rfl) := by
    rw [Finset.range_orderEmbOfFin]
    exact mem_coe.mpr (mem_image_of_mem a (mem_range.mpr hr))
  obtain ⟨⟨k, hk⟩, hka⟩ := hmem
  have hk' : k < nUniq a m := hk
  exact ⟨k, hk', by rw [uniqAt, dite_eq_left hk']; exact hka⟩

lemma lt_nUniq_index {k l : ℕ} (h : uniqAt a m k < uniqAt a m l)
    (hk : k < nUniq a m) : k < l := by
  by_contra hc
  push Not at hc
  exact absurd (uniqAt_mono hc hk) (not_le.mpr h)

lemma nUniq_pos (hm : 0 < m) : 0 < nUniq a m :=
  card_pos.mpr ⟨a 0, mem_image_of_mem a (mem_range.mpr hm)⟩

lemma uniqAt_zero (hs : VaRSorted a m) (hm : 0 < m) : uniqAt a m 0 = a 0 := by
  obtain ⟨r, hr, hra⟩ := exists_eq_uniqAt (a := a) (nUniq_pos hm)
  obtain ⟨k, hk, hka⟩ := exists_uniqAt_eq (a := a) hm
  refine le_antisymm ?_ (hra ▸ hs 0 r (Nat.zero_le r) hr)
  rw [← hka]
  exact uniqAt_mono (Nat.zero_le k) hk

lemma uniqAt_last (hs : VaRSorted a m) (hm : 0 < m) :
    uniqAt a m (nUniq a m - 1) = a (m - 1) := by
  have hN := nUniq_pos (a := a) hm
  obtain ⟨r, hr, hra⟩ := exists_eq_uniqAt (a := a) (show nUniq a m - 1 < nUniq a m by omega)
  obtain ⟨k, hk, hka⟩ := exists_uniqAt_eq (a := a) (show m - 1 < m by omega)
  refine le_antisymm (hra ▸ hs r (m - 1) (by omega) (by omega)) ?_
  rw [← hka]
  exact uniqAt_mono (by omega) (by omega)

lemma two_le_nUniq (hs : VaRSorted a m) (hm : 0 < m) (hlt : a 0 < a (m - 1)) :
    2 ≤ nUniq a m := by
  have hN := nUniq_pos (a := a) hm
  by_contra hc
  have h1 : nUniq a m - 1 = 0 := by omega
  have := uniqAt_last hs hm
  rw [h1, uniqAt_zero hs hm] at this
  exact hlt.ne this

lemma uniqAt_succ_le {k r : ℕ} (hk : k + 1 < nUniq a m) (hr : r < m)
    (h : uniqAt a m k < a r) : uniqAt a m (k + 1) ≤ a r := by
  obtain ⟨l, hl, hla⟩ := exists_uniqAt_eq (a := a) hr
  rw [← hla] at h ⊢
  exact uniqAt_mono (lt_nUniq_index h (by omega)) hl

end Uniq

theorem lpmRaw_prefixPower (d : ℕ) (t : ℝ) (a : ℕ → ℝ) (K : ℕ) :
    lpmRaw d t (prefixPower a K) = ∑ r ∈ range K, (t - a r) ^ d := by
  unfold lpmRaw prefixPower
  simp_rw [mul_sum]
  rw [sum_comm]
  refine sum_congr rfl fun r _ => ?_
  rw [sub_eq_neg_add, add_pow]
  refine sum_congr rfl fun j _ => ?_
  rw [neg_pow]
  ring

theorem upmRaw_prefixPower (d : ℕ) (t : ℝ) (a : ℕ → ℝ) {K m : ℕ} (hK : K ≤ m) :
    upmRaw d t (prefixPower a K) (totalPower a m) = ∑ r ∈ Ico K m, (a r - t) ^ d := by
  have hsuf : ∀ j, totalPower a m j - prefixPower a K j = ∑ r ∈ Ico K m, a r ^ j := by
    intro j
    unfold totalPower prefixPower
    rw [← sum_range_add_sum_Ico _ hK]
    ring
  unfold upmRaw
  simp_rw [hsuf, mul_sum]
  rw [sum_comm]
  refine sum_congr rfl fun r _ => ?_
  rw [sub_eq_add_neg, add_pow]
  refine sum_congr rfl fun j hj => ?_
  rw [neg_pow]
  ring

private lemma lowerDev_eq_max_pow' {d : ℕ} (hd : d ≠ 0) (t v : ℝ) :
    lowerDev d t v = (max (t - v) 0) ^ d := by
  unfold lowerDev
  split_ifs with h
  · rw [max_eq_left (sub_nonneg.mpr h)]
  · rw [max_eq_right (by linarith [not_le.mp h]), zero_pow hd]

private lemma upperDev_eq_max_pow' {d : ℕ} (hd : d ≠ 0) (t v : ℝ) :
    upperDev d t v = (max (v - t) 0) ^ d := by
  unfold upperDev
  split_ifs with h
  · rw [max_eq_left (sub_nonneg.mpr h.le)]
  · rw [max_eq_right (by linarith [not_lt.mp h]), zero_pow hd]

def arrSample (a : ℕ → ℝ) (m : ℕ) : Fin m → ℝ := fun i => a i

section Interval

variable {d : ℕ} {a : ℕ → ℝ} {m : ℕ}

lemma kBreak_le (a : ℕ → ℝ) (m k : ℕ) : kBreak a m k ≤ m := countLe_le a m _

theorem lpmRaw_eq_sum_lowerDev (hd : d ≠ 0) (hs : VaRSorted a m) (k : ℕ) {t : ℝ}
    (hlo : uniqAt a m k ≤ t) (hhi : ∀ r < m, uniqAt a m k < a r → t ≤ a r) :
    lpmRaw d t (prefixPower a (kBreak a m k)) = ∑ r ∈ range m, lowerDev d t (a r) := by
  rw [lpmRaw_prefixPower, ← sum_range_add_sum_Ico _ (kBreak_le a m k)]
  have h2 : ∑ r ∈ Ico (kBreak a m k) m, lowerDev d t (a r) = 0 := by
    refine sum_eq_zero fun r hr => ?_
    have hr := mem_Ico.mp hr
    have hgt : uniqAt a m k < a r := not_le.mp fun h =>
      absurd ((lt_countLe_iff hs _ hr.2).mpr h) (not_lt.mpr hr.1)
    rw [lowerDev_eq_max_pow' hd, max_eq_right (by linarith [hhi r hr.2 hgt]), zero_pow hd]
  rw [h2, add_zero]
  refine sum_congr rfl fun r hr => ?_
  have hr := mem_range.mp hr
  have hle : a r ≤ uniqAt a m k :=
    (lt_countLe_iff hs _ (lt_of_lt_of_le hr (countLe_le a m _))).mp hr
  simp [lowerDev, hle.trans hlo]

theorem upmRaw_eq_sum_upperDev (hd : d ≠ 0) (hs : VaRSorted a m) (k : ℕ) {t : ℝ}
    (hlo : uniqAt a m k ≤ t) (hhi : ∀ r < m, uniqAt a m k < a r → t ≤ a r) :
    upmRaw d t (prefixPower a (kBreak a m k)) (totalPower a m) =
      ∑ r ∈ range m, upperDev d t (a r) := by
  rw [upmRaw_prefixPower d t a (kBreak_le a m k),
    ← sum_range_add_sum_Ico _ (kBreak_le a m k)]
  have h1 : ∑ r ∈ range (kBreak a m k), upperDev d t (a r) = 0 := by
    refine sum_eq_zero fun r hr => ?_
    have hr := mem_range.mp hr
    have hle : a r ≤ uniqAt a m k :=
      (lt_countLe_iff hs _ (lt_of_lt_of_le hr (countLe_le a m _))).mp hr
    simp [upperDev, not_lt.mpr (hle.trans hlo)]
  rw [h1, zero_add]
  refine sum_congr rfl fun r hr => ?_
  have hr := mem_Ico.mp hr
  have hgt : uniqAt a m k < a r := not_le.mp fun h =>
    absurd ((lt_countLe_iff hs _ hr.2).mpr h) (not_lt.mpr hr.1)
  rw [upperDev_eq_max_pow' hd, max_eq_left (by linarith [hhi r hr.2 hgt])]

lemma admissible_of_mem {k : ℕ} (hk : k + 1 < nUniq a m) {t : ℝ}
    (ht : t ≤ uniqAt a m (k + 1)) : ∀ r < m, uniqAt a m k < a r → t ≤ a r :=
  fun _ hr h => ht.trans (uniqAt_succ_le hk hr h)

lemma admissible_self (k : ℕ) : ∀ r < m, uniqAt a m k < a r → uniqAt a m k ≤ a r :=
  fun _ _ h => h.le

end Interval

lemma continuous_lpmRaw (d : ℕ) (P : ℕ → ℝ) : Continuous fun t => lpmRaw d t P := by
  unfold lpmRaw
  fun_prop

lemma continuous_upmRaw (d : ℕ) (P T : ℕ → ℝ) : Continuous fun t => upmRaw d t P T := by
  unfold upmRaw
  fun_prop

section Core

variable {d : ℕ} {a : ℕ → ℝ} {m : ℕ}

lemma sum_arrSample (f : ℝ → ℝ) : ∑ i, f (arrSample a m i) = ∑ r ∈ range m, f (a r) :=
  Fin.sum_univ_eq_sum_range (fun r => f (a r)) m

lemma arrSample_ne (hm : 0 < m) (hlt : a 0 < a (m - 1)) (t : ℝ) : ∃ i, arrSample a m i ≠ t := by
  by_cases h : a 0 = t
  · exact ⟨⟨m - 1, by omega⟩, fun e => hlt.ne (h.trans e.symm)⟩
  · exact ⟨⟨0, hm⟩, h⟩

lemma uniqAt_mem_range {k : ℕ} (hk : k < nUniq a m) :
    (∃ i, arrSample a m i ≤ uniqAt a m k) ∧ ∃ j, uniqAt a m k ≤ arrSample a m j := by
  obtain ⟨r, hr, hra⟩ := exists_eq_uniqAt hk
  exact ⟨⟨⟨r, hr⟩, hra.le⟩, ⟨⟨r, hr⟩, hra.ge⟩⟩

lemma raw_ratio_eq (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m) (hlt : a 0 < a (m - 1))
    (k : ℕ) {t : ℝ} (hlo : uniqAt a m k ≤ t) (hhi : ∀ r < m, uniqAt a m k < a r → t ≤ a r) :
    lpmRaw d t (prefixPower a (kBreak a m k)) /
        (lpmRaw d t (prefixPower a (kBreak a m k)) +
          upmRaw d t (prefixPower a (kBreak a m k)) (totalPower a m)) =
      LPMratio d t (arrSample a m) ∧
    0 < lpmRaw d t (prefixPower a (kBreak a m k)) +
          upmRaw d t (prefixPower a (kBreak a m k)) (totalPower a m) := by
  rw [lpmRaw_eq_sum_lowerDev hd hs k hlo hhi, upmRaw_eq_sum_upperDev hd hs k hlo hhi,
    LPMratio_eq_sum_div hd (arrSample_ne hm hlt t), sum_arrSample, sum_arrSample]
  refine ⟨rfl, ?_⟩
  have hpos := ratio_denom_pos (d := d) (arrSample_ne hm hlt t)
  unfold LPM UPM at hpos
  rw [← add_div, sum_arrSample, sum_arrSample] at hpos
  exact (div_pos_iff_of_pos_right (Nat.cast_pos.mpr hm)).mp hpos

lemma ratioBreakRaw_eq (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m) (hlt : a 0 < a (m - 1))
    (k : ℕ) : ratioBreakRaw d a m k = LPMratio d (uniqAt a m k) (arrSample a m) :=
  (raw_ratio_eq hd hs hm hlt k le_rfl (admissible_self k)).1

lemma ratioBreakRaw_lt (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m)
    (hlt : a 0 < a (m - 1)) {k l : ℕ} (hkl : k < l) (hl : l < nUniq a m) :
    ratioBreakRaw d a m k < ratioBreakRaw d a m l := by
  rw [ratioBreakRaw_eq hd hs hm hlt, ratioBreakRaw_eq hd hs hm hlt]
  exact LPMratio_lt_LPMratio hd (uniqAt_mem_range (hkl.trans hl)).1 (uniqAt_mem_range hl).2
    (uniqAt_strictMono hkl hl)

lemma ratioBreakRaw_zero (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m)
    (hlt : a 0 < a (m - 1)) : ratioBreakRaw d a m 0 = 0 := by
  rw [ratioBreakRaw_eq hd hs hm hlt, uniqAt_zero hs hm]
  refine LPMratio_eq_zero_of_forall_ge hd (fun i => hs 0 i (Nat.zero_le _) i.isLt)
    ⟨⟨m - 1, by omega⟩, hlt⟩

lemma ratioBreakRaw_last (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m)
    (hlt : a 0 < a (m - 1)) : ratioBreakRaw d a m (nUniq a m - 1) = 1 := by
  rw [ratioBreakRaw_eq hd hs hm hlt, uniqAt_last hs hm]
  exact LPMratio_eq_one_of_forall_le hm d fun i => hs i (m - 1) (by omega) (by omega)

lemma cummaxSeq_eq {r : ℕ → ℝ} {K : ℕ} (hr : ∀ k, k + 1 ≤ K → r k ≤ r (k + 1)) :
    ∀ k, k ≤ K → cummaxSeq r k = r k
  | 0, _ => rfl
  | k + 1, hk => by
    rw [cummaxSeq, cummaxSeq_eq hr k (by omega), max_eq_right (hr k hk)]

theorem ratioBreak_eq (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m) (hlt : a 0 < a (m - 1))
    {k : ℕ} (hk : k < nUniq a m) : ratioBreak d a m k = ratioBreakRaw d a m k := by
  unfold ratioBreak
  split_ifs with h1 h2
  · rw [h1, ratioBreakRaw_last hd hs hm hlt]
  · rw [h2, ratioBreakRaw_zero hd hs hm hlt]
  · have hmono : ∀ l, l + 1 ≤ nUniq a m - 1 →
        clampP (ratioBreakRaw d a m l) ≤ clampP (ratioBreakRaw d a m (l + 1)) := by
      intro l hl
      have := (ratioBreakRaw_lt hd hs hm hlt (Nat.lt_succ_self l) (by omega)).le
      unfold clampP
      exact min_le_min_right _ (max_le_max_right _ this)
    rw [cummaxSeq_eq hmono k (by omega), ratioBreakRaw_eq hd hs hm hlt,
      clampP_of_mem (LPMratio_mem_Icc _ _ _)]

theorem ratioBreak_lt (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m) (hlt : a 0 < a (m - 1))
    {k l : ℕ} (hkl : k < l) (hl : l < nUniq a m) :
    ratioBreak d a m k < ratioBreak d a m l := by
  rw [ratioBreak_eq hd hs hm hlt (hkl.trans hl), ratioBreak_eq hd hs hm hlt hl]
  exact ratioBreakRaw_lt hd hs hm hlt hkl hl

lemma ratioBreak_sorted (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m)
    (hlt : a 0 < a (m - 1)) : VaRSorted (ratioBreak d a m) (nUniq a m) := by
  intro k l hkl hl
  rcases hkl.lt_or_eq with h | h
  · exact (ratioBreak_lt hd hs hm hlt h hl).le
  · rw [h]

private lemma countLe_ratioBreak_bounds (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m)
    (hlt : a 0 < a (m - 1)) {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) :
    0 < countLe (ratioBreak d a m) (nUniq a m) p ∧
      countLe (ratioBreak d a m) (nUniq a m) p < nUniq a m := by
  have hN := two_le_nUniq hs hm hlt
  have hsorted := ratioBreak_sorted hd hs hm hlt
  have h0 : ratioBreak d a m 0 = 0 := by
    rw [ratioBreak_eq hd hs hm hlt (by omega), ratioBreakRaw_zero hd hs hm hlt]
  have hlast : ratioBreak d a m (nUniq a m - 1) = 1 := by
    rw [ratioBreak_eq hd hs hm hlt (by omega), ratioBreakRaw_last hd hs hm hlt]
  refine ⟨(lt_countLe_iff hsorted p (by omega)).mpr (h0 ▸ hp0.le), ?_⟩
  by_contra hc
  have := (lt_countLe_iff hsorted p (show nUniq a m - 1 < nUniq a m by omega)).mp (by omega)
  rw [hlast] at this
  linarith

theorem locate_eq_countLe (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m)
    (hlt : a 0 < a (m - 1)) {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) :
    locate d a m p = countLe (ratioBreak d a m) (nUniq a m) p := by
  obtain ⟨hpos, hlt'⟩ := countLe_ratioBreak_bounds hd hs hm hlt hp0 hp1
  simp only [locate]
  rw [ite_eq_right (by omega)]
  omega

theorem locate_spec (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m) (hlt : a 0 < a (m - 1))
    {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) :
    1 ≤ locate d a m p ∧ locate d a m p < nUniq a m ∧
      ratioBreak d a m (locate d a m p - 1) ≤ p ∧ p < ratioBreak d a m (locate d a m p) := by
  have hsorted := ratioBreak_sorted hd hs hm hlt
  obtain ⟨hfi_pos, hfi_lt⟩ := countLe_ratioBreak_bounds hd hs hm hlt hp0 hp1
  rw [locate_eq_countLe hd hs hm hlt hp0 hp1]
  refine ⟨hfi_pos, hfi_lt, (lt_countLe_iff hsorted p (by omega)).mp (by omega), ?_⟩
  exact not_le.mp fun h => absurd ((lt_countLe_iff hsorted p hfi_lt).mpr h) (lt_irrefl _)

theorem solveInterval_spec (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m)
    (hlt : a 0 < a (m - 1)) {tol p : ℝ} (htol : 0 ≤ tol) {i : ℕ} (hi1 : 1 ≤ i)
    (hiN : i < nUniq a m) (hlo : ratioBreak d a m (i - 1) ≤ p)
    (hhi : p < ratioBreak d a m i) :
    solveInterval tol d a m p i ∈ Set.Icc (uniqAt a m (i - 1)) (uniqAt a m i) ∧
      |LPMratio d (solveInterval tol d a m p i) (arrSample a m) - p| ≤ tol := by
  have hk1 : i - 1 + 1 = i := by omega
  have hlu : uniqAt a m (i - 1) < uniqAt a m i := uniqAt_strictMono (by omega) hiN
  set lower := uniqAt a m (i - 1) with hlower
  set upper := uniqAt a m i with hupper
  set P := prefixPower a (kBreak a m (i - 1)) with hP
  set T := totalPower a m with hT
  have hadm : ∀ t, t ≤ upper → ∀ r < m, lower < a r → t ≤ a r := fun t ht =>
    admissible_of_mem (k := i - 1) (by omega) (by rw [hk1]; exact ht)
  have hR : ∀ t, lower ≤ t → t ≤ upper →
      lpmRaw d t P / (lpmRaw d t P + upmRaw d t P T) = LPMratio d t (arrSample a m) ∧
        0 < lpmRaw d t P + upmRaw d t P T := fun t h1 h2 =>
    raw_ratio_eq hd hs hm hlt (i - 1) h1 (hadm t h2)
  have hf : ∀ t, lower ≤ t → t ≤ upper →
      (1 - p) * lpmRaw d t P - p * upmRaw d t P T =
        (lpmRaw d t P + upmRaw d t P T) * (LPMratio d t (arrSample a m) - p) := by
    intro t h1 h2
    obtain ⟨hr, hpos⟩ := hR t h1 h2
    rw [← hr, mul_sub, mul_div_cancel₀ _ hpos.ne']
    ring
  have hrl : LPMratio d lower (arrSample a m) = ratioBreak d a m (i - 1) := by
    rw [ratioBreak_eq hd hs hm hlt (by omega), ratioBreakRaw_eq hd hs hm hlt]
  have hru : LPMratio d upper (arrSample a m) = ratioBreak d a m i := by
    rw [ratioBreak_eq hd hs hm hlt hiN, ratioBreakRaw_eq hd hs hm hlt]
  have hfl : (1 - p) * lpmRaw d lower P - p * upmRaw d lower P T ≤ 0 := by
    rw [hf lower le_rfl hlu.le]
    exact mul_nonpos_of_nonneg_of_nonpos (hR lower le_rfl hlu.le).2.le (by linarith)
  have hfu : 0 ≤ (1 - p) * lpmRaw d upper P - p * upmRaw d upper P T := by
    rw [hf upper hlu.le le_rfl]
    exact mul_nonneg (hR upper hlu.le le_rfl).2.le (by linarith)
  have hgl : lpmRaw d lower P / (lpmRaw d lower P + upmRaw d lower P T) - p =
      LPMratio d lower (arrSample a m) - p := by rw [(hR lower le_rfl hlu.le).1]
  have hgu : lpmRaw d upper P / (lpmRaw d upper P + upmRaw d upper P T) - p =
      LPMratio d upper (arrSample a m) - p := by rw [(hR upper hlu.le le_rfl).1]
  simp only [solveInterval, ← hlower, ← hupper, ← hP, ← hT]
  split_ifs with h1 h2 h3 h4
  · exact absurd h1 hlu.ne
  · exact ⟨⟨le_rfl, hlu.le⟩, hgl ▸ h2⟩
  · exact ⟨⟨hlu.le, le_rfl⟩, hgu ▸ h3⟩
  · set f : ℝ → ℝ := fun t => (1 - p) * lpmRaw d t P - p * upmRaw d t P T with hfdef
    have hcont : Continuous f :=
      ((continuous_const.mul (continuous_lpmRaw d P)).sub
        (continuous_const.mul (continuous_upmRaw d P T)))
    have hex : ∃ t ∈ Set.Icc lower upper, f t = 0 :=
      intermediate_value_Icc hlu.le hcont.continuousOn ⟨hfl, hfu⟩
    have hroot : unirootIdeal f lower upper ∈ Set.Icc lower upper ∧
        f (unirootIdeal f lower upper) = 0 := by
      rw [unirootIdeal, dite_eq_left hex]
      exact hex.choose_spec
    set t := unirootIdeal f lower upper
    have hclamp : min (max t lower) upper = t := by
      rw [max_eq_left hroot.1.1, min_eq_left hroot.1.2]
    rw [hclamp]
    refine ⟨hroot.1, ?_⟩
    have h0 := hroot.2
    simp only [hfdef] at h0
    rw [hf t hroot.1.1 hroot.1.2] at h0
    have hpos := (hR t hroot.1.1 hroot.1.2).2
    have : LPMratio d t (arrSample a m) - p = 0 := by
      rcases mul_eq_zero.mp h0 with h | h
      · exact absurd h hpos.ne'
      · exact h
    rw [this, abs_zero]
    exact htol
  all_goals exact absurd (mul_nonpos_of_nonpos_of_nonneg hfl hfu) h4

theorem intVaRCore_spec (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m)
    (hlt : a 0 < a (m - 1)) {tol p : ℝ} (htol : 0 ≤ tol) (hp : p ∈ Set.Icc (0 : ℝ) 1) :
    intVaRCore tol d a m p ∈ Set.Icc (a 0) (a (m - 1)) ∧
      |LPMratio d (intVaRCore tol d a m p) (arrSample a m) - p| ≤ tol := by
  have hle : a 0 ≤ a (m - 1) := hlt.le
  unfold intVaRCore
  split_ifs with h0 h1
  · have hp0 : p = 0 := le_antisymm h0 hp.1
    refine ⟨⟨le_rfl, hle⟩, ?_⟩
    rw [hp0, LPMratio_eq_zero_of_forall_ge (x := arrSample a m) hd
      (fun i => hs 0 i (Nat.zero_le _) i.isLt)
      ⟨⟨m - 1, by omega⟩, hlt⟩]
    simpa using htol
  · have hp1 : p = 1 := le_antisymm hp.2 h1
    refine ⟨⟨hle, le_rfl⟩, ?_⟩
    rw [hp1, LPMratio_eq_one_of_forall_le (x := arrSample a m) hm d
      fun i => hs i (m - 1) (by omega) (by omega)]
    simpa using htol
  · push Not at h0 h1
    obtain ⟨hi1, hiN, hlo, hhi⟩ := locate_spec hd hs hm hlt h0 h1
    obtain ⟨hmem, hres⟩ := solveInterval_spec hd hs hm hlt htol hi1 hiN hlo hhi
    refine ⟨⟨?_, ?_⟩, hres⟩
    · rw [← uniqAt_zero hs hm]
      exact (uniqAt_mono (Nat.zero_le _) (by omega)).trans hmem.1
    · rw [← uniqAt_last hs hm]
      exact hmem.2.trans (uniqAt_mono (by omega) (by omega))

theorem intVaRCore_same_interval (hd : d ≠ 0) (hs : VaRSorted a m) (hm : 0 < m)
    (hlt : a 0 < a (m - 1)) {tol p : ℝ} (htol : 0 ≤ tol) (hp0 : 0 < p) (hp1 : p < 1) {t : ℝ}
    (ht : t ∈ Set.Icc (a 0) (a (m - 1))) (htp : LPMratio d t (arrSample a m) = p) :
    intVaRCore tol d a m p ∈
        Set.Icc (uniqAt a m (locate d a m p - 1)) (uniqAt a m (locate d a m p)) ∧
      t ∈ Set.Icc (uniqAt a m (locate d a m p - 1)) (uniqAt a m (locate d a m p)) := by
  obtain ⟨hi1, hiN, hlo, hhi⟩ := locate_spec hd hs hm hlt hp0 hp1
  have hres : intVaRCore tol d a m p = solveInterval tol d a m p (locate d a m p) := by
    rw [intVaRCore, ite_eq_right (not_le.mpr hp0), ite_eq_right (not_le.mpr hp1)]
  refine ⟨hres ▸ (solveInterval_spec hd hs hm hlt htol hi1 hiN hlo hhi).1, ⟨?_, ?_⟩⟩
  · by_contra hc
    push Not at hc
    have h := LPMratio_lt_LPMratio hd (x := arrSample a m) ⟨⟨0, hm⟩, ht.1⟩
      (uniqAt_mem_range (show locate d a m p - 1 < nUniq a m by omega)).2 hc
    rw [htp, ← ratioBreakRaw_eq hd hs hm hlt, ← ratioBreak_eq hd hs hm hlt (by omega)] at h
    linarith
  · by_contra hc
    push Not at hc
    have h := LPMratio_lt_LPMratio hd (x := arrSample a m) (uniqAt_mem_range hiN).1
      ⟨⟨m - 1, by omega⟩, ht.2⟩ hc
    rw [htp, ← ratioBreakRaw_eq hd hs hm hlt, ← ratioBreak_eq hd hs hm hlt hiN] at h
    linarith

end Core

theorem LPMratio_comp_perm (d : ℕ) (t : ℝ) (x : Fin n → ℝ) (σ : Equiv.Perm (Fin n)) :
    LPMratio d t (x ∘ σ) = LPMratio d t x := by
  have hL : ∀ e, LPM e t (x ∘ σ) = LPM e t x := fun e => by
    unfold LPM; rw [← Equiv.sum_comp σ (fun i => lowerDev e t (x i))]; rfl
  have hU : ∀ e, UPM e t (x ∘ σ) = UPM e t x := fun e => by
    unfold UPM; rw [← Equiv.sum_comp σ (fun i => upperDev e t (x i))]; rfl
  have hall : (∀ i, (x ∘ σ) i = t) ↔ ∀ i, x i = t :=
    ⟨fun h i => by simpa using h (σ.symm i), fun h i => h (σ i)⟩
  unfold LPMratio
  simp only [hL, hU, hall]

theorem LPMratio_sub_const (d : ℕ) (t c : ℝ) (x : Fin n → ℝ) :
    LPMratio d t (fun i => x i - c) = LPMratio d (t + c) x := by
  have h := LPMratio_affine d one_pos c t (fun i => x i - c)
  simp only [one_mul, sub_add_cancel] at h
  exact h.symm

theorem LPMratio_shiftedArr (d : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    LPMratio d t (arrSample (shiftedArr x) n) = LPMratio d (t + vshift x) x := by
  have : arrSample (shiftedArr x) n = fun i => (x ∘ Tuple.sort x) i - vshift x := by
    funext i
    simp [arrSample, sortedAt, shiftedArr]
  rw [this, LPMratio_sub_const, LPMratio_comp_perm]

lemma shiftedArr_sorted (x : Fin n → ℝ) : VaRSorted (shiftedArr x) n := fun k l hkl hl => by
  simp only [shiftedArr]
  linarith [sortedAt_sorted x k l hkl hl]

lemma shiftedArr_zero (x : Fin n → ℝ) : shiftedArr x 0 + vshift x = varMin x := by
  simp [shiftedArr, varMin]

lemma shiftedArr_last (x : Fin n → ℝ) : shiftedArr x (n - 1) + vshift x = varMax x := by
  simp [shiftedArr, varMax]

lemma shiftedArr_lt {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j) :
    shiftedArr x 0 < shiftedArr x (n - 1) := by
  obtain ⟨i, _, _⟩ := id hx
  have h := (varMin_lt_varMax_iff (Fin.pos i) x).mpr hx
  rw [← shiftedArr_zero, ← shiftedArr_last] at h
  linarith

theorem intLPMVaR_of_const (hn : 0 < n) {x : Fin n → ℝ} (hx : ¬ ∃ i j, x i ≠ x j)
    (tol : ℝ) (d : ℕ) (p : ℝ) : intLPMVaR tol d p x = varMin x := by
  push Not at hx
  have hc : shiftedArr x 0 = shiftedArr x (n - 1) := by
    obtain ⟨i, hi⟩ := exists_eq_sortedAt x hn
    obtain ⟨j, hj⟩ := exists_eq_sortedAt x (show n - 1 < n by omega)
    simp only [shiftedArr]
    rw [← hi, ← hj, hx i j]
  rw [intLPMVaR, ite_eq_left (Or.inr hc), shiftedArr_zero]

theorem intLPMVaR_eq_core {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j) (tol : ℝ) (d : ℕ)
    (p : ℝ) : intLPMVaR tol d p x = intVaRCore tol d (shiftedArr x) n (clampP p) + vshift x := by
  obtain ⟨i, j, hij⟩ := id hx
  have hn1 : n ≠ 1 := by
    rintro rfl
    exact hij (by rw [Subsingleton.elim i j])
  have hgate : ¬ (n = 1 ∨ shiftedArr x 0 = shiftedArr x (n - 1)) := by
    rintro (h | h)
    · exact hn1 h
    · exact (shiftedArr_lt hx).ne h
  rw [intLPMVaR, ite_eq_right hgate]

theorem intLPMVaR_spec {d : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    {tol : ℝ} (htol : 0 ≤ tol) (p : ℝ) :
    intLPMVaR tol d p x ∈ Set.Icc (varMin x) (varMax x) ∧
      |LPMratio d (intLPMVaR tol d p x) x - clampP p| ≤ tol := by
  obtain ⟨i, _, _⟩ := id hx
  have hcore := intVaRCore_spec hd (shiftedArr_sorted x) (Fin.pos i) (shiftedArr_lt hx) htol
    (clampP_mem p)
  rw [intLPMVaR_eq_core hx, ← shiftedArr_zero, ← shiftedArr_last, ← LPMratio_shiftedArr]
  exact ⟨⟨by linarith [hcore.1.1], by linarith [hcore.1.2]⟩, hcore.2⟩

theorem intLPMVaR_same_interval {d : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    {tol : ℝ} (htol : 0 ≤ tol) {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) {t : ℝ}
    (ht : t ∈ Set.Icc (varMin x) (varMax x)) (htp : LPMratio d t x = p) :
    intLPMVaR tol d p x - vshift x ∈
        Set.Icc (uniqAt (shiftedArr x) n (locate d (shiftedArr x) n p - 1))
          (uniqAt (shiftedArr x) n (locate d (shiftedArr x) n p)) ∧
      t - vshift x ∈
        Set.Icc (uniqAt (shiftedArr x) n (locate d (shiftedArr x) n p - 1))
          (uniqAt (shiftedArr x) n (locate d (shiftedArr x) n p)) := by
  obtain ⟨i, _, _⟩ := id hx
  have hcl : clampP p = p := clampP_of_mem ⟨hp0.le, hp1.le⟩
  have h0 := shiftedArr_zero x
  have h1 := shiftedArr_last x
  have h := intVaRCore_same_interval hd (shiftedArr_sorted x) (Fin.pos i) (shiftedArr_lt hx)
    htol hp0 hp1 (t := t - vshift x) ⟨by linarith [ht.1], by linarith [ht.2]⟩
    (by rw [LPMratio_shiftedArr, sub_add_cancel, htp])
  rw [intLPMVaR_eq_core hx, hcl, add_sub_cancel_right]
  exact h

theorem intLPMVaR_dist_root_le {d : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    {tol : ℝ} (htol : 0 ≤ tol) {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) {t : ℝ}
    (ht : t ∈ Set.Icc (varMin x) (varMax x)) (htp : LPMratio d t x = p) :
    |intLPMVaR tol d p x - t| ≤
      uniqAt (shiftedArr x) n (locate d (shiftedArr x) n p) -
        uniqAt (shiftedArr x) n (locate d (shiftedArr x) n p - 1) := by
  obtain ⟨⟨h1, h2⟩, ⟨h3, h4⟩⟩ := intLPMVaR_same_interval hd hx htol hp0 hp1 ht htp
  rw [abs_le]
  constructor <;> linarith

theorem intLPMVaR_snap_min {d : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    {tol p : ℝ} (hp0 : 0 < p) (hptol : p ≤ tol) (hp1 : p < 1)
    (hbelow : ∀ i, varMin x < x i → p < LPMratio d (x i) x) :
    intLPMVaR tol d p x = varMin x := by
  obtain ⟨i0, _, _⟩ := id hx
  have hm : 0 < n := Fin.pos i0
  set a := shiftedArr x with ha
  have hs : VaRSorted a n := shiftedArr_sorted x
  have hlt : a 0 < a (n - 1) := shiftedArr_lt hx
  have hN := two_le_nUniq hs hm hlt
  rw [intLPMVaR_eq_core hx, clampP_of_mem ⟨hp0.le, hp1.le⟩, intVaRCore,
    ite_eq_right (not_le.mpr hp0), ite_eq_right (not_le.mpr hp1)]
  obtain ⟨hi1, hiN, hlo, -⟩ := locate_spec hd hs hm hlt hp0 hp1
  have hloc : locate d a n p = 1 := by
    by_contra hne
    have hrb : ratioBreak d a n 1 ≤ p :=
      (ratioBreak_sorted hd hs hm hlt 1 _ (by omega) (by omega)).trans hlo
    rw [ratioBreak_eq hd hs hm hlt (by omega), ratioBreakRaw_eq hd hs hm hlt,
      LPMratio_shiftedArr] at hrb
    obtain ⟨r, hr, hra⟩ := exists_eq_uniqAt (a := a) (m := n) (k := 1) (by omega)
    obtain ⟨j, hj⟩ := exists_eq_sortedAt x hr
    have hxj : uniqAt a n 1 + vshift x = x j := by
      rw [← hra, hj, ha, shiftedArr, sub_add_cancel]
    have hmin : varMin x < x j := by
      have h01 := uniqAt_strictMono (a := a) (m := n) zero_lt_one (by omega)
      rw [uniqAt_zero hs hm] at h01
      rw [← shiftedArr_zero, ← hxj]
      linarith
    rw [hxj] at hrb
    exact absurd (hbelow j hmin) (not_lt.mpr hrb)
  rw [hloc]
  have hne : uniqAt a n 0 ≠ uniqAt a n 1 := (uniqAt_strictMono zero_lt_one (by omega)).ne
  have hg : lpmRaw d (uniqAt a n 0) (prefixPower a (kBreak a n 0)) /
      (lpmRaw d (uniqAt a n 0) (prefixPower a (kBreak a n 0)) +
        upmRaw d (uniqAt a n 0) (prefixPower a (kBreak a n 0)) (totalPower a n)) = 0 :=
    ratioBreakRaw_zero hd hs hm hlt
  simp only [solveInterval, Nat.sub_self]
  rw [ite_eq_right hne, hg, zero_sub, abs_neg, abs_of_pos hp0, ite_eq_left hptol,
    uniqAt_zero hs hm, shiftedArr_zero]

private lemma sortedAt_fin2 (k : ℕ) (hk : k < 2) : sortedAt ![(0 : ℝ), 1] k = k := by
  have hmono : Monotone ![(0 : ℝ), 1] := by
    refine Fin.monotone_iff_le_succ.mpr fun i => ?_
    fin_cases i
    simp
  have hs := Tuple.sort_eq_refl_iff_monotone.mpr hmono
  simp only [sortedAt, dite_eq_left hk, hs, Equiv.refl_apply]
  interval_cases k <;> simp

private lemma varMin_fin2 : varMin ![(0 : ℝ), 1] = 0 := by
  rw [varMin, sortedAt_fin2 0 (by norm_num)]; simp

private lemma varMax_fin2 : varMax ![(0 : ℝ), 1] = 1 := by
  rw [varMax, sortedAt_fin2 1 (by norm_num)]; simp

private lemma nonconst_fin2 : ∃ i j, ![(0 : ℝ), 1] i ≠ ![(0 : ℝ), 1] j := ⟨0, 1, by simp⟩

private lemma LPMratio_four_fin2 {t : ℝ} (h0 : 0 < t) (h1 : t < 1) :
    LPMratio 4 t ![(0 : ℝ), 1] = t ^ 4 / (t ^ 4 + (1 - t) ^ 4) := by
  rw [LPMratio_eq_sum_div (by norm_num) ⟨0, by simpa using h0.ne⟩]
  simp [Fin.sum_univ_two, lowerDev, upperDev, h0.le, h1, not_le.mpr h1, not_lt.mpr h0.le]

lemma lt_rTol : (1 : ℝ) / 10 ^ 12 < rTol := by
  rw [rTol]; norm_num

theorem intLPMVaR_snap_example :
    intLPMVaR rTol 4 (1 / 10 ^ 12) ![(0 : ℝ), 1] = 0 ∧
      LPMratio 4 (intLPMVaR rTol 4 (1 / 10 ^ 12) ![(0 : ℝ), 1]) ![(0 : ℝ), 1] = 0 ∧
      (∃! t, t ∈ Set.Icc (0 : ℝ) 1 ∧ LPMratio 4 t ![(0 : ℝ), 1] = 1 / 10 ^ 12) ∧
      ∀ t, LPMratio 4 t ![(0 : ℝ), 1] = 1 / 10 ^ 12 → 1 / 2000 < t ∧ t < 1 / 500 := by
  have hres : intLPMVaR rTol 4 (1 / 10 ^ 12) ![(0 : ℝ), 1] = 0 := by
    rw [intLPMVaR_snap_min (by norm_num) nonconst_fin2 (by norm_num) lt_rTol.le (by norm_num),
      varMin_fin2]
    intro i hi
    rw [varMin_fin2] at hi
    rw [LPMratio_eq_one_of_forall_le (by norm_num) 4 (x := ![(0 : ℝ), 1])]
    · norm_num
    · intro j
      fin_cases i <;> fin_cases j <;> simp_all
  have hzero : LPMratio 4 0 ![(0 : ℝ), 1] = 0 := by
    have := LPMratio_varMin (d := 4) (by norm_num) nonconst_fin2
    rwa [varMin_fin2] at this
  refine ⟨hres, by rw [hres, hzero], ?_, ?_⟩
  · have := existsUnique_LPMratio_eq (d := 4) (by norm_num) nonconst_fin2
      (p := 1 / 10 ^ 12) ⟨by norm_num, by norm_num⟩
    rwa [varMin_fin2, varMax_fin2] at this
  · intro t ht
    have hmono := LPMratio_monotone 4 ![(0 : ℝ), 1]
    have hlo : LPMratio 4 (1 / 2000) ![(0 : ℝ), 1] < 1 / 10 ^ 12 := by
      rw [LPMratio_four_fin2 (by norm_num) (by norm_num)]; norm_num
    have hhi : 1 / 10 ^ 12 < LPMratio 4 (1 / 500) ![(0 : ℝ), 1] := by
      rw [LPMratio_four_fin2 (by norm_num) (by norm_num)]; norm_num
    constructor
    · by_contra hc
      have := hmono (not_lt.mp hc)
      simp only at this
      linarith
    · by_contra hc
      have := hmono (not_lt.mp hc)
      simp only at this
      linarith

open Classical in
def fallbackLPMVaR (d : ℕ) (p : ℝ) (x : Fin n → ℝ) : ℝ :=
  let q := clampP p
  let g : ℝ → ℝ := fun b => LPMratio d b x - q
  if varMin x = varMax x then varMin x
  else if g (varMin x) = 0 then varMin x
  else if g (varMax x) = 0 then varMax x
  else if Real.sign (g (varMin x)) = Real.sign (g (varMax x)) then
    (if |g (varMin x)| ≤ |g (varMax x)| then varMin x else varMax x)
  else unirootIdeal g (varMin x) (varMax x)

open Classical in
def fallbackUPMVaR (d : ℕ) (p : ℝ) (x : Fin n → ℝ) : ℝ :=
  let q := clampP p
  let g : ℝ → ℝ := fun b => UPMratio d b x - q
  if varMin x = varMax x then varMin x
  else if g (varMin x) = 0 then varMin x
  else if g (varMax x) = 0 then varMax x
  else if Real.sign (g (varMin x)) = Real.sign (g (varMax x)) then
    (if |g (varMin x)| ≤ |g (varMax x)| then varMin x else varMax x)
  else unirootIdeal g (varMin x) (varMax x)

theorem fallbackLPMVaR_spec {d : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    (p : ℝ) : fallbackLPMVaR d p x ∈ Set.Icc (varMin x) (varMax x) ∧
      LPMratio d (fallbackLPMVaR d p x) x = clampP p := by
  obtain ⟨i, _, _⟩ := id hx
  have hn : 0 < n := Fin.pos i
  have hlt := (varMin_lt_varMax_iff hn x).mpr hx
  have hq := clampP_mem p
  have hlo := LPMratio_varMin hd hx
  have hhi := LPMratio_varMax hn d x
  simp only [fallbackLPMVaR]
  split_ifs with h1 h2 h3 h4
  · exact absurd h1 hlt.ne
  · exact ⟨⟨le_rfl, hlt.le⟩, sub_eq_zero.mp h2⟩
  · exact ⟨⟨hlt.le, le_rfl⟩, sub_eq_zero.mp h3⟩
  · rw [hlo, zero_sub] at h2 h4
    rw [hhi] at h3 h4
    have hq0 : 0 < clampP p := lt_of_le_of_ne hq.1 (fun e => h2 (by rw [← e, neg_zero]))
    have hq1 : clampP p < 1 := lt_of_le_of_ne hq.2 (fun e => h3 (by rw [e, sub_self]))
    rw [Real.sign_of_neg (by linarith), Real.sign_of_pos (by linarith)] at h4
    norm_num at h4
  · rw [hlo, zero_sub] at h2 h4
    rw [hhi] at h3 h4
    have hq0 : 0 < clampP p := lt_of_le_of_ne hq.1 (fun e => h2 (by rw [← e, neg_zero]))
    have hq1 : clampP p < 1 := lt_of_le_of_ne hq.2 (fun e => h3 (by rw [e, sub_self]))
    rw [Real.sign_of_neg (by linarith), Real.sign_of_pos (by linarith)] at h4
    norm_num at h4
  · have hcont : Continuous fun b => LPMratio d b x - clampP p :=
      (LPMratio_continuous hd hx).sub continuous_const
    have hex : ∃ t ∈ Set.Icc (varMin x) (varMax x), LPMratio d t x - clampP p = 0 :=
      intermediate_value_Icc hlt.le hcont.continuousOn
        ⟨by rw [hlo]; linarith [hq.1], by rw [hhi]; linarith [hq.2]⟩
    rw [unirootIdeal, dite_eq_left hex]
    exact ⟨hex.choose_spec.1, sub_eq_zero.mp hex.choose_spec.2⟩

theorem fallbackUPMVaR_spec {d : ℕ} (hd : d ≠ 0) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    (p : ℝ) : fallbackUPMVaR d p x ∈ Set.Icc (varMin x) (varMax x) ∧
      UPMratio d (fallbackUPMVaR d p x) x = clampP p := by
  obtain ⟨i, _, _⟩ := id hx
  have hn : 0 < n := Fin.pos i
  have hlt := (varMin_lt_varMax_iff hn x).mpr hx
  have hq := clampP_mem p
  have hlo : UPMratio d (varMin x) x = 1 := by
    rw [UPMratio_eq_one_sub_LPMratio hn, LPMratio_varMin hd hx, sub_zero]
  have hhi : UPMratio d (varMax x) x = 0 := by
    rw [UPMratio_eq_one_sub_LPMratio hn, LPMratio_varMax hn, sub_self]
  simp only [fallbackUPMVaR]
  split_ifs with h1 h2 h3 h4
  · exact absurd h1 hlt.ne
  · exact ⟨⟨le_rfl, hlt.le⟩, sub_eq_zero.mp h2⟩
  · exact ⟨⟨hlt.le, le_rfl⟩, sub_eq_zero.mp h3⟩
  · rw [hlo] at h2 h4
    rw [hhi, zero_sub] at h3 h4
    have hq1 : clampP p < 1 := lt_of_le_of_ne hq.2 (fun e => h2 (by rw [e, sub_self]))
    have hq0 : 0 < clampP p := lt_of_le_of_ne hq.1 (fun e => h3 (by rw [← e, neg_zero]))
    rw [Real.sign_of_pos (by linarith), Real.sign_of_neg (by linarith)] at h4
    norm_num at h4
  · rw [hlo] at h2 h4
    rw [hhi, zero_sub] at h3 h4
    have hq1 : clampP p < 1 := lt_of_le_of_ne hq.2 (fun e => h2 (by rw [e, sub_self]))
    have hq0 : 0 < clampP p := lt_of_le_of_ne hq.1 (fun e => h3 (by rw [← e, neg_zero]))
    rw [Real.sign_of_pos (by linarith), Real.sign_of_neg (by linarith)] at h4
    norm_num at h4
  · have hcont : Continuous fun b => UPMratio d b x - clampP p := by
      have : (fun b => UPMratio d b x - clampP p) = fun b => 1 - LPMratio d b x - clampP p :=
        funext fun b => by rw [UPMratio_eq_one_sub_LPMratio hn]
      rw [this]
      exact (continuous_const.sub (LPMratio_continuous hd hx)).sub continuous_const
    have hex : ∃ t ∈ Set.Icc (varMin x) (varMax x), UPMratio d t x - clampP p = 0 :=
      intermediate_value_Icc' hlt.le hcont.continuousOn
        ⟨by rw [hhi]; linarith [hq.1], by rw [hlo]; linarith [hq.2]⟩
    rw [unirootIdeal, dite_eq_left hex]
    exact ⟨hex.choose_spec.1, sub_eq_zero.mp hex.choose_spec.2⟩

theorem fallbackUPMVaR_eq {d : ℕ} (hd : d ≠ 0) (hn : 0 < n) (x : Fin n → ℝ) (p : ℝ) :
    fallbackUPMVaR d p x = fallbackLPMVaR d (1 - clampP p) x := by
  by_cases hx : ∃ i j, x i ≠ x j
  · obtain ⟨hU, hUr⟩ := fallbackUPMVaR_spec hd hx p
    obtain ⟨hL, hLr⟩ := fallbackLPMVaR_spec hd hx (1 - clampP p)
    have hq := clampP_mem p
    rw [clampP_of_mem (p := 1 - clampP p) ⟨by linarith [hq.2], by linarith [hq.1]⟩] at hLr
    refine eq_of_LPMratio_eq hd hn hU hL ?_
    rw [hLr, ← hUr, UPMratio_eq_one_sub_LPMratio hn]
    ring
  · have hc : varMin x = varMax x := by
      by_contra h
      exact hx ((varMin_lt_varMax_iff hn x).mp
        (lt_of_le_of_ne (varMin_le x ⟨0, hn⟩ |>.trans (le_varMax x ⟨0, hn⟩)) h))
    simp only [fallbackUPMVaR, fallbackLPMVaR, hc, ↓reduceIte]

def LPMVaR (d : ℕ) (p : ℝ) (x : Fin n → ℝ) : ℝ :=
  if d = 0 then LPMVaR0 p x
  else if d ≤ 4 then intLPMVaR rTol d (clampP p) x
  else fallbackLPMVaR d (clampP p) x

def UPMVaR (d : ℕ) (p : ℝ) (x : Fin n → ℝ) : ℝ :=
  if d = 0 then UPMVaR0 p x
  else if d ≤ 4 then intUPMVaR rTol d (clampP p) x
  else fallbackUPMVaR d (clampP p) x

lemma clampP_clampP (p : ℝ) : clampP (clampP p) = clampP p := clampP_of_mem (clampP_mem p)

lemma clampP_one_sub_clampP (p : ℝ) : clampP (1 - clampP p) = 1 - clampP p := by
  have h := clampP_mem p
  exact clampP_of_mem ⟨by linarith [h.2], by linarith [h.1]⟩

theorem UPMVaR_eq_LPMVaR (hn : 0 < n) (d : ℕ) (p : ℝ) (x : Fin n → ℝ) :
    UPMVaR d p x = LPMVaR d (1 - clampP p) x := by
  unfold UPMVaR LPMVaR
  split_ifs with h0 h4
  · exact UPMVaR0_eq_LPMVaR0 p x
  · rw [intUPMVaR, clampP_clampP, clampP_one_sub_clampP]
  · rw [fallbackUPMVaR_eq h0 hn, clampP_clampP, clampP_one_sub_clampP]

theorem LPMVaR_of_const (hn : 0 < n) {x : Fin n → ℝ} (hx : ¬ ∃ i j, x i ≠ x j) {d : ℕ}
    (hd : d ≠ 0) (p : ℝ) : LPMVaR d p x = varMin x := by
  have hc : varMin x = varMax x := by
    by_contra h
    exact hx ((varMin_lt_varMax_iff hn x).mp
      (lt_of_le_of_ne ((varMin_le x ⟨0, hn⟩).trans (le_varMax x ⟨0, hn⟩)) h))
  unfold LPMVaR
  split_ifs with h0 h4
  · exact absurd h0 hd
  · exact intLPMVaR_of_const hn hx rTol d _
  · simp only [fallbackLPMVaR, hc, ↓reduceIte]

theorem LPMVaR_spec_of_five_le {d : ℕ} (hd : 5 ≤ d) {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    {p : ℝ} (hp : p ∈ Set.Icc (0 : ℝ) 1) :
    LPMVaR d p x ∈ Set.Icc (varMin x) (varMax x) ∧ LPMratio d (LPMVaR d p x) x = p := by
  have h := fallbackLPMVaR_spec (d := d) (by omega) hx p
  rw [clampP_of_mem hp] at h
  rw [LPMVaR, ite_eq_right (by omega), ite_eq_right (by omega), clampP_of_mem hp]
  exact h

theorem LPMVaR_spec_of_le_four {d : ℕ} (hd1 : 1 ≤ d) (hd4 : d ≤ 4) {x : Fin n → ℝ}
    (hx : ∃ i j, x i ≠ x j) (p : ℝ) :
    LPMVaR d p x ∈ Set.Icc (varMin x) (varMax x) ∧
      |LPMratio d (LPMVaR d p x) x - clampP p| ≤ rTol := by
  have h := intLPMVaR_spec (d := d) (tol := rTol) (by omega) hx (by norm_num [rTol]) (clampP p)
  rw [clampP_clampP] at h
  rw [LPMVaR, ite_eq_right (by omega), ite_eq_left hd4]
  exact h

theorem LPMVaR_dist_root_le_of_le_four {d : ℕ} (hd1 : 1 ≤ d) (hd4 : d ≤ 4) {x : Fin n → ℝ}
    (hx : ∃ i j, x i ≠ x j) {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) {t : ℝ}
    (ht : t ∈ Set.Icc (varMin x) (varMax x)) (htp : LPMratio d t x = p) :
    |LPMVaR d p x - t| ≤
      uniqAt (shiftedArr x) n (locate d (shiftedArr x) n p) -
        uniqAt (shiftedArr x) n (locate d (shiftedArr x) n p - 1) := by
  rw [LPMVaR, ite_eq_right (by omega), ite_eq_left hd4, clampP_of_mem ⟨hp0.le, hp1.le⟩]
  exact intLPMVaR_dist_root_le (by omega) hx (by norm_num [rTol]) hp0 hp1 ht htp

theorem UPMVaR_spec_of_le_four {d : ℕ} (hd1 : 1 ≤ d) (hd4 : d ≤ 4) {x : Fin n → ℝ}
    (hx : ∃ i j, x i ≠ x j) (p : ℝ) :
    UPMVaR d p x ∈ Set.Icc (varMin x) (varMax x) ∧
      |UPMratio d (UPMVaR d p x) x - clampP p| ≤ rTol := by
  obtain ⟨i, _, _⟩ := id hx
  have hn : 0 < n := Fin.pos i
  rw [UPMVaR_eq_LPMVaR hn]
  obtain ⟨hmem, habs⟩ := LPMVaR_spec_of_le_four hd1 hd4 hx (1 - clampP p)
  rw [clampP_one_sub_clampP] at habs
  refine ⟨hmem, ?_⟩
  rw [UPMratio_eq_one_sub_LPMratio hn, abs_sub_comm]
  convert habs using 2
  ring

lemma varMin_neg (hn : 0 < n) (x : Fin n → ℝ) :
    varMin (fun i => -x i) = -varMax x := by
  obtain ⟨i, hi⟩ := exists_eq_varMin hn (fun i => -x i)
  obtain ⟨j, hj⟩ := exists_eq_varMax hn x
  refine le_antisymm ?_ ?_
  · rw [← hj]; exact varMin_le (fun i => -x i) j
  · rw [← hi]; exact neg_le_neg (le_varMax x i)

lemma varMax_neg (hn : 0 < n) (x : Fin n → ℝ) :
    varMax (fun i => -x i) = -varMin x := by
  obtain ⟨i, hi⟩ := exists_eq_varMax hn (fun i => -x i)
  obtain ⟨j, hj⟩ := exists_eq_varMin hn x
  refine le_antisymm ?_ ?_
  · rw [← hi]; exact neg_le_neg (varMin_le x i)
  · rw [← hj]; exact le_varMax (fun i => -x i) j

end

end NNS
