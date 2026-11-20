-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.VaRRoot
public import NNSProofs.MeasurePM

namespace NNS.GapCertainty

open Finset NNS MeasureTheory ProbabilityTheory Set

@[expose] public noncomputable section

variable {n m : ℕ}

def constVal (x : Fin n → ℝ) : Option ℝ :=
  if h : 0 < n then (if ∀ i, x i = x ⟨0, h⟩ then some (x ⟨0, h⟩) else none) else none

def pooledMean (x : Fin n → ℝ) (y : Fin m → ℝ) : ℝ :=
  ((n : ℝ) * mean x + m * mean y) / (n + m)

def rMedian (z : Fin n → ℝ) : ℝ :=
  if n % 2 = 1 then sortedAt z (n / 2) else (sortedAt z (n / 2 - 1) + sortedAt z (n / 2)) / 2

def grandStat (x : Fin n → ℝ) (y : Fin m → ℝ) (medians : Bool) : ℝ :=
  if medians then rMedian (Fin.append x y)
  else if (constVal x).isSome ∧ constVal x = constVal y then (constVal x).getD 0
  else pooledMean x y

def areaDist (cx cy : Option ℝ) (r1 r2 target p : ℝ) : ℝ :=
  min p (max (if cx = some target then 0 else |r1 - p|) (if cy = some target then 0 else |r2 - p|))

def ltFrac (t : ℝ) (x : Fin n → ℝ) : ℝ := (∑ i, if x i < t then (1 : ℝ) else 0) / n

def medianDist (t : ℝ) (x : Fin n → ℝ) : ℝ :=
  max (max 0 (ltFrac t x - 1 / 2)) (1 / 2 - LPMratio 0 t x)

structure Targets where
  u25 : ℝ
  l25 : ℝ
  u125 : ℝ
  l125 : ℝ

def defaultTargets (x : Fin n → ℝ) (y : Fin m → ℝ) : Targets where
  u25 := (UPMVaR 1 (1 / 4) x + UPMVaR 1 (1 / 4) y) / 2
  l25 := (LPMVaR 1 (1 / 4) x + LPMVaR 1 (1 / 4) y) / 2
  u125 := (UPMVaR 1 (1 / 8) x + UPMVaR 1 (1 / 8) y) / 2
  l125 := (LPMVaR 1 (1 / 8) x + LPMVaR 1 (1 / 8) y) / 2

def centralDist (x : Fin n → ℝ) (y : Fin m → ℝ) (medians : Bool) (z : ℝ) : ℝ :=
  if medians then max (medianDist z x) (medianDist z y)
  else areaDist (constVal x) (constVal y) (LPMratio 1 z x) (LPMratio 1 z y) z (1 / 2)

def dU25 (x : Fin n → ℝ) (y : Fin m → ℝ) (T : Targets) : ℝ :=
  areaDist (constVal x) (constVal y) (UPMratio 1 T.u25 x) (UPMratio 1 T.u25 y) T.u25 (1 / 4)

def dL25 (x : Fin n → ℝ) (y : Fin m → ℝ) (T : Targets) : ℝ :=
  areaDist (constVal x) (constVal y) (LPMratio 1 T.l25 x) (LPMratio 1 T.l25 y) T.l25 (1 / 4)

def dU125 (x : Fin n → ℝ) (y : Fin m → ℝ) (T : Targets) : ℝ :=
  areaDist (constVal x) (constVal y) (UPMratio 1 T.u125 x) (UPMratio 1 T.u125 y) T.u125 (1 / 8)

def dL125 (x : Fin n → ℝ) (y : Fin m → ℝ) (T : Targets) : ℝ :=
  areaDist (constVal x) (constVal y) (LPMratio 1 T.l125 x) (LPMratio 1 T.l125 y) T.l125 (1 / 8)

def rawScore (meansOnly : Bool) (a b c d e : ℝ) : ℝ :=
  if meansOnly then (1 / 2 - a) ^ 2 / (1 / 4)
  else ((1 / 2 - a) ^ 2 / (1 / 4) + 1 / 2 * ((1 / 4 - b) ^ 2 / (1 / 4) ^ 2) +
    1 / 2 * ((1 / 4 - c) ^ 2 / (1 / 4) ^ 2) + 1 / 4 * ((1 / 8 - d) ^ 2 / (1 / 8) ^ 2) +
    1 / 4 * ((1 / 8 - e) ^ 2 / (1 / 8) ^ 2)) / (5 / 2)

def popAdj (n m : ℕ) : ℝ := (((n : ℝ) + m - 2) / (n + m)) ^ 2

def certaintyAt (x : Fin n → ℝ) (y : Fin m → ℝ) (meansOnly medians : Bool) (z : ℝ)
    (T : Targets) : ℝ :=
  min 1 (rawScore meansOnly (centralDist x y medians z) (dU25 x y T) (dL25 x y T)
    (dU125 x y T) (dL125 x y T) * popAdj n m)

def certainty (x : Fin n → ℝ) (y : Fin m → ℝ) (meansOnly medians : Bool) : ℝ :=
  certaintyAt x y meansOnly medians (grandStat x y medians) (defaultTargets x y)

def controlCDF (x : Fin n → ℝ) (medians : Bool) (z : ℝ) : ℝ :=
  if medians then LPMratio 0 z x else LPMratio 1 z x

lemma constVal_eq_some {x : Fin n → ℝ} {c : ℝ} :
    constVal x = some c ↔ 0 < n ∧ ∀ i, x i = c := by
  unfold constVal
  constructor
  · intro h
    split_ifs at h with hn hall
    cases h
    exact ⟨hn, hall⟩
  · rintro ⟨hn, hall⟩
    simp [hn, hall]

lemma mean_const {x : Fin n → ℝ} {c : ℝ} (hn : 0 < n) (h : ∀ i, x i = c) : mean x = c := by
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  simp [mean, h, hn']

lemma mul_mean (x : Fin n → ℝ) : (n : ℝ) * mean x = ∑ i, x i := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp
  · rw [mean, mul_div_cancel₀ _ (Nat.cast_ne_zero.mpr hn.ne')]

theorem pooledMean_eq_mean_append (x : Fin n → ℝ) (y : Fin m → ℝ) :
    pooledMean x y = mean (Fin.append x y) := by
  rw [pooledMean, mul_mean, mul_mean, mean, Fin.sum_univ_add]
  simp [Fin.append_left, Fin.append_right, Nat.cast_add]

theorem grandStat_eq_mean_append (x : Fin n → ℝ) (y : Fin m → ℝ) :
    grandStat x y false = mean (Fin.append x y) := by
  simp only [grandStat, Bool.false_eq_true, ↓reduceIte]
  split_ifs with hc
  · obtain ⟨hs, he⟩ := hc
    obtain ⟨c, hcx⟩ := Option.isSome_iff_exists.mp hs
    obtain ⟨hn, hx⟩ := constVal_eq_some.mp hcx
    obtain ⟨hm, hy⟩ := constVal_eq_some.mp (he ▸ hcx)
    rw [hcx, Option.getD_some, mean_const (by omega)]
    intro i
    refine Fin.addCases (fun i => ?_) (fun j => ?_) i
    · simp [hx]
    · simp [hy]
  · exact pooledMean_eq_mean_append x y

theorem areaDist_mem_Icc (cx cy : Option ℝ) (r1 r2 t : ℝ) {p : ℝ} (hp : 0 ≤ p) :
    areaDist cx cy r1 r2 t p ∈ Icc 0 p := by
  refine ⟨le_min hp (le_max_of_le_left ?_), min_le_left _ _⟩
  split_ifs <;> simp

theorem areaDist_eq_book {cx cy : Option ℝ} {t : ℝ} (hx : cx ≠ some t) (hy : cy ≠ some t)
    (r1 r2 p : ℝ) : areaDist cx cy r1 r2 t p = min p (max |r1 - p| |r2 - p|) := by
  simp [areaDist, hx, hy]

lemma abs_sub_half_le {r : ℝ} (hr : r ∈ Icc (0 : ℝ) 1) : |r - 1 / 2| ≤ 1 / 2 := by
  rw [abs_le]; constructor <;> linarith [hr.1, hr.2]

theorem centralDist_eq_book (x : Fin n → ℝ) (y : Fin m → ℝ) {z : ℝ}
    (hx : constVal x ≠ some z) (hy : constVal y ≠ some z) :
    centralDist x y false z =
        min (1 / 2) (max |LPMratio 1 z x - 1 / 2| |LPMratio 1 z y - 1 / 2|) ∧
      centralDist x y false z = max |LPMratio 1 z x - 1 / 2| |LPMratio 1 z y - 1 / 2| := by
  have hb := areaDist_eq_book hx hy (LPMratio 1 z x) (LPMratio 1 z y) (1 / 2)
  simp only [centralDist, Bool.false_eq_true, ↓reduceIte]
  refine ⟨hb, hb.trans (min_eq_right (max_le ?_ ?_))⟩
  · exact abs_sub_half_le (LPMratio_mem_Icc 1 z x)
  · exact abs_sub_half_le (LPMratio_mem_Icc 1 z y)

theorem centralDist_const_witness :
    centralDist ![(3 : ℝ), 3] ![(3 : ℝ), 3] false (grandStat ![(3 : ℝ), 3] ![(3 : ℝ), 3] false)
        = 0 ∧
      min (1 / 2) (max |LPMratio 1 3 ![(3 : ℝ), 3] - 1 / 2| |LPMratio 1 3 ![(3 : ℝ), 3] - 1 / 2|)
        = 1 / 2 := by
  have hc : constVal ![(3 : ℝ), 3] = some 3 :=
    constVal_eq_some.mpr ⟨by norm_num, fun i => by fin_cases i <;> rfl⟩
  have hg : grandStat ![(3 : ℝ), 3] ![(3 : ℝ), 3] false = 3 := by
    simp [grandStat, hc]
  have hr : LPMratio 1 3 ![(3 : ℝ), 3] = 1 := by
    have h3 : 0 < 2 ∧ ∀ i, ![(3 : ℝ), 3] i = 3 := ⟨by norm_num, fun i => by fin_cases i <;> rfl⟩
    simp only [LPMratio, one_ne_zero, h3, ↓reduceIte]
    simp
  refine ⟨?_, ?_⟩
  · rw [hg]; simp [centralDist, areaDist, hc]
  · rw [hr]; norm_num

theorem UPMVaR_one_eq (hn : 0 < n) (x : Fin n → ℝ) {p : ℝ} (hp : p ∈ Icc (0 : ℝ) 1) :
    UPMVaR 1 p x = LPMVaR 1 (1 - p) x := by
  rw [UPMVaR_eq_LPMVaR hn, clampP_of_mem hp]

theorem LPMratio_one_eq_iff_balance {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j)
    (t p : ℝ) : LPMratio 1 t x = p ↔ (1 - p) * LPM 1 t x = p * UPM 1 t x := by
  obtain ⟨i, _, _⟩ := id hx
  have hn' : (0 : ℝ) < n := Nat.cast_pos.mpr (Fin.pos i)
  rw [LPMratio_one_eq_iff_expectile hx, LPM, UPM, mul_div_assoc', mul_div_assoc',
    div_left_inj' hn'.ne']

theorem UPMratio_eq_iff_LPMratio (hn : 0 < n) (x : Fin n → ℝ) (t p : ℝ) :
    UPMratio 1 t x = p ↔ LPMratio 1 t x = 1 - p := by
  rw [UPMratio_eq_one_sub_LPMratio hn]
  constructor <;> intro h <;> linarith

theorem lowerTarget_ratio {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j) {p : ℝ}
    (hp : p ∈ Icc (0 : ℝ) 1) : |LPMratio 1 (LPMVaR 1 p x) x - p| ≤ rTol := by
  have h := (LPMVaR_spec_of_le_four le_rfl (by norm_num) hx p).2
  rwa [clampP_of_mem hp] at h

theorem upperTarget_ratio {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j) {p : ℝ}
    (hp : p ∈ Icc (0 : ℝ) 1) : |UPMratio 1 (UPMVaR 1 p x) x - p| ≤ rTol := by
  have h := (UPMVaR_spec_of_le_four le_rfl (by norm_num) hx p).2
  rwa [clampP_of_mem hp] at h

theorem rawScore_eq_book (a b c d e : ℝ) :
    rawScore false a b c d e =
      1 / (25 / 10) * ((1 / 2 - a) ^ 2 / (25 / 100) +
        5 / 10 * ((25 / 100 - b) ^ 2 / (625 / 10000)) +
        5 / 10 * ((25 / 100 - c) ^ 2 / (625 / 10000)) +
        25 / 100 * ((125 / 1000 - d) ^ 2 / (15625 / 1000000)) +
        25 / 100 * ((125 / 1000 - e) ^ 2 / (15625 / 1000000))) := by
  simp only [rawScore, Bool.false_eq_true, ↓reduceIte]
  ring

theorem rawScore_zero (meansOnly : Bool) : rawScore meansOnly 0 0 0 0 0 = 1 ∧
    (1 : ℝ) + 1 / 2 + 1 / 2 + 1 / 4 + 1 / 4 = 5 / 2 := by
  cases meansOnly <;> norm_num [rawScore]

theorem rawScore_meansOnly (a b c d e : ℝ) :
    rawScore true a b c d e = ((1 / 2 - a) ^ 2 / (1 / 4)) / 1 := by
  simp [rawScore]

lemma sum_ind_le (x : Fin n → ℝ) (P : ℝ → Prop) [DecidablePred P] :
    ∑ i, (if P (x i) then (1 : ℝ) else 0) ≤ n := by
  calc ∑ i, (if P (x i) then (1 : ℝ) else 0) ≤ ∑ _i : Fin n, (1 : ℝ) :=
        sum_le_sum fun i _ => by split_ifs <;> norm_num
    _ = n := by simp

lemma sum_ind_nonneg (x : Fin n → ℝ) (P : ℝ → Prop) [DecidablePred P] :
    0 ≤ ∑ i, (if P (x i) then (1 : ℝ) else 0) :=
  sum_nonneg fun i _ => by split_ifs <;> norm_num

lemma ltFrac_mem_Icc (t : ℝ) (x : Fin n → ℝ) : ltFrac t x ∈ Icc (0 : ℝ) 1 := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [ltFrac]
  have hn' : (0 : ℝ) < n := Nat.cast_pos.mpr hn
  exact ⟨div_nonneg (sum_ind_nonneg x (· < t)) hn'.le, (div_le_one hn').mpr (sum_ind_le x (· < t))⟩

lemma medianDist_mem_Icc (t : ℝ) (x : Fin n → ℝ) : medianDist t x ∈ Icc (0 : ℝ) (1 / 2) := by
  have h1 := ltFrac_mem_Icc t x
  have h2 := LPMratio_mem_Icc 0 t x
  refine ⟨le_max_of_le_left (le_max_left _ _), max_le (max_le (by norm_num) ?_) ?_⟩ <;>
    linarith [h1.2, h2.1]

theorem centralDist_mem_Icc (x : Fin n → ℝ) (y : Fin m → ℝ) (medians : Bool) (z : ℝ) :
    centralDist x y medians z ∈ Icc (0 : ℝ) (1 / 2) := by
  unfold centralDist
  split_ifs
  · exact ⟨le_max_of_le_left (medianDist_mem_Icc z x).1,
      max_le (medianDist_mem_Icc z x).2 (medianDist_mem_Icc z y).2⟩
  · exact areaDist_mem_Icc _ _ _ _ _ (by norm_num)

lemma cap_le {δ p q : ℝ} (hq : q * p ^ 2 = 1) (hq0 : 0 ≤ q) (h0 : 0 ≤ δ) (h1 : δ ≤ p) :
    q * (p - δ) ^ 2 ≤ 1 := by
  have : (p - δ) ^ 2 ≤ p ^ 2 := by nlinarith
  nlinarith

lemma cap_eq {δ p q : ℝ} (hq : q * p ^ 2 = 1) (h0 : 0 ≤ δ) (h1 : δ ≤ p)
    (h : q * (p - δ) ^ 2 = 1) : δ = 0 := by
  have hq' : q ≠ 0 := by rintro rfl; simp at hq
  have : (p - δ) ^ 2 = p ^ 2 := mul_left_cancel₀ hq' (h.trans hq.symm)
  nlinarith

lemma rawScore_poly (meansOnly : Bool) (a b c d e : ℝ) :
    rawScore meansOnly a b c d e =
      if meansOnly then 4 * (1 / 2 - a) ^ 2
      else 2 / 5 * (4 * (1 / 2 - a) ^ 2 + 1 / 2 * (16 * (1 / 4 - b) ^ 2) +
        1 / 2 * (16 * (1 / 4 - c) ^ 2) + 1 / 4 * (64 * (1 / 8 - d) ^ 2) +
        1 / 4 * (64 * (1 / 8 - e) ^ 2)) := by
  cases meansOnly <;> simp only [rawScore, Bool.false_eq_true, ↓reduceIte] <;> ring

theorem rawScore_mem_Icc (meansOnly : Bool) {a b c d e : ℝ} (ha : a ∈ Icc (0 : ℝ) (1 / 2))
    (hb : b ∈ Icc (0 : ℝ) (1 / 4)) (hc : c ∈ Icc (0 : ℝ) (1 / 4))
    (hd : d ∈ Icc (0 : ℝ) (1 / 8)) (he : e ∈ Icc (0 : ℝ) (1 / 8)) :
    rawScore meansOnly a b c d e ∈ Icc (0 : ℝ) 1 := by
  obtain ⟨ha0, ha1⟩ := ha; obtain ⟨hb0, hb1⟩ := hb; obtain ⟨hc0, hc1⟩ := hc
  obtain ⟨hd0, hd1⟩ := hd; obtain ⟨he0, he1⟩ := he
  have ta : 4 * (1 / 2 - a) ^ 2 ∈ Icc (0 : ℝ) 1 :=
    ⟨by positivity, cap_le (by norm_num) (by norm_num) ha0 ha1⟩
  have tb : 16 * (1 / 4 - b) ^ 2 ∈ Icc (0 : ℝ) 1 :=
    ⟨by positivity, cap_le (by norm_num) (by norm_num) hb0 hb1⟩
  have tc : 16 * (1 / 4 - c) ^ 2 ∈ Icc (0 : ℝ) 1 :=
    ⟨by positivity, cap_le (by norm_num) (by norm_num) hc0 hc1⟩
  have td : 64 * (1 / 8 - d) ^ 2 ∈ Icc (0 : ℝ) 1 :=
    ⟨by positivity, cap_le (by norm_num) (by norm_num) hd0 hd1⟩
  have te : 64 * (1 / 8 - e) ^ 2 ∈ Icc (0 : ℝ) 1 :=
    ⟨by positivity, cap_le (by norm_num) (by norm_num) he0 he1⟩
  rw [rawScore_poly]
  split_ifs
  · exact ta
  · constructor <;> linarith [ta.1, ta.2, tb.1, tb.2, tc.1, tc.2, td.1, td.2, te.1, te.2]

theorem rawScore_eq_one_iff (meansOnly : Bool) {a b c d e : ℝ}
    (ha : a ∈ Icc (0 : ℝ) (1 / 2)) (hb : b ∈ Icc (0 : ℝ) (1 / 4))
    (hc : c ∈ Icc (0 : ℝ) (1 / 4)) (hd : d ∈ Icc (0 : ℝ) (1 / 8))
    (he : e ∈ Icc (0 : ℝ) (1 / 8)) :
    rawScore meansOnly a b c d e = 1 ↔
      a = 0 ∧ (meansOnly = true ∨ (b = 0 ∧ c = 0 ∧ d = 0 ∧ e = 0)) := by
  obtain ⟨ha0, ha1⟩ := ha; obtain ⟨hb0, hb1⟩ := hb; obtain ⟨hc0, hc1⟩ := hc
  obtain ⟨hd0, hd1⟩ := hd; obtain ⟨he0, he1⟩ := he
  have ta : 4 * (1 / 2 - a) ^ 2 ≤ 1 := cap_le (by norm_num) (by norm_num) ha0 ha1
  have tb : 16 * (1 / 4 - b) ^ 2 ≤ 1 := cap_le (by norm_num) (by norm_num) hb0 hb1
  have tc : 16 * (1 / 4 - c) ^ 2 ≤ 1 := cap_le (by norm_num) (by norm_num) hc0 hc1
  have td : 64 * (1 / 8 - d) ^ 2 ≤ 1 := cap_le (by norm_num) (by norm_num) hd0 hd1
  have te : 64 * (1 / 8 - e) ^ 2 ≤ 1 := cap_le (by norm_num) (by norm_num) he0 he1
  rw [rawScore_poly]
  cases meansOnly
  · simp only [Bool.false_eq_true, ↓reduceIte, false_or]
    constructor
    · intro h
      have ea : 4 * (1 / 2 - a) ^ 2 = 1 := by linarith
      have eb : 16 * (1 / 4 - b) ^ 2 = 1 := by linarith
      have ec : 16 * (1 / 4 - c) ^ 2 = 1 := by linarith
      have ed : 64 * (1 / 8 - d) ^ 2 = 1 := by linarith
      have ee : 64 * (1 / 8 - e) ^ 2 = 1 := by linarith
      exact ⟨cap_eq (by norm_num) ha0 ha1 ea, cap_eq (by norm_num) hb0 hb1 eb,
        cap_eq (by norm_num) hc0 hc1 ec, cap_eq (by norm_num) hd0 hd1 ed,
        cap_eq (by norm_num) he0 he1 ee⟩
    · rintro ⟨rfl, rfl, rfl, rfl, rfl⟩; norm_num
  · simp only [↓reduceIte, true_or, and_true]
    constructor
    · exact cap_eq (by norm_num) ha0 ha1
    · rintro rfl; norm_num

theorem popAdj_mem (hn : 0 < n) (hm : 0 < m) : 0 ≤ popAdj n m ∧ popAdj n m < 1 := by
  have h1 : (1 : ℝ) ≤ n := Nat.one_le_cast.mpr hn
  have h2 : (1 : ℝ) ≤ m := Nat.one_le_cast.mpr hm
  have hN : (0 : ℝ) < n + m := by linarith
  have hq0 : 0 ≤ ((n : ℝ) + m - 2) / (n + m) := div_nonneg (by linarith) hN.le
  have hq1 : ((n : ℝ) + m - 2) / (n + m) < 1 := (div_lt_one hN).mpr (by linarith)
  exact ⟨sq_nonneg _, by unfold popAdj; nlinarith⟩

lemma tails_mem (x : Fin n → ℝ) (y : Fin m → ℝ) (T : Targets) :
    dU25 x y T ∈ Icc (0 : ℝ) (1 / 4) ∧ dL25 x y T ∈ Icc (0 : ℝ) (1 / 4) ∧
      dU125 x y T ∈ Icc (0 : ℝ) (1 / 8) ∧ dL125 x y T ∈ Icc (0 : ℝ) (1 / 8) :=
  ⟨areaDist_mem_Icc _ _ _ _ _ (by norm_num), areaDist_mem_Icc _ _ _ _ _ (by norm_num),
    areaDist_mem_Icc _ _ _ _ _ (by norm_num), areaDist_mem_Icc _ _ _ _ _ (by norm_num)⟩

lemma rawScore_at_mem (x : Fin n → ℝ) (y : Fin m → ℝ) (meansOnly medians : Bool) (z : ℝ)
    (T : Targets) :
    rawScore meansOnly (centralDist x y medians z) (dU25 x y T) (dL25 x y T) (dU125 x y T)
      (dL125 x y T) ∈ Icc (0 : ℝ) 1 := by
  obtain ⟨h1, h2, h3, h4⟩ := tails_mem x y T
  exact rawScore_mem_Icc meansOnly (centralDist_mem_Icc x y medians z) h1 h2 h3 h4

theorem certaintyAt_eq_mul (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    (meansOnly medians : Bool) (z : ℝ) (T : Targets) :
    certaintyAt x y meansOnly medians z T =
      rawScore meansOnly (centralDist x y medians z) (dU25 x y T) (dL25 x y T) (dU125 x y T)
        (dL125 x y T) * popAdj n m := by
  have hr := rawScore_at_mem x y meansOnly medians z T
  have hp := popAdj_mem hn hm
  refine min_eq_right ?_
  nlinarith [hr.1, hr.2]

theorem certaintyAt_le_popAdj (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    (meansOnly medians : Bool) (z : ℝ) (T : Targets) :
    0 ≤ certaintyAt x y meansOnly medians z T ∧ certaintyAt x y meansOnly medians z T ≤ popAdj n m ∧
      popAdj n m < 1 := by
  have hr := rawScore_at_mem x y meansOnly medians z T
  have hp := popAdj_mem hn hm
  rw [certaintyAt_eq_mul hn hm]
  exact ⟨mul_nonneg hr.1 hp.1, by nlinarith [hr.1, hr.2], hp.2⟩

theorem isMedian_iff (μ : Measure ℝ) [IsProbabilityMeasure μ] (t : ℝ) :
    (1 / 2 ≤ μ.real (Iic t) ∧ 1 / 2 ≤ μ.real (Ici t)) ↔
      μ.real (Iio t) ≤ 1 / 2 ∧ 1 / 2 ≤ cdf μ t := by
  rw [cdf_eq_real, ← compl_Iio, probReal_compl_eq_one_sub measurableSet_Iio]
  constructor <;> rintro ⟨h1, h2⟩ <;> constructor <;> linarith

def medDist (a b : ℝ) : ℝ := max (max 0 (a - 1 / 2)) (1 / 2 - b)

theorem medianDist_eq (t : ℝ) (x : Fin n → ℝ) :
    medianDist t x = medDist (ltFrac t x) (LPM 0 t x) := by
  simp [medianDist, medDist, LPMratio]

lemma ltFrac_le_LPM_zero (t : ℝ) (x : Fin n → ℝ) : ltFrac t x ≤ LPM 0 t x := by
  unfold ltFrac LPM lowerDev
  refine div_le_div_of_nonneg_right (sum_le_sum fun i _ => ?_) (Nat.cast_nonneg n)
  by_cases h : x i < t
  · simp [h, h.le]
  · simp only [h, ↓reduceIte]; split_ifs <;> norm_num

theorem medDist_le {a b : ℝ} (s : ℝ) (hs : s ∈ Icc a b) : medDist a b ≤ |s - 1 / 2| := by
  unfold medDist
  refine max_le (max_le (abs_nonneg _) ?_) ?_
  · exact (sub_le_sub_right hs.1 _).trans (le_abs_self _)
  · exact (by linarith [hs.2] : 1 / 2 - b ≤ -(s - 1 / 2)).trans (neg_le_abs _)

theorem exists_medDist_eq {a b : ℝ} (hab : a ≤ b) :
    ∃ s ∈ Icc a b, medDist a b = |s - 1 / 2| := by
  unfold medDist
  rcases le_or_gt a (1 / 2) with h1 | h1
  · rcases le_or_gt (1 / 2) b with h2 | h2
    · refine ⟨1 / 2, ⟨h1, h2⟩, ?_⟩
      simp only [sub_self, abs_zero]
      exact le_antisymm (max_le (max_le le_rfl (by linarith)) (by linarith))
        (le_max_of_le_left (le_max_left _ _))
    · refine ⟨b, ⟨hab, le_rfl⟩, ?_⟩
      rw [abs_of_neg (by linarith), max_eq_right (max_le (by linarith) (by linarith))]
      ring
  · refine ⟨a, ⟨le_rfl, hab⟩, ?_⟩
    rw [abs_of_pos (by linarith), max_eq_right (by linarith : (0 : ℝ) ≤ a - 1 / 2),
      max_eq_left (by linarith : 1 / 2 - b ≤ a - 1 / 2)]

theorem medDist_eq_zero_iff (a b : ℝ) : medDist a b = 0 ↔ a ≤ 1 / 2 ∧ 1 / 2 ≤ b := by
  unfold medDist
  constructor
  · intro h
    have h1 := le_max_right (max 0 (a - 1 / 2)) (1 / 2 - b)
    have h2 := (le_max_right 0 (a - 1 / 2)).trans (le_max_left _ (1 / 2 - b))
    constructor <;> linarith
  · rintro ⟨h1, h2⟩
    rw [max_eq_left (show a - 1 / 2 ≤ 0 by linarith), max_eq_left (show 1 / 2 - b ≤ 0 by linarith)]

theorem medDist_self (a : ℝ) : medDist a a = |a - 1 / 2| := by
  obtain ⟨s, hs, h⟩ := exists_medDist_eq (le_refl a)
  rwa [le_antisymm hs.2 hs.1] at h

theorem medDist_eq_abs_of_no_atom (μ : Measure ℝ) [IsProbabilityMeasure μ] {t : ℝ}
    (ht : μ {t} = 0) : medDist (μ.real (Iio t)) (cdf μ t) = |cdf μ t - 1 / 2| := by
  have : μ.real (Iio t) = cdf μ t := by
    have hd : Disjoint (Iio t) {t} := by simp
    rw [cdf_eq_real, ← Iio_union_right, measureReal_union (μ := μ) hd (measurableSet_singleton t)]
    simp [measureReal_def, ht]
  rw [this, medDist_self]

theorem medianDist_eq_abs_of_no_atom {t : ℝ} {x : Fin n → ℝ} (ht : ∀ i, x i ≠ t) :
    medianDist t x = |LPM 0 t x - 1 / 2| := by
  have : ltFrac t x = LPM 0 t x := by
    unfold ltFrac LPM lowerDev
    congr 1
    refine sum_congr rfl fun i _ => ?_
    by_cases h : x i < t
    · simp [h, h.le]
    · have : ¬ x i ≤ t := fun h' => h (lt_of_le_of_ne h' (ht i))
      simp [h, this]
  rw [medianDist_eq, this, medDist_self]

theorem sample_isMedian_iff (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    ((n : ℝ) / 2 ≤ ∑ i, (if x i ≤ t then (1 : ℝ) else 0) ∧
        (n : ℝ) / 2 ≤ ∑ i, (if t ≤ x i then (1 : ℝ) else 0)) ↔ medianDist t x = 0 := by
  have hn' : (0 : ℝ) < n := Nat.cast_pos.mpr hn
  have hsplit : ∑ i, (if t ≤ x i then (1 : ℝ) else 0) =
      n - ∑ i, (if x i < t then (1 : ℝ) else 0) := by
    rw [eq_sub_iff_add_eq, ← sum_add_distrib]
    rw [show (n : ℝ) = ∑ _i : Fin n, (1 : ℝ) by simp]
    refine sum_congr rfl fun i _ => ?_
    by_cases h : x i < t
    · simp [h, not_le.mpr h]
    · simp [h, not_lt.mp h]
  have hL : LPM 0 t x = (∑ i, (if x i ≤ t then (1 : ℝ) else 0)) / n := by
    simp [LPM, lowerDev]
  rw [medianDist_eq, medDist_eq_zero_iff, hsplit, hL, ltFrac, div_le_iff₀ hn', le_div_iff₀ hn']
  constructor <;> rintro ⟨h1, h2⟩ <;> constructor <;> linarith

lemma sum_range_ind_le {a : ℕ → ℝ} {N j : ℕ} (hj : j ≤ N) {t : ℝ}
    (h : ∀ k, j ≤ k → k < N → t ≤ a k) :
    ∑ k ∈ range N, (if a k < t then (1 : ℝ) else 0) ≤ j := by
  rw [← sum_range_add_sum_Ico _ hj]
  have h2 : ∑ k ∈ Ico j N, (if a k < t then (1 : ℝ) else 0) = 0 :=
    sum_eq_zero fun k hk => by
      rw [Finset.mem_Ico] at hk
      simp [not_lt.mpr (h k hk.1 hk.2)]
  rw [h2, add_zero]
  calc ∑ k ∈ range j, (if a k < t then (1 : ℝ) else 0) ≤ ∑ _k ∈ range j, (1 : ℝ) :=
        sum_le_sum fun k _ => by split_ifs <;> norm_num
    _ = j := by simp

lemma le_sum_range_ind {a : ℕ → ℝ} {N j : ℕ} (hj : j ≤ N) {t : ℝ}
    (h : ∀ k, k < j → a k ≤ t) :
    (j : ℝ) ≤ ∑ k ∈ range N, (if a k ≤ t then (1 : ℝ) else 0) := by
  rw [← sum_range_add_sum_Ico _ hj]
  have h1 : ∑ k ∈ range j, (if a k ≤ t then (1 : ℝ) else 0) = j := by
    rw [sum_congr rfl fun k hk => (show (if a k ≤ t then (1 : ℝ) else 0) = 1 by
      simp [h k (mem_range.mp hk)])]
    simp
  rw [h1]
  refine le_add_of_nonneg_right (sum_nonneg fun k _ => ?_)
  split_ifs <;> norm_num

theorem medianDist_rMedian (hN : 0 < n) (z : Fin n → ℝ) : medianDist (rMedian z) z = 0 := by
  rw [← sample_isMedian_iff hN]
  have hs := sortedAt_sorted z
  set s := sortedAt z
  have hsum : ∀ (P : ℝ → Prop) [DecidablePred P],
      ∑ i, (if P (z i) then (1 : ℝ) else 0) = ∑ k ∈ range n, (if P (s k) then (1 : ℝ) else 0) :=
    fun P _ => (sum_sortedAt z fun v => if P v then (1 : ℝ) else 0).symm
  have hcompl : ∀ t, ∑ i, (if t ≤ z i then (1 : ℝ) else 0) =
      n - ∑ i, (if z i < t then (1 : ℝ) else 0) := fun t => by
    rw [eq_sub_iff_add_eq, ← sum_add_distrib, show (n : ℝ) = ∑ _i : Fin n, (1 : ℝ) by simp]
    refine sum_congr rfl fun i _ => ?_
    by_cases h : z i < t
    · simp [h, not_le.mpr h]
    · simp [h, not_lt.mp h]
  obtain ⟨j₁, j₂, hj₁, hj₂, hlo, hhi, hc1, hc2⟩ : ∃ j₁ j₂ : ℕ, j₁ ≤ n ∧ j₂ ≤ n ∧
      (∀ k, k < j₁ → s k ≤ rMedian z) ∧ (∀ k, j₂ ≤ k → k < n → rMedian z ≤ s k) ∧
      (n : ℝ) / 2 ≤ j₁ ∧ (j₂ : ℝ) ≤ n / 2 := by
    unfold rMedian
    split_ifs with hodd
    · have hc : n = 2 * (n / 2) + 1 := by omega
      have hr : (n : ℝ) = 2 * ((n / 2 : ℕ) : ℝ) + 1 := by exact_mod_cast hc
      refine ⟨n / 2 + 1, n / 2, by omega, by omega, fun k hk => hs _ _ (by omega) (by omega),
        fun k hk hkn => hs _ _ hk hkn, ?_, ?_⟩
      · push_cast; linarith
      · linarith
    · have hc : n = 2 * (n / 2) := by omega
      have hr : (n : ℝ) = 2 * ((n / 2 : ℕ) : ℝ) := by exact_mod_cast hc
      have hle := hs (n / 2 - 1) (n / 2) (by omega) (by omega)
      refine ⟨n / 2, n / 2, by omega, by omega, fun k hk => ?_, fun k hk hkn => ?_, ?_, ?_⟩
      · have := hs k (n / 2 - 1) (by omega) (by omega); linarith
      · have := hs (n / 2) k hk hkn; linarith
      · linarith
      · linarith
  constructor
  · rw [hsum (· ≤ rMedian z)]
    exact hc1.trans (le_sum_range_ind hj₁ hlo)
  · rw [hcompl, hsum (· < rMedian z)]
    have := sum_range_ind_le hj₂ hhi
    linarith

lemma ltFrac_append_self (t : ℝ) (x : Fin n → ℝ) : ltFrac t (Fin.append x x) = ltFrac t x := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [ltFrac]
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  simp only [ltFrac, Fin.sum_univ_add, Fin.append_left, Fin.append_right, Nat.cast_add]
  rw [← two_mul, ← two_mul, mul_div_mul_left _ _ two_ne_zero]

lemma LPM_zero_append_self (t : ℝ) (x : Fin n → ℝ) : LPM 0 t (Fin.append x x) = LPM 0 t x := by
  simp only [LPM, Fin.sum_univ_add, Fin.append_left, Fin.append_right, Nat.cast_add]
  rw [← two_mul, ← two_mul, mul_div_mul_left _ _ two_ne_zero]

theorem centralDist_median_self (hn : 0 < n) (x : Fin n → ℝ) :
    centralDist x x true (grandStat x x true) = 0 := by
  simp only [centralDist, grandStat, ↓reduceIte, max_self]
  have h := medianDist_rMedian (by omega : 0 < n + n) (Fin.append x x)
  rwa [medianDist_eq, ltFrac_append_self, LPM_zero_append_self, ← medianDist_eq] at h

theorem LPM_one_const {x : Fin n → ℝ} {c : ℝ} (h : ∀ i, x i = c) :
    LPM 1 c x = 0 ∧ UPM 1 c x = 0 := by
  simp [LPM, UPM, lowerDev, upperDev, h]

theorem balance_const {x : Fin n → ℝ} {c : ℝ} (h : ∀ i, x i = c) (p : ℝ) :
    (1 - p) * LPM 1 c x = p * UPM 1 c x := by
  rw [(LPM_one_const h).1, (LPM_one_const h).2, mul_zero, mul_zero]

theorem areaDist_const_target (t : ℝ) (cy : Option ℝ) (r1 r2 p : ℝ) :
    areaDist (some t) cy r1 r2 t p = min p (if cy = some t then 0 else |r2 - p|) := by
  have h : 0 ≤ (if cy = some t then 0 else |r2 - p|) := by split_ifs <;> simp
  simp only [areaDist, ↓reduceIte, max_eq_right h]

theorem areaDist_const_self (t r1 r2 : ℝ) {p : ℝ} (hp : 0 ≤ p) :
    areaDist (some t) (some t) r1 r2 t p = 0 := by
  simp [areaDist, hp]

lemma constVal_ne_of_nonconst {y : Fin m → ℝ} (hy : ∃ i j, y i ≠ y j) (c : ℝ) :
    constVal y ≠ some c := fun h => by
  obtain ⟨i, j, hij⟩ := hy
  obtain ⟨-, h'⟩ := constVal_eq_some.mp h
  exact hij ((h' i).trans (h' j).symm)

lemma not_nonconst {x : Fin n → ℝ} {c : ℝ} (h : ∀ i, x i = c) : ¬ ∃ i j, x i ≠ x j := by
  rintro ⟨i, j, hij⟩; exact hij ((h i).trans (h j).symm)

lemma sortedAt_const {x : Fin n → ℝ} {c : ℝ} (h : ∀ i, x i = c) {k : ℕ} (hk : k < n) :
    sortedAt x k = c := by
  obtain ⟨i, hi⟩ := exists_eq_sortedAt x hk
  rw [← hi, h]

lemma LPMVaR_one_const (hn : 0 < n) {x : Fin n → ℝ} {c : ℝ} (h : ∀ i, x i = c) (p : ℝ) :
    LPMVaR 1 p x = c := by
  rw [LPMVaR_of_const hn (not_nonconst h) one_ne_zero, varMin, sortedAt_const h hn]

lemma UPMVaR_one_const (hn : 0 < n) {x : Fin n → ℝ} {c : ℝ} (h : ∀ i, x i = c) (p : ℝ) :
    UPMVaR 1 p x = c := by
  rw [UPMVaR_eq_LPMVaR hn, LPMVaR_one_const hn h]

lemma rMedian_const {z : Fin n → ℝ} {c : ℝ} (hn : 0 < n) (h : ∀ i, z i = c) : rMedian z = c := by
  unfold rMedian
  split_ifs
  · exact sortedAt_const h (by omega)
  · rw [sortedAt_const h (by omega), sortedAt_const h (by omega)]; ring

lemma append_const {x : Fin n → ℝ} {y : Fin m → ℝ} {c : ℝ} (hx : ∀ i, x i = c)
    (hy : ∀ j, y j = c) (i : Fin (n + m)) : Fin.append x y i = c := by
  refine Fin.addCases (fun i => ?_) (fun j => ?_) i
  · simp [hx]
  · simp [hy]

theorem certainty_const_const (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ} {y : Fin m → ℝ}
    {c : ℝ} (hx : ∀ i, x i = c) (hy : ∀ j, y j = c) (meansOnly medians : Bool) :
    certainty x y meansOnly medians = popAdj n m := by
  have hcx : constVal x = some c := constVal_eq_some.mpr ⟨hn, hx⟩
  have hcy : constVal y = some c := constVal_eq_some.mpr ⟨hm, hy⟩
  have hz : grandStat x y medians = c := by
    unfold grandStat
    split_ifs
    · exact rMedian_const (by omega) (append_const hx hy)
    · simp [hcx]
    · simp_all
  have hT : defaultTargets x y = ⟨c, c, c, c⟩ := by
    simp only [defaultTargets, UPMVaR_one_const hn hx, UPMVaR_one_const hm hy,
      LPMVaR_one_const hn hx, LPMVaR_one_const hm hy, add_self_div_two]
  have hmed : ∀ {k : ℕ} {w : Fin k → ℝ}, 0 < k → (∀ i, w i = c) → medianDist c w = 0 :=
    fun {k w} hk hw => by
      rw [medianDist_eq, medDist_eq_zero_iff]
      have h1 : ltFrac c w = 0 := by simp [ltFrac, hw]
      have h2 : LPM 0 c w = 1 := by
        have hk' : (k : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hk.ne'
        simp [LPM, lowerDev, hw, hk']
      rw [h1, h2]; norm_num
  have hc0 : centralDist x y medians c = 0 := by
    unfold centralDist
    split_ifs
    · rw [hmed hn hx, hmed hm hy, max_self]
    · rw [hcx, hcy, areaDist_const_self _ _ _ (by norm_num)]
  rw [certainty, certaintyAt_eq_mul hn hm, hz, hT, hc0]
  simp only [dU25, dL25, dU125, dL125, hcx, hcy]
  rw [areaDist_const_self _ _ _ (by norm_num), areaDist_const_self _ _ _ (by norm_num),
    areaDist_const_self _ _ _ (by norm_num), areaDist_const_self _ _ _ (by norm_num),
    (rawScore_zero meansOnly).1, one_mul]

theorem certainty_const_var_meansOnly (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ}
    {y : Fin m → ℝ} {c : ℝ} (hx : ∀ i, x i = c) (hy : ∃ i j, y i ≠ y j) (hmean : mean y = c) :
    certainty x y true false = popAdj n m := by
  have hcx : constVal x = some c := constVal_eq_some.mpr ⟨hn, hx⟩
  have hz : grandStat x y false = c := by
    have hne : ¬ ((constVal x).isSome ∧ constVal x = constVal y) := fun h =>
      constVal_ne_of_nonconst hy c (h.2 ▸ hcx)
    simp only [grandStat, Bool.false_eq_true, ↓reduceIte, hne]
    have hN : (n : ℝ) + m ≠ 0 := by positivity
    rw [pooledMean, mean_const hn hx, hmean]
    field_simp
  have hc0 : centralDist x y false c = 0 := by
    simp only [centralDist, Bool.false_eq_true, ↓reduceIte, areaDist, hcx,
      constVal_ne_of_nonconst hy c]
    rw [← hmean, LPMratio_one_mean hy]
    norm_num
  rw [certainty, certaintyAt_eq_mul hn hm, hz, hc0]
  simp only [rawScore, ↓reduceIte]
  norm_num

theorem centralDist_eq_piecewise (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    (z : ℝ) :
    centralDist x y false z =
        min (1 / 2) (max (if ∀ i, x i = z then 0 else |LPMratio 1 z x - 1 / 2|)
          (if ∀ j, y j = z then 0 else |LPMratio 1 z y - 1 / 2|)) ∧
      centralDist x y false z =
        max (if ∀ i, x i = z then 0 else |LPMratio 1 z x - 1 / 2|)
          (if ∀ j, y j = z then 0 else |LPMratio 1 z y - 1 / 2|) := by
  have hx : (constVal x = some z) ↔ ∀ i, x i = z := by rw [constVal_eq_some]; simp [hn]
  have hy : (constVal y = some z) ↔ ∀ j, y j = z := by rw [constVal_eq_some]; simp [hm]
  have h1 : centralDist x y false z =
      min (1 / 2) (max (if ∀ i, x i = z then 0 else |LPMratio 1 z x - 1 / 2|)
        (if ∀ j, y j = z then 0 else |LPMratio 1 z y - 1 / 2|)) := by
    simp only [centralDist, Bool.false_eq_true, ↓reduceIte, areaDist]
    congr 2
    · exact if_congr hx rfl rfl
    · exact if_congr hy rfl rfl
  refine ⟨h1, h1.trans (min_eq_right (max_le ?_ ?_))⟩
  · split_ifs
    · norm_num
    · exact abs_sub_half_le (LPMratio_mem_Icc 1 z x)
  · split_ifs
    · norm_num
    · exact abs_sub_half_le (LPMratio_mem_Icc 1 z y)

lemma areaDist_const_off {x : Fin n → ℝ} {c w p : ℝ} (hn : 0 < n) (hx : ∀ i, x i = c)
    (hw : w ≠ c) (hp0 : 0 ≤ p) (hp : p ≤ 1 / 2) (cy : Option ℝ) (r1 r2 : ℝ)
    (hr : r1 = 0 ∨ r1 = 1) : areaDist (constVal x) cy r1 r2 w p = p := by
  have hcx : constVal x = some c := constVal_eq_some.mpr ⟨hn, hx⟩
  have hne : constVal x ≠ some w := by
    rw [hcx]; intro h; exact hw (Option.some_injective _ h).symm
  have hd : p ≤ |r1 - p| := by
    rcases hr with rfl | rfl
    · rw [abs_of_nonpos (by linarith)]; linarith
    · rw [abs_of_nonneg (by linarith)]; linarith
  simp only [areaDist, hne, ↓reduceIte]
  exact min_eq_left (le_max_of_le_left hd)

lemma UPMratio_const_off {x : Fin n → ℝ} {c w : ℝ} (hn : 0 < n) (hx : ∀ i, x i = c)
    (hw : w ≠ c) : UPMratio 1 w x = 0 ∨ UPMratio 1 w x = 1 := by
  rcases lt_or_gt_of_ne hw with h | h
  · exact Or.inr (UPMratio_eq_one_of_forall_lt hn 1 fun i => by rw [hx]; exact h)
  · exact Or.inl (UPMratio_eq_zero_of_forall_le 1 fun i => by rw [hx]; exact h.le)

lemma LPMratio_const_off {x : Fin n → ℝ} {c w : ℝ} (hn : 0 < n) (hx : ∀ i, x i = c)
    (hw : w ≠ c) : LPMratio 1 w x = 0 ∨ LPMratio 1 w x = 1 := by
  rcases lt_or_gt_of_ne hw with h | h
  · exact Or.inl (LPMratio_eq_zero_of_forall_lt 1 fun i => by rw [hx]; exact h)
  · exact Or.inr (LPMratio_eq_one_of_forall_le hn 1 fun i => by rw [hx]; exact h.le)

lemma UPMVaR_ne_mean {y : Fin m → ℝ} (hy : ∃ i j, y i ≠ y j) {p : ℝ} (hp : p ∈ Icc (0 : ℝ) 1)
    (hgap : rTol < |1 / 2 - p|) : UPMVaR 1 p y ≠ mean y := fun h => by
  have h1 := upperTarget_ratio hy hp
  rw [h, UPMratio_one_mean hy] at h1
  linarith

lemma LPMVaR_ne_mean {y : Fin m → ℝ} (hy : ∃ i j, y i ≠ y j) {p : ℝ} (hp : p ∈ Icc (0 : ℝ) 1)
    (hgap : rTol < |1 / 2 - p|) : LPMVaR 1 p y ≠ mean y := fun h => by
  have h1 := lowerTarget_ratio hy hp
  rw [h, LPMratio_one_mean hy] at h1
  linarith

theorem certainty_const_var_full (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ} {y : Fin m → ℝ}
    {c : ℝ} (hx : ∀ i, x i = c) (hy : ∃ i j, y i ≠ y j) (hmean : mean y = c) :
    certainty x y false false = 2 / 5 * popAdj n m := by
  have hcx : constVal x = some c := constVal_eq_some.mpr ⟨hn, hx⟩
  have hz : grandStat x y false = c := by
    have hne : ¬ ((constVal x).isSome ∧ constVal x = constVal y) := fun h =>
      constVal_ne_of_nonconst hy c (h.2 ▸ hcx)
    simp only [grandStat, Bool.false_eq_true, ↓reduceIte, hne]
    have hN : (n : ℝ) + m ≠ 0 := by positivity
    rw [pooledMean, mean_const hn hx, hmean]
    field_simp
  have hc0 : centralDist x y false c = 0 := by
    simp only [centralDist, Bool.false_eq_true, ↓reduceIte, areaDist, hcx,
      constVal_ne_of_nonconst hy c]
    rw [← hmean, LPMratio_one_mean hy]
    norm_num
  have g14 : rTol < |1 / 2 - (1 / 4 : ℝ)| := by norm_num [rTol, abs_of_pos]
  have g18 : rTol < |1 / 2 - (1 / 8 : ℝ)| := by norm_num [rTol, abs_of_pos]
  have m14 : (1 / 4 : ℝ) ∈ Icc (0 : ℝ) 1 := ⟨by norm_num, by norm_num⟩
  have m18 : (1 / 8 : ℝ) ∈ Icc (0 : ℝ) 1 := ⟨by norm_num, by norm_num⟩
  have avg : ∀ {v : ℝ}, v ≠ c → (c + v) / 2 ≠ c := fun hv h => hv (by linarith)
  have tU25 : (defaultTargets x y).u25 ≠ c := by
    simp only [defaultTargets, UPMVaR_one_const hn hx]
    exact avg (hmean ▸ UPMVaR_ne_mean hy m14 g14)
  have tL25 : (defaultTargets x y).l25 ≠ c := by
    simp only [defaultTargets, LPMVaR_one_const hn hx]
    exact avg (hmean ▸ LPMVaR_ne_mean hy m14 g14)
  have tU125 : (defaultTargets x y).u125 ≠ c := by
    simp only [defaultTargets, UPMVaR_one_const hn hx]
    exact avg (hmean ▸ UPMVaR_ne_mean hy m18 g18)
  have tL125 : (defaultTargets x y).l125 ≠ c := by
    simp only [defaultTargets, LPMVaR_one_const hn hx]
    exact avg (hmean ▸ LPMVaR_ne_mean hy m18 g18)
  have hU25 : dU25 x y (defaultTargets x y) = 1 / 4 :=
    areaDist_const_off hn hx tU25 (by norm_num) (by norm_num) _ _ _ (UPMratio_const_off hn hx tU25)
  have hL25 : dL25 x y (defaultTargets x y) = 1 / 4 :=
    areaDist_const_off hn hx tL25 (by norm_num) (by norm_num) _ _ _ (LPMratio_const_off hn hx tL25)
  have hU125 : dU125 x y (defaultTargets x y) = 1 / 8 :=
    areaDist_const_off hn hx tU125 (by norm_num) (by norm_num) _ _ _
      (UPMratio_const_off hn hx tU125)
  have hL125 : dL125 x y (defaultTargets x y) = 1 / 8 :=
    areaDist_const_off hn hx tL125 (by norm_num) (by norm_num) _ _ _
      (LPMratio_const_off hn hx tL125)
  rw [certainty, certaintyAt_eq_mul hn hm, hz, hc0, hU25, hL25, hU125, hL125]
  simp only [rawScore, Bool.false_eq_true, ↓reduceIte]
  norm_num

theorem certainty_chapter_example :
    certainty (fun _ : Fin 10 => (0 : ℝ)) (fun j : Fin 10 => (-1 : ℝ) ^ (j : ℕ)) true false
        = 81 / 100 ∧
      certainty (fun _ : Fin 10 => (0 : ℝ)) (fun j : Fin 10 => (-1 : ℝ) ^ (j : ℕ)) false false
        = 81 / 250 := by
  have hy : ∃ i j, (fun j : Fin 10 => (-1 : ℝ) ^ (j : ℕ)) i ≠
      (fun j : Fin 10 => (-1 : ℝ) ^ (j : ℕ)) j := ⟨0, 1, by norm_num⟩
  have hmean : mean (fun j : Fin 10 => (-1 : ℝ) ^ (j : ℕ)) = 0 := by
    norm_num [mean, Fin.sum_univ_succ]
  have hp : popAdj 10 10 = 81 / 100 := by norm_num [popAdj]
  refine ⟨?_, ?_⟩
  · rw [certainty_const_var_meansOnly (x := fun _ : Fin 10 => (0 : ℝ)) (c := 0) (by norm_num)
      (by norm_num) (fun _ => rfl) hy hmean, hp]
  · rw [certainty_const_var_full (x := fun _ : Fin 10 => (0 : ℝ)) (c := 0) (by norm_num)
      (by norm_num) (fun _ => rfl) hy hmean, hp]
    norm_num

def c01 : ℝ := 3602879701896397 / 36028797018963968

def c01ctrl : Fin 2 → ℝ := fun _ => c01

def c01trt : Fin 4 → ℝ := ![c01 - 1 / 16, c01 + 1 / 16, c01, c01]

theorem certainty_c01 :
    mean c01ctrl = c01 ∧ mean c01trt = c01 ∧ grandStat c01ctrl c01trt false = c01 ∧
      certainty c01ctrl c01trt true false = 4 / 9 := by
  have hy : ∃ i j, c01trt i ≠ c01trt j := ⟨0, 1, by norm_num [c01trt, c01]⟩
  have hmy : mean c01trt = c01 := by norm_num [mean, Fin.sum_univ_succ, c01trt]; ring
  have hmx : mean c01ctrl = c01 := mean_const (by norm_num) fun _ => rfl
  refine ⟨hmx, hmy, ?_, ?_⟩
  · rw [grandStat_eq_mean_append, ← pooledMean_eq_mean_append, pooledMean, hmx, hmy]; ring
  · rw [certainty_const_var_meansOnly (x := c01ctrl) (c := c01) (by norm_num) (by norm_num)
      (fun _ => rfl) hy hmy]
    norm_num [popAdj]

def sepCtrl : Fin 3 → ℝ := fun _ => 0

def sepTrt : Fin 3 → ℝ := ![-1, -1, 2]

lemma sep_pooled_median : grandStat sepCtrl sepTrt true = 0 := by
  have h := medianDist_rMedian (by norm_num : 0 < 3 + 3) (Fin.append sepCtrl sepTrt)
  rw [medianDist_eq, medDist_eq_zero_iff] at h
  set t := rMedian (Fin.append sepCtrl sepTrt)
  change t = 0
  obtain ⟨h1, h2⟩ := h
  simp only [ltFrac, LPM, lowerDev, pow_zero] at h1 h2
  rw [Fin.sum_univ_add] at h1 h2
  simp only [Fin.append_left, Fin.append_right, Fin.sum_univ_three, sepCtrl, sepTrt,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
    Matrix.tail_cons] at h1 h2
  by_contra ht
  rcases lt_or_gt_of_ne ht with ht | ht
  · have a : ¬ (0 : ℝ) ≤ t := not_le.mpr ht
    have b : ¬ (2 : ℝ) ≤ t := by intro h; linarith
    simp only [a, b, ↓reduceIte] at h2
    by_cases hb : (-1 : ℝ) ≤ t <;> simp only [hb, ↓reduceIte] at h2 <;> norm_num at h2
  · have a : (0 : ℝ) < t := ht
    have b : (-1 : ℝ) < t := by linarith
    simp only [a, b, ↓reduceIte] at h1
    by_cases hb : (2 : ℝ) < t <;> simp only [hb, ↓reduceIte] at h1 <;> norm_num at h1

theorem certainty_modes_distinct :
    certainty sepCtrl sepTrt true false = 4 / 9 ∧
      certainty sepCtrl sepTrt true true = 16 / 81 ∧
      certainty sepCtrl sepTrt false false = 8 / 45 := by
  have hy : ∃ i j, sepTrt i ≠ sepTrt j := ⟨0, 2, by norm_num [sepTrt]⟩
  have hmean : mean sepTrt = 0 := by norm_num [mean, Fin.sum_univ_succ, sepTrt]
  have hp : popAdj 3 3 = 4 / 9 := by norm_num [popAdj]
  refine ⟨?_, ?_, ?_⟩
  · rw [certainty_const_var_meansOnly (x := sepCtrl) (c := 0) (by norm_num) (by norm_num)
      (fun _ => rfl) hy hmean, hp]
  · have hx0 : medianDist 0 sepCtrl = 0 := by
      rw [medianDist_eq, medDist_eq_zero_iff]
      simp [ltFrac, LPM, lowerDev, sepCtrl]
      norm_num
    have hy0 : medianDist 0 sepTrt = 1 / 6 := by
      rw [medianDist_eq]
      simp only [medDist, ltFrac, LPM, lowerDev, Fin.sum_univ_three, sepTrt]
      norm_num
    rw [certainty, certaintyAt_eq_mul (by norm_num) (by norm_num), sep_pooled_median, hp]
    simp only [centralDist, ↓reduceIte, hx0, hy0, rawScore]
    norm_num
  · rw [certainty_const_var_full (x := sepCtrl) (c := 0) (by norm_num) (by norm_num)
      (fun _ => rfl) hy hmean, hp]
    norm_num

theorem controlCDF_medians (x : Fin n → ℝ) (z : ℝ) : controlCDF x true z = LPM 0 z x := by
  simp [controlCDF, LPMratio]

end

end NNS.GapCertainty
