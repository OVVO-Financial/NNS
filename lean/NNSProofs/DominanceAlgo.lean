-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.Integration
public import NNSProofs.Ratios

namespace NNS

open Finset Set

@[expose] public noncomputable section

variable {n m : ℕ}

def sortedVals (x : Fin n → ℝ) : List ℝ := (univ.val.map x).sort (· ≤ ·)

def prefix1 (x : Fin n → ℝ) (k : ℕ) : ℝ := ((sortedVals x).take k).sum

def prefix2 (x : Fin n → ℝ) (k : ℕ) : ℝ := (((sortedVals x).take k).map (· ^ 2)).sum

def total1 (x : Fin n → ℝ) : ℝ := prefix1 x n

def sampleMin (x : Fin n → ℝ) : ℝ := (sortedVals x).headD 0

def sampleMean (x : Fin n → ℝ) : ℝ := total1 x / n

def walkCount (x : Fin n → ℝ) (t : ℝ) : ℕ :=
  ((sortedVals x).takeWhile fun v => decide (v ≤ t)).length

abbrev identicalSamples (x : Fin n → ℝ) (y : Fin m → ℝ) : Prop := sortedVals x = sortedVals y

def lpmDeg1 (x : Fin n → ℝ) (k : ℕ) (t : ℝ) : ℝ := ((k : ℝ) * t - prefix1 x k) / n

def upmDeg1 (x : Fin n → ℝ) (k : ℕ) (t : ℝ) : ℝ :=
  ((total1 x - prefix1 x k) - ((n : ℝ) - k) * t) / n

def lpmDeg2 (x : Fin n → ℝ) (k : ℕ) (t : ℝ) : ℝ :=
  ((k : ℝ) * t ^ 2 - 2 * t * prefix1 x k + prefix2 x k) / n

lemma coe_sortedVals (x : Fin n → ℝ) : ((sortedVals x : List ℝ) : Multiset ℝ) = univ.val.map x :=
  Multiset.sort_eq _ _

lemma sortedVals_pairwise (x : Fin n → ℝ) : (sortedVals x).Pairwise (· ≤ ·) :=
  Multiset.pairwise_sort _ _

lemma length_sortedVals (x : Fin n → ℝ) : (sortedVals x).length = n := by
  rw [← Multiset.coe_card, coe_sortedVals, Multiset.card_map, card_val, card_univ,
    Fintype.card_fin]

lemma mem_sortedVals {x : Fin n → ℝ} {v : ℝ} : v ∈ sortedVals x ↔ ∃ i, x i = v := by
  rw [← Multiset.mem_coe, coe_sortedVals, Multiset.mem_map]
  simp

lemma takeWhile_le_eq_filter {l : List ℝ} (hl : l.Pairwise (· ≤ ·)) (t : ℝ) :
    l.takeWhile (fun v => decide (v ≤ t)) = l.filter (fun v => decide (v ≤ t)) := by
  induction l with
  | nil => rfl
  | cons a l ih =>
    rw [List.pairwise_cons] at hl
    by_cases ha : a ≤ t
    · rw [List.takeWhile_cons_of_pos (by simpa using ha),
        List.filter_cons_of_pos (by simpa using ha), ih hl.2]
    · rw [List.takeWhile_cons_of_neg (by simpa using ha),
        List.filter_cons_of_neg (by simpa using ha), eq_comm, List.filter_eq_nil_iff]
      intro b hb
      have := hl.1 b hb
      simp only [decide_eq_true_eq]
      intro hbt
      exact ha (this.trans hbt)

lemma take_walkCount (x : Fin n → ℝ) (t : ℝ) :
    (sortedVals x).take (walkCount x t) = (sortedVals x).filter fun v => decide (v ≤ t) := by
  rw [walkCount, ← List.prefix_iff_eq_take.mp (List.takeWhile_prefix _),
    takeWhile_le_eq_filter (sortedVals_pairwise x)]

lemma sum_map_filter_sortedVals (x : Fin n → ℝ) (t : ℝ) (f : ℝ → ℝ) :
    (((sortedVals x).filter fun v => decide (v ≤ t)).map f).sum
      = ∑ i ∈ univ.filter (fun i => x i ≤ t), f (x i) := by
  rw [sum_eq_multiset_sum, ← Multiset.sum_coe, ← Multiset.map_coe, ← Multiset.filter_coe,
    coe_sortedVals, Multiset.filter_map, Multiset.map_map, filter_val]
  rfl

theorem walkCount_eq_lowerCount (x : Fin n → ℝ) (t : ℝ) : walkCount x t = lowerCount t x := by
  rw [walkCount, takeWhile_le_eq_filter (sortedVals_pairwise x), lowerCount, ← Multiset.coe_card,
    ← Multiset.filter_coe, coe_sortedVals, Multiset.filter_map, Multiset.card_map, card_def,
    filter_val]
  rfl

theorem prefix1_walkCount (x : Fin n → ℝ) (t : ℝ) :
    prefix1 x (walkCount x t) = lowerSum1 t x := by
  have h := sum_map_filter_sortedVals x t id
  rw [List.map_id] at h
  rw [prefix1, take_walkCount, lowerSum1, h]
  rfl

theorem prefix2_walkCount (x : Fin n → ℝ) (t : ℝ) :
    prefix2 x (walkCount x t) = lowerSum2 t x := by
  rw [prefix2, take_walkCount, lowerSum2, ← sum_map_filter_sortedVals x t (· ^ 2)]

theorem total1_eq_sum (x : Fin n → ℝ) : total1 x = ∑ i, x i := by
  rw [total1, prefix1, List.take_of_length_le (length_sortedVals x).le, sum_eq_multiset_sum,
    ← Multiset.sum_coe, coe_sortedVals]

theorem sampleMean_eq_mean (x : Fin n → ℝ) : sampleMean x = mean x := by
  rw [sampleMean, total1_eq_sum, mean]

theorem sampleMin_le (x : Fin n → ℝ) (i : Fin n) : sampleMin x ≤ x i := by
  have hmem : x i ∈ sortedVals x := mem_sortedVals.mpr ⟨i, rfl⟩
  have hp := sortedVals_pairwise x
  unfold sampleMin
  cases h : sortedVals x with
  | nil => rw [h] at hmem; simp at hmem
  | cons a l =>
    rw [h] at hmem hp
    rw [List.pairwise_cons] at hp
    rcases List.mem_cons.mp hmem with h' | h'
    · simp [h']
    · exact hp.1 _ h'

theorem exists_eq_sampleMin (hn : 0 < n) (x : Fin n → ℝ) : ∃ i, x i = sampleMin x := by
  unfold sampleMin
  cases h : sortedVals x with
  | nil =>
    have := length_sortedVals x
    rw [h] at this
    simp at this
    omega
  | cons a l =>
    obtain ⟨i, hi⟩ := mem_sortedVals.mp (show a ∈ sortedVals x by rw [h]; simp)
    exact ⟨i, hi⟩

theorem sampleMin_le_iff (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    sampleMin y ≤ sampleMin x ↔ ∀ i, ∃ j, y j ≤ x i := by
  constructor
  · intro h i
    obtain ⟨j, hj⟩ := exists_eq_sampleMin hm y
    exact ⟨j, hj ▸ h.trans (sampleMin_le x i)⟩
  · intro h
    obtain ⟨i, hi⟩ := exists_eq_sampleMin hn x
    obtain ⟨j, hj⟩ := h i
    exact (sampleMin_le y j).trans (hj.trans hi.le)

lemma LPM_eq_sortedVals (d : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    LPM d t x = ((sortedVals x).map (lowerDev d t)).sum / (sortedVals x).length := by
  rw [LPM, length_sortedVals, sum_eq_multiset_sum, ← Multiset.sum_coe, ← Multiset.map_coe,
    coe_sortedVals, Multiset.map_map]
  rfl

theorem LPM_eq_of_identical {x : Fin n → ℝ} {y : Fin m → ℝ} (h : identicalSamples x y)
    (d : ℕ) (t : ℝ) : LPM d t x = LPM d t y := by
  unfold identicalSamples at h
  rw [LPM_eq_sortedVals, LPM_eq_sortedVals, h]

theorem lpmDeg1_walkCount {x : Fin n → ℝ} {u w s : ℝ} (hx : ∀ i, x i ∉ Ioo u w)
    (hs : s ∈ Icc u w) : lpmDeg1 x (walkCount x u) s = LPM 1 s x := by
  rw [LPM_one_eq_on_Icc hx hs, lpmDeg1, walkCount_eq_lowerCount, ← walkCount_eq_lowerCount,
    prefix1_walkCount, walkCount_eq_lowerCount]

theorem lpmDeg2_walkCount {x : Fin n → ℝ} {u w s : ℝ} (hx : ∀ i, x i ∉ Ioo u w)
    (hs : s ∈ Icc u w) : lpmDeg2 x (walkCount x u) s = LPM 2 s x := by
  rw [LPM_two_eq_on_Icc hx hs, lpmDeg2, prefix1_walkCount, prefix2_walkCount,
    walkCount_eq_lowerCount]

lemma not_mem_Ioo_self (u v : ℝ) : v ∉ Ioo u u := fun h => lt_irrefl u (h.1.trans h.2)

theorem lpmDeg1_walkCount_self (x : Fin n → ℝ) (t : ℝ) :
    lpmDeg1 x (walkCount x t) t = LPM 1 t x :=
  lpmDeg1_walkCount (w := t) (fun _ => not_mem_Ioo_self t _) ⟨le_rfl, le_rfl⟩

theorem lpmDeg2_walkCount_self (x : Fin n → ℝ) (t : ℝ) :
    lpmDeg2 x (walkCount x t) t = LPM 2 t x :=
  lpmDeg2_walkCount (w := t) (fun _ => not_mem_Ioo_self t _) ⟨le_rfl, le_rfl⟩

theorem walkCount_div_eq_LPM_zero (x : Fin n → ℝ) (t : ℝ) :
    (walkCount x t : ℝ) / n = LPM 0 t x := by
  rw [walkCount_eq_lowerCount, LPM_zero_eq_ecdf, lowerCount]

def knotList (x : Fin n → ℝ) (y : Fin m → ℝ) : List ℝ :=
  (Finset.univ.image x ∪ Finset.univ.image y).sort (· ≤ ·)

def walkSteps : List ℝ → List (ℝ × Option ℝ)
  | [] => []
  | [t] => [(t, none)]
  | t :: t' :: rest => (t, some t') :: walkSteps (t' :: rest)

lemma mem_knotList {x : Fin n → ℝ} {y : Fin m → ℝ} {t : ℝ} :
    t ∈ knotList x y ↔ IsKnot x y t := by
  rw [knotList, mem_sort, isKnot_iff_mem]

lemma knotList_pairwise (x : Fin n → ℝ) (y : Fin m → ℝ) : (knotList x y).Pairwise (· < ·) :=
  ((pairwise_sort _ (· ≤ ·)).and (sort_nodup _ (· ≤ ·))).imp fun h => lt_of_le_of_ne h.1 h.2

lemma fst_mem_of_mem_walkSteps {l : List ℝ} {t : ℝ} {o : Option ℝ}
    (h : (t, o) ∈ walkSteps l) : t ∈ l := by
  induction l using walkSteps.induct with
  | case1 => simp [walkSteps] at h
  | case2 a => simp_all [walkSteps]
  | case3 a b rest ih =>
    simp only [walkSteps, List.mem_cons, Prod.mk.injEq] at h
    rcases h with ⟨rfl, -⟩ | h
    · simp
    · exact List.mem_cons_of_mem _ (ih h)

lemma adjacent_of_mem_walkSteps {l : List ℝ} (hl : l.Pairwise (· < ·)) {t w : ℝ}
    (h : (t, some w) ∈ walkSteps l) : t < w ∧ ∀ z ∈ l, z ∉ Ioo t w := by
  induction l using walkSteps.induct with
  | case1 => simp [walkSteps] at h
  | case2 a => simp [walkSteps] at h
  | case3 a b rest ih =>
    have hl' := List.pairwise_cons.mp hl
    simp only [walkSteps, List.mem_cons, Prod.mk.injEq, Option.some.injEq] at h
    rcases h with ⟨rfl, rfl⟩ | h
    · refine ⟨hl'.1 _ (by simp), fun z hz hzi => ?_⟩
      rcases List.mem_cons.mp hz with rfl | hz
      · exact lt_irrefl _ hzi.1
      rcases List.mem_cons.mp hz with rfl | hz
      · exact lt_irrefl _ hzi.2
      · exact lt_asymm hzi.2 ((List.pairwise_cons.mp hl'.2).1 z hz)
    · obtain ⟨h1, h2⟩ := ih hl'.2 h
      refine ⟨h1, fun z hz hzi => ?_⟩
      rcases List.mem_cons.mp hz with rfl | hz
      · exact lt_asymm hzi.1 (hl'.1 t (fst_mem_of_mem_walkSteps h))
      · exact h2 z hz hzi

lemma exists_mem_walkSteps {l : List ℝ} {t : ℝ} (ht : t ∈ l) :
    ∃ o, (t, o) ∈ walkSteps l := by
  induction l using walkSteps.induct with
  | case1 => simp at ht
  | case2 a => exact ⟨none, by simp_all [walkSteps]⟩
  | case3 a b rest ih =>
    rcases List.mem_cons.mp ht with rfl | ht
    · exact ⟨some b, by simp [walkSteps]⟩
    · obtain ⟨o, ho⟩ := ih ht
      exact ⟨o, by simp only [walkSteps, List.mem_cons]; exact Or.inr ho⟩

lemma mem_walkSteps_of_adjacent {l : List ℝ} (hl : l.Pairwise (· < ·)) {u w : ℝ}
    (hu : u ∈ l) (hw : w ∈ l) (huw : u < w) (hgap : ∀ z ∈ l, z ∉ Ioo u w) :
    (u, some w) ∈ walkSteps l := by
  induction l using walkSteps.induct with
  | case1 => simp at hu
  | case2 a =>
    simp only [List.mem_singleton] at hu hw
    exact absurd huw (by rw [hu, hw]; exact lt_irrefl a)
  | case3 a b rest ih =>
    have hl' := List.pairwise_cons.mp hl
    simp only [walkSteps, List.mem_cons, Prod.mk.injEq, Option.some.injEq]
    rcases List.mem_cons.mp hu with rfl | hu'
    · left
      refine ⟨rfl, ?_⟩
      rcases List.mem_cons.mp hw with rfl | hw'
      · exact absurd huw (lt_irrefl _)
      rcases List.mem_cons.mp hw' with rfl | hw''
      · rfl
      · exact absurd ⟨hl'.1 b (by simp), (List.pairwise_cons.mp hl'.2).1 w hw''⟩
          (hgap b (by simp))
    · right
      have hw' : w ∈ b :: rest := by
        rcases List.mem_cons.mp hw with rfl | hw'
        · exact absurd (hl'.1 u hu') (not_lt.mpr huw.le)
        · exact hw'
      exact ih hl'.2 hu' hw' fun z hz => hgap z (List.mem_cons_of_mem _ hz)

def flagStep (kv kd vv : Bool) (s : Bool × Bool) : Bool × Bool :=
  if !(s.1 || kv) && vv then (true, true) else (s.1 || kv, s.2 || kd)

def runFlags (steps : List (ℝ × Option ℝ)) (kv kd vv : ℝ × Option ℝ → Bool) : Bool × Bool :=
  steps.foldl (fun s e => flagStep (kv e) (kd e) (vv e) s) (false, false)

def verdict (s : Bool × Bool) : ℕ := if s.2 = false then 0 else if s.1 = true then 0 else 1

lemma flagStep_fst (kv kd vv : Bool) (s : Bool × Bool) :
    (flagStep kv kd vv s).1 = (s.1 || kv || vv) := by
  obtain ⟨g, d⟩ := s
  cases g <;> cases kv <;> cases vv <;> rfl

lemma flagStep_snd (kv kd vv : Bool) (s : Bool × Bool) (h : (flagStep kv kd vv s).1 = false) :
    (flagStep kv kd vv s).2 = (s.2 || kd) := by
  obtain ⟨g, d⟩ := s
  cases g <;> cases kv <;> cases vv <;> simp_all [flagStep]

lemma foldl_flagStep_fst (kv kd vv : ℝ × Option ℝ → Bool) (L : List (ℝ × Option ℝ))
    (s : Bool × Bool) :
    (L.foldl (fun s e => flagStep (kv e) (kd e) (vv e) s) s).1 = false ↔
      s.1 = false ∧ ∀ e ∈ L, kv e = false ∧ vv e = false := by
  induction L generalizing s with
  | nil => simp
  | cons e L ih =>
    rw [List.foldl_cons, ih, flagStep_fst]
    simp only [Bool.or_eq_false_iff, List.mem_cons, forall_eq_or_imp]
    tauto

lemma foldl_flagStep_snd (kv kd vv : ℝ × Option ℝ → Bool) (L : List (ℝ × Option ℝ))
    (s : Bool × Bool) (h : (L.foldl (fun s e => flagStep (kv e) (kd e) (vv e) s) s).1 = false) :
    (L.foldl (fun s e => flagStep (kv e) (kd e) (vv e) s) s).2 = true ↔
      s.2 = true ∨ ∃ e ∈ L, kd e = true := by
  induction L generalizing s with
  | nil => simp
  | cons e L ih =>
    rw [List.foldl_cons] at h ⊢
    have h1 := ((foldl_flagStep_fst kv kd vv L _).mp h).1
    rw [ih _ h, flagStep_snd _ _ _ _ h1]
    simp only [Bool.or_eq_true, List.mem_cons, exists_eq_or_imp]
    tauto

theorem verdict_runFlags_eq_one (L : List (ℝ × Option ℝ)) (kv kd vv : ℝ × Option ℝ → Bool) :
    verdict (runFlags L kv kd vv) = 1 ↔
      (∀ e ∈ L, kv e = false ∧ vv e = false) ∧ ∃ e ∈ L, kd e = true := by
  have key : verdict (runFlags L kv kd vv) = 1 ↔
      (runFlags L kv kd vv).1 = false ∧ (runFlags L kv kd vv).2 = true := by
    unfold verdict
    cases (runFlags L kv kd vv).1 <;> cases (runFlags L kv kd vv).2 <;> simp
  rw [key]
  constructor
  · rintro ⟨h1, h2⟩
    refine ⟨((foldl_flagStep_fst kv kd vv L _).mp h1).2, ?_⟩
    simpa using (foldl_flagStep_snd kv kd vv L _ h1).mp h2
  · rintro ⟨h1, h2⟩
    have h1' : (runFlags L kv kd vv).1 = false := (foldl_flagStep_fst kv kd vv L _).mpr ⟨rfl, h1⟩
    exact ⟨h1', (foldl_flagStep_snd kv kd vv L _ h1').mpr (Or.inr h2)⟩

def fsdRatio (x : Fin n → ℝ) (t : ℝ) : ℝ := (walkCount x t : ℝ) / n

def fsdKv (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) : Bool :=
  decide (fsdRatio x e.1 - fsdRatio y e.1 > 0)

def fsdKd (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) : Bool :=
  decide (|fsdRatio x e.1 - fsdRatio y e.1| > 0)

def ssdDiff (x : Fin n → ℝ) (y : Fin m → ℝ) (t : ℝ) : ℝ :=
  lpmDeg1 x (walkCount x t) t - lpmDeg1 y (walkCount y t) t

def ssdKv (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) : Bool :=
  decide (ssdDiff x y e.1 > 0)

def ssdKd (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) : Bool :=
  decide (|ssdDiff x y e.1| > 0)

def tsdDiff (x : Fin n → ℝ) (y : Fin m → ℝ) (t : ℝ) : ℝ :=
  lpmDeg2 x (walkCount x t) t - lpmDeg2 y (walkCount y t) t

def tsdKv (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) : Bool :=
  decide (tsdDiff x y e.1 > 0)

def tsdKd (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) : Bool :=
  decide (|tsdDiff x y e.1| > 0)

def tsdVv (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) : Bool :=
  match e.2 with
  | none => false
  | some tn =>
    let kx := walkCount x e.1
    let ky := walkCount y e.1
    let a := (kx : ℝ) / n - (ky : ℝ) / m
    let b := prefix1 x kx / n - prefix1 y ky / m
    decide (a < 0 ∧ (b / a > e.1 ∧ b / a < tn) ∧
      lpmDeg2 x kx (b / a) - lpmDeg2 y ky (b / a) > 0)

def omegaRatio (x : Fin n → ℝ) (t : ℝ) : ℝ :=
  let k := walkCount x t
  if lpmDeg1 x k t + upmDeg1 x k t > 0 then
    lpmDeg1 x k t / (lpmDeg1 x k t + upmDeg1 x k t)
  else (k : ℝ) / n

def ratioCoefs (x : Fin n → ℝ) (k : ℕ) : ℝ × ℝ × ℝ × ℝ :=
  ((k : ℝ) / n, -prefix1 x k / n, (2 * (k : ℝ) - n) / n, (total1 x - 2 * prefix1 x k) / n)

def omegaKv (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) : Bool :=
  decide (omegaRatio x e.1 - omegaRatio y e.1 > 0)

def omegaKd (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) : Bool :=
  decide (|omegaRatio x e.1 - omegaRatio y e.1| > 0)

def stationaryFires (ax bx gx dx ay by' gy dy t tn : ℝ) : Prop :=
  let Ax := ax * dx - bx * gx
  let Ay := ay * dy - by' * gy
  let r := √(Ay / Ax)
  let den := gy - r * gx
  let ts := (r * dx - dy) / den
  let Dx := gx * ts + dx
  let Dy := gy * ts + dy
  (Ax ≠ 0 ∧ Ay ≠ 0 ∧ (Ax > 0 ↔ Ay > 0)) ∧ den ≠ 0 ∧ (ts > t ∧ ts < tn) ∧
    (Dx > 0 ∧ Dy > 0) ∧ (ax * ts + bx) / Dx - (ay * ts + by') / Dy > 0

open scoped Classical in
def omegaSv (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) : Bool :=
  match e.2 with
  | none => false
  | some tn =>
    let cx := ratioCoefs x (walkCount x e.1)
    let cy := ratioCoefs y (walkCount y e.1)
    decide (stationaryFires cx.1 cx.2.1 cx.2.2.1 cx.2.2.2 cy.1 cy.2.1 cy.2.2.1 cy.2.2.2 e.1 tn)

def omegaMid (x : Fin n → ℝ) (y : Fin m → ℝ) (t tn : ℝ) : ℝ :=
  let s := (t + tn) / 2
  let kx := walkCount x t
  let ky := walkCount y t
  lpmDeg1 x kx s / (lpmDeg1 x kx s + upmDeg1 x kx s) -
    lpmDeg1 y ky s / (lpmDeg1 y ky s + upmDeg1 y ky s)

def omegaStep (x : Fin n → ℝ) (y : Fin m → ℝ) (s : Bool × Bool) (e : ℝ × Option ℝ) :
    Bool × Bool :=
  let s1 := flagStep (omegaKv x y e) (omegaKd x y e) (omegaSv x y e) s
  match e.2 with
  | none => s1
  | some tn =>
    if s1.2 = false then
      (s1.1 || decide (omegaMid x y e.1 tn > 0), s1.2 || decide (|omegaMid x y e.1 tn| > 0))
    else s1

def sdDomPair (x : Fin n → ℝ) (y : Fin m → ℝ) (degree : ℕ) (discrete : Bool) : ℕ :=
  let steps := walkSteps (knotList x y)
  if degree = 1 then
    if ¬ sampleMin x ≥ sampleMin y then 0
    else if identicalSamples x y then 0
    else if discrete then verdict (runFlags steps (fsdKv x y) (fsdKd x y) fun _ => false)
    else verdict (steps.foldl (omegaStep x y) (false, false))
  else if ¬ sampleMin x ≥ sampleMin y ∨ sampleMean y - sampleMean x > 0 then 0
  else if identicalSamples x y then 0
  else if degree = 2 then verdict (runFlags steps (ssdKv x y) (ssdKd x y) fun _ => false)
  else verdict (runFlags steps (tsdKv x y) (tsdKd x y) (tsdVv x y))

def tsdKnotOnly (x : Fin n → ℝ) (y : Fin m → ℝ) : ℕ :=
  if ¬ sampleMin x ≥ sampleMin y ∨ sampleMean y - sampleMean x > 0 then 0
  else if identicalSamples x y then 0
  else verdict (runFlags (walkSteps (knotList x y)) (tsdKv x y) (tsdKd x y) fun _ => false)

lemma snd_mem_of_mem_walkSteps {l : List ℝ} {t w : ℝ} (h : (t, some w) ∈ walkSteps l) :
    w ∈ l := by
  induction l using walkSteps.induct with
  | case1 => simp [walkSteps] at h
  | case2 a => simp [walkSteps] at h
  | case3 a b rest ih =>
    simp only [walkSteps, List.mem_cons, Prod.mk.injEq, Option.some.injEq] at h
    rcases h with ⟨-, rfl⟩ | h
    · simp
    · exact List.mem_cons_of_mem _ (ih h)

lemma gap_of_mem_walkSteps {x : Fin n → ℝ} {y : Fin m → ℝ} {u w : ℝ}
    (h : (u, some w) ∈ walkSteps (knotList x y)) :
    IsKnot x y u ∧ IsKnot x y w ∧ u < w ∧ (∀ i, x i ∉ Ioo u w) ∧ ∀ j, y j ∉ Ioo u w := by
  obtain ⟨huw, hgap⟩ := adjacent_of_mem_walkSteps (knotList_pairwise x y) h
  exact ⟨mem_knotList.mp (fst_mem_of_mem_walkSteps h), mem_knotList.mp
    (snd_mem_of_mem_walkSteps h), huw, fun i => hgap _ (mem_knotList.mpr (Or.inl ⟨i, rfl⟩)),
    fun j => hgap _ (mem_knotList.mpr (Or.inr ⟨j, rfl⟩))⟩

lemma forall_mem_of_walkSteps_chain {l : List ℝ} (f : ℝ → Prop)
    (hstep : ∀ u w, (u, some w) ∈ walkSteps l → f u → f w) (h0 : ∀ a ∈ l.head?, f a) :
    ∀ t ∈ l, f t := by
  induction l using walkSteps.induct with
  | case1 => simp
  | case2 a => simpa using h0
  | case3 a b rest ih =>
    have ha : f a := h0 a (by simp)
    have hb : f b := hstep a b (by simp [walkSteps]) ha
    have ih' := ih (fun u w h => hstep u w (by
      simp only [walkSteps, List.mem_cons]; exact Or.inr h)) (by simpa using hb)
    intro t ht
    rcases List.mem_cons.mp ht with rfl | ht
    · exact ha
    · exact ih' t ht

lemma isKnot_symm {x : Fin n → ℝ} {y : Fin m → ℝ} {t : ℝ} (h : IsKnot y x t) : IsKnot x y t :=
  Or.symm h

lemma LPM_succ_eq_zero_of_le {d : ℕ} {t : ℝ} {x : Fin n → ℝ} (h : ∀ i, t ≤ x i) :
    LPM (d + 1) t x = 0 := by
  unfold LPM
  rw [sum_eq_zero fun i _ => ?_, zero_div]
  unfold lowerDev
  split_ifs with hi
  · rw [le_antisymm hi (h i), sub_self, zero_pow (Nat.succ_ne_zero d)]
  · rfl

lemma LPM_two_of_forall_le (hn : 0 < n) {t : ℝ} {x : Fin n → ℝ} (h : ∀ i, x i ≤ t) :
    LPM 2 t x = t ^ 2 - 2 * t * mean x + (∑ i, x i ^ 2) / n := by
  have hn' : (n : ℝ) ≠ 0 := by positivity
  unfold LPM mean
  rw [sum_congr rfl fun i _ => (show lowerDev 2 t (x i) = t ^ 2 - 2 * t * x i + x i ^ 2 by
    simp only [lowerDev, h i, ↓reduceIte]; ring)]
  rw [sum_add_distrib, sum_sub_distrib, sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul,
    ← mul_sum]
  field_simp

theorem mean_le_of_LPM_two_le (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ} {y : Fin m → ℝ}
    (h : ∀ t, LPM 2 t x ≤ LPM 2 t y) : mean y ≤ mean x := by
  by_contra hlt
  push Not at hlt
  obtain ⟨bx, hbx⟩ := exists_sample_upperBound x
  obtain ⟨by', hby⟩ := exists_sample_upperBound y
  set Qx := (∑ i, x i ^ 2) / n
  set Qy := (∑ j, y j ^ 2) / m
  set δ := mean y - mean x with hδ
  have hδpos : 0 < δ := sub_pos.mpr hlt
  set t := max (max bx by') ((Qy - Qx) / (2 * δ) + 1)
  have ht := h t
  rw [LPM_two_of_forall_le hn fun i => (hbx i).trans ((le_max_left _ _).trans (le_max_left _ _)),
    LPM_two_of_forall_le hm fun j => (hby j).trans ((le_max_right _ _).trans (le_max_left _ _))]
    at ht
  have htt : (Qy - Qx) / (2 * δ) + 1 ≤ t := le_max_right _ _
  have key : 2 * δ * ((Qy - Qx) / (2 * δ) + 1) = Qy - Qx + 2 * δ := by
    field_simp
  nlinarith [mul_le_mul_of_nonneg_left htt (by positivity : (0 : ℝ) ≤ 2 * δ)]

lemma LPM_two_sub_eq {u w : ℝ} {x : Fin n → ℝ} (huw : u ≤ w) (hx : ∀ i, x i ∉ Ioo u w) :
    LPM 2 w x - LPM 2 u x = (w - u) * (LPM 1 u x + LPM 1 w x) := by
  rw [LPM_two_eq_on_Icc hx ⟨huw, le_rfl⟩, LPM_two_eq_on_Icc hx ⟨le_rfl, huw⟩,
    LPM_one_eq_on_Icc hx ⟨le_rfl, huw⟩, LPM_one_eq_on_Icc hx ⟨huw, le_rfl⟩]
  ring

theorem LPM_one_eq_of_eq_at_knots (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ} {y : Fin m → ℝ}
    (h : ∀ t, IsKnot x y t → LPM 1 t x = LPM 1 t y) (t : ℝ) : LPM 1 t x = LPM 1 t y :=
  le_antisymm (LPM_one_le_of_le_at_knots hm (fun t ht => (h t ht).le) t)
    (LPM_one_le_of_le_at_knots hn (fun t ht => (h t (isKnot_symm ht)).ge) t)

theorem LPM_two_eq_of_eq_at_knots (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ} {y : Fin m → ℝ}
    (h : ∀ t, IsKnot x y t → LPM 2 t x = LPM 2 t y) (t : ℝ) : LPM 2 t x = LPM 2 t y := by
  have h1 : ∀ t ∈ knotList x y, LPM 1 t x = LPM 1 t y := by
    refine forall_mem_of_walkSteps_chain (fun t => LPM 1 t x = LPM 1 t y) ?_ ?_
    · intro u w hstep hu
      obtain ⟨ku, kw, huw, hx, hy⟩ := gap_of_mem_walkSteps hstep
      have ex := LPM_two_sub_eq huw.le hx
      have ey := LPM_two_sub_eq huw.le hy
      rw [h u ku, h w kw] at ex
      have hne : w - u ≠ 0 := sub_ne_zero.mpr huw.ne'
      have : (w - u) * (LPM 1 w x - LPM 1 w y) = 0 := by
        rw [hu] at ex
        linear_combination ey - ex
      simpa [hne, sub_eq_zero] using this
    · intro a ha
      have hmin : ∀ v ∈ knotList x y, a ≤ v := by
        intro v hv
        cases hl : knotList x y with
        | nil => rw [hl] at hv; simp at hv
        | cons b rest =>
          rw [hl] at ha hv
          have hp := knotList_pairwise x y
          rw [hl, List.pairwise_cons] at hp
          simp only [List.head?_cons, Option.mem_def, Option.some.injEq] at ha
          subst ha
          rcases List.mem_cons.mp hv with rfl | hv
          · exact le_rfl
          · exact (hp.1 v hv).le
      change LPM (0 + 1) a x = LPM (0 + 1) a y
      rw [LPM_succ_eq_zero_of_le fun i => hmin _ (mem_knotList.mpr (Or.inl ⟨i, rfl⟩)),
        LPM_succ_eq_zero_of_le fun j => hmin _ (mem_knotList.mpr (Or.inr ⟨j, rfl⟩))]
  have h1' : ∀ s, LPM 1 s x = LPM 1 s y :=
    LPM_one_eq_of_eq_at_knots hn hm fun t ht => h1 t (mem_knotList.mpr ht)
  obtain ⟨ax, hax⟩ := exists_sample_lowerBound x
  obtain ⟨ay, hay⟩ := exists_sample_lowerBound y
  rw [LPM_succ_eq_integral 1 (a := min ax ay) (fun i => (min_le_left _ _).trans (hax i)) t,
    LPM_succ_eq_integral 1 (a := min ax ay) (fun j => (min_le_right _ _).trans (hay j)) t]
  simp_rw [h1']

lemma LPM_zero_eq_of_iff {t u : ℝ} {x : Fin n → ℝ} (h : ∀ i, x i ≤ t ↔ x i ≤ u) :
    LPM 0 t x = LPM 0 u x := by
  rw [LPM_zero_eq_ecdf, LPM_zero_eq_ecdf, filter_congr fun i _ => h i]

lemma LPM_zero_cases (x : Fin n → ℝ) (y : Fin m → ℝ) (t : ℝ) :
    (LPM 0 t x = 0 ∧ LPM 0 t y = 0) ∨
      ∃ u, IsKnot x y u ∧ LPM 0 t x = LPM 0 u x ∧ LPM 0 t y = LPM 0 u y := by
  set K := Finset.univ.image x ∪ Finset.univ.image y
  by_cases hne : (K.filter (· ≤ t)).Nonempty
  · right
    set u := (K.filter (· ≤ t)).max' hne
    have hu := mem_filter.mp ((K.filter (· ≤ t)).max'_mem hne)
    have hiff : ∀ v ∈ K, v ≤ t ↔ v ≤ u := fun v hv =>
      ⟨fun hvt => le_max' _ v (mem_filter.mpr ⟨hv, hvt⟩), fun hvu => hvu.trans hu.2⟩
    refine ⟨u, isKnot_iff_mem.mpr hu.1, LPM_zero_eq_of_iff fun i => hiff _ ?_,
      LPM_zero_eq_of_iff fun j => hiff _ ?_⟩
    · exact isKnot_iff_mem.mp (Or.inl ⟨i, rfl⟩)
    · exact isKnot_iff_mem.mp (Or.inr ⟨j, rfl⟩)
  · left
    have hlt : ∀ v ∈ K, t < v := fun v hv => not_le.mp fun hvt => hne ⟨v, mem_filter.mpr ⟨hv, hvt⟩⟩
    exact ⟨LPM_eq_zero_of_lt fun i => hlt _ (isKnot_iff_mem.mp (Or.inl ⟨i, rfl⟩)),
      LPM_eq_zero_of_lt fun j => hlt _ (isKnot_iff_mem.mp (Or.inr ⟨j, rfl⟩))⟩

theorem LPM_zero_le_of_le_at_knots {x : Fin n → ℝ} {y : Fin m → ℝ}
    (h : ∀ t, IsKnot x y t → LPM 0 t x ≤ LPM 0 t y) (t : ℝ) : LPM 0 t x ≤ LPM 0 t y := by
  rcases LPM_zero_cases x y t with ⟨hx, hy⟩ | ⟨u, hu, hx, hy⟩
  · rw [hx, hy]
  · rw [hx, hy]; exact h u hu

theorem LPM_zero_eq_of_eq_at_knots {x : Fin n → ℝ} {y : Fin m → ℝ}
    (h : ∀ t, IsKnot x y t → LPM 0 t x = LPM 0 t y) (t : ℝ) : LPM 0 t x = LPM 0 t y := by
  rcases LPM_zero_cases x y t with ⟨hx, hy⟩ | ⟨u, hu, hx, hy⟩
  · rw [hx, hy]
  · rw [hx, hy]; exact h u hu

lemma fsdKv_eq_false (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) :
    fsdKv x y e = false ↔ LPM 0 e.1 x ≤ LPM 0 e.1 y := by
  simp [fsdKv, fsdRatio, walkCount_div_eq_LPM_zero]

lemma fsdKd_eq_true (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) :
    fsdKd x y e = true ↔ LPM 0 e.1 x ≠ LPM 0 e.1 y := by
  simp [fsdKd, fsdRatio, walkCount_div_eq_LPM_zero, sub_eq_zero]

lemma ssdDiff_eq (x : Fin n → ℝ) (y : Fin m → ℝ) (t : ℝ) :
    ssdDiff x y t = LPM 1 t x - LPM 1 t y := by
  rw [ssdDiff, lpmDeg1_walkCount_self, lpmDeg1_walkCount_self]

lemma ssdKv_eq_false (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) :
    ssdKv x y e = false ↔ LPM 1 e.1 x ≤ LPM 1 e.1 y := by
  simp [ssdKv, ssdDiff_eq]

lemma ssdKd_eq_true (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) :
    ssdKd x y e = true ↔ LPM 1 e.1 x ≠ LPM 1 e.1 y := by
  simp [ssdKd, ssdDiff_eq, sub_eq_zero]

lemma tsdDiff_eq (x : Fin n → ℝ) (y : Fin m → ℝ) (t : ℝ) :
    tsdDiff x y t = LPM 2 t x - LPM 2 t y := by
  rw [tsdDiff, lpmDeg2_walkCount_self, lpmDeg2_walkCount_self]

lemma tsdKv_eq_false (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) :
    tsdKv x y e = false ↔ LPM 2 e.1 x ≤ LPM 2 e.1 y := by
  simp [tsdKv, tsdDiff_eq]

lemma tsdKd_eq_true (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) :
    tsdKd x y e = true ↔ LPM 2 e.1 x ≠ LPM 2 e.1 y := by
  simp [tsdKd, tsdDiff_eq, sub_eq_zero]

lemma tsdVv_eq_false {x : Fin n → ℝ} {y : Fin m → ℝ} {e : ℝ × Option ℝ}
    (he : e ∈ walkSteps (knotList x y)) :
    tsdVv x y e = false ↔ ∀ tn, e.2 = some tn → tsdCoefA e.1 x y < 0 →
      tsdCoefB e.1 x y / tsdCoefA e.1 x y ∈ Ioo e.1 tn →
      LPM 2 (tsdCoefB e.1 x y / tsdCoefA e.1 x y) x ≤
        LPM 2 (tsdCoefB e.1 x y / tsdCoefA e.1 x y) y := by
  obtain ⟨t, o⟩ := e
  cases o with
  | none => simp [tsdVv]
  | some tn =>
    obtain ⟨-, -, -, hx, hy⟩ := gap_of_mem_walkSteps he
    have ha : (walkCount x t : ℝ) / n - (walkCount y t : ℝ) / m = tsdCoefA t x y := by
      simp [tsdCoefA, walkCount_eq_lowerCount]
    have hb : prefix1 x (walkCount x t) / n - prefix1 y (walkCount y t) / m =
        tsdCoefB t x y := by
      rw [prefix1_walkCount, prefix1_walkCount]; rfl
    simp only [tsdVv, ha, hb, decide_eq_false_iff_not, Option.some.injEq, forall_eq', Set.mem_Ioo]
    constructor
    · intro h hA hin
      rw [← lpmDeg2_walkCount hx ⟨hin.1.le, hin.2.le⟩, ← lpmDeg2_walkCount hy ⟨hin.1.le, hin.2.le⟩]
      by_contra hc
      exact h ⟨hA, hin, by linarith [not_le.mp hc]⟩
    · rintro h ⟨hA, hin, hD⟩
      have := h hA hin
      rw [← lpmDeg2_walkCount hx ⟨hin.1.le, hin.2.le⟩,
        ← lpmDeg2_walkCount hy ⟨hin.1.le, hin.2.le⟩] at this
      linarith

lemma sdDomPair_fsd_branch (x : Fin n → ℝ) (y : Fin m → ℝ) :
    sdDomPair x y 1 true = 1 ↔ sampleMin y ≤ sampleMin x ∧ ¬ identicalSamples x y ∧
      verdict (runFlags (walkSteps (knotList x y)) (fsdKv x y) (fsdKd x y) fun _ => false) = 1 := by
  unfold sdDomPair
  split_ifs <;> simp_all

lemma sdDomPair_ssd_tsd_branch (x : Fin n → ℝ) (y : Fin m → ℝ) {d : ℕ} (hd : d ≠ 1)
    (b : Bool) :
    sdDomPair x y d b = 1 ↔ sampleMin y ≤ sampleMin x ∧ ¬ sampleMean y - sampleMean x > 0 ∧
      ¬ identicalSamples x y ∧
      (if d = 2 then verdict (runFlags (walkSteps (knotList x y)) (ssdKv x y) (ssdKd x y)
          fun _ => false)
        else verdict (runFlags (walkSteps (knotList x y)) (tsdKv x y) (tsdKd x y)
          (tsdVv x y))) = 1 := by
  unfold sdDomPair
  simp only [hd, ↓reduceIte]
  by_cases hg : sampleMin y ≤ sampleMin x ∧ ¬ sampleMean y - sampleMean x > 0
  · rw [ite_eq_right (by tauto)]
    simp only [hg, not_false_eq_true, true_and]
    split_ifs <;> simp_all
  · rw [ite_eq_left (by tauto)]
    simp only [zero_ne_one, false_iff]
    tauto

theorem sdDomPair_fsd_iff (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    sdDomPair x y 1 true = 1 ↔
      (∀ t, LPM 0 t x ≤ LPM 0 t y) ∧ ∃ t, LPM 0 t x ≠ LPM 0 t y := by
  rw [sdDomPair_fsd_branch, verdict_runFlags_eq_one]
  constructor
  · rintro ⟨-, -, hall, e, -, hkd⟩
    refine ⟨LPM_zero_le_of_le_at_knots fun t ht => ?_, e.1, (fsdKd_eq_true x y e).mp hkd⟩
    obtain ⟨o, ho⟩ := exists_mem_walkSteps (mem_knotList.mpr ht)
    exact (fsdKv_eq_false x y _).mp (hall _ ho).1
  · rintro ⟨hle, t, hne⟩
    refine ⟨(sampleMin_le_iff hn hm x y).mpr (exists_le_of_forall_LPM_le hle),
      fun hid => hne (LPM_eq_of_identical hid 0 t),
      fun e _ => ⟨(fsdKv_eq_false x y e).mpr (hle _), rfl⟩, ?_⟩
    by_contra hcon
    push Not at hcon
    refine hne (LPM_zero_eq_of_eq_at_knots (fun u hu => ?_) t)
    obtain ⟨o, ho⟩ := exists_mem_walkSteps (mem_knotList.mpr hu)
    have := hcon _ ho
    rw [ne_eq, fsdKd_eq_true] at this
    exact not_not.mp this

theorem sdDomPair_ssd_iff (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    (b : Bool) :
    sdDomPair x y 2 b = 1 ↔
      (∀ t, LPM 1 t x ≤ LPM 1 t y) ∧ ∃ t, LPM 1 t x ≠ LPM 1 t y := by
  rw [sdDomPair_ssd_tsd_branch x y (by norm_num) b]
  simp only [↓reduceIte]
  rw [verdict_runFlags_eq_one]
  constructor
  · rintro ⟨-, -, -, hall, e, -, hkd⟩
    refine ⟨LPM_one_le_of_le_at_knots hm fun t ht => ?_, e.1, (ssdKd_eq_true x y e).mp hkd⟩
    obtain ⟨o, ho⟩ := exists_mem_walkSteps (mem_knotList.mpr ht)
    exact (ssdKv_eq_false x y _).mp (hall _ ho).1
  · rintro ⟨hle, t, hne⟩
    refine ⟨(sampleMin_le_iff hn hm x y).mpr (exists_le_of_forall_LPM_le hle), ?_,
      fun hid => hne (LPM_eq_of_identical hid 1 t),
      fun e _ => ⟨(ssdKv_eq_false x y e).mpr (hle _), rfl⟩, ?_⟩
    · rw [sampleMean_eq_mean, sampleMean_eq_mean]
      linarith [mean_le_of_LPM_one_le hn hm hle]
    by_contra hcon
    push Not at hcon
    refine hne (LPM_one_eq_of_eq_at_knots hn hm (fun u hu => ?_) t)
    obtain ⟨o, ho⟩ := exists_mem_walkSteps (mem_knotList.mpr hu)
    have := hcon _ ho
    rw [ne_eq, ssdKd_eq_true] at this
    exact not_not.mp this

theorem sdDomPair_tsd_iff (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    (b : Bool) :
    sdDomPair x y 3 b = 1 ↔
      (∀ t, LPM 2 t x ≤ LPM 2 t y) ∧ ∃ t, LPM 2 t x ≠ LPM 2 t y := by
  rw [sdDomPair_ssd_tsd_branch x y (by norm_num) b]
  simp only [show (3 : ℕ) ≠ 2 by norm_num, ↓reduceIte]
  rw [verdict_runFlags_eq_one]
  constructor
  · rintro ⟨-, hmean, -, hall, e, -, hkd⟩
    have hknot : ∀ t, IsKnot x y t → LPM 2 t x ≤ LPM 2 t y := fun t ht => by
      obtain ⟨o, ho⟩ := exists_mem_walkSteps (mem_knotList.mpr ht)
      exact (tsdKv_eq_false x y _).mp (hall _ ho).1
    have hmean' : mean y ≤ mean x := by
      rw [← sampleMean_eq_mean, ← sampleMean_eq_mean]
      linarith [not_lt.mp hmean]
    refine ⟨LPM_two_le_of_le_at_knots hm hmean' hknot ?_, e.1, (tsdKd_eq_true x y e).mp hkd⟩
    intro u w hu hw hx hy hA hin
    have hstep : (u, some w) ∈ walkSteps (knotList x y) :=
      mem_walkSteps_of_adjacent (knotList_pairwise x y) (mem_knotList.mpr hu)
        (mem_knotList.mpr hw) (hin.1.trans hin.2) fun z hz hzi => by
          rcases mem_knotList.mp hz with ⟨i, rfl⟩ | ⟨j, rfl⟩
          · exact hx i hzi
          · exact hy j hzi
    exact (tsdVv_eq_false hstep).mp (hall _ hstep).2 w rfl hA hin
  · rintro ⟨hle, t, hne⟩
    refine ⟨(sampleMin_le_iff hn hm x y).mpr (exists_le_of_forall_LPM_le hle), ?_,
      fun hid => hne (LPM_eq_of_identical hid 2 t),
      fun e he => ⟨(tsdKv_eq_false x y e).mpr (hle _),
        (tsdVv_eq_false he).mpr fun _ _ _ _ => hle _⟩, ?_⟩
    · rw [sampleMean_eq_mean, sampleMean_eq_mean]
      linarith [mean_le_of_LPM_two_le hn hm hle]
    by_contra hcon
    push Not at hcon
    refine hne (LPM_two_eq_of_eq_at_knots hn hm (fun u hu => ?_) t)
    obtain ⟨o, ho⟩ := exists_mem_walkSteps (mem_knotList.mpr hu)
    have := hcon _ ho
    rw [ne_eq, tsdKd_eq_true] at this
    exact not_not.mp this

theorem sdDomPair_eq_zero_or_one (x : Fin n → ℝ) (y : Fin m → ℝ) (d : ℕ) (b : Bool) :
    sdDomPair x y d b = 0 ∨ sdDomPair x y d b = 1 := by
  have hv : ∀ s : Bool × Bool, verdict s = 0 ∨ verdict s = 1 := fun s => by
    unfold verdict; split_ifs <;> simp
  unfold sdDomPair
  split_ifs <;> first | exact Or.inl rfl | exact hv _

def sampleMax (x : Fin n → ℝ) : ℝ := (sortedVals x).getLastD 0

def shiftSample (x : Fin n → ℝ) (c : ℝ) : Fin n → ℝ := fun i => x i - c

def sdDomPairCentred (x : Fin n → ℝ) (y : Fin m → ℝ) (degree : ℕ) (discrete : Bool) : ℕ :=
  if degree = 1 ∧ discrete = true then sdDomPair x y degree discrete
  else if ¬ sampleMin x ≥ sampleMin y ∨ identicalSamples x y then 0
  else if degree ≥ 2 ∧ sampleMean y - sampleMean x > 0 then 0
  else
    let c := min (sampleMin x) (sampleMin y) / 2 + max (sampleMax x) (sampleMax y) / 2
    sdDomPair (shiftSample x c) (shiftSample y c) degree discrete

def LPMDominates (d : ℕ) (x : Fin n → ℝ) (y : Fin m → ℝ) : Prop :=
  (∀ t, LPM d t x ≤ LPM d t y) ∧ ∃ t, LPM d t x ≠ LPM d t y

lemma LPM_shiftSample (d : ℕ) (t c : ℝ) (x : Fin n → ℝ) :
    LPM d t (shiftSample x c) = LPM d (t + c) x := by
  unfold LPM lowerDev shiftSample
  congr 1
  refine sum_congr rfl fun i _ => ?_
  by_cases h : x i ≤ t + c
  · rw [ite_eq_left (by linarith), ite_eq_left h]; ring_nf
  · rw [ite_eq_right (by intro h'; exact h (by linarith)), ite_eq_right h]

theorem LPMDominates_shift (d : ℕ) (c : ℝ) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    LPMDominates d (shiftSample x c) (shiftSample y c) ↔ LPMDominates d x y := by
  unfold LPMDominates
  simp only [LPM_shiftSample]
  constructor
  · rintro ⟨h1, t, h2⟩
    exact ⟨fun s => by simpa using h1 (s - c), t + c, h2⟩
  · rintro ⟨h1, t, h2⟩
    exact ⟨fun s => h1 (s + c), t - c, by simpa using h2⟩

lemma sdDomPair_eq_one_iff_dominates (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ)
    (y : Fin m → ℝ) {d : ℕ} (hd : d = 2 ∨ d = 3) (b : Bool) :
    sdDomPair x y d b = 1 ↔ LPMDominates (d - 1) x y := by
  rcases hd with rfl | rfl
  · exact sdDomPair_ssd_iff hn hm x y b
  · exact sdDomPair_tsd_iff hn hm x y b

theorem sdDomPairCentred_eq (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    {d : ℕ} (hd : d = 1 ∨ d = 2 ∨ d = 3) (b : Bool) (hb : d = 1 → b = true) :
    sdDomPairCentred x y d b = sdDomPair x y d b := by
  rcases hd with rfl | hd
  · simp [sdDomPairCentred, hb rfl]
  have hd1 : d ≠ 1 := by rcases hd with rfl | rfl <;> norm_num
  have hd2 : 2 ≤ d := by rcases hd with rfl | rfl <;> norm_num
  have hzero : sdDomPair x y d b ≠ 1 → sdDomPair x y d b = 0 := fun h =>
    (sdDomPair_eq_zero_or_one x y d b).resolve_right h
  unfold sdDomPairCentred
  simp only [hd1, false_and, ↓reduceIte]
  split_ifs with h1 h2
  · refine (hzero fun h => ?_).symm
    rw [sdDomPair_ssd_tsd_branch x y hd1 b] at h
    tauto
  · refine (hzero fun h => ?_).symm
    rw [sdDomPair_ssd_tsd_branch x y hd1 b] at h
    exact h.2.1 h2.2
  · rcases sdDomPair_eq_zero_or_one x y d b with h | h
    · rcases sdDomPair_eq_zero_or_one (shiftSample x _) (shiftSample y _) d b with h' | h'
      · rw [h, h']
      · rw [sdDomPair_eq_one_iff_dominates hn hm _ _ hd b, LPMDominates_shift,
          ← sdDomPair_eq_one_iff_dominates hn hm _ _ hd b] at h'
        omega
    · rw [h]
      rw [sdDomPair_eq_one_iff_dominates hn hm _ _ hd b, ← LPMDominates_shift (d - 1)
        (min (sampleMin x) (sampleMin y) / 2 + max (sampleMax x) (sampleMax y) / 2),
        ← sdDomPair_eq_one_iff_dominates hn hm _ _ hd b] at h
      exact h

theorem tsdKnotOnly_eq_one_iff (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    tsdKnotOnly x y = 1 ↔ (∀ i, ∃ j, y j ≤ x i) ∧ mean y ≤ mean x ∧
      ¬ identicalSamples x y ∧ (∀ t, IsKnot x y t → LPM 2 t x ≤ LPM 2 t y) ∧
      ∃ t, IsKnot x y t ∧ LPM 2 t x ≠ LPM 2 t y := by
  have hbr : tsdKnotOnly x y = 1 ↔ sampleMin y ≤ sampleMin x ∧
      ¬ sampleMean y - sampleMean x > 0 ∧ ¬ identicalSamples x y ∧
      verdict (runFlags (walkSteps (knotList x y)) (tsdKv x y) (tsdKd x y)
        fun _ => false) = 1 := by
    unfold tsdKnotOnly
    by_cases hg : sampleMin y ≤ sampleMin x ∧ ¬ sampleMean y - sampleMean x > 0
    · rw [ite_eq_right (by tauto)]
      simp only [hg, not_false_eq_true, true_and]
      split_ifs <;> simp_all
    · rw [ite_eq_left (by tauto)]
      simp only [zero_ne_one, false_iff]
      tauto
  rw [hbr, verdict_runFlags_eq_one, sampleMin_le_iff hn hm, sampleMean_eq_mean,
    sampleMean_eq_mean, gt_iff_lt, sub_pos, not_lt]
  constructor
  · rintro ⟨h1, h2, h3, hall, e, he, hkd⟩
    refine ⟨h1, h2, h3, fun t ht => ?_, e.1, mem_knotList.mp (fst_mem_of_mem_walkSteps he),
      (tsdKd_eq_true x y e).mp hkd⟩
    obtain ⟨o, ho⟩ := exists_mem_walkSteps (mem_knotList.mpr ht)
    exact (tsdKv_eq_false x y _).mp (hall _ ho).1
  · rintro ⟨h1, h2, h3, hk, t, ht, hne⟩
    refine ⟨h1, h2, h3, fun e he => ⟨(tsdKv_eq_false x y e).mpr
      (hk _ (mem_knotList.mp (fst_mem_of_mem_walkSteps he))), rfl⟩, ?_⟩
    obtain ⟨o, ho⟩ := exists_mem_walkSteps (mem_knotList.mpr ht)
    exact ⟨_, ho, (tsdKd_eq_true x y _).mpr hne⟩

theorem tsdKnotOnly_vertexCex : tsdKnotOnly tsdVertexCexX tsdVertexCexY = 1 := by
  rw [tsdKnotOnly_eq_one_iff (by norm_num) (by norm_num)]
  refine ⟨fun i => ⟨3, ?_⟩, tsdVertexCex_mean_le, fun h => ?_, tsdVertexCex_le_at_knots,
    1, Or.inl ⟨3, rfl⟩, ?_⟩
  · fin_cases i <;> norm_num [tsdVertexCexX, tsdVertexCexY]
  · have := congrArg List.length h
    rw [length_sortedVals, length_sortedVals] at this
    norm_num at this
  · norm_num [LPM, lowerDev, Fin.sum_univ_five, Fin.sum_univ_four, tsdVertexCexX,
      tsdVertexCexY]

theorem sdDomPair_tsd_vertexCex (b : Bool) :
    sdDomPair tsdVertexCexX tsdVertexCexY 3 b = 0 := by
  refine (sdDomPair_eq_zero_or_one _ _ 3 b).resolve_right fun h => ?_
  rw [sdDomPair_tsd_iff (by norm_num) (by norm_num)] at h
  exact absurd (h.1 3) (not_le.mpr tsdVertexCex_not_le_at_three)

theorem upmDeg1_walkCount (hn : 0 < n) {x : Fin n → ℝ} {u w s : ℝ} (hx : ∀ i, x i ∉ Ioo u w)
    (hs : s ∈ Icc u w) : upmDeg1 x (walkCount x u) s = UPM 1 s x := by
  have h1 := UPM_one_sub_LPM_one hn s x
  have h2 := lpmDeg1_walkCount hx hs
  have h3 := sampleMean_eq_mean x
  have hn' : (n : ℝ) ≠ 0 := by positivity
  unfold sampleMean at h3
  unfold lpmDeg1 at h2
  have : UPM 1 s x =
      ((walkCount x u : ℝ) * s - prefix1 x (walkCount x u)) / n + total1 x / n - s := by
    rw [h2, h3]; linarith
  rw [this, upmDeg1]
  field_simp
  ring

theorem upmDeg1_walkCount_self (hn : 0 < n) (x : Fin n → ℝ) (t : ℝ) :
    upmDeg1 x (walkCount x t) t = UPM 1 t x :=
  upmDeg1_walkCount hn (w := t) (fun _ => not_mem_Ioo_self t _) ⟨le_rfl, le_rfl⟩

lemma ratioCoefs_num (x : Fin n → ℝ) (k : ℕ) (s : ℝ) :
    (ratioCoefs x k).1 * s + (ratioCoefs x k).2.1 = lpmDeg1 x k s := by
  simp only [ratioCoefs, lpmDeg1]; ring

lemma ratioCoefs_den (x : Fin n → ℝ) (k : ℕ) (s : ℝ) :
    (ratioCoefs x k).2.2.1 * s + (ratioCoefs x k).2.2.2 = lpmDeg1 x k s + upmDeg1 x k s := by
  simp only [ratioCoefs, lpmDeg1, upmDeg1]; ring

lemma ratio_denom_nonneg (t : ℝ) (x : Fin n → ℝ) : 0 ≤ LPM 1 t x + UPM 1 t x :=
  add_nonneg (LPM_nonneg _ _ _) (UPM_nonneg _ _ _)

lemma LPMratio_one_of_pos {t : ℝ} {x : Fin n → ℝ} (hD : 0 < LPM 1 t x + UPM 1 t x) :
    LPMratio 1 t x = LPM 1 t x / (LPM 1 t x + UPM 1 t x) := by
  unfold LPMratio
  simp only [one_ne_zero, ↓reduceIte]
  refine ite_eq_right ?_
  rintro ⟨-, h⟩
  exact hD.ne' ((ratio_denom_eq_zero_iff one_ne_zero t x).mpr h)

theorem omegaRatio_eq (hn : 0 < n) (x : Fin n → ℝ) (t : ℝ) :
    omegaRatio x t = LPMratio 1 t x := by
  unfold omegaRatio
  simp only
  rw [lpmDeg1_walkCount_self, upmDeg1_walkCount_self hn]
  split_ifs with h
  · exact (LPMratio_one_of_pos h).symm
  · have hD : LPM 1 t x + UPM 1 t x = 0 := le_antisymm (not_lt.mp h) (ratio_denom_nonneg t x)
    have hall := (ratio_denom_eq_zero_iff one_ne_zero t x).mp hD
    rw [LPMratio_of_forall_eq one_ne_zero hn hall, walkCount_div_eq_LPM_zero, LPM_zero_eq_ecdf,
      filter_true_of_mem fun i _ => (hall i).le, card_univ, Fintype.card_fin]
    exact div_self (by positivity)

theorem omega_formula_eq (hn : 0 < n) {x : Fin n → ℝ} {u w s : ℝ} (hx : ∀ i, x i ∉ Ioo u w)
    (hs : s ∈ Icc u w) (hD : 0 < (ratioCoefs x (walkCount x u)).2.2.1 * s +
      (ratioCoefs x (walkCount x u)).2.2.2) :
    ((ratioCoefs x (walkCount x u)).1 * s + (ratioCoefs x (walkCount x u)).2.1) /
      ((ratioCoefs x (walkCount x u)).2.2.1 * s + (ratioCoefs x (walkCount x u)).2.2.2) =
      LPMratio 1 s x := by
  rw [ratioCoefs_den, lpmDeg1_walkCount hx hs, upmDeg1_walkCount hn hx hs] at hD ⊢
  rw [ratioCoefs_num, lpmDeg1_walkCount hx hs, LPMratio_one_of_pos hD]

lemma ratio_denom_pos_of_gap (hn : 0 < n) {x : Fin n → ℝ} {u w s : ℝ}
    (hx : ∀ i, x i ∉ Ioo u w) (hs : s ∈ Ioo u w) : 0 < LPM 1 s x + UPM 1 s x :=
  ratio_denom_pos ⟨⟨0, hn⟩, fun h => hx _ (h ▸ hs)⟩

theorem omegaMid_eq (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ} {y : Fin m → ℝ} {t tn : ℝ}
    (hx : ∀ i, x i ∉ Ioo t tn) (hy : ∀ j, y j ∉ Ioo t tn) (htn : t < tn) :
    omegaMid x y t tn = LPMratio 1 ((t + tn) / 2) x - LPMratio 1 ((t + tn) / 2) y := by
  have hs : (t + tn) / 2 ∈ Ioo t tn := ⟨by linarith, by linarith⟩
  unfold omegaMid
  simp only
  have hsI := Ioo_subset_Icc_self hs
  rw [lpmDeg1_walkCount hx hsI, upmDeg1_walkCount hn hx hsI,
    lpmDeg1_walkCount hy hsI, upmDeg1_walkCount hm hy hsI,
    LPMratio_one_of_pos (ratio_denom_pos_of_gap hn hx hs),
    LPMratio_one_of_pos (ratio_denom_pos_of_gap hm hy hs)]

def omegaCross (x : Fin n → ℝ) (y : Fin m → ℝ) (s : ℝ) : ℝ :=
  LPM 1 s x * (LPM 1 s y + UPM 1 s y) - LPM 1 s y * (LPM 1 s x + UPM 1 s x)

lemma omegaCross_eq_zero {x : Fin n → ℝ} {y : Fin m → ℝ} {s : ℝ}
    (h : LPMratio 1 s x = LPMratio 1 s y) : omegaCross x y s = 0 := by
  unfold omegaCross
  rcases (ratio_denom_nonneg s x).eq_or_lt with hx | hx
  · have hL : LPM 1 s x = 0 := le_antisymm (by linarith [UPM_nonneg 1 s x]) (LPM_nonneg _ _ _)
    rw [← hx, hL]; ring
  rcases (ratio_denom_nonneg s y).eq_or_lt with hy | hy
  · have hL : LPM 1 s y = 0 := le_antisymm (by linarith [UPM_nonneg 1 s y]) (LPM_nonneg _ _ _)
    rw [← hy, hL]; ring
  rw [LPMratio_one_of_pos hx, LPMratio_one_of_pos hy, div_eq_div_iff hx.ne' hy.ne'] at h
  linarith

lemma quad_eq_zero_of_three {A B C u v w : ℝ} (huv : u ≠ v) (hvw : v ≠ w) (huw : u ≠ w)
    (hu : A * u ^ 2 + B * u + C = 0) (hv : A * v ^ 2 + B * v + C = 0)
    (hw : A * w ^ 2 + B * w + C = 0) (s : ℝ) : A * s ^ 2 + B * s + C = 0 := by
  have e1 : (u - v) * (A * (u + v) + B) = 0 := by linear_combination hu - hv
  have e2 : (v - w) * (A * (v + w) + B) = 0 := by linear_combination hv - hw
  have f1 : A * (u + v) + B = 0 := (mul_eq_zero.mp e1).resolve_left (sub_ne_zero.mpr huv)
  have f2 : A * (v + w) + B = 0 := (mul_eq_zero.mp e2).resolve_left (sub_ne_zero.mpr hvw)
  have g : (u - w) * A = 0 := by linear_combination f1 - f2
  have hA : A = 0 := (mul_eq_zero.mp g).resolve_left (sub_ne_zero.mpr huw)
  have hB : B = 0 := by rw [hA] at f1; linarith
  have hC : C = 0 := by rw [hA, hB] at hu; linarith
  rw [hA, hB, hC]; ring

lemma omegaCross_quadratic (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ} {y : Fin m → ℝ}
    {u w : ℝ} (hx : ∀ i, x i ∉ Ioo u w) (hy : ∀ j, y j ∉ Ioo u w) :
    ∃ A B C : ℝ, ∀ s ∈ Icc u w, omegaCross x y s = A * s ^ 2 + B * s + C := by
  set cx := ratioCoefs x (walkCount x u)
  set cy := ratioCoefs y (walkCount y u)
  refine ⟨cx.1 * cy.2.2.1 - cy.1 * cx.2.2.1,
    cx.1 * cy.2.2.2 + cx.2.1 * cy.2.2.1 - cy.1 * cx.2.2.2 - cy.2.1 * cx.2.2.1,
    cx.2.1 * cy.2.2.2 - cy.2.1 * cx.2.2.2, fun s hs => ?_⟩
  have ex1 := ratioCoefs_num x (walkCount x u) s
  have ex2 := ratioCoefs_den x (walkCount x u) s
  have ey1 := ratioCoefs_num y (walkCount y u) s
  have ey2 := ratioCoefs_den y (walkCount y u) s
  rw [lpmDeg1_walkCount hx hs] at ex1 ex2
  rw [upmDeg1_walkCount hn hx hs] at ex2
  rw [lpmDeg1_walkCount hy hs] at ey1 ey2
  rw [upmDeg1_walkCount hm hy hs] at ey2
  unfold omegaCross
  rw [← ex2, ← ey2, ← ex1, ← ey1]
  ring

theorem LPMratio_one_eq_of_eq_at_knots_mid (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ}
    {y : Fin m → ℝ} (hk : ∀ t, IsKnot x y t → LPMratio 1 t x = LPMratio 1 t y)
    (hmid : ∀ u w, (u, some w) ∈ walkSteps (knotList x y) →
      LPMratio 1 ((u + w) / 2) x = LPMratio 1 ((u + w) / 2) y) (t : ℝ) :
    LPMratio 1 t x = LPMratio 1 t y := by
  set K := Finset.univ.image x ∪ Finset.univ.image y
  have hxK : ∀ i, x i ∈ K := fun i => isKnot_iff_mem.mp (Or.inl ⟨i, rfl⟩)
  have hyK : ∀ j, y j ∈ K := fun j => isKnot_iff_mem.mp (Or.inr ⟨j, rfl⟩)
  by_cases hlo : ∃ z ∈ K, z ≤ t
  swap
  · push Not at hlo
    have h0 : ∀ {k : ℕ} (z : Fin k → ℝ), 0 < k → (∀ i, t < z i) → LPMratio 1 t z = 0 :=
      fun z hk hz => by
        rw [LPMratio_one_of_pos (ratio_denom_pos ⟨⟨0, hk⟩, (hz _).ne'⟩),
          LPM_eq_zero_of_lt hz, zero_div]
    rw [h0 x hn fun i => hlo _ (hxK i), h0 y hm fun j => hlo _ (hyK j)]
  by_cases hhi : ∃ z ∈ K, t ≤ z
  swap
  · push Not at hhi
    have h1 : ∀ {k : ℕ} (z : Fin k → ℝ), 0 < k → (∀ i, z i < t) → LPMratio 1 t z = 1 :=
      fun z hk hz => by
        have hL := LPM_pos_of_lt (d := 1) (hz ⟨0, hk⟩)
        rw [LPMratio_one_of_pos (ratio_denom_pos ⟨⟨0, hk⟩, (hz _).ne⟩),
          UPM_eq_zero_of_le fun i => (hz i).le, add_zero, div_self hL.ne']
    rw [h1 x hn fun i => hhi _ (hxK i), h1 y hm fun j => hhi _ (hyK j)]
  obtain ⟨u, hu, w, hw, hs, hgap⟩ := exists_bracket K hlo hhi
  have hku := isKnot_iff_mem.mpr hu
  have hkw := isKnot_iff_mem.mpr hw
  rcases hs.1.eq_or_lt with rfl | hut
  · exact hk _ hku
  rcases hs.2.eq_or_lt with rfl | htw
  · exact hk _ hkw
  have hx : ∀ i, x i ∉ Ioo u w := fun i => hgap _ (hxK i)
  have hy : ∀ j, y j ∉ Ioo u w := fun j => hgap _ (hyK j)
  have hstep : (u, some w) ∈ walkSteps (knotList x y) :=
    mem_walkSteps_of_adjacent (knotList_pairwise x y) (mem_knotList.mpr hku)
      (mem_knotList.mpr hkw) (hut.trans htw) fun z hz => hgap z (isKnot_iff_mem.mp
        (mem_knotList.mp hz))
  obtain ⟨A, B, C, hq⟩ := omegaCross_quadratic hn hm hx hy
  have huw : u < w := hut.trans htw
  have hzero := quad_eq_zero_of_three (A := A) (B := B) (C := C) (v := (u + w) / 2)
    (by intro h; linarith) (by intro h; linarith) huw.ne
    (by rw [← hq u ⟨le_rfl, huw.le⟩]; exact omegaCross_eq_zero (hk _ hku))
    (by rw [← hq _ ⟨by linarith, by linarith⟩]; exact omegaCross_eq_zero (hmid u w hstep))
    (by rw [← hq w ⟨huw.le, le_rfl⟩]; exact omegaCross_eq_zero (hk _ hkw)) t
  rw [← hq t ⟨hut.le, htw.le⟩, omegaCross] at hzero
  have hDx := ratio_denom_pos_of_gap hn hx ⟨hut, htw⟩
  have hDy := ratio_denom_pos_of_gap hm hy ⟨hut, htw⟩
  rw [LPMratio_one_of_pos hDx, LPMratio_one_of_pos hDy, div_eq_div_iff hDx.ne' hDy.ne']
  linarith

def omegaMv (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) : Bool :=
  match e.2 with
  | none => false
  | some tn => decide (omegaMid x y e.1 tn > 0)

def omegaMd (x : Fin n → ℝ) (y : Fin m → ℝ) (e : ℝ × Option ℝ) : Bool :=
  match e.2 with
  | none => false
  | some tn => decide (|omegaMid x y e.1 tn| > 0)

lemma omegaStep_fst {x : Fin n → ℝ} {y : Fin m → ℝ} {s : Bool × Bool} {e : ℝ × Option ℝ}
    (hs : s.1 = false) (hkv : omegaKv x y e = false) (hsv : omegaSv x y e = false)
    (hmv : omegaMv x y e = false) : (omegaStep x y s e).1 = false := by
  obtain ⟨g, d⟩ := s
  obtain ⟨t, o⟩ := e
  simp only at hs
  subst hs
  have h1 : (flagStep (omegaKv x y (t, o)) (omegaKd x y (t, o)) (omegaSv x y (t, o))
      (false, d)).1 = false := by
    rw [flagStep_fst, hkv, hsv]; rfl
  cases o with
  | none => exact h1
  | some tn =>
    simp only [omegaMv] at hmv
    simp only [omegaStep]
    split_ifs
    · simp only [h1, Bool.false_or]; exact hmv
    · exact h1

lemma flagStep_snd_of (kv kd vv : Bool) (s : Bool × Bool) (h : s.2 = true ∨ kd = true) :
    (flagStep kv kd vv s).2 = true := by
  obtain ⟨g, d⟩ := s
  cases g <;> cases kv <;> cases vv <;> cases d <;> cases kd <;> simp_all [flagStep]

lemma omegaStep_snd {x : Fin n → ℝ} {y : Fin m → ℝ} {s : Bool × Bool} {e : ℝ × Option ℝ}
    (h : s.2 = true ∨ omegaKd x y e = true ∨ omegaMd x y e = true) :
    (omegaStep x y s e).2 = true := by
  obtain ⟨t, o⟩ := e
  cases o with
  | none =>
    simp only [omegaMd, Bool.false_eq_true, or_false] at h
    exact flagStep_snd_of _ _ _ _ h
  | some tn =>
    simp only [omegaStep]
    by_cases hc : (flagStep (omegaKv x y (t, some tn)) (omegaKd x y (t, some tn))
        (omegaSv x y (t, some tn)) s).2 = false
    · rw [ite_eq_left hc]
      rcases h with h | h | h
      · have := flagStep_snd_of (omegaKv x y (t, some tn)) (omegaKd x y (t, some tn))
          (omegaSv x y (t, some tn)) s (Or.inl h)
        rw [hc] at this; exact absurd this Bool.false_ne_true
      · have := flagStep_snd_of (omegaKv x y (t, some tn)) (omegaKd x y (t, some tn))
          (omegaSv x y (t, some tn)) s (Or.inr h)
        rw [hc] at this; exact absurd this Bool.false_ne_true
      · simp only [omegaMd] at h
        simp only [Bool.or_eq_true]
        exact Or.inr h
    · rw [ite_eq_right hc]
      simpa using hc

lemma foldl_omegaStep_fst {x : Fin n → ℝ} {y : Fin m → ℝ} (L : List (ℝ × Option ℝ))
    (hL : ∀ e ∈ L, omegaKv x y e = false ∧ omegaSv x y e = false ∧ omegaMv x y e = false)
    (s : Bool × Bool) (hs : s.1 = false) : (L.foldl (omegaStep x y) s).1 = false := by
  induction L generalizing s with
  | nil => exact hs
  | cons e L ih =>
    have he := hL e (by simp)
    exact ih (fun e' he' => hL e' (by simp [he'])) _ (omegaStep_fst hs he.1 he.2.1 he.2.2)

lemma foldl_omegaStep_snd {x : Fin n → ℝ} {y : Fin m → ℝ} (L : List (ℝ × Option ℝ))
    (s : Bool × Bool) (h : s.2 = true ∨ ∃ e ∈ L, omegaKd x y e = true ∨ omegaMd x y e = true) :
    (L.foldl (omegaStep x y) s).2 = true := by
  induction L generalizing s with
  | nil => simpa using h
  | cons e L ih =>
    rw [List.foldl_cons]
    refine ih _ ?_
    rcases h with h | ⟨e', he', h'⟩
    · exact Or.inl (omegaStep_snd (Or.inl h))
    · rcases List.mem_cons.mp he' with rfl | he'
      · exact Or.inl (omegaStep_snd (Or.inr h'))
      · exact Or.inr ⟨e', he', h'⟩

lemma UPM_eq_sortedVals (d : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    UPM d t x = ((sortedVals x).map (upperDev d t)).sum / (sortedVals x).length := by
  rw [UPM, length_sortedVals, sum_eq_multiset_sum, ← Multiset.sum_coe, ← Multiset.map_coe,
    coe_sortedVals, Multiset.map_map]
  rfl

theorem LPMratio_eq_of_identical {x : Fin n → ℝ} {y : Fin m → ℝ} (h : identicalSamples x y)
    (d : ℕ) (t : ℝ) : LPMratio d t x = LPMratio d t y := by
  have hlen : n = m := by
    have := congrArg List.length h
    rwa [length_sortedVals, length_sortedVals] at this
  have hU : UPM d t x = UPM d t y := by
    unfold identicalSamples at h
    rw [UPM_eq_sortedVals, UPM_eq_sortedVals, h]
  have hall : (∀ i, x i = t) ↔ ∀ j, y j = t := by
    have hx : (∀ i, x i = t) ↔ ∀ v ∈ sortedVals x, v = t := by
      simp only [mem_sortedVals]; constructor
      · rintro h v ⟨i, rfl⟩; exact h i
      · exact fun h i => h _ ⟨i, rfl⟩
    have hy : (∀ j, y j = t) ↔ ∀ v ∈ sortedVals y, v = t := by
      simp only [mem_sortedVals]; constructor
      · rintro h v ⟨i, rfl⟩; exact h i
      · exact fun h i => h _ ⟨i, rfl⟩
    unfold identicalSamples at h
    rw [hx, hy, h]
  subst hlen
  simp only [LPMratio, LPM_eq_of_identical h, hU, hall]

theorem sdDomPair_omega_of_dominates (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ)
    (y : Fin m → ℝ) (hle : ∀ t, LPMratio 1 t x ≤ LPMratio 1 t y)
    (hne : ∃ t, LPMratio 1 t x ≠ LPMratio 1 t y) : sdDomPair x y 1 false = 1 := by
  have hmin : sampleMin y ≤ sampleMin x := by
    refine (sampleMin_le_iff hn hm x y).mpr fun i => ?_
    by_contra hc
    push Not at hc
    obtain ⟨j0, hj0⟩ := exists_eq_sampleMin hm y
    set t := (x i + sampleMin y) / 2
    have hxt : x i < t := by have := hc j0; rw [hj0] at this; simp only [t]; linarith
    have hyt : ∀ j, t < y j := fun j => by
      have h1 := sampleMin_le y j
      have h2 := hc j0
      rw [hj0] at h2
      simp only [t]; linarith
    have hRy : LPMratio 1 t y = 0 := by
      rw [LPMratio_one_of_pos (ratio_denom_pos ⟨⟨0, hm⟩, (hyt _).ne'⟩), LPM_eq_zero_of_lt hyt,
        zero_div]
    have hRx : 0 < LPMratio 1 t x := by
      rw [LPMratio_one_of_pos (ratio_denom_pos ⟨i, hxt.ne⟩)]
      exact div_pos (LPM_pos_of_lt hxt) (ratio_denom_pos ⟨i, hxt.ne⟩)
    linarith [hle t]
  have hid : ¬ identicalSamples x y := fun h => by
    obtain ⟨t, ht⟩ := hne
    exact ht (LPMratio_eq_of_identical h 1 t)
  unfold sdDomPair
  simp only [↓reduceIte, ge_iff_le, hmin, not_true_eq_false, hid, Bool.false_eq_true]
  set L := walkSteps (knotList x y)
  have hfst : (L.foldl (omegaStep x y) (false, false)).1 = false := by
    refine foldl_omegaStep_fst L (fun e he => ?_) _ rfl
    obtain ⟨t, o⟩ := e
    refine ⟨?_, ?_, ?_⟩
    · simp only [omegaKv, omegaRatio_eq hn, omegaRatio_eq hm, decide_eq_false_iff_not]
      linarith [hle t]
    · cases o with
      | none => rfl
      | some tn =>
        obtain ⟨-, -, -, gx, gy⟩ := gap_of_mem_walkSteps he
        simp only [omegaSv, decide_eq_false_iff_not]
        intro hfire
        simp only [stationaryFires] at hfire
        obtain ⟨-, -, ⟨h1, h2⟩, ⟨hDx, hDy⟩, hpos⟩ := hfire
        rw [omega_formula_eq hn gx ⟨h1.le, h2.le⟩ hDx,
          omega_formula_eq hm gy ⟨h1.le, h2.le⟩ hDy] at hpos
        exact absurd hpos (not_lt.mpr (sub_nonpos.mpr (hle _)))
    · cases o with
      | none => rfl
      | some tn =>
        obtain ⟨-, -, htn, gx, gy⟩ := gap_of_mem_walkSteps he
        simp only [omegaMv, omegaMid_eq hn hm gx gy htn, decide_eq_false_iff_not]
        linarith [hle ((t + tn) / 2)]
  have hsnd : (L.foldl (omegaStep x y) (false, false)).2 = true := by
    refine foldl_omegaStep_snd L _ (Or.inr ?_)
    by_contra hc
    push Not at hc
    obtain ⟨t, ht⟩ := hne
    refine ht (LPMratio_one_eq_of_eq_at_knots_mid hn hm (fun u hu => ?_) (fun u w hw => ?_) t)
    · obtain ⟨o, ho⟩ := exists_mem_walkSteps (mem_knotList.mpr hu)
      have := (hc _ ho).1
      simp only [omegaKd, omegaRatio_eq hn, omegaRatio_eq hm, ne_eq, decide_eq_true_eq,
        abs_pos, sub_ne_zero, not_not] at this
      exact this
    · obtain ⟨-, -, huw, gx, gy⟩ := gap_of_mem_walkSteps hw
      have := (hc _ hw).2
      simp only [omegaMd, omegaMid_eq hn hm gx gy huw, ne_eq, decide_eq_true_eq,
        abs_pos, sub_ne_zero, not_not] at this
      exact this
  unfold verdict
  simp [hfst, hsnd]

theorem omega_core {ax bx gx dx ay by' gy dy u w s0 : ℝ}
    (hDx : ∀ s ∈ Icc u w, 0 < gx * s + dx) (hDy : ∀ s ∈ Icc u w, 0 < gy * s + dy)
    (hu : (ax * u + bx) / (gx * u + dx) - (ay * u + by') / (gy * u + dy) ≤ 0)
    (hw : (ax * w + bx) / (gx * w + dx) - (ay * w + by') / (gy * w + dy) ≤ 0)
    (hs0 : s0 ∈ Ioo u w)
    (h0 : 0 < (ax * s0 + bx) / (gx * s0 + dx) - (ay * s0 + by') / (gy * s0 + dy)) :
    stationaryFires ax bx gx dx ay by' gy dy u w := by
  simp only [stationaryFires]
  set Ax := ax * dx - bx * gx with hAx
  set Ay := ay * dy - by' * gy with hAy
  set r := √(Ay / Ax) with hr
  set den := gy - r * gx with hden
  set ts := (r * dx - dy) / den with hts
  set q : ℝ → ℝ := fun s => (ax * s + bx) * (gy * s + dy) - (ay * s + by') * (gx * s + dx)
    with hq
  have hfq : ∀ s ∈ Icc u w, (ax * s + bx) / (gx * s + dx) - (ay * s + by') / (gy * s + dy)
      = q s / ((gx * s + dx) * (gy * s + dy)) := fun s hs => by
    rw [div_sub_div _ _ (hDx s hs).ne' (hDy s hs).ne']
    simp only [hq]; ring
  have huw : u ≤ w := hs0.1.le.trans hs0.2.le
  have hu' : q u ≤ 0 := by
    have := hfq u ⟨le_rfl, huw⟩
    rw [this, div_nonpos_iff] at hu
    have := mul_pos (hDx u ⟨le_rfl, huw⟩) (hDy u ⟨le_rfl, huw⟩)
    rcases hu with h | h <;> [linarith [h.2]; exact h.1]
  have hw' : q w ≤ 0 := by
    have := hfq w ⟨huw, le_rfl⟩
    rw [this, div_nonpos_iff] at hw
    have := mul_pos (hDx w ⟨huw, le_rfl⟩) (hDy w ⟨huw, le_rfl⟩)
    rcases hw with h | h <;> [linarith [h.2]; exact h.1]
  have hs0I : s0 ∈ Icc u w := Ioo_subset_Icc_self hs0
  have h0' : 0 < q s0 := by
    rw [hfq s0 hs0I] at h0
    exact (div_pos_iff_of_pos_right (mul_pos (hDx s0 hs0I) (hDy s0 hs0I))).mp h0
  have hqc : Continuous q := by simp only [hq]; fun_prop
  obtain ⟨u1, hu1I, hu1⟩ := intermediate_value_Icc hs0.1.le hqc.continuousOn ⟨hu', h0'.le⟩
  obtain ⟨u2, hu2I, hu2⟩ := intermediate_value_Icc' hs0.2.le hqc.continuousOn ⟨hw', h0'.le⟩
  have hu1s : u1 < s0 := lt_of_le_of_ne hu1I.2 (by rintro rfl; linarith)
  have hs2 : s0 < u2 := lt_of_le_of_ne hu2I.1 (by rintro rfl; linarith)
  have hu1w : u1 ∈ Icc u w := ⟨hu1I.1, hu1I.2.trans hs0.2.le⟩
  have hu2w : u2 ∈ Icc u w := ⟨hs0.1.le.trans hu2I.1, hu2I.2⟩
  set A := ax * gy - ay * gx with hA
  have hfac : ∀ s, q s = A * (s - u1) * (s - u2) := by
    have e1 : (u1 - u2) * (A * (u1 + u2) + (ax * dy + bx * gy - ay * dx - by' * gx)) = 0 := by
      simp only [hq] at hu1 hu2; linear_combination hu1 - hu2
    have e2 : A * (u1 + u2) + (ax * dy + bx * gy - ay * dx - by' * gx) = 0 :=
      (mul_eq_zero.mp e1).resolve_left (by linarith)
    intro s
    simp only [hq] at hu1 ⊢
    linear_combination hu1 + (s - u1) * e2
  have hAneg : A < 0 := by
    by_contra hA
    push Not at hA
    have h1 := mul_nonneg hA (mul_pos (sub_pos.mpr hu1s) (sub_pos.mpr hs2)).le
    have h2 : q s0 = -(A * ((s0 - u1) * (u2 - s0))) := by rw [hfac s0]; ring
    linarith
  set X1 := gx * u1 + dx
  set X2 := gx * u2 + dx
  set Y1 := gy * u1 + dy
  set Y2 := gy * u2 + dy
  have hX1 : 0 < X1 := hDx u1 hu1w
  have hX2 : 0 < X2 := hDx u2 hu2w
  have hY1 : 0 < Y1 := hDy u1 hu1w
  have hY2 : 0 < Y2 := hDy u2 hu2w
  have hq1 : (ax * u1 + bx) * Y1 = (ay * u1 + by') * X1 := by
    have := hu1; simp only [hq] at this; linarith
  have hq2 : (ax * u2 + bx) * Y2 = (ay * u2 + by') * X2 := by
    have := hu2; simp only [hq] at this; linarith
  have hne12 : u2 - u1 ≠ 0 := by linarith
  have hkey : Ax * (Y1 * Y2) = Ay * (X1 * X2) := by
    have : (u2 - u1) * (Ax * (Y1 * Y2) - Ay * (X1 * X2)) = 0 := by
      simp only [hAx, hAy, X1, X2, Y1, Y2]
      linear_combination (gx * u1 + dx) * (gy * u1 + dy) * hq2 -
        (gx * u2 + dx) * (gy * u2 + dy) * hq1
    linarith [(mul_eq_zero.mp this).resolve_left hne12]
  have hP := mul_pos hX1 hX2
  have hQ := mul_pos hY1 hY2
  have hs0D := mul_pos (hDx s0 hs0I) (hDy s0 hs0I)
  have hAx0 : Ax ≠ 0 := by
    intro hAx0
    have hAy0 : Ay = 0 := by
      rw [hAx0, zero_mul] at hkey
      exact (mul_eq_zero.mp hkey.symm).resolve_right hP.ne'
    have ex : (ax * s0 + bx) * X1 = (ax * u1 + bx) * (gx * s0 + dx) := by
      have : Ax * (s0 - u1) = 0 := by rw [hAx0, zero_mul]
      simp only [hAx, X1] at this ⊢; linear_combination this
    have ey : (ay * s0 + by') * Y1 = (ay * u1 + by') * (gy * s0 + dy) := by
      have : Ay * (s0 - u1) = 0 := by rw [hAy0, zero_mul]
      simp only [hAy, Y1] at this ⊢; linear_combination this
    have : q s0 * (X1 * Y1) = 0 := by
      simp only [hq]
      linear_combination (gy * s0 + dy) * Y1 * ex - (gx * s0 + dx) * X1 * ey +
        (gx * s0 + dx) * (gy * s0 + dy) * hq1
    have := mul_pos h0' (mul_pos hX1 hY1)
    linarith
  have hAy0 : Ay ≠ 0 := by
    intro h
    rw [h, zero_mul] at hkey
    exact hAx0 ((mul_eq_zero.mp hkey).resolve_right hQ.ne')
  have hsign : Ax > 0 ↔ Ay > 0 := by
    constructor
    · intro h; by_contra h'; push Not at h'
      linarith [mul_pos h hQ, mul_nonpos_of_nonpos_of_nonneg h' hP.le]
    · intro h; by_contra h'; push Not at h'
      linarith [mul_pos h hP, mul_nonpos_of_nonpos_of_nonneg h' hQ.le]
  have hratio : Ay / Ax = (Y1 * Y2) / (X1 * X2) := by
    rw [div_eq_div_iff hAx0 hP.ne']; linarith
  have hrpos : 0 < r := Real.sqrt_pos.mpr (by rw [hratio]; positivity)
  have hr2 : r ^ 2 * (X1 * X2) = Y1 * Y2 := by
    rw [hr, Real.sq_sqrt (by rw [hratio]; positivity), hratio]; field_simp
  set c0 := dy - r * dx with hc0
  have hphi : ∀ s, gy * s + dy - r * (gx * s + dx) = den * s + c0 := fun s => by
    simp only [hden, hc0]; ring
  by_cases hab : X1 * Y2 = Y1 * X2
  · exfalso
    have h1 : r * X1 = Y1 := by
      have e : (r * X1) ^ 2 * X2 = Y1 ^ 2 * X2 := by
        linear_combination X1 * hr2 + Y1 * hab
      have e' : (r * X1) ^ 2 = Y1 ^ 2 := (mul_left_inj' hX2.ne').mp e
      exact (pow_left_inj₀ (mul_pos hrpos hX1).le hY1.le two_ne_zero).mp e'
    have h2 : r * X2 = Y2 := by
      have e : (r * X2) ^ 2 * X1 = Y2 ^ 2 * X1 := by
        linear_combination X2 * hr2 - Y2 * hab
      have e' : (r * X2) ^ 2 = Y2 ^ 2 := (mul_left_inj' hX1.ne').mp e
      exact (pow_left_inj₀ (mul_pos hrpos hX2).le hY2.le two_ne_zero).mp e'
    have hd0 : den = 0 ∧ c0 = 0 := by
      have a1 := hphi u1; have a2 := hphi u2
      have z1 : den * u1 + c0 = 0 := by rw [← a1]; simp only [X1, Y1] at h1; linarith
      have z2 : den * u2 + c0 = 0 := by rw [← a2]; simp only [X2, Y2] at h2; linarith
      have : den * (u2 - u1) = 0 := by linear_combination z2 - z1
      have hd : den = 0 := (mul_eq_zero.mp this).resolve_right hne12
      exact ⟨hd, by rw [hd] at z1; linarith⟩
    have hDy_eq : ∀ s, gy * s + dy = r * (gx * s + dx) := fun s => by
      have := hphi s; rw [hd0.1, hd0.2] at this; linarith
    -- `psi = r Nx - Ny` vanishes at `u1`, `u2`, so it is identically zero
    have p1 : r * (ax * u1 + bx) - (ay * u1 + by') = 0 := by
      have := hq1; rw [show Y1 = r * X1 from h1.symm] at this
      have h' : X1 * (r * (ax * u1 + bx) - (ay * u1 + by')) = 0 := by linear_combination this
      exact (mul_eq_zero.mp h').resolve_left hX1.ne'
    have p2 : r * (ax * u2 + bx) - (ay * u2 + by') = 0 := by
      have := hq2; rw [show Y2 = r * X2 from h2.symm] at this
      have h' : X2 * (r * (ax * u2 + bx) - (ay * u2 + by')) = 0 := by linear_combination this
      exact (mul_eq_zero.mp h').resolve_left hX2.ne'
    have pslope : r * ax - ay = 0 := by
      have : (u2 - u1) * (r * ax - ay) = 0 := by linear_combination p2 - p1
      exact (mul_eq_zero.mp this).resolve_left hne12
    have p0 : r * (ax * s0 + bx) - (ay * s0 + by') = 0 := by
      linear_combination p1 + (s0 - u1) * pslope
    have : q s0 = 0 := by
      simp only [hq]; rw [hDy_eq s0]; linear_combination (gx * s0 + dx) * p0
    linarith
  have hprod : (den * u1 + c0) * (den * u2 + c0) < 0 := by
    rw [← hphi u1, ← hphi u2]
    change (Y1 - r * X1) * (Y2 - r * X2) < 0
    have hlt : 2 * (Y1 * Y2) < r * (X1 * Y2 + Y1 * X2) := by
      by_contra hc
      push Not at hc
      have hnn : 0 ≤ r * (X1 * Y2 + Y1 * X2) := by positivity
      have hsq : (r * (X1 * Y2 + Y1 * X2)) ^ 2 ≤ (2 * (Y1 * Y2)) ^ 2 :=
        pow_le_pow_left₀ hnn hc 2
      have e : (r * (X1 * Y2 + Y1 * X2)) ^ 2 * (X1 * X2) =
          (Y1 * Y2) * (X1 * Y2 + Y1 * X2) ^ 2 := by
        linear_combination (X1 * Y2 + Y1 * X2) ^ 2 * hr2
      have h3 : (Y1 * Y2) * (X1 * Y2 + Y1 * X2) ^ 2 ≤
          (Y1 * Y2) * (4 * (X1 * Y2) * (Y1 * X2)) :=
        calc (Y1 * Y2) * (X1 * Y2 + Y1 * X2) ^ 2
            = (r * (X1 * Y2 + Y1 * X2)) ^ 2 * (X1 * X2) := e.symm
          _ ≤ (2 * (Y1 * Y2)) ^ 2 * (X1 * X2) := mul_le_mul_of_nonneg_right hsq hP.le
          _ = (Y1 * Y2) * (4 * (X1 * Y2) * (Y1 * X2)) := by ring
      have h4 : (X1 * Y2 + Y1 * X2) ^ 2 ≤ 4 * (X1 * Y2) * (Y1 * X2) :=
        le_of_mul_le_mul_left h3 hQ
      have h5 : (X1 * Y2 - Y1 * X2) ^ 2 = 0 :=
        le_antisymm (by linarith [show (X1 * Y2 - Y1 * X2) ^ 2 = (X1 * Y2 + Y1 * X2) ^ 2 -
          4 * (X1 * Y2) * (Y1 * X2) by ring]) (sq_nonneg _)
      exact hab (sub_eq_zero.mp (pow_eq_zero_iff two_ne_zero |>.mp h5))
    have e : (Y1 - r * X1) * (Y2 - r * X2) = 2 * (Y1 * Y2) - r * (X1 * Y2 + Y1 * X2) := by
      linear_combination hr2
    rw [e]
    linarith
  have hden0 : den ≠ 0 := by
    intro h
    rw [h, zero_mul, zero_mul, zero_add] at hprod
    linarith [mul_self_nonneg c0]
  have hroot : den * ts + c0 = 0 := by
    simp only [hts, hc0]; field_simp; ring
  have hts_in : u1 < ts ∧ ts < u2 := by
    have e1 : den * u1 + c0 = den * (u1 - ts) := by linear_combination hroot
    have e2 : den * u2 + c0 = den * (u2 - ts) := by linear_combination hroot
    rw [e1, e2, show den * (u1 - ts) * (den * (u2 - ts)) = den ^ 2 * ((u1 - ts) * (u2 - ts))
      by ring] at hprod
    have hneg : (u1 - ts) * (u2 - ts) < 0 := by
      by_contra h
      push Not at h
      linarith [mul_nonneg (sq_nonneg den) h]
    constructor
    · by_contra h
      push Not at h
      linarith [mul_nonneg (sub_nonneg.mpr h) (sub_nonneg.mpr (h.trans (hu1s.trans hs2).le))]
    · by_contra h
      push Not at h
      linarith [mul_nonneg_of_nonpos_of_nonpos (sub_nonpos.mpr ((hu1s.trans hs2).le.trans h))
        (sub_nonpos.mpr h)]
  have htsI : ts ∈ Icc u w := ⟨hu1I.1.trans hts_in.1.le, hts_in.2.le.trans hu2I.2⟩
  refine ⟨⟨hAx0, hAy0, hsign⟩, hden0, ⟨hu1I.1.trans_lt hts_in.1, hts_in.2.trans_le hu2I.2⟩,
    ⟨hDx ts htsI, hDy ts htsI⟩, ?_⟩
  rw [hfq ts htsI]
  have hqts : 0 < q ts := by
    rw [hfac ts, mul_assoc]
    exact mul_pos_of_neg_of_neg hAneg
      (mul_neg_of_pos_of_neg (sub_pos.mpr hts_in.1) (sub_neg.mpr hts_in.2))
  exact div_pos hqts (mul_pos (hDx ts htsI) (hDy ts htsI))

lemma affine_pos_on_Icc {g d u w : ℝ} (hu : 0 < g * u + d) (hw : 0 < g * w + d) :
    ∀ s ∈ Icc u w, 0 < g * s + d := by
  intro s hs
  rcases hs.2.eq_or_lt with rfl | hsw
  · exact hw
  have huw : u < w := hs.1.trans_lt hsw
  have key : (w - u) * (g * s + d) = (w - s) * (g * u + d) + (s - u) * (g * w + d) := by ring
  have : 0 < (w - u) * (g * s + d) := by
    rw [key]
    exact add_pos_of_pos_of_nonneg (mul_pos (sub_pos.mpr hsw) hu)
      (mul_nonneg (sub_nonneg.mpr hs.1) hw.le)
  exact pos_of_mul_pos_right this (sub_pos.mpr huw).le

lemma LPMratio_one_of_forall_lt {k : ℕ} (hk : 0 < k) {t : ℝ} {z : Fin k → ℝ}
    (hz : ∀ i, t < z i) : LPMratio 1 t z = 0 := by
  rw [LPMratio_one_of_pos (ratio_denom_pos ⟨⟨0, hk⟩, (hz _).ne'⟩), LPM_eq_zero_of_lt hz,
    zero_div]

lemma LPMratio_one_of_forall_gt {k : ℕ} (hk : 0 < k) {t : ℝ} {z : Fin k → ℝ}
    (hz : ∀ i, z i < t) : LPMratio 1 t z = 1 := by
  have hL := LPM_pos_of_lt (d := 1) (hz ⟨0, hk⟩)
  rw [LPMratio_one_of_pos (ratio_denom_pos ⟨⟨0, hk⟩, (hz _).ne⟩),
    UPM_eq_zero_of_le fun i => (hz i).le, add_zero, div_self hL.ne']

theorem omega_interval_le (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ} {y : Fin m → ℝ}
    (hmin : ∀ i, ∃ j, y j ≤ x i) {u w : ℝ} (hstep : (u, some w) ∈ walkSteps (knotList x y))
    (hu : LPMratio 1 u x ≤ LPMratio 1 u y) (hw : LPMratio 1 w x ≤ LPMratio 1 w y)
    (hsv : omegaSv x y (u, some w) = false) (s : ℝ) (hs : s ∈ Ioo u w) :
    LPMratio 1 s x ≤ LPMratio 1 s y := by
  obtain ⟨-, -, huw, gx, gy⟩ := gap_of_mem_walkSteps hstep
  by_contra h0
  push Not at h0
  have hsI := Ioo_subset_Icc_self hs
  have hmono_y := LPMratio_monotone 1 y hs.1.le
  simp only at hmono_y
  have hDxu : 0 < LPM 1 u x + UPM 1 u x := by
    rcases (ratio_denom_nonneg u x).eq_or_lt with h | h
    · have hall := (ratio_denom_eq_zero_iff one_ne_zero u x).mp h.symm
      have := LPMratio_of_forall_eq one_ne_zero hn hall
      linarith [LPMratio_le_one 1 s x]
    · exact h
  have hDyu : 0 < LPM 1 u y + UPM 1 u y := by
    rcases (ratio_denom_nonneg u y).eq_or_lt with h | h
    · have hall := (ratio_denom_eq_zero_iff one_ne_zero u y).mp h.symm
      have := LPMratio_of_forall_eq one_ne_zero hm hall
      linarith [LPMratio_le_one 1 s x]
    · exact h
  have hxs : (∀ i, x i = w) ∨ (∀ j, y j = w) → ∀ i, s < x i := by
    rintro (h | h) i
    · rw [h i]; exact hs.2
    · obtain ⟨j, hj⟩ := hmin i
      rw [h j] at hj
      exact hs.2.trans_le hj
  have hDxw : 0 < LPM 1 w x + UPM 1 w x := by
    rcases (ratio_denom_nonneg w x).eq_or_lt with h | h
    · have hall := (ratio_denom_eq_zero_iff one_ne_zero w x).mp h.symm
      rw [LPMratio_one_of_forall_lt hn (hxs (Or.inl hall))] at h0
      linarith [LPMratio_nonneg 1 s y]
    · exact h
  have hDyw : 0 < LPM 1 w y + UPM 1 w y := by
    rcases (ratio_denom_nonneg w y).eq_or_lt with h | h
    · have hall := (ratio_denom_eq_zero_iff one_ne_zero w y).mp h.symm
      rw [LPMratio_one_of_forall_lt hn (hxs (Or.inr hall))] at h0
      linarith [LPMratio_nonneg 1 s y]
    · exact h
  have huI : u ∈ Icc u w := ⟨le_rfl, huw.le⟩
  have hwI : w ∈ Icc u w := ⟨huw.le, le_rfl⟩
  have hDenx : ∀ r ∈ Icc u w, (ratioCoefs x (walkCount x u)).2.2.1 * r +
      (ratioCoefs x (walkCount x u)).2.2.2 = LPM 1 r x + UPM 1 r x := fun r hr => by
    rw [ratioCoefs_den, lpmDeg1_walkCount gx hr, upmDeg1_walkCount hn gx hr]
  have hDeny : ∀ r ∈ Icc u w, (ratioCoefs y (walkCount y u)).2.2.1 * r +
      (ratioCoefs y (walkCount y u)).2.2.2 = LPM 1 r y + UPM 1 r y := fun r hr => by
    rw [ratioCoefs_den, lpmDeg1_walkCount gy hr, upmDeg1_walkCount hm gy hr]
  have hDx := affine_pos_on_Icc (by rw [hDenx u huI]; exact hDxu) (by rw [hDenx w hwI]; exact hDxw)
  have hDy := affine_pos_on_Icc (by rw [hDeny u huI]; exact hDyu) (by rw [hDeny w hwI]; exact hDyw)
  have hRx : ∀ r ∈ Icc u w, _ = LPMratio 1 r x := fun r hr => omega_formula_eq hn gx hr (hDx r hr)
  have hRy : ∀ r ∈ Icc u w, _ = LPMratio 1 r y := fun r hr => omega_formula_eq hm gy hr (hDy r hr)
  have hfire := omega_core hDx hDy (by rw [hRx u huI, hRy u huI]; linarith)
    (by rw [hRx w hwI, hRy w hwI]; linarith) hs (by rw [hRx s hsI, hRy s hsI]; linarith)
  simp only [omegaSv, decide_eq_false_iff_not] at hsv
  exact hsv hfire

theorem LPMratio_one_le_of_omega_walk (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ}
    {y : Fin m → ℝ} (hmin : ∀ i, ∃ j, y j ≤ x i)
    (hk : ∀ t, IsKnot x y t → LPMratio 1 t x ≤ LPMratio 1 t y)
    (hsv : ∀ u w, (u, some w) ∈ walkSteps (knotList x y) → omegaSv x y (u, some w) = false)
    (t : ℝ) : LPMratio 1 t x ≤ LPMratio 1 t y := by
  set K := Finset.univ.image x ∪ Finset.univ.image y
  have hxK : ∀ i, x i ∈ K := fun i => isKnot_iff_mem.mp (Or.inl ⟨i, rfl⟩)
  have hyK : ∀ j, y j ∈ K := fun j => isKnot_iff_mem.mp (Or.inr ⟨j, rfl⟩)
  by_cases hlo : ∃ z ∈ K, z ≤ t
  swap
  · push Not at hlo
    rw [LPMratio_one_of_forall_lt hn fun i => hlo _ (hxK i),
      LPMratio_one_of_forall_lt hm fun j => hlo _ (hyK j)]
  by_cases hhi : ∃ z ∈ K, t ≤ z
  swap
  · push Not at hhi
    rw [LPMratio_one_of_forall_gt hn fun i => hhi _ (hxK i),
      LPMratio_one_of_forall_gt hm fun j => hhi _ (hyK j)]
  obtain ⟨u, hu, w, hw, hs, hgap⟩ := exists_bracket K hlo hhi
  have hku := isKnot_iff_mem.mpr hu
  have hkw := isKnot_iff_mem.mpr hw
  rcases hs.1.eq_or_lt with rfl | hut
  · exact hk _ hku
  rcases hs.2.eq_or_lt with rfl | htw
  · exact hk _ hkw
  have hstep : (u, some w) ∈ walkSteps (knotList x y) :=
    mem_walkSteps_of_adjacent (knotList_pairwise x y) (mem_knotList.mpr hku)
      (mem_knotList.mpr hkw) (hut.trans htw) fun z hz => hgap z (isKnot_iff_mem.mp
        (mem_knotList.mp hz))
  exact omega_interval_le hn hm hmin hstep (hk _ hku) (hk _ hkw) (hsv u w hstep) t ⟨hut, htw⟩

lemma omegaStep_flag_fst {x : Fin n → ℝ} {y : Fin m → ℝ} {s : Bool × Bool} {e : ℝ × Option ℝ}
    (h : (omegaStep x y s e).1 = false) :
    (flagStep (omegaKv x y e) (omegaKd x y e) (omegaSv x y e) s).1 = false := by
  obtain ⟨t, o⟩ := e
  cases o with
  | none => exact h
  | some tn =>
    simp only [omegaStep] at h
    split_ifs at h
    · simp only [Bool.or_eq_false_iff] at h; exact h.1
    · exact h

lemma omegaStep_fst_eq_false {x : Fin n → ℝ} {y : Fin m → ℝ} {s : Bool × Bool}
    {e : ℝ × Option ℝ} (h : (omegaStep x y s e).1 = false) :
    s.1 = false ∧ omegaKv x y e = false ∧ omegaSv x y e = false := by
  have h1 := omegaStep_flag_fst h
  rw [flagStep_fst] at h1
  simpa [Bool.or_eq_false_iff, and_assoc] using h1

lemma omegaStep_snd_eq_true {x : Fin n → ℝ} {y : Fin m → ℝ} {s : Bool × Bool}
    {e : ℝ × Option ℝ} (h1 : (omegaStep x y s e).1 = false) (h2 : (omegaStep x y s e).2 = true) :
    s.2 = true ∨ omegaKd x y e = true ∨ omegaMd x y e = true := by
  have hf := omegaStep_flag_fst h1
  have hs2 := flagStep_snd (omegaKv x y e) (omegaKd x y e) (omegaSv x y e) s hf
  obtain ⟨t, o⟩ := e
  cases o with
  | none =>
    change (flagStep _ _ _ s).2 = true at h2
    rw [hs2, Bool.or_eq_true] at h2
    tauto
  | some tn =>
    simp only [omegaStep] at h2
    split_ifs at h2
    · rw [hs2, Bool.or_eq_true, Bool.or_eq_true] at h2
      simp only [omegaMd]
      tauto
    · rw [hs2, Bool.or_eq_true] at h2
      tauto

lemma foldl_omegaStep_fst_false {x : Fin n → ℝ} {y : Fin m → ℝ} (L : List (ℝ × Option ℝ))
    (s : Bool × Bool) (h : (L.foldl (omegaStep x y) s).1 = false) :
    s.1 = false ∧ ∀ e ∈ L, omegaKv x y e = false ∧ omegaSv x y e = false := by
  induction L generalizing s with
  | nil => exact ⟨h, by simp⟩
  | cons e L ih =>
    rw [List.foldl_cons] at h
    obtain ⟨h1, h2⟩ := ih _ h
    obtain ⟨g, k, v⟩ := omegaStep_fst_eq_false h1
    refine ⟨g, fun e' he' => ?_⟩
    rcases List.mem_cons.mp he' with rfl | he'
    · exact ⟨k, v⟩
    · exact h2 e' he'

lemma foldl_omegaStep_snd_true {x : Fin n → ℝ} {y : Fin m → ℝ} (L : List (ℝ × Option ℝ))
    (s : Bool × Bool) (h1 : (L.foldl (omegaStep x y) s).1 = false)
    (h2 : (L.foldl (omegaStep x y) s).2 = true) :
    s.2 = true ∨ ∃ e ∈ L, omegaKd x y e = true ∨ omegaMd x y e = true := by
  induction L generalizing s with
  | nil => exact Or.inl h2
  | cons e L ih =>
    rw [List.foldl_cons] at h1 h2
    rcases ih _ h1 h2 with h | ⟨e', he', h'⟩
    · have hfst := (foldl_omegaStep_fst_false L _ h1).1
      rcases omegaStep_snd_eq_true hfst h with h | h
      · exact Or.inl h
      · exact Or.inr ⟨e, by simp, h⟩
    · exact Or.inr ⟨e', List.mem_cons_of_mem _ he', h'⟩

theorem sdDomPair_omega_iff (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    sdDomPair x y 1 false = 1 ↔
      (∀ t, LPMratio 1 t x ≤ LPMratio 1 t y) ∧ ∃ t, LPMratio 1 t x ≠ LPMratio 1 t y := by
  refine ⟨fun h => ?_, fun h => sdDomPair_omega_of_dominates hn hm x y h.1 h.2⟩
  set L := walkSteps (knotList x y)
  have hbr : sampleMin y ≤ sampleMin x ∧
      verdict (L.foldl (omegaStep x y) (false, false)) = 1 := by
    unfold sdDomPair at h
    simp only [↓reduceIte, Bool.false_eq_true] at h
    split_ifs at h with h1 h2
    · exact ⟨h1, h⟩
  have hv : (L.foldl (omegaStep x y) (false, false)).1 = false ∧
      (L.foldl (omegaStep x y) (false, false)).2 = true := by
    have := hbr.2
    unfold verdict at this
    split_ifs at this with h1 h2
    exact ⟨by simpa using h2, by simpa using h1⟩
  obtain ⟨-, hall⟩ := foldl_omegaStep_fst_false L _ hv.1
  have hle := LPMratio_one_le_of_omega_walk hn hm ((sampleMin_le_iff hn hm x y).mp hbr.1)
    (fun t ht => by
      obtain ⟨o, ho⟩ := exists_mem_walkSteps (mem_knotList.mpr ht)
      have := (hall _ ho).1
      simp only [omegaKv, omegaRatio_eq hn, omegaRatio_eq hm, decide_eq_false_iff_not,
        not_lt, sub_nonpos] at this
      exact this)
    (fun u w hw => (hall _ hw).2)
  refine ⟨hle, ?_⟩
  rcases foldl_omegaStep_snd_true L _ hv.1 hv.2 with h | ⟨⟨t, o⟩, he, hkd | hmd⟩
  · simp at h
  · refine ⟨t, ?_⟩
    simpa [omegaKd, omegaRatio_eq hn, omegaRatio_eq hm, sub_eq_zero] using hkd
  · cases o with
    | none => simp [omegaMd] at hmd
    | some tn =>
      obtain ⟨-, -, htn, gx, gy⟩ := gap_of_mem_walkSteps he
      refine ⟨(t + tn) / 2, ?_⟩
      simpa [omegaMd, omegaMid_eq hn hm gx gy htn, sub_eq_zero] using hmd

lemma UPM_shiftSample (d : ℕ) (t c : ℝ) (x : Fin n → ℝ) :
    UPM d t (shiftSample x c) = UPM d (t + c) x := by
  unfold UPM upperDev shiftSample
  congr 1
  refine sum_congr rfl fun i _ => ?_
  by_cases h : t + c < x i
  · rw [ite_eq_left (by linarith), ite_eq_left h]; ring_nf
  · rw [ite_eq_right (by intro h'; exact h (by linarith)), ite_eq_right h]

lemma LPMratio_shiftSample (d : ℕ) (t c : ℝ) (x : Fin n → ℝ) :
    LPMratio d t (shiftSample x c) = LPMratio d (t + c) x := by
  have hall : (∀ i, shiftSample x c i = t) ↔ ∀ i, x i = t + c := by
    simp only [shiftSample, sub_eq_iff_eq_add]
  simp only [LPMratio, LPM_shiftSample, UPM_shiftSample, hall]

theorem sdDomPairCentred_omega_eq (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ)
    (y : Fin m → ℝ) : sdDomPairCentred x y 1 false = sdDomPair x y 1 false := by
  have hzero : sdDomPair x y 1 false ≠ 1 → sdDomPair x y 1 false = 0 := fun h =>
    (sdDomPair_eq_zero_or_one x y 1 false).resolve_right h
  have hshift : ∀ c : ℝ, sdDomPair (shiftSample x c) (shiftSample y c) 1 false = 1 ↔
      sdDomPair x y 1 false = 1 := fun c => by
    rw [sdDomPair_omega_iff hn hm, sdDomPair_omega_iff hn hm]
    simp only [LPMratio_shiftSample]
    constructor
    · rintro ⟨h1, t, h2⟩
      exact ⟨fun s => by simpa using h1 (s - c), t + c, h2⟩
    · rintro ⟨h1, t, h2⟩
      exact ⟨fun s => h1 (s + c), t - c, by simpa using h2⟩
  unfold sdDomPairCentred
  simp only [Bool.false_eq_true, and_false, ↓reduceIte, ge_iff_le, Nat.not_ofNat_le_one,
    false_and]
  split_ifs with h1
  · refine (hzero fun h => ?_).symm
    unfold sdDomPair at h
    simp only [↓reduceIte] at h
    split_ifs at h <;> tauto
  · set c := min (sampleMin x) (sampleMin y) / 2 + max (sampleMax x) (sampleMax y) / 2
    rcases sdDomPair_eq_zero_or_one x y 1 false with h | h
    · rcases sdDomPair_eq_zero_or_one (shiftSample x c) (shiftSample y c) 1 false with h' | h'
      · rw [h, h']
      · rw [hshift] at h'; omega
    · rw [h]; exact (hshift c).mpr h

theorem sdDomPairCentred_eq_one_iff (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ)
    (y : Fin m → ℝ) {d : ℕ} (hd : d = 1 ∨ d = 2 ∨ d = 3) (b : Bool) :
    sdDomPairCentred x y d b = 1 ↔
      if d = 1 then
        if b then LPMDominates 0 x y
        else (∀ t, LPMratio 1 t x ≤ LPMratio 1 t y) ∧ ∃ t, LPMratio 1 t x ≠ LPMratio 1 t y
      else LPMDominates (d - 1) x y := by
  rcases hd with rfl | hd
  · cases b with
    | false =>
      rw [sdDomPairCentred_omega_eq hn hm, sdDomPair_omega_iff hn hm]
      simp
    | true =>
      rw [sdDomPairCentred_eq hn hm x y (Or.inl rfl) true fun _ => rfl,
        sdDomPair_fsd_iff hn hm]
      simp [LPMDominates]
  · have hd1 : d ≠ 1 := by rcases hd with rfl | rfl <;> norm_num
    simp only [hd1, ↓reduceIte]
    rw [sdDomPairCentred_eq hn hm x y (Or.inr hd) b fun h => absurd h hd1]
    exact sdDomPair_eq_one_iff_dominates hn hm x y hd b

end

end NNS
