-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.DominanceAlgo
public import NNSProofs.MeasurePM

namespace NNS.GapSuperiority

open MeasureTheory ProbabilityTheory Set NNS

@[expose] public noncomputable section

def pGt (μ ν : Measure ℝ) : ℝ := (μ.prod ν).real {p | p.2 < p.1}

def pTie (μ ν : Measure ℝ) : ℝ := (μ.prod ν).real {p | p.1 = p.2}

def pStar (μ ν : Measure ℝ) : ℝ := pGt μ ν + pTie μ ν / 2

def cSS (μ ν : Measure ℝ) : ℝ := |2 * pStar μ ν - 1|

lemma measurableSet_gt : MeasurableSet {p : ℝ × ℝ | p.2 < p.1} :=
  measurableSet_lt measurable_snd measurable_fst

lemma measurableSet_lt' : MeasurableSet {p : ℝ × ℝ | p.1 < p.2} :=
  measurableSet_lt measurable_fst measurable_snd

lemma measurableSet_diag : MeasurableSet {p : ℝ × ℝ | p.1 = p.2} :=
  measurableSet_eq_fun measurable_fst measurable_snd

variable {μ ν : Measure ℝ}

lemma prod_swap_real [SFinite μ] [SFinite ν] {s : Set (ℝ × ℝ)} (hs : MeasurableSet s) :
    (ν.prod μ).real s = (μ.prod ν).real (Prod.swap ⁻¹' s) := by
  rw [← Measure.prod_swap, map_measureReal_apply measurable_swap hs]

lemma pGt_swap [SFinite μ] [SFinite ν] : pGt ν μ = (μ.prod ν).real {p | p.1 < p.2} := by
  rw [pGt, prod_swap_real measurableSet_gt]
  rfl

lemma pTie_comm [SFinite μ] [SFinite ν] : pTie ν μ = pTie μ ν := by
  rw [pTie, prod_swap_real measurableSet_diag, pTie]
  congr 1
  ext p
  exact eq_comm

theorem pGt_add_pGt_add_pTie [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    pGt μ ν + pGt ν μ + pTie μ ν = 1 := by
  rw [pGt_swap (μ := μ) (ν := ν), pGt, pTie]
  have h1 : Disjoint {p : ℝ × ℝ | p.2 < p.1} {p | p.1 < p.2} :=
    Set.disjoint_left.mpr fun p (h : p.2 < p.1) (h' : p.1 < p.2) => lt_asymm h h'
  have h2 : Disjoint ({p : ℝ × ℝ | p.2 < p.1} ∪ {p | p.1 < p.2}) {p | p.1 = p.2} :=
    Set.disjoint_left.mpr fun p h h' => by
      rcases h with h | h
      · exact (ne_of_lt h) (h'.symm)
      · exact (ne_of_lt h) h'
  rw [← measureReal_union h1 measurableSet_lt',
    ← measureReal_union h2 measurableSet_diag]
  have : {p : ℝ × ℝ | p.2 < p.1} ∪ {p | p.1 < p.2} ∪ {p | p.1 = p.2} = univ := by
    ext p
    simp only [mem_union, mem_ofPred_eq, mem_univ, iff_true]
    rcases lt_trichotomy p.1 p.2 with h | h | h
    · exact Or.inl (Or.inr h)
    · exact Or.inr h
    · exact Or.inl (Or.inl h)
  rw [this, probReal_univ]

theorem pStar_add_pStar [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    pStar μ ν + pStar ν μ = 1 := by
  have h := pGt_add_pGt_add_pTie (μ := μ) (ν := ν)
  have hc : pTie ν μ = pTie μ ν := pTie_comm
  rw [pStar, pStar, hc]
  linarith

theorem pStar_self [IsProbabilityMeasure μ] : pStar μ μ = 1 / 2 := by
  have h := pStar_add_pStar (μ := μ) (ν := μ)
  linarith

theorem pStar_sub_half [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    pStar μ ν - 1 / 2 = (pGt μ ν - pGt ν μ) / 2 := by
  have h := pGt_add_pGt_add_pTie (μ := μ) (ν := ν)
  rw [pStar]
  linarith

theorem half_lt_pStar_iff [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    1 / 2 < pStar μ ν ↔ pGt ν μ < pGt μ ν := by
  have h := pStar_sub_half (μ := μ) (ν := ν)
  constructor <;> intro h' <;> linarith

theorem pStar_eq_half_iff [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    pStar μ ν = 1 / 2 ↔ pGt μ ν = pGt ν μ := by
  have h := pStar_sub_half (μ := μ) (ν := ν)
  constructor <;> intro h' <;> linarith

lemma pGt_nonneg : 0 ≤ pGt μ ν := measureReal_nonneg

lemma pTie_nonneg : 0 ≤ pTie μ ν := measureReal_nonneg

theorem pStar_mem_Icc [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    pStar μ ν ∈ Icc (0 : ℝ) 1 := by
  have h := pGt_add_pGt_add_pTie (μ := μ) (ν := ν)
  have h1 := pGt_nonneg (μ := μ) (ν := ν)
  have h2 := pGt_nonneg (μ := ν) (ν := μ)
  have h3 := pTie_nonneg (μ := μ) (ν := ν)
  rw [pStar]
  constructor <;> linarith

theorem cSS_eq [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    cSS μ ν = |pGt μ ν - pGt ν μ| := by
  have h := pStar_sub_half (μ := μ) (ν := ν)
  rw [cSS, show 2 * pStar μ ν - 1 = pGt μ ν - pGt ν μ by linarith]

theorem cSS_mem_Icc [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    cSS μ ν ∈ Icc (0 : ℝ) 1 := by
  have h := pStar_mem_Icc (μ := μ) (ν := ν)
  refine ⟨abs_nonneg _, abs_le.mpr ⟨?_, ?_⟩⟩ <;> linarith [h.1, h.2]

theorem cSS_eq_zero_iff [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    cSS μ ν = 0 ↔ pGt μ ν = pGt ν μ := by
  rw [cSS_eq, abs_eq_zero, sub_eq_zero]

theorem cSS_eq_one_iff [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    cSS μ ν = 1 ↔ pGt μ ν = 1 ∨ pGt ν μ = 1 := by
  have h := pGt_add_pGt_add_pTie (μ := μ) (ν := ν)
  have h1 := pGt_nonneg (μ := μ) (ν := ν)
  have h2 := pGt_nonneg (μ := ν) (ν := μ)
  have h3 := pTie_nonneg (μ := μ) (ν := ν)
  rw [cSS_eq]
  constructor
  · intro habs
    rcases abs_eq (zero_le_one) |>.mp habs with h' | h'
    · left; linarith
    · right; linarith
  · rintro (h' | h')
    · rw [abs_eq zero_le_one]; left; linarith
    · rw [abs_eq zero_le_one]; right; linarith

theorem pTie_eq_zero_left [NullSingletonClass μ] [SFinite μ] [SFinite ν] : pTie μ ν = 0 := by
  rw [pTie, measureReal_def, Measure.prod_apply_symm measurableSet_diag]
  simp

theorem pTie_eq_zero_right [NullSingletonClass ν] [SFinite μ] [SFinite ν] : pTie μ ν = 0 := by
  rw [← pTie_comm]
  exact pTie_eq_zero_left

theorem pGt_add_pGt_of_noAtoms [IsProbabilityMeasure μ] [IsProbabilityMeasure ν]
    (h : NullSingletonClass μ ∨ NullSingletonClass ν) : pGt μ ν + pGt ν μ = 1 := by
  have h0 : pTie μ ν = 0 := h.elim (fun _ => pTie_eq_zero_left) fun _ => pTie_eq_zero_right
  linarith [pGt_add_pGt_add_pTie (μ := μ) (ν := ν)]

section RV

variable {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P]
  {X Y : Ω → ℝ}

lemma map_pair (hX : Measurable X) (hY : Measurable Y) (hXY : IndepFun X Y P) :
    P.map (fun ω => (X ω, Y ω)) = (P.map X).prod (P.map Y) :=
  (indepFun_iff_map_prod_eq_prod_map_map hX.aemeasurable hY.aemeasurable).mp hXY

theorem measureReal_gt_eq_pGt (hX : Measurable X) (hY : Measurable Y) (hXY : IndepFun X Y P) :
    P.real {ω | Y ω < X ω} = pGt (P.map X) (P.map Y) := by
  rw [pGt, ← map_pair hX hY hXY, map_measureReal_apply (hX.prodMk hY) measurableSet_gt]
  rfl

theorem measureReal_eq_eq_pTie (hX : Measurable X) (hY : Measurable Y) (hXY : IndepFun X Y P) :
    P.real {ω | X ω = Y ω} = pTie (P.map X) (P.map Y) := by
  rw [pTie, ← map_pair hX hY hXY, map_measureReal_apply (hX.prodMk hY) measurableSet_diag]
  rfl

end RV

section Sample

open Finset

variable {n m : ℕ}

def leftPtr (y : Fin m → ℝ) (v : ℝ) : ℕ :=
  ((sortedVals y).takeWhile fun w => decide (w < v)).length

def lessCount (x : Fin n → ℝ) (y : Fin m → ℝ) : ℕ :=
  ((sortedVals x).map fun v => leftPtr y v).sum

def tieCount (x : Fin n → ℝ) (y : Fin m → ℝ) : ℕ :=
  ((sortedVals x).map fun v => walkCount y v - leftPtr y v).sum

def ssGt (x : Fin n → ℝ) (y : Fin m → ℝ) : ℝ := (lessCount x y : ℝ) / ((n : ℝ) * m)

def ssTie (x : Fin n → ℝ) (y : Fin m → ℝ) : ℝ := (tieCount x y : ℝ) / ((n : ℝ) * m)

def ssStar (x : Fin n → ℝ) (y : Fin m → ℝ) : ℝ := ssGt x y + 1 / 2 * ssTie x y

lemma takeWhile_lt_eq_filter {l : List ℝ} (hl : l.Pairwise (· ≤ ·)) (t : ℝ) :
    l.takeWhile (fun v => decide (v < t)) = l.filter (fun v => decide (v < t)) := by
  induction l with
  | nil => rfl
  | cons a l ih =>
    rw [List.pairwise_cons] at hl
    by_cases ha : a < t
    · rw [List.takeWhile_cons_of_pos (by simpa using ha),
        List.filter_cons_of_pos (by simpa using ha), ih hl.2]
    · rw [List.takeWhile_cons_of_neg (by simpa using ha),
        List.filter_cons_of_neg (by simpa using ha), eq_comm, List.filter_eq_nil_iff]
      intro b hb
      simp only [decide_eq_true_eq]
      intro hbt
      exact ha ((hl.1 b hb).trans_lt hbt)

lemma leftPtr_eq (y : Fin m → ℝ) (v : ℝ) : leftPtr y v = #{j | y j < v} := by
  rw [leftPtr, takeWhile_lt_eq_filter (sortedVals_pairwise y), ← Multiset.coe_card,
    ← Multiset.filter_coe, coe_sortedVals, Multiset.filter_map, Multiset.card_map, card_def,
    filter_val]
  rfl

lemma walkCount_eq (y : Fin m → ℝ) (v : ℝ) : walkCount y v = #{j | y j ≤ v} :=
  walkCount_eq_lowerCount y v

lemma tie_eq (y : Fin m → ℝ) (v : ℝ) : walkCount y v - leftPtr y v = #{j | y j = v} := by
  rw [walkCount_eq, leftPtr_eq]
  have h : ({j | y j ≤ v} : Finset (Fin m)) =
      ({j | y j < v} : Finset (Fin m)) ∪ ({j | y j = v} : Finset (Fin m)) := by
    ext j; simp [le_iff_lt_or_eq]
  have hd : Disjoint ({j | y j < v} : Finset (Fin m)) ({j | y j = v} : Finset (Fin m)) :=
    disjoint_filter.mpr fun j _ h1 h2 => (ne_of_lt h1) h2
  rw [h, card_union_of_disjoint hd, Nat.add_sub_cancel_left]

lemma sum_sortedVals (x : Fin n → ℝ) (f : ℝ → ℕ) :
    ((sortedVals x).map f).sum = ∑ i, f (x i) := by
  rw [sum_eq_multiset_sum, ← Multiset.sum_coe, ← Multiset.map_coe, coe_sortedVals,
    Multiset.map_map]
  rfl

theorem ssGt_eq_sum (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ssGt x y = (∑ i, ∑ j, if y j < x i then (1 : ℝ) else 0) / ((n : ℝ) * m) := by
  rw [ssGt, lessCount, sum_sortedVals]
  simp only [leftPtr_eq, Finset.card_filter]
  push_cast
  rfl

theorem ssTie_eq_sum (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ssTie x y = (∑ i, ∑ j, if x i = y j then (1 : ℝ) else 0) / ((n : ℝ) * m) := by
  rw [ssTie, tieCount, sum_sortedVals]
  simp only [tie_eq, Finset.card_filter]
  push_cast
  congr 1
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
  by_cases h : x i = y j
  · simp [h]
  · simp [h, Ne.symm h]

theorem ssStar_eq_sum (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ssStar x y = (∑ i, ∑ j, ((if y j < x i then (1 : ℝ) else 0) +
      1 / 2 * (if x i = y j then (1 : ℝ) else 0))) / ((n : ℝ) * m) := by
  rw [ssStar, ssGt_eq_sum, ssTie_eq_sum, mul_div_assoc', ← add_div]
  congr 1
  rw [Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.mul_sum, ← Finset.sum_add_distrib]

theorem ssGt_eq_UPM_zero (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ssGt x y = (∑ j, UPM 0 (y j) x) / m := by
  rw [ssGt_eq_sum, Finset.sum_comm]
  simp only [UPM, upperDev, pow_zero, ← Finset.sum_div, div_div]

lemma measureReal_prod_emp (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    {s : Set (ℝ × ℝ)} (hs : MeasurableSet s) :
    ((empMeasure x).prod (empMeasure y)).real s =
      (∑ i, ∑ j, s.indicator (fun _ => (1 : ℝ)) (x i, y j)) / ((n : ℝ) * m) := by
  have := isProbabilityMeasure_empMeasure hn x
  have := isProbabilityMeasure_empMeasure hm y
  rw [← integral_indicator_one hs, integral_prod]
  · simp only [integral_empMeasure, ← Finset.sum_div, div_div]
    rw [mul_comm]
    rfl
  · exact (integrable_const (1 : ℝ)).indicator hs

theorem ssGt_eq_pGt (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ssGt x y = pGt (empMeasure x) (empMeasure y) := by
  rw [ssGt_eq_sum, pGt, measureReal_prod_emp hn hm x y measurableSet_gt]
  rfl

theorem ssTie_eq_pTie (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ssTie x y = pTie (empMeasure x) (empMeasure y) := by
  rw [ssTie_eq_sum, pTie, measureReal_prod_emp hn hm x y measurableSet_diag]
  rfl

theorem ssStar_eq_pStar (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ssStar x y = pStar (empMeasure x) (empMeasure y) := by
  rw [ssStar, pStar, ssGt_eq_pGt hn hm, ssTie_eq_pTie hn hm]
  ring

end Sample

section Reflection

open Finset

variable {n : ℕ}

lemma measurable_mean : Measurable fun x : Fin n → ℝ => mean x := by
  unfold mean
  fun_prop

def lowerSet (i : Fin n) : Set (Fin n → ℝ) := {x | x i ≤ mean x}

lemma measurableSet_lowerSet (i : Fin n) : MeasurableSet (lowerSet i) :=
  measurableSet_le (measurable_pi_apply i) measurable_mean

lemma LPM_zero_mean_eq (x : Fin n → ℝ) :
    LPM 0 (mean x) x = (∑ i, (lowerSet i).indicator (1 : (Fin n → ℝ) → ℝ) x) / n := by
  unfold LPM lowerDev lowerSet
  congr 1

lemma tie_null (μ : Measure ℝ) [IsProbabilityMeasure μ] [NullSingletonClass μ]
    (hn : 2 ≤ n) (i : Fin n) :
    (Measure.pi fun _ : Fin n => μ) {x | x i = mean x} = 0 := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
  have hmp := measurePreserving_piFinSuccAbove (fun _ : Fin (k + 1) => μ) i
  set e := MeasurableEquiv.piFinSuccAbove (fun _ : Fin (k + 1) => ℝ) i
  have hE : MeasurableSet {x : Fin (k + 1) → ℝ | x i = mean x} :=
    measurableSet_eq_fun (measurable_pi_apply i) measurable_mean
  have hpre : {x : Fin (k + 1) → ℝ | x i = mean x} = e ⁻¹' (e.symm ⁻¹' {x | x i = mean x}) := by
    ext x; simp
  rw [hpre, hmp.measure_preimage_equiv, Measure.prod_apply_symm (e.symm.measurable hE)]
  have hk : (1 : ℝ) ≤ k := by exact_mod_cast (show 1 ≤ k by omega)
  have hslice : ∀ z : Fin k → ℝ,
      μ ((fun a => (a, z)) ⁻¹' (e.symm ⁻¹' {x | x i = mean x})) = 0 := fun z => by
    refine Set.Subsingleton.measure_zero (fun a ha b hb => ?_) μ
    simp only [Set.mem_preimage, Set.mem_ofPred_eq, e, MeasurableEquiv.piFinSuccAbove_symm_apply,
      Fin.insertNthEquiv_apply, mean, Fin.sum_univ_succAbove _ i, Fin.insertNth_apply_same,
      Fin.insertNth_apply_succAbove] at ha hb
    have hk1 : ((k + 1 : ℕ) : ℝ) ≠ 0 := by positivity
    field_simp at ha hb
    push_cast at ha hb
    nlinarith
  exact (lintegral_congr fun z => hslice z).trans lintegral_zero

theorem integral_ecdf_mean_eq_half (μ : Measure ℝ) [IsProbabilityMeasure μ]
    [NullSingletonClass μ] {c : ℝ} (hsym : μ.map (fun v => 2 * c - v) = μ) (hn : 2 ≤ n) :
    ∫ x, LPM 0 (mean x) x ∂(Measure.pi fun _ : Fin n => μ) = 1 / 2 := by
  set P := Measure.pi fun _ : Fin n => μ
  have hn0 : (n : ℝ) ≠ 0 := by positivity
  -- reflection `x ↦ 2c - x` preserves `P` and swaps `x_i ≤ x̄` with `x̄ ≤ x_i`
  have hR : MeasurePreserving (fun v : ℝ => 2 * c - v) μ μ := ⟨by fun_prop, hsym⟩
  have hRpi := measurePreserving_pi (fun _ : Fin n => μ) (fun _ => μ) (fun _ => hR)
  have hmeanR : ∀ x : Fin n → ℝ, mean (fun j => 2 * c - x j) = 2 * c - mean x := fun x => by
    simp only [mean, Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ,
      Fintype.card_fin, nsmul_eq_mul]
    field_simp
  have hhalf : ∀ i : Fin n, P.real (lowerSet i) = 1 / 2 := fun i => by
    have hB : MeasurableSet {x : Fin n → ℝ | mean x ≤ x i} :=
      measurableSet_le measurable_mean (measurable_pi_apply i)
    have hAB : P.real (lowerSet i) = P.real {x | mean x ≤ x i} := by
      rw [← hRpi.measureReal_preimage (measurableSet_lowerSet i).nullMeasurableSet]
      congr 1
      ext x
      simp only [Set.mem_preimage, lowerSet, Set.mem_ofPred_eq, hmeanR]
      constructor <;> intro h <;> linarith
    have hunion : lowerSet i ∪ {x : Fin n → ℝ | mean x ≤ x i} = univ := by
      ext x; simp only [lowerSet, Set.mem_union, Set.mem_ofPred_eq, mem_univ, iff_true]
      exact le_total _ _
    have hinter : lowerSet i ∩ {x : Fin n → ℝ | mean x ≤ x i} = {x | x i = mean x} := by
      ext x; simp only [lowerSet, Set.mem_inter_iff, Set.mem_ofPred_eq]
      exact ⟨fun h => le_antisymm h.1 h.2, fun h => ⟨h.le, h.ge⟩⟩
    have h := measureReal_union_add_inter (μ := P) (s := lowerSet i) hB
    rw [hunion, hinter, probReal_univ, measureReal_def, tie_null μ hn i] at h
    simp only [ENNReal.toReal_zero, add_zero] at h
    linarith
  simp_rw [LPM_zero_mean_eq]
  have hint : ∀ i ∈ (Finset.univ : Finset (Fin n)),
      Integrable ((lowerSet i).indicator (1 : (Fin n → ℝ) → ℝ)) P := fun i _ =>
    (integrable_const (1 : ℝ)).indicator (measurableSet_lowerSet i)
  rw [integral_div, integral_finsetSum _ hint]
  simp only [integral_indicator_one (measurableSet_lowerSet _), hhalf, Finset.sum_const,
    Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  field_simp

theorem integral_ecdf_mean_one (μ : Measure ℝ) [IsProbabilityMeasure μ] :
    ∫ x, LPM 0 (mean x) x ∂(Measure.pi fun _ : Fin 1 => μ) = 1 := by
  have : ∀ x : Fin 1 → ℝ, LPM 0 (mean x) x = 1 := fun x => by
    simp [LPM, lowerDev, mean, Finset.filter_singleton]
  simp [this]

end Reflection

end

end NNS.GapSuperiority
