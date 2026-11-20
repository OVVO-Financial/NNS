-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.DominanceAlgo
public import NNSProofs.GapQuantile

namespace NNS.GapDominance

open Finset NNS

@[expose] public noncomputable section

variable {n m : ℕ}

def CurveDom (f g : ℝ → ℝ) : Prop := (∀ t, f t ≤ g t) ∧ ∃ t, f t ≠ g t

theorem curveDom_irrefl (f : ℝ → ℝ) : ¬ CurveDom f f := fun ⟨_, _, h⟩ => h rfl

theorem curveDom_trans {f g h : ℝ → ℝ} (hfg : CurveDom f g) (hgh : CurveDom g h) :
    CurveDom f h := by
  refine ⟨fun t => (hfg.1 t).trans (hgh.1 t), ?_⟩
  obtain ⟨t, ht⟩ := hfg.2
  exact ⟨t, fun he => ht (le_antisymm (hfg.1 t) (he ▸ hgh.1 t))⟩

theorem curveDom_asymm {f g : ℝ → ℝ} (hfg : CurveDom f g) : ¬ CurveDom g f :=
  fun hgf => curveDom_irrefl f (curveDom_trans hfg hgf)

def domCurve (d : ℕ) (b : Bool) (x : Fin n → ℝ) : ℝ → ℝ :=
  if d = 1 then (if b then fun t => LPM 0 t x else fun t => LPMratio 1 t x)
  else fun t => LPM (d - 1) t x

theorem fsdUni_curveDom (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    {d : ℕ} (hd : d = 1 ∨ d = 2 ∨ d = 3) (b : Bool) :
    sdDomPairCentred x y d b = 1 ↔ CurveDom (domCurve d b x) (domCurve d b y) := by
  rw [sdDomPairCentred_eq_one_iff hn hm x y hd b]
  unfold domCurve CurveDom LPMDominates
  split_ifs <;> rfl

theorem fsd_iff_lt (x : Fin n → ℝ) (y : Fin m → ℝ) :
    LPMDominates 0 x y ↔ (∀ t, LPM 0 t x ≤ LPM 0 t y) ∧ ∃ t, LPM 0 t x < LPM 0 t y := by
  refine ⟨fun ⟨h, t, ht⟩ => ⟨h, t, lt_of_le_of_ne (h t) ht⟩,
    fun ⟨h, t, ht⟩ => ⟨h, t, ht.ne⟩⟩

theorem fsdUni_eq_one_iff (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    sdDomPairCentred x y 1 true = 1 ↔
      (∀ t, (#{i | x i ≤ t} : ℝ) / n ≤ (#{j | y j ≤ t} : ℝ) / m) ∧
        ∃ t, (#{i | x i ≤ t} : ℝ) / n < (#{j | y j ≤ t} : ℝ) / m := by
  rw [sdDomPairCentred_eq_one_iff hn hm x y (Or.inl rfl) true]
  simp only [↓reduceIte]
  rw [fsd_iff_lt]
  simp only [LPM_zero_eq_ecdf]

section Efficient

variable {k r : ℕ}

def domMatrix (X : Fin k → Fin r → ℝ) (d : ℕ) (b : Bool) (i j : Fin k) : ℕ :=
  if i = j then 0 else sdDomPairCentred (X i) (X j) d b

def edge (D : Fin k → Fin k → ℕ) (R : Finset (Fin k)) (a b : Fin k) : Prop :=
  a ∈ R ∧ b ∈ R ∧ D a b = 1

def sameComp (D : Fin k → Fin k → ℕ) (R : Finset (Fin k)) (i j : Fin k) : Prop :=
  i = j ∨ (Relation.TransGen (edge D R) i j ∧ Relation.TransGen (edge D R) j i)

def sccKeep (D : Fin k → Fin k → ℕ) (R : Finset (Fin k)) (v : Fin k) : Prop :=
  ¬ ∃ i j, edge D R i j ∧ sameComp D R j v ∧ ¬ sameComp D R i j

open Classical in
def effSet (D : Fin k → Fin k → ℕ) (R : Finset (Fin k)) : Finset (Fin k) :=
  R.filter (sccKeep D R)

lemma domMatrix_trans (hr : 0 < r) (X : Fin k → Fin r → ℝ) {d : ℕ}
    (hd : d = 1 ∨ d = 2 ∨ d = 3) (b : Bool) {a c e : Fin k}
    (hac : domMatrix X d b a c = 1) (hce : domMatrix X d b c e = 1) :
    domMatrix X d b a e = 1 := by
  unfold domMatrix at *
  have hac' : sdDomPairCentred (X a) (X c) d b = 1 := by
    by_cases h : a = c
    · simp [h] at hac
    · simpa [h] using hac
  have hce' : sdDomPairCentred (X c) (X e) d b = 1 := by
    by_cases h : c = e
    · simp [h] at hce
    · simpa [h] using hce
  rw [fsdUni_curveDom hr hr _ _ hd] at hac' hce'
  by_cases hae : a = e
  · subst hae
    exact absurd hce' (curveDom_asymm hac')
  simp only [hae, ↓reduceIte]
  rw [fsdUni_curveDom hr hr _ _ hd]
  exact curveDom_trans hac' hce'

lemma domMatrix_irrefl (X : Fin k → Fin r → ℝ) (d : ℕ) (b : Bool) (a : Fin k) :
    domMatrix X d b a a ≠ 1 := by
  simp [domMatrix]

lemma transGen_edge (hr : 0 < r) (X : Fin k → Fin r → ℝ) {d : ℕ}
    (hd : d = 1 ∨ d = 2 ∨ d = 3) (b : Bool) {R : Finset (Fin k)} {a c : Fin k}
    (h : Relation.TransGen (edge (domMatrix X d b) R) a c) : edge (domMatrix X d b) R a c := by
  induction h with
  | single h => exact h
  | tail _ hce ih => exact ⟨ih.1, hce.2.1, domMatrix_trans hr X hd b ih.2.2 hce.2.2⟩

lemma sameComp_iff (hr : 0 < r) (X : Fin k → Fin r → ℝ) {d : ℕ}
    (hd : d = 1 ∨ d = 2 ∨ d = 3) (b : Bool) (R : Finset (Fin k)) (i j : Fin k) :
    sameComp (domMatrix X d b) R i j ↔ i = j := by
  refine ⟨fun h => h.elim id fun ⟨h1, h2⟩ => ?_, Or.inl⟩
  have e1 := transGen_edge hr X hd b h1
  have e2 := transGen_edge hr X hd b h2
  exact absurd (domMatrix_trans hr X hd b e1.2.2 e2.2.2) (domMatrix_irrefl X d b i)

theorem effSet_eq (hr : 0 < r) (X : Fin k → Fin r → ℝ) {d : ℕ}
    (hd : d = 1 ∨ d = 2 ∨ d = 3) (b : Bool) (R : Finset (Fin k)) (v : Fin k) :
    v ∈ effSet (domMatrix X d b) R ↔
      v ∈ R ∧ ∀ j ∈ R, sdDomPairCentred (X j) (X v) d b ≠ 1 := by
  simp only [effSet, mem_filter, sccKeep, sameComp_iff hr X hd b, edge]
  refine and_congr_right fun hv => ⟨fun h j hj hjv => ?_, fun h => ?_⟩
  · have hne : j ≠ v := by
      rintro rfl
      exact curveDom_irrefl _ ((fsdUni_curveDom hr hr _ _ hd b).mp hjv)
    exact h ⟨j, v, ⟨hj, hv, by simp [domMatrix, hne, hjv]⟩, rfl, hne⟩
  · rintro ⟨i, j, ⟨hi, -, hij⟩, rfl, hne⟩
    exact h i hi (by simpa [domMatrix, hne] using hij)

theorem exists_eff_dominates (hr : 0 < r) (X : Fin k → Fin r → ℝ) {d : ℕ}
    (hd : d = 1 ∨ d = 2 ∨ d = 3) (b : Bool) {R : Finset (Fin k)} {v : Fin k} (hv : v ∈ R)
    (hnot : v ∉ effSet (domMatrix X d b) R) :
    ∃ e ∈ effSet (domMatrix X d b) R, sdDomPairCentred (X e) (X v) d b = 1 := by
  let rel : Fin k → Fin k → Prop := fun a c => domMatrix X d b a c = 1
  have : IsTrans (Fin k) rel := ⟨fun _ _ _ h1 h2 => domMatrix_trans hr X hd b h1 h2⟩
  have : Std.Irrefl rel := ⟨fun a => domMatrix_irrefl X d b a⟩
  have wf := Finite.wellFounded_of_trans_of_irrefl rel
  rw [effSet_eq hr X hd b] at hnot
  push Not at hnot
  obtain ⟨j, hj, hjv⟩ := hnot hv
  have hS : ({a | a ∈ R ∧ rel a v} : Set (Fin k)).Nonempty := by
    have hne : j ≠ v := by
      rintro rfl
      exact curveDom_irrefl _ ((fsdUni_curveDom hr hr _ _ hd b).mp hjv)
    exact ⟨j, hj, by simp [rel, domMatrix, hne, hjv]⟩
  obtain ⟨e, ⟨he, hev⟩, hmin⟩ := wf.has_min _ hS
  refine ⟨e, (effSet_eq hr X hd b R e).mpr ⟨he, fun c hc hce => ?_⟩, ?_⟩
  · have hne : c ≠ e := by
      rintro rfl
      exact curveDom_irrefl _ ((fsdUni_curveDom hr hr _ _ hd b).mp hce)
    have hce' : rel c e := by simp [rel, domMatrix, hne, hce]
    exact hmin c ⟨hc, domMatrix_trans hr X hd b hce' hev⟩ hce'
  · have hne : e ≠ v := by rintro rfl; exact domMatrix_irrefl X d b e hev
    simpa [rel, domMatrix, hne] using hev

def Acyclic (D : Fin k → Fin k → ℕ) : Prop :=
  ∀ i, ¬ Relation.TransGen (fun a b => D a b = 1) i i

lemma transGen_edge_le {D : Fin k → Fin k → ℕ} {R : Finset (Fin k)} {a b : Fin k}
    (h : Relation.TransGen (edge D R) a b) : Relation.TransGen (fun a b => D a b = 1) a b := by
  induction h with
  | single h => exact .single h.2.2
  | tail _ h ih => exact ih.tail h.2.2

lemma sameComp_iff_of_acyclic {D : Fin k → Fin k → ℕ} (hD : Acyclic D) (R : Finset (Fin k))
    (i j : Fin k) : sameComp D R i j ↔ i = j :=
  ⟨fun h => h.elim id fun ⟨h1, h2⟩ => absurd (transGen_edge_le (h1.trans h2)) (hD i), Or.inl⟩

theorem effSet_eq_of_acyclic {D : Fin k → Fin k → ℕ} (hD : Acyclic D) (R : Finset (Fin k))
    (v : Fin k) : v ∈ effSet D R ↔ v ∈ R ∧ ∀ j ∈ R, D j v ≠ 1 := by
  simp only [effSet, mem_filter, sccKeep, sameComp_iff_of_acyclic hD, edge]
  refine and_congr_right fun hv => ⟨fun h j hj hjv => ?_, fun h => ?_⟩
  · have hne : j ≠ v := by
      rintro rfl
      exact hD j (Relation.TransGen.single hjv)
    exact h ⟨j, v, ⟨hj, hv, hjv⟩, rfl, hne⟩
  · rintro ⟨i, j, ⟨hi, -, hij⟩, rfl, -⟩
    exact h i hi hij

theorem effSet_antichain_of_acyclic {D : Fin k → Fin k → ℕ} (hD : Acyclic D) (R : Finset (Fin k))
    {i j : Fin k} (hi : i ∈ effSet D R) (hj : j ∈ effSet D R) : D i j ≠ 1 :=
  ((effSet_eq_of_acyclic hD R j).mp hj).2 i ((effSet_eq_of_acyclic hD R i).mp hi).1

theorem acyclic_domMatrix (hr : 0 < r) (X : Fin k → Fin r → ℝ) {d : ℕ}
    (hd : d = 1 ∨ d = 2 ∨ d = 3) (b : Bool) : Acyclic (domMatrix X d b) := by
  intro i h
  have key : ∀ {a c : Fin k}, Relation.TransGen (fun u w => domMatrix X d b u w = 1) a c →
      domMatrix X d b a c = 1 := fun h => by
    induction h with
    | single h => exact h
    | tail _ hce ih => exact domMatrix_trans hr X hd b ih hce
  exact domMatrix_irrefl X d b i (key h)

end Efficient

section Utility

def stepSum (u : ℝ → ℝ) : List ℝ → ℝ → ℝ
  | a :: b :: l, v => (u b - u a) * (if b ≤ v then 1 else 0) + stepSum u (b :: l) v
  | _, _ => 0

def stepW (u : ℝ → ℝ) : List ℝ → (ℝ → ℝ) → ℝ
  | a :: b :: l, g => (u b - u a) * g b + stepW u (b :: l) g
  | _, _ => 0

lemma stepSum_eq_zero (u : ℝ → ℝ) {l : List ℝ} {v : ℝ} (h : ∀ s ∈ l, v < s) :
    stepSum u l v = 0 := by
  induction l with
  | nil => rfl
  | cons a l ih =>
    cases l with
    | nil => rfl
    | cons b l =>
      have hb : ¬ b ≤ v := not_le.mpr (h b (by simp))
      simp only [stepSum, hb, ↓reduceIte, mul_zero, zero_add]
      exact ih fun s hs => h s (List.mem_cons_of_mem a hs)

lemma eq_add_stepSum (u : ℝ → ℝ) {a : ℝ} {l : List ℝ} (hl : (a :: l).Pairwise (· < ·))
    {v : ℝ} (hv : v ∈ a :: l) : u v = u a + stepSum u (a :: l) v := by
  induction l generalizing a with
  | nil =>
    rw [List.mem_singleton] at hv
    subst hv
    simp [stepSum]
  | cons b l ih =>
    have hab : a < b := List.rel_of_pairwise_cons hl (by simp)
    have hl' := hl.of_cons
    rcases List.mem_cons.mp hv with rfl | hv'
    · have hb : ¬ b ≤ v := not_le.mpr hab
      have hz : stepSum u (b :: l) v = 0 := stepSum_eq_zero u fun s hs => by
        rcases List.mem_cons.mp hs with rfl | hs
        · exact hab
        · exact List.rel_of_pairwise_cons hl (List.mem_cons_of_mem b hs)
      simp [stepSum, hb, hz]
    · have hbv : b ≤ v := by
        rcases List.mem_cons.mp hv' with rfl | hv''
        · exact le_rfl
        · exact (List.rel_of_pairwise_cons hl' hv'').le
      simp only [stepSum, hbv, ↓reduceIte]
      rw [ih hl' hv']
      ring

lemma sum_stepSum (u : ℝ → ℝ) (l : List ℝ) (x : Fin n → ℝ) :
    ∑ i, stepSum u l (x i) = stepW u l fun b => (#{i | b ≤ x i} : ℝ) := by
  induction l with
  | nil => simp [stepSum, stepW]
  | cons a l ih =>
    cases l with
    | nil => simp [stepSum, stepW]
    | cons b l =>
      simp only [stepSum, stepW, sum_add_distrib, ← mul_sum, ih]
      congr 2
      rw [sum_boole]

lemma stepW_mono (u : ℝ → ℝ) (hu : Monotone u) {l : List ℝ} (hl : l.Pairwise (· < ·))
    {g h : ℝ → ℝ} (hgh : ∀ b, g b ≤ h b) : stepW u l g ≤ stepW u l h := by
  induction l with
  | nil => exact le_rfl
  | cons a l ih =>
    cases l with
    | nil => exact le_rfl
    | cons b l =>
      have hab : u a ≤ u b := hu (List.rel_of_pairwise_cons hl (by simp)).le
      simp only [stepW]
      exact add_le_add (mul_le_mul_of_nonneg_left (hgh b) (by linarith)) (ih hl.of_cons)

lemma stepW_div (u : ℝ → ℝ) (l : List ℝ) (g : ℝ → ℝ) (c : ℝ) :
    stepW u l g / c = stepW u l fun b => g b / c := by
  induction l with
  | nil => simp [stepW]
  | cons a l ih =>
    cases l with
    | nil => simp [stepW]
    | cons b l => simp only [stepW, add_div, mul_div_assoc, ih]

lemma strict_count_le {x : Fin n → ℝ} {y : Fin m → ℝ}
    (h : ∀ t, LPM 0 t x ≤ LPM 0 t y) (b : ℝ) :
    (#{i | x i < b} : ℝ) / n ≤ (#{j | y j < b} : ℝ) / m := by
  set K := (univ.image x ∪ univ.image y).filter (· < b)
  have key : ∃ t, (∀ i, x i < b ↔ x i ≤ t) ∧ ∀ j, y j < b ↔ y j ≤ t := by
    by_cases hK : K.Nonempty
    · have hm := mem_filter.mp (K.max'_mem hK)
      refine ⟨K.max' hK, fun i => ⟨fun hi => K.le_max' _ ?_, fun hi => hi.trans_lt hm.2⟩,
        fun j => ⟨fun hj => K.le_max' _ ?_, fun hj => hj.trans_lt hm.2⟩⟩
      · exact mem_filter.mpr ⟨mem_union_left _ (mem_image_of_mem x (mem_univ i)), hi⟩
      · exact mem_filter.mpr ⟨mem_union_right _ (mem_image_of_mem y (mem_univ j)), hj⟩
    · refine ⟨b - 1, fun i => ⟨fun hi => absurd ⟨x i, ?_⟩ hK, fun hi => by linarith⟩,
        fun j => ⟨fun hj => absurd ⟨y j, ?_⟩ hK, fun hj => by linarith⟩⟩
      · exact mem_filter.mpr ⟨mem_union_left _ (mem_image_of_mem x (mem_univ i)), hi⟩
      · exact mem_filter.mpr ⟨mem_union_right _ (mem_image_of_mem y (mem_univ j)), hj⟩
  obtain ⟨t, hx, hy⟩ := key
  have := h t
  rwa [LPM_zero_eq_ecdf, LPM_zero_eq_ecdf, ← filter_congr fun i _ => hx i,
    ← filter_congr fun j _ => hy j] at this

lemma count_ge_div (hn : 0 < n) (x : Fin n → ℝ) (b : ℝ) :
    (#{i | b ≤ x i} : ℝ) / n = 1 - (#{i | x i < b} : ℝ) / n := by
  have hc : #{i | b ≤ x i} + #{i | x i < b} = n := by
    have := card_filter_add_card_filter_not (s := (univ : Finset (Fin n)))
      (fun i => b ≤ x i)
    simpa [not_le] using this
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  rw [eq_sub_iff_add_eq, ← add_div, ← Nat.cast_add, hc, div_self hn']

lemma mean_comp_eq (hn : 0 < n) (u : ℝ → ℝ) (x : Fin n → ℝ) {a : ℝ} {l : List ℝ}
    (hl : (a :: l).Pairwise (· < ·)) (hx : ∀ i, x i ∈ a :: l) :
    mean (fun i => u (x i)) = u a + stepW u (a :: l) fun b => (#{i | b ≤ x i} : ℝ) / n := by
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  rw [mean, sum_congr rfl fun i _ => eq_add_stepSum u hl (hx i), sum_add_distrib,
    sum_stepSum, ← stepW_div, add_div]
  simp [hn']

theorem mean_comp_le_of_LPM_zero_le (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ}
    {y : Fin m → ℝ} (h : ∀ t, LPM 0 t x ≤ LPM 0 t y) {u : ℝ → ℝ} (hu : Monotone u) :
    mean (fun j => u (y j)) ≤ mean (fun i => u (x i)) := by
  obtain ⟨a, l, hal⟩ : ∃ a l, knotList x y = a :: l := by
    rcases hk : knotList x y with _ | ⟨a, l⟩
    · have : x ⟨0, hn⟩ ∈ knotList x y := mem_knotList.mpr (Or.inl ⟨_, rfl⟩)
      simp [hk] at this
    · exact ⟨a, l, rfl⟩
  have hl : (a :: l).Pairwise (· < ·) := hal ▸ knotList_pairwise x y
  have hx : ∀ i, x i ∈ a :: l := fun i => hal ▸ mem_knotList.mpr (Or.inl ⟨i, rfl⟩)
  have hy : ∀ j, y j ∈ a :: l := fun j => hal ▸ mem_knotList.mpr (Or.inr ⟨j, rfl⟩)
  rw [mean_comp_eq hn u x hl hx, mean_comp_eq hm u y hl hy]
  refine add_le_add_right (stepW_mono u hu hl fun b => ?_) _
  rw [count_ge_div hn, count_ge_div hm]
  linarith [strict_count_le h b]

theorem mean_comp_le_of_fsdUni (hn : 0 < n) (hm : 0 < m) {x : Fin n → ℝ} {y : Fin m → ℝ}
    (h : sdDomPairCentred x y 1 true = 1) {u : ℝ → ℝ} (hu : Monotone u) :
    mean (fun j => u (y j)) ≤ mean (fun i => u (x i)) :=
  mean_comp_le_of_LPM_zero_le hn hm ((sdDomPairCentred_eq_one_iff hn hm x y (Or.inl rfl) true).mp
    h).1 hu

end Utility

section Laws

open MeasureTheory ProbabilityTheory GapQuantile

def qf (μ : Measure ℝ) (p : ℝ) : ℝ := if p ∈ Set.Ioo (0 : ℝ) 1 then quantile μ p else 0

lemma qf_preimage_Iic (μ : Measure ℝ) (t : ℝ) :
    qf μ ⁻¹' Set.Iic t =
      (Set.Ioo 0 1 ∩ {p | p ≤ cdf μ t}) ∪ ((Set.Ioo 0 1)ᶜ ∩ {_p | (0 : ℝ) ≤ t}) := by
  ext p
  by_cases hp : p ∈ Set.Ioo (0 : ℝ) 1
  · simp [qf, hp, quantile_le_iff μ hp.1 hp.2]
  · simp [qf, hp]

lemma measurable_qf (μ : Measure ℝ) : Measurable (qf μ) := by
  refine measurable_of_Iic fun t => ?_
  rw [qf_preimage_Iic]
  exact (measurableSet_Ioo.inter (measurableSet_le measurable_id measurable_const)).union
    (measurableSet_Ioo.compl.inter (MeasurableSet.const _))

lemma volume_Ioo_inter_Iic {a : ℝ} (h0 : 0 ≤ a) (h1 : a ≤ 1) :
    volume (Set.Ioo (0 : ℝ) 1 ∩ Set.Iic a) = ENNReal.ofReal a := by
  rcases h1.lt_or_eq with h1 | rfl
  · have : Set.Ioo (0 : ℝ) 1 ∩ Set.Iic a = Set.Ioc 0 a := by
      ext p
      simp only [Set.mem_inter_iff, Set.mem_Ioo, Set.mem_Iic, Set.mem_Ioc]
      exact ⟨fun h => ⟨h.1.1, h.2⟩, fun h => ⟨⟨h.1, h.2.trans_lt h1⟩, h.2⟩⟩
    rw [this, Real.volume_Ioc, sub_zero]
  · rw [Set.inter_eq_left.mpr fun _ h => Set.mem_Iic.mpr h.2.le, Real.volume_Ioo, sub_zero]

theorem map_qf (μ : Measure ℝ) [IsProbabilityMeasure μ] :
    (volume.restrict (Set.Ioo (0 : ℝ) 1)).map (qf μ) = μ := by
  refine (Measure.ext_of_Iic μ _ fun t => ?_).symm
  rw [Measure.map_apply (measurable_qf μ) measurableSet_Iic,
    Measure.restrict_apply (measurable_qf μ measurableSet_Iic), ← ofReal_cdf]
  have : qf μ ⁻¹' Set.Iic t ∩ Set.Ioo 0 1 = Set.Ioo (0 : ℝ) 1 ∩ Set.Iic (cdf μ t) := by
    rw [qf_preimage_Iic]
    ext p
    simp only [Set.mem_inter_iff, Set.mem_union, Set.mem_ofPred_eq, Set.mem_compl_iff, Set.mem_Iic]
    tauto
  rw [this, volume_Ioo_inter_Iic (cdf_nonneg μ t) (cdf_le_one μ t)]

theorem integral_le_of_cdf_le {μ ν : Measure ℝ} [IsProbabilityMeasure μ] [IsProbabilityMeasure ν]
    (h : ∀ t, cdf μ t ≤ cdf ν t) {u : ℝ → ℝ} (hu : Monotone u) (hiμ : Integrable u μ)
    (hiν : Integrable u ν) : ∫ x, u x ∂ν ≤ ∫ x, u x ∂μ := by
  have hm := hu.measurable
  have rep : ∀ (κ : Measure ℝ) [IsProbabilityMeasure κ], Integrable u κ →
      (∫ x, u x ∂κ = ∫ p, u (qf κ p) ∂(volume.restrict (Set.Ioo (0 : ℝ) 1))) ∧
        Integrable (fun p => u (qf κ p)) (volume.restrict (Set.Ioo (0 : ℝ) 1)) := by
    intro κ _ hi
    have hmap := map_qf κ
    rw [← hmap] at hi
    refine ⟨?_, (integrable_map_measure hm.aestronglyMeasurable
      (measurable_qf κ).aemeasurable).mp hi⟩
    conv_lhs => rw [← hmap]
    exact integral_map (measurable_qf κ).aemeasurable hm.aestronglyMeasurable
  obtain ⟨eμ, iμ⟩ := rep μ hiμ
  obtain ⟨eν, iν⟩ := rep ν hiν
  rw [eμ, eν]
  refine integral_mono_ae iν iμ ?_
  filter_upwards [ae_restrict_mem measurableSet_Ioo] with p hp
  have := (forall_cdf_le_iff_forall_quantile_le μ ν).mp h p hp
  simp only [qf, hp, ↓reduceIte]
  exact hu this

def LawFSD (μ ν : Measure ℝ) : Prop := CurveDom (cdf μ) (cdf ν)

open Classical in
def effLaws {k : ℕ} (M : Fin k → Measure ℝ) : Finset (Fin k) :=
  Finset.univ.filter fun i => ∀ j, ¬ LawFSD (M j) (M i)

theorem exists_effLaws_dominates {k : ℕ} (M : Fin k → Measure ℝ) {i : Fin k}
    (hi : i ∉ effLaws M) : ∃ e ∈ effLaws M, LawFSD (M e) (M i) := by
  let rel : Fin k → Fin k → Prop := fun a b => LawFSD (M a) (M b)
  have : IsTrans (Fin k) rel := ⟨fun _ _ _ h1 h2 => curveDom_trans h1 h2⟩
  have : Std.Irrefl rel := ⟨fun a => curveDom_irrefl _⟩
  have wf := Finite.wellFounded_of_trans_of_irrefl rel
  simp only [effLaws, Finset.mem_filter, Finset.mem_univ, true_and, not_forall,
    not_not] at hi
  obtain ⟨j, hj⟩ := hi
  obtain ⟨e, hev, hmin⟩ := wf.has_min {a | rel a i} ⟨j, hj⟩
  refine ⟨e, ?_, hev⟩
  simp only [effLaws, Finset.mem_filter, Finset.mem_univ, true_and]
  exact fun c hce => hmin c (curveDom_trans hce hev) hce

end Laws

section KS

lemma LPM_zero_mem_Icc (t : ℝ) (x : Fin n → ℝ) : LPM 0 t x ∈ Set.Icc (0 : ℝ) 1 := by
  have h := LPMratio_mem_Icc 0 t x
  simpa [LPMratio] using h

lemma abs_ecdf_sub_le_one (t : ℝ) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    |LPM 0 t x - LPM 0 t y| ≤ 1 := by
  have hx := LPM_zero_mem_Icc t x
  have hy := LPM_zero_mem_Icc t y
  rw [abs_le]
  constructor <;> linarith [hx.1, hx.2, hy.1, hy.2]

def ksStat (x : Fin n → ℝ) (y : Fin m → ℝ) : ℝ :=
  sSup (Set.range fun t => |LPM 0 t x - LPM 0 t y|)

theorem ksStat_bddAbove (x : Fin n → ℝ) (y : Fin m → ℝ) :
    BddAbove (Set.range fun t => |LPM 0 t x - LPM 0 t y|) :=
  ⟨1, by rintro _ ⟨t, rfl⟩; exact abs_ecdf_sub_le_one t x y⟩

theorem le_ksStat (x : Fin n → ℝ) (y : Fin m → ℝ) (t : ℝ) :
    |LPM 0 t x - LPM 0 t y| ≤ ksStat x y :=
  le_csSup (ksStat_bddAbove x y) ⟨t, rfl⟩

theorem ksStat_attained (hn : 0 < n) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ∃ u, IsKnot x y u ∧ ksStat x y = |LPM 0 u x - LPM 0 u y| := by
  set K := univ.image x ∪ univ.image y
  have hK : K.Nonempty := ⟨x ⟨0, hn⟩, mem_union_left _ (mem_image_of_mem x (mem_univ _))⟩
  obtain ⟨u, hu, hmax⟩ := K.exists_max_image (fun t => |LPM 0 t x - LPM 0 t y|) hK
  have hall : ∀ t, |LPM 0 t x - LPM 0 t y| ≤ |LPM 0 u x - LPM 0 u y| := fun t => by
    rcases LPM_zero_cases x y t with ⟨hx, hy⟩ | ⟨w, hw, hx, hy⟩
    · rw [hx, hy, sub_self, abs_zero]; exact abs_nonneg _
    · rw [hx, hy]; exact hmax w (isKnot_iff_mem.mp hw)
  refine ⟨u, isKnot_iff_mem.mpr hu, le_antisymm ?_ (le_ksStat x y u)⟩
  exact csSup_le (Set.range_nonempty _) (by rintro _ ⟨t, rfl⟩; exact hall t)

end KS

section Shift

inductive Verdict
  | xDom
  | yDom
  | noDom

def nnsFSD (x : Fin n → ℝ) (y : Fin m → ℝ) (discrete : Bool) : Verdict :=
  if sdDomPairCentred x y 1 discrete = 1 then .xDom
  else if sdDomPairCentred y x 1 discrete = 1 then .yDom
  else .noDom

theorem nnsFSD_eq_yDom_iff (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    (b : Bool) :
    nnsFSD x y b = .yDom ↔ CurveDom (domCurve 1 b y) (domCurve 1 b x) := by
  have hxy := fsdUni_curveDom hn hm x y (Or.inl rfl) b
  have hyx := fsdUni_curveDom hm hn y x (Or.inl rfl) b
  unfold nnsFSD
  constructor
  · intro h
    split_ifs at h with h1 h2
    exact hyx.mp h2
  · intro h
    have h1 : sdDomPairCentred x y 1 b ≠ 1 := fun h' => curveDom_asymm h (hxy.mp h')
    simp [h1, hyx.mpr h]

theorem nnsFSD_eq_xDom_iff (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    (b : Bool) :
    nnsFSD x y b = .xDom ↔ CurveDom (domCurve 1 b x) (domCurve 1 b y) := by
  have hxy := fsdUni_curveDom hn hm x y (Or.inl rfl) b
  unfold nnsFSD
  constructor
  · intro h
    split_ifs at h with h1 h2
    exact hxy.mp h1
  · intro h
    simp [hxy.mpr h]

lemma LPM_add_const (d : ℕ) (t c : ℝ) (x : Fin n → ℝ) :
    LPM d t (fun i => x i + c) = LPM d (t - c) x := by
  have : (fun i => x i + c) = shiftSample x (-c) := funext fun i => by simp [shiftSample]
  rw [this, LPM_shiftSample, sub_eq_add_neg]

theorem LPMDominates_add_const (hn : 0 < n) (x : Fin n → ℝ) {c : ℝ} (hc : 0 < c) :
    LPMDominates 0 (fun i => x i + c) x := by
  refine ⟨fun t => ?_, ?_⟩
  · rw [LPM_add_const]
    exact LPM_monotone_target 0 x (by linarith)
  · obtain ⟨i, hi⟩ := exists_eq_sampleMin hn x
    refine ⟨sampleMin x + c / 2, ne_of_lt ?_⟩
    rw [LPM_add_const]
    have h0 : LPM 0 (sampleMin x + c / 2 - c) x = 0 :=
      LPM_eq_zero_of_lt fun j => by linarith [sampleMin_le x j]
    rw [h0]
    exact LPM_pos_of_lt (i := i) (by rw [hi]; linarith)

theorem fsdUni_shift (hn : 0 < n) (x : Fin n → ℝ) {c : ℝ} (hc : 0 < c) :
    sdDomPairCentred (fun i => x i + c) x 1 true = 1 :=
  (sdDomPairCentred_eq_one_iff hn hn _ x (Or.inl rfl) true).mpr
    (by simpa using LPMDominates_add_const hn x hc)

lemma sdDomPairCentred_fsd_zero_or_one (x : Fin n → ℝ) (y : Fin m → ℝ) :
    sdDomPairCentred x y 1 true = 0 ∨ sdDomPairCentred x y 1 true = 1 := by
  have : sdDomPairCentred x y 1 true = sdDomPair x y 1 true := by simp [sdDomPairCentred]
  rw [this]
  exact sdDomPair_eq_zero_or_one x y 1 true

theorem fsdUni_self (hn : 0 < n) (x : Fin n → ℝ) : sdDomPairCentred x x 1 true = 0 :=
  (sdDomPairCentred_fsd_zero_or_one _ _).resolve_right fun h =>
    curveDom_irrefl _ ((fsdUni_curveDom hn hn _ _ (Or.inl rfl) true).mp h)

theorem effSet_eq_of_cover {k r : ℕ} (hr : 0 < r) (X : Fin k → Fin r → ℝ) {d : ℕ}
    (hd : d = 1 ∨ d = 2 ∨ d = 3) (b : Bool) (S : Finset (Fin k))
    (hanti : ∀ i ∈ S, ∀ j ∈ S, sdDomPairCentred (X i) (X j) d b ≠ 1)
    (hcover : ∀ v ∉ S, ∃ s ∈ S, sdDomPairCentred (X s) (X v) d b = 1) :
    effSet (domMatrix X d b) univ = S := by
  have hc := fun x y => fsdUni_curveDom hr hr (X x) (X y) hd b
  ext v
  rw [effSet_eq hr X hd b]
  refine ⟨fun h => ?_, fun hv => ⟨mem_univ _, fun j _ hjv => ?_⟩⟩
  · by_contra hv
    obtain ⟨s, -, hs⟩ := hcover v hv
    exact h.2 s (mem_univ _) hs
  · by_cases hj : j ∈ S
    · exact hanti j hj v hv hjv
    · obtain ⟨s, hsS, hs⟩ := hcover j hj
      by_cases hsv : s = v
      · subst hsv
        exact curveDom_asymm ((hc _ _).mp hs) ((hc _ _).mp hjv)
      · exact hanti s hsS v hv ((hc _ _).mpr (curveDom_trans ((hc _ _).mp hs) ((hc _ _).mp hjv)))

end Shift

section Peel

variable {k r : ℕ}

def peelRemaining (D : Fin k → Fin k → ℕ) : ℕ → Finset (Fin k)
  | 0 => univ
  | t + 1 => peelRemaining D t \ effSet D (peelRemaining D t)

def peelLevel (D : Fin k → Fin k → ℕ) (t : ℕ) : Finset (Fin k) :=
  effSet D (peelRemaining D t)

def peelFamily : Fin 3 → Fin 2 → ℝ := ![![2, 2], ![1, 1], ![-5, 5]]

lemma not_fsd_of_gt {x y : Fin 2 → ℝ} (t : ℝ) (h : LPM 0 t y < LPM 0 t x) :
    sdDomPairCentred x y 1 true ≠ 1 := by
  rw [Ne, sdDomPairCentred_eq_one_iff (by norm_num) (by norm_num) x y (Or.inl rfl) true]
  simp only [↓reduceIte]
  exact fun hd => absurd (hd.1 t) (not_le.mpr h)

lemma peelFamily_dom (i j : Fin 3) :
    sdDomPairCentred (peelFamily i) (peelFamily j) 1 true = 1 ↔ i = 0 ∧ j = 1 := by
  have hAB : sdDomPairCentred (peelFamily 0) (peelFamily 1) 1 true = 1 := by
    have h := fsdUni_shift (n := 2) (by norm_num) (peelFamily 1) (c := 1) one_pos
    convert h using 2
    funext i; fin_cases i <;> norm_num [peelFamily]
  fin_cases i <;> fin_cases j <;> simp only [Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk,
    and_true, Fin.isValue, iff_true, iff_false, Fin.reduceEq, and_self, and_false]
  · exact (fsdUni_self (n := 2) (by norm_num) _).symm ▸ zero_ne_one
  · exact hAB
  · exact not_fsd_of_gt 2 (by simp only [LPM, lowerDev, Fin.sum_univ_two]; norm_num [peelFamily])
  · exact not_fsd_of_gt 1 (by simp only [LPM, lowerDev, Fin.sum_univ_two]; norm_num [peelFamily])
  · exact (fsdUni_self (n := 2) (by norm_num) _).symm ▸ zero_ne_one
  · exact not_fsd_of_gt 1 (by simp only [LPM, lowerDev, Fin.sum_univ_two]; norm_num [peelFamily])
  · exact not_fsd_of_gt (-5) (by simp only [LPM, lowerDev, Fin.sum_univ_two]; norm_num [peelFamily])
  · exact not_fsd_of_gt (-5) (by simp only [LPM, lowerDev, Fin.sum_univ_two]; norm_num [peelFamily])
  · exact (fsdUni_self (n := 2) (by norm_num) _).symm ▸ zero_ne_one

theorem peel_witness :
    (0 : Fin 3) ∈ peelLevel (domMatrix peelFamily 1 true) 0 ∧
      (2 : Fin 3) ∈ peelLevel (domMatrix peelFamily 1 true) 0 ∧
      (1 : Fin 3) ∈ peelLevel (domMatrix peelFamily 1 true) 1 ∧
      sdDomPairCentred (peelFamily 0) (peelFamily 1) 1 true = 1 ∧
      sdDomPairCentred (peelFamily 2) (peelFamily 1) 1 true ≠ 1 := by
  have hd : (1 = 1 ∨ 1 = 2 ∨ 1 = 3) := Or.inl rfl
  have hr : 0 < 2 := by norm_num
  have h0 : (0 : Fin 3) ∈ peelLevel (domMatrix peelFamily 1 true) 0 := by
    rw [peelLevel, peelRemaining, effSet_eq hr _ hd]
    exact ⟨mem_univ _, fun j _ h => absurd ((peelFamily_dom j _).mp h).2 (by decide)⟩
  have h2 : (2 : Fin 3) ∈ peelLevel (domMatrix peelFamily 1 true) 0 := by
    rw [peelLevel, peelRemaining, effSet_eq hr _ hd]
    exact ⟨mem_univ _, fun j _ h => absurd ((peelFamily_dom j _).mp h).2 (by decide)⟩
  have h1 : (1 : Fin 3) ∉ peelLevel (domMatrix peelFamily 1 true) 0 := by
    rw [peelLevel, peelRemaining, effSet_eq hr _ hd]
    exact fun h => h.2 0 (mem_univ _) ((peelFamily_dom 0 1).mpr ⟨rfl, rfl⟩)
  refine ⟨h0, h2, ?_, (peelFamily_dom 0 1).mpr ⟨rfl, rfl⟩,
    fun h => absurd ((peelFamily_dom 2 1).mp h).1 (by decide)⟩
  rw [peelLevel, effSet_eq hr _ hd]
  refine ⟨?_, fun j hj h => ?_⟩
  · rw [peelRemaining, mem_sdiff]
    exact ⟨mem_univ _, h1⟩
  · obtain ⟨rfl, -⟩ := (peelFamily_dom j 1).mp h
    rw [peelRemaining, mem_sdiff] at hj
    exact hj.2 h0

end Peel

section Omega

def c591X : Fin 3 → ℝ := ![27 / 10, 23 / 5, 14 / 5]

def c591Y : Fin 3 → ℝ := ![37 / 10, 37 / 10, 12 / 5]

lemma lowerDev_one_of_le {t v : ℝ} (h : v ≤ t) : lowerDev 1 t v = t - v := by
  simp [lowerDev, h]

lemma lowerDev_one_of_ge {t v : ℝ} (h : t ≤ v) : lowerDev 1 t v = 0 := by
  unfold lowerDev
  split_ifs with h'
  · rw [le_antisymm h h', sub_self, pow_one]
  · rfl

lemma upperDev_one_of_ge {t v : ℝ} (h : t ≤ v) : upperDev 1 t v = v - t := by
  unfold upperDev
  split_ifs with h'
  · rw [pow_one]
  · rw [le_antisymm h (not_lt.mp h'), sub_self]

lemma upperDev_one_of_le {t v : ℝ} (h : v ≤ t) : upperDev 1 t v = 0 := by
  simp [upperDev, not_lt.mpr h]

lemma c591_ratio (s : Fin 3 → ℝ) (hs : ∃ i j, s i ≠ s j) (t : ℝ) :
    LPMratio 1 t s = LPM 1 t s / (LPM 1 t s + UPM 1 t s) := by
  refine LPMratio_one_of_pos (ratio_denom_pos ?_)
  obtain ⟨i, j, hij⟩ := hs
  by_contra h
  push Not at h
  exact hij ((h i).trans (h j).symm)

lemma c591X_nonconst : ∃ i j, c591X i ≠ c591X j := ⟨0, 1, by norm_num [c591X]⟩

lemma c591Y_nonconst : ∃ i j, c591Y i ≠ c591Y j := ⟨0, 2, by norm_num [c591Y]⟩

lemma c591_cross (t : ℝ) :
    LPM 1 t c591X * UPM 1 t c591Y ≤ LPM 1 t c591Y * UPM 1 t c591X := by
  simp only [LPM, UPM, Fin.sum_univ_three, c591X, c591Y, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons,
    Nat.cast_ofNat]
  rcases le_total t (12 / 5) with h1 | h1
  · rw [lowerDev_one_of_ge (v := 12 / 5) (by linarith),
      upperDev_one_of_ge (v := 12 / 5) (by linarith),
      lowerDev_one_of_ge (v := 27 / 10) (by linarith),
      upperDev_one_of_ge (v := 27 / 10) (by linarith),
      lowerDev_one_of_ge (v := 14 / 5) (by linarith),
      upperDev_one_of_ge (v := 14 / 5) (by linarith),
      lowerDev_one_of_ge (v := 37 / 10) (by linarith),
      upperDev_one_of_ge (v := 37 / 10) (by linarith),
      lowerDev_one_of_ge (v := 23 / 5) (by linarith),
      upperDev_one_of_ge (v := 23 / 5) (by linarith)]
    nlinarith [sq_nonneg (t - 47 / 15)]
  rcases le_total t (27 / 10) with h2 | h2
  · rw [lowerDev_one_of_le (v := 12 / 5) (by linarith),
      upperDev_one_of_le (v := 12 / 5) (by linarith),
      lowerDev_one_of_ge (v := 27 / 10) (by linarith),
      upperDev_one_of_ge (v := 27 / 10) (by linarith),
      lowerDev_one_of_ge (v := 14 / 5) (by linarith),
      upperDev_one_of_ge (v := 14 / 5) (by linarith),
      lowerDev_one_of_ge (v := 37 / 10) (by linarith),
      upperDev_one_of_ge (v := 37 / 10) (by linarith),
      lowerDev_one_of_ge (v := 23 / 5) (by linarith),
      upperDev_one_of_ge (v := 23 / 5) (by linarith)]
    nlinarith [sq_nonneg (t - 47 / 15)]
  rcases le_total t (14 / 5) with h3 | h3
  · rw [lowerDev_one_of_le (v := 12 / 5) (by linarith),
      upperDev_one_of_le (v := 12 / 5) (by linarith),
      lowerDev_one_of_le (v := 27 / 10) (by linarith),
      upperDev_one_of_le (v := 27 / 10) (by linarith),
      lowerDev_one_of_ge (v := 14 / 5) (by linarith),
      upperDev_one_of_ge (v := 14 / 5) (by linarith),
      lowerDev_one_of_ge (v := 37 / 10) (by linarith),
      upperDev_one_of_ge (v := 37 / 10) (by linarith),
      lowerDev_one_of_ge (v := 23 / 5) (by linarith),
      upperDev_one_of_ge (v := 23 / 5) (by linarith)]
    nlinarith [sq_nonneg (t - 47 / 15)]
  rcases le_total t (37 / 10) with h4 | h4
  · rw [lowerDev_one_of_le (v := 12 / 5) (by linarith),
      upperDev_one_of_le (v := 12 / 5) (by linarith),
      lowerDev_one_of_le (v := 27 / 10) (by linarith),
      upperDev_one_of_le (v := 27 / 10) (by linarith),
      lowerDev_one_of_le (v := 14 / 5) (by linarith),
      upperDev_one_of_le (v := 14 / 5) (by linarith),
      lowerDev_one_of_ge (v := 37 / 10) (by linarith),
      upperDev_one_of_ge (v := 37 / 10) (by linarith),
      lowerDev_one_of_ge (v := 23 / 5) (by linarith),
      upperDev_one_of_ge (v := 23 / 5) (by linarith)]
    nlinarith [sq_nonneg (t - 47 / 15)]
  rcases le_total t (23 / 5) with h5 | h5
  · rw [lowerDev_one_of_le (v := 12 / 5) (by linarith),
      upperDev_one_of_le (v := 12 / 5) (by linarith),
      lowerDev_one_of_le (v := 27 / 10) (by linarith),
      upperDev_one_of_le (v := 27 / 10) (by linarith),
      lowerDev_one_of_le (v := 14 / 5) (by linarith),
      upperDev_one_of_le (v := 14 / 5) (by linarith),
      lowerDev_one_of_le (v := 37 / 10) (by linarith),
      upperDev_one_of_le (v := 37 / 10) (by linarith),
      lowerDev_one_of_ge (v := 23 / 5) (by linarith),
      upperDev_one_of_ge (v := 23 / 5) (by linarith)]
    nlinarith [sq_nonneg (t - 47 / 15)]
  rw [lowerDev_one_of_le (v := 12 / 5) (by linarith),
    upperDev_one_of_le (v := 12 / 5) (by linarith),
    lowerDev_one_of_le (v := 27 / 10) (by linarith),
    upperDev_one_of_le (v := 27 / 10) (by linarith),
    lowerDev_one_of_le (v := 14 / 5) (by linarith), upperDev_one_of_le (v := 14 / 5) (by linarith),
    lowerDev_one_of_le (v := 37 / 10) (by linarith),
    upperDev_one_of_le (v := 37 / 10) (by linarith),
    lowerDev_one_of_le (v := 23 / 5) (by linarith), upperDev_one_of_le (v := 23 / 5) (by linarith)]
  nlinarith [sq_nonneg (t - 47 / 15)]

theorem c591_omega : sdDomPairCentred c591X c591Y 1 false = 1 := by
  rw [fsdUni_curveDom (by norm_num) (by norm_num) _ _ (Or.inl rfl)]
  simp only [domCurve, ↓reduceIte, Bool.false_eq_true]
  refine ⟨fun t => ?_, 3, ?_⟩
  · change LPMratio 1 t c591X ≤ LPMratio 1 t c591Y
    rw [c591_ratio _ c591X_nonconst, c591_ratio _ c591Y_nonconst,
      div_le_div_iff₀ (ratio_denom_pos ?_) (ratio_denom_pos ?_)]
    · nlinarith [c591_cross t]
    · obtain ⟨i, j, hij⟩ := c591X_nonconst
      by_contra h; push Not at h; exact hij ((h i).trans (h j).symm)
    · obtain ⟨i, j, hij⟩ := c591Y_nonconst
      by_contra h; push Not at h; exact hij ((h i).trans (h j).symm)
  · change LPMratio 1 3 c591X ≠ LPMratio 1 3 c591Y
    rw [c591_ratio _ c591X_nonconst, c591_ratio _ c591Y_nonconst]
    simp only [LPM, UPM, lowerDev, upperDev, Fin.sum_univ_three]
    norm_num [c591X, c591Y]

theorem c591_not_fsd : sdDomPairCentred c591X c591Y 1 true = 0 := by
  refine (sdDomPairCentred_fsd_zero_or_one _ _).resolve_right fun h => ?_
  rw [sdDomPairCentred_eq_one_iff (by norm_num) (by norm_num) _ _ (Or.inl rfl)] at h
  have := h.1 3
  simp only [LPM, lowerDev, Fin.sum_univ_three] at this
  norm_num [c591X, c591Y] at this

theorem c591_utility :
    Monotone (fun v : ℝ => if 3 < v then (1 : ℝ) else 0) ∧
      mean (fun i => if 3 < c591X i then (1 : ℝ) else 0) <
        mean (fun j => if 3 < c591Y j then (1 : ℝ) else 0) := by
  refine ⟨fun a b hab => ?_, ?_⟩
  · simp only
    split_ifs with h1 h2
    · exact le_rfl
    · exact absurd (h1.trans_le hab) h2
    · norm_num
    · exact le_rfl
  · simp only [mean, Fin.sum_univ_three]
    norm_num [c591X, c591Y]

end Omega

section Cluster

variable {k r : ℕ}

def clusterLoop (D : Fin k → Fin k → ℕ) (minC : ℕ) : ℕ → Finset (Fin k) → List (Finset (Fin k))
  | 0, R => if R.Nonempty then [R] else []
  | fuel + 1, R =>
    if minC < R.card then
      if effSet D R = ∅ then [R]
      else if (R \ effSet D R).card ≤ minC then [effSet D R, R \ effSet D R]
      else effSet D R :: clusterLoop D minC fuel (R \ effSet D R)
    else if R.Nonempty then [R] else []

def mergeLast (minC : ℕ) (l : List (Finset (Fin k))) : List (Finset (Fin k)) :=
  match l.reverse with
  | last :: prev :: rest => if last.card < minC then ((prev ∪ last) :: rest).reverse else l
  | _ => l

def sdCluster (D : Fin k → Fin k → ℕ) (minC : ℕ) : List (Finset (Fin k)) :=
  mergeLast minC (clusterLoop D minC k univ)

lemma clusterLoop_dropLast (D : Fin k → Fin k → ℕ) (minC fuel : ℕ) (R : Finset (Fin k)) :
    ∀ S ∈ (clusterLoop D minC fuel R).dropLast, ∃ R', S = effSet D R' := by
  induction fuel generalizing R with
  | zero =>
    intro S hS
    simp only [clusterLoop] at hS
    split_ifs at hS <;> simp at hS
  | succ fuel ih =>
    intro S hS
    simp only [clusterLoop] at hS
    split_ifs at hS
    · simp at hS
    · simp only [List.dropLast_cons_cons, List.dropLast_singleton, List.mem_cons,
        List.not_mem_nil, or_false] at hS
      exact ⟨R, hS⟩
    · rcases List.eq_nil_or_concat (clusterLoop D minC fuel (R \ effSet D R)) with h | ⟨L, b, h⟩
      · rw [h] at hS; simp at hS
      · rw [h, List.concat_eq_append, ← List.cons_append, List.dropLast_concat] at hS
        rcases List.mem_cons.mp hS with rfl | hS
        · exact ⟨R, rfl⟩
        · refine ih (R \ effSet D R) S ?_
          rw [h, List.concat_eq_append, List.dropLast_concat]
          exact hS
    · simp at hS
    · simp at hS

lemma mergeLast_dropLast (minC : ℕ) (l : List (Finset (Fin k))) :
    ∀ S ∈ (mergeLast minC l).dropLast, S ∈ l.dropLast.dropLast ∨ S ∈ l.dropLast := by
  intro S hS
  unfold mergeLast at hS
  rcases h : l.reverse with _ | ⟨last, _ | ⟨prev, rest⟩⟩ <;> simp only [h] at hS
  · exact Or.inr hS
  · exact Or.inr hS
  · split_ifs at hS
    · left
      have hl : l = (rest.reverse ++ [prev]) ++ [last] := by
        rw [← List.reverse_reverse l, h]; simp
      rw [hl, List.dropLast_concat, List.dropLast_concat]
      simpa [List.reverse_cons, List.dropLast_concat] using hS
    · exact Or.inr hS

theorem sdCluster_antichain_of_acyclic {D : Fin k → Fin k → ℕ} (hD : Acyclic D) (minC : ℕ)
    {S : Finset (Fin k)} (hS : S ∈ (sdCluster D minC).dropLast) {i j : Fin k} (hi : i ∈ S)
    (hj : j ∈ S) : D i j ≠ 1 := by
  have hS' : S ∈ (clusterLoop D minC k univ).dropLast := by
    rcases mergeLast_dropLast minC _ S hS with h | h
    · exact List.mem_of_mem_dropLast h
    · exact h
  obtain ⟨R', rfl⟩ := clusterLoop_dropLast _ minC k univ S hS'
  exact effSet_antichain_of_acyclic hD R' hi hj

def chainFamily : Fin 3 → Fin 1 → ℝ := ![![2], ![1], ![0]]

lemma chainFamily_dom (i j : Fin 3) :
    sdDomPairCentred (chainFamily i) (chainFamily j) 1 true = 1 ↔ i < j := by
  have hs : ∀ i j : Fin 3, i < j → sdDomPairCentred (chainFamily i) (chainFamily j) 1 true = 1 := by
    intro i j hij
    have h := fsdUni_shift (n := 1) one_pos (chainFamily j) (c := (j : ℝ) - i)
      (by have : (i : ℝ) < j := by exact_mod_cast hij
          linarith)
    convert h using 2
    funext l
    fin_cases i <;> fin_cases j <;> fin_cases l <;> simp_all [chainFamily] <;> norm_num
  have hc := fun x y => fsdUni_curveDom (n := 1) (m := 1) one_pos one_pos (chainFamily x)
    (chainFamily y) (d := 1) (Or.inl rfl) true
  constructor
  · intro h
    by_contra hij
    rcases lt_or_eq_of_le (not_lt.mp hij) with hji | rfl
    · exact curveDom_asymm ((hc _ _).mp h) ((hc _ _).mp (hs j i hji))
    · exact curveDom_irrefl _ ((hc _ _).mp h)
  · exact hs i j

theorem cluster_witness :
    sdCluster (domMatrix chainFamily 1 true) 2 = [{0}, {1, 2}] ∧
      sdDomPairCentred (chainFamily 1) (chainFamily 2) 1 true = 1 := by
  have he : effSet (domMatrix chainFamily 1 true) univ = {0} := by
    ext v
    rw [effSet_eq one_pos _ (Or.inl rfl), mem_singleton]
    constructor
    · rintro ⟨-, h⟩
      by_contra hv
      exact h 0 (mem_univ _) ((chainFamily_dom 0 v).mpr (Fin.pos_iff_ne_zero.mpr hv))
    · rintro rfl
      exact ⟨mem_univ _, fun j _ hj => absurd ((chainFamily_dom j 0).mp hj) (Fin.not_lt_zero j)⟩
  refine ⟨?_, (chainFamily_dom 1 2).mpr (by decide)⟩
  have h1 : (univ : Finset (Fin 3)) \ {0} = {1, 2} := by decide
  simp only [sdCluster, clusterLoop, he, h1, card_univ, Fintype.card_fin]
  decide

end Cluster

section Vignette

def ofChunks76 : List ℕ → ℕ
  | [] => 0
  | c :: cs => c + 2 ^ 304 * ofChunks76 cs

def digits68 : ℕ → ℕ → List ℕ
  | 0, _ => []
  | k + 1, N => N % 2 ^ 68 :: digits68 k (N / 2 ^ 68)

lemma length_digits68 (k N : ℕ) : (digits68 k N).length = k := by
  induction k generalizing N with
  | zero => rfl
  | succ k ih => simp [digits68, ih]

def dval (d : ℕ) : ℝ := ((d : ℝ) - 2 ^ 67) / 2 ^ 64

lemma dval_le_iff {a b : ℕ} : dval a ≤ dval b ↔ a ≤ b := by
  unfold dval
  rw [div_le_div_iff_of_pos_right (by positivity), sub_le_sub_iff_right, Nat.cast_le]

lemma dval_lt_iff {a b : ℕ} : dval a < dval b ↔ a < b := by
  unfold dval
  rw [div_lt_div_iff_of_pos_right (by positivity), sub_lt_sub_iff_right, Nat.cast_lt]

def allLe : List ℕ → List ℕ → Bool
  | a :: l, b :: m => Nat.ble a b && allLe l m
  | [], [] => true
  | _, _ => false

def allLt : List ℕ → List ℕ → Bool
  | a :: l, b :: m => Nat.blt a b && allLt l m
  | [], [] => true
  | _, _ => false

def cntLe (c : ℕ) : List ℕ → ℕ
  | [] => 0
  | a :: l => (if Nat.ble a c then 1 else 0) + cntLe c l

lemma allLe_getD {l m : List ℕ} (h : allLe l m = true) {i : ℕ} (hi : i < l.length) :
    l.getD i 0 ≤ m.getD i 0 := by
  induction l generalizing m i with
  | nil => simp at hi
  | cons a l ih =>
    cases m with
    | nil => simp [allLe] at h
    | cons b m =>
      simp only [allLe, Bool.and_eq_true, Nat.ble_eq] at h
      cases i with
      | zero => simpa using h.1
      | succ i => simpa using ih h.2 (by simpa using hi)

lemma allLt_getD {l m : List ℕ} (h : allLt l m = true) {i : ℕ} (hi : i < l.length) :
    l.getD i 0 < m.getD i 0 := by
  induction l generalizing m i with
  | nil => simp at hi
  | cons a l ih =>
    cases m with
    | nil => simp [allLt] at h
    | cons b m =>
      simp only [allLt, Bool.and_eq_true, Nat.blt_eq] at h
      cases i with
      | zero => simpa using h.1
      | succ i => simpa using ih h.2 (by simpa using hi)

lemma card_filter_getD_le (c : ℕ) : ∀ {k : ℕ} (l : List ℕ), l.length = k →
    #{i : Fin k | l.getD i 0 ≤ c} = cntLe c l
  | 0, [], _ => rfl
  | k + 1, a :: l, hl => by
    have hl' : l.length = k := by simpa using hl
    rw [card_filter, Fin.sum_univ_succ]
    simp only [Fin.val_zero, List.getD_cons_zero, Fin.val_succ, List.getD_cons_succ, cntLe,
      Nat.ble_eq]
    rw [← card_filter_getD_le c l hl', card_filter]
  | 0, _ :: _, hl => by simp at hl
  | _ + 1, [], hl => by simp at hl

def colOf {k : ℕ} (l : List ℕ) : Fin k → ℝ := fun i => dval (l.getD i 0)

lemma LPM_zero_colOf {k : ℕ} (l : List ℕ) (hl : l.length = k) (c : ℕ) :
    LPM 0 (dval c) (colOf (k := k) l) = (cntLe c l : ℝ) / k := by
  rw [LPM_zero_eq_ecdf]
  simp only [colOf, dval_le_iff]
  rw [card_filter_getD_le c l hl]

def col (N : ℕ) : Fin 1000 → ℝ := colOf (digits68 1000 N)

def cnt (N c : ℕ) : ℕ := cntLe c (digits68 1000 N)

lemma LPM_zero_col (N c : ℕ) : LPM 0 (dval c) (col N) = (cnt N c : ℝ) / 1000 := by
  rw [col, LPM_zero_colOf _ (length_digits68 1000 N), cnt, Nat.cast_ofNat]

lemma col_le {N M : ℕ} (h : allLe (digits68 1000 N) (digits68 1000 M) = true) (i : Fin 1000) :
    col N i ≤ col M i :=
  dval_le_iff.mpr (allLe_getD h (by rw [length_digits68]; exact i.2))

lemma col_lt {N M : ℕ} (h : allLt (digits68 1000 N) (digits68 1000 M) = true) (i : Fin 1000) :
    col N i < col M i :=
  dval_lt_iff.mpr (allLt_getD h (by rw [length_digits68]; exact i.2))

lemma not_fsd_col {N M c : ℕ} (h : cnt M c < cnt N c) :
    sdDomPairCentred (col N) (col M) 1 true ≠ 1 := by
  rw [Ne, sdDomPairCentred_eq_one_iff (by norm_num) (by norm_num) _ _ (Or.inl rfl)]
  simp only [↓reduceIte]
  intro hd
  have := hd.1 (dval c)
  rw [LPM_zero_col, LPM_zero_col] at this
  have : (cnt N c : ℝ) ≤ cnt M c := by linarith
  exact absurd (by exact_mod_cast this) (not_le.mpr h)

lemma LPM_zero_le_of_forall_le {x y : Fin n → ℝ} (h : ∀ i, x i ≤ y i) (t : ℝ) :
    LPM 0 t y ≤ LPM 0 t x := by
  rw [LPM_zero_eq_ecdf, LPM_zero_eq_ecdf]
  refine div_le_div_of_nonneg_right ?_ (Nat.cast_nonneg n)
  exact_mod_cast card_le_card fun i hi => by
    simp only [mem_filter, mem_univ, true_and] at hi ⊢
    exact (h i).trans hi

theorem fsdUni_of_forall_lt (hn : 0 < n) {x y : Fin n → ℝ} (h : ∀ i, x i < y i) :
    sdDomPairCentred y x 1 true = 1 := by
  rw [sdDomPairCentred_eq_one_iff hn hn _ _ (Or.inl rfl)]
  simp only [↓reduceIte]
  rw [fsd_iff_lt]
  refine ⟨LPM_zero_le_of_forall_le fun i => (h i).le, ?_⟩
  obtain ⟨i₀, hi₀⟩ := exists_eq_sampleMin hn x
  refine ⟨x i₀, ?_⟩
  rw [LPM_eq_zero_of_lt fun j => by rw [hi₀]; exact (sampleMin_le x j).trans_lt (h j),
    LPM_zero_eq_ecdf]
  refine div_pos ?_ (Nat.cast_pos.mpr hn)
  exact_mod_cast card_pos.mpr ⟨i₀, by simp⟩

lemma LPM_comp_perm (d : ℕ) (t : ℝ) (x : Fin n → ℝ) (σ : Equiv.Perm (Fin n)) :
    LPM d t (x ∘ σ) = LPM d t x := by
  unfold LPM
  rw [show (∑ i, lowerDev d t ((x ∘ σ) i)) = ∑ i, lowerDev d t (x (σ i)) from rfl,
    Equiv.sum_comp σ fun i => lowerDev d t (x i)]

lemma fsdUni_perm (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    (σ : Equiv.Perm (Fin n)) (τ : Equiv.Perm (Fin m)) :
    sdDomPairCentred (x ∘ σ) (y ∘ τ) 1 true = 1 ↔ sdDomPairCentred x y 1 true = 1 := by
  rw [sdDomPairCentred_eq_one_iff hn hm _ _ (Or.inl rfl),
    sdDomPairCentred_eq_one_iff hn hm _ _ (Or.inl rfl)]
  simp only [↓reduceIte, LPMDominates, LPM_comp_perm]

theorem nnsFSD_perm (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    (σ : Equiv.Perm (Fin n)) (τ : Equiv.Perm (Fin m)) :
    nnsFSD (x ∘ σ) (y ∘ τ) true = nnsFSD x y true := by
  unfold nnsFSD
  have h1 := fsdUni_perm hn hm x y σ τ
  have h2 := fsdUni_perm hm hn y x τ σ
  by_cases a : sdDomPairCentred x y 1 true = 1
  · simp [a, h1.mpr a]
  · by_cases b : sdDomPairCentred y x 1 true = 1
    · simp [a, b, mt h1.mp a, h2.mpr b]
    · simp [a, b, mt h1.mp a, mt h2.mp b]

def c519xN : ℕ := ofChunks76 [
    0x81b04000565f702aa3aba000055b5a971c94308000556cdc37368bae000530b29b4cbc4ea000,
    0xc66734d1958d60005c03e300ffd3ae0005b0da4efe2e7ea000588bae55e8903a00057df91952,
    0xd82c0005ebe8209528c4a0005e387d919cdfe00005d7a47d62c66ea0005ca13f2ce9f7a20005,
    0x65a70a6bebc60005f2a1867621db20005f29a0b1dc513e0005f25e63dce17040005ede34201d,
    0x93f000609cd87e27912b0006088bc72c85f01000601db3c2bedda50005fc5cb3682bae80005f,
    0xbe27e81d90e00060ef8a6d65290700060e9a68921cb0e00060e6ce4fff8fc700060ab310cffd,
    0x5d000624f9418f09923000623c8202d9d4f900061d19cb5d9f28300061329ae474f54700060f,
    0x438798d48c0006385368908e6ac000635b020fb4ef64000626ca7cfa60dc9000625bae75fdcc,
    0xc000645f4ff1ecb8040006458dd5d52b14500063fba09ff1b22a00063f2bd6205bc9a00063e5,
    0xc1daba1fc00065431045fabbf700065067c77bf7c4700065034ddfffab7f00064ed26f4a5afd,
    0x00065975c81f769ba00065954dc0cec1d3000655205a1cf9a4b00065502aae74bf0000065501,
    0xf675e058000661d994195c611000661d2705722568000660cde5aa5062a00065d163f8fd585c,
    0x0067384efac3e2f500066f49b1ba2b30500066d87f5dbdf6a000066601b9af22da3000662ce4,
    0x040aa4c0006787f1ae326d5d00067793bbfd9caaf000675a59107e867b000673c4ea9b8f8b20,
    0x0680554faa3c3c100067e039dcee011d00067da74305f4a6800067c3ebd9ab2c5900067be2fe,
    0x66b08d00068a88d843331bd00068a80e4f31cd6e00068a137e082698d000689ca62cfb27fc00,
    0x68fbe074d0370c00068ed6d1423108300068e5d048eb543500068c96c71d50eb100068b250cb,
    0xa9e33000696857a476c910000695b01587f84e400069455d1af818c1000694323c9de61dd000,
    0xa2ce7108c97ae0006a2a266a56781f00069ca8d14f4e21800069c4aceca23d8700069acf4cda,
    0xb66d0006a61e48fdca8bb0006a56b1cebde9b30006a454ee469eeac0006a3e17ac9514720006,
    0xd53136c827150006ac1837973202e0006a946169b999f10006a862c129a07320006a62b3dd76,
    0x3d80006b06f4f8acf48f0006aef92ca00fe970006aed807df71cfd0006ae89161ea76540006a,
    0x43e308d8a270006bc24f269ed7420006bc0efc15895db0006b71cded3fcf0a0006b6852bb479,
    0xb10006c24ae86953a630006bfacae17664410006bdc93419538800006bd667743e18b50006bd,
    0xb5e00f20260006c77f0b44380550006c6ee4f56338010006c56129702b9a00006c3839ad978e,
    0xb0006cd7995e3eaa240006cd2715c78b9c80006caa37f749b9610006c93aef69294860006c84,
    0x55a1824040006d090e99a368fe0006d084608e986190006d054c3e5c71140006cfc2c97d960b,
    0x0006d80569e4f3d010006d4d3429fac6f90006d49418239e6d30006d2c44d81015800006d185,
    0x2f2070c50006dcd9ab7f9e44b0006dca30ec52a00b0006dc5d58cbbb6180006d815dac7dc132,
    0x006e23288b6071c20006e14e4079cc35a0006e07bf4b548e370006e019f94afbfda0006dd4a1,
    0xf36c02e0006e76cc9e4605110006e6c95be090ad20006e683d45919c650006e6052e6cac2830,
    0x06ed00cbaeeb8ac0006ec2e28fad95d20006e9944f71477e50006e8f64aeb168270006e83f18,
    0xc3f1b50006efc9dbc59c9320006efbe99e8b982f0006eea31b08e49a60006ed9f171625cd500,
    0x6f553e1020851f0006f36922bca77160006f28e7d5ea268f0006f1ea67314c8110006f10035e,
    0x406250006f63d0b4f94bad0006f606f362e9a7d0006f5d78226d9da40006f554e623e928b000,
    0xfc021ff0d66450006fb3ea4c90b5cb0006f9d2b21175b9c0006f957c5c0f84a40006f93c7ad9,
    0xefd18006fe9f9326ecba30006fe8e9f41baac00006fddb07d5661af0006fc63037a024550006,
    0x26e1788b52a4800701eb0d280521b8007019cf8192655b000700d2fa2eefa34000700cf29134,
    0x3b080070c62b9bf5493e000709c3c37cf87018007072b6329cf3f3000706a74ed7dda6a80070,
    0x2089146538c80071417864d22c5a80070fbbed71fa68480070df9afd4f929c00070d724dbe49,
    0xba80071b56bfaa1308780071ac9778b7320b80071aa329da56a39000718855d3468556800716,
    0x8eadb774398007226dbfc0a2ec58007224b59b1881bc8007207e6b3e1ef9580071fef41cc830,
    0xb00072a72a3daa3140000727c96163b1cc7800726a3edf0cc70e80072679c6a6844640007257,
    0xa84f542fa00072bfcc56141d1080072ba75af14efa380072b4fef8cf374f00072a893a239f64,
    0x0007361cdeda466648007360681ac772168007356c74219976d000733530db4b171800072e6b,
    0x800fc385000737c19eae319e50007376b825e4a11100073710802d3297c8007361e1fd9e681c,
    0x0073c42c7c09a7eb80073b03093d6bd23000739d1a34cfb063000737ed398b7817f000737dfe,
    0x3bdf71f000742efd92cf41e580074237cc2f29f5000073f4fc69971df600073d9815c18ead28,
    0x0745ca316aeb833000745936ead1c79b800745676267da72e00074463b16534691800743e432,
    0xbf13ff800748129904273078007477b61178370400074744ee6fadcea800746e5260a41fa300,
    0x74b9f90ef9de9b80074a71d19e690d900074a22cba20c3a1800748e5e3552898c0007486463f,
    0xca6ea00074e83982e5e7b800074e27aef53bc2800074d344d99d4e4300074c8b06fafe8d3000,
    0x524d9e2c31d690007502cd97260b0d0007502a6957f3f6180074fdeaad4e80fc00074e9159f3,
    0x1df38007578ca49c6c3e88007565f5e9cc5fa0000755d1ab981e1568007531397b095bd38007,
    0xad9cfcd1812f80075a5b37ae77e6880075919cfad966c5800758a1e2a9dee94800757dfb215b,
    0x82880075ed3ed71e824a80075d8123410603780075bfab0e7fa4ad00075b3f2aeff245580075,
    0x8aa1676d7fd0007623142b6c6a800007601c5896aa39c80075ffd6d31df25000075f418c2273,
    0x7b000765ea59b6cc4910007657ef9e4a29088007653e9ef5b73510007643598c562c23000763,
    0x1531be5715000767c729aa9fe03800767156008ae19d00076655634017b6780076603a808c94,
    0xc80076c2b0e9cf674300076b3f2b0bb1caf80076917c549f748f80076879c25a3c3b58007683,
    0x12f5a4d0f00076cf4dce24643700076cf47a875e15000076cb60fbbbe78380076cb16df132e9,
    0x00077084ab03da92180077028dee22884c00076d9b3c0f1f73080076d2b2772c3ec080076d10,
    0x970d246d800771d94728eecb0800771b4653875b028007715a53803c64d80077108f7c13d587,
    0x007791b6e8119fd1000779072cfbd976180077816829a0ba29800776109a571c6b0800773248,
    0x1773a3580077b6580e286cf680077a980e33bb90e00077a35cac93333a00077a226a1d0f69c8,
    0x077f6fe7a09fd1f80077f67baece1bad80077e60eccb6debe00077be33da37196100077bb5a7,
    0x87e5974007826ad3ce7f0b70007824bc810786e7800781ceda5c3ff4ac007802e65fb175d700,
    0x78614aefacdf5d4007846f4ab3d8fecc0078426f3f2247cd000783d7da55cc96980078297d99,
    0x4a1eec007893821cdb3237800788894655af598800788034a8e318b7000786f7246dfa8d8000,
    0x8bfb3a1d8dd3480078ad034a8551e980078aa891ce62b9700078aa6f19fca0f8000789f96434,
    0x143d00078de918dd59893c0078dd733b9fba0980078dae6d1fa6fc400078cb69ba21e22d8007,
    0x07989db05d03c007902be91f3ee5c800790165015db05fc0078f0fed6b3192ac0078e22c4ad3,
    0x1c8c00795dd8362fa4c6c0079548d233c2f23c00793d73f62700e2c007929dd31b31ca440079,
    0x5c3533ad968c00798dc8a18560c1c00796194bc6ffa09400796030e703545a000795edd13243,
    0xf9c0079e799e7f8a9a400079c78ec3e77f2000079bfc11b9b506380079aacfc4baf73800079a,
    0x7da8ac56908007a01b780a2edc3c0079fd87d7de507000079ea9798a4acc4c0079e9973d5934,
    0x1c007a3d8abb10f9288007a2e75afa062dac007a1b3a5486f2c8c007a11c6c95846c90007a0a,
    0x8194f34304007a62a3fd09f54b4007a5e0488e4e12f4007a5c5a08bf202d8007a53cb585c084,
    0xc007a7f50a46e332b4007a70772bd9eb234007a6ca82742dba28007a67d4feab4ad58007a676,
    0x254acc730007ac8fbf6368f2c4007ac81d577e63e70007ab616c2b1e5e9c007aab2ebc05aa75,
    0x007aeff5ea2289160007adbda7b77e913c007ad98728d851f34007ad607e930ee488007ace16,
    0xd20bb6a8007b342c3d8564270007b2c0bb0b044c74007b1ccbddeaa9f4c007b1ac6e4d70a4a8,
    0x07b6d2d595f6b7a0007b670d555e0b3bc007b4b62b2bf5450c007b47631fe90e3c0007b47089,
    0xd578224007b88eb9e93eaaa4007b83802e57e2a3c007b7f0c9c17facb4007b7191dbfeeadcc0,
    0x7bb77978a7bae84007ba4cacc3e91ee8007ba0c697d8c8280007b95ba84e1f8328007b8b73ca,
    0x928e70007be70bdf9bf6184007bce0a015ec4fd8007bc1f741ccb087c007bbfe414d4d007400,
    0xc3832ea2af5d04007c2cefdb17dcb38007c2b9db0ee53bba007c0d8cd11d0b768007c0350a2c,
    0x52fcc007c61650f4118b44007c5b67f43740ee2007c513168a64d882007c3a9242a020c04007,
    0x7ec507a4de784007c78e3153c7ee8e007c7572b52661512007c668f029ffb9e6007c633df683,
    0x8f88007c93497e2447a4e007c9163b69837e02007c8dcd2b37ba916007c832cbc631b926007c,
    0x4afea6a640de007ca0d5321c78b46007c9b080b2e10802007c982c67434096c007c964b3fbcf,
    0x640007cc6d3d273ea132007cc39f17d62e1ba007cb71820347ef94007cac5eeeaa1ea72007ca,
    0xa33c2a34c54007cd85e16e36a280007cd800caeaae61c007ccd4052d805410007cc999d88631,
    0x4c007d05a64d3150a52007cf3a4a75e321aa007cee7d89aad1888007ce4c0717876eb6007cdd,
    0xacb9d723d2007d1f1aa4d4ddfb2007d12bcd3f54bee0007d0d1224bf58eaa007d09697334644,
    0x6007d77a0ae394e758007d66b83fce5e9f0007d632f56c254666007d450f644436410007d229,
    0x21e8d8f84007dde9d182c6bcf8007dc719d9ca09ea2007dbd644fc4bce3e007d83b498e6dd5a,
    0x007e12ddecc9b24d0007e0962a4a2a8c2a007df73305bf3aa60007de422c125dde26007de3d7,
    0x8557094f007e3ea063e83af97007e2e87eed47bfac007e239ff8c5c7ed5007e16b8dd1197032,
    0x07e5eafa9cb16ff2007e5cc9e6f91c4e3007e4f2b934543fa1007e4eb542b827dff007e447c1,
    0xd9939e4007e712cec46555e5007e70ff978b5a55d007e70d6f53d79114007e6268289cba07f0,
    0x7e83504752c79f6007e8076b487323f0007e7a46aaa4df6b9007e78e2011dd6b1c007e73f94d,
    0x26d5f3007e953cbfbcc02ff007e8f762af5dd775007e8e0d0a9ba8a03007e897c368c7591b00,
    0xebb8c80b462b9f007eafc5a48da77d4007eac5c23777816a007eaa8532a0164b9007e9cd1a3a,
    0x28731007ecbb8d13819f46007eca3d5af4b119e007ec2668017755b5007ec1d9bef9a1b29007,
    0x1a65148623670007f06df8eddb2f08007f0268ddbd56464007edbec0d628a69d007ed2b6ec76,
    0x88c2007f2de2448bcd610007f2b15054f823f5007f22b3679b19f3a007f1c6b0e72badfd807f,
    0x475e535429d6807f6a57f218b4064807f5067113f81fa1807f4791022884b57007f39106414a,
    0x306007f89dd70968338a007f865f9d573dcb0007f78af13852fcb3807f777a9a5f5df6e007f7,
    0x74224fd0953807f956d7be8f229b807f8ffca099bb0ef807f8bb1a087c1121807f8b128eb5ed,
    0x41f08006e579af9f0fbf07fef5e18a4a5ffd707fd2c19ba3628a5607fc5659705b9552407fbe,
    0x5cc7250ff780801ddc36b7ef7aac080179c2fded9690f0800bb2c87221ec148800b06694f48d,
    0x900806511984448393808055f2e2924ed08c08053c41f44b87f340804ee2d02a9206000802d9,
    0xda37a6b1000809ac80229ba46200808d856b761bd0780808121696b469d980806ae811e31c6c,
    0x8080e03368104b1048080d91af5c5b76fc0080c12501a927b128080bf6ba7c51080d0080a8e3,
    0x99a2470500810724406b469a10080f4bc59b0141a28080e6418aba61b890080e571cd571d8f0,
    0x08120cd683431a5900811dc332014fc5b00811134eeb4285b900810b70b53a5da3200810aded,
    0xf02fe24008145a3e25e6f85d00813f53e033e98720081315c933d6e8ed008127985317adb140,
    0x8175d4c9b460195008172ac9cbdc50bc00816d648d9aab1e000815b157dabd11ed00814729c4,
    0x8cd73700819d5e245fa2693008196d88803294bc00818b9d966a4150e0081836a0483eebb100,
    0x1c55b3d5377cb90081b70d49c944d610081b0d9855d5e0960081a88052ec884c400819e052c1,
    0x906ea0081e6cea9be719780081e3eb7fb25bc490081e1e165f07dddc0081ddbf0e67b6478008,
    0x0563c710c0d5c0081fb4e966c88e930081f0eb89fba07520081ee776434869f90081e86db16a,
    0x145e008239beacb7417ac00822505b4e41f81a008211900424a48ce00820c425a1bb3b180082,
    0x7b98a1f394a400824fc5c1acc8af400824e40b0c79e90600824c47e07d27de6008247827e2ad,
    0xe86008274375ed3181d200826d488964eed52008266e45ddec9f46008258eafed72d20c00825,
    0x27a7b936ff60082a3302d8c10a9000828484f1835a3ba008281dafddf6d94200827f9d207d4c,
    0x760082cbfdea871f13a0082cbade80c8f0b00082baba9189801560082b06553f5fd1700082a5,
    0x3ce721648600830521faf9706c00082e8d89a12bed6e0082e69e7a42c70f20082d9ddf47b936,
    0x400834a209bde5363e008349fa8cb9e6eec008323f472a3a203c00830b2ef189b5b1a008307e,
    0x3476efea200836c5f33d4820a600836445984b9d4980083516ee78c6ca2200834cfafbff5ee0,
    0x0083880890093ae340083755c5e1fbcb0e0083747f2579a303e00836ec043c183e1e00836e5e,
    0x2cf174600083b09b46177dad20083a2fcdaa3af39c00839f06e71cbb66c0083894bb5d3fb1a2,
    0x083d1d85e2d22b680083cc83af76de8300083c424b5de189400083bcf699bc324640083b2a59,
    0x09619e80083ff48349d123560083e624ca4930c0c0083e58290694da4a0083db52faced34300,
    0x84490d5e792bfb400842a70cf7d9f06c00841c32179b3d51400840d97b6d142b9800840103dc,
    0xf946f800846bc9d5de66ad400846b9ede9120e8000846158f8dbe1e5c0084564a18746c0c000,
    0x4a1273451036fc00849531c5b14e29800848ac49bacebe940084795b10c31dd7400846e11d0b,
    0x623b00084ccf537b32c4f80084cb22bf60b68480084c90f5947d9db40084c58a4b86ea5e0008,
    0xf6b400402e1680084f5d45f28fa2340084db40d943b08980084d22e437f808a40084d1863019,
    0x409400850661541d14c34008502fc2c3e1c4b00084f7baa98e556040084f79a80ef789140084,
    0xefa9ac6baa8000854c82884f9a3bc00853059c081fb5000085205e2d58d742c0085182062bd1,
    0x14c008589a3f2f028d34008557aa9ca1284e40085557bbff48f60c008550b3a6c3cd43c00854,
    0x69ec22eea200085c1cc24bdc277c0085bddf7faaae2a40085b3562a24abf480085b07ca2ace8,
    0xd40085f038d7d2514f00085e74759280c7480085e65e127d15c340085d71de6fd4e66c0085d3,
    0x9a6609662c008609c7f8a0aa02c0085fe8a4103cf8fc0085f8cc02941c0fc0085f7d3badf016,
    0x400861debff4544b8000861a51f8757de0800861300d9ab95b4c00860e063922514c800860cf,
    0xc31f3ac2400864f70f8dc999c0008636be6393703a800862d2a31e72ef7400862a1bc82a09a1,
    0x0086907c9e9ca2490008673e08ff62ca180086698f53178d0a800865d8fcabc40f6000865321,
    0x5c2832240086baac0020011f40086b426e4ec0001c00869ab43033c16d8008695378a5a3f3d8,
    0x086f680ddef4a9380086e5c7fd59833e80086e27234ad5df100086d9c564cc84df40086d2cc2,
    0x5fc96bc00872bde044f5f60400870566c207119280086fc000b1dafa480086f6f18b7de63c00,
    0x875155ebdeaa78800874cbe8aee6e730008746bc2418198a4008745f361763ebd80087395c50,
    0x4379200087708246a1cb05800875fe9abeeb86a400875e9508ed436e4008751c3c1e001bd000,
    0x7bce5bd39b87c80087b0a4781159dd80087a1f1cc9682d9c00879e3dedad8d8cc0087993d068,
    0x6a3280087fa80db623094c0087f7321024ae41c0087e02b10a9565540087c471aece9ce20008,
    0x4bca9ad57c5f800884511228dee32800882efbf5449f48000882980bed781370008810e6b108,
    0x52b80088883af83500f800088871fc6409c0d00088660c83ee137e800884f7deda3e75000088,
    0x0ec40fc887180088a9198fe199de00088a2fb79577619800889d9f55d5d4d98008890fe5ffd8,
    0xb680088dcd8c1e9efa080088d5a353ccd85a00088ccfae58c6e9e80088c63be3e5020900088b,
    0xb5cf0ed2de000890599d00a506e8008902015012c19a80088ff0c0b41d6a700088fe1f1cd794,
    0xd000896344984cf6f9000895a93f277a95d8008953eecbd3eb4b00089427264870ecc8008917,
    0x8cc11761080089af99279036fb000899c80db017775000897091ad916ccd800896ff7caf1d24,
    0x80089e34518a1f40580089dd76373d4e8c00089dd00ab52964400089bd8f0f9e435200089ba3,
    0x69a1f77d8008a1791a0a895da8008a134c9fdf3ca100089fcac00ae6d0780089f3d91f5175cb,
    0x008a44ef9025abe58008a3b5bf61095790008a2f63aeb99cb70008a2d906574e8fc8008a209b,
    0x3130c050008a6b49a57cb4128008a658c282a19220008a5c36f783ca540008a4f5dc8a329f68,
    0x08aa47eae8f0bb70008aa172875827b88008a92b84df8a48d0008a8ef61195ddb20008a73be5,
    0xbb029a8008af09f24f22a0a0008acfec40253ee08008ac5a5d7d03c190008abb2bb234571e00,
    0x8b3a82361478bf0008b38c0f74129788008b04aba4a11ab58008b01b504a4d10b0008af2ee71,
    0xacb578008b62b8b28d96190008b52482a6c11b90008b4444e9f642aa8008b41a247a4b8eb000,
    0xbd6d3337f98998008bc8affe25680c0008b715ee9d58edc8008b67cd593bd07a0008b673b9a1,
    0xc3dd0008bf823f3faa2b10008be072d82cc9960008bdab9f86076f18008bda60547cdf8b0008,
    0x4dff3ca5769d0008c1bbfdacf11878008c1a3da34096440008c109ab3acae330008c0458db3f,
    0x2a98008c8bf6c382f8d50008c7abcb43467720008c7635ef414f540008c4e09c77eb7db8008c,
    0x56b1a33cd6d0008ca3181e00383b8008ca1a0a793d21c8008c9a940cbb1d2c0008c98196ab19,
    0x248008d03b01c53947a0008cd08b8df5d1df0008cc2f8e839d75a8008cc2120e4546430008ca,
    0xc183d53f238008d2532344a53340008d24e44b37814b0008d19eb88372d890008d1522125713,
    0x78008d72d9e3bc87230008d67936351dd1a8008d59b57553f6378008d36faeeb97e9d8008d34,
    0x34e644e378008d8b00c64b66608008d881dcba4f0948008d852c5b12a1300008d84ba4f5196f,
    0x0008dd3918eb18e8c8008dd181053f2fc28008db5f4e73009ee0008daee881e37c108008d982,
    0xf36435498008e0cd5b327a1230008df39186106beb0008ddca2bbbc95000008dda3b8132abfa,
    0x008e27874817ef6c8008e25e43017b7a98008e220223424ee78008e1e93cacf14e98008e194e,
    0x25276d50008e74cda3e9f47a0008e63c1136084880008e526f48fffae30008e301fc221649e8,
    0x08eacbf9b3fac740008ea2e94c1947448008e950f0901f7348008e94509918813a8008e78915,
    0xe496e00008f492d7b25d3e70008f0f03e07fe66b8008ec19b884b54cd0008eb435c9d0d06f00,
    0x8f7c59f05dae3a8008f79ed3c0e9a790008f6683ba2201138008f615684ed6f970008f581605,
    0x059f30008fc827916245540008fc6e4fc367dc00008fb6bc1facfa4b8008fa1aed88e089d000,
    0x02bbca846069f000901781476e48730008fe5644c9cea480008fd81d3993935e8008fd6ec783,
    0x36f80009047eb90b78b910009043668bde11cc000903d350d63b845000903491c6cb78a40009,
    0x6db62abdcfd40009068bd8648fade00090650fc30a19fb0009052f4d28e4475000904e77a4a4,
    0xb9d00090d7e7fa4c997000090c90f0efbe47700090befde3d64ef700090b9926e90bc4200090,
    0xf277ad6a0820009120d56c7ce78100090fa3f042b334600090de8107d941b200090ddecebb77,
    0xa8000915b3d1fe77891000915a975d4ed4d3000914af4699810910009132f99bfe18f2000912,
    0xb1eab8fc2800091c23bc717ca7800091c1f01f43ac76000918ca710a740f9000918bb55328f5,
    0xf00092004b626fb9a0000920002e113e01400091d599b68dd19e00091c579560798d100091c3,
    0x418ea4a1d0009252ff694661c9000923af79b2135d70009219f50cc5c9580009218ccc33aad5,
    0x00092e61def36379d00092b1b349c016d000092a4e0bd15f50900092726d41725df100092618,
    0xa8e6402e000935c33fb2283e80009353cff5d87deb0009310cc7de1fbbf00092feb798dfdb0d,
    0x0093b838a073411c0009395d6cb7a1a260009389113dcf608d000937d5f744534b5000936ef4,
    0xbd1b4d900093f15a324a77dc00093c96c9051ae9c00093c96bc602371900093c56682075cab0,
    0x0946cddc47aca4a000945b346919ac2b000943601a8c18ab2000940fa03c13b6ea00093ffa66,
    0xd86d1800094ffa4bfca0c1100094d19deacd785000094cd9d20ee1a6100094b4915372cd9f00,
    0x95e3df0f3df71e00095c53b807a8a0f000955e47b7fe979300095405e0d98ba1d0009515cb59,
    0x90ea600096169666ec388c000960e9f8f57b70100095eb894ed7ead000095e5cb82312b01000,
    0x71ce15c65646400097103a33cb0ef000096e2ed90b866ba0009676fd69dc5bcb0009664a9ff8,
    0x269f00097ce2f0adaaf0100097be9737d4961f00097bd8f8fa9b0db00097216ebae976de0009,
    0x438ba62f554a000984376ae35ec720009832c58c0dd87a000980ffa15c1be4e00097efdd501d,
    0x8cd00098bee93d916328000989d68fbf70c4e000988592b84f7ac1000984eb73b1ea8b900098,
    0x9460d86174200099937deee70c2b000997e717d79f34700098f0782114616900098c2e45b43e,
    0xa10009a398e7e888ca10009a2318a41ba1e30009a07b49897162700099fa5a3b4141f900099d,
    0xe612ee12b80009a7b99a2dbdf740009a6e3ed56c7d070009a6a1df2ceeeaf0009a40e7e8ca3e,
    0x60009b005dd46534870009af372db67a9450009acfa7962ade190009abca3fcfca0910009aa7,
    0x80b09f00c0009bd171a409567a0009b9694006e79a20009b70e7fbd2da910009b6599e765d45,
    0x0009d8075722cf2f10009d7e04c193b6e60009d41a92756e68d0009c97323ab5cabb0009c78d,
    0x2f2774960009e7012f31451d20009e0bbb24cc114e0009ddea00a37e4550009dc231e43cd956,
    0x009f75c07aa34f410009f48e2533db3580009f36273ae1cf7a0009f2234e2f09cd60009e8bb0,
    0x6b85a04000a01b5d364a33d6000a00a2b54d699900009ff49605aa04050009f9b51028e83e30,
    0x0a20e2394162332000a19a0bd51984fe000a16331d21241e6000a0cd25996cc8d0000a099e73,
    0xfbff22000a32e53c2a8fda6000a2ff50e1677ce6000a2b40b2b997d04000a29528c84ed38000,
    0xa5fee9d4202b66000a5223e548a4bce000a4e188ac423a1a000a4cb3dabd68212000a4b07393,
    0x26a9a000a6e235579afee0000a6ab1a8945d56a000a65bf721b871f8000a62077be0c794e000,
    0x9350acfbaed20000a924b14bf87ae8000a8d931e8f25b10000a7ced6dcb2f008000a754373d0,
    0x000000000000000000000000b3db4cb0d9b8b6000ab1142b3e20970000aaf52eb5c0ce66000a]

def c519yN : ℕ := ofChunks76 [
    0x272b00006655dc123b4f2200065ee48721ff03600064dfee4fba971a0005f3bf63838e990000,
    0xafaf0327844a80006a9ff5b0c7bb7600069f62427e27bc000068994af07d145a0006735e459c,
    0x42080006d36c912621f3c0006cd135d98ddf920006c055d9082ece00006bfbb80dda43860006,
    0x7ba3a66699420006dae8075372c560006d65f5004b8e9e0006d5e03e01dc2500006d4e293131,
    0x2ab0006fe9a51012af560006fb3738546f5740006f49d268bf0f600006e8073eb8cf9c60006e,
    0x05cbd53117a0007022a8dd928e4f000701e5e61e0d8db000701840d76c48a000070035e8c0d2,
    0x5300071e4a4ec6d3276000710c845fc43c5700070f6014fefa8a5000708f6eaf7c9a3e000708,
    0x3830d6797400072a42ae08f67b9000725f787e34d7cb000722e6b87b8bdf200071e5a81d5805,
    0x9000738a9cb4569dee0007389286a70c045000733ad95b47f2f80007336bcb229fc4f00072cd,
    0x5596a7e0a000745db654f601160007424c618def4580007406aac9c2744300073faaaaabb0eb,
    0x00074ee21ad2bdc9000074c65441c0a17800074c0a7f1ecff4800074a248ee51778f00074917,
    0x999707b1000759fc3e5b6d4ab000758cb39bc70b57000758c43e4ff9e6600074fbbae2db8e5f,
    0x0076e7bb95b9dbfb00076a569046c89d400076a3e032297d2b000769316af53c17b0007662f5,
    0xa5c6777000771306bc7f48720007712fd6929b4ad000770e07334c049a00076ecf404997a4f0,
    0x077c13b78666dad00077af7f4ae7afde0007743ece6909f600007743389b2ff828000771328f,
    0x8802930007902b07300490d00078e7740020388d0007831744f66f301000781117e717782700,
    0x794d8cbfd559700007946c7ed226178000792e58fafc061b000792507889da6cb0007910211a,
    0x04a2600079dc755c6c59d000079965fea657301000798b4dfb1e38150007986faeb4d5ecd000,
    0xae14b527e7b750007a833efc3b69730007a45a4110040ec0007a29a63cfeb9af00079f56494f,
    0xa9a10007b3864ccbdd98f0007afc1ae44a6f820007ae4dc6cb597de0007ae3cc39c5d2b30007,
    0x851558ce19af0007b616b15b58d340007b602d1dbe95db0007b5ec3a82f29a80007b53ee4a52,
    0x9360007bcf4175d2b8f80007bce0d1692ea510007ba918649077e30007b9790c3e5799f0007b,
    0x84a8214a8290007c5027a97514e70007c2576b5ddc8080007c14d1425df7690007bfa2a56093,
    0xed0007c89d767c06d2b0007c802809dc8f080007c7ba879f4d5b70007c71092694fe260007c5,
    0x40d30a4f2b0007d07d10d0fb85b0007d000e99b1a1820007cd95a768380230007cc93d3e4945,
    0x80007d4f7f28c464340007d4aedff6443450007d305719da64f90007d20b112ed07350007d13,
    0xb23378dd30007d9498efff5b5f0007d68b429e207bb0007d574eb44c29ea0007d525889bf22b,
    0x0007ddb0c8a0cecfd0007dd0610efb7bef0007dcf492b23d2640007db4fa037e099a0007da2b,
    0xf397f3d80007e3303adb799ce0007e11b36c0368730007df662bd54e97f0007ddd5916d7601f,
    0x007e61f81c917daa0007e5708000feb0d0007e512bea22a1f00007e4e956cc68af20007e4449,
    0xe08e29d0007e8bcd327d46a20007e81d70ee0af5f0007e721db1c652f90007e65b1d0d334390,
    0x07ed25bba324e6a0007ec714425b79830007ebdbb3b9ef8e50007e9f427acfc3590007e92dc7,
    0x9b6c550007f2ecdd73d957e0007ee973ba51de290007edef16e495a630007eda1fe74ea36000,
    0x7f7636ba6c7cea0007f716d27b2df340007f6cfb9533f4320007f5d0cfb5b85310007f5c5819,
    0x9dc820007fa1dbb1da35e40007f9f3cdf258f500007f976102e09e100007f8059a9681023000,
    0x011355b35d7208008005fc97a69c560007ff6be59c51d480007fac706e2318510007fa47720a,
    0xdb71800802f5679666c65000802ae0480d496980080233d188d2195800801cf729a86d4d0008,
    0x87098e3d92eb00080629ab13095c58008061a24eb05460000804835ba618a04000803b3b56e9,
    0xa0080080c5ebabe7819280080ba93c0daa91500080af6e0415145400080aca61da382a680080,
    0x8262608ce7f00081192a6cc765cc8008103d0303f0c2300080fa4845a83ffa00080e523c30d3,
    0x5680081bc53dfa7848b00081a1b484bfe16f000819101adb89496000815dcf75b30229800815,
    0xa861894e01000821bc119f5e8cf00082121d982f66bd0008200e53e1cda3c80081d106eb353b,
    0xe000825c0e0aeffad1000824cabdf035c57000823c1b134ae55e00082312275b1647b000822e,
    0x012522eb180082a02d8e9c0ed800082731cc5aea8e18008269c7e40ccea48008260ceb903877,
    0x00082d708c2313f0e80082d34700ee9d8d80082c6766a94a5bb00082c48c627d89ca00082b7e,
    0xc3dea67880082f69bd777a20780082ecd9f17c668900082dd3d0d96b26580082dc20005fb04b,
    0x0083517672191327000831bdc541079d68008318c4eccc5b198008318a5a48f6d52000830f56,
    0x3c7dd6180083786ea7cd75748008374256a7d8aba000836b93e24aea6f0008359aba685ba850,
    0x083b4f205526b2f00083b49058f47c4800083a7ee1ada7a05800838edcd1633ef30008388dd6,
    0x43c32700083e13d0696c72980083ddd13211fdd200083be778387026280083b65090be5c2800,
    0x84087a953375b200084026c4a7af35380084015b83c859b700084011f8888081d80083e36784,
    0x7262a00084362f0bddca85000842621b1ebc730800841b9f72c029d5800840d9a61f982bd800,
    0x4c568b6dc38d000084b2e83c1bb0d10008461dfc2a494bd800844fd516de0f3d0008436f9be5,
    0x448880084e7c606aa974880084e189ec6d47ef80084d44f39125fe080084c6a271dcbef98008,
    0x5562377a3dc2000855291181ada10800855069cf3d042d80085225b2c01236f00085060bb2b5,
    0xcc20008586e309250fc400085800734939589800856a059e85bc38800855aeb74b21bac80085,
    0x57960d25ce6800859780148ba6390008594a2dce54e2d000858fc07061e261800858e3c6f5bd,
    0x9700086004b1752a5a500085e78e3cfddb9c00085d5b7780f47f080085cc84b1d0742b80085b,
    0x51a692e381000865f17c7b8f306800864315fa52fd3f000863133c07e0ceb00086160e86c833,
    0x7800869b4c1d0471c3800868b3195f4b503800867f87df60b822000867e3f95d35a5c000866d,
    0xabb2d183780086afe1f083db5080086ae7b94e6f4fa00086adf64eeb431380086a969b0a2380,
    0x00086f137571405de00086dd93e2f0d4d100086c9477a4b9ee400086c62b3c241d5600086b09,
    0xd2da53e90008738d28425ed59800872b89edf5d04d800871c0bb71c87df8008713e4469d82f0,
    0x00879af25f51fea50008791ac44b93c0a8008779196cc4b172800876e5f5decaefc000873cea,
    0x327df2a80087d6620145edab00087c0e02f919fde00087a06478c3d0d2000879f399c9443710,
    0x087fd0df89c200600087f18f77fe884680087f0293a111f8200087e2492cade4f600087d97bc,
    0xe13b06000882adf947d361a8008823343228501900088169b38b4d1ec0008801a33b3aa4e900,
    0x887bb6ae4882a10008879348b33b3330008863e180186d860008857f2cd9e70c700088403aa0,
    0xe6c5800088abf53f08be9000088a1d9e31d99800008883600dd0acfb800887bbff633d175000,
    0x8bddd4cb79e4500088bcc4f3fec2a900088bc6ce834c33880088b82ec021754680088b52966e,
    0xd90500088d40863027cef00088c8f8c064225600088c55fcc6dd90880088c1bcf785e7d08008,
    0x17a8d7a27da8000890ebf8e2c0e0000088fba4e113c9e100088e6671f1e928a80088d9413215,
    0x99c000895bc311fb5b77000893e8ee5dd0000000892a341f3b54238008917ea87ae77cc00089,
    0xde9b333a0ea00089d0564b09607700089c1071b49334a00089bff5402403fc000898a670cb4b,
    0x820008a3d83431209d40008a2128510162020008a0c8387d420b500089ff4cf3d3019500089e,
    0xe751db0bce0008a61800ff937ea8008a5fe12efb99d80008a5480f0dc346e0008a41cbbcaf33,
    0x58008a93eeb7a303d58008a8c4603da99490008a7c35294be2750008a76dc4574e1d10008a62,
    0x5622248350008aab896afa10c50008aa68d2481f5888008a98dc2f3d76d88008a94ad3a9124c,
    0x0008abe6fe464c8d40008abc41c78398f60008abc3e7c993be60008ab1fd133665a18008aae6,
    0x445babe58008affd2000179648008af1ec1dcb71de0008ae1f1c4872dcf0008adab397fa2d8a,
    0x008b1660d5628c4e0008b1579708dab620008b143f34b9ff370008b08268b79e7b00008b0433,
    0xcfa73160008b593a30dc942b0008b473d02e531300008b39012712103b0008b206544a6d4800,
    0x08b8216b773a6120008b772bfd0a37230008b6029229948d90008b5b821c5cc6b50008b59d57,
    0x0018f70008b95f3de9ab2c20008b86a2ca688e8e0008b859a0918bfb28008b825996dddc6b00,
    0x8be96f498ba8930008bd3a191afdb210008bcc1a03a6ff358008bc165ccecc2c48008b9dd9a1,
    0x51d338008c1e6101ea67348008c1a07d28ee4cb0008c14f4bf17a0a28008bfe4dd5add554800,
    0xc7813b09d12280008c69fadb3e1fc50008c302222ac85fb8008c2cc748f9d31d0008c21a0229,
    0xa1550008c9ebb8d11819b8008c9296ed241f890008c80a3c517eac30008c7e8511a742440008,
    0xc545786e1c738008cbd00ff24a90b0008ca9b255580dda0008ca61967a5eb780008ca219d4ef,
    0x2eb0008ceb45fbf103fc0008cdd99293df00b0008cce5c412376f40008cc570dc3c16040008c,
    0xdf3b991c0138008d0eea9691b5d40008cf6635f3eb5e38008cf44a250be16c0008cf2e0f518e,
    0xd28008d5e667041c2f20008d4f6e7ca537ab0008d4ed8e0a3503e8008d431f612973258008d1,
    0x78166535990008d64a6f61e7c0f0008d63d2db369da58008d5fe328ea3f8f8008d5f258c474b,
    0x20008d98be0e265a8a0008d955fe934ee868008d7d5edd0b10f40008d78b4cb97ebff0008d66,
    0xb9e0ae5a90008dafac4d720de30008da83ba68752eb0008d9e72bc2d982a0008d9a7618435d8,
    0x8008de29c53f092ec0008dddb702e2e2cd8008dd26224a2f9398008dc95fd890ec4f8008dc83,
    0x5521d5ec0008e1f93c582cddc0008e15c04a6e10988008dff2e4ea12ecc0008df8b776dd8519,
    0x008e33a47c9b8f840008e32f4c1771bb68008e31b650dbb6018008e280d47e3e56c0008e2365,
    0x3536de68008e65ec36b807300008e4fa069bf58010008e4161e2814b288008e403cb03057f60,
    0x08e93bd8a92f45f8008e801c1101f2738008e7e5ed61369cf0008e761d934486a30008e6657c,
    0xf631a90008edcaa36ccc2630008ed717f745e44e8008ec801430ca5e80008ec5284380c5a380,
    0x8ef28c6c09ebad0008ef1d2ef49a5b80008ee655ff7377018008ee07c0563d8a98008edfbc72,
    0xd25308008f16a610c7aa178008f1533520b38700008f0a2f8af9b66e0008f025d632a4835000,
    0xf4d8a153560e68008f330c1f2d76900008f2c2985f0247f0008f27bbd09366ba0008f205d290,
    0xab038008f627b210b59978008f56ab8dbd90998008f4e73012b56430008f4e6d53d20d618008,
    0x767d721e39930008f724e147f06bc0008f71b24eded6d98008f6cca25b6e7d18008f64868b12,
    0xa6d8008fb6d4a2b8b4300008fb65a59ebfbfe8008fa6bac4ba2fa38008f7a3b21731b520008f,
    0x0b82934d6c18008fde9bd945eac48008fd6a3d2d32cd28008fd4129011f2b18008fc7aa7f607,
    0xab8008fefa949e78a438008fe82baf525dfd8008fe3f5217c66000008fe144c2d567928008fe,
    0xe6534535c100090165ad7534c0a0009003214eb6f5f40008ffe6caaa581f80008ff155c64836,
    0xf000909dcaebe271f8000908d6ce33220ac0009059d462d924c200090527a719d1dfe000902c,
    0xb5fa9ec2000090e06e6615dd8300090df73396a735400090c83fb7a9af1e00090c4d34f87ace,
    0x000911b97b5811a5b000911b968073de7d000910a84c37bc05900090e4699ab3ac6100090e0e,
    0x0abe3cca0009139eeb4c95838000912f2016beba37000912d75bc2b84ab0009127fb7c240b1a,
    0x009185a9dee6a82a000916ea56c1cce8200091655a2e6228ce0009157c7b23a5eb600091540f,
    0xd36f6ba00091b73c15741d3a00091b70b9eda70dc00091a18aee299b4f000918ceb742970860,
    0x091e11d4d2931cc00091d9383116ed0a00091cca6dace4d5d00091c93293551fba00091c298e,
    0x8f58e70009204b2ca042a6a0009200238793bc7c00091f25c8efd792d00091ecea3031616f00,
    0x9244d7e5fae40e0009242aa9a3b97d7000922e43dff359e8000921456b88e93260009205c906,
    0x330de00092907f2e88dd4f000928cc8050bc0c700092816a1907ed4f000926de7bcd59bcf000,
    0x2a5b6add978c500092a1ee86e3134400092a1de5e89cdc7000929ec447ec0bf0000929462152,
    0xe32800092bf00add343fd00092bc0844d76b8000092afe2e1bebd7900092ab22db310f600009,
    0x126efb2dc8f5000930586dc00a4c200092fe2957d20f8200092edc3b548226100092c5099da6,
    0x2c7000933ab9ede966b10009333cd35bb3a0b0009318ad822420ff0009313f231d8275d00093,
    0x1e66cb76fab0009351d6f9c270690009350938f5e9ec20009343b5fd1bd8060009342352558c,
    0x1f000936b35f8cb163b000935b4873290709000935aa7f522a3ba0009358ea93020578000935,
    0x2b5b56f10100093a44f9031fa9800093a403d81a54a00009378923e769f76000936e81ece841,
    0x000093f2e3bafaa1d700093cf240d65f2e100093b0ace674d74b00093ab9a9c7eae6200093a8,
    0xe1d42e78e000941527a55d4f6d000940f4351a3a35c00093fee069d06a2d00093fbbce2d8780,
    0x000944a263c29aae900094454c126348f200094388ca78eb9ea0009425c1b866974f000941ad,
    0x553149260009473e55df4a29600094706c12d8b11c000946ea359a80f67000945f2c668c5842,
    0x0094a3fb1922b23500094a3918cb9d58e0009491c6613385d0000948457b72ccadb000947771,
    0x2040b3400094cc3f8128aa6600094b5abe60cdb6400094aed04c9d2ab200094ae52040012ea0,
    0x094ec72ed85a45900094eb11db43b6f900094dd5b2e50e1d000094d9467162fd8400094cd56c,
    0x29aaf4000954d8e27d2880800095327c071beee0000951c989814ed4a00094f1dae5321f2600,
    0x957c788a1ac0a8000957561f3fa520f000956ec5552db12d000955ad821121dbe00095526b2a,
    0xae00e00095c405f379222200095ab6d3570407e00095a52d4d4babc20009594d3fa49c0ec000,
    0x608097a6d0b5600095f772261b883f00095e4dfa02a360200095db9ac154975c00095cc31c6a,
    0x16cc00096317e11b6a0df000962d66e037f2020009629981aba2214000961affd13e50750009,
    0xa10198f4f97700096892b3e1149380009686f297369fdd000965fb195392134000964c0c40cd,
    0x21900096bd04d16ff89e00096b81f2c448ed000096b69cb47bbf0e00096a29d3553dd8e00096,
    0x0f46fb0f51600096cbfadb4b571c00096ca416b55be0300096c74006fd863b00096c00a2a04b,
    0x2c00096da72b255bf8800096d88a87b567b800096d82b4f497cd000096d0fadb9b65b200096d,
    0x2169d45f2500096f39507ba118100096f0564458599300096f01bef11b98f00096eb7c21320e,
    0xd000970fe7d41eec90000970b76e743470a0009706ad960cf01e000970501192d03a200096f6,
    0x78d83fbfc000972d7ea2d6763c0009724191e23d5ea000971bccbd3d9514000971455f0dff62,
    0x00097610b56c006970009753186ace8b2c000973b32248bded50009732f0eaa0f6d6000972e2,
    0x630bf53e000978127a6f9b826000977f1432e9f121000976ceaedfc19ab000976bb580c0b6ff,
    0x0097a2428c059734000979b35183327080009790972d2cc0d2000978e8b08c80778000978723,
    0x995769800097d8479576b3d500097bedefff450b000097b94c89c9e61c00097a76343049f2f0,
    0x09825b54f510fb60009822c0f2e750100009818023a1c8582000981759058398dc000980116c,
    0x2371d2000983a495cd024b600098394710177d65000982a1448daaddf00098295c71f7fabe00,
    0x9896f6572252880009894735979005c000989410f2cc898b000988fecd272de92000986816e2,
    0x4f4c500098ddc6073a98de00098d204e8561fd500098c7c840bf8cdc00098a5e61e408001000,
    0x918290ed2dd7600098edfdb834cc4600098ec77527e3f6a00098ea7a840885c000098e21da78,
    0x63da0009940dc237734e7000992dc5f3484dfa000992a3c6c384c5e000991dd08d867d620009,
    0xb2b5a51fb59400099789c5f596c360009974fd552d56fe000995af3e72892880009941e6d854,
    0x19000099e041da7a329f00099d2a0e2d79aa900099cb1806efd01100099b5fc5fdc53ac00099,
    0x17098aa044300099f76bb7d3cee200099f563f5e34d1800099ee3ce9243cd200099eae7e96b2,
    0x050009a1ab7fec7c7cb0009a15765c1001de0009a103d989339c00009a08793c005ce60009a0,
    0x262ee7860f0009a4ba5c5f087160009a3acddf77fa800009a3511bc64e0050009a25e398126a,
    0x20009a7a09e2705e930009a74b7f9dbb3cc0009a699692ab08180009a5fa178df4a5f0009a5f,
    0x91d8cdfeb0009a9897f62f84200009a8e8f08c830b40009a8e81f36bae660009a811434bd2e9,
    0x0009aef82091c98340009aeac363f6d8720009ae8c63fdb975a0009ab334aadbf29e0009ab14,
    0xd4666dae0009b2b9074336d4e0009aff71832fcf7a0009afa68f62242a20009af6926d7b4cf0,
    0x009b984612e8fc4c0009b8a1ed43be9e60009b7788ee1055b60009b761824ed6e620009b7004,
    0x48282e60009bbede5114a6360009bacce92d1a0800009bac477c2c57090009bab1f064a27d40,
    0x09c100a9504f0560009c0d7f1bac02340009bdd300a1771a40009bdc45aac373fb0009bc976c,
    0x3ba89c0009c2face88af4f00009c2e5348ed68a00009c2d79cd2fa1120009c271427d4641600,
    0x9c884dffa70b390009c4bce72eaa4900009c4a0f9f4528080009c49c274d7dbd20009c33d1b0,
    0xafb1a0009ce17f5a4a6a480009ccf59635a7e900009c9a52e8f3cffd0009c8d0f1e13591a000,
    0xd2d82a036409c0009d1d3264d75c040009d0330b997c1a20009ce622174ad3520009ce5d269a,
    0xcf7c0009d45ef5d6cbfaa0009d3ee24f73040e0009d367072756ab20009d33d0337116da0009,
    0x98cedaefbe4a0009d985ab7fbed160009d90c39448c2d20009d6958bb73bf9c0009d66eea55a,
    0xb6f0009db564af59099e0009da868e6c2fbf80009da1342e107bc20009d9fd8582ce4dc0009d,
    0x47d5aade0620009dfa8c026728860009ded4c3bb825d40009de693c5aaefc00009de4b81d4af,
    0xae0009e20739ae898e20009e1da1fe82dd8e0009e0d504ef365a60009e05dcca4447df0009e0,
    0xd4bc89525c0009e47f89767f4e80009e358cc4ddadf60009e296bad7c14740009e20c70801c8,
    0x40009e72ffe8ef97c40009e722e8a09a78c0009e64890a0cd6060009e63c2c6f46cc80009e5f,
    0x33ef0caf20009ea064e899283a0009e8900744706dd0009e8554f00c69290009e79e47632fbe,
    0x0009ec62a25a690120009eb1868301bbb20009ea7e749c3fbac0009ea3b2ac7220ea0009ea2d,
    0x8a4817c50009ef83b917d9bba0009edbcd75c5c3c20009ecf096cffe7c90009ecee4deb33242,
    0x009f4c1b3195fe360009f3dabe90ddda80009f2ba916bdbb2e0009f2a2a0568649a0009f1772,
    0x40721120009fd5d1fb8b1d310009fd32e65d17e0a0009f82cda41305510009f508b7c68e6880,
    0x0a046a7260b2eee000a0335a3353e1c4000a026e578f2de460009ff8fef58eeaac0009fe54b4,
    0xc9c772000a07366a1428c5c000a0643fc1021474000a05f08cdc1cb60000a04cd528dbfd0400,
    0xa0a5625f994092000a095b35db6e718000a08a3d19c9bad2000a082c13a048c12000a07f1442,
    0xe3ac8000a0c3e62ad17a3c000a0b3c76a1908b0000a0a6647693a6dc000a0a6305a166824000,
    0x0f087a90f0f60000a0ec3144059dd4000a0ea74b699aaa4000a0d3c4b60909c8000a0cf3a712,
    0x4346000a153c95600b348000a14ed08a4a4d6c000a1021ef10164f2000a0f282e7dccb06000a,
    0xa41b6814923a000a18bed118f52c6000a1861ef912abb8000a1836470598e6c000a16d00b749,
    0x6a2000a1c836d259611a000a1be8a7356174c000a1ba5902123368000a1b59b27bc05a4000a1,
    0x050c5c1de14000a1ecf964effd0c000a1dc40edb8da88000a1dbb2115663ce000a1d0556cba2,
    0x36000a24d2660eb0758000a23f06416b74ee000a210f9b6bbe128000a202bc08e96fb6000a1f,
    0xa40ae98086000a286cd2019bc10000a2683d67f3741a000a25714b29e30c4000a24ee669cba8,
    0x8000a2f09f8dc03fe8000a2d762c11956d0000a2b220cbfd6750000a29c86da84acdc000a294,
    0x954a9a228000a397cc8dc568c0000a38e48124b9e62000a35fee0d6b9c14000a3273ae888fbf,
    0x000a46cf5d1a88a40000a43e2141d2c392000a42d32f833aa74000a41ba39a3c0b74000a3e63,
    0x1527f250000a4b069a11ee1a6000a4a942a0135c96000a48b168197553c000a47714d021825c,
    0x00a4e259ebb570b2000a4d26c26db2962000a4c8e89e2c6a1c000a4ba1064e03f16000a4b819,
    0xa35535a000a542a942c52e2c000a536ebaa8f5134000a5231987277934000a507751b494b6c0,
    0x0a5bd185130fde2000a5b6c3173fa460000a57b6338340e10000a5775bd05fe95c000a559a29,
    0x8bc918000a5f8056bce9832000a5e8dd50cb205a000a5e7d53d861aa8000a5d5b32062b92c00,
    0xa618fca917e4e4000a60f610e3f3104000a60af5d428558e000a601ce8a301a94000a60120ad,
    0x6af78000a66243d3398dde000a6608d0a4520dc000a653342f3f1546000a6508a49dc5d38000,
    0x6baad2bbf60bc000a6a4f8d05b478e000a69057006443dc000a679963d2495f6000a665833b8,
    0x871c000a6cd33d9f4f09c000a6ca07076f6e04000a6c67d144bed7a000a6bbaa693f430c000a,
    0x5c05401e9f84000a73ebd55a816a0000a733b95ddfacfa000a7094319ebef00000a6ed99489e,
    0x78e000a7b7ebe9b47fd6000a7b600715056a0000a796ca111d259a000a78459e21e72f4000a7,
    0x4cea1a80380000a862e392426cac000a7f2fb522cb8ea000a7e415c9961072000a7cf56d78ee,
    0x5e000a9243caf065c24000a90d9355946254000a8cafa86e83392000a898e4c2032038000a88,
    0x364fe9e36c000a94972e4d69674000a9442e22f96dc0000a93ee1b76d8cd0000a93475daa00f,
    0x8000a9f66a1e34f3d0000a9dd80832fdb94000a9aada647eb5d8000a99367adf1e34c000a966,
    0x90f4d50d0000aacf24bfb5d4c4000aa959a21eca400000aa74edf1f5504c000aa30287962749,
    0x000ac0732814bbdbc000abfceacd94ee98000ab8131fd2f333c000ab0aad9b305dda000aad30,
    0x78c21f6c000ac5186f8d07e08000ac4f52a2c49be8000ac3b9495394564000ac313cd4b27cde,
    0x00ad9bf1131e91ec000ad27af774eb7d4000aca593dac5edbc000ac77345ba5a75c000ac5f2f,
    0xb6a3324000ae6f738364aa98000ae66f55b8506b2000ae201c9e3e8464000ade6034f56c7c60,
    0x0af90319da6ab90000af5853a818bd34000af423940368610000af093d94397134000ae7b0ef,
    0x7761b6000b2091113f91530000b109075c04e870000b02e418c31c816000b016d4006da53400,
    0xb5d7447832d0b6000b48bbdb8efa6ba000b46025221a7794000b36e732fcb1296000b318f1e6,
    0xed0f2000b72afff54cdb8a000b6cab5202cbc0e000b6a8c910006b12000b69e3debd81722000,
    0xd50cc3c398e34000bd0eae62c4bc36000bcc21d3824e116000b93c7d78a31c7e000b7a8f15b7,
    0x000000000000000000000000c63ef5796ac870000c4a5f59e5f87b0000c2f1d89e51fd8c000b]

def c520N1 : ℕ := ofChunks76 [
    0x4a48ce008120cd683431a590098f078211461690007c513168a64d8820077084ab03da921800,
    0x502a6957f3f618006bc24f269ed742000875fe9abeeb86a4009b70e7fbd2da91000821190042,
    0x77cb90086698f53178d0a80085c1cc24bdc277c009395d6cb7a1a2600078de918dd59893c007,
    0x88bc72c85f0100087f7321024ae41c009c97323ab5cabb000771b4653875b0280081c55b3d53,
    0x4a40007c832cbc631b926006eea31b08e49a6000786f7246dfa8d80008b38c0f741297880060,
    0x7936351dd1a80065034ddfffab7f00075ffd6d31df250000745676267da72e0006f957c5c0f8,
    0x3c00086d2cc25c28322400940fa03c13b6ea0006dca30ec52a00b0008274375ed3181d2008d6,
    0xaba4a11ab58008d2532344a53340008e0cd5b327a1230008e526f48fffae30007b47631fe90e,
    0x800079e9973d5934f9c007b1ac6e4d70a4a8007f0268ddbd564640088dcd8c1e9efa08008b04,
    0xff5d87deb000a2b40b2b997d040006bc0efc15895db0007cac5eeeaa1ea720074e27aef53bc2,
    0x9008c7abcb4346772000788894655af598800798dc8a18560c1c006e07bf4b548e370009353c,
    0x82312b010007f5067113f81fa1807f8b128eb5ed30600840d97b6d142b98007eaa8532a0164b,
    0x00895a93f277a95d80067384efac3e2f5000984376ae35ec720007c633df68352fcc0095e5cb,
    0xc05aa75c0077f67baece1bad800861300d9ab95b4c0083747f2579a303e0081fb4e966c88e93,
    0x0872bde044f5f6040084db40d943b0898006ed9f171625cd50006fb3ea4c90b5cb0007aab2eb,
    0xe2e7ea0007824bc810786e7800a0cd25996cc8d00008ec19b884b54cd00080d91af5c5b76fc0,
    0x9068bd8648fade00074fdeaad4e80fc00074a71d19e690d9000901781476e48730005b0da4ef,
    0xed9690f07dc719d9ca09ea20082e69e7a42c70f2006c77f0b44380550007b7191dbfeeadcc00,
    0xc78e3153c7ee8e008a4f5dc8a329f68007a11c6c95846c9000862a1bc82a09a140080179c2fd,
    0x4a1d0007ac8fbf6368f2c40086f680ddef4a93800918ca710a740f9000854efa9ac6baa80007,
    0xf418c227382880083d1d85e2d22b680088c63be3e502090008fe5644c9cea4800092618418ea,
    0x0c0400988592b84f7ac1000a2ff50e1677ce600076655634017b6780095c53b807a8a0f00075,
    0xd8cd11d0b76800841c32179b3d5140074a22cba20c3a18006f93c7ad9406250007c3a9242a02,
    0xfc00073710802d3297c8007f4791022884b570070c62b9bf5493e0007a70772bd9eb234007c0,
    0x38cc11761080076cb60fbbbe7838008eb435c9d0d06f00079ea9798a4acc4c0065501c1daba1,
    0x960084d1863019623b000884f7deda3e7500007f1c6b0e72badfd80661d270572256800089ba,
    0x165f07dddc006f9d2b21175b9c00072679c6a68446400075bfab0e7fa4ad00081b0d9855d5e0,
    0x8009d8075722cf2f10007be70bdf9bf6184007826ad3ce7f0b700070d724dbe493b080081e1e,
    0x0d628a69d00709c3c37cf8701800813f53e033e98720083c424b5de189400075919cfad966c5,
    0x00793d73f62700e2c0080a8e3da37a6b100087395c505fc96bc00971ce15c6564640007edbec,
    0xf2774960008bd6d3337f989980068a137e082698d0009219f50cc5c9580005f25e63dce17040,
    0x066d87f5dbdf6a00007bce0a015ec4fd8008b3a82361478bf00068e5d048eb54350009e8bb02,
    0xa4d10b000689ca62cfb27fc00077816829a0ba2980066601b9af22da300067c3ebd9ab2c5900,
    0x8c4dff3ca5769d0008c9a940cbb1d2c0006b6852bb4793d8000a19a0bd51984fe0008b01b504,
    0x2c19a8007b83802e57e2a3c007e16b8dd1197032006fddb07d5661af0008550b3a6c3cd43c00,
    0x0d7e7fa4c99700007a01b780a2edc3c008fa1aed88e089d0007a0a7da8ac5690800890201501,
    0xea5e00079548d233c2f23c00b3db4cb0d9b8b60006bd667743e18b50006f36922bca77160009,
    0xe74759280c74800884511228dee3280078426f3f2247cd0008a2f63aeb99cb700084c58a4b86,
    0xf50000a20e23941623320007f7475e535429d680810b70b53a5da32007c8dcd2b37ba9160085,
    0xc000b1dafa480084f7baa98e5560400809ac80229ba462006e76cc9e46051100074237cc2f29,
    0xc58007a67d4feab4ad5800943601a8c18ab20006efc9dbc59c93200078aa891ce62b9700086f,
    0x157dabd11ed0091c23bc717ca780007cd85e16e36a280007c3832ea2af5d04007226dbfc0a2e,
    0xb1007ace16254acc73000836e5e3476efea2007802e65fb175d70008c109ab3acae33000815b,
    0x0db0177750009ff49605aa04050006b06f4f8acf48f00071ac9778b7320b80081836a0483eeb,
    0x000a32e53c2a8fda60006d084608e986190007638aa1676d7fd0006bfacae1766441000899c8,
    0xefacdf5d400795edd132431c8c0088b0ec40fc88718007bc1f741ccb087c0094ffa4bfca0c11,
    0x007f22b3679b19f3a009a6a1df2ceeeaf000767c729aa9fe038007361e1fd9e681c00078614a,
    0x280521b80077be33da37196100093b838a073411c00083e624ca4930c0c0081e86db16a906ea,
    0x06c3839ad978eb1000746e5260a41fa300078f0fed6b3192ac009acfa7962ade19000701eb0d,
    0x6c11b900091c1f01f43ac7600089e34518a1f40580076d1012f5a4d0f0006b71cded3fcf0a00,
    0x7486463fbf13ff80074b9f90ef9de9b80080f4bc59b0141a2807a2e75afa062dac008b52482a,
    0x183e1e007e8e0d0a9ba8a03009f48e2533db3580006fc021ff0d66450008e27874817ef6c800,
    0xd12bcd3f54bee0006aed807df71cfd00076cf4dce246437000742efd92cf41e5800836ec043c,
    0xfd1f8007361cdeda46664800737ed398b7817f000853059c081fb5000086b426e4ec0001c007,
    0x7012f31451d20007d229acb9d723d2006dcd9ab7f9e44b00097efdd501d269f00077f6fe7a09,
    0x07000087c471aece9ce2000755d1ab981e1568006a3e17ac9514720007e6268289cba07f009e,
    0x37deee70c2b0008172ac9cbdc50bc007a7f50a46e332b40077028dee22884c00079fd87d7de5,
    0xaf0007e2e87eed47bfac008a1791a0a895da800914af4699810910007e953cbfbcc02ff00999,
    0x9d20ee1a6100080c12501a927b128078297d9987e59740077a980e33bb90e00067793bbfd9ca,
    0x000063e5438798d48c0007dde9d182c6bcf80098c2e45b43e8cd000a4b07393fbff2200094cd,
    0xe71bb029a8008f6683ba2201138008d852c5b12a130000816d648d9aab1e00079c78ec3e77f2,
    0xd0082cbfdea871f13a0078dae6d1fa6fc40008d98234e644e3780069acf4cdaa9e330008af2e,
    0xfe25680c000655205a1cf9a4b0008065119844483938086d9c564cc84df40081315c933d6e8e,
    0x00822505b4e41f81a0081e3eb7fb25bc49007bbfe414d4d007400862d2a31e72ef74008bc8af,
    0x9c016d000082b06553f5fd170007c7ec507a4de784009a40e7e8ca3ea100083894bb5d3fb1a2,
    0x0a00a2b54d6999000076c2b0e9cf67430009252ff694661c900090ddecebb77b9d00092b1b34,
    0xd1b4d90008055f2e2924ed08c06a62b3dd76b66d0009ddea00a37e455000811134eeb4285b90,
    0x6f28e7d5ea268f00070fbbed71fa68480073f4fc69971df6000748e5e3552898c00093ffa66b,
    0x23719000836445984b9d498005fc5cb3682bae8000854c82884f9a3bc00790165015db05fc00,
    0x45f4ff1ecb8040008c1bbfdacf118780094d19deacd7850000a099e736b85a0400093c96bc60,
    0x38800007e4f2b934543fa1008b41a247a4b8eb0007a5e0488e4e12f40076603a808c947b0006,
    0x7d5f744534b500083cc83af76de830008e950f0901f7348009af372db67a9450006bdc934195,
    0x2f08008af09f24f22a0a00077a226a1d0f69c8008a92b84df8a48d0006a946169b999f100093,
    0x7eb90b78b91000801ddc36b7ef7aac0955e47b7fe97930008a209b69a1f77d8007f06df8eddb,
    0xb20008609c7f8a0aa02c00984eb73b1ea8b900074744ee6fadcea8006cfc2c97d960bb000904,
    0xafa9cb16ff2008dda3b8132abfa0007cc999d88631640006a2ce7108c97ae0005f2a1867621d,
    0xf50080e571cd571d8f0809abca3fcfca0910008f581605e496e000089fcac00ae6d078007e5e,
    0x4ed7dda6a80076d9b3c0f1f7308008196d88803294bc0063f2bd6205bc9a0007f2b15054f823,
    0x8007929dd31b31ca4400601db3c2bedda5000903d350d63b8450007d1f1aa4d4ddfb200706a7,
    0x13852fcb3808695378a5a3f3d8008557aa9ca1284e40071b56bfaa1308780081ddbf0e67b647,
    0x008a6b49a57cb4128007cb71820347ef9400a924b14bf87ae8000588bae55e8903a0007f78af,
    0x6609662c007ca4afea6a640de008d15221257132480090650fc30a19fb000846158f8dbe1e5c,
    0x0a6ab1a8945d56a000789f964344a1eec008db5f4e73009ee00070df9afd4f929c000860cf9a,
    0x9f748f800882980bed781370008d34c183d53f238007893821cdb323780065954dc0cec1d300,
    0x635b020fb4ef640007fc5659705b955240824fc5c1acc8af400700d2fa2eefa3400076917c54,
    0xd8d8cc006f1ea67314c81100082cbade80c8f0b000830b2ef189b5b1a00808d856b761bd0780,
    0x907989db05d03c006dd4a12f2070c5000874cbe8aee6e73000960e9f8f57b701000879e3deda,
    0xefa38008e25e43017b7a98005d7a47d62c66ea00075a5b37ae77e688008589a3f2f028d34007,
    0x881dcba4f094800739d1a34cfb063000980ffa15c1be4e00076d2b2772c3ec080072ba75af14,
    0x0378006d4d3429fac6f90007ed2b6ec7628731007a53cb585c0841c006bd43e308d8a270008d,
    0x768194f343040065975c81f769ba0008abb2bb234571e0007f89dd70968338a0075d81234106,
    0xfb40087e02b10a9565540083a2fcdaa3af39c00776109a571c6b08008c1a3da34096440007a6,
    0xa971c94308000796194bc6ffa094007e0962a4a2a8c2a008a73be53130c0500084490d5e792b,
    0xc000771d94728eecb08008890fe5ffd852b80086e27234ad5df10007e83504752c79f60055b5,
    0x990427307800945b346919ac2b000820563c710c0d5c00849531c5b14e298009c78d80b09f00,
    0x0009a2318a41ba1e3000802d95cc7250ff780a65bf721b871f800078cb69ba21e22d80074812,
    0xf41baac000084ccf537b32c4f800860e063922514c8007cf3a4a75e321aa0068fbe074d0370c,
    0x00915b3d1fe778910008b673b9a1acb578006ec2e28fad95d2000804ee2d02a92060006fe8e9,
    0x195ddb20008746bc2418198a4006c24ae86953a6300093c56682075cab0005c66734d1958d60,
    0x0870566c207119280082a527a7b936ff60075ad9cfcd1812f8007ccd4052d805410008a8ef61,
    0xd64ef700096e2ed90b866ba00065d163f8fd585c0005f29a0b1dc513e0008e220223424ee780,
    0x556cdc37368bae0008eacbf9b3fac740008b715ee9d58edc80086f6f18b7de63c00090befde3,
    0x06beb0007b47089d20bb6a80083b09b46177dad200783d7da55cc96980091c3b1eab8fc28000,
    0xe83f18f36c02e00079bfc11b9b5063800884bca9ad57c5f8007a6ca82742dba28008df391861,
    0xe014000810aded99a24705009b9694006e79a20008bdab9f86076f1800935c33fb2283e80006,
    0x2b9db0ee53bba006ad53136c827150007b7f0c9c17facb4009f9b51028e83e3000920002e113,
    0x16d4007ac81d577e63e70009b6599e765d456000826d488964eed52007c93497e2447a4e007c,
    0x05e0d98ba1d000850661541d14c34008053c41f44b87f3407c5b67f43740ee20085f7d3badf0,
    0x500008ea2e94c1947448008c7635ef414f540008b67cd593bd07a00081f0eb89fba075200954,
    0xfd69dc5bcb0007e4eb542b827dff0079e799e7f8a9a40009a07b49897162700076cf47a875e1,
    0x6e006a454ee469eeac0007207e6b3e1ef958006e8f64aeb1682700094b4915372cd9f0009676,
    0x59b6cc49100088d5a353ccd85a00085d369ec22eea200082a3302d8c10a900082e8d89a12bed,
    0x90067da74305f4a680008c0458db3fc3dd000906db62abdcfd40007019cf8192655b000765ea,
    0x80b462b9f008266e45ddec9f46005ede34201dd82c00071aa329da56a390007e7a46aaa4df6b,
    0x00834a209bde5363e008e1e93cacf14e980083755c5e1fbcb0e007e712cec46555e5007ebb8c,
    0x26d9da400084f5d45f28fa234007de422c125dde2600743e4323bdf71f0007623142b6c6a800,
    0x074d344d99d4e430007e447c18557094f008f79ed3c0e9a790007d0d1224bf58eaa006f5d782,
    0xb7417ac0093c96c9051ae9c00088ccfae58c6e9e80091d599b68dd19e0007b95ba84e1f83280,
    0x77e60eccb6debe00089af99279036fb00077108f7c13d58700086907c9e9ca2490008239beac,
    0x5f5090008cd08b8df5d1df0009f2234e2f09cd6000820c425a1bb3b180069455d1af818c1000,
    0xba0c697d8c8280007cc39f17d62e1ba007643598c562c2300085bddf7faaae2a40092a4e0bd1,
    0x8e8c8008ddca2bbbc95000006cd7995e3eaa240008b4444e9f642aa800788034a8e318b70007,
    0x15a53803c64d800a6e235579afee00008a3b5bf61095790006cd2715c78b9c80008dd3918eb1,
    0xd0700083ff48349d123560091c579560798d1000737c19eae319e50008d84ba4f5196f780077,
    0x800caeaae61c007791b6e8119fd10007f2de2448bcd6100068a80e4f31cd6e0009a6e3ed56c7,
    0x1f0006e9944f71477e500089427264870ecc80072a893a239f64b00075ed3ed71e824a8007cd,
    0x44984cf6f900088871fc6409c0d000819e052c18cd737006d054c3e5c711400097be9737d496,
    0x500068b250cb66b08d0008f615684ed6f970008145a3e25e6f85d007b2c0bb0b044c74008963,
    0x5d7d03c1900095eb894ed7ead000078e22c4ad3143d00085205e2d58d742c00737dfe800fc38,
    0x3007aeff5ea228916000806ae811e31c6c90067e039dcee011d0008127985317adb14008ac5a,
    0x0f4118b44006f606f362e9a7d00087fa80db623094c006d18555a1824040007e5cc9e6f91c4e,
    0x006aef92ca00fe970008953eecbd3eb4b0007376b825e4a111000861a51f8757de08007c6165,
    0x1086a328007b6d2d595f6b7a0008d72d9e3bc872300087708246a1cb05800530b29b4cbc4ea0,
    0x092e61def36379d00060ef8a6d6529070007df73305bf3aa60006d815dac7dc1320008810e6b,
    0x40aa4c0007f777a9a5f5df6e007f8ffca099bb0ef80912f277ad6a0820009dc231e43cd95600,
    0x696857a476c910000757dfb215b1df38007ca0d5321c78b46008ca56b1a33cd6d00067be2fe0,
    0x534870006c84b5e00f202600079a5c3533ad968c00726a3edf0cc70e8007b342c3d856427000,
    0x8bfb3a1d8dd34800565f702aa3aba00009132f99bfe18f20007fbe74224fd0953809b005dd46,
    0x1409400695b01587f84e400098bee93d9163280006c6ee4f56338010007531397b095bd38007,
    0x0f03e07fe66b80071fef41cc830ba80082d9ddf47b93676008d8b00c64b666080085182062bd,
    0xd69000a16331d21241e600069ca8d14f4e2180006efbe99e8b982f00082baba918980156008f,
    0x6b1cebde9b300084f6b400402e1680088883af83500f8000624f9418f099230007524d9e2c31,
    0x4ac008a2d906574e8fc800923af79b2135d70007e239ff8c5c7ed50060e9a68921cb0e0006a5,
    0x4c9fdf3ca1000828484f1835a3ba0084564a18746c0c00072a72a3daa3140000781ceda5c3ff,
    0xb0007ab616c2b1e5e9c0072b4fef8cf374f0008e63c11360848800079aacfc4baf7380008a13,
    0xf0f9e435200081b70d49c944d610060fbe27e81d90e0008fd6ec783059f30008bda60547cdf8,
    0x8009052f4d28e447500072bfcc56141d1080087b0a4781159dd80068c96c71d50eb100089bd8,
    0x66ec388c0008dd181053f2fc280081ee776434869f9008c4e09c77eb7db800889d9f55d5d4d9,
    0x00808121696b469d9805c03e300ffd3ae0007f8bb1a087c1121809f75c07aa34f41000961696,
    0xe286cf6800902bbca846069f0008917b5cf0ed2de0007d83b498e6dd5a6008349fa8cb9e6eec,
    0x08175d4c9b46019500773248970d246d800865d8fcabc40f60007b4b62b2bf5450c0077b6580,
    0xe5e7b800072578eadb774398006ac1837973202e0006e14e4079cc35a000609cd87e27912b00,
    0x700cf29134efd1800745ca316aeb8330008fb6bc1facfa4b800861debff4544b800074e83982,
    0x7e913c0087bce5bd39b87c8007c2cefdb17dcb3800795dd8362fa4c6c006f553e1020851f000,
    0x310cc7de1fbbf00092feb798dfdb0d0007e897c368c7591b005ebe8209528c4a0007adbda7b7,
    0x79114007d66b83fce5e9f000a754373d026a9a000673c4ea9b8f8b20007360681ac772168009,
    0x5a59107e867b000745936ead1c79b800662ce4f675e05800086baac0020011f4007e70d6f53d,
    0x36fc00897091ad916ccd8006a2a266a56781f00081e6cea9be719780074e9159f3ca6ea00067,
    0xa3f042b33460008bf823f3faa2b1000839f06e71cbb66c00718855d346855680084a12734510,
    0xa000097103a33cb0ef00007e9cd1a3a26d5f3007e8076b487323f0007c982c67434096c0090f,
    0xe08ff62ca180085f8cc02941c0fc007b670d555e0b3bc008d59b57553f63780092004b626fb9,
    0xfbf06787f1ae326d5d0008a44ef9025abe58006458dd5d52b1450006f554e623e928b0008673,
    0xb165346918008307e3ce721648600890599d00a506e800840103dc09619e8008006e579af9f0,
    0x00060e6ce4fff8fc70008e194ef36435498009bd171a409567a0008fc8279162455400074463,
    0x4dd9939e40089f3d91f5175cb80088660c83ee137e8007f1a65148623670009664a9ff890ea6,
    0x008fd81d3993935e8008b62b8b28d9619000904e77a4a436f80007ecbb8d13819f46007e73f9,
    0x3e209700005ca13f2ce9f7a200083516ee78c6ca22008aa172875827b8800a62077be0c794e0,
    0x0989d68fbf70c4e00085fe8a4103cf8fc00a5fee9d4202b660007846f4ab3d8fecc00ab1142b,
    0xc9b24d0007d05a64d3150a5200836c5f33d4820a600882efbf5449f480007e3ea063e83af970,
    0x7cdda33c2a34c54005f65a70a6bebc60007cc6d3d273ea13200903491c6cb78a40007e12ddec,
    0xd15c340064ed26f4a5afdc00089dd76373d4e8c00089dd00ab52964400088a2fb79577619800,
    0xad98728d851f34007c668f029ffb9e600946cddc47aca4a0008f7c59f05dae3a80085e65e127,
    0xc6110008751c3c1e001bd0007902be91f3ee5c80065502aae74bf0000097ce2f0adaaf010007,
    0x8b73cad578224007fef5e18a4a5ffd709e0bbb24cc114e00084795b10c31dd7400661d994195,
    0x8d41f08e7891525276d50008d03b01c53947a0007b88eb9e93eaaa40087993d068437920007b,
    0x2120e45464300076831531be57150006ae89161ea76540006d2c44d8101580000800b06694f4,
    0x8c800864f70f8dc999c000758a1e2a9dee9480061d19cb5d9f28300060ab310cffd93f0008cc,
    0xeb88372d890007d450f6444364100085557bbff48f60c008e301fc221649e800716208914653,
    0x48008a5c36f783ca5400083b2a592cf174600078dd733b9fba098008636be6393703a8008d19,
    0xede9120e800092726d41725df10008daee881e37c108007578ca49c6c3e880085b3562a24abf,
    0x800846bc9d5de66ad400a29528c84ed3800007eca3d5af4b119e00824e40b0c79e90600846b9,
    0xfbff5ee04007ec1d9bef9a1b290066f49b1ba2b30500057df9195281b040007d77a0ae394e75,
    0x00936ef4a8e6402e0007ce4c0717876eb6008d24e44b37814b000846e11d0bf946f800834cfa,
    0x3777816a008be072d82cc9960006c56129702b9a00006caa37f749b96100071417864d22c5a8,
    0x07a1b3a5486f2c8c0076879c25a3c3b58007bb77978a7bae84008ca3181e00383b8007eac5c2,
    0x8ae19d00090de8107d941b200068ed6d142310830006d49418239e6d3000625bae75fdcc5d00,
    0x848ac49bacebe94007cee7d89aad18880098438ba62f554a0008ca1a0a793d21c80076715600,
    0x8ead280084c90f5947d9db40080e6418aba61b890072e6ba84f542fa00063fba09ff1b22a000,
    0xa658c282a1922000810724406b469a10078aa6f19fca0f8000aaf52eb5c0ce6600073d9815c1,
    0x3dcb00099d9460d86174200090b9926e90bc4200075b3f2aeff24558007f956d7be8f229b808,
    0x80569e4f3d010009043668bde11cc0007e70ff978b5a55d0088ff0c0b41d6a70007f865f9d57,
    0xd20c007a62a3fd09f54b40068a88d843331bd0007657ef9e4a2908800a5223e548a4bce0006d,
    0x4aceca23d870008247827e2ad145e008e94509918813a80099fa5a3b4141f90008258eafed72,
    0x5f000623c8202d9d4f9000a8d931e8f25b100007d632f56c254666007224b59b1881bc80069c,
    0x47e07d27de6006dc5d58cbbb6180009aa7e612ee12b8000779072cfbd97618009218ccc33aad,
    0xa8009f36273ae1cf7a00097bd8f8fa9b0db0008e74cda3e9f47a0006e683d45919c65000824c,
    0x1cc9682d9c007de3d721e8d8f840093f15a324a77dc0009d7e04c193b6e60008cc2f8e839d75,
    0x80088a9198fe199de0009389113dcf608d0007d096973346444c007072b6329cf3f300087a1f,
    0xc31f3ac240073c42c7c09a7eb8007653e9ef5b73510006f63d0b4f94bad000875155ebdeaa78,
    0x009d41a92756e68d0006e23288b6071c200088fe1f1cd794b68007026e1788b52a4800865321,
    0x7d2514f0007a5c5a08bf202d80083db52faced34300074c8b06fafe8d3000875e9508ed436e4,
    0x084cb22bf60b6848007356c74219976d0006fc63037a0245500083e58290694da4a0085f038d,
    0x14a88c2008d36faeeb97e9d8007601c5896aa39c800915a975d4ed4d30009a398e7e888ca100,
    0x6d090e99a368fe000a9350acfbaed20000842a70cf7d9f06c0084d22e437f808a4007f391064,
    0x976de000918bb55328f5a8000896ff7caf1d24d0006385368908e6ac000819d5e245fa269300,
    0x96030e703545a000997e717d79f347000869ab43033c16d80061329ae474f54700097216ebae,
    0x6d180008aa47eae8f0bb700085d71de6fd4e66c007f6a57f218b4064807c9b080b2e10802007,
    0xe5c7fd59833e800a7ced6dcb2f0080008323f472a3a203c007e78e2011dd6b1c009515cb59d8,
    0x28a5607c7572b5266151200800bb2c87221ec1486a862c129a0732000830521faf9706c00086,
    0x163b69837e02007477b6117837040007502cd97260b0d00076cb16df132e9c8007fd2c19ba36,
    0xc780065067c77bf7c470007a3d8abb10f92880090c90f0efbe47700095e3df0f3df71e0007c9,
    0x5d364a33d60007565f5e9cc5fa000081a88052ec884c40078ad034a8551e9800727c96163b1c,
    0x86006fe9f9326ecba30007dbd644fc4bce3e006c93aef69294860007ba4cacc3e91ee800a01b,
    0x890093ae3400660cde5aa5062a00085b07ca2ace814c0083bcf699bc3246400827f9d207d4ce,
    0x0008c8bf6c382f8d50008f492d7b25d3e7000694323c9de61dd00084f79a80ef789140083880,
    0x32014fc5b00818b9d966a4150e0080bf6ba7c51080d008281dafddf6d94200a4cb3dabd68212,
    0x8076b3f2b0bb1caf8007ec2668017755b50065431045fabbf7000626ca7cfa60dc9000811dc3,
    0xc367dc00006e6052e6cac283000680554faa3c3c10005e387d919cdfe000080e03368104b104,
    0x07ad607e930ee48800814729c4f02fe2400733530db4b17180006e6c95be090ad20008fc6e4f,
    0x773a358007b1ccbddeaa9f4c006e019f94afbfda000a4e188ac423a1a0008257b98a1f394a40,
    0x8c98196ab192a98007eafc5a48da77d40073b03093d6bd230009832c58c0dd87a00077bb5a71,
    0x63ebd8006ed00cbaeeb8ac0008acfec40253ee08009a7b99a2dbdf740006f10035ec3f1b5000,
    0x120d56c7ce7810007e8f762af5dd775008502fc2c3e1c4b0007c964b3fbcf8f88008745f3617,
    0x0000000000000000000000007c0350a2c928e700077a35cac93333a0006a61e48fdca8bb0009]

def c520N2 : ℕ := ofChunks76 [
    0x4a48d0009120cd683431a6000a8f078211461680008c513168a64d8800087084ab03da921800,
    0x502a6957f3f618007bc24f269ed742000975fe9abeeb86a000ab70e7fbd2da90000921190042,
    0x77cc00096698f53178d0a00095c1cc24bdc278000a395d6cb7a1a2600088de918dd598940008,
    0x88bc72c85f0100097f7321024ae42000ac97323ab5cabc000871b4653875b0280091c55b3d53,
    0x4a40008c832cbc631b928007eea31b08e49a6000886f7246dfa8d80009b38c0f741297800070,
    0x7936351dd1a00075034ddfffab7f00085ffd6d31df250000845676267da72e0007f957c5c0f8,
    0x3c00096d2cc25c28322000a40fa03c13b6ea0007dca30ec52a00b0009274375ed3181d0009d6,
    0xaba4a11ab60009d2532344a53340009e0cd5b327a1230009e526f48fffae30008b47631fe90e,
    0x800089e9973d5934fa0008b1ac6e4d70a4a8008f0268ddbd564600098dcd8c1e9efa00009b04,
    0xff5d87dec000b2b40b2b997d040007bc0efc15895db0008cac5eeeaa1ea700084e27aef53bc2,
    0x8009c7abcb4346772000888894655af598800898dc8a18560c20007e07bf4b548e37000a353c,
    0x82312b000008f5067113f81fa0008f8b128eb5ed30800940d97b6d142ba0008eaa8532a0164b,
    0x00995a93f277a95e00077384efac3e2f5000a84376ae35ec720008c633df68352fd000a5e5cb,
    0xc05aa7600087f67baece1bad800961300d9ab95b500093747f2579a30400091fb4e966c88e90,
    0x0972bde044f5f6000094db40d943b08a0007ed9f171625cd50007fb3ea4c90b5cb0008aab2eb,
    0xe2e7ea0008824bc810786e7800b0cd25996cc8d00009ec19b884b54cd00090d91af5c5b77000,
    0xa068bd8648fade00084fdeaad4e80fc00084a71d19e690d9000a01781476e48740006b0da4ef,
    0xed9690008dc719d9ca09ea00092e69e7a42c70f0007c77f0b44380550008b7191dbfeeadd000,
    0xc78e3153c7ee90009a4f5dc8a329f60008a11c6c95846c9000962a1bc82a09a100090179c2fd,
    0x4a1c0008ac8fbf6368f2c00096f680ddef4a94000a18ca710a740f8000954efa9ac6baa80008,
    0xf418c227382880093d1d85e2d22b600098c63be3e502090009fe5644c9cea48000a2618418ea,
    0x0c0000a88592b84f7ac0000b2ff50e1677ce600086655634017b67800a5c53b807a8a1000085,
    0xd8cd11d0b76800941c32179b3d5100084a22cba20c3a18007f93c7ad9406250008c3a9242a02,
    0xfc00083710802d3297c8008f4791022884b580080c62b9bf5493e0008a70772bd9eb230008c0,
    0x38cc11761000086cb60fbbbe7838009eb435c9d0d06f00089ea9798a4acc500075501c1daba1,
    0x900094d1863019623b000984f7deda3e7500008f1c6b0e72bae0000761d270572256800099ba,
    0x165f07dde0007f9d2b21175b9c00082679c6a68446400085bfab0e7fa4ad00091b0d9855d5e0,
    0x800ad8075722cf2f00008be70bdf9bf6180008826ad3ce7f0b700080d724dbe493b080091e1e,
    0x0d628a6a000809c3c37cf8701800913f53e033e98700093c424b5de189400085919cfad966c5,
    0x00893d73f62700e300090a8e3da37a6b100097395c505fc96c000a71ce15c6564640008edbec,
    0xf2774960009bd6d3337f989a00078a137e082698d000a219f50cc5c9580006f25e63dce17040,
    0x076d87f5dbdf6a00008bce0a015ec4fd8009b3a82361478bf00078e5d048eb5435000ae8bb02,
    0xa4d10b000789ca62cfb27fc00087816829a0ba2980076601b9af22da300077c3ebd9ab2c5900,
    0x9c4dff3ca5769d0009c9a940cbb1d2c0007b6852bb4793d8000b19a0bd51984fe0009b01b504,
    0x2c19a0008b83802e57e2a40008e16b8dd1197030007fddb07d5661af0009550b3a6c3cd44000,
    0x0d7e7fa4c99700008a01b780a2edc40009fa1aed88e089d0008a0a7da8ac5690800990201501,
    0xea5e00089548d233c2f24000c3db4cb0d9b8b80007bd667743e18b50007f36922bca7716000a,
    0xe74759280c74000984511228dee3200088426f3f2247cd0009a2f63aeb99cb700094c58a4b86,
    0xf50000b20e23941623320008f7475e535429d800910b70b53a5da30008c8dcd2b37ba9180095,
    0xc000b1dafa400094f7baa98e5560000909ac80229ba460007e76cc9e46051100084237cc2f29,
    0xc58008a67d4feab4ad5800a43601a8c18ab20007efc9dbc59c93200088aa891ce62b9700096f,
    0x157dabd11f000a1c23bc717ca780008cd85e16e36a280008c3832ea2af5d00008226dbfc0a2e,
    0xb0008ace16254acc73000936e5e3476efea0008802e65fb175d70009c109ab3acae33000915b,
    0x0db017775000aff49605aa04040007b06f4f8acf48f00081ac9778b7320b80091836a0483eeb,
    0x000b32e53c2a8fda60007d084608e986190008638aa1676d7fd0007bfacae1766441000999c8,
    0xefacdf5d000895edd132431c900098b0ec40fc88720008bc1f741ccb088000a4ffa4bfca0c10,
    0x008f22b3679b19f3800aa6a1df2ceeeb0000867c729aa9fe038008361e1fd9e681c00088614a,
    0x280521b80087be33da371961000a3b838a073411c00093e624ca4930c100091e86db16a906f0,
    0x07c3839ad978eb1000846e5260a41fa300088f0fed6b3192b000aacfa7962ade18000801eb0d,
    0x6c11b9000a1c1f01f43ac7600099e34518a1f40600086d1012f5a4d0f0007b71cded3fcf0a00,
    0x8486463fbf13ff80084b9f90ef9de9b80090f4bc59b0141a0008a2e75afa062db0009b52482a,
    0x183e20008e8e0d0a9ba8a0000af48e2533db3580007fc021ff0d66450009e27874817ef6c000,
    0xd12bcd3f54bee0007aed807df71cfd00086cf4dce246437000842efd92cf41e5800936ec043c,
    0xfd1f8008361cdeda46664800837ed398b7817f000953059c081fb5000096b426e4ec00020008,
    0x7012f31451d20008d229acb9d723d0007dcd9ab7f9e44b000a7efdd501d26a000087f6fe7a09,
    0x07000097c471aece9ce2000855d1ab981e1568007a3e17ac9514720008e6268289cba08000ae,
    0x37deee70c2c0009172ac9cbdc50c0008a7f50a46e332b00087028dee22884c00089fd87d7de5,
    0xaf0008e2e87eed47bfb0009a1791a0a895da000a14af4699810900008e953cbfbcc030000a99,
    0x9d20ee1a6000090c12501a927b100088297d9987e59700087a980e33bb90e00077793bbfd9ca,
    0x000073e5438798d48c0008dde9d182c6bcf800a8c2e45b43e8cc000b4b07393fbff22000a4cd,
    0xe71bb029a0009f6683ba2201140009d852c5b12a130000916d648d9aab1e00089c78ec3e77f2,
    0x00092cbfdea871f1400088dae6d1fa6fc40009d98234e644e3800079acf4cdaa9e330009af2e,
    0xfe25680c000755205a1cf9a4b0009065119844483900096d9c564cc84df00091315c933d6e8f,
    0x00922505b4e41f8200091e3eb7fb25bc50008bbfe414d4d007000962d2a31e72ef70009bc8af,
    0x9c016d000092b06553f5fd170008c7ec507a4de78000aa40e7e8ca3ea000093894bb5d3fb1a0,
    0x0b00a2b54d6999000086c2b0e9cf6743000a252ff694661c8000a0ddecebb77b9c000a2b1b34,
    0xd1b4d80009055f2e2924ed10007a62b3dd76b66d000addea00a37e454000911134eeb4285c00,
    0x7f28e7d5ea268f00080fbbed71fa68480083f4fc69971df6000848e5e3552898c000a3ffa66b,
    0x23718000936445984b9d4a0006fc5cb3682bae8000954c82884f9a3c000890165015db060000,
    0x45f4ff1ecb8040009c1bbfdacf1188000a4d19deacd7850000b099e736b85a04000a3c96bc60,
    0x38800008e4f2b934543fa0009b41a247a4b8eb0008a5e0488e4e12f00086603a808c947b0007,
    0x7d5f744534b400093cc83af76de830009e950f0901f734000aaf372db67a9440007bdc934195,
    0x2f08009af09f24f22a0a00087a226a1d0f69c8009a92b84df8a48d0007a946169b999f1000a3,
    0x7eb90b78b90000901ddc36b7ef7b000a55e47b7fe97940009a209b69a1f77e0008f06df8eddb,
    0xb20009609c7f8a0aa03000a84eb73b1ea8b800084744ee6fadcea8007cfc2c97d960bb000a04,
    0xafa9cb16ff0009dda3b8132abfa0008cc999d88631640007a2ce7108c97ae0006f2a1867621d,
    0xf80090e571cd571d8f000aabca3fcfca0900009f581605e496e000099fcac00ae6d080008e5e,
    0x4ed7dda6a80086d9b3c0f1f7308009196d88803294c00073f2bd6205bc9a0008f2b15054f823,
    0x0008929dd31b31ca4000701db3c2bedda5000a03d350d63b8440008d1f1aa4d4ddfb000806a7,
    0x13852fcb0009695378a5a3f3e0009557aa9ca1284e00081b56bfaa1308780091ddbf0e67b648,
    0x009a6b49a57cb4120008cb71820347ef9000b924b14bf87ae8000688bae55e8903a0008f78af,
    0x66096630008ca4afea6a640e0009d152212571324000a0650fc30a19fc000946158f8dbe1e60,
    0x0b6ab1a8945d56a000889f964344a1ef0009db5f4e73009ee00080df9afd4f929c000960cf9a,
    0x9f748f800982980bed781370009d34c183d53f240008893821cdb323780075954dc0cec1d300,
    0x735b020fb4ef640008fc5659705b955000924fc5c1acc8af000800d2fa2eefa3400086917c54,
    0xd8d8d0007f1ea67314c81100092cbade80c8f0b000930b2ef189b5b2000908d856b761bd0000,
    0x907989db05d040007dd4a12f2070c5000974cbe8aee6e73000a60e9f8f57b700000979e3deda,
    0xefa38009e25e43017b7aa0006d7a47d62c66ea00085a5b37ae77e688009589a3f2f028d30008,
    0x881dcba4f094000839d1a34cfb063000a80ffa15c1be4e00086d2b2772c3ec080082ba75af14,
    0x0378007d4d3429fac6f90008ed2b6ec7628730008a53cb585c08420007bd43e308d8a270009d,
    0x768194f343000075975c81f769ba0009abb2bb234571e0008f89dd7096833880085d81234106,
    0xfb00097e02b10a9565500093a2fcdaa3af3a000876109a571c6b08009c1a3da34096440008a6,
    0xa971c94308000896194bc6ffa090008e0962a4a2a8c28009a73be53130c0500094490d5e792b,
    0xc000871d94728eecb08009890fe5ffd852c00096e27234ad5df10008e83504752c79f80065b5,
    0x990427307800a45b346919ac2c000920563c710c0d6000949531c5b14e2a000ac78d80b09f00,
    0x000aa2318a41ba1e4000902d95cc7250ff000b65bf721b871f800088cb69ba21e22d80084812,
    0xf41baac000094ccf537b32c50000960e063922514c0008cf3a4a75e321a80078fbe074d0370c,
    0x00a15b3d1fe778900009b673b9a1acb580007ec2e28fad95d2000904ee2d02a92060007fe8e9,
    0x195ddb20009746bc2418198a0007c24ae86953a63000a3c56682075cac0006c66734d1958d60,
    0x0970566c207119200092a527a7b936ff00085ad9cfcd1812f8008ccd4052d805410009a8ef61,
    0xd64ef8000a6e2ed90b866ba00075d163f8fd585c0006f29a0b1dc513e0009e220223424ee800,
    0x656cdc37368bae0009eacbf9b3fac740009b715ee9d58edc00096f6f18b7de63c000a0befde3,
    0x06beb0008b47089d20bb6a80093b09b46177dad000883d7da55cc969800a1c3b1eab8fc28000,
    0xe83f18f36c02e00089bfc11b9b5063800984bca9ad57c600008a6ca82742dba28009df391861,
    0xe014000910aded99a2470000ab9694006e79a20009bdab9f86076f2000a35c33fb2283e80007,
    0x2b9db0ee53bb8007ad53136c827150008b7f0c9c17facb000af9b51028e83e4000a20002e113,
    0x16d0008ac81d577e63e7000ab6599e765d456000926d488964eed50008c93497e2447a50008c,
    0x05e0d98ba1c000950661541d14c30009053c41f44b87f0008c5b67f43740ee00095f7d3badf0,
    0x500009ea2e94c1947440009c7635ef414f540009b67cd593bd07a00091f0eb89fba075000a54,
    0xfd69dc5bcc0008e4eb542b827e000089e799e7f8a9a4000aa07b49897162800086cf47a875e1,
    0x70007a454ee469eeac0008207e6b3e1ef958007e8f64aeb16827000a4b4915372cda0000a676,
    0x59b6cc49100098d5a353ccd85a00095d369ec22eea200092a3302d8c10a900092e8d89a12bed,
    0x80077da74305f4a680009c0458db3fc3dd000a06db62abdcfd40008019cf8192655b000865ea,
    0x80b462ba0009266e45ddec9f40006ede34201dd82c00081aa329da56a390008e7a46aaa4df6b,
    0x00934a209bde53640009e1e93cacf14ea00093755c5e1fbcb10008e712cec46555e8008ebb8c,
    0x26d9da400094f5d45f28fa230008de422c125dde2800843e4323bdf71f0008623142b6c6a800,
    0x084d344d99d4e430008e447c185570950009f79ed3c0e9a790008d0d1224bf58ea8007f5d782,
    0xb7417b000a3c96c9051ae9c00098ccfae58c6e9e000a1d599b68dd19e0008b95ba84e1f83280,
    0x87e60eccb6debe00099af99279036fb00087108f7c13d58700096907c9e9ca2490009239beac,
    0x5f5080009cd08b8df5d1df000af2234e2f09cd6000920c425a1bb3b200079455d1af818c1000,
    0xba0c697d8c8280008cc39f17d62e1b8008643598c562c2300095bddf7faaae2a000a2a4e0bd1,
    0x8e8c0009ddca2bbbc95000007cd7995e3eaa240009b4444e9f642aa000888034a8e318b70008,
    0x15a53803c64d800b6e235579afee00009a3b5bf61095790007cd2715c78b9c80009dd3918eb1,
    0xd0800093ff48349d1235000a1c579560798d0000837c19eae319e50009d84ba4f5196f800087,
    0x800caeaae620008791b6e8119fd10008f2de2448bcd6100078a80e4f31cd6e000aa6e3ed56c7,
    0x200007e9944f71477e500099427264870ecc00082a893a239f64b00085ed3ed71e824a8008cd,
    0x44984cf6f900098871fc6409c0d000919e052c18cd730007d054c3e5c7114000a7be9737d496,
    0x500078b250cb66b08d0009f615684ed6f970009145a3e25e6f860008b2c0bb0b044c70009963,
    0x5d7d03c19000a5eb894ed7ead000088e22c4ad3143d00095205e2d58d743000837dfe800fc38,
    0x0008aeff5ea228916000906ae811e31c6d00077e039dcee011d0009127985317adb10009ac5a,
    0x0f4118b40007f606f362e9a7d00097fa80db6230950007d18555a1824040008e5cc9e6f91c4e,
    0x007aef92ca00fe970009953eecbd3eb4b0008376b825e4a111000961a51f8757de00008c6165,
    0x1086a320008b6d2d595f6b7a0009d72d9e3bc872300097708246a1cb06000630b29b4cbc4ea0,
    0x0a2e61def36379c00070ef8a6d6529070008df73305bf3aa60007d815dac7dc1320009810e6b,
    0x40aa4c0008f777a9a5f5df70008f8ffca099bb0f000a12f277ad6a082000adc231e43cd95600,
    0x796857a476c910000857dfb215b1df38008ca0d5321c78b48009ca56b1a33cd6d00077be2fe0,
    0x534880007c84b5e00f202600089a5c3533ad969000826a3edf0cc70e8008b342c3d856427000,
    0x8bfb3a1d8dd34800665f702aa3aba0000a132f99bfe18f20008fbe74224fd095000ab005dd46,
    0x1409000795b01587f84e4000a8bee93d9163280007c6ee4f56338010008531397b095bd38008,
    0x0f03e07fe66c00081fef41cc830ba80092d9ddf47b93670009d8b00c64b666000095182062bd,
    0xd69000b16331d21241e600079ca8d14f4e2180007efbe99e8b982f00092baba918980150009f,
    0x6b1cebde9b300094f6b400402e1600098883af83500f8000724f9418f099230008524d9e2c31,
    0x4b0009a2d906574e8fc000a23af79b2135d80008e239ff8c5c7ed80070e9a68921cb0e0007a5,
    0x4c9fdf3ca1000928484f1835a3c00094564a18746c0c00082a72a3daa3140000881ceda5c3ff,
    0xb0008ab616c2b1e5ea00082b4fef8cf374f0009e63c11360848800089aacfc4baf7380009a13,
    0xf0f9e435200091b70d49c944d600070fbe27e81d90e0009fd6ec783059f30009bda60547cdf8,
    0x000a052f4d28e447400082bfcc56141d1080097b0a4781159de00078c96c71d50eb100099bd8,
    0x66ec388c0009dd181053f2fc200091ee77643486a00009c4e09c77eb7dc000989d9f55d5d4da,
    0x00908121696b469e0006c03e300ffd3ae0008f8bb1a087c112000af75c07aa34f40000a61696,
    0xe286cf6800a02bbca84606a00009917b5cf0ed2de0008d83b498e6dd5a8009349fa8cb9e6ef0,
    0x09175d4c9b46019000873248970d246d800965d8fcabc40f60008b4b62b2bf545100087b6580,
    0xe5e7b800082578eadb774398007ac1837973202e0007e14e4079cc35a000709cd87e27912b00,
    0x800cf29134efd1800845ca316aeb8330009fb6bc1facfa4c000961debff4544b800084e83982,
    0x7e91400097bce5bd39b87c0008c2cefdb17dcb3800895dd8362fa4c70007f553e1020851f000,
    0x310cc7de1fbc0000a2feb798dfdb0c0008e897c368c75918006ebe8209528c4a0008adbda7b7,
    0x79110008d66b83fce5e9f000b754373d026a9a000773c4ea9b8f8b20008360681ac77216800a,
    0x5a59107e867b000845936ead1c79b800762ce4f675e05800096baac0020011f0008e70d6f53d,
    0x370000997091ad916cce0007a2a266a56781f00091e6cea9be719800084e9159f3ca6ea00077,
    0xa3f042b33460009bf823f3faa2b1000939f06e71cbb67000818855d346855680094a12734510,
    0xa0000a7103a33cb0ef00008e9cd1a3a26d5f0008e8076b487323f0008c982c67434097000a0f,
    0xe08ff62ca200095f8cc02941c100008b670d555e0b3c0009d59b57553f638000a2004b626fb9,
    0x00007787f1ae326d5d0009a44ef9025abe60007458dd5d52b1450007f554e623e928b0009673,
    0xb165346918009307e3ce721648000990599d00a506e000940103dc09619e0009006e579af9f1,
    0x00070e6ce4fff8fc70009e194ef364354a000abd171a409567a0009fc8279162455400084463,
    0x4dd9939e00099f3d91f5175cc00098660c83ee137e0008f1a6514862367000a664a9ff890ea6,
    0x009fd81d3993935e0009b62b8b28d9619000a04e77a4a436f80008ecbb8d13819f48008e73f9,
    0x3e209700006ca13f2ce9f7a200093516ee78c6ca20009aa172875827b8000b62077be0c794e0,
    0x0a89d68fbf70c4e00095fe8a4103cf90000b5fee9d4202b660008846f4ab3d8fed000bb1142b,
    0xc9b24d0008d05a64d3150a5000936c5f33d4820a000982efbf5449f480008e3ea063e83af980,
    0x8cdda33c2a34c50006f65a70a6bebc60008cc6d3d273ea13000a03491c6cb78a40008e12ddec,
    0xd15c300074ed26f4a5afdc00099dd76373d4e8c00099dd00ab52964400098a2fb7957761a000,
    0xad98728d851f30008c668f029ffb9e800a46cddc47aca4a0009f7c59f05dae3a00095e65e127,
    0xc6110009751c3c1e001bd0008902be91f3ee5c80075502aae74bf00000a7ce2f0adaaf000008,
    0x8b73cad578220008fef5e18a4a600000ae0bbb24cc114e00094795b10c31dd7000761d994195,
    0x8d40009e7891525276d50009d03b01c53947a0008b88eb9e93eaaa00097993d068437920008b,
    0x2120e45464300086831531be57150007ae89161ea76540007d2c44d8101580000900b06694f4,
    0x8c800964f70f8dc999c000858a1e2a9dee9480071d19cb5d9f28300070ab310cffd93f0009cc,
    0xeb88372d890008d450f6444364100095557bbff48f610009e301fc221649e000816208914653,
    0x40009a5c36f783ca5400093b2a592cf174600088dd733b9fba098009636be6393703a0009d19,
    0xede9120e8000a2726d41725df00009daee881e37c100008578ca49c6c3e880095b3562a24abf,
    0x800946bc9d5de66ad000b29528c84ed3800008eca3d5af4b11a000924e40b0c79e90000946b9,
    0xfbff5ee00008ec1d9bef9a1b280076f49b1ba2b30500067df9195281b040008d77a0ae394e75,
    0x00a36ef4a8e6402e0008ce4c0717876eb8009d24e44b37814b000946e11d0bf9470000934cfa,
    0x37778168009be072d82cc9960007c56129702b9a00007caa37f749b96100081417864d22c5a8,
    0x08a1b3a5486f2c900086879c25a3c3b58008bb77978a7bae80009ca3181e00383c0008eac5c2,
    0x8ae19d000a0de8107d941b200078ed6d142310830007d49418239e6d3000725bae75fdcc5d00,
    0x948ac49bacebe90008cee7d89aad188800a8438ba62f554a0009ca1a0a793d21c00086715600,
    0x8ead280094c90f5947d9db00090e6418aba61b900082e6ba84f542fa00073fba09ff1b22a000,
    0xa658c282a1922000910724406b469a00088aa6f19fca0f8000baf52eb5c0ce6600083d9815c1,
    0x3dcb000a9d9460d861742000a0b9926e90bc4200085b3f2aeff24558008f956d7be8f2298009,
    0x80569e4f3d01000a043668bde11cc0008e70ff978b5a5600098ff0c0b41d6a70008f865f9d57,
    0xd210008a62a3fd09f54b00078a88d843331bd0008657ef9e4a2908800b5223e548a4bce0007d,
    0x4aceca23d870009247827e2ad1460009e94509918813a000a9fa5a3b4141f80009258eafed72,
    0x60000723c8202d9d4f9000b8d931e8f25b100008d632f56c254668008224b59b1881bc80079c,
    0x47e07d27de0007dc5d58cbbb618000aaa7e612ee12b8000879072cfbd9761800a218ccc33aad,
    0xa000af36273ae1cf7a000a7bd8f8fa9b0dc0009e74cda3e9f47a0007e683d45919c65000924c,
    0x1cc9682da0008de3d721e8d8f8000a3f15a324a77dc000ad7e04c193b6e60009cc2f8e839d75,
    0x00098a9198fe199de000a389113dcf608c0008d0969733464450008072b6329cf3f300097a1f,
    0xc31f3ac200083c42c7c09a7eb8008653e9ef5b73510007f63d0b4f94bad000975155ebdeaa78,
    0x00ad41a92756e68c0007e23288b6071c200098fe1f1cd794b60008026e1788b52a4800965321,
    0x7d2514f0008a5c5a08bf202d80093db52faced34300084c8b06fafe8d3000975e9508ed436e0,
    0x094cb22bf60b6840008356c74219976d0007fc63037a0245500093e58290694da500095f038d,
    0x14a88c0009d36faeeb97e9e0008601c5896aa39c800a15a975d4ed4d4000aa398e7e888ca000,
    0x7d090e99a368fe000b9350acfbaed20000942a70cf7d9f0700094d22e437f808a0008f391064,
    0x976de000a18bb55328f5a8000996ff7caf1d24d0007385368908e6ac000919d5e245fa269000,
    0x96030e703545a000a97e717d79f348000969ab43033c16e00071329ae474f547000a7216ebae,
    0x6d180009aa47eae8f0bb700095d71de6fd4e670008f6a57f218b4068008c9b080b2e10800008,
    0xe5c7fd59833e000b7ced6dcb2f0080009323f472a3a2040008e78e2011dd6b2000a515cb59d8,
    0x28a8008c7572b5266151000900bb2c87221ec0007a862c129a0732000930521faf9706c00096,
    0x163b69837e00008477b6117837040008502cd97260b0d00086cb16df132e9c8008fd2c19ba36,
    0xc780075067c77bf7c470008a3d8abb10f928800a0c90f0efbe478000a5e3df0f3df71e0008c9,
    0x5d364a33d60008565f5e9cc5fa000091a88052ec884c00088ad034a8551e9800827c96163b1c,
    0x80007fe9f9326ecba30008dbd644fc4bce40007c93aef69294860008ba4cacc3e91ee800b01b,
    0x890093ae3000760cde5aa5062a00095b07ca2ace81500093bcf699bc3246000927f9d207d4ce,
    0x0009c8bf6c382f8d50009f492d7b25d3e7000794323c9de61dd00094f79a80ef789100093880,
    0x32014fc6000918b9d966a415100090bf6ba7c510810009281dafddf6d94000b4cb3dabd68212,
    0x0086b3f2b0bb1caf8008ec2668017755b80075431045fabbf7000726ca7cfa60dc9000911dc3,
    0xc367dc00007e6052e6cac283000780554faa3c3c10006e387d919cdfe000090e03368104b100,
    0x08ad607e930ee48800914729c4f02fe2000833530db4b17180007e6c95be090ad20009fc6e4f,
    0x773a358008b1ccbddeaa9f50007e019f94afbfda000b4e188ac423a1a0009257b98a1f394a00,
    0x9c98196ab192aa0008eafc5a48da77d00083b03093d6bd23000a832c58c0dd87a00087bb5a71,
    0x63ebe0007ed00cbaeeb8ac0009acfec40253ee0000aa7b99a2dbdf740007f10035ec3f1b5000,
    0x120d56c7ce7800008e8f762af5dd778009502fc2c3e1c4b0008c964b3fbcf8f88009745f3617,
    0x0000000000000000000000008c0350a2c928e700087a35cac93333a0007a61e48fdca8bb000a]

def c520N3 : ℕ := ofChunks76 [
    0x272b00007de29c53f092ec2007fb65a59ebfbfe6806f5c58199b6c550007011355b35d720800,
    0xaf6926d7b4cf0800a6a8c910006b1200083fee069d06a2d60090a6305a1668250005735e459c,
    0xb8268006c802809dc8f08000ad50cc3c398e34000acc21d3824e11600078d9413215d904c008,
    0xcf4175d2b8f800083a403d81a54a000082fe2957d20f826007c9ebb8d11819ba00878127a6f9,
    0x0dd0007d64a6f61e7c0ee007d5fe328ea3f8f6009bfceacd94ee970008491c6613385d08006b,
    0x4edf1f5504b00083a82b5b56f101400659fc3e5b6d4ab0008e5fd4bc89525c00096608d0a452,
    0x340006cd95a768380230008b8a1ed43be9e580086b81f2c448ecfc0096a4f8d05b478e0009a7,
    0x492b23d2640008ef83b917d9bb980078abf53f08be8fc0070ba93c0daa91500084cd56c2040b,
    0x20009d27af774eb7d300080e0eb5fa9ec2050086da72b255bf87c0084454c126348f20006dcf,
    0xaedfc19aac007bd3a191afdb210007e93bd8a92f45f50089b2b5a51fb5940006fa47720a9dc8,
    0x2009c0732814bbdbd00078169b38b4d1ebc006afc1ae44a6f820006fac706e231851000876ce,
    0xabe78192800878723630bf53e8007cf6635f3eb5e3600854d8e27d2880800080e4699ab3ac61,
    0x0080df73396a73538807318c4eccc5b19800895af3e72892878009286cd2019bc1000070c5eb,
    0x15663ce000763133c07e0ceb00097b600715056a000088ec77527e3f6a8008b761824ed6e620,
    0x06e92dc7e08e29d00070629ab13095c58007d667816653598e009095b35db6e71700091dbb21,
    0xd420b4c00888fecd272de928009c77345ba5a75d0007edcaa36ccc2634008753186ace8b2bc0,
    0x78c8f8c0642256000857c788a1ac0a880076ae7b94e6f4fa0006f976102e09e100007a0c8387,
    0xa3323000744fd516de0f3d0007e4fa069bf5800e007aa68d2481f588800883a495cd024b6000,
    0xa8e81f36bae6600064ee21ad2bdc9000083b0ace674d74b40085526b2a29aaf48009e7b0efb6,
    0x744300078c1bcf785e7d08008eb1868301bbb280071bc53dfa7848b0006f9f3cdf258f500008,
    0x46a7260b2eee000892a3c6c384c5e00059f62427e27bc00009c5186f8d07e090006406aac9c2,
    0x0af80090a6647693a6dd0008706ad960cf01e800871bccbd3d9513c0075e78e3cfddb9c00090,
    0x3d0303f0c230008ccf59635a7e908008ea7e749c3fbac8007c14f4bf17a0a2a0087bedefff45,
    0xb6800703b3b56e9db718007b9dd9a10018f700082907f2e88dd4f20069965fea657301000710,
    0x8b7c68e68880076004b1752a5a5000829462152330de8006ae14b527e7b750008157c7b23a5e,
    0x9b008d33d0337116da80090ec3144059dd30007158262608ce7f000a6cab5202cbc0e0008f50,
    0x960d25ce68005fe9a51012af56000993475daa00f5e00078bddd4cb79e44c007ee07c0563d8a,
    0x000779af25f51fea50007cf2e0f518e2eb4009653342f3f15460006902b07300490d00075b57,
    0x95b9dbfb00084cc3f8128aa6640078e6671f1e928a8007f330c1f2d768fc00b2f1d89e51fd8c,
    0x00704835ba618a04000877f1432e9f1210007e761d934486a340087d8479576b3d500066e7bb,
    0x2f96dc000094c8e89e2c6a1c00073b65090be5c2800074e7c606aa97488006fa1dbb1da35e40,
    0x0811b97b5811a5ae007abc41c78398f5c007a41cbbcaf33820007d7d5edd0b10f3e0099442e2,
    0x88fbf800072731cc5aea8e1800740d9a61f982bd80084a3fb1922b235400818ceb74297085b0,
    0x70aca61da382a68007f16a610c7aa17b8084d9467162fd848007738d28425ed5980093273ae8,
    0x665a180081e11d4d2931cc5009665833b86af7900082c5099da6e3280008974fd552d56fd800,
    0x82a1448daaddf00072a02d8e9c0ed80006d07d10d0fb85b0008473e55df4a2964007ab1fd133,
    0xf574000759780148ba6390007e280d47e3e56c1007ef1d2ef49a5b84007aab896afa10c54008,
    0xb7ebe9b47fd60006b616b15b58d340008c2face88af4f08008594d3fa49c0ec0005fb3738546,
    0xde650082ab22db310f6080070af6e0415145400095775bd05fe95c0008629981aba221400097,
    0xa9cb4569dee0007fd6a3d2d32cd2a2076c9477a4b9ee40008c4bce72eaa4908007e6657c3536,
    0xfb8007f4d8a153560e678074b2e83c1bb0d10008200238793bc7be00738edcd1633ef3000638,
    0x60d5628c4dc0072d708c2313f0e80092b220cbfd675100089b5fc5fdc53ab8007883600dd0ac,
    0x6000a109075c04e8700007ae1f1c4872dcf4005cd135d98ddf9200097094319ebeeff0007b16,
    0xc6f5bdcc200082a5b6add978c5400722ea861894e0100082816a1907ed4f400a318f1e67761b,
    0x80089eae7e96b218f80085f772261b883f0007318a5a48f6d52000973ebd55a816a0000758e3,
    0x118f52c60008a15765c1001de00074087a953375b20008d9fd8582ce4db80073e13d0696c729,
    0x0084b5abe60cdb64400782adf947d361a80091ba5902123369000702f5679666c65000918bed,
    0x4a2f9398006c5027a97514e7000844a263c29aae9000792a341f3b5423800833ab9ede966b16,
    0x06d68b429e207bb0006831744f66f3010005d36c912621f3c0008d367072756ab28007dd2622,
    0xd0742b8005d5e03e01dc250000701cf729a86d4d00086b69cb47bbf0e000669316af53c17b00,
    0x9294a40ae980860006e5708000feb0d0008db564af59099e80079bff5402403fc40075cc84b1,
    0x74bd28009f42394036860f000725c0e0aeffad100085cc31c6aae00e400865fb195392134000,
    0xb08268b79e7afc00610c845fc43c570007cbd00ff24a90b00062cd3830d679740007d5f258c4,
    0xfb198007461dfc2a494bd8009026e578f2de460008ae8c63fdb975a800793e8ee5dcffffc007,
    0x699692ab0817800758fc07061e2618007374256a7d8aba00096ca07076f6e050008ce5d269aa,
    0x9c800086317e11b6a0df4009679963d2495f600088c7c840bf8cdc80082edc3b54822614008a,
    0x942a0135c96000862d66e037f201c0086f01bef11b98f00089f76bb7d3cee280090d3c4b6090,
    0x369007e15c04a6e1098a0089789c5f596c358006e512bea22a1f00006ff6be59c51d4800094a,
    0xd2000179648007bcc1a03a6ff358008c0d7f1bac02340008bdc45aac373fb000812f2016beba,
    0x198096ed99489e871c00086d0fadb9b65b200091021ef10164f2000812d75bc2b84ab1007aff,
    0xd3553dd8e00089cb1806efd0110008e20c70801c8ae8009202bc08e96fb60007fe0b82934d6c,
    0x00080c83fb7a9af1e20083ab9a9c7eae61c007e4161e2814b285007b859a0918bfb280086a29,
    0xd890ec4fa0078403aa0e13b060007edfbc72f631a8e007e403cb03057f60008f17728a4817c5,
    0x007ed717f745e44e9006a29a63cfeb9af0005afaf0327844a80007ee655ff7377016007dc95f,
    0x9dbb3cc8007a5480f0dc346e0006022a8dd928e4f0006cc93d3e4945ed0007bc165ccecc2c48,
    0x081540f0abe3cc980097e415c99610720007f627b210b59975809c5f2f78c21f6d0008a74b7f,
    0x4706dd0007f2c2985f0247f20072121d982f66bd00085327c071beee000081d9383116ed09e0,
    0x85ab6d3570407e4006da2bb23378dd30007adab397fa2d8a40074011f8888081d8008e890074,
    0xd3a0080081b73c15741d398006c7ba879f4d5b70008a5fa178df4a5f00086cbfadb4b571c400,
    0xb7788ee1055b600077f18f77fe8846800957b6338340e100007ffe6caaa581f804070e523c30,
    0xff4800075225b2c01236f000870501192d03a20007857f2cd9e70c6c0074015b83c859b70008,
    0x1840d76c48a000076a969b0a2380780085c405f3792221c00943e2141d2c39200064c0a7f1ec,
    0x8bf0008d98cedaefbe498008351d6f9c27068e00b63ef5796ac86e0009e6f738364aa9900060,
    0x6f9be57262a000645db654f601160007a61800ff937ea80074c568b6dc38d00009397cc8dc56,
    0xf20007e32f4c1771bb65007c9296ed241f88c007b593a30dc942ac0089f563f5e34d18000743,
    0xa4e113c9e0c007ca9b255580dda40068e7740020388d0008d66eea55acf7c000622e6b87b8bd,
    0xe00073786ea7cd757480076afe1f083db50800975c05401e9f830007c80a3c517eac340078fb,
    0x275b1647b0008686f297369fdcc006d574eb44c29ea000731bdc541079d680067af7f4ae7afd,
    0x8007359aba685ba850006910211a8802930007da83ba68752eb40084dd5b2e50e1d080072312,
    0xe71778270007ff155c64836ab7e8a016d4006da5340008e72ffe8ef97c48008e296bad7c1474,
    0x0081c298ed36f6ba5009e66f55b8506b200086892b3e11493780073b49058f47c48000681117,
    0xbb570b20007d431f612973258006d4f7f28c464340008c4a0f9f452808000923f06416b74ee0,
    0x07823343228501940087a2428c059734000876bb580c0b6ff4008e047d5aade06200094e259e,
    0x9339c000070fa4845a83ffa0006946c7ed22617800094ba1064e03f160006ae3cc39c5d2b300,
    0x89ee3ce9243cd28006d1340d30a4f2b000772b89edf5d04d8006bce0d1692ea510008a103d98,
    0x7fa8000088ea7a840885bf800871455f0dff62d00086c00a2a04b218c008724191e23d5e9800,
    0xb984612e8fc4c00081cca6dace4d5d3007d9e72bc2d9829c0060805cbd53117a0008a3acddf7,
    0x825d00090ea74b699aaa30006ed25bba324e6a0007801a33b3aa4e94006d000e99b1a1820008,
    0x683d67f3741a000873b32248bded50006c71092694fe260008c9a52e8f3cffd000947714d021,
    0xdae000755aeb74b21bac8006e9f427acfc35900076f137571405de00082afe2e1bebd78e0092,
    0x0dc237734e700094b8191527f24f00083a44f9031fa986007917ea87ae77cbc008b7004d4666,
    0xda000810a84c37bc0591007b825996dddc6b00064c65441c0a17800095d5b32062b92b000894,
    0x8ca78eb9ea80084eb11db43b6f8c008ce622174ad3520006d4aedff6443450008941e6d85463,
    0xa0008a5f262ee7860f0008ecf096cffe7c90007e23655521d5ebc008822c0f2e750108008438,
    0xf95d35a5c0008ab334aadbf29e80084ae52040012ea40090a5625f9940920007d98be0e265a8,
    0x4007e65ec36b807303008d90c39448c2d20007fefa949e78a43bf0994972e4d69674000767e3,
    0x40cd16cc80087a76343049f2f40066a569046c89d40008c8d0f1e13591a8007b86a2ca688e8e,
    0x008a8e8f08c830b38005d4e293131420800085e4dfa02a360200054dfee4fba971a000864c0c,
    0x80f47f080095f8056bce98320007fe3f5217c66000f074e189ec6d47ef80078bcc4f3fec2a90,
    0x06424c618def4580006d20b112ed0735000857561f3fa520f00088edfdb834cc4580075d5b77,
    0x484dfa0009d9bf1131e91ed0007f71b24eded6d9b0088d204e8561fd5000601e5e61e0d8db00,
    0x88295c71f7fabd80077fd0df89c20060009559a29a35535a0008d985ab7fbed16000892dc5f3,
    0x5e8cf00086bd04d16ff89e4008c33d1b03ba89b80095e8dd50cb205a0008de693c5aaefbf800,
    0x6c74006fd863b0008a25e398126a05000880116c99576978008bac477c2c5709000721bc119f,
    0x2eb180082bf00add343fd0008aff71832fcf7a8007ec5284380c5a35007cc545786e1c736008,
    0xafac4d720de30007087098e3d92eb0006ae4dc6cb597de0007b5b821c5cc6b500072b7e01252,
    0x7cb00076c62b3c241d5600091a41b6814923a00064fbbae2db8e5f000633ad95b47f2f80007d,
    0xfe7d41eec8f8008a9897f62f841f800952319872779350007e31b650dbb601600625f787e34d,
    0xf3bca769b4c1d0471c38008bacce92d1a08080073a7ee1ada7a0580092f09f8dc03fe7000870,
    0x49f397f3d80008c49c274d7dbd28008f82cda4130551000a36e732fcb12960008003214eb6f5,
    0x8000836e81ece8411f4006e81d70ee0af5f000a48bbdb8efa6ba000736b93e24aea6f0006e44,
    0xd0d96b2658007d63d2db369da5a00916d00b7494346000904cd528dbfd040008e47f89767f4e,
    0x00091ecf964effd0c0008a08793c005ce5800719101adb894960007b1579708dab6200072dd3,
    0x7632fbe400090643fc10214740007200e53e1cda3c80085db9ac154975c800a2091113f91530,
    0x009210f9b6bbe1290008c2e5348ed689f800671306bc7f48720007c302222ac85fb8008e79e4,
    0x582cddbc008aeac363f6d8718007b95f3de9ab2c20008894735979005c8007b473d02e531300,
    0x0768b3195f4b5038007b772bfd0a3722c0081c93293551fb9e007a7c35294be2754007e1f93c,
    0x0040ec0009884cea1a803800008351e66cb76fab6007fde9bd945eac45007af1ec1dcb71de00,
    0x86d0f46fb0f5168009c313cd4b27cde0008a1ab7fec7c7cb0007cce5c412376f42006a45a411,
    0x1b5d44008fd32e65d17e0a80086a10198f4f9770007b143f34b9ff36c007fc7aa7f607a6da80,
    0x1c836d259611a0008bab1f064a27d40008378923e769f75800829ec447ec0bf02007d0eea969,
    0x8cd000082a1ee86e31343e0096508a49dc5d380007ec801430ca5e8300847771553149258009,
    0x33a47c9b8f844006743389b2ff82800089d2a0e2d79aa90007f025d632a4834c80993ee1b76d,
    0x83a0008f2ba916bdbb2e8006f6cfb9533f43200073ddd13211fdd20007fa6bac4ba2fa39807e,
    0xecdd73d957e0006e3303adb799ce00086ca416b55be030006b3864ccbdd98f0008ea064e8992,
    0x8b000809dcaebe271f83808fe54b4407211180075060bb2b54488800886816e22371d28006f2,
    0x240d65f2e0c008ea2d33ef0caf18006d9498efff5b5f000779f399c944371000889410f2cc89,
    0x380074d44f39125fe08009536ebaa8f5133000811b968073de7ca007c2cc748f9d31cc0083cf,
    0x231d8275cc008d1d3264d75c0380066a3e032297d2b00076160e86c8339700074026c4a7af35,
    0x00072d34700ee9d8d80077d6620145edab0007f4e73012b5642c808350938f5e9ec26008313f,
    0xb72ccadac0073517672191327000870b76e743470a800692e58fafc061b0008da1342e107bc2,
    0x008ab1491d8cdfeb0009c4f52a2c49be800077a06478c3d0d2000742621b1ebc730800848457,
    0x3290708c006dd0610efb7bef00078c55fcc6dd9088008205c9068f58e76005dae8075372c560,
    0x05c055d9082ece00009f90319da6ab8f0007abe6fe464c8d4400821456b88e9326800835b487,
    0x60b82200078bc6ce834c3388008a811434bd2e920007b6029229948d90008d6958bb73bf9b80,
    0x8fd5d1fb8b1d31000a7a8f15b7ed0f20007ca61967a5eb77c0064a248ee51778f000767f87df,
    0x327080008e722e8a09a78c800924ee669cba8360008ea3b2ac7220ea8009acf24bfb5d4c3000,
    0xfe82baf525dfd7a06662f5999707b10007779196cc4b17280076adf64eeb4313800879b35183,
    0x87ae000773cead2da53e90008185a9dee6a82a0009082c13a048c120008e63c2c6f46cc88007,
    0x64868b12ab037807a98dc2f3d76d880096baad2bbf60bc0008bc976c48282e6000b4a5f59e5f,
    0x9420008da868e6c2fbf78007b0433445babe58009ad3090f4d50d1000872d7ea2d6763b8007f,
    0xe667041c2f220091d0556cba26a20006b53ee4a52a9a100077e2492cade4f60005e7ba3a6669,
    0xb780071a1b484bfe16f00091dc40edb8da8700095e7d53d861aa900079c1071b493349c007d5,
    0xb95ddfacfa00099dd80832fdb940007863e180186d85c006c584a8214a82900086d88a87b567,
    0x80006d305719da64f900061e5a81d58055300076dd93e2f0d4d100091be8a7356174c0009733,
    0x7006443db0007d1df3b991c013a005d65f5004b8e9e0007b59d57cfa73164008f3dabe90ddda,
    0xf8072c6766a94a5bb0006743ece6909f6000072dc20005fb04b0008a4ba5c5f0871580096905,
    0xd7534c09920905f08cdc1cb6000060f6014fefa8a50008e358cc4ddadf5800808d6ce33220ab,
    0x0078d40863027ceec00755562377a3dc2000658c43e4ff9e66000891dd08d867d6280080165a,
    0xf919fde00099aada647eb5d900086f62169d45f25000907366a1428c5b00090b3c76a1908af0,
    0x08ec62a25a69012000755291181ada108007f6cca25b6e7d188090335a3353e1c300077c0e02,
    0x28558e000787bbff633d174c0077d97bc327df2a8008e0d504ef365a68009618fca917e4e300,
    0x942d32f833aa750007a94ad3a9124c58008d45ef5d6cbfaa80085a52d4d4babc1c00960af5d4,
    0xb697300058994af07d145a00072ecd9f17c66890008a017098aa0443000730f56c3dea678800,
    0xbfa2a5609393600063faaaaabb0eb90008e20739ae898e2000822e43dff359e86006a833efc3,
    0x103b4007791ac44b93c0a800881759058398dc000790ebf8e2c0dffc009c3b94953945640006,
    0x7bb6ae4882a10007ceb45fbf103fc000802ce6534535c09e081655a2e6228ce1007b39012712,
    0x7aac006e61f81c917daa0009ca593dac5edbb0009de6034f56c7c60004f3bf63838e99000078,
    0xa1fe74ea360000946cf5d1a88a3f0008ff8fef58eeaac000960120ad8bc9170007d4f6e7ca53,
    0x9c0008242aa9a3b97d70008c100a9504f0558006ddd5916d7601f000a93c7d78a31c7e0006ed,
    0xa56c1cce81d0099f66a1e34f3d100070233d188d219580079ede9b333a0ea000798a670cb4b9,
    0xb8006ddb0c8a0cecfd00079ff4cf3d30195400835aa7f522a3b9e008425c1b866974ec00816e,
    0xb3b9ef8e50006ec714425b7983000694d8cbfd559700005e8073eb8cf9c6000840f4351a3a35,
    0x80086f05644585992c0067c13b78666dad000999367adf1e34b00081b70b9eda70dc5006ebdb,
    0x498ba892c008d2d82a036409c8007f27bbd09366b9f80702ae0480d49698008e64890a0cd605,
    0x00771c0bb71c87df8007f0a2f8af9b66e20084ec72ed85a4590007c7813b09d12284007be96f,
    0x25df769000715dcf75b302298008dfa8c026728860009e201c9e3e84650007dc83b9e0ae5a8e,
    0x074c6a271dcbef980084aed04c9d2ab28009542a942c52e2b0007a3d83431209d3c006c14d14,
    0x426cac00087610b56c00696c006e11b36c0368730007dddb702e2e2cda008e1da1fe82dd8e00,
    0x8918290ed2dd758007580073493958980083126efb2dc8f4a0086f39507ba1180c009862e392,
    0x2d0b60007f724e147f06bbf007f4e6d53d20d61a80658cb39bc70b570006edef16e495a63000,
    0x0cf3a712e3ac90008ce17f5a4a6a4780082bc0844d76b806006c89d767c06d2b000a5d744783,
    0xe61c0009f093d943971340007ef28c6c09ebacd00741b9f72c029d58007fd4129011f2b1b609,
    0x767d721e399310088a5e61e40800100096c67d144bed7a0005f49d268bf0f6000087b94c89c9,
    0xf66c0077f0293a111f82000841ade1d42e78e8007d9a7618435d81e007fb6d4a2b8b4300c07f,
    0x5b54f510fb680091b59b27bc05a40008333cd35bb3a0b60071192a6cc765cc800846ea359a80,
    0x2900089e041da7a329f00080e06e6615dd82e006e8bcd327d46a2000908a3d19c9bad2000882,
    0xe7bcd59bcf20082a1de5e89cdc70006f716d27b2df34000990d93559462530006ee973ba51de,
    0x80006f8059a96810230008732f0eaa0f6d58007ca219d4efa154e0081f25c8efd792d700826d,
    0x64fe9e36b000941ba39a3c0b740007f1533520b386ff007586e309250fc40006b5ec3a82f29a,
    0x00080c4d34f87acef7807c1a07d28ee4cb000861affd13e5075400851c989814ed4a40099663,
    0x4574e1d1000855ad821121dbe400776e5f5decaefc00084a3918cb9d58e4006986faeb4d5ecd,
    0x0091836470598e6b0007a2128510162020007c7e8511a74244400756a059e85bc388007a76dc,
    0x50bc0c720097f2fb522cb8ea0006c2576b5ddc80800086eb7c21320e2c4008358ea930205782,
    0x069dc755c6c59d00008ded4c3bb825d380084f1dae5321f25800724cabdf035c57000828cc80,
    0x2248350008127fb7c240b1a70091f050c5c1de150007d78b4cb97ebff2009b0aad9b305dda00,
    0x8d3ee24f73040e8007f56ab8dbd90997006e65b1d0d3343900062a42ae08f67b90007aae6562,
    0x807778006db4fa037e099a000966243d3398dde0009601ce8a301a9300061e4a4ec6d3276000,
    0xde4b81d4afb6f0008edbcd75c5c3c18007dff2e4ea12ecbc006f7636ba6c7cea000878e8b08c,
    0x049a000856ec5552db12cc007061a24eb0546000079d0564b0960774007269c7e40ccea48008,
    0xb82ec0217546800a46025221a7794000978459e21e72f300055ee48721ff036000670e07334c,
    0x2c740088ddc6073a98de0009796ca111d259a0007fe144c2d567928208c2d79cd2fa11280078,
    0xb4dfb1e38150007a62e751db0bcdc009153c95600b3480006336bcb229fc4f0008342352558c,
    0xf00006d525889bf22b800090f282e7dccb060007df8b776dd8519a007cc570dc3c1604000698,
    0x9a21eca400000929c86da84acdc00095bd185130fde200097cf56d78ee78e0007713e4469d82,
    0x200090f087a90f0f6000088e21da784f4c500076b09abb2d18378008343b5fd1bd805c009a95,
    0xd290d25306008a7a09e2705e93000845f2c668c5842400608f6eaf7c9a3e000a69e3debd8172,
    0x0008818023a1c85820007388dd63c7dd618008bdd300a1771a4800671328fa5c67770007f205,
    0x8443c327000914ed08a4a4d6d0006712fd6929b4ad000948b168197553c0008e8554f00c6929,
    0x007c1e6101ea6734600925714b29e30c40009507751b494b6c0006b602d1dbe95db00073e367,
    0x773a6124007e801c1101f2735007d955fe934ee86800723c1b134ae55e0008c271427d464158,
    0x07f7a3b21731b5210078b52966ee6c584009898e4c20320390008bbede5114a6368007b8216b,
    0x353b5680088394710177d650005655dc123b4f220007abc3e7c993be60009a30287962749900,
    0x6f5d0cfb5b8531000a72afff54cdb8a00096bbaa693f430c0007005fc97a69c5600071d106eb,
    0xb163ac00938e48124b9e620008244d7e5fae40da007e7e5ed61369cec00907f1442c9c772000,
    0xa9ff5b0c7bb760008afa68f62242a28007879348b33b332c007260ceb903877e000836b35f8c,
    0xe1a60007cf44a250be16be008059d462d924c1e8066ecf404997a4f00091861ef912abb90005,
    0xc48c627d89ca0008b2b9074336d4d80081a18aee299b4f300830586dc00a4c1e0094b069a11e,
    0x0050006b9790c3e5799f0008aef82091c983400092d762c11956d00007cdd99293df00b20072,
    0xbb80dda4386000935fee0d6b9c13000766d51a692e38100069f56494f04a260008a3511bc64e,
    0x3e60093e63954a9a22900072f69bd777a20780094d26c26db29620006389286a70c0450005bf,
    0x6544a6d4804007594a2dce54e2d0008d0330b997c1a200086d82b4f497ccf8007d4ed8e0a350,
    0x94008790972d2cc0d20008204b2ca042a6a800ad0eae62c4bc36000795bc311fb5b770007b20,
    0xa719d1dfe3807a93eeb7a303d5800960f610e3f31050007c21a022951d336007a8c4603da994,
    0x000692507889da6cb0009f5853a818bd340007917a8d7a27da7c008608097a6d0b5600080527,
    0x3195fe3600060035e8c0d22ab00095b6c3173fa4600007a5fe12efb99d7c00a02e418c31c816,
    0x008318ad822420fee0073b4f205526b2f0006ba918649077e30008e05dcca4447df0008f4c1b,
    0xd54e97f00099243caf065c240007c69fadb3e1fc500073be7783870262800924d2660eb07580,
    0x0755069cf3d042d8006e4e956cc68af2000872e278d83fbfc8008f2a2a0568649a0006df662b,
    0x6a7e0a000765f17c7b8f30680096cd33d9f4f09b00083fbbce2d878000008c884dffa70b3900,
    0x9b8131fd2f333b00081ecea3031616f400764315fa52fd3f00090c3e62ad17a3d00064917559,
    0xaa1d6a008ecee4deb332418006e721db1c652f900098cafa86e833920007bfe4dd5add554800,
    0x41527a55d4f6d4008139eeb4c958378006b851558ce19af00074362f0bddca8500083f2e3baf,
    0x00000000000000000000000078a1d9e31d997fc008896f65722528780084706c12d8b11b8008]

def c520N4 : ℕ := ofChunks76 [
    0x272b00008de29c53f092ec0008fb65a59ebfbfe8007f5c58199b6c550008011355b35d720800,
    0xaf6926d7b4cf0000b6a8c910006b1200093fee069d06a2d000a0a6305a1668240006735e459c,
    0xb8260007c802809dc8f08000bd50cc3c398e34000bcc21d3824e11600088d9413215d9050009,
    0xcf4175d2b8f800093a403d81a54a000092fe2957d20f820008c9ebb8d11819b800978127a6f9,
    0x0dc0008d64a6f61e7c0f0008d5fe328ea3f8f800abfceacd94ee980009491c6613385d00007b,
    0x4edf1f5504c00093a82b5b56f101000759fc3e5b6d4ab0009e5fd4bc89525c000a6608d0a452,
    0x340007cd95a768380230009b8a1ed43be9e600096b81f2c448ed0000a6a4f8d05b478e000aa7,
    0x492b23d2640009ef83b917d9bba00088abf53f08be9000080ba93c0daa91500094cd56c2040b,
    0x2000ad27af774eb7d400090e0eb5fa9ec2000096da72b255bf8800094454c126348f20007dcf,
    0xaedfc19ab0008bd3a191afdb210008e93bd8a92f45f80099b2b5a51fb5940007fa47720a9dc8,
    0x000ac0732814bbdbc00088169b38b4d1ec0007afc1ae44a6f820007fac706e231851000976ce,
    0xabe78192800978723630bf53e0008cf6635f3eb5e3800954d8e27d2880800090e4699ab3ac61,
    0x0090df73396a73540008318c4eccc5b19800995af3e7289288000a286cd2019bc1000080c5eb,
    0x15663ce000863133c07e0ceb000a7b600715056a000098ec77527e3f6a0009b761824ed6e620,
    0x07e92dc7e08e29d00080629ab13095c58008d667816653599000a095b35db6e718000a1dbb21,
    0xd420b5000988fecd272de92000ac77345ba5a75c0008edcaa36ccc2630009753186ace8b2c00,
    0x88c8f8c0642256000957c788a1ac0a800086ae7b94e6f4fa0007f976102e09e100008a0c8387,
    0xa3324000844fd516de0f3d0008e4fa069bf58010008aa68d2481f588800983a495cd024b6000,
    0xa8e81f36bae6600074ee21ad2bdc9000093b0ace674d74b00095526b2a29aaf4000ae7b0efb6,
    0x744300088c1bcf785e7d08009eb1868301bbb200081bc53dfa7848b0007f9f3cdf258f500009,
    0x46a7260b2eee000992a3c6c384c5e00069f62427e27bc0000ac5186f8d07e080007406aac9c2,
    0x0b0000a0a6647693a6dc0009706ad960cf01e000971bccbd3d951400085e78e3cfddb9c000a0,
    0x3d0303f0c230009ccf59635a7e900009ea7e749c3fbac0008c14f4bf17a0a280097bedefff45,
    0xb6000803b3b56e9db718008b9dd9a10018f700092907f2e88dd4f00079965fea657301000810,
    0x8b7c68e68800086004b1752a5a5000929462152330de0007ae14b527e7b750009157c7b23a5e,
    0x98009d33d0337116da000a0ec3144059dd40008158262608ce7f000b6cab5202cbc0e0009f50,
    0x960d25ce68006fe9a51012af56000a93475daa00f5e00088bddd4cb79e450008ee07c0563d8a,
    0x000879af25f51fea50008cf2e0f518e2eb000a653342f3f15460007902b07300490d00085b57,
    0x95b9dbfb00094cc3f8128aa6600088e6671f1e928a8008f330c1f2d7690000c2f1d89e51fd8c,
    0x00804835ba618a04000977f1432e9f1210008e761d934486a300097d8479576b3d500076e7bb,
    0x2f96dc0000a4c8e89e2c6a1c00083b65090be5c2800084e7c606aa97488007fa1dbb1da35e40,
    0x0911b97b5811a5b0008abc41c78398f60008a41cbbcaf33820008d7d5edd0b10f4000a9442e2,
    0x88fbf800082731cc5aea8e1800840d9a61f982bd80094a3fb1922b235000918ceb7429708600,
    0x80aca61da382a68008f16a610c7aa1780094d9467162fd840008738d28425ed59800a3273ae8,
    0x665a180091e11d4d2931cc000a665833b86af7800092c5099da6e3280009974fd552d56fe000,
    0x82a1448daaddf00082a02d8e9c0ed80007d07d10d0fb85b0009473e55df4a2960008ab1fd133,
    0xf574000859780148ba6390008e280d47e3e56c0008ef1d2ef49a5b80008aab896afa10c50009,
    0xb7ebe9b47fd60007b616b15b58d340009c2face88af4f00009594d3fa49c0ec0006fb3738546,
    0xde680092ab22db310f6000080af6e04151454000a5775bd05fe95c0009629981aba2214000a7,
    0xa9cb4569dee0008fd6a3d2d32cd280086c9477a4b9ee40009c4bce72eaa4900008e6657c3536,
    0xfb8008f4d8a153560e680084b2e83c1bb0d10009200238793bc7c000838edcd1633ef3000738,
    0x60d5628c4e00082d708c2313f0e800a2b220cbfd675000099b5fc5fdc53ac0008883600dd0ac,
    0x6000b109075c04e8700008ae1f1c4872dcf0006cd135d98ddf92000a7094319ebef000008b16,
    0xc6f5bdcc200092a5b6add978c5000822ea861894e0100092816a1907ed4f000b318f1e67761b,
    0x80099eae7e96b219000095f772261b883f0008318a5a48f6d52000a73ebd55a816a0000858e3,
    0x118f52c60009a15765c1001de00084087a953375b20009d9fd8582ce4dc00083e13d0696c729,
    0x0094b5abe60cdb64000882adf947d361a800a1ba5902123368000802f5679666c65000a18bed,
    0x4a2f9398007c5027a97514e7000944a263c29aae9000892a341f3b5423800933ab9ede966b10,
    0x07d68b429e207bb0007831744f66f3010006d36c912621f3c0009d367072756ab20008dd2622,
    0xd0742b8006d5e03e01dc250000801cf729a86d4d00096b69cb47bbf0e000769316af53c17b00,
    0xa294a40ae980860007e5708000feb0d0009db564af59099e00089bff5402403fc00085cc84b1,
    0x74bd2800af423940368610000825c0e0aeffad100095cc31c6aae00e000965fb195392134000,
    0xb08268b79e7b0000710c845fc43c570008cbd00ff24a90b00072cd3830d679740008d5f258c4,
    0xfb1a0008461dfc2a494bd800a026e578f2de460009ae8c63fdb975a000893e8ee5dd00000008,
    0x699692ab0818000858fc07061e2618008374256a7d8aba000a6ca07076f6e040009ce5d269aa,
    0x9c800096317e11b6a0df000a679963d2495f600098c7c840bf8cdc00092edc3b54822610009a,
    0x942a0135c96000962d66e037f20200096f01bef11b98f00099f76bb7d3cee2000a0d3c4b6090,
    0x370008e15c04a6e109880099789c5f596c360007e512bea22a1f00007ff6be59c51d48000a4a,
    0xd2000179648008bcc1a03a6ff358009c0d7f1bac02340009bdc45aac373fb000912f2016beba,
    0x1800a6ed99489e871c00096d0fadb9b65b2000a1021ef10164f2000912d75bc2b84ab0008aff,
    0xd3553dd8e00099cb1806efd0110009e20c70801c8ae000a202bc08e96fb60008fe0b82934d6c,
    0x00090c83fb7a9af1e00093ab9a9c7eae620008e4161e2814b288008b859a0918bfb280096a29,
    0xd890ec4f80088403aa0e13b060008edfbc72f631a90008e403cb03057f60009f17728a4817c5,
    0x008ed717f745e44e8007a29a63cfeb9af0006afaf0327844a80008ee655ff7377018008dc95f,
    0x9dbb3cc0008a5480f0dc346e0007022a8dd928e4f0007cc93d3e4945ed0008bc165ccecc2c48,
    0x091540f0abe3cca000a7e415c99610720008f627b210b5997800ac5f2f78c21f6c0009a74b7f,
    0x4706dd0008f2c2985f0247f00082121d982f66bd00095327c071beee000091d9383116ed0a00,
    0x95ab6d3570407e0007da2bb23378dd30008adab397fa2d8a00084011f8888081d8009e890074,
    0xd3a0080091b73c15741d3a0007c7ba879f4d5b70009a5fa178df4a5f00096cbfadb4b571c000,
    0xb7788ee1055b600087f18f77fe8846800a57b6338340e100008ffe6caaa581f800080e523c30,
    0xff4800085225b2c01236f000970501192d03a20008857f2cd9e70c700084015b83c859b70009,
    0x1840d76c48a000086a969b0a2380780095c405f3792222000a43e2141d2c39200074c0a7f1ec,
    0x8c00009d98cedaefbe4a0009351d6f9c27069000c63ef5796ac870000ae6f738364aa9800070,
    0x6f9be57262a000745db654f601160008a61800ff937ea80084c568b6dc38d0000a397cc8dc56,
    0xf20008e32f4c1771bb68008c9296ed241f890008b593a30dc942b00099f563f5e34d18000843,
    0xa4e113c9e10008ca9b255580dda00078e7740020388d0009d66eea55acf7c000722e6b87b8bd,
    0xe00083786ea7cd757480086afe1f083db50800a75c05401e9f840008c80a3c517eac300088fb,
    0x275b1647b0009686f297369fdd0007d574eb44c29ea000831bdc541079d680077af7f4ae7afd,
    0x0008359aba685ba850007910211a8802930008da83ba68752eb00094dd5b2e50e1d000082312,
    0xe71778270008ff155c64836ab800b016d4006da5340009e72ffe8ef97c40009e296bad7c1474,
    0x0091c298ed36f6ba000ae66f55b8506b200096892b3e11493800083b49058f47c48000781117,
    0xbb570b20008d431f612973258007d4f7f28c464340009c4a0f9f452808000a23f06416b74ee0,
    0x08823343228501900097a2428c059734000976bb580c0b6ff0009e047d5aade062000a4e259e,
    0x9339c000080fa4845a83ffa0007946c7ed226178000a4ba1064e03f160007ae3cc39c5d2b300,
    0x99ee3ce9243cd20007d1340d30a4f2b000872b89edf5d04d8007bce0d1692ea510009a103d98,
    0x7fa8000098ea7a840885c0000971455f0dff62d00096c00a2a04b2190009724191e23d5ea000,
    0xb984612e8fc4c00091cca6dace4d5d0008d9e72bc2d982a00070805cbd53117a0009a3acddf7,
    0x825c000a0ea74b699aaa40007ed25bba324e6a0008801a33b3aa4e90007d000e99b1a1820009,
    0x683d67f3741a000973b32248bded50007c71092694fe260009c9a52e8f3cffd000a47714d021,
    0xdae000855aeb74b21bac8007e9f427acfc35900086f137571405de00092afe2e1bebd79000a2,
    0x0dc237734e7000a4b8191527f25000093a44f9031fa980008917ea87ae77cc0009b7004d4666,
    0xda000910a84c37bc0590008b825996dddc6b00074c65441c0a178000a5d5b32062b92c000994,
    0x8ca78eb9ea00094eb11db43b6f90009ce622174ad3520007d4aedff6443450009941e6d85463,
    0xa0009a5f262ee7860f0009ecf096cffe7c90008e23655521d5ec0009822c0f2e750100009438,
    0xf95d35a5c0009ab334aadbf29e00094ae52040012ea000a0a5625f9940920008d98be0e265a8,
    0x0008e65ec36b807300009d90c39448c2d20008fefa949e78a43800a94972e4d69674000867e3,
    0x40cd16cc00097a76343049f2f00076a569046c89d40009c8d0f1e13591a0008b86a2ca688e8e,
    0x009a8e8f08c830b40006d4e293131420800095e4dfa02a360200064dfee4fba971a000964c0c,
    0x80f47f0800a5f8056bce98320008fe3f5217c660000084e189ec6d47ef80088bcc4f3fec2a90,
    0x07424c618def4580007d20b112ed0735000957561f3fa520f00098edfdb834cc4600085d5b77,
    0x484dfa000ad9bf1131e91ec0008f71b24eded6d980098d204e8561fd5000701e5e61e0d8db00,
    0x98295c71f7fabe00087fd0df89c2006000a559a29a35535a0009d985ab7fbed16000992dc5f3,
    0x5e8cf00096bd04d16ff89e0009c33d1b03ba89c000a5e8dd50cb205a0009de693c5aaefc0000,
    0x6c74006fd863b0009a25e398126a05000980116c99576980009bac477c2c5709000821bc119f,
    0x2eb180092bf00add343fd0009aff71832fcf7a0008ec5284380c5a38008cc545786e1c738009,
    0xafac4d720de30008087098e3d92eb0007ae4dc6cb597de0008b5b821c5cc6b500082b7e01252,
    0x7cb00086c62b3c241d56000a1a41b6814923a00074fbbae2db8e5f000733ad95b47f2f80008d,
    0xfe7d41eec900009a9897f62f8420000a52319872779340008e31b650dbb601800725f787e34d,
    0xf4000869b4c1d0471c38009bacce92d1a08000083a7ee1ada7a05800a2f09f8dc03fe8000970,
    0x49f397f3d80009c49c274d7dbd20009f82cda4130551000b36e732fcb12960009003214eb6f5,
    0x8000936e81ece8411f0007e81d70ee0af5f000b48bbdb8efa6ba000836b93e24aea6f0007e44,
    0xd0d96b2658008d63d2db369da5800a16d00b7494346000a04cd528dbfd040009e47f89767f4e,
    0x000a1ecf964effd0c0009a08793c005ce6000819101adb894960008b1579708dab6200082dd3,
    0x7632fbe4000a0643fc10214740008200e53e1cda3c80095db9ac154975c000b2091113f91530,
    0x00a210f9b6bbe1280009c2e5348ed68a0000771306bc7f48720008c302222ac85fb8009e79e4,
    0x582cddc0009aeac363f6d8720008b95f3de9ab2c20009894735979005c0008b473d02e531300,
    0x0868b3195f4b5038008b772bfd0a372300091c93293551fba0008a7c35294be2750008e1f93c,
    0x0040ec000a884cea1a803800009351e66cb76fab0008fde9bd945eac48008af1ec1dcb71de00,
    0x96d0f46fb0f516000ac313cd4b27cde0009a1ab7fec7c7cb0008cce5c412376f40007a45a411,
    0x1b5d40009fd32e65d17e0a00096a10198f4f9770008b143f34b9ff370008fc7aa7f607a6d800,
    0x1c836d259611a0009bab1f064a27d40009378923e769f76000929ec447ec0bf00008d0eea969,
    0x8cd000092a1ee86e31344000a6508a49dc5d380008ec801430ca5e800094777155314926000a,
    0x33a47c9b8f840007743389b2ff82800099d2a0e2d79aa90008f025d632a4835000a93ee1b76d,
    0x83a0009f2ba916bdbb2e0007f6cfb9533f43200083ddd13211fdd20008fa6bac4ba2fa38008e,
    0xecdd73d957e0007e3303adb799ce00096ca416b55be030007b3864ccbdd98f0009ea064e8992,
    0x8b000909dcaebe271f80009fe54b4407211200085060bb2b54488800986816e22371d20007f2,
    0x240d65f2e10009ea2d33ef0caf20007d9498efff5b5f000879f399c944371000989410f2cc89,
    0x380084d44f39125fe0800a536ebaa8f5134000911b968073de7d0008c2cc748f9d31d00093cf,
    0x231d8275d0009d1d3264d75c0400076a3e032297d2b00086160e86c8339700084026c4a7af35,
    0x00082d34700ee9d8d80087d6620145edab0008f4e73012b56430009350938f5e9ec20009313f,
    0xb72ccadb00083517672191327000970b76e743470a000792e58fafc061b0009da1342e107bc2,
    0x009ab1491d8cdfeb000ac4f52a2c49be800087a06478c3d0d2000842621b1ebc730800948457,
    0x32907090007dd0610efb7bef00088c55fcc6dd9088009205c9068f58e70006dae8075372c560,
    0x06c055d9082ece0000af90319da6ab900008abe6fe464c8d4000921456b88e9326000935b487,
    0x60b82200088bc6ce834c3388009a811434bd2e920008b6029229948d90009d6958bb73bf9c00,
    0x9fd5d1fb8b1d31000b7a8f15b7ed0f20008ca61967a5eb7800074a248ee51778f000867f87df,
    0x327080009e722e8a09a78c000a24ee669cba8360009ea3b2ac7220ea000aacf24bfb5d4c4000,
    0xfe82baf525dfd8007662f5999707b10008779196cc4b17280086adf64eeb4313800979b35183,
    0x87b0000873cead2da53e90009185a9dee6a82a000a082c13a048c120009e63c2c6f46cc80008,
    0x64868b12ab038008a98dc2f3d76d8800a6baad2bbf60bc0009bc976c48282e6000c4a5f59e5f,
    0x9420009da868e6c2fbf80008b0433445babe5800aad3090f4d50d0000972d7ea2d6763c0008f,
    0xe667041c2f2000a1d0556cba26a20007b53ee4a52a9a100087e2492cade4f60006e7ba3a6669,
    0xb800081a1b484bfe16f000a1dc40edb8da88000a5e7d53d861aa800089c1071b49334a0008d5,
    0xb95ddfacfa000a9dd80832fdb940008863e180186d860007c584a8214a82900096d88a87b567,
    0x80007d305719da64f900071e5a81d58055300086dd93e2f0d4d1000a1be8a7356174c000a733,
    0x7006443dc0008d1df3b991c0138006d65f5004b8e9e0008b59d57cfa73160009f3dabe90ddda,
    0x00082c6766a94a5bb0007743ece6909f6000082dc20005fb04b0009a4ba5c5f08716000a6905,
    0xd7534c0a000a05f08cdc1cb6000070f6014fefa8a50009e358cc4ddadf6000908d6ce33220ac,
    0x0088d40863027cef000855562377a3dc2000758c43e4ff9e66000991dd08d867d6200090165a,
    0xf919fde000a9aada647eb5d800096f62169d45f25000a07366a1428c5c000a0b3c76a1908b00,
    0x09ec62a25a69012000855291181ada108008f6cca25b6e7d1800a0335a3353e1c400087c0e02,
    0x28558e000887bbff633d17500087d97bc327df2a8009e0d504ef365a6000a618fca917e4e400,
    0xa42d32f833aa740008a94ad3a9124c58009d45ef5d6cbfaa00095a52d4d4babc2000a60af5d4,
    0xb697300068994af07d145a00082ecd9f17c66890009a017098aa0443000830f56c3dea678800,
    0xbfa2a5609393600073faaaaabb0eb90009e20739ae898e2000922e43dff359e80007a833efc3,
    0x103b0008791ac44b93c0a800981759058398dc000890ebf8e2c0e00000ac3b94953945640007,
    0x7bb6ae4882a10008ceb45fbf103fc000902ce6534535c100091655a2e6228ce0008b39012712,
    0x7ab0007e61f81c917daa000aca593dac5edbc000ade6034f56c7c60005f3bf63838e99000088,
    0xa1fe74ea360000a46cf5d1a88a400009ff8fef58eeaac000a60120ad8bc9180008d4f6e7ca53,
    0x9c0009242aa9a3b97d70009c100a9504f0560007ddd5916d7601f000b93c7d78a31c7e0007ed,
    0xa56c1cce82000a9f66a1e34f3d000080233d188d219580089ede9b333a0ea000898a670cb4b9,
    0xc0007ddb0c8a0cecfd00089ff4cf3d30195000935aa7f522a3ba0009425c1b866974f000916e,
    0xb3b9ef8e50007ec714425b7983000794d8cbfd559700006e8073eb8cf9c6000940f4351a3a35,
    0x00096f0564458599300077c13b78666dad000a99367adf1e34c00091b70b9eda70dc0007ebdb,
    0x498ba8930009d2d82a036409c0008f27bbd09366ba000802ae0480d49698009e64890a0cd606,
    0x00871c0bb71c87df8008f0a2f8af9b66e00094ec72ed85a4590008c7813b09d12280008be96f,
    0x25df769000815dcf75b302298009dfa8c02672886000ae201c9e3e84640008dc83b9e0ae5a90,
    0x084c6a271dcbef980094aed04c9d2ab2000a542a942c52e2c0008a3d83431209d40007c14d14,
    0x426cac00097610b56c006970007e11b36c0368730008dddb702e2e2cd8009e1da1fe82dd8e00,
    0x9918290ed2dd760008580073493958980093126efb2dc8f500096f39507ba1181000a862e392,
    0x2d0b60008f724e147f06bc0008f4e6d53d20d61800758cb39bc70b570007edef16e495a63000,
    0x0cf3a712e3ac80009ce17f5a4a6a4800092bc0844d76b800007c89d767c06d2b000b5d744783,
    0xe61c000af093d943971340008ef28c6c09ebad000841b9f72c029d58008fd4129011f2b1800a,
    0x767d721e399300098a5e61e408001000a6c67d144bed7a0006f49d268bf0f6000097b94c89c9,
    0xf6700087f0293a111f82000941ade1d42e78e0008d9a7618435d820008fb6d4a2b8b4300008f,
    0x5b54f510fb6000a1b59b27bc05a40009333cd35bb3a0b00081192a6cc765cc800946ea359a80,
    0x2900099e041da7a329f00090e06e6615dd830007e8bcd327d46a2000a08a3d19c9bad2000982,
    0xe7bcd59bcf00092a1de5e89cdc70007f716d27b2df34000a90d93559462540007ee973ba51de,
    0x80007f8059a96810230009732f0eaa0f6d60008ca219d4efa15500091f25c8efd792d000926d,
    0x64fe9e36c000a41ba39a3c0b740008f1533520b38700008586e309250fc40007b5ec3a82f29a,
    0x00090c4d34f87acef0008c1a07d28ee4cb000961affd13e5075000951c989814ed4a000a9663,
    0x4574e1d1000955ad821121dbe000876e5f5decaefc00094a3918cb9d58e0007986faeb4d5ecd,
    0x00a1836470598e6c0008a2128510162020008c7e8511a74244000856a059e85bc388008a76dc,
    0x50bc0c7000a7f2fb522cb8ea0007c2576b5ddc80800096eb7c21320e2c0009358ea930205780,
    0x079dc755c6c59d00009ded4c3bb825d400094f1dae5321f26000824cabdf035c57000928cc80,
    0x2248350009127fb7c240b1a000a1f050c5c1de140008d78b4cb97ebff000ab0aad9b305dda00,
    0x9d3ee24f73040e0008f56ab8dbd90998007e65b1d0d3343900072a42ae08f67b90008aae6562,
    0x807780007db4fa037e099a000a66243d3398dde000a601ce8a301a9400071e4a4ec6d3276000,
    0xde4b81d4afb6f0009edbcd75c5c3c20008dff2e4ea12ecc0007f7636ba6c7cea000978e8b08c,
    0x049a000956ec5552db12d0008061a24eb0546000089d0564b0960770008269c7e40ccea48009,
    0xb82ec0217546800b46025221a7794000a78459e21e72f400065ee48721ff036000770e07334c,
    0x2c700098ddc6073a98de000a796ca111d259a0008fe144c2d567928009c2d79cd2fa11200088,
    0xb4dfb1e38150008a62e751db0bce000a153c95600b3480007336bcb229fc4f0009342352558c,
    0xf00007d525889bf22b8000a0f282e7dccb060008df8b776dd85198008cc570dc3c1604000798,
    0x9a21eca400000a29c86da84acdc000a5bd185130fde2000a7cf56d78ee78e0008713e4469d82,
    0x2000a0f087a90f0f6000098e21da784f4c500086b09abb2d18378009343b5fd1bd806000aa95,
    0xd290d25308009a7a09e2705e93000945f2c668c5842000708f6eaf7c9a3e000b69e3debd8172,
    0x0009818023a1c85820008388dd63c7dd618009bdd300a1771a4000771328fa5c67770008f205,
    0x8443c327000a14ed08a4a4d6c0007712fd6929b4ad000a48b168197553c0009e8554f00c6929,
    0x008c1e6101ea6734800a25714b29e30c4000a507751b494b6c0007b602d1dbe95db00083e367,
    0x773a6120008e801c1101f2738008d955fe934ee86800823c1b134ae55e0009c271427d464160,
    0x08f7a3b21731b5200088b52966ee6c58000a898e4c20320380009bbede5114a6360008b8216b,
    0x353b5680098394710177d650006655dc123b4f220008abc3e7c993be6000aa30287962749800,
    0x7f5d0cfb5b8531000b72afff54cdb8a000a6bbaa693f430c0008005fc97a69c5600081d106eb,
    0xb163b000a38e48124b9e620009244d7e5fae40e0008e7e5ed61369cf000a07f1442c9c772000,
    0xa9ff5b0c7bb760009afa68f62242a20008879348b33b3330008260ceb903877e000936b35f8c,
    0xe1a60008cf44a250be16c0009059d462d924c200076ecf404997a4f000a1861ef912abb80006,
    0xc48c627d89ca0009b2b9074336d4e00091a18aee299b4f000930586dc00a4c2000a4b069a11e,
    0x0050007b9790c3e5799f0009aef82091c9834000a2d762c11956d00008cdd99293df00b00082,
    0xbb80dda4386000a35fee0d6b9c14000866d51a692e38100079f56494f04a260009a3511bc64e,
    0x3e800a3e63954a9a22800082f69bd777a207800a4d26c26db29620007389286a70c0450006bf,
    0x6544a6d4800008594a2dce54e2d0009d0330b997c1a200096d82b4f497cd00008d4ed8e0a350,
    0x90009790972d2cc0d20009204b2ca042a6a000bd0eae62c4bc36000895bc311fb5b770008b20,
    0xa719d1dfe0008a93eeb7a303d5800a60f610e3f31040008c21a022951d338008a8c4603da994,
    0x000792507889da6cb000af5853a818bd340008917a8d7a27da80009608097a6d0b5600090527,
    0x3195fe3600070035e8c0d22ab000a5b6c3173fa4600008a5fe12efb99d8000b02e418c31c816,
    0x009318ad822420ff00083b4f205526b2f0007ba918649077e30009e05dcca4447df0009f4c1b,
    0xd54e97f000a9243caf065c240008c69fadb3e1fc500083be7783870262800a24d2660eb07580,
    0x0855069cf3d042d8007e4e956cc68af2000972e278d83fbfc0009f2a2a0568649a0007df662b,
    0x6a7e0a000865f17c7b8f306800a6cd33d9f4f09c00093fbbce2d878000009c884dffa70b3900,
    0xab8131fd2f333c00091ecea3031616f000864315fa52fd3f000a0c3e62ad17a3c00074917559,
    0xaa1d70009ecee4deb332420007e721db1c652f9000a8cafa86e833920008bfe4dd5add554800,
    0x41527a55d4f6d0009139eeb4c958380007b851558ce19af00084362f0bddca8500093f2e3baf,
    0x00000000000000000000000088a1d9e31d99800009896f65722528800094706c12d8b11c0009]

def c520N5 : ℕ := ofChunks76 [
    0xa42910009381f4d966dd4e0007755a697c045f180083ca7f5fa096a100077d0789d8af12d000,
    0x0bc7e438e8efa00075b37b70d1f9708006316f7673e18a20007627dc69019c3580082c942b30,
    0x8698800870f89613266968008a727914a6f8b9000729fbf311e21e380077071dfbc6a927000a,
    0x37ccfe804da33806d3cfe7de7320600091f6644478d42f0007b8d69ff427aed8008c06ac308c,
    0xfd600091c15d42f4f43e000821481f847b2618008ae4632b6c0d8980074bbd7577dabe28007f,
    0xa9ce10952e2c0072d2e01ffe4a83000777235cf2e933700082402003cc2a2d000a0ef9e9b70f,
    0xfc4006e822e6eeb56bc0006def2ad5056fb30009fde6baf3d5c5d0007ee2e92f51733c1007bc,
    0x59a2399d49a0071791976fda48c0008d2346307ee0e20008966a4c082e3918007978f7eb73a4,
    0x2a807b00e5f563b47c000916d54d6ad55ad000931f90d1e723df0007f504f73eb3acfa008208,
    0xde8c25da20007b09214fe5d676000658b4b4e1a3a3d0008a718d445eb55d80080bc63f5ba1ff,
    0x800746bea9b6ddbc98007fa9f400a57b910007f0a93364c63e3d8077a0d7520b8c5c8008932a,
    0xc5565d740008d5874114c55ab0007a7292f423dd1e8007ccd5803a241ebe0070169c4fb15165,
    0x00899f2b7fe98dd20007edf94311a3275f008660f5e219514e0007fcf2a7aa879d6020987a65,
    0x496de78000703385e17f7c1380091ce8017c559b200083c9185d6decc3200837ef5190cfd488,
    0x07d60c1e3f9b151e0077093fd738b2d880082e4a59ee67e4fe008595f6e58f0a66000972a737,
    0x659f9c800896038bda9571800076602a855ea0980005efb06c242f52800085fa9ceae0a79240,
    0x65d7b9ff4a365600083c5a0c70482e920075d431f8a8777480073310a074991288007b49cd77,
    0x35f690007367b3584d31d700077f25241cded360005d4c0b27a231f40008df08e1279e971800,
    0xb2fb5266b30fa0007e1567d350488c1008400511c734cf5c006ec0e43ba18b520006c05bbc96,
    0x866a80071f391695120fa8007de62a8363d9146006acf9515eb5cc50005ade26cda226b20007,
    0x80389b7ccca70006f4c6994ff4bba00082cf35ba670cd280090caeb33e2549a00073406f38dd,
    0x2658007f2bcf7ba024e218070933c85755caa800751d82cccf805c8007b71090c038f4300078,
    0x3e3bee2e0e480086348300402415c007ccc32a55fe63c40083a12f9105a8612008bdbeefe84e,
    0x9b8007c016e47ae40bec006dabbbbe4e1d3d000874afa646d1bde4007dc0746f8b117c4008cb,
    0x9627b634340007e05a92770242090085b794aaa2da0e4008321711a3c748900079467550c301,
    0xcc007367e481a426bb80095b59a8f287cf30008474d2071cc41f400711fe3fc4488b900073bd,
    0xbdaa294e650081f1f5cfc40bb7a007720464e9498f6800854a77b2fb47e400079d7dd58ef077,
    0x0009062e01b7a25e60008d8a7a163389d30008d4d955762f32e8007394793650c5158007f3d0,
    0xca41ef80800842df2ed13d10b40078917b6e4e9ab6c0077e5f12e8edaef0009449b6c9d17cb9,
    0x0097b719546bda02000761b8f5091ad4a0007939facb4fd1178007a90e429e25d99c008a14cf,
    0x7d4942bc006f0bd6cd09d758000809c851b13b22ae807d5199b215ec8ec0057de8cf590de2a0,
    0x05cc199120398c00006562e002e867de0008f478d6db2c8ee8006086cf109abd41000862ab45,
    0x6bcd73000770d29146345b880075932e4f030c2b8008297e19c72d1c940073c8714c3d089600,
    0x9a24d256990d4600087f7889603b2c80007ea05c3a0951295006f0ebde8480b9f0008554cebd,
    0x1ff280007ff2b96ea486b70e0874024bc3dd48d8009b6fab51bbf51500087abd5bc5eeac9000,
    0xad6fc01eaf0414008b6575b49cc8bb8007ae58772a65e1d0006c67cc6f2b0f050005be53bb4e,
    0x4485c008378db5148629ac007b9d2a2851e1d40006efaa0f2abd33000088b2c8e58ce1590007,
    0xb3ffc8f7bbd1c0085fd516ab503a300075449db7c3ea4d80081ce289b7e7a60e007a3ff4c2dd,
    0x953000682a6af6caea330007538a74b310a9c00060d1c4a1cd8c020007f3f492f2d39f10007e,
    0x5cfd3262dbd0008489a063ef3ae1c007d6d5eb0d484ed80078190e089d2988800664efa94fec,
    0xff0007b13162ebc259e80077943ba9f0fb32000903170ab6cf299000871c437c59d7d2000917,
    0x18da0b0dcd0007274823334042e800651db7b78dcf370008004f1d2f88cf3d74915cda0c53f5,
    0xd800780eeb5354ff8d000803a4193ccf3052e08555c983d52215c00763a6b91bd39af800901f,
    0x9c50c217ec0093c1357a036f41000835a6f54806b3480086447cc3e38c294007a632dcecbb96,
    0x000a482efc85666e40008b490ce26ec2aa00082a78a02b96b36800713a8a344175938007ccd8,
    0x03769bed800948025e1a27baa0008e1eb133599ead0007ac7417a1b8a144007894c11ef023bd,
    0x0089f94947f44446000709aeaf546342b00076f51d650197878008aa214a11dfc4d80075760a,
    0xa147e628009cb83ace1bf12b00076b641f725d967800863ab303cf19cc00081b437020b3513b,
    0x07dcbc685dbbb474006c91595d8a7b250008b3e75a9bfd2730007a02729a2480cf8008aa5050,
    0x5571500007c7301eaaebe23e008e915e5a30b33c0007cfa83a1475cee8009674276976a87500,
    0x76816cc5fbe98a00097251474966d5b0008390534c40ea63c0071eaa63ff43c3c8005b3b48ab,
    0x6451e00073901c6c5e36498006fbc98615a7fff000709019a6122a34000846a68f7c453d6800,
    0x9358e2276cb8180091d8ef432ec1d1000897fa7093744780008f0b8297c4421b80095e7ada05,
    0x9619b8072b3f81965bdc58006f36e1f2db9c08000835aa43acbe0df8008c273a256c076d8007,
    0x8bb4e48210c6000700e5f3292fbbd0007adbb74658e8c0000834a42bf533be9a0080cd044d9a,
    0x760000674406f22efb8300088c95cec1818e08007b46912c1d278200086cdc3b7d4c6400005e,
    0x99b569f0e210006045d22edad62400069651da251d9ee0007e37fc117481230006bf25604aee,
    0x6e8006ddb3f3a1bd4650007ba5a7470de0d84006f8544d81858820008e69f5f80317fd8008c8,
    0xb7268877c8000739cce108f9f7f0007216b0a35f6cff0007f86700bbcbac95c0851015abb073,
    0x900072838617b768a080086060a89f480240008922da2c71694800086fb63940eedfb8006f01,
    0x8a24514ce0007f5ea620928eaf4807d2ac43cb7a2aec00887e777ab7da31800636bd96559efb,
    0xe007496c85dda706d800796391509133b9800924de3039f09e900071b428bd8c977500090839,
    0x8a4d5824800719785bd9ef61f0007caf3ea5c68399c00863d8553c2cc078007cb5438810682e,
    0x007d4ab6433db39fa00874ba33b86933700077dd22dacb95118005754400dbed5700007b7e77,
    0x087157f0008a7f28fbf972f40008a5491aec5db6000087c69584459cd14008b0363ee00ac968,
    0x0873eadc61bec84c00784feeb8fa211e4009a14c17c46a09b000801109bd9a41b09f0642ec5f,
    0x1249b130075a090eab71a790007ba5c3bf739c0800088e985a6aa79d5000824b57399029fa20,
    0x8694e5c3e5926fc008080c7ed40d967d807d08e1024ea3da600835335a4d7eb42c008126c06a,
    0x7e0ce00097a2a20e064a06000a3365066d3154c000923c1ea39a5c780007394c0336c81bd800,
    0xc73200246ea738008e15995f60a00980080ddfead54bb46e0062d57428b35fca00078bab4590,
    0x6ed60008a046ddf7077e18008223d16a0d09e1a006e308b2fec414a0006c52f50bc4067e0008,
    0xccb80fb1c97c0006aff91d730e8bc00070832c216aca6c0008d0ff1092427150008df787313e,
    0x8170006f0dd3865fcbf3000855fa4f8cd7b714008d3ff1e9c6b4630008814344ef6c0a20005f,
    0xb832ccb915900088d28ae8ac517f0009c76f0edcf413f0007a13e1221a26aa40099155e28a6c,
    0x3cc007710536cac6b430006fcc9a658b585b0007ef72f39ab28ee800854a258f5eb54e800930,
    0x8e8268e494c007d2706d7f3d5d66007f77ef21970664c008483601910ba1f4007f3f0bda7ceb,
    0xe9808d8a53f4056c3c00068b56efc9c1022000776cccd8852044800779c2c0fa1b61180082ef,
    0x9948e8d130007c091e45dbc50d00074253285c8b4ae0007ea979a2c46c6e30080c03c279b0bc,
    0x0008203f2e9d48e2c4008fb509c5b0c621000614bf93a1eb3e10006f3613fabdd9500006e9c8,
    0x66e94cab2008be0cc3c678a2000089fe53d956c87e800929edeb3a9c8f60006e21af1860c5dc,
    0x008ccc734231dc5c80086799ce8d898560007e6e601489aae2f00793bea9015219b4007c633e,
    0xb0ed4b8000819d1df9e50a362006fc03f2dbb80a800077dad8fdede714800996e7734110e0a0,
    0x07857e30bc150aac007c1c450d41da9de009254b5d4206e94000792463b62d16e5400b05a89e,
    0x69916c00084d45d3a1398df4009a44bfcb17d9b30007c9ebdf407f486e0074a221a7d99a0d80,
    0x692a8d4935235a0008c68ef4d8e499d000a125da6db9c66c0008d6b7638e7b30480087cc199f,
    0xb3e0640094a6d89e7b4f52000526c5a38e2e04000066919783ca9c3e0006f76bcd1620abb000,
    0xdd72c29afa4d600077cf88853c0d9b80092e7e81f39f7270008651ccc5fa054dc007819b4fea,
    0x64a80008bb16bb5862751800739c679ce399ad000839d3387e8d78f800859aa6394478934007,
    0xd71ed9c8d5660008a22226886cd7c8006c7120b6c9b1bc0007b6282fdb76854c0088baec51d5,
    0x181000777f9a5f5f52ee800784514298967874008df64fa9dbebf080096169982da9c6000096,
    0xc58d783ccc4000841a337b73c78bc008b58d305a13a490007969af7ecdda3fc007698d7de61b,
    0x180005b1495c8c996be0008266262bff9ab6c007fa278b7a2d5a47c0623b087729ad20000624,
    0xb342c1df870005fcb4a290b6f5a000976379a0b6fb790008070e3f652d66f54072efd508a3a2,
    0x30007ade74f7f62abd000854b1e1143554580082ce6c1e5a4b9d0006eb7102ad8d8850006be9,
    0x673a1fa320007cc23a208a5b0d00077083fd897259b00075e941de7276be8006633e6d56da2b,
    0x0006c2806d7cc497f00074b6b03a80f4de000813861065eb691e00752bf2aa0df4ac00099f75,
    0x742b9972a80817b23466ca211e00abc36ab61255f60006f309ffb4a6dbc0008a2cd7b3f4e3b0,
    0x007bfbf13c1802000008a9f7208a1c7c18009d7f92773e36850007eef6e49f00c05b0080f6cc,
    0x4884761a0093706dc4160dba0007370af0ddddd298006f615effe4571a0006d57495b7ad8aa0,
    0x074e332f45df1138009585db7bceed1a0008cac468f1fb8c700073e283b27038f58007cf9ee7,
    0xf0a2a700086807a21137d11c007fec283d782ae7580858648b61cdfb580078e18e69a6973d80,
    0x7f35e282a87c5e98087c441f428c12bc007b823e56d0702bc007f08a81bceeacc7009b6d5b6e,
    0xb534156066faeb061926ee0008718ab46162caa4007c6a9d39d670c7c008a0b53608dad05800,
    0x5f3b0f8cde62c0007964118ad5af784006cf389a0da817500077bd944e0f1b8b8008038f30ad,
    0x1e588008d41637f997f148008f0dec5d62cf678008b18826e055d41800817a3c9da7bccfc008,
    0x0e6ecc97b9d4b0092b6e89e7670e00007b83bf7a84d1f98006abfdf4955d80d0007acee46ff2,
    0xe5a800799720efac81d2000848254f72cf2db00087dcb4a78075d600065cda7450bf53e0007e,
    0xbaeeb0e58eac6078b302074f4b35c0094c981d7dc6cda00095e9030b83e2c200082c616b6126,
    0xa220097f7d90105c5bf000a021c0aed9e760000745b1f54168b5780068bc0d05509f0d000802,
    0x85fa754aa5b00a2044ef2d83e740006a224f9c0b5c4100074eba9aaf4e00d0007c54f8ff4749,
    0x70008ba001c7b2a08a8006e6818af8b9ddd000913fc389762ec5000767f2f6dac5fb300081d5,
    0xa9f4a28bc8008a43cb6f09ace80009738725df497e30007ca3249d5d4712200970dce049e475,
    0x4008c63e11cff2dbb8007774cd2219cae400074ba4ef9e0838f0008930628d449e2d00075239,
    0xdcc4d3f20006cb5f154c3498e000904ddb6efe68c50007f9bf20968dd181c07863dad1dbdf0d,
    0x00780ea82743b57880086c03ae72149e340089d4ea4a8f0f32800a0abd97d08eff8000998652,
    0x9d0b27b8007503947260575b800849c57f4915bfb80077f7e18b660b6700087dafa5abaf8ef4,
    0x075a82c9b58b6628009fb721880b263e00093e352ed88a77e0008b0fd8c616996a0008c9d4d4,
    0xb719450008f20a1ac9d9d4c00078ae2582a27d428006929e99c014a670008f7582e73dca2880,
    0x7f58e55271ea25f007a130b23dc67694007197c74dfdcf2b8008587be8ae2463e800742f8edb,
    0x7aab5200879930d7aa1356c0083e3b8fa1193c420060ddca32f62acb0007633f13cb8f191000,
    0xa0571b55cc1a6000913d960843cd100008e673bf33ac551800782e16dd62ca0dc00822ffb9f0,
    0x328c8006e40d3b16cd29c0007a7cb8c2583cfe0006386f8f104e21100091fe422f5179140005,
    0xc08b02209db600087e701addd1817000888d5de91f44958009511f4b6f4dc8b00075eb1632d5,
    0x4ec800706d942faca4d48008dc67083075d000008a0f93e55829260008324260ae73755c005e,
    0xd8a361d119e8007338cb401f099a0009bd1e433238dfa00069ad3cb4b0816f0008f002d24b10,
    0xd50008088b627f87444b00924ee22a1fac2b0006efb7e060d80f500085162a34e62520c008c3,
    0xf9ff1fe6560008058cdac837ad29c08de55d5ffc4354800837db6fc0da36d0008d1832c03c37,
    0x38008d80731c10c219000860ea10333c136000886fc5e02710238007a5d7798eee02080091d7,
    0x1f16588b20008e89c7d116f06300076b7a4d26d8f2e80082d7d30a41ea73000761035beabbe8,
    0x600934c94acd8f3ae0007fc7da24720e547c08ca55cb1671f410009bb73fec5cdc6300075d02,
    0x59dfa0b500086c4c3f282748b40088737bf89fec5580089031f4f6b27580007ea583f86917d1,
    0x008aecb88de70e7d80078b30296cce401c00651829daa974c80006bfc8ab5f5069e0007673a7,
    0x56714a40008b3b6e3b97a446000a280eea86ee9e80006cbf1317724e510007c2ce7081d655ea,
    0x068f2f880a5cbfe0007968c54ff86ad8c007637ca8792ac7b000936898cec178e90009f25590,
    0xc859b50006cc3cc667955ed0007145b1a603734900099fa821fc2a0a90008bd8fab10a98a480,
    0x648fca343b0f52000892777d868ca2b8009ebeebc2c3b43900067bc807a9067640008e1f7309,
    0xc987b8007f548cac8af038680588160cf9f7d0400081aed32e1ca5cf90085adbf6f9e7376c00,
    0xba1d403a8eb630008d9d3f264d5e0600083d1c8ba5d594fc0071cda1ab244b8e000769a9473f,
    0x29de8009280540f3b003c00070c32b7b2817350008be5d4fdde101d00055e13fa998cde60009,
    0xb8f8ff7f7bef80087ea0281ee4d7300068469c14c003af000971e8b6961ec4c0008588a64037,
    0x24200078e170b47d4e16000924de49bb564600007c4d223ff8c2e880090c40b33cefe2000084,
    0x2935fdc44b100078b1e1ef655eee4006ee457d5870eb30008067f4205bf779b007921d7f5c0f,
    0x23a0085187e5d4d394bc007933c092f5b11240091d3fcb706b32f0009a6de3ef7964ff00091c,
    0xaee6246c8b00072cfb1ccb56ac9000841e58128923c64006bdb06022f9a83000819369e1d365,
    0xa080706ba9371b03f58007db1f9eff8813040074ac4f9d1152bf0008b9bab9d5a2cec0006f07,
    0x0e90d9b3e0006ecbbe4a90962e000880a6ce498deaf80087e0a6c781465900080e525fb76ec7,
    0x800813ccbe98577baa00828509463f3bf180085f770ec5c586780097e3364e960f5b0009420f,
    0x2a37680c8008a760f490721f980082306df68b3e77a008605d4bdb01b60400841ee87f707520,
    0x007e678b9d2c0d4860063e8ad64131bd90007f65a6b53d2e7ab009211cbfcdd1fe90007aa0b1,
    0x94a1a06f8060f579ff2d8a1d00080f69f36ef5bbee0093ae35b0fdd6b30008726feeb0caa018,
    0x073365dd4122c928006f970816d0d8a80007061defc562a078006cd00728e2d2850008013ca7,
    0xce329280082f25fdf181550e0078f999500fb20e80072b20ec0ac7c380009231ae45f456e700,
    0x6d40a3e4aae55500077b3c0001e66f080072c7c97c580ece00075552d209a60a00008e3efeb4,
    0x9830b00068f117290266c4000778e92c096703c8005d16b3114cdc920008eb874e3b33c22000,
    0x6fede5fe287240009a277ca50964a4000989c9cff2ec1970008075624c337be9bc060b102702,
    0x8f0960057034a0b9bd0060008b5c9c8be110e100077636cd0cbe18980074d944f89275498007,
    0xb4ef5baf829a0005fe6eeae9a12b80007418d3f25ca0328007982e816e5d7c400083f694c5cd,
    0xa3d00080bdb9c0487354e0071ae7b941f1bf50007fc14c881448168a0659d76a550e0320007c,
    0x25b30a2181580076364dc7cbb90a0006fa4418431245e00077f0af9e528501000876045d0446,
    0xc84007c77777dc8b4cfc007ce7d0046e8f9940073829679c9c677800681b8acde4caa0000755,
    0x186234c0a8000979cd2cd86141d0008556b76f89e6f740094ab85e0d4fdc40008545cd8f37b9,
    0xb000866f4adf205f374006ceec2a9f91351000897e9a56c3b50c800959c567a8b115400075a2,
    0xe58b8cd6d0007c38e48646cdf0a00775d33256f6e6a0008f27db5a5dd2790006c3268fdf42c0,
    0x8008acdac3a3aa1758008fde6e6b5e46d28009b3bf6d90391ef000725d7a4b7ef29000075209,
    0x4f7b8b5a00099a805370828dd0008a96602303308c8008bc805f74705ea000812d8c1f44138a,
    0x007a2e20799007b9c007882c657c057f94007c8fcc6ac8ebc80007cb3732807861f0008d99c3,
    0x780d1f1c0074c53d4f024c018005abf09f843c3180009390eb67a6eb990007c368ae6c60d634,
    0x06e211a74220d0e0007c41bf6d23305b80099e7d60a9563fa0006b77c6093901670007ddfb87,
    0xb6940e9007e02ff078341d5e009de7765696232f0008df554a142f0f700071bbf32e8fc45480,
    0x5abb60b88eea180008b9000b58e241f8009510fa0b50881200070d2cce69f435100081b663d6,
    0x7eaad0007ba0dfafd4feac0006e44f95051e07c0007c410e03a61b50c0075ad20a35661eb000,
    0x6a9cde8c0d9f8c008cf9cabaa4d258800855ecd0b3dfd9100085bfb9b26772a58006ddd6ccee,
    0xa568000806919a2904eaf4805d8c11c4a7e6da00079d26116d4cfda40099817db8ef8c8a0008,
    0x999292caf8ec80091842c297a44cc0009669fd7412b459000799489fae19ffac008a65c7e70b,
    0x9a40008781c0a17bfffd0007d2b3d6e621e8a000772d0dac0b59b5000811dcc22eba70310088,
    0xc029a7c21d9800a79c14e46803440007459131c7dd73c0007b431d4711dbe0c008f7c187c2b2,
    0xe43808d3a875cf15b7780080fcaafb313781d808ee54a4b841551000718ad75d883ce6800884,
    0x16be74338a800869ac0097c1ec08008173252bb98c253007f10acb2fc3f8eb8080d4290a947f,
    0x5c007aec06aa818f8180082bda691ea3a7c4007d0b24f58a68df4007583a1d86f9d7b000763d,
    0x46446bc21c009e45acf0e4e762000691acb4f2fff6c0006e5f4ab020783b00085125e67917f3,
    0x2007ca286a453781b4008c8fd872e50ce78007b854aeb479b1980079a5c5d3f7a75b400865db,
    0xcc4ad04820089a291c930dc8400077fa4d5eba79fb0007bc3af0fd7e789c008312e723eec2df,
    0x007ca7fb34397f5b2007dd5c7e6c1c7b4a006dfe5fdd5440c300094d43984e27f9200082438f,
    0xcab57aa80069d4c125062a770007e481733620c11d008ff73ab634eb6f0006ecd5d274237790,
    0x0753875c2e02e22000a8afbc957134ce00090a01b46b540560006e508dc19da9210007686978,
    0xb97bb60005abc66dd33df28000ab57e5d4ed34d80008479cb717d2e57c007e8639d3b7e32e40,
    0x7452b4005ca4508007b95eab3b79e264007c5ba61757440f0006cc15a48de815c00077c6d659,
    0x6b3c00008d2d88a0a0cb7f000a64c202d6523b4000841cf08afa11970009d58b3aa688f89000,
    0x5ab2baa3e5f9f00076c0ed0a3a58ee800845e5d888db3e5c0073f1b3321757710006bc1480bf,
    0xbd5c00075b41084d5b8728007dd8a03679d44900086bceab65d55a84008bb89e17327a540008,
    0x0ed162423d258007eb9ac91aa72033008ec1e30d429c1600085744e358da5dd8006b5ef144df,
    0x20000088bf54ea720b400007814ddfeaa29d48009ef0ee95a51305000a3a3f62a59ad520008c,
    0x179e44539100008470ffac1927498007f10f1cfe14186d808b0979ccfffcba80074c3bec73fa,
    0x6cb007ec628aae7171cf0078dd8308bf305dc007ea732e8500e35f0093493a5a7d2641000724,
    0x8330e8cb5d800a75807290253c00006e3a797bb394ed0009a3098516397400008138ed81f8cb,
    0x08008bb557fcb044d8000758754651b605f8005f11da5c3db6e600082cb56d9feaa2a20073a0,
    0x5d806fb850006280ba28679bde0007888211ae6b6acc008691b42d912c0f00088156f4277477,
    0xc007b275edc688e364008d662d3bac071700089f239d28901960007c529f12940d03a006ed5b,
    0xaf2421f08009898253c4ce62600080ecef70043b861807d64f5cd4a81f16007d933da6512a7e,
    0x00714376d9424b8f0006a967680b6afe600081427680360c497007cb94aba01526d800745f62,
    0xc3d5ec1000788c77e938b516800923dcf3b3c73af00086d546700c549dc00833445df7e8cab2,
    0x091aa309f00276000060dd21015c725300066759fdfa612e0000722ac71990f22f80097432d4,
    0x1ac01e00099614d66812daa0007b4d5c0c91dba4c007de9cabc7fd12f40084e0c194384f7ac0,
    0x88f84415c98924800915d5f36b13c8800069be64f9e9fabf0006ff0aea0eb3e690009337e459,
    0x4249c000755b949d6899ec0005fa0d86c07362c0008657de0429012600077d3059a4296cc000,
    0x29be41c12e170000903e87b2bfe7a10006b7d6361b60ce60007ed55a9aee4f2940086f5ca1d9,
    0xe2488008dd1d1240382a98008bb0983dffbc628007d085424d6ad2760070f598c94e3dad8009,
    0xcc7b711161570006157fd47f05b9300085856dca30e09740071f180358d269d800778a3ff0f4,
    0xa10c007c25be8fa69384a006cc99dec6af08100064bd3de2156d090008c609ce53517cc80070,
    0x897311033e12006fa81d665eeebb00088266a6db3223880063ef6ec3a034c20007b2ab4f9fa0,
    0xf10007dd9b5ed3c967e0007ef1bd6cc634062008b25642db52d7f8007cdad719fa0dc10007d7,
    0x6a59ab3d5e80087992f4d72735040089fa05906a4a2580085f3480c3d3b7ec009a79d1890584,
    0x00009f839f14985c9c00083e2b0557ea79d0007f1619f282ed69700955f2d33a00b530008925,
    0xcfcce1de40009b5cbb6b917aaa0007649685f9c134500094a7f9081f35210007784f68ba3182,
    0x00088d3cebe3bd64c0008f705c0bf739cb80090951aeb17fd1600072f854c68bd29a00075570,
    0x68f662cd000691456b58a57d9000a0e04cc31935b600083dea862992090000845efaef728aa8,
    0x007deb7c40db844ea007c11d9c31370336008f5f412629e6c38006b88a45dd1dede0007037f7,
    0xbb8c20f0006972ee307c0c21000a05d5ddf2dd82200076bc4ec21c582e8009c886cbaf801a00,
    0x0956c70b38486d00007d3c0f5f0691d820079c84ad236d688000852b80edab5a960009515dac,
    0xc47866000a93f40cbf84162000b6bccdb336a5c00006524feed13233a0009b42c9f24976b100,
    0x8bd676d9ae71d80007e83114a60343350075621c2b09f69c8007ef84aadeae9aa50077984821,
    0x9006a8005cbcf467fcf5a400097a932292a83bc000831275ee1e408440081a75fc45bfcfd400,
    0x8df3579ce704bc0090f9b7bc83162900073786b7753e691000824f887c10d3a400078230c9d7,
    0x0000000000000000000000008d5944bdf0bb938009306313cdd5733000792213fdb742720007]

def c520N6 : ℕ := ofChunks76 [
    0xa4291000a381f4d966dd4e0008755a697c045f180093ca7f5fa096a100087d0789d8af12d000,
    0x0bc7e438e8efa00085b37b70d1f9708007316f7673e18a20008627dc69019c3580092c942b30,
    0x8698000970f89613266960009a727914a6f8b9000829fbf311e21e380087071dfbc6a927000b,
    0x37ccfe804da30007d3cfe7de73206000a1f6644478d4300008b8d69ff427aed8009c06ac308c,
    0xfd6000a1c15d42f4f43e000921481f847b2620009ae4632b6c0d8a00084bbd7577dabe28008f,
    0xa9ce10952e300082d2e01ffe4a83000877235cf2e933700092402003cc2a2d000b0ef9e9b70f,
    0xfc0007e822e6eeb56bc0007def2ad5056fb3000afde6baf3d5c5c0008ee2e92f51733c0008bc,
    0x59a2399d4a00081791976fda48c0009d2346307ee0e20009966a4c082e3920008978f7eb73a4,
    0x30008b00e5f563b47c000a16d54d6ad55ac000a31f90d1e723e00008f504f73eb3acf8009208,
    0xde8c25da20008b09214fe5d676000758b4b4e1a3a3d0009a718d445eb55e00090bc63f5ba1ff,
    0x800846bea9b6ddbc98008fa9f400a57b910008f0a93364c63e400087a0d7520b8c5c8009932a,
    0xc5565d740009d5874114c55ab0008a7292f423dd1e8008ccd5803a241ec00080169c4fb15165,
    0x00999f2b7fe98dd20008edf94311a32760009660f5e219514e0008fcf2a7aa879d6000a87a65,
    0x496de78000803385e17f7c13800a1ce8017c559b200093c9185d6decc3000937ef5190cfd480,
    0x08d60c1e3f9b15200087093fd738b2d880092e4a59ee67e500009595f6e58f0a66000a72a737,
    0x659f9c800996038bda9571800086602a855ea0980006efb06c242f52800095fa9ceae0a79200,
    0x75d7b9ff4a365600093c5a0c70482e900085d431f8a8777480083310a074991288008b49cd77,
    0x35f690008367b3584d31d700087f25241cded360006d4c0b27a231f40009df08e1279e972000,
    0xb2fb5266b30fa0008e1567d350488c0009400511c734cf60007ec0e43ba18b520007c05bbc96,
    0x866a80081f391695120fa8008de62a8363d9148007acf9515eb5cc50006ade26cda226b20008,
    0x80389b7ccca70007f4c6994ff4bba00092cf35ba670cd2000a0caeb33e2549a00083406f38dd,
    0x2660008f2bcf7ba024e200080933c85755caa800851d82cccf805c8008b71090c038f4300088,
    0x3e3bee2e0e4000963483004024160008ccc32a55fe63c00093a12f9105a8610009bdbeefe84e,
    0x9b8008c016e47ae40bf0007dabbbbe4e1d3d000974afa646d1bde0008dc0746f8b117c0009cb,
    0x9627b634340008e05a92770242080095b794aaa2da0e0009321711a3c748900089467550c301,
    0xd0008367e481a426bb800a5b59a8f287cf40009474d2071cc41f000811fe3fc4488b900083bd,
    0xbdaa294e680091f1f5cfc40bb80008720464e9498f6800954a77b2fb47e400089d7dd58ef077,
    0x000a062e01b7a25e60009d8a7a163389d30009d4d955762f32e0008394793650c5158008f3d0,
    0xca41ef80000942df2ed13d10b00088917b6e4e9ab700087e5f12e8edaef000a449b6c9d17cb8,
    0x00a7b719546bda02000861b8f5091ad4a0008939facb4fd1178008a90e429e25d9a0009a14cf,
    0x7d4942c0007f0bd6cd09d758000909c851b13b22b0008d5199b215ec8f00067de8cf590de2a0,
    0x06cc199120398c00007562e002e867de0009f478d6db2c8ee0007086cf109abd41000962ab45,
    0x6bcd73000870d29146345b880085932e4f030c2b8009297e19c72d1c900083c8714c3d089600,
    0xaa24d256990d4600097f7889603b2c80008ea05c3a0951298007f0ebde8480b9f0009554cebd,
    0x1ff280008ff2b96ea486b7000974024bc3dd48e000ab6fab51bbf51400097abd5bc5eeac9000,
    0xad6fc01eaf0410009b6575b49cc8bc0008ae58772a65e1d0007c67cc6f2b0f050006be53bb4e,
    0x44860009378db5148629b0008b9d2a2851e1d40007efaa0f2abd33000098b2c8e58ce1590008,
    0xb3ffc8f7bbd200095fd516ab503a300085449db7c3ea4d80091ce289b7e7a610008a3ff4c2dd,
    0x953000782a6af6caea330008538a74b310a9c00070d1c4a1cd8c020008f3f492f2d39f10008e,
    0x5cfd3262dbc0009489a063ef3ae20008d6d5eb0d484ed80088190e089d2988800764efa94fec,
    0x000008b13162ebc259e80087943ba9f0fb32000a03170ab6cf298000971c437c59d7d2000a17,
    0x18da0b0dcc0008274823334042e800751db7b78dcf370009004f1d2f88cf4000a15cda0c53f6,
    0xd800880eeb5354ff8d000903a4193ccf3050009555c983d52216000863a6b91bd39af800a01f,
    0x9c50c217f000a3c1357a036f40000935a6f54806b3400096447cc3e38c290008a632dcecbb96,
    0x000b482efc85666e40009b490ce26ec2aa00092a78a02b96b36000813a8a344175938008ccd8,
    0x03769bed800a48025e1a27baa0009e1eb133599ead0008ac7417a1b8a140008894c11ef023bd,
    0x0099f94947f44446000809aeaf546342b00086f51d650197878009aa214a11dfc4e00085760a,
    0xa147e62000acb83ace1bf12c00086b641f725d967800963ab303cf19cc00091b437020b35140,
    0x08dcbc685dbbb470007c91595d8a7b250009b3e75a9bfd2730008a02729a2480cf8009aa5050,
    0x5571500008c7301eaaebe240009e915e5a30b33c0008cfa83a1475cee800a674276976a87400,
    0x86816cc5fbe98a000a7251474966d5c0009390534c40ea6400081eaa63ff43c3c8006b3b48ab,
    0x6451e00083901c6c5e36498007fbc98615a7fff000809019a6122a34000946a68f7c453d6000,
    0x9358e2276cb81800a1d8ef432ec1d0000997fa7093744780009f0b8297c4421c000a5e7ada05,
    0x961a00082b3f81965bdc58007f36e1f2db9c08000935aa43acbe0e00009c273a256c076e0008,
    0x8bb4e48210c6000800e5f3292fbbd0008adbb74658e8c0000934a42bf533bea00090cd044d9a,
    0x760000774406f22efb8300098c95cec1818e00008b46912c1d278200096cdc3b7d4c6400006e,
    0x99b569f0e210007045d22edad62400079651da251d9ee0008e37fc117481230007bf25604aee,
    0x6e0007ddb3f3a1bd4650008ba5a7470de0d80007f8544d81858820009e69f5f80317fe0009c8,
    0xb7268877c8000839cce108f9f7f0008216b0a35f6cff0008f86700bbcbac9800951015abb073,
    0x900082838617b768a080096060a89f480240009922da2c71694800096fb63940eedfc0007f01,
    0x8a24514ce0008f5ea620928eaf8008d2ac43cb7a2af000987e777ab7da32000736bd96559efb,
    0x0008496c85dda706d800896391509133b9800a24de3039f09e800081b428bd8c9775000a0839,
    0x8a4d5824800819785bd9ef61f0008caf3ea5c6839a000963d8553c2cc080008cb5438810682f,
    0x008d4ab6433db39f800974ba33b86933700087dd22dacb95118006754400dbed5700008b7e77,
    0x087157f0009a7f28fbf972f40009a5491aec5db6000097c69584459cd10009b0363ee00ac960,
    0x0973eadc61bec85000884feeb8fa211e000aa14c17c46a09c000901109bd9a41b1000742ec5f,
    0x1249b100085a090eab71a790008ba5c3bf739c0800098e985a6aa79d5000924b57399029fa00,
    0x9694e5c3e592700009080c7ed40d9680008d08e1024ea3da800935335a4d7eb430009126c06a,
    0x7e0ce000a7a2a20e064a06000b3365066d3154c000a23c1ea39a5c780008394c0336c81bd800,
    0xc73200246ea740009e15995f60a00a00090ddfead54bb4700072d57428b35fca00088bab4590,
    0x6ed60009a046ddf7077e20009223d16a0d09e20007e308b2fec414a0007c52f50bc4067e0009,
    0xccb80fb1c97c0007aff91d730e8bc00080832c216aca6c0009d0ff1092427150009df787313e,
    0x8180007f0dd3865fcbf3000955fa4f8cd7b710009d3ff1e9c6b4630009814344ef6c0a20006f,
    0xb832ccb915800098d28ae8ac517f000ac76f0edcf41400008a13e1221a26aa000a9155e28a6c,
    0x3d0008710536cac6b430007fcc9a658b585b0008ef72f39ab28ee800954a258f5eb54e000a30,
    0x8e8268e4950008d2706d7f3d5d68008f77ef219706650009483601910ba1f0008f3f0bda7ceb,
    0xf0009d8a53f4056c3c00078b56efc9c1022000876cccd8852044800879c2c0fa1b61180092ef,
    0x9948e8d130008c091e45dbc50d00084253285c8b4ae0008ea979a2c46c6e00090c03c279b0bc,
    0x0009203f2e9d48e2c0009fb509c5b0c621000714bf93a1eb3e10007f3613fabdd9500007e9c8,
    0x66e94cab0009be0cc3c678a2000099fe53d956c87e000a29edeb3a9c8f60007e21af1860c5dc,
    0x009ccc734231dc5c00096799ce8d898560008e6e601489aae3000893bea9015219b0008c633e,
    0xb0ed4b8000919d1df9e50a360007fc03f2dbb80a800087dad8fdede714800a96e7734110e0a0,
    0x08857e30bc150ab0008c1c450d41da9e000a254b5d4206e94000892463b62d16e5000c05a89e,
    0x69916c00094d45d3a1398df000aa44bfcb17d9b40008c9ebdf407f48700084a221a7d99a0d80,
    0x792a8d4935235a0009c68ef4d8e499d000b125da6db9c66c0009d6b7638e7b30400097cc199f,
    0xb3e06000a4a6d89e7b4f52000626c5a38e2e04000076919783ca9c3e0007f76bcd1620abb000,
    0xdd72c29afa4d600087cf88853c0d9b800a2e7e81f39f7280009651ccc5fa054e0008819b4fea,
    0x64a80009bb16bb5862752000839c679ce399ad000939d3387e8d790000959aa6394478930008,
    0xd71ed9c8d5660009a22226886cd7c0007c7120b6c9b1bc0008b6282fdb7685500098baec51d5,
    0x181000877f9a5f5f52ee800884514298967870009df64fa9dbebf0000a6169982da9c60000a6,
    0xc58d783ccc4000941a337b73c78c0009b58d305a13a490008969af7ecdda400008698d7de61b,
    0x180006b1495c8c996be0009266262bff9ab70008fa278b7a2d5a4800723b087729ad20000724,
    0xb342c1df870006fcb4a290b6f5a000a76379a0b6fb780009070e3f652d66f00082efd508a3a2,
    0x30008ade74f7f62abd000954b1e1143554600092ce6c1e5a4b9d0007eb7102ad8d8850007be9,
    0x673a1fa320008cc23a208a5b0d00087083fd897259b00085e941de7276be8007633e6d56da2b,
    0x0007c2806d7cc497f00084b6b03a80f4de000913861065eb692000852bf2aa0df4ac000a9f75,
    0x742b9973000917b23466ca212000bbc36ab61255f60007f309ffb4a6dbc0009a2cd7b3f4e3b0,
    0x008bfbf13c1802000009a9f7208a1c7c2000ad7f92773e36840008eef6e49f00c0580090f6cc,
    0x4884761800a3706dc4160dba0008370af0ddddd298007f615effe4571a0007d57495b7ad8aa0,
    0x084e332f45df113800a585db7bceed1a0009cac468f1fb8c700083e283b27038f58008cf9ee7,
    0xf0a2a800096807a21137d120008fec283d782ae7800958648b61cdfb600088e18e69a6973d80,
    0x8f35e282a87c5e80097c441f428c12c0008b823e56d0702c0008f08a81bceeacc800ab6d5b6e,
    0xb534100076faeb061926ee0009718ab46162caa0008c6a9d39d670c80009a0b53608dad06000,
    0x5f3b0f8cde62c0008964118ad5af780007cf389a0da817500087bd944e0f1b8b8009038f30ad,
    0x1e588009d41637f997f140009f0dec5d62cf680009b18826e055d42000917a3c9da7bcd00009,
    0x0e6ecc97b9d4800a2b6e89e7670e00008b83bf7a84d1f98007abfdf4955d80d0008acee46ff2,
    0xe5a000899720efac81d2000948254f72cf2db00097dcb4a78075d600075cda7450bf53e0008e,
    0xbaeeb0e58eb00088b302074f4b36000a4c981d7dc6cda000a5e9030b83e2c200092c616b6126,
    0xa2000a7f7d90105c5c0000b021c0aed9e760000845b1f54168b5780078bc0d05509f0d000902,
    0x85fa754aa6000b2044ef2d83e740007a224f9c0b5c4100084eba9aaf4e00d0008c54f8ff4749,
    0x80009ba001c7b2a08a0007e6818af8b9ddd000a13fc389762ec4000867f2f6dac5fb300091d5,
    0xa9f4a28bc8009a43cb6f09ace8000a738725df497e40008ca3249d5d4712000a70dce049e475,
    0x0009c63e11cff2dbc0008774cd2219cae400084ba4ef9e0838f0009930628d449e2d00085239,
    0xdcc4d3f20007cb5f154c3498e000a04ddb6efe68c40008f9bf20968dd180008863dad1dbdf0d,
    0x00880ea82743b57880096c03ae72149e300099d4ea4a8f0f32000b0abd97d08eff8000a98652,
    0x9d0b27c0008503947260575b800949c57f4915bfc00087f7e18b660b6700097dafa5abaf8ef0,
    0x085a82c9b58b662800afb721880b263e000a3e352ed88a77e0009b0fd8c616996a0009c9d4d4,
    0xb719450009f20a1ac9d9d4c00088ae2582a27d428007929e99c014a670009f7582e73dca2800,
    0x8f58e55271ea260008a130b23dc67690008197c74dfdcf2b8009587be8ae2463e000842f8edb,
    0x7aab5000979930d7aa135700093e3b8fa1193c400070ddca32f62acb0008633f13cb8f191000,
    0xa0571b55cc1a6000a13d960843cd100009e673bf33ac552000882e16dd62ca0e000922ffb9f0,
    0x328c8007e40d3b16cd29c0008a7cb8c2583cfe0007386f8f104e211000a1fe422f5179140006,
    0xc08b02209db600097e701addd1817000988d5de91f4496000a511f4b6f4dc8c00085eb1632d5,
    0x4ec000806d942faca4d48009dc67083075d000009a0f93e55829260009324260ae737560006e,
    0xd8a361d119e0008338cb401f099a000abd1e433238dfa00079ad3cb4b0816f0009f002d24b10,
    0xd50009088b627f87445000a24ee22a1fac2c0007efb7e060d80f500095162a34e625210009c3,
    0xf9ff1fe6560009058cdac837ad30009de55d5ffc4354000937db6fc0da36d0009d1832c03c37,
    0x38009d80731c10c219000960ea10333c136000986fc5e02710240008a5d7798eee020800a1d7,
    0x1f16588b20009e89c7d116f06300086b7a4d26d8f2e80092d7d30a41ea73000861035beabbe8,
    0x800a34c94acd8f3ae0008fc7da24720e548009ca55cb1671f41000abb73fec5cdc6400085d02,
    0x59dfa0b500096c4c3f282748b00098737bf89fec5600099031f4f6b27580008ea583f86917d1,
    0x009aecb88de70e7e00088b30296cce402000751829daa974c80007bfc8ab5f5069e0008673a7,
    0x56714a40009b3b6e3b97a446000b280eea86ee9e80007cbf1317724e510008c2ce7081d655e8,
    0x078f2f880a5cbfe0008968c54ff86ad90008637ca8792ac7b000a36898cec178e8000af25590,
    0xc859b50007cc3cc667955ed0008145b1a6037349000a9fa821fc2a0a80009bd8fab10a98a400,
    0x748fca343b0f52000992777d868ca2c000aebeebc2c3b43800077bc807a9067640009e1f7309,
    0xc987b8008f548cac8af038800688160cf9f7d0400091aed32e1ca5d000095adbf6f9e7377000,
    0xba1d403a8eb640009d9d3f264d5e0600093d1c8ba5d595000081cda1ab244b8e000869a9473f,
    0x29de000a280540f3b003c00080c32b7b2817350009be5d4fdde101d00065e13fa998cde6000a,
    0xb8f8ff7f7bf000097ea0281ee4d7300078469c14c003af000a71e8b6961ec4c0009588a64037,
    0x24200088e170b47d4e16000a24de49bb564600008c4d223ff8c2e8800a0c40b33cefe2000094,
    0x2935fdc44b000088b1e1ef655eee0007ee457d5870eb30009067f4205bf77a0008921d7f5c0f,
    0x2400095187e5d4d394c0008933c092f5b112000a1d3fcb706b330000aa6de3ef796500000a1c,
    0xaee6246c8b00082cfb1ccb56ac9000941e58128923c60007bdb06022f9a83000919369e1d365,
    0xa000806ba9371b03f58008db1f9eff8813000084ac4f9d1152bf0009b9bab9d5a2cec0007f07,
    0x0e90d9b3e0007ecbbe4a90962e000980a6ce498deb000097e0a6c781465900090e525fb76ec7,
    0x000913ccbe98577bb000928509463f3bf200095f770ec5c5868000a7e3364e960f5c000a420f,
    0x2a37680c8009a760f490721fa00092306df68b3e780009605d4bdb01b60000941ee87f707520,
    0x008e678b9d2c0d4880073e8ad64131bd90008f65a6b53d2e7a800a211cbfcdd1fe80008aa0b1,
    0x94a1a0700070f579ff2d8a1d00090f69f36ef5bbf000a3ae35b0fdd6b40009726feeb0caa020,
    0x083365dd4122c928007f970816d0d8a80008061defc562a078007cd00728e2d2850009013ca7,
    0xce329200092f25fdf18155100088f999500fb20e80082b20ec0ac7c38000a231ae45f456e800,
    0x7d40a3e4aae55500087b3c0001e66f080082c7c97c580ece00085552d209a60a00009e3efeb4,
    0x9830b00078f117290266c4000878e92c096703c8006d16b3114cdc920009eb874e3b33c22000,
    0x6fede5fe28724000aa277ca50964a4000a89c9cff2ec1980009075624c337bea00070b102702,
    0x8f0900067034a0b9bd0060009b5c9c8be110e100087636cd0cbe18980084d944f89275498008,
    0xb4ef5baf829a0006fe6eeae9a12b80008418d3f25ca0328008982e816e5d7c400093f694c5cd,
    0xa3d00090bdb9c048735500081ae7b941f1bf50008fc14c88144816800759d76a550e0320008c,
    0x25b30a2181580086364dc7cbb90a0007fa4418431245e00087f0af9e528501000976045d0446,
    0xc80008c77777dc8b4d00008ce7d0046e8f9900083829679c9c677800781b8acde4caa0000855,
    0x186234c0a8000a79cd2cd86141c0009556b76f89e6f7000a4ab85e0d4fdc40009545cd8f37b9,
    0xb000966f4adf205f370007ceec2a9f91351000997e9a56c3b50c000a59c567a8b115400085a2,
    0xe58b8cd6d0008c38e48646cdf0800875d33256f6e6a0009f27db5a5dd2790007c3268fdf42c0,
    0x0009acdac3a3aa1760009fde6e6b5e46d2000ab3bf6d90391f0000825d7a4b7ef29000085209,
    0x4f7b8b5a000a9a805370828dc0009a96602303308c0009bc805f74705ea000912d8c1f44138a,
    0x008a2e20799007ba0008882c657c057f90008c8fcc6ac8ebc80008cb3732807861f0009d99c3,
    0x780d1f200084c53d4f024c018006abf09f843c318000a390eb67a6eb980008c368ae6c60d630,
    0x07e211a74220d0e0008c41bf6d23305b800a9e7d60a9563fa0007b77c6093901670008ddfb87,
    0xb6940f0008e02ff078341d6000ade776569623300009df554a142f0f700081bbf32e8fc45480,
    0x6abb60b88eea180009b9000b58e2420000a510fa0b50881200080d2cce69f435100091b663d6,
    0x7eaad0008ba0dfafd4feac0007e44f95051e07c0008c410e03a61b5100085ad20a35661eb000,
    0x6a9cde8c0d9f90009cf9cabaa4d258000955ecd0b3dfd9100095bfb9b26772a60007ddd6ccee,
    0xa568000906919a2904eaf0006d8c11c4a7e6da00089d26116d4cfda000a9817db8ef8c8a0009,
    0x999292caf8ec000a1842c297a44cc000a669fd7412b458000899489fae19ffb0009a65c7e70b,
    0x9a40009781c0a17bfffd0008d2b3d6e621e8a000872d0dac0b59b5000911dcc22eba70300098,
    0xc029a7c21da000b79c14e46803440008459131c7dd73c0008b431d4711dbe10009f7c187c2b2,
    0xe40009d3a875cf15b7800090fcaafb3137820009ee54a4b841551000818ad75d883ce6800984,
    0x16be74338a800969ac0097c1ec00009173252bb98c250008f10acb2fc3f8e80090d4290a947f,
    0x60008aec06aa818f8180092bda691ea3a7c0008d0b24f58a68df0008583a1d86f9d7b000863d,
    0x46446bc22000ae45acf0e4e762000791acb4f2fff6c0007e5f4ab020783b00095125e67917f3,
    0x0008ca286a453781b0009c8fd872e50ce80008b854aeb479b1980089a5c5d3f7a75b000965db,
    0xcc4ad04800099a291c930dc8400087fa4d5eba79fb0008bc3af0fd7e78a0009312e723eec2df,
    0x008ca7fb34397f5b0008dd5c7e6c1c7b48007dfe5fdd5440c3000a4d43984e27f9200092438f,
    0xcab57aa80079d4c125062a770008e481733620c120009ff73ab634eb6f0007ecd5d274237790,
    0x0853875c2e02e22000b8afbc957134ce000a0a01b46b540560007e508dc19da9210008686978,
    0xb97bb60006abc66dd33df28000bb57e5d4ed34d80009479cb717d2e580008e8639d3b7e32e00,
    0x8452b4005ca4508008b95eab3b79e260008c5ba61757440f0007cc15a48de815c00087c6d659,
    0x6b3c00009d2d88a0a0cb7f000b64c202d6523b4000941cf08afa1197000ad58b3aa688f88000,
    0x5ab2baa3e5f9f00086c0ed0a3a58ee800945e5d888db3e600083f1b3321757710007bc1480bf,
    0xbd5c00085b41084d5b8728008dd8a03679d44900096bceab65d55a80009bb89e17327a540009,
    0x0ed162423d260008eb9ac91aa72030009ec1e30d429c1600095744e358da5de0007b5ef144df,
    0x20000098bf54ea720b400008814ddfeaa29d4800aef0ee95a51304000b3a3f62a59ad520009c,
    0x179e44539100009470ffac19274a0008f10f1cfe141870009b0979ccfffcba00084c3bec73fa,
    0x6d0008ec628aae7171d00088dd8308bf305e0008ea732e8500e36000a3493a5a7d2640000824,
    0x8330e8cb5d800b75807290253c00007e3a797bb394ed000aa3098516397400009138ed81f8cb,
    0x00009bb557fcb044d8000858754651b605f8006f11da5c3db6e600092cb56d9feaa2a00083a0,
    0x5d806fb850007280ba28679bde0008888211ae6b6ad0009691b42d912c0f00098156f4277477,
    0x0008b275edc688e360009d662d3bac071700099f239d28901960008c529f12940d038007ed5b,
    0xaf2421f0800a898253c4ce62600090ecef70043b860008d64f5cd4a81f18008d933da6512a7f,
    0x00814376d9424b8f0007a967680b6afe600091427680360c490008cb94aba01526d800845f62,
    0xc3d5ec0000888c77e938b516800a23dcf3b3c73b000096d546700c549e000933445df7e8cab0,
    0x0a1aa309f00276000070dd21015c725300076759fdfa612e0000822ac71990f22f800a7432d4,
    0x1ac01e000a9614d66812daa0008b4d5c0c91dba50008de9cabc7fd12f00094e0c194384f7b00,
    0x98f84415c98924000a15d5f36b13c8800079be64f9e9fabf0007ff0aea0eb3e69000a337e459,
    0x4249c000855b949d6899ec0006fa0d86c07362c0009657de0429012600087d3059a4296cc000,
    0x29be41c12e170000a03e87b2bfe7a00007b7d6361b60ce60008ed55a9aee4f2900096f5ca1d9,
    0xe2488009dd1d1240382aa0009bb0983dffbc620008d085424d6ad2780080f598c94e3dad800a,
    0xcc7b711161570007157fd47f05b9300095856dca30e09700081f180358d269d800878a3ff0f4,
    0xa110008c25be8fa693848007cc99dec6af08100074bd3de2156d090009c609ce53517cc00080,
    0x897311033e10007fa81d665eeebb00098266a6db3223800073ef6ec3a034c20008b2ab4f9fa0,
    0xf00008dd9b5ed3c967e0008ef1bd6cc634060009b25642db52d800008cdad719fa0dc10008d7,
    0x6a59ab3d5e00097992f4d72735000099fa05906a4a2600095f3480c3d3b7f000aa79d1890584,
    0x0000af839f14985c9c00093e2b0557ea79d0008f1619f282ed69800a55f2d33a00b540009925,
    0xcfcce1de4000ab5cbb6b917aaa0008649685f9c1345000a4a7f9081f35200008784f68ba3182,
    0x00098d3cebe3bd64c0009f705c0bf739cc000a0951aeb17fd1600082f854c68bd29a00085570,
    0x68f662cd000791456b58a57d9000b0e04cc31935b600093dea862992090000945efaef728aa8,
    0x008deb7c40db844e8008c11d9c31370338009f5f412629e6c40007b88a45dd1dede0008037f7,
    0xbb8c2100007972ee307c0c21000b05d5ddf2dd82200086bc4ec21c582e800ac886cbaf801a00,
    0x0a56c70b38486d00008d3c0f5f0691d800089c84ad236d688000952b80edab5a96000a515dac,
    0xc47866000b93f40cbf84162000c6bccdb336a5c00007524feed13233a000ab42c9f24976b000,
    0x9bd676d9ae71d80008e83114a60343380085621c2b09f69c8008ef84aadeae9aa80087984821,
    0x9006a8006cbcf467fcf5a4000a7a932292a83bc000931275ee1e408400091a75fc45bfcfd000,
    0x8df3579ce704c000a0f9b7bc83162800083786b7753e691000924f887c10d3a400088230c9d7,
    0x0000000000000000000000009d5944bdf0bb94000a306313cdd5734000892213fdb742720008]

def c520N7 : ℕ := ofChunks76 [
    0xf928000074d7ec1e95888d80068d450a1f751870007ac181bc9a6c320007d98572fa3ac2ee00,
    0x9c1d71f051d348007d0c2df14b524a2008e9da9ac29d8298007f66bf83c0ff71600a9936a7e5,
    0xd58800064a7b08649da8f00077f11e1603c2f500095d29fd41f2be90007f2807eefc632ab008,
    0x4285c251f64f00081fd7410adc9ec000976a1e113188e40007c323dce30f7db0008e6075a1ed,
    0xa8300078a6f59abb59aa4008b9c0d4929bfab8007a6a6442363fc1c00859a87959430c080097,
    0x1cdf978a89e0007871fdc93cb7ef0005f85b491a505340008258687b649920c009af4c6187de,
    0x1c00073d13209ee32fb0008a5da31d6fb69f80080f284c540fc202005d2574e81416fc00065e,
    0x2a10264d6c0008b3f7ecf732e7800067310c697a2184000765e1501d77bc88006abd0fd43c1d,
    0x37007d4fa11e879a3940090c1659f474bd40007603d05e91a13f8006679af9c2d39ee0008931,
    0xd81037bb1800934cc59afc71c70008a07834b6c222480086ad0415c8fad38008026fa22eedce,
    0x0007e413575b6bbf830084d636f2913c134007a337fd3094da240069f57c3cfcb7390008c5c6,
    0x7770fb3100082875c7604bab4600903b99fcf0d42400073dec3181887838005d2223960a6864,
    0x00a1b476cada0bdc0008abfe9d8df2dec8007d75ef270f106260083746f09a1aef2e00687064,
    0x625eaf78405dd21ab49cccfe0007dc7a082d5524360077eb42832a59d700067bedf8111539f0,
    0x097cfc0ffeaa82b000514ff82c6dd762000a3c507689e58860007682612ed2a9430007f801f7,
    0x3c9a7b00074a08fdea1c5dd80082595788369e376008bc2ef70cc85370009020de7954b6dd00,
    0x9138458b065236000867135c9503b300008b1dbe91639ec800071d3fae2b24e5f80091b0aa8d,
    0x328fe00084e4161e1ce87600084d54e3c27625340076b4465127543900076748bbaa361fc000,
    0x74d61e4172ede8008d7b5ef2093d50800825143ca02645cc0087c43230e6a4e440095fbe5b2b,
    0xf01920064b0c27f1f799200085690d37d465cf4006fa3ab6a49d765000917cd0195f73ff0007,
    0xc41dff4b638280063e17664dc99fc0007c4bd6713db877a0066e99f7de6944f00083f174c701,
    0x1f700073b3d59f59774100076afa0b490108980088e05b62f3baaa0007394d20199849200075,
    0xd2bc86158f20006823dee620329d000917ee54990022c0009ab4ec080fe12c000938b5408ab7,
    0x122007c93eb240b01b6e008bc0573eb81d268007ed949e29aca309006c57e8fcfc4040000946,
    0xd30598fadc0008694e4731ce99a00086d881a74475cd80082019100c274250007d97fdb72e86,
    0x18007aec179f678e48c008973abcc66dfdd8008496e3fdc65f72c007d308b7a31cd362006460,
    0x22e9de4670007e0e3be12b6fe6c007456ce35ccd34e8006a1e5c40dbc9a100074c1f05a7dabf,
    0x0008513ca8e7635f180066429b7cb3448200085e148f87f21438008983e6f745238400075d63,
    0xb3ea6cb340079158d48f8e2870007f93792e0df328a807bd19aca4c14a900088aad280f230d7,
    0x00768d74eefed1640006c6726b785227900086cda67be76cb0c006fe329e4bc9a390007b4e49,
    0x5c7769c800783c67db9147370007a1a6566c548b48007210b896354bef0007161bd0c4ab6ab0,
    0x08a1c0ab7cd2a9680081c6d02ad1af6170066a44267dca92a0007bf1ef4c59eebf80074ccf4c,
    0x96f44dd0083d5058b672211000a3b6429a0b0e180006963f09a44cff40006b786d8182c8b000,
    0x7743f1d4221e0200080988a15376b99a0087e1e29e6761220007cd2e9d8a48d7260080d9e0d3,
    0xc335700073f90a9a17601e8007ae2a257fa01860007a7cfd9ac53424800791493f9bd66eb000,
    0x59a892d60ce788008975c13bc67a2e80073a184ad1b9bfa00081804ba92329e1400912f85143,
    0x52010007cd18ce221e171600717f78bd178d56000866a031ac7f1c2c007dcc795a638d77a008,
    0x6454500a5ae380082b03ac624775c4007ab29754387f650006d0f1cca37cb290007003765b54,
    0x00200092bb0d9bb722e900078da2a0a10912780062b67e20db41bd0007d32f63f81e88b40071,
    0x71a678320d3000888e02e604a76b0007c7b728d20fa9f0006ce9c01d38c58e0006d93ccec75f,
    0x03800a072608f1916160008abe460d954f428008cd0756d91e5c8800aa75b8d696aa7200093c,
    0x0d2c70ea31000694751c21f8273000840d65001a37480007c961eef0b5eb64007b2f69ff44f2,
    0x9e00701cda6a2b363e80078ac851f69c00d80076ccd143a942418009ce89bf984132400062a2,
    0x7c6f43d8d800803c3afd1e23838a07c579bbe1d5ec3a0089754fc3b339e9000824e1f112ee24,
    0xb007285d4cf32a2830007f07bdb92611cf6807704bcb512d3ae0007fd3e4da5628c20208f38e,
    0x66562298800583732b5df6276000964d84c48868ab00071a21892ad901e80081a155663157a3,
    0x0087b59b178e806d4007c3fa3461170d14009fcb5c8c749ecd00078310a5ddd3400c00866bb3,
    0xe075bd00009ff5681a79d08c000884a380643d6e580083f677582391570008630b7bb0c09494,
    0x08a67c4343a1afc0008954c729bdf26d8009ca4ee30daf43100083248da5f92e73c008ed7137,
    0x23cde2e007b9bbaa17b913fc008b9b9c2095cec98008293e48a9ecb6e20093eba4bbad9b9c00,
    0x770ec1640cfe1b800a08df088a013ea00075c3a74a5a0c9900080a2604ebe621b88081ff541a,
    0xad9c880092b0a55625f0040007fab88067e41394009530d04a1d6c1d0007a4eba3fa92d1e000,
    0xb25a25d7cd856c00707e82b20588bb0007bfc299c9cc5688007fcb668b541bb7340764843e64,
    0xc7828007eaced9da3030ec007b58219933d995c008db81aafab0b4f8007584069a8fe9c10007,
    0x14ea6395433100050cd533e71c41400066358238fdd9330008371748333b0cb40072c125ec63,
    0x7818008ad585574cd1728005c37431bd2478c0008a27d93651fd610008c1dbe710f1a5e8009e,
    0xb9e2e8bafc520091edd9cabcb2500008d335bd8787b7e8007ec1ecc98d1544a00885128c4c20,
    0xc6a007af5473b5e67a040093971520252ac60008fd05567ef01cd8009f2019c583daa00007e0,
    0xe5be6d793fa008439640dc6f301c008684957238a204400713d17116cb4d70007de6ebb872cf,
    0xa000712229cc475ee90009fc304929e0cf40007c49da5acbe27a000751be77fa3927f0008175,
    0x986b7c030c007b6968e46347c84006c25f946f9096b000698b8a4e73bd5a0007980a7873a53c,
    0x8009224f2b99b1c5c0007cda9dba342ec2a00739e556de5c4e580076212bd34c515500079bd4,
    0xa30cc96a0005fee3128d6ccd00006ee529ef5f09490006b77ff39ba34ee0007a4672439b1613,
    0x009fd3f1274e0df40006d9d96cc0b7d7100070609a2e6b0898000740a91226a9a7f80064a1b9,
    0xfc21e170007d6b86a5563453a008c6c86f76991bd800a907d20bfb0dd00009c04a206b32fa50,
    0x08aa9dd094fe43f800947d0382d934cb000780ccde7560eaa4006e6342ea795b6d0006f4ecdd,
    0x6a725d6008538f20209ecd680067244e37d419510007d018cf637d8ccc008e3b648f111f5580,
    0x7c5872bc8f59de2008dea416830a3ad800a2f25fa08be5ac0007b22fc9d3723348007c90f309,
    0xa0957c0085f713252f84f08008fd8a20d626aba0007b69eaa21906a9400618b7bc4675844000,
    0xbefc135199a848007a77984bacdaf04007d4bbc2f03707020090b8b41b950119000845a8d416,
    0x31e10006c595c2215f3480008c2a3187de31ff0007bb0c003a44fe80008bc090a09eed520007,
    0xc2e1d5f908e1800797a40ef990caec008c2a8787bfcdf00008209c364fb8ff60006ba9b7b961,
    0x7b10008670f26283f1d4c00718de7f0144d2500085fd37e28fd39dc005291ea4be159f400076,
    0x5cdf97891ab0007d8f2c29c8b2b62008b25c9838112998007810b730c4752640096d5dbf5a8d,
    0xad0009cbf54ee43a66800094e11f5b3f22ae00086af67cb3b69f28006d81d880f7cc4500085e,
    0x0313b3b561200885d104e2a66b780077df220c0fbcb40008889f94a9427dd8006d007d2477c5,
    0xe0008477197be5b1004008bf0710318e4378008ae9e12971d7a78007c080e453ebe3c0007d25,
    0x6e878b0930008175adc1a117659008474384c0ba0d2000860e2c816a940a800876cbf0b79e26,
    0xa0088ca15f17c4f5680091b62499abbffb0008cb3f2eb8e32548008b5bbf5c3b0b348008ea65,
    0x0f55db9ec007214aa27132de1800783104e1b891850006e918776233cb00007d8d519fde9eb7,
    0x00760b0cef9155230006f9e1b1e8a29230007cda5ff9538aeaa00779864b7dbe44f8007b394f,
    0xe0f8c2e000903485308d14c70007988ee3ac57af1800836aee138864ede0085b07c28a86fed8,
    0x086338c3e36884f400797974a109344a400991e06a820b4fa00087b961173a1238c0065bdef0,
    0x286dde00071603259e582018008d0cb637451f020007fa5d38fb541d7ec0771140993f91c080,
    0x940b41b17cb21d0008ce8e52c65c77700086a6bf341ff31c0008883fca0f2f816800771f743a,
    0x2b4e88e0799d56825c92358008275e37340a5d3c006ec2044faea87d0008c51511bcb002d800,
    0xe14028f8bfc88a00827ba6ca3d1291e00846b21f8807fb1c008ead84262fa5778007fe801f7d,
    0x83ca80099dd7ed564d214000794156dbf322dc0006d2459ea0e6a39000865a30136447e48007,
    0x48528e57d51f00070a69551ec54640006c57d4f8ecfd5f000709c7dddd8e99c00089dba3a064,
    0x3620008fa6d279ea44600007a19de72d4ab0dc0076e7163e2ecae500086ee38182d9fed80089,
    0xd54f3f64bc68007d152f9f5cea01c0085a3fa6d2c86964007f10c16e54837c1007b1639fb7bb,
    0x970009550094db2c33900087824df215a690c006eb438e33dc217000796497b0808e10800743,
    0x7c0f2d7a872007dc14b95a9b2f9e007594f370aa07d70008d00f05c70184a80083e87ff00b14,
    0x30008aef3b9b86a18680087aaf4b024bed5c00816c1b01904283f00691c18936deb160007dd2,
    0xbfa737b54e008b7922a9ba023a80089e1e67f0d1f8400068700724e915ca000950ceb405fa43,
    0xc004def41387ac3cc000836e291c1031e3e00695e13321dabd800089a3ac65b534880007ca3b,
    0x4e7cf7d80005571fe853dc99e00098abfa96edcf5600099a001005b661200084fe30d8e747e1,
    0x008629203a82ec2e0005c8ab7a54239960007f999126856d8b0407a678af568e1b8400888f78,
    0xe6ba2f0c006f97ceca395fae0007cf0c52f3117e8c0085fb738b72387440077928d540048718,
    0x07893699749523e4008912d91d38e66980080de7de6b0e97d880898567258247698008554c5a,
    0xa9f65b00090aa5d3188106a00076b9bef0a0c0cd8007a45bedb4c52b8000948e8b4ca8e89c00,
    0x6fd46a36ca42ce0009984a72a8e98f4000838731fdceeccb4008ae3f72a52feb900091d68eb2,
    0xfe9f2800895f1e0937e3ab80080dfae49ef73b20009b6e283ab93ad400092344e4ee0dcd5000,
    0x4a8271f03595b0008e54d0047422ad80072424fff8333130009318e972ac97b200089c157165,
    0x17f9000842464c87214928007af0e1bdff1a5c80076dbcad72ebeaa800a13aff364659b00006,
    0x6bb35258db660007a8f7e4049263b00080fbda411ea136e007bcb2ec6e03aed8006bfeb8c7fa,
    0xcca200729f11a76c34b680092dc8ef977c8820007e50543c970d3260063696da9e5cc7e000a0,
    0x3518b07b1b500073247690a02b88800621f5a38fc31330009f17b1e76738c80007dc6db30419,
    0x5d80075aac6ed4dfec8800877c21ba5a0c5dc007224566afc062c8008a128503da86948006cd,
    0x01434eaf3b00060f2ded680bb2300082312ad9143866a008db30e5f9457178007d8f415bd767,
    0x57007c579177ed15c5a006e29ac0f69b3df0009dc32ae268ee0a0008b6ae2b171b68b800776d,
    0x208bbf2db0009d6bb5fbe7945400069e0bc67a996ea0007447685cea24518007eb1cecf656b5,
    0x60086e30844295e2d4007b8fe86e3696ec40061db40099ef7b4000697dcd363b8e3b00097a74,
    0x185fdd240006d59ea889413ff0007c0e38238fc265e00921177c09d53c80007df07d6408b527,
    0x0080e564c4d13a8f7807bca1a8a0076294007b74d81cabd7090007c416d2c4ef737400985170,
    0xaacb4630008a116367bb2ad78008b8a2a759cc6010008aeec8f937789c8005ae6c2ca5664a00,
    0x072ed575ca6d96a00085c1a90b755be5800690a9902aabc0a0008196c46edaba7c600796c078,
    0xffe6b100082c33646c799a02005f0dfc7c49edc8000712fd3a059d0d00009511fe50c9a00000,
    0x7bcc120f202b08c009df51222bab15f0008cdd50f67170fa80071278a4db38abc00065e93df8,
    0xb77c48007feb29b42057934709bb02225688c2a00069f6514fb573b30008f8be23b2a76fa000,
    0x1a0dbb86d7a3840080ce1f859274efe00730d50fca85220000865c5703d3427b80073859c4b3,
    0x134a00081fbe0a62e7512d0079fe9b2e8c772e00076ef3e2d581e1280092b4bae43533c40008,
    0x8cff27d90d7e0008f5d34b7c4ec1c0009950a2beb475b40009e545d8341364b0008e5ba55722,
    0xc5800076f60bbd315e3b00074d5d2ae032d0d0008c0e028529f8c38008f283e00495dc300061,
    0x2d74f6702b800063f175127a55350008173387595fd1f3007fcc9b8b1a867eb406a17f5ced2a,
    0x96800901f9b8827a87d0007b3668f20fe6244006f91e225bec45d00082cab0dd525ccec0066b,
    0x014d0bec770006e0166ea7f851b00089da32b7c68b8d800820874974ab9eac00792f77bc77ce,
    0xe8006892259f3c521c000792e47c0986517800945e10dae5ab0c0008eac184b97236000093e2,
    0x7171f094d4009b830b03b16800000a8a55991f3f79200089e582130ef413800882db93383c5a,
    0x0006d11e21e75ae740009fe4523d678ad300089ef921f92b21b8006fe43d7f044089000788fd,
    0x68d35051000831ce72fd464ef8009d651e2d50b4b20009bf68f0506c8cf000738310beca5476,
    0x007fbef83db2f726180a275b1564652d400092adbf60c930f20007751a1910f1962000661389,
    0xa124d5a80071735cd26d042780082516c3a19b4cbe0084e274f3f2ce0640075e599c986800d0,
    0x084088a50d64c510008a0312a46dec588008035bdc3e19be19c0963f4287cdf69300089785b8,
    0x8cc65800083b57733d3e4b2600932a8015d667e40007c94fff4504ea460087629320ddc2da80,
    0x55263d132368620007eceeb0eab1b211008f511735dc55240008a165a374f8d0980077f45c85,
    0x6137400099e5e958bb83e20007a69e813f28823c0073fd6ad263a6b08009a2dbab5345bdd000,
    0x1a711a857cbc70006afe701abc5e5600074023ff93f2e120006e4b98240600b400061a1a0ed0,
    0x199d80079389b588c69a580078f83198cd4d07800880396b0a6c7b2800943b604c3eacff0009,
    0xb1db9862bceb8008e8c57be1d012800087b9a754757e42400779acd7f4d109a00076e92c35d2,
    0x3c400071f1e6f3b6eaad8008136a0d478a19e30059d3cafafdf65a000851a95ad34cf7000078,
    0x2c6c8a903128008ec1d003a5d8bf8005cdb27a8c57eb600092f92303995cef000919397870d3,
    0x8200072e547e1a96e2f0008b06d44f8acfaa0008e604267375d3f0007d607f6116bc7260076f,
    0xb6d8271c1180066db3f8fe968ec00094079975fa48b10007886dc1f528000800851efe51dbfe,
    0x57007a313c525e327ac008a5d2b8920a6a0800563f62b4373f620009641edf8c6591e0008d1d,
    0x5fd386f3a8008eda3c73185eba0008b12d4a02ee1360007badc6bb0bc42300081c377e1fee2b,
    0xc007c5d03f42721416007dc1465dfe670ba007936c45fe3a2660007e4c9ff421b236000897d7,
    0x091203790007006f7fa6725ae800895a712b8c698980072ce74f60468b38008347ef06e5db14,
    0x008b9e01f0620e8600090681e5aa5686c00086a507219c9cecc00841fdb4838186a800684062,
    0x5515ee90006dc0100a06a2c7000704fe19b4c76778008b10cac129724780083843fc489f4c0e,
    0x097bd3dbc1e0c8d0008e8095282d35af8008a9561f5a870a300080b735450b4e5a9006c8d5c7,
    0x7653763008410a9d6655ed8800756b68de442cef8007df2ed85d709e5a00854467d5d7539c40,
    0x7045d95e5c5f880006fe7fa356c5bab000752ae9dcbf945f8007acc45339f70dc8007ef8669a,
    0x9d94120066bba222727e33000791e6d13e53a3f800840287d85ca3f9c0097805443f77379000,
    0xc6ab90e25ca048007caf84c969f8ed8007ea0302e3fb95a100780b0f18eb37794008315f8492,
    0x346fa009165303a8ab2c1000a5f7257affb5e00007d26e7fc736218c007ee6deb6f03fc66008,
    0x3b210fbe0c258807113de0f185ed380092ddb71649a712000744d035459508480082148dff35,
    0xb5e400870b8199274bfcc007a1fe89778acf2000796b43fcc39ffa4008da430f42eeeb08007f,
    0x7b87aaa6012c008403aee816a877c007ca6676974b10d4007c030fb16fab87400823b4615162,
    0x64800688e794b1cbaa50009511098033b70d000794ad317ece35680094b54be0bfada3000877,
    0xc287dabc8f00086bbd924a5b47200084f378febb65a44006a5a83971af4ab00075adb9b47380,
    0xe800671f12aa7a0bf500087e6c7c3ffe30a0007160cdc1f59f8d80080c40d8968dbf90806a24,
    0x3eef66f668007eac59dd0ca2095006f6cc75a9d9aa70007f4a1ae86cb7cf4008dd5e1d2fbe80,
    0x0007d3a0cd087c6d440076a2faaced376e000680adbadd077230007a2b2a8cd96ba34007bbfd,
    0x816a6c0d8008c836dde41b1dd00095f848e33b99c400084d678e44eee2700068449955600df8,
    0x007f5e20c2ff3c55b0089b1bd4b3f2d0c800788c8647501ac4000757d79d978e8ec000733935,
    0xe4d6655c006f4becfb26bab500088f2e194fbd1a180067cee93997c02b00070e10c6eba29eb8,
    0x0693a2fa98870950009181b9f9d394070008b8c5fc3861a60800812b4c688b5fa5f0087f3c2e,
    0x601f9e8007faf1b9d432c847807b14aa2e7e9392c008b014a625a93a7800905b5b147bc17f00,
    0x69794295ec407d00097577a6af65a720006493474dd41d020008a91694c373c9280087d07ad3,
    0xa8d230007c69e8e39b4edbe0067d084a3e3c8e10008f4490849ec5860005eb38d40a21a42000,
    0x50dec8226ea064008a5a92236123ea00078ed19dfddbf7c400634961f696fdaa0006a749e1de,
    0x1b010007f9fa4c92f934d8407e6bbea31396c5f007fbf7589e26375a8075b86d6025f6ef8008,
    0x91c93a14066d800860a4fae485d90c006ead4021bfea8800085c67367b34afb00078dd06f2d6,
    0xf03000703c015fe8a6d40007928aa72fbdb108008b4af93bb8a4380007aa1635eea7b370007a,
    0x73ed0296a0f00079545055e5a8dd8005aaadadae28bfe0008ab45bad8b459a0006cc05ade35e,
    0xbb80062cf23f839dcfc0008ab94d991caffc000801c8531000781af06a002983db158300082e,
    0xc8e6c1c5cf000708ae0891c24a08007586a88ae9c9888009c323a1abe12bc00071cdf8bc6283,
    0xa800833a3e5496307da0074bdc76455525a00082f98721722253a00790d0719ec93400007246,
    0x0522ecb1f300732f960126bca50008bff493c0c398d8008560ffdd881560c008b3b32cdce2ac,
    0x20080d05349c0cd3bb00836732e865ed0aa0080e11b847269afb808aa0ed12f9aaeb0007e094,
    0x623d2b8c80082c8b385c7b6f9c006e68e4fc12ad5d0007f9af08f38be5578082ed3a383835a3,
    0x00606d8588a8b958000826f364d62b5cec007786cc7553ac7400096f8f3417b9c7800078b6d1,
    0xbaffec62007d9526e5dd9281e0093b50ec63d729000072bea86f6c39a500075753b8b6c26568,
    0x0981b21cf56666200084cfbaae2f1eb48007780913fe0732480084ce04228572268008251bda,
    0x81729b00082dd60f3d9f87060072a5cb64ba14df8008b6e2ec163e1a380089df33bf72287400,
    0x70684fd77dc4a78009f939230e00b0000081e0f45e6fb311d007584b4273334df00084b18c08,
    0xa5a300006f54af8b407b1100085a1c75d45b6e60006ed1a0ef5f5d2900088d8dfa2e566ac800,
    0xc70bf02e327c1e009576080990763b0008053e0634372c3d806bb4785d5b7c280008972140c9,
    0x1a60400a7bacc6d39b05e000a094abb1cc93f200096779042a7f3f00007e4c3209485b58a007,
    0x90368272796500086f832edb5ea4d8008a6b20d0d54e0b8006f1e291824e79200079b3b5ee8b,
    0x89a000829783a234bdcba0098385af94ed9300008292b7521659e460082674621d9f9df80094,
    0x9676eb3c63c800777169d9863e610007f49d7cfc5411348087c0957821b6c04008aea6f61668,
    0x09000768c7da66740ba000667a6d1130e56300094748d3f54527b00060d2d562c89a330007a9,
    0xd38b4fce1f8008f1b8bfa1ede7a00080ecf4dc2cbe4be809349d81a25ffdc0007a131fc0de5b,
    0x50008e0a5dedd5fab780088c1b794ce2c5d80087cabe39710efac009abdecd799dd400007b4a,
    0x9b7139b42ba073bb7a1c981c20800836e481eef6d996007fd4e4c5eea4a99606848c6d5573e8,
    0x000848cf111a05805800637157f6048c2e000838415e17daa02a00747563d75b86db8008031d,
    0x73640af980080ae226f91d05d880606e23afb099f40008f078945b364ee00075d1d2e29fba29,
    0x008b63c8e7059a7b000773865d3bfebda00066b846458550f70008f0c627cf809220007c7727,
    0x4fbe64a000654da73423d569000849dfbe605106780085b38b64a9632a00089ce54925a44000,
    0x08f42a1fbea471f8007d398b1918c8ae60062f6f5b40d49740008373d199fc71908008d9df36,
    0x109620c0079a609c8a433984008e857ef229483c80075e641c07905f10008b42b6c52a9be380,
    0x8d1ea5d48c3e70800758e87c5b1d76e800863fccb67dbdf94008edb451b0ac9aa000861c8101,
    0x07a8a800900efc2f831e3c0007d15bb651731d6400712202eb94d02200097a38dbac90d9b000,
    0xa35df41b3b58b00091105795bb3ac30007815b676e00d5b400647d6fcc4bd498000725e3dc10,
    0x36d0000891219b6dd91f780080c03700961cd6c808f6410cb2afe4f8007c06f2900192480006,
    0x772026c3cc8500087db14905ee167000897584541c69d90008e698f12e60d0d00099fbd29ecc,
    0x36e00080ea4ae85d558268086753b03dafe7940088ac77fe3e2a7b0008d6d59b8c3aa5d0006d,
    0xfb8b310fa6c00076366c78ae9e2f00094144cd7ae02b400082751eb8d9c27dc006ba9c9b1e6e,
    0x7800089fc1f7dea2d36800923ee69019468400076f8b7b2ba60a38009d688fe9795ca4000632,
    0x256a8597280009bb754248b10cb0007b31611413c2078007f5398d2be65856009a861649ff40,
    0x2400753d4a2ad54bb70008be0a9de15cb2480082b7b1891d86c20008518a740e44655c00a281,
    0x15d6ecd610008cb7d7cbdef9998007b3d86040df5be000702a8c9af1398b00087f7188e398fa,
    0x800a2a68e3c16fa9400096f09c034ecf3e0008808ece8e167a58008d9d20077f382300073caf,
    0x60da1e35000832f4a42280fddc0098f0496ed4439a00086400e8129440900076cb42ca252c40,
    0x006f3e420ef6067f000a997849d622de20009668384922747c0007009f331f4c12b8007589fb,
    0x69bb1eb5a07b15ca90ee62a1400811d5f8e09797470066e93053db83ef00081c4c2c39c34f34,
    0x09d5290e5e098c800086e9efbe07e163c00699f0137f764b7000ae5152d67302580007fcea13,
    0x003ed0c0094cf0ef40bb19f00080c5e0240b3c0b98072e3f8e52c2d6a0007639406cee8fda00,
    0x6864892bca50990005fb70324a39c920008120f781fc78e5a006f965be477c1740007a828f87,
    0x4254c00089e1c195a401b58005d374f84faf19800087574efff201fd000916226d0a1f102000,
    0x447ab711af0170006924f7d4ed9cfb00084cab6bc6fe9e300078d9375604a3e54005b1f8ee89,
    0x00000000000000000000000091efb4ab90b69a0007533ca446a8e05000936e9d9f89e2510009]

def c520N8 : ℕ := ofChunks76 [
    0xf928000084d7ec1e95888d80078d450a1f751870008ac181bc9a6c320008d98572fa3ac2f000,
    0x9c1d71f051d340008d0c2df14b524a0009e9da9ac29d82a0008f66bf83c0ff71800b9936a7e5,
    0xd58800074a7b08649da8f00087f11e1603c2f5000a5d29fd41f2be80008f2807eefc632a8009,
    0x4285c251f65000091fd7410adc9ec000a76a1e113188e40008c323dce30f7db0009e6075a1ed,
    0xa8400088a6f59abb59aa0009b9c0d4929bfac0008a6a6442363fc2000959a87959430c0000a7,
    0x1cdf978a89e0008871fdc93cb7ef0006f85b491a505340009258687b649921000aaf4c6187de,
    0x1c00083d13209ee32fb0009a5da31d6fb6a000090f284c540fc200006d2574e81416fc00075e,
    0x2a10264d6c0009b3f7ecf732e7800077310c697a2184000865e1501d77bc88007abd0fd43c1d,
    0x30008d4fa11e879a39000a0c1659f474bd40008603d05e91a13f8007679af9c2d39ee0009931,
    0xd81037bb2000a34cc59afc71c80009a07834b6c222400096ad0415c8fad40009026fa22eedce,
    0x0008e413575b6bbf800094d636f2913c130008a337fd3094da200079f57c3cfcb7390009c5c6,
    0x7770fb3100092875c7604bab4000a03b99fcf0d42400083dec3181887838006d2223960a6864,
    0x00b1b476cada0bdc0009abfe9d8df2dec0008d75ef270f106280093746f09a1aef3000787064,
    0x625eaf78006dd21ab49cccfe0008dc7a082d5524380087eb42832a59d700077bedf8111539f0,
    0x0a7cfc0ffeaa82c000614ff82c6dd762000b3c507689e58860008682612ed2a9430008f801f7,
    0x3c9a7c00084a08fdea1c5dd80092595788369e370009bc2ef70cc8537000a020de7954b6dc00,
    0xa138458b065236000967135c9503b300009b1dbe91639ec800081d3fae2b24e5f800a1b0aa8d,
    0x328fe00094e4161e1ce87600094d54e3c27625300086b4465127543900086748bbaa361fc000,
    0x74d61e4172ede8009d7b5ef2093d50000925143ca02645d00097c43230e6a4e4000a5fbe5b2b,
    0xf01900074b0c27f1f799200095690d37d465cf0007fa3ab6a49d765000a17cd0195f74000008,
    0xc41dff4b638280073e17664dc99fc0008c4bd6713db87780076e99f7de6944f00093f174c701,
    0x1f800083b3d59f59774100086afa0b490108980098e05b62f3baaa0008394d20199849200085,
    0xd2bc86158f20007823dee620329d000a17ee54990022c000aab4ec080fe12c000a38b5408ab7,
    0x120008c93eb240b01b70009bc0573eb81d260008ed949e29aca308007c57e8fcfc4040000a46,
    0xd30598fadc0009694e4731ce99a00096d881a74475ce00092019100c274250008d97fdb72e86,
    0x18008aec179f678e490009973abcc66dfde0009496e3fdc65f730008d308b7a31cd360007460,
    0x22e9de4670008e0e3be12b6fe70008456ce35ccd34e8007a1e5c40dbc9a100084c1f05a7dabf,
    0x0009513ca8e7635f200076429b7cb3448200095e148f87f21440009983e6f745238400085d63,
    0xb3ea6cb300089158d48f8e2870008f93792e0df3288008bd19aca4c14a900098aad280f230d7,
    0x00868d74eefed1640007c6726b785227900096cda67be76cb10007fe329e4bc9a390008b4e49,
    0x5c7769c800883c67db9147370008a1a6566c548b48008210b896354bef0008161bd0c4ab6ab0,
    0x09a1c0ab7cd2a9600091c6d02ad1af6100076a44267dca92a0008bf1ef4c59eebf80084ccf4c,
    0x96f44e00093d5058b672211000b3b6429a0b0e180007963f09a44cff40007b786d8182c8b000,
    0x8743f1d4221e0200090988a15376b9a00097e1e29e6761220008cd2e9d8a48d7280090d9e0d3,
    0xc335800083f90a9a17601e8008ae2a257fa01860008a7cfd9ac53424800891493f9bd66eb000,
    0x59a892d60ce780009975c13bc67a2e00083a184ad1b9bfa00091804ba92329e1000a12f85143,
    0x52010008cd18ce221e171800817f78bd178d56000966a031ac7f1c30008dcc795a638d778009,
    0x6454500a5ae380092b03ac624775c0008ab29754387f650007d0f1cca37cb290008003765b54,
    0x002000a2bb0d9bb722e800088da2a0a10912780072b67e20db41bd0008d32f63f81e88b00081,
    0x71a678320d4000988e02e604a76b0008c7b728d20fa9f0007ce9c01d38c58e0007d93ccec75f,
    0x03800b072608f1916160009abe460d954f420009cd0756d91e5c8000ba75b8d696aa72000a3c,
    0x0d2c70ea31000794751c21f8273000940d65001a37480008c961eef0b5eb60008b2f69ff44f2,
    0xa000801cda6a2b363e80088ac851f69c00d80086ccd143a94241800ace89bf984132400072a2,
    0x7c6f43d8e000903c3afd1e23840008c579bbe1d5ec380099754fc3b339e9000924e1f112ee24,
    0x0008285d4cf32a2830008f07bdb92611cf8008704bcb512d3ae0008fd3e4da5628c20009f38e,
    0x66562298000683732b5df6276000a64d84c48868ac00081a21892ad901e80091a155663157a4,
    0x0097b59b178e806d0008c3fa3461170d1000afcb5c8c749ecc00088310a5ddd3401000966bb3,
    0xe075bd0000aff5681a79d08c000984a380643d6e600093f677582391570009630b7bb0c09490,
    0x09a67c4343a1afc0009954c729bdf26e000aca4ee30daf43000093248da5f92e740009ed7137,
    0x23cde30008b9bbaa17b91400009b9b9c2095ceca0009293e48a9ecb6e000a3eba4bbad9b9c00,
    0x870ec1640cfe1b800b08df088a013ea00085c3a74a5a0c9900090a2604ebe621c00091ff541a,
    0xad9c8800a2b0a55625f0040008fab88067e4139000a530d04a1d6c1c0008a4eba3fa92d1e000,
    0xb25a25d7cd857000807e82b20588bb0008bfc299c9cc5688008fcb668b541bb7000864843e64,
    0xc7828008eaced9da3030f0008b58219933d9960009db81aafab0b500008584069a8fe9c10008,
    0x14ea6395433000060cd533e71c41400076358238fdd9330009371748333b0cb00082c125ec63,
    0x7820009ad585574cd1720006c37431bd2478c0009a27d93651fd610009c1dbe710f1a5e000ae,
    0xb9e2e8bafc5000a1edd9cabcb2500009d335bd8787b7e0008ec1ecc98d1544800985128c4c20,
    0xc68008af5473b5e67a0000a3971520252ac60009fd05567ef01ce000af2019c583daa00008e0,
    0xe5be6d79400009439640dc6f3020009684957238a204000813d17116cb4d70008de6ebb872cf,
    0xa000812229cc475ee9000afc304929e0cf40008c49da5acbe27a000851be77fa3927f0009175,
    0x986b7c0310008b6968e46347c80007c25f946f9096b000798b8a4e73bd5a0008980a7873a53c,
    0x800a224f2b99b1c5c0008cda9dba342ec2800839e556de5c4e580086212bd34c515500089bd4,
    0xa30cc96a0006fee3128d6ccd00007ee529ef5f09490007b77ff39ba34ee0008a4672439b1613,
    0x00afd3f1274e0df40007d9d96cc0b7d7100080609a2e6b0898000840a91226a9a7f80074a1b9,
    0xfc21e170008d6b86a55634538009c6c86f76991be000b907d20bfb0dd0000ac04a206b32fa40,
    0x09aa9dd094fe440000a47d0382d934cc000880ccde7560eaa0007e6342ea795b6d0007f4ecdd,
    0x6a725d8009538f20209ecd600077244e37d419510008d018cf637d8cd0009e3b648f111f5600,
    0x8c5872bc8f59de0009dea416830a3ae000b2f25fa08be5ac0008b22fc9d3723348008c90f309,
    0xa095800095f713252f84f00009fd8a20d626aba0008b69eaa21906a9000718b7bc4675844000,
    0xbefc135199a848008a77984bacdaf00008d4bbc2f037070000a0b8b41b950118000945a8d416,
    0x31e10007c595c2215f3480009c2a3187de31ff0008bb0c003a44fe80009bc090a09eed520008,
    0xc2e1d5f908e1800897a40ef990caf0009c2a8787bfcdf00009209c364fb8ff60007ba9b7b961,
    0x7b00009670f26283f1d5000818de7f0144d2500095fd37e28fd39e0006291ea4be159f400086,
    0x5cdf97891ab0008d8f2c29c8b2b60009b25c98381129a0008810b730c47526000a6d5dbf5a8d,
    0xad000acbf54ee43a668000a4e11f5b3f22ae00096af67cb3b69f20007d81d880f7cc4500095e,
    0x0313b3b561000985d104e2a66b800087df220c0fbcb40009889f94a9427de0007d007d2477c5,
    0xe0009477197be5b1000009bf0710318e4380009ae9e12971d7a80008c080e453ebe3c0008d25,
    0x6e878b0930009175adc1a117660009474384c0ba0d2000960e2c816a940a000976cbf0b79e26,
    0x80098ca15f17c4f56000a1b62499abbffc0009cb3f2eb8e32540009b5bbf5c3b0b340009ea65,
    0x0f55db9f0008214aa27132de1800883104e1b891850007e918776233cb00008d8d519fde9eb7,
    0x00860b0cef9155230007f9e1b1e8a29230008cda5ff9538aea800879864b7dbe44f8008b394f,
    0xe0f8c2e000a03485308d14c80008988ee3ac57af1800936aee138864ee00095b07c28a86fee0,
    0x096338c3e36884f000897974a109344a000a91e06a820b4fa00097b961173a123900075bdef0,
    0x286dde00081603259e582018009d0cb637451f020008fa5d38fb541d8000871140993f91c080,
    0xa40b41b17cb21c0009ce8e52c65c77700096a6bf341ff31c0009883fca0f2f816000871f743a,
    0x2b4e8800899d56825c92358009275e37340a5d40007ec2044faea87d0009c51511bcb002e000,
    0xe14028f8bfc88800927ba6ca3d1292000946b21f8807fb20009ead84262fa5780008fe801f7d,
    0x83ca000a9dd7ed564d214000894156dbf322dc0007d2459ea0e6a39000965a30136447e40008,
    0x48528e57d51f00080a69551ec54640007c57d4f8ecfd5f000809c7dddd8e99c00099dba3a064,
    0x3620009fa6d279ea44600008a19de72d4ab0e00086e7163e2ecae500096ee38182d9fee00099,
    0xd54f3f64bc68008d152f9f5cea0200095a3fa6d2c86960008f10c16e54837c0008b1639fb7bb,
    0x97000a550094db2c33800097824df215a6910007eb438e33dc217000896497b0808e10800843,
    0x7c0f2d7a870008dc14b95a9b2fa0008594f370aa07d70009d00f05c70184a00093e87ff00b14,
    0x40009aef3b9b86a18600097aaf4b024bed6000916c1b01904284000791c18936deb160008dd2,
    0xbfa737b550009b7922a9ba023a00099e1e67f0d1f8400078700724e915ca000a50ceb405fa43,
    0x0005def41387ac3cc000936e291c1031e4000795e13321dabd800099a3ac65b534880008ca3b,
    0x4e7cf7d80006571fe853dc99e000a8abfa96edcf56000a9a001005b661200094fe30d8e747e2,
    0x009629203a82ec2e0006c8ab7a54239960008f999126856d8b0008a678af568e1b8000988f78,
    0xe6ba2f10007f97ceca395fae0008cf0c52f3117e900095fb738b72387400087928d540048718,
    0x08893699749523e0009912d91d38e66a00090de7de6b0e97e0009985672582476a0009554c5a,
    0xa9f65c000a0aa5d3188106a00086b9bef0a0c0cd8008a45bedb4c52b8000a48e8b4ca8e89c00,
    0x7fd46a36ca42ce000a984a72a8e98f4000938731fdceeccb0009ae3f72a52feb9000a1d68eb2,
    0xfe9f2000995f1e0937e3ac00090dfae49ef73b2000ab6e283ab93ad4000a2344e4ee0dcd4000,
    0x4a8271f03595b0009e54d0047422ae00082424fff833313000a318e972ac97b200099c157165,
    0x17f9000942464c87214920008af0e1bdff1a5c80086dbcad72ebeaa800b13aff364659b00007,
    0x6bb35258db660008a8f7e4049263b00090fbda411ea1370008bcb2ec6e03aed8007bfeb8c7fa,
    0xcca000829f11a76c34b6800a2dc8ef977c8820008e50543c970d3280073696da9e5cc7e000b0,
    0x3518b07b1b500083247690a02b88800721f5a38fc3133000af17b1e76738c80008dc6db30419,
    0x5d80085aac6ed4dfec8800977c21ba5a0c5e0008224566afc062c8009a128503da86940007cd,
    0x01434eaf3b00070f2ded680bb2300092312ad91438670009db30e5f9457180008d8f415bd767,
    0x58008c579177ed15c58007e29ac0f69b3df000adc32ae268ee0a0009b6ae2b171b68c000876d,
    0x208bbf2dc000ad6bb5fbe7945400079e0bc67a996ea0008447685cea24518008eb1cecf656b5,
    0x80096e30844295e2d0008b8fe86e3696ec00071db40099ef7b4000797dcd363b8e3b000a7a74,
    0x185fdd240007d59ea889413ff0008c0e38238fc266000a21177c09d53c80008df07d6408b527,
    0x0090e564c4d13a8f0008bca1a8a0076290008b74d81cabd7090008c416d2c4ef737000a85170,
    0xaacb4630009a116367bb2ad80009b8a2a759cc6010009aeec8f937789c0006ae6c2ca5664a00,
    0x082ed575ca6d96a00095c1a90b755be6000790a9902aabc0a0009196c46edaba7c000896c078,
    0xffe6b100092c33646c799a00006f0dfc7c49edc8000812fd3a059d0d0000a511fe50c9a00000,
    0x8bcc120f202b09000adf51222bab1600009cdd50f67170fa00081278a4db38abc00075e93df8,
    0xb77c48008feb29b4205793800abb02225688c2a00079f6514fb573b30009f8be23b2a76fa000,
    0x1a0dbb86d7a3800090ce1f859274f0000830d50fca85220000965c5703d3427c00083859c4b3,
    0x134a00091fbe0a62e751300089fe9b2e8c772e00086ef3e2d581e12800a2b4bae43533c40009,
    0x8cff27d90d7e0009f5d34b7c4ec1c000a950a2beb475b4000ae545d8341364c0009e5ba55722,
    0xc5800086f60bbd315e3b00084d5d2ae032d0d0009c0e028529f8c40009f283e00495dc300071,
    0x2d74f6702b800073f175127a55350009173387595fd1f0008fcc9b8b1a867e8007a17f5ced2a,
    0x96800a01f9b8827a87c0008b3668f20fe6240007f91e225bec45d00092cab0dd525ccf00076b,
    0x014d0bec780007e0166ea7f851b00099da32b7c68b8e000920874974ab9eb000892f77bc77ce,
    0xe0007892259f3c521c000892e47c0986517800a45e10dae5ab0c0009eac184b972360000a3e2,
    0x7171f094d000ab830b03b16800000b8a55991f3f79200099e582130ef414000982db93383c5a,
    0x0007d11e21e75ae74000afe4523d678ad400099ef921f92b21c0007fe43d7f044089000888fd,
    0x68d35051000931ce72fd464f0000ad651e2d50b4b2000abf68f0506c8d0000838310beca5476,
    0x008fbef83db2f726000b275b1564652d4000a2adbf60c930f20008751a1910f1962000761389,
    0xa124d5a00081735cd26d042780092516c3a19b4cc00094e274f3f2ce0600085e599c986800d0,
    0x094088a50d64c510009a0312a46dec580009035bdc3e19be2000a63f4287cdf69400099785b8,
    0x8cc65800093b57733d3e4b2000a32a8015d667e40008c94fff4504ea480097629320ddc2da00,
    0x65263d132368620008eceeb0eab1b210009f511735dc55240009a165a374f8d0a00087f45c85,
    0x61374000a9e5e958bb83e20008a69e813f2882400083fd6ad263a6b0800aa2dbab5345bdc000,
    0x1a711a857cbc80007afe701abc5e5600084023ff93f2e120007e4b98240600b400071a1a0ed0,
    0x199d80089389b588c69a580088f83198cd4d07800980396b0a6c7b2000a43b604c3ead00000a,
    0xb1db9862bceb8009e8c57be1d012800097b9a754757e42000879acd7f4d109a00086e92c35d2,
    0x3c400081f1e6f3b6eaad8009136a0d478a19e00069d3cafafdf65a000951a95ad34cf7000088,
    0x2c6c8a903128009ec1d003a5d8c00006cdb27a8c57eb6000a2f92303995cf0000a19397870d3,
    0x8200082e547e1a96e2f0009b06d44f8acfaa0009e604267375d3f0008d607f6116bc7280086f,
    0xb6d8271c1200076db3f8fe968ec000a4079975fa48b00008886dc1f528000800951efe51dbfe,
    0x50008a313c525e327b0009a5d2b8920a6a0000663f62b4373f62000a641edf8c6591e0009d1d,
    0x5fd386f3a0009eda3c73185eba0009b12d4a02ee1360008badc6bb0bc42300091c377e1fee2b,
    0x0008c5d03f42721418008dc1465dfe670b8008936c45fe3a2660008e4c9ff421b236000997d7,
    0x091203790008006f7fa6725ae800995a712b8c698a00082ce74f60468b38009347ef06e5db15,
    0x009b9e01f0620e86000a0681e5aa5686c00096a507219c9ced000941fdb4838186a000784062,
    0x5515ee90007dc0100a06a2c7000804fe19b4c76778009b10cac129724800093843fc489f4c10,
    0x0a7bd3dbc1e0c8c0009e8095282d35b00009a9561f5a870a300090b735450b4e5b0007c8d5c7,
    0x7653760009410a9d6655ed8000856b68de442cef8008df2ed85d709e5800954467d5d7539c00,
    0x8045d95e5c5f880007fe7fa356c5bab000852ae9dcbf945f8008acc45339f70dc8008ef8669a,
    0x9d94100076bba222727e33000891e6d13e53a3f800940287d85ca3fa000a7805443f77378000,
    0xc6ab90e25ca040008caf84c969f8ed8008ea0302e3fb95a000880b0f18eb37790009315f8492,
    0x3470000a165303a8ab2c0000b5f7257affb5e00008d26e7fc7362190008ee6deb6f03fc68009,
    0x3b210fbe0c258008113de0f185ed3800a2ddb71649a712000844d035459508480092148dff35,
    0xb5e000970b8199274bfd0008a1fe89778acf2000896b43fcc39ffa0009da430f42eeeb00008f,
    0x7b87aaa60130009403aee816a8780008ca6676974b10d0008c030fb16fab87000923b4615162,
    0x64800788e794b1cbaa5000a511098033b70c000894ad317ece356800a4b54be0bfada4000977,
    0xc287dabc8f00096bbd924a5b47200094f378febb65a40007a5a83971af4ab00085adb9b47380,
    0xe000771f12aa7a0bf500097e6c7c3ffe30a0008160cdc1f59f8d80090c40d8968dbf90007a24,
    0x3eef66f668008eac59dd0ca2098007f6cc75a9d9aa70008f4a1ae86cb7cf0009dd5e1d2fbe80,
    0x0008d3a0cd087c6d400086a2faaced376e000780adbadd077230008a2b2a8cd96ba30008bbfd,
    0x816a6c0d8009c836dde41b1dd000a5f848e33b99c400094d678e44eee2700078449955600df8,
    0x008f5e20c2ff3c5580099b1bd4b3f2d0c000888c8647501ac4000857d79d978e8ec000833935,
    0xe4d66560007f4becfb26bab500098f2e194fbd1a200077cee93997c02b00080e10c6eba29eb8,
    0x0793a2fa9887095000a181b9f9d394080009b8c5fc3861a60000912b4c688b5fa600097f3c2e,
    0x601f9e0008faf1b9d432c848008b14aa2e7e93930009b014a625a93a8000a05b5b147bc18000,
    0x79794295ec407d000a7577a6af65a720007493474dd41d020009a91694c373c9200097d07ad3,
    0xa8d230008c69e8e39b4edc00077d084a3e3c8e10009f4490849ec5860006eb38d40a21a42000,
    0x50dec8226ea060009a5a92236123ea00088ed19dfddbf7c000734961f696fdaa0007a749e1de,
    0x1b010008f9fa4c92f934d8008e6bbea31396c60008fbf7589e2637580085b86d6025f6ef8009,
    0x91c93a14066d800960a4fae485d910007ead4021bfea8800095c67367b34afb00088dd06f2d6,
    0xf03000803c015fe8a6d40008928aa72fbdb108009b4af93bb8a4380008aa1635eea7b370008a,
    0x73ed0296a0f00089545055e5a8dd8006aaadadae28bfe0009ab45bad8b459a0007cc05ade35e,
    0xbb80072cf23f839dcfc0009ab94d991caffc000901c85310007820007a002983db158300092e,
    0xc8e6c1c5cf000808ae0891c24a08008586a88ae9c988800ac323a1abe12bc00081cdf8bc6283,
    0xa000933a3e5496307e00084bdc76455525a00092f98721722254000890d0719ec93400008246,
    0x0522ecb1f000832f960126bca50009bff493c0c398e0009560ffdd8815610009b3b32cdce2ac,
    0x00090d05349c0cd3c000936732e865ed0b00090e11b847269b00009aa0ed12f9aaeb0008e094,
    0x623d2b8c80092c8b385c7b6fa0007e68e4fc12ad5d0008f9af08f38be5580092ed3a383835a3,
    0x00706d8588a8b958000926f364d62b5cf0008786cc7553ac74000a6f8f3417b9c7800088b6d1,
    0xbaffec60008d9526e5dd9282000a3b50ec63d729000082bea86f6c39a500085753b8b6c26568,
    0x0a81b21cf56666200094cfbaae2f1eb40008780913fe0732480094ce04228572260009251bda,
    0x81729b00092dd60f3d9f87000082a5cb64ba14df8009b6e2ec163e1a400099df33bf72287400,
    0x80684fd77dc4a7800af939230e00b0000091e0f45e6fb3120008584b4273334df00094b18c08,
    0xa5a300007f54af8b407b1100095a1c75d45b6e60007ed1a0ef5f5d2900098d8dfa2e566ac000,
    0xc70bf02e327c2000a576080990763c0009053e0634372c40007bb4785d5b7c280009972140c9,
    0x1a60000b7bacc6d39b05e000b094abb1cc93f2000a6779042a7f3f00008e4c3209485b588008,
    0x90368272796400096f832edb5ea4e0009a6b20d0d54e0c0007f1e291824e79200089b3b5ee8b,
    0x89a000929783a234bdcc000a8385af94ed9300009292b7521659e400092674621d9f9e0000a4,
    0x9676eb3c63c800877169d9863e610008f49d7cfc5411380097c0957821b6c00009aea6f61668,
    0x09000868c7da66740ba000767a6d1130e563000a4748d3f54527c00070d2d562c89a330008a9,
    0xd38b4fce1f8009f1b8bfa1ede7a00090ecf4dc2cbe4c000a349d81a25ffdc0008a131fc0de5b,
    0x50009e0a5dedd5fab800098c1b794ce2c5e00097cabe39710efb000aabdecd799dd400008b4a,
    0x9b7139b4300083bb7a1c981c20800936e481eef6d990008fd4e4c5eea4a98007848c6d5573e8,
    0x000948cf111a05806000737157f6048c2e000938415e17daa03000847563d75b86db8009031d,
    0x73640af980090ae226f91d05e000706e23afb099f40009f078945b364ee00085d1d2e29fba29,
    0x009b63c8e7059a7b000873865d3bfebda00076b846458550f70009f0c627cf809220008c7727,
    0x4fbe64a000754da73423d569000949dfbe605106800095b38b64a9632a00099ce54925a44000,
    0x09f42a1fbea47200008d398b1918c8ae80072f6f5b40d49740009373d199fc71900009d9df36,
    0x10962100089a609c8a433980009e857ef229483c00085e641c07905f10009b42b6c52a9be400,
    0x9d1ea5d48c3e70000858e87c5b1d76e800963fccb67dbdf90009edb451b0ac9aa000961c8101,
    0x07a8a800a00efc2f831e3c0008d15bb651731d6000812202eb94d022000a7a38dbac90d9c000,
    0xa35df41b3b58b000a1105795bb3ac40008815b676e00d5b000747d6fcc4bd498000825e3dc10,
    0x36d0000991219b6dd91f800090c03700961cd70009f6410cb2afe500008c06f2900192480007,
    0x772026c3cc8500097db14905ee167000997584541c69d90009e698f12e60d0d000a9fbd29ecc,
    0x36e00090ea4ae85d558200096753b03dafe7900098ac77fe3e2a7b0009d6d59b8c3aa5d0007d,
    0xfb8b310fa6c00086366c78ae9e2f000a4144cd7ae02b400092751eb8d9c27e0007ba9c9b1e6e,
    0x7800099fc1f7dea2d36000a23ee69019468400086f8b7b2ba60a3800ad688fe9795ca4000732,
    0x256a859728000abb754248b10cc0008b31611413c2078008f5398d2be6585800aa861649ff40,
    0x2000853d4a2ad54bb70009be0a9de15cb2400092b7b1891d86c20009518a740e44656000b281,
    0x15d6ecd610009cb7d7cbdef99a0008b3d86040df5be000802a8c9af1398b00097f7188e398fa,
    0x800b2a68e3c16fa94000a6f09c034ecf3e0009808ece8e167a60009d9d20077f382300083caf,
    0x60da1e35000932f4a42280fde000a8f0496ed4439a00096400e8129440900086cb42ca252c40,
    0x007f3e420ef6067f000b997849d622de2000a668384922747c0008009f331f4c12b8008589fb,
    0x69bb1eb8008b15ca90ee62a1000911d5f8e09797400076e93053db83ef00091c4c2c39c34f30,
    0x0ad5290e5e098c800096e9efbe07e164000799f0137f764b7000be5152d67302580008fcea13,
    0x003ed1000a4cf0ef40bb1a000090c5e0240b3c0c00082e3f8e52c2d6a0008639406cee8fda00,
    0x7864892bca50990006fb70324a39c920009120f781fc78e60007f965be477c1740008a828f87,
    0x4254c00099e1c195a401b60006d374f84faf19800097574efff201fd000a16226d0a1f102000,
    0x447ab711af0180007924f7d4ed9cfb00094cab6bc6fe9e300088d9375604a3e50006b1f8ee89,
    0x000000000000000000000000a1efb4ab90b69a0008533ca446a8e05000a36e9d9f89e250000a]

theorem c519_verdicts (σ τ : Equiv.Perm (Fin 1000)) :
    nnsFSD (col c519xN ∘ σ) (col c519yN ∘ τ) true = .yDom ∧
      nnsFSD (col c519yN ∘ τ) (col c519xN ∘ σ) true = .xDom := by
  have hpair : allLe (digits68 1000 c519xN) (digits68 1000 c519yN) = true := by
    decide +kernel
  have hmin : cnt c519yN (95742758182101884928) = 0 := by decide +kernel
  have hle := col_le hpair
  have hdom : sdDomPairCentred (col c519yN) (col c519xN) 1 true = 1 := by
    rw [sdDomPairCentred_eq_one_iff (by norm_num) (by norm_num) _ _ (Or.inl rfl)]
    simp only [↓reduceIte]
    rw [fsd_iff_lt]
    refine ⟨LPM_zero_le_of_forall_le hle, dval (95742758182101884928), ?_⟩
    rw [LPM_zero_col, LPM_zero_col, hmin]
    refine div_lt_div_of_pos_right ?_ (by norm_num)
    have : 0 < cnt c519xN (95742758182101884928) := by decide +kernel
    exact_mod_cast this
  have hcd := fsdUni_curveDom (by norm_num : 0 < 1000) (by norm_num : 0 < 1000)
  rw [nnsFSD_perm (by norm_num) (by norm_num), nnsFSD_perm (by norm_num) (by norm_num)]
  refine ⟨(nnsFSD_eq_yDom_iff (by norm_num) (by norm_num) _ _ true).mpr
    ((hcd _ _ (Or.inl rfl) true).mp hdom), ?_⟩
  exact (nnsFSD_eq_xDom_iff (by norm_num) (by norm_num) _ _ true).mpr
    ((hcd _ _ (Or.inl rfl) true).mp hdom)

def c520Cols : Fin 8 → Fin 1000 → ℝ :=
  ![col c520N1, col c520N2, col c520N3, col c520N4, col c520N5, col c520N6, col c520N7,
    col c520N8]

theorem c520_effSet :
    (sdDomPairCentred (c520Cols 1) (c520Cols 0) 1 true = 1 ∧
      sdDomPairCentred (c520Cols 3) (c520Cols 2) 1 true = 1 ∧
      sdDomPairCentred (c520Cols 5) (c520Cols 4) 1 true = 1 ∧
      sdDomPairCentred (c520Cols 7) (c520Cols 6) 1 true = 1) ∧
    (∀ i ∈ ({1, 3, 5, 7} : Finset (Fin 8)), ∀ j ∈ ({1, 3, 5, 7} : Finset (Fin 8)),
      sdDomPairCentred (c520Cols i) (c520Cols j) 1 true ≠ 1) ∧
    effSet (domMatrix c520Cols 1 true) univ = {1, 3, 5, 7} := by
  have s1 : sdDomPairCentred (c520Cols 1) (c520Cols 0) 1 true = 1 :=
    fsdUni_of_forall_lt (n := 1000) (by norm_num) (col_lt (by decide +kernel))
  have s3 : sdDomPairCentred (c520Cols 3) (c520Cols 2) 1 true = 1 :=
    fsdUni_of_forall_lt (n := 1000) (by norm_num) (col_lt (by decide +kernel))
  have s5 : sdDomPairCentred (c520Cols 5) (c520Cols 4) 1 true = 1 :=
    fsdUni_of_forall_lt (n := 1000) (by norm_num) (col_lt (by decide +kernel))
  have s7 : sdDomPairCentred (c520Cols 7) (c520Cols 6) 1 true = 1 :=
    fsdUni_of_forall_lt (n := 1000) (by norm_num) (col_lt (by decide +kernel))
  have hanti : ∀ i ∈ ({1, 3, 5, 7} : Finset (Fin 8)), ∀ j ∈ ({1, 3, 5, 7} : Finset (Fin 8)),
      sdDomPairCentred (c520Cols i) (c520Cols j) 1 true ≠ 1 := by
    intro i hi j hj
    simp only [mem_insert, mem_singleton] at hi hj
    rcases hi with rfl | rfl | rfl | rfl <;> rcases hj with rfl | rfl | rfl | rfl
    · exact (fsdUni_self (n := 1000) (by norm_num) _).symm ▸ zero_ne_one
    · exact not_fsd_col (N := c520N2) (M := c520N4) (c := 173634949939558014976) (by decide +kernel)
    · exact not_fsd_col (N := c520N2) (M := c520N6) (c := 145072786777335156736) (by decide +kernel)
    · exact not_fsd_col (N := c520N2) (M := c520N8) (c := 174523099954693971968) (by decide +kernel)
    · exact not_fsd_col (N := c520N4) (M := c520N2) (c := 152433226522401349632) (by decide +kernel)
    · exact (fsdUni_self (n := 1000) (by norm_num) _).symm ▸ zero_ne_one
    · exact not_fsd_col (N := c520N4) (M := c520N6) (c := 152545171685730668544) (by decide +kernel)
    · exact not_fsd_col (N := c520N4) (M := c520N8) (c := 152433226522401349632) (by decide +kernel)
    · exact not_fsd_col (N := c520N6) (M := c520N2) (c := 163740347441912875008) (by decide +kernel)
    · exact not_fsd_col (N := c520N6) (M := c520N4) (c := 173635158385156161536) (by decide +kernel)
    · exact (fsdUni_self (n := 1000) (by norm_num) _).symm ▸ zero_ne_one
    · exact not_fsd_col (N := c520N6) (M := c520N8) (c := 157813139893371787264) (by decide +kernel)
    · exact not_fsd_col (N := c520N8) (M := c520N2) (c := 163778957461474406400) (by decide +kernel)
    · exact not_fsd_col (N := c520N8) (M := c520N4) (c := 179414436066368184320) (by decide +kernel)
    · exact not_fsd_col (N := c520N8) (M := c520N6) (c := 141181622740032040960) (by decide +kernel)
    · exact (fsdUni_self (n := 1000) (by norm_num) _).symm ▸ zero_ne_one
  refine ⟨⟨s1, s3, s5, s7⟩, hanti, effSet_eq_of_cover (by norm_num) _ (Or.inl rfl) true _ hanti ?_⟩
  intro v hv
  fin_cases v <;> simp only [mem_insert, mem_singleton] at hv
  · exact ⟨1, by simp, s1⟩
  · exact absurd (Or.inl rfl) hv
  · exact ⟨3, by simp, s3⟩
  · exact absurd (Or.inr (Or.inl rfl)) hv
  · exact ⟨5, by simp, s5⟩
  · exact absurd (Or.inr (Or.inr (Or.inl rfl))) hv
  · exact ⟨7, by simp, s7⟩
  · exact absurd (Or.inr (Or.inr (Or.inr rfl))) hv

end Vignette

end

end NNS.GapDominance
