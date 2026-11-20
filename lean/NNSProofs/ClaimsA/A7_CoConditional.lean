-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.GapCoMoments
public import NNSProofs.MeasurePM

namespace NNS.ClaimsA

open MeasureTheory ProbabilityTheory Set Filter Topology Finset
open scoped ENNReal

@[expose] public noncomputable section

namespace A7_CoConditional

inductive Tail where
  | lower
  | upper

open Tail

def Tail.region : Tail → ℝ → Set ℝ
  | lower => Iic
  | upper => Ioi

def Tail.xEvent (s : Tail) (a : ℝ) := s.region a ×ˢ (univ : Set ℝ)
def Tail.yEvent (s : Tail) (b : ℝ) := (univ : Set ℝ) ×ˢ s.region b

def Tail.xMoment : Tail → ℝ → Measure (ℝ × ℝ) → ℝ
  | lower => mLPMx 0
  | upper => mUPMx 0

def Tail.yMoment : Tail → ℝ → Measure (ℝ × ℝ) → ℝ
  | lower => mLPMy 0
  | upper => mUPMy 0

def cell : Tail → Tail → ℝ → ℝ → Measure (ℝ × ℝ) → ℝ
  | lower, lower => mCoLPM 0 0
  | lower, upper => mDUPM 0 0
  | upper, lower => mDLPM 0 0
  | upper, upper => mCoUPM 0 0

theorem cell_eq (s t : Tail) (a b : ℝ) (μ : Measure (ℝ × ℝ)) :
    cell s t a b μ = μ.real (s.region a ×ˢ t.region b) := by
  cases s <;> cases t
  exacts [mCoLPM_zero .., mDUPM_zero .., mDLPM_zero .., mCoUPM_zero ..]

theorem unique_region (a b : ℝ) (p : ℝ × ℝ) :
    ∃! st : Tail × Tail, p ∈ st.1.region a ×ˢ st.2.region b := by
  by_cases hx : p.1 ≤ a <;> by_cases hy : p.2 ≤ b
  all_goals
    first | refine ⟨(lower, lower), And.intro hx hy, ?_⟩
          | refine ⟨(lower, upper), And.intro hx (lt_of_not_ge hy), ?_⟩
          | refine ⟨(upper, lower), And.intro (lt_of_not_ge hx) hy, ?_⟩
          | refine ⟨(upper, upper), And.intro (lt_of_not_ge hx) (lt_of_not_ge hy), ?_⟩
    rintro ⟨s, t⟩ h
    cases s <;> cases t <;> simp_all [Tail.region, mem_prod, Prod.le_def] <;> grind

end A7_CoConditional

open A7_CoConditional A7_CoConditional.Tail

variable {Ω : Type*} [MeasurableSpace Ω]

theorem claim_C224 (P : Measure Ω) [IsProbabilityMeasure P] {A B : Set Ω}
    (hB : MeasurableSet B) (_h : 0 < P.real B) :
    (P[|B]).real A = P.real (A ∩ B) / P.real B ∧
      P.real (A ∩ B) = (P[|B]).real A * P.real B := by
  rw [inter_comm A B]
  exact ⟨cond_real_eq hB A, (cond_real_mul hB A).symm⟩

theorem claim_C226 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ]
    (a b : ℝ) :
    mLPMx 0 a μ = μ.real (xLe a) ∧ mLPMy 0 b μ = μ.real (yLe b) ∧
    mUPMx 0 a μ = μ.real (xGt a) ∧ mUPMy 0 b μ = μ.real (yGt b) :=
  ⟨mLPMx_zero .., mLPMy_zero .., mUPMx_zero .., mUPMy_zero ..⟩

theorem claim_C227 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ]
    (a b : ℝ) :
    mCoLPM 0 0 a b μ = μ.real (Iic a ×ˢ Iic b) ∧
    mCoUPM 0 0 a b μ = μ.real (Ioi a ×ˢ Ioi b) ∧
    mDUPM 0 0 a b μ = μ.real (Iic a ×ˢ Ioi b) ∧
    mDLPM 0 0 a b μ = μ.real (Ioi a ×ˢ Iic b) :=
  ⟨mCoLPM_zero .., mCoUPM_zero .., mDUPM_zero .., mDLPM_zero ..⟩

theorem claim_C228 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] (a b : ℝ) :
    (∀ p : ℝ × ℝ, ∃! st : Tail × Tail, p ∈ st.1.region a ×ˢ st.2.region b) ∧
    mCoLPM 0 0 a b μ + mCoUPM 0 0 a b μ + mDUPM 0 0 a b μ + mDLPM 0 0 a b μ = 1 :=
  ⟨unique_region a b, quadrants_sum_eq_one_measure ..⟩

theorem claim_C229 {n : ℕ} (hn : 0 < n) (x y : Fin n → ℝ) (a b : ℝ) :
    1 = CoUPM 0 0 a b x y + DUPM 0 0 a b x y + DLPM 0 0 a b x y + CoLPM 0 0 a b x y ∧
    CoUPM 0 0 a b x y = (#{i | a < x i ∧ b < y i} : ℝ) / n ∧
    DUPM 0 0 a b x y = (#{i | x i ≤ a ∧ b < y i} : ℝ) / n ∧
    DLPM 0 0 a b x y = (#{i | a < x i ∧ y i ≤ b} : ℝ) / n ∧
    CoLPM 0 0 a b x y = (#{i | x i ≤ a ∧ y i ≤ b} : ℝ) / n :=
  ⟨by linarith [quadrants_sum_eq_one hn a b x y], GapCoMoments.CoUPM_zero_eq_card ..,
    GapCoMoments.DUPM_zero_eq_card .., GapCoMoments.DLPM_zero_eq_card ..,
    GapCoMoments.CoLPM_zero_eq_card ..⟩

theorem claim_C230 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ]
    (a b : ℝ) (s t : Tail) :
    (0 < s.xMoment a μ → (μ[|s.xEvent a]).real (t.yEvent b) = cell s t a b μ / s.xMoment a μ) ∧
    (0 < t.yMoment b μ → (μ[|t.yEvent b]).real (s.xEvent a) = cell s t a b μ / t.yMoment b μ) := by
  cases s <;> cases t
  exacts [⟨fun _ => cond_yLe_xLe .., fun _ => cond_xLe_yLe ..⟩,
    ⟨fun _ => cond_yGt_xLe .., fun _ => cond_xLe_yGt ..⟩,
    ⟨fun _ => cond_yLe_xGt .., fun _ => cond_xGt_yLe ..⟩,
    ⟨fun _ => cond_yGt_xGt .., fun _ => cond_xGt_yGt ..⟩]

theorem claim_C232 {n : ℕ} (x y : Fin n → ℝ) (a b : ℝ) :
    (pairEmp x y).real (xGt a) = UPM 0 a x ∧ (pairEmp x y).real (yGt b) = UPM 0 b y ∧
    (0 < UPM 0 a x → ((pairEmp x y)[|xGt a]).real (yGt b) = CoUPM 0 0 a b x y / UPM 0 a x) ∧
    (0 < UPM 0 b y → ((pairEmp x y)[|yGt b]).real (xGt a) = CoUPM 0 0 a b x y / UPM 0 b y) := by
  refine ⟨?_, ?_, fun _ => ?_, fun _ => ?_⟩
  · rw [← mUPMx_zero, mUPMx_pairEmp]
  · rw [← mUPMy_zero, mUPMy_pairEmp]
  · rw [cond_yGt_xGt, mCoUPM_pairEmp, mUPMx_pairEmp]
  · rw [cond_xGt_yGt, mCoUPM_pairEmp, mUPMy_pairEmp]

theorem claim_C233 (P : Measure Ω) [IsProbabilityMeasure P] {A B : Set Ω}
    (hA : MeasurableSet A) (hB : MeasurableSet B) (ha : 0 < P.real A) (hb : 0 < P.real B) :
    P.real (A ∩ B) = (P[|B]).real A * P.real B ∧
    P.real (A ∩ B) = (P[|A]).real B * P.real A ∧
    (P[|B]).real A = (P[|A]).real B * P.real A / P.real B := by
  have ab := claim_C224 (A := A) P hB hb
  have ba := (claim_C224 (A := B) P hA ha).2
  rw [inter_comm B A] at ba
  exact ⟨ab.2, ba, by rw [ab.1, ← ba]⟩

theorem claim_C236 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ]
    (a b : ℝ) (s t : Tail)
    (hx : 0 < s.xMoment a μ) (hy : 0 < t.yMoment b μ) :
    (μ[|s.xEvent a]).real (t.yEvent b) = cell s t a b μ / s.xMoment a μ ∧
    (μ[|t.yEvent b]).real (s.xEvent a) = cell s t a b μ / t.yMoment b μ ∧
    cell s t a b μ = (μ[|s.xEvent a]).real (t.yEvent b) * s.xMoment a μ ∧
    (μ[|t.yEvent b]).real (s.xEvent a) =
      (μ[|s.xEvent a]).real (t.yEvent b) * s.xMoment a μ / t.yMoment b μ := by
  obtain ⟨h, k⟩ := claim_C230 μ a b s t
  rw [h hx, k hy, div_mul_cancel₀ _ hx.ne']
  exact ⟨rfl, rfl, rfl, rfl⟩

theorem claim_C234 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ]
    (a b : ℝ) (s t : Tail)
    (hx : 0 < s.xMoment a μ) (hy : 0 < t.yMoment b μ) :
    cell s t a b μ = cell t s b a (μ.map Prod.swap) ∧
    (μ[|t.yEvent b]).real (s.xEvent a) =
      (μ[|s.xEvent a]).real (t.yEvent b) * s.xMoment a μ / t.yMoment b μ := by
  have symmCell : cell s t a b μ = cell t s b a (μ.map Prod.swap) := by
    rw [cell_eq, cell_eq, measureReal_def, measureReal_def, Measure.map_apply measurable_swap]
    · rw [preimage_swap_prod]
    · cases s <;> cases t <;> apply MeasurableSet.prod <;>
        first | exact measurableSet_Iic | exact measurableSet_Ioi
  refine ⟨symmCell, ?_⟩
  obtain ⟨forward, reverse⟩ := claim_C230 μ a b s t
  rw [forward hx, reverse hy, symmCell, div_mul_cancel₀ _ hx.ne']

theorem claim_C235 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ]
    (a b : ℝ) (hx : 0 < mLPMx 0 a μ) (hy : 0 < mLPMy 0 b μ) :
    mCoLPM 0 0 a b μ = (μ[|xLe a]).real (yLe b) * mLPMx 0 a μ ∧
    (μ[|yLe b]).real (xLe a) = (μ[|xLe a]).real (yLe b) * mLPMx 0 a μ / mLPMy 0 b μ ∧
    (μ[|yLe b]).real (xLe a) = (μ[|xLe a]).real (yLe b) * μ.real (xLe a) / μ.real (yLe b) := by
  obtain ⟨_, _, h, k⟩ := claim_C236 μ a b lower lower hx hy
  exact ⟨h, k, by rw [← mLPMx_zero a μ, ← mLPMy_zero b μ]; exact k⟩

theorem claim_C237 (P : Measure Ω) [IsProbabilityMeasure P] {A B : Set Ω}
    (hA : MeasurableSet A) (hB : MeasurableSet B) (ha : 0 < P.real A) :
    (P[|A]).real (Aᶜ ∩ B) = 0 ∧ (P[|A]).real (Aᶜ ∩ Bᶜ) = 0 ∧
    (P[|A]).real (A ∩ B) = P.real (A ∩ B) / P.real A ∧
    (P[|A]).real (A ∩ Bᶜ) = P.real (A ∩ Bᶜ) / P.real A ∧
    (P[|A]).real (A ∩ B) + (P[|A]).real (A ∩ Bᶜ) = 1 := by
  have split : P.real (A ∩ B) + P.real (A ∩ Bᶜ) = P.real A := by
    rw [← measureReal_union (Set.disjoint_left.mpr fun _ h k => k.2 h.2) (hA.inter hB.compl)]
    congr 1
    ext
    simp only [mem_union, mem_inter_iff, mem_compl_iff]
    tauto
  simp only [cond_real_eq hA, ← inter_assoc, inter_compl_self, empty_inter, measureReal_empty,
    zero_div, inter_self, ← add_div, split, div_self ha.ne', and_self]

theorem claim_C238 {n : ℕ} (x y : Fin n → ℝ) (hx : LPM 0 0 x = 4 / 10)
    (hy : LPM 0 0 y = 5 / 10) (hxy : CoLPM 0 0 0 0 x y = 3 / 10) :
    ((pairEmp x y)[|xLe 0]).real (yLe 0) = 3 / 4 ∧
    ((pairEmp x y)[|yLe 0]).real (xLe 0) = 3 / 5 ∧
    ((pairEmp x y)[|xLe 0]).real (yLe 0) * LPM 0 0 x / LPM 0 0 y =
      ((pairEmp x y)[|yLe 0]).real (xLe 0) := by
  norm_num [cond_yLe_xLe, cond_xLe_yLe, mCoLPM_pairEmp, mLPMx_pairEmp, mLPMy_pairEmp, hx, hy, hxy]

namespace A7_CoConditional

def magnitude (a b : ℝ) (μ : Measure (ℝ × ℝ)) :=
  mCoLPM 1 1 a b μ + mCoUPM 1 1 a b μ + mDLPM 1 1 a b μ + mDUPM 1 1 a b μ

end A7_CoConditional

theorem claim_C241 {μ : Measure (ℝ × ℝ)} {a b : ℝ}
    (h : Integrable (fun p : ℝ × ℝ => (p.1 - a) * (p.2 - b)) μ) (hp : 0 < magnitude a b μ) :
    (let S := magnitude a b μ
     0 ≤ mCoLPM 1 1 a b μ / S ∧ 0 ≤ mCoUPM 1 1 a b μ / S ∧
     0 ≤ mDLPM 1 1 a b μ / S ∧ 0 ≤ mDUPM 1 1 a b μ / S ∧
     mCoLPM 1 1 a b μ / S + mCoUPM 1 1 a b μ / S +
       mDLPM 1 1 a b μ / S + mDUPM 1 1 a b μ / S = 1) ∧
    (let ν := pairEmp (![-1, 2] : Fin 2 → ℝ) ![-1, 2]
     magnitude 0 0 ν > 0 ∧ mCoLPM 1 1 0 0 ν / magnitude 0 0 ν ≠ mCoLPM 0 0 0 0 ν) := by
  obtain ⟨_, hs, hll, huu, hul, hlu⟩ := hinge_weights h rfl hp.ne'
  refine ⟨⟨div_nonneg hll hp.le, div_nonneg huu hp.le, div_nonneg hul hp.le,
    div_nonneg hlu hp.le, hs⟩, ?_⟩
  norm_num [magnitude, mCoLPM_pairEmp, mCoUPM_pairEmp, mDLPM_pairEmp, mDUPM_pairEmp,
    CoLPM, CoUPM, DLPM, DUPM, Fin.sum_univ_two, lowerDev, upperDev]

private theorem pair_axis_iff {n : ℕ} (x y : Fin n → ℝ) (a b : ℝ) :
    (∀ᵐ p ∂pairEmp x y, p.1 = a ∨ p.2 = b) ↔ ∀ i, x i = a ∨ y i = b := by
  simp [pairEmp]

theorem claim_C242 {μ : Measure (ℝ × ℝ)} {a b : ℝ}
    (h : Integrable (fun p : ℝ × ℝ => (p.1 - a) * (p.2 - b)) μ) :
    (magnitude a b μ = 0 ↔ ∀ᵐ p ∂μ, p.1 = a ∨ p.2 = b) ∧
    (∀ {n : ℕ}, 0 < n → ∀ (x y : Fin n → ℝ),
      magnitude a b (pairEmp x y) = 0 ↔ ∀ i, x i = a ∨ y i = b) ∧
    (let ν := pairEmp (![0, 1] : Fin 2 → ℝ) ![1, 0]
     magnitude 0 0 ν = 0 ∧ ¬ IndepFun Prod.fst Prod.snd ν) := by
  have sample {n : ℕ} (hn : 0 < n) (x y : Fin n → ℝ) (u v : ℝ) :
      magnitude u v (pairEmp x y) = 0 ↔ ∀ i, x i = u ∨ y i = v := by
    apply Iff.trans (hinge_total_eq_zero_iff ?_) (pair_axis_iff x y u v)
    apply Integrable.smul_measure _ (by simp [hn.ne'])
    exact integrable_finsetSum_measure.mpr fun i _ => integrable_dirac enorm_lt_top
  refine ⟨hinge_total_eq_zero_iff h, fun hn x y => sample hn x y a b,
    (sample (by decide) _ _ 0 0).2 ?_, ?_⟩
  · intro i
    fin_cases i <;> norm_num
  · intro hind
    let ν := pairEmp (![0, 1] : Fin 2 → ℝ) ![1, 0]
    have : IsProbabilityMeasure ν := isProbabilityMeasure_pairEmp (by decide) _ _
    have he := mCoLPM_map_of_indepFun (P := ν) measurable_fst measurable_snd hind 0 0 0 0
    change mCoLPM 0 0 0 0 (ν.map id) = mLPMx 0 0 ν * mLPMy 0 0 ν at he
    norm_num [Measure.map_id, ν, mCoLPM_pairEmp, mLPMx_pairEmp, mLPMy_pairEmp,
      CoLPM, LPM, lowerDev, Fin.sum_univ_two] at he

private theorem rectangle_mass (μ : Measure (ℝ × ℝ)) [IsFiniteMeasure μ]
    {a b c d : ℝ} (hab : a ≤ b) (hcd : c ≤ d) :
    μ.real (Ioc a b ×ˢ Ioc c d) = mCoLPM 0 0 b d μ - mCoLPM 0 0 a d μ -
      mCoLPM 0 0 b c μ + mCoLPM 0 0 a c μ := by
  have xsplit (T : Set ℝ) (hT : MeasurableSet T) :
      μ.real (Ioc a b ×ˢ T) = μ.real (Iic b ×ˢ T) - μ.real (Iic a ×ˢ T) := by
    rw [← measureReal_sdiff (μ := μ)
      (Set.prod_mono (Set.Iic_subset_Iic.mpr hab) (Set.Subset.rfl (s := T)))
      (measurableSet_Iic.prod hT)]
    congr 1
    ext p
    simp only [mem_prod, Set.mem_Ioc, Set.mem_sdiff, Set.mem_Iic]
    grind
  have ysplit (T : Set ℝ) (hT : MeasurableSet T) :
      μ.real (T ×ˢ Ioc c d) = μ.real (T ×ˢ Iic d) - μ.real (T ×ˢ Iic c) := by
    rw [← measureReal_sdiff (μ := μ)
      (Set.prod_mono (Set.Subset.rfl (s := T)) (Set.Iic_subset_Iic.mpr hcd))
      (hT.prod measurableSet_Iic)]
    congr 1
    ext p
    simp only [mem_prod, Set.mem_Ioc, Set.mem_sdiff, Set.mem_Iic]
    grind
  rw [xsplit _ measurableSet_Ioc, ysplit _ measurableSet_Iic, ysplit _ measurableSet_Iic]
  simp only [mCoLPM_zero]
  ring

theorem claim_C244 (μ ν : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    (∀ a b, mCoLPM 0 0 a b μ = μ.real (Iic a ×ˢ Iic b)) ∧
    (∀ a b c d : ℝ, a ≤ b → c ≤ d → μ.real (Ioc a b ×ˢ Ioc c d) =
      mCoLPM 0 0 b d μ - mCoLPM 0 0 a d μ - mCoLPM 0 0 b c μ + mCoLPM 0 0 a c μ) ∧
    ((∀ a b, mCoLPM 0 0 a b μ = mCoLPM 0 0 a b ν) → μ = ν) ∧
    (HingeFinite μ → HingeFinite ν → (∀ a b, mCoLPM 1 1 a b μ = mCoLPM 1 1 a b ν) → μ = ν) :=
  ⟨fun a b => mCoLPM_zero a b μ, fun _ _ _ _ => rectangle_mass μ,
    ext_of_mCoLPM_zero, ext_of_hinge⟩

theorem claim_C245 (μ : Measure (ℝ × ℝ)) :
    (Integrable (fun p : ℝ × ℝ => (1 + |p.1|) * (1 + |p.2|)) μ → HingeFinite μ) ∧
    (HingeFinite μ → ∀ a b,
      HasDerivWithinAt (fun x => mCoLPM 1 1 x b μ) (mCoLPM 0 1 a b μ) (Ici a) a ∧
      HasDerivWithinAt (fun y => mCoLPM 0 1 a y μ) (μ.real (Iic a ×ˢ Iic b)) (Ici b) b ∧
      derivWithin (fun y => derivWithin (fun x => mCoLPM 1 1 x y μ) (Ici a) a)
        (Ici b) b = μ.real (Iic a ×ˢ Iic b)) :=
  ⟨hingeFinite_of_integrable μ, fun h a b => ⟨hasDerivWithinAt_hinge_x h a b,
    by rw [← mCoLPM_zero]; exact hasDerivWithinAt_hinge_y h a b,
    hinge_recovers_joint_cdf h a b⟩⟩

theorem claim_C250 {Ω' : Type*} [MeasurableSpace Ω'] {P : Measure Ω} {Q : Measure Ω'}
    {X : Ω → ℝ} {Y : Ω' → ℝ} (hX : AEMeasurable X P) (hY : AEMeasurable Y Q)
    (h : P.map X = Q.map Y) (d t : ℝ) :
    (∫ ω, rLowerDev d t (X ω) ∂P) = (∫ ω, rLowerDev d t (Y ω) ∂Q) ∧
    (∫ ω, rUpperDev d t (X ω) ∂P) = (∫ ω, rUpperDev d t (Y ω) ∂Q) := by
  constructor
  · rw [← integral_map hX (GapCoMoments.measurable_rLowerDev' d t).aestronglyMeasurable, h,
      integral_map hY (GapCoMoments.measurable_rLowerDev' d t).aestronglyMeasurable]
  · rw [← integral_map hX (GapCoMoments.measurable_rUpperDev' d t).aestronglyMeasurable, h,
      integral_map hY (GapCoMoments.measurable_rUpperDev' d t).aestronglyMeasurable]

theorem claim_C252 {Ω' : Type*} [MeasurableSpace Ω'] {P : Measure Ω} {Q : Measure Ω'}
    {X : Ω → ℝ} {Y : Ω' → ℝ} (hX : AEMeasurable X P) (hY : AEMeasurable Y Q)
    (h : P.map X = Q.map Y) (d t : ℝ) :
    let x : Fin 2 → ℝ := ![-1, 1]
    let y : Fin 2 → ℝ := ![1, -1]
    let μ := pairEmp x x
    let ν := pairEmp x y
    ((∫ ω, rLowerDev d t (X ω) ∂P) = (∫ ω, rLowerDev d t (Y ω) ∂Q) ∧
      (∫ ω, rUpperDev d t (X ω) ∂P) = (∫ ω, rUpperDev d t (Y ω) ∂Q)) ∧
    μ.map Prod.fst = ν.map Prod.fst ∧ μ.map Prod.snd = ν.map Prod.snd ∧
    magnitude 0 0 μ = 1 ∧ magnitude 0 0 ν = 1 ∧
    mCoLPM 1 1 0 0 μ / magnitude 0 0 μ = 1 / 2 ∧
    mCoLPM 1 1 0 0 ν / magnitude 0 0 ν = 0 := by
  dsimp only
  refine ⟨claim_C250 hX hY h d t, ?_, ?_, ?_⟩
  · simp [pairEmp, Fin.sum_univ_two, Measure.map_add, measurable_fst, measurable_fst.aemeasurable]
  · simp [pairEmp, Fin.sum_univ_two, Measure.map_add, measurable_snd,
      measurable_snd.aemeasurable, add_comm]
  · norm_num [magnitude, mCoLPM_pairEmp, mCoUPM_pairEmp, mDLPM_pairEmp, mDUPM_pairEmp,
      CoLPM, CoUPM, DLPM, DUPM, Fin.sum_univ_two, lowerDev, upperDev]

namespace A7_CoConditional

end A7_CoConditional

theorem claim_C255 {μ : Measure (ℝ × ℝ)} (hμ : HingeFinite μ) :
    (∀ a b, derivWithin (fun y => derivWithin (fun x => mCoLPM 1 1 x y μ) (Ici a) a)
      (Ici b) b = μ.real (Iic a ×ˢ Iic b)) ∧
    (∀ a b, (∀ᶠ y in 𝓝 b, DifferentiableAt ℝ (fun x => mCoLPM 1 1 x y μ) a) →
      DifferentiableAt ℝ (fun y => deriv (fun x => mCoLPM 1 1 x y μ) a) b →
      deriv (fun y => deriv (fun x => mCoLPM 1 1 x y μ) a) b = μ.real (Iic a ×ˢ Iic b)) ∧
    ∀ f : ℝ × ℝ → ℝ, μ = volume.withDensity (fun p => ENNReal.ofReal (f p)) →
      Continuous f → (∀ p, 0 ≤ f p) →
      (∀ b x₀ : ℝ, ∃ φ : ℝ → ℝ, Integrable φ ∧
        ∀ᶠ x in 𝓝 x₀, ∀ y, max (b - y) 0 * f (x, y) ≤ φ y) →
      (∃ H₁ H₂ H₃ : ℝ → ℝ → ℝ,
        (∀ a b, HasDerivAt (fun x => mCoLPM 1 1 x b μ) (H₁ a b) a) ∧
        (∀ a b, HasDerivAt (fun x => H₁ x b) (H₂ a b) a) ∧
        (∀ a b, HasDerivAt (fun y => H₂ a y) (H₃ a b) b) ∧
        (∀ a b, HasDerivAt (fun y => H₃ a y) (f (a, b)) b)) ∧
      ∀ a b, deriv (fun v => deriv (fun y => deriv (fun u => deriv (fun x =>
        mCoLPM 1 1 x y μ) u) a) v) b = f (a, b) := by
  refine ⟨hinge_recovers_joint_cdf hμ, hinge_mixed_deriv_eq hμ, ?_⟩
  rintro f rfl hf hf₀ hdom
  have hH (a b : ℝ) : Integrable fun p : ℝ × ℝ => f p * (max (a - p.1) 0 * max (b - p.2) 0) := by
    simpa only [ENNReal.toReal_ofReal (hf₀ _), smul_eq_mul, GapCoMoments.lowerDev_one_eq_hinge]
      using (integrable_withDensity_iff_integrable_smul' hf.measurable.ennreal_ofReal
        (Eventually.of_forall fun _ => ENNReal.ofReal_lt_top)).mp (hμ a b)
  exact ⟨GapCoMoments.hinge_fourth_deriv_eq_density_of_dominated hf hf₀ hH hdom,
    GapCoMoments.deriv_hinge_fourth_eq_density_of_dominated hf hf₀ hH hdom⟩

theorem claim_C256 {μ : Measure (ℝ × ℝ)} [IsProbabilityMeasure μ] (hμ : HingeFinite μ) :
    (∀ a b, derivWithin (fun y => derivWithin (fun x => mCoLPM 1 1 x y μ) (Ici a) a)
      (Ici b) b = μ.real (Iic a ×ˢ Iic b)) ∧
    (∀ (ν : Measure (ℝ × ℝ)) [IsProbabilityMeasure ν], HingeFinite ν →
      (∀ a b, mCoLPM 1 1 a b μ = mCoLPM 1 1 a b ν) → μ = ν) :=
  ⟨hinge_recovers_joint_cdf hμ, fun _ _ => ext_of_hinge hμ⟩

theorem claim_C257 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] (t : ℝ) :
    mDLPM 0 0 t t μ = μ.real (Ioi t ×ˢ Iic t) ∧ mDUPM 0 0 t t μ = μ.real (Iic t ×ˢ Ioi t) :=
  ⟨mDLPM_zero .., mDUPM_zero ..⟩

theorem claim_C258 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] (t : ℝ) :
    (0 < mDLPM 0 0 t t μ - mDUPM 0 0 t t μ ↔ μ.real (yGt t) < μ.real (xGt t)) ∧
    (mDLPM 0 0 t t μ - mDUPM 0 0 t t μ < 0 ↔ μ.real (xGt t) < μ.real (yGt t)) := by
  rw [mDLPM_zero_sub_mDUPM_zero_eq_upper, mUPMx_zero, mUPMy_zero, sub_pos, sub_neg]
  exact ⟨Iff.rfl, Iff.rfl⟩

theorem claim_C259 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] (t : ℝ) :
    μ.real (Ioi t ×ˢ Iic t) - μ.real (Iic t ×ˢ Ioi t) =
      cdf (μ.map Prod.snd) t - cdf (μ.map Prod.fst) t := by
  rw [← mDLPM_zero, ← mDUPM_zero, mDLPM_zero_sub_mDUPM_zero]

theorem claim_C385 (a b : ℝ) (p : ℝ × ℝ) :
    ∃! st : Tail × Tail, p ∈ st.1.region a ×ˢ st.2.region b :=
  unique_region a b p

namespace A7_CoConditional

def signedMoment (a b : ℝ) (μ : Measure (ℝ × ℝ)) :=
  mCoUPM 1 1 a b μ + mCoLPM 1 1 a b μ - mDLPM 1 1 a b μ - mDUPM 1 1 a b μ

end A7_CoConditional

theorem claim_C389 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ]
    (_hx : Integrable Prod.fst μ) (_hy : Integrable Prod.snd μ)
    (h : Integrable (fun p : ℝ × ℝ => (p.1 - ∫ q, q.1 ∂μ) * (p.2 - ∫ q, q.2 ∂μ)) μ) :
    cov[Prod.fst, Prod.snd; μ] = signedMoment (∫ q, q.1 ∂μ) (∫ q, q.2 ∂μ) μ :=
  (covariance_decomposition_measure h).symm

theorem claim_C403 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ]
    (a b : ℝ) :
    (∀ s t, cell s t a b μ = μ.real (s.region a ×ˢ t.region b)) ∧
    (∀ p : ℝ × ℝ, ∃! st : Tail × Tail, p ∈ st.1.region a ×ˢ st.2.region b) :=
  ⟨fun _ _ => cell_eq .., unique_region a b⟩

theorem claim_C408 (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ]
    (hx : Integrable Prod.fst μ) (hy : Integrable Prod.snd μ)
    (h : Integrable (fun p : ℝ × ℝ => (p.1 - ∫ q, q.1 ∂μ) * (p.2 - ∫ q, q.2 ∂μ)) μ) :
    cov[Prod.fst, Prod.snd; μ] = signedMoment (∫ q, q.1 ∂μ) (∫ q, q.2 ∂μ) μ :=
  claim_C389 μ hx hy h

private theorem shifted_integral (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ]
    (hx : Integrable Prod.fst μ) (hy : Integrable Prod.snd μ)
    (hxy : Integrable (fun p : ℝ × ℝ => p.1 * p.2) μ) (a b : ℝ) :
    Integrable (fun p : ℝ × ℝ => (p.1 - a) * (p.2 - b)) μ ∧
    (∫ p, (p.1 - a) * (p.2 - b) ∂μ) =
      (∫ p, p.1 * p.2 ∂μ) - (∫ p, p.1 ∂μ) * b - a * (∫ p, p.2 ∂μ) + a * b := by
  have e : (fun p : ℝ × ℝ => (p.1 - a) * (p.2 - b)) =
      fun p => p.1 * p.2 - p.1 * b - a * p.2 + a * b := funext fun _ => by ring
  have h₁ : Integrable (fun p : ℝ × ℝ => p.1 * p.2 - p.1 * b) μ := hxy.sub (hx.mul_const b)
  have h₂ : Integrable (fun p : ℝ × ℝ => p.1 * p.2 - p.1 * b - a * p.2) μ :=
    h₁.sub (hy.const_mul a)
  rw [e]
  refine ⟨h₂.add (integrable_const _), ?_⟩
  rw [integral_add h₂ (integrable_const _), integral_sub h₁ (hy.const_mul a),
    integral_sub hxy (hx.mul_const b), integral_mul_const, integral_const_mul, integral_const]
  simp

theorem claim_C409 :
    (∀ (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ], Integrable Prod.fst μ →
      Integrable Prod.snd μ → Integrable (fun p : ℝ × ℝ => p.1 * p.2) μ → ∀ a b : ℝ,
      signedMoment a b μ = (∫ p, (p.1 - a) * (p.2 - b) ∂μ) ∧
      signedMoment a b μ = cov[Prod.fst, Prod.snd; μ] +
        ((∫ p, p.1 ∂μ) - a) * ((∫ p, p.2 ∂μ) - b)) ∧
    (let ν := empMeasure (![9, 11] : Fin 2 → ℝ)
     let μ := ν.prod ν
     IndepFun Prod.fst Prod.snd μ ∧ cov[Prod.fst, Prod.snd; μ] = 0 ∧
       signedMoment 0 0 μ = 100) := by
  constructor
  · intro μ _ hx hy hxy a b
    have h := shifted_integral μ hx hy hxy a b
    have hc := (shifted_integral μ hx hy hxy (∫ p, p.1 ∂μ) (∫ p, p.2 ∂μ)).2
    change cov[Prod.fst, Prod.snd; μ] = _ at hc
    have e := coMoment_decomposition_measure h.1
    exact ⟨e, by rw [signedMoment, e, h.2, hc]; ring⟩
  let x : Fin 2 → ℝ := ![9, 11]
  let ν := empMeasure x
  have : IsProbabilityMeasure ν := isProbabilityMeasure_empMeasure (by decide) x
  have hx := memLp_empMeasure (by decide) x 2
  refine ⟨indepFun_prod (X := id) (Y := id) measurable_id measurable_id,
    covariance_fst_snd_prod hx hx, ?_⟩
  norm_num [signedMoment, mCoUPM_prod, mCoLPM_prod, mDLPM_prod, mDUPM_prod, ν,
    integral_empMeasure, x, Fin.sum_univ_two, upperDev, lowerDev]

theorem claim_C417 {P : Measure Ω} [IsProbabilityMeasure P] {X Y : Ω → ℝ}
    (hX : AEMeasurable X P) (hY : AEMeasurable Y P) (h : IndepFun X Y P) (d a b : ℝ) :
    let μ := P.map fun ω => (X ω, Y ω)
    (Integrable (fun ω => rLowerDev d a (X ω)) P →
      Integrable (fun ω => rLowerDev d b (Y ω)) P → GapCoMoments.quadLL d a b μ =
        (∫ ω, rLowerDev d a (X ω) ∂P) * ∫ ω, rLowerDev d b (Y ω) ∂P) ∧
    (Integrable (fun ω => rUpperDev d a (X ω)) P →
      Integrable (fun ω => rUpperDev d b (Y ω)) P → GapCoMoments.quadUU d a b μ =
        (∫ ω, rUpperDev d a (X ω) ∂P) * ∫ ω, rUpperDev d b (Y ω) ∂P) := by
  dsimp only
  rw [(indepFun_iff_map_prod_eq_prod_map_map hX hY).mp h]
  constructor <;> intro _ _
  · rw [GapCoMoments.quadLL, integral_prod_mul,
      integral_map hX (GapCoMoments.measurable_rLowerDev' d a).aestronglyMeasurable,
      integral_map hY (GapCoMoments.measurable_rLowerDev' d b).aestronglyMeasurable]
  · rw [GapCoMoments.quadUU, integral_prod_mul,
      integral_map hX (GapCoMoments.measurable_rUpperDev' d a).aestronglyMeasurable,
      integral_map hY (GapCoMoments.measurable_rUpperDev' d b).aestronglyMeasurable]

theorem claim_C480 {n : ℕ} (x y : Fin n → ℝ) (h : 0 < UPM 0 0 y) :
    ((pairEmp x y)[|yGt 0]).real (xGt 0) = CoUPM 0 0 0 0 x y / UPM 0 0 y :=
  (claim_C232 x y 0 0).2.2.2 h

theorem claim_C561 {n : ℕ} (x y : Fin n → ℝ) (a b : ℝ) :
    CoLPM 0 0 a b x y = (#{i | x i ≤ a ∧ y i ≤ b} : ℝ) / n ∧
    CoLPM 1 1 a b x y = (∑ i, max (a - x i) 0 * max (b - y i) 0) / n :=
  ⟨GapCoMoments.CoLPM_zero_eq_card .., GapCoMoments.CoLPM_one_one_eq_hinge ..⟩

theorem claim_C564 {n : ℕ} (x y : Fin n → ℝ) (a b : ℝ) :
    CoUPM 0 0 a b x y = (#{i | a < x i ∧ b < y i} : ℝ) / n ∧
    CoUPM 1 1 a b x y = (∑ i, max (x i - a) 0 * max (y i - b) 0) / n :=
  ⟨GapCoMoments.CoUPM_zero_eq_card .., GapCoMoments.CoUPM_one_one_eq_hinge ..⟩

end

end NNS.ClaimsA
