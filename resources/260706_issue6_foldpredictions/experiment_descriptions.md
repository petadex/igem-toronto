# Structure Figures Planning Good - Petadex

## Alex Leonardos

**Quick Background**: We have the ESMC6B base and will have a finetuned ESMC6B on the PETadex. We want to use ESMFold2, which takes 6B ESMC embeddings, and draw some conclusions about its use.

### Experiment Flow Diagram: Doesn't show any results, just a diagram of the movement of data and how things are analyzed

- [ ] Experiment Flow Completion

#### *Important Note:* Since the ESMAtlas is so new, the amount of points that have ground truth PDB structures that were introduced after ESMC training was finished is likely very very small. For this reason, I will use prediction quality metrics as a proxy for quality metrics, which is a standard assumption in the field

## *Letters in the following list represent Graphs and Charts:*

- [ ] **A.** *What ESMFold2 Parameters should we use?*

- We can't run the full centroid list on all parameter combinations, so run a *grid-search and see where the predicted quality metrics plateau on a subset*.
- Could also graph the latency depending on the # of loops and sampling steps.
- I'll compare some curves with the paper's defaults: *"We use N = 100 which reduces to 68 sampling steps; we observed no benefit of increased sampling steps for structural quality (Figure S13)"*, with specific focus on the *active site* to see if their defaults would be good or if more steps improves the active site specifically significantly
- This test should be run over a small subset that is indicative of the complexity of the full dataset (not sure how to do this, random, stratified by difficulty?)
  
- [ ] ***FOR NOW: I'll do this analysis on the non-data leakage set from the current petadex version.***

### Group B: This set of figures compares Base vs Finetune (Does finetuning do anything? This should be run over the full centroids)

- [ ] **Ba.** *Does finetune affect TM-Score?*: TM-score(base_pred, finetuned_pred)

- [ ] **Bb**. *Does finetune affect molecular feasibility distributions?*: Overlaid distributions of Molprobity(base_pred) and Molprobity(finetuned_pred)

- [ ] **Bc**. *Does finetune affect molecular feasibility (pointwise)?*: Molprobity score of the predicted structure itself as the coordinate
  - Concretely, x coord: Molprobity(base_pred), y coord: Molprobity(finetune_pred)
  - want points below y=x, as lower molprobity is better
  - Same data as Bb, but it is a pointwise comparison

- [ ] ***FOR NOW: I'll have the TM-Score and Molprobity distributions for the base model.***

### Group C: Compares the predictions to ground truth data when the crystal structures exist.  Note that we treat pLDDT/pTM as proxies for the rest of the analysis, but this ground truth comparison is valid as these graphs compares base to finetune realtive to each other (points relative to y=x on a scatterplot), not absolutely

#### *Ground Truth Note 1: Determining which structures have ground truth data:*

1. Source 1: Cores file
    The cores file has various types of accessions which can be mapped to PDB accessions. This can then by reverse mapped by ORFid.
2. Source 2: ORFs file
    The ORFs file has GenBank accessions which can be mapped to PDB accessions. This can then be reverse mapped by ORFid.

#### *Ground Truth Note 2: Determining which PDBs (if any) minimize data leakage*

1. Method 1: ESMFold2 Training cutoff is 2023, so anything after that is guaranteed not in the training set. (this probably isn't very many)
    Another issue: it was trained on AFDB.

2. Method 2: No near-neighbor advantage: Exclude anything that clusters with a training chain at 40pid (paper does this). Could also do FoldSeek and take the lowest?

    Rationale for this is that you just want to directly minimize leakage.

    Doesn't make as much sense, since <40pid would be much too low for most PETases. Also, we don't need a fully novel fold for this initial benchmark. This is why it wasn't really considered for this first test set.

**This has been completed!**

![alt text](image.png)

- [ ] **Ca**. *Difference in Quality (lDDT)*:
  - Scatterplot with x axis corresponding to base and y corresponding to finetune. Coord value along an axis is a mean lDDT comparison metric of the axis' model's prediction vs the crystal structure.
  - Concretely, x coord: mean-lDDT(base_pred, crystal structure), y coord: mean-lDDT(finetune_pred, crystal structure)
  - Want points above y=x

- [ ] **Cb**. *Difference in Quality*: Same as Ca but with TM-score vs the crystal structure as the comparison function.
  - Concretely, x coord: TM-score(base_pred, crystal structure), y coord: TM-score(finetune_pred, crystal structure)
  - Want points above y=x

- [ ] ***FOR NOW: I can just compute the distributions of these values for just the base model now to ensure that it's not horrible, but they are still most meaningful when relative to the finetune, not absolute.***

### Group D: Metrics for the final distribution from the best model

- [ ] **Da**. *Quality for best model (structural feasibility)*: Distribution of mean-pLDDT or pTM and its corresponding fraction above 70 or 80. (Indicative of usable structures). Use these results to conclude about the *use of the generated structures*.

- [ ] **Db**. *Novelty of new structures (novelty)*: Distribution of FoldSeek best-hit similarity score

  - (searches through PDB or AFDB quickly and finds the closest structure and a similarity score for that structure), and its corresponding fraction above a threshold.

  - Make the other dimension confidence to isolate for confidently novel structures - most interesting.

- [ ] ***FOR NOW: This was supposed to just use the better model between finetune and base, but I can just make the graphs using the base for now.***

**Extra (In the future/downstream plan)**: Compare binding affinity models (bindcraft2, nesso1 and boltzgen)

**Extra (Not in the Current Plan)**. *Accuracy of Confidence*: pLDDT vs LDDT / pTM vs TM - one graph for each of base or finetune (Exp. 2)

- Concretely, x coord for base graph: Mean-pLDDT(base_pred)/100
  y coord for base graph: Mean-LDDT(base_pred, crystal structure)
- Could also do this per-residue instead of over the full protein
