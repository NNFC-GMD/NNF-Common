## The papers behind each step

This is a reading guide to 14 papers that scSLIDE builds on or is compared with: for each, what it did and why it matters here, in plain language. Eleven were asked for by a reader, and three more were added because they are the sources for three methods the reader asked about: diffusion maps, weighted nearest neighbours (WNN) and optimal transport. Each entry was written from the text of the paper, not from memory, and says how much of it was read.

### Testing abundance between two groups

These five papers ask where cells become more or less common between groups of samples, which the scSLIDE paper describes as the prevailing approach.

#### cydar: abundance testing without clusters, in mass cytometry

Lun, Richard and Marioni (2017). *Nature Methods* 14:707-709. [doi:10.1038/nmeth.4295](https://doi.org/10.1038/nmeth.4295). Reference 19 in the scSLIDE paper.

**In one sentence.** Cydar counts each sample's cells inside many small, overlapping neighbourhoods of cells with similar marker levels and tests each one for a change in abundance, with no clustering step.

**What they did.** Mass cytometry measures more than 30 protein markers on each of millions of cells. Earlier analyses clustered the cells first, so results depended on how the clustering was set up. The authors instead centre small balls ("hyperspheres") on cells, count each sample's cells in every ball, and test the counts with a model that allows for variation between replicates. Because balls are denser in some regions than in others, they control a "spatial" false discovery rate: the share of the flagged space, not of the flagged balls, that is flagged wrongly. The method recovered known changes in a reprogramming time course, and in a simulation caught opposite changes in two neighbouring populations that a cluster-based method missed.

**Why scSLIDE cites it.** The scSLIDE paper cites it as an example of the prevailing approach, testing cell abundance between groups.

**For village and differentiation data.** Our reading: this is the template that Milo later adapted to single-cell RNA data. It was built for protein markers, and it relies on replicate samples to judge variability.

**Read:** full text.

#### DA-seq: a per-cell abundance score at several scales

Zhao et al. (2021). *PNAS* 118:e2100293118. [doi:10.1073/pnas.2100293118](https://doi.org/10.1073/pnas.2100293118). Reference 20 in the scSLIDE paper.

**In one sentence.** DA-seq measures how far each cell's most similar cells lean towards one of two biological states, then groups cells with extreme values into regions that need not match any cluster.

**What they did.** Comparing two states, such as healthy and diseased, usually starts by clustering cells, which can miss differences within a cluster or across several. For each cell, DA-seq compares how many nearest neighbours (the most similar cells) come from each state, at several neighbourhood sizes. A statistical model (logistic regression) that predicts the cell's state from these numbers gives each cell one differential abundance (DA) measure. Cells with values beyond what shuffled state labels give by chance, or beyond cut-offs set by the user, are clustered into regions ("DA subpopulations"), each described by marker genes. On melanoma immunotherapy, embryonic mouse skin and COVID-19 data it outlined five regions each: some overlapped published differences, others the authors call new. In a simulation it located planted differences better than cydar.

**Why scSLIDE cites it.** It is another of the scSLIDE paper's examples of abundance testing between groups.

**For village and differentiation data.** Our reading: DA-seq compares exactly two states, such as two time points, and labels each cell only by its state, so one donor line alone could produce a region. When the authors split the melanoma samples into two halves and analysed each, two regions seen in only one half had an excess of cells from single patients. The regions also depend on a cut-off: the figures use 0.8 on a scale from -1 to 1, and in the melanoma data the shuffling test alone flagged roughly 70% of cells. For a time series the authors suggest, as a possible extension, comparing states in pairs. They also warn that uncorrected batch effects (technical differences between samples) can appear as regions.

**Read:** full text of the main paper. Its supplement, which holds the ageing brain analysis and the full simulations, was not read.

#### propeller: testing cell-type proportions with replicates

Phipson et al. (2022). *Bioinformatics* 38:4720-4726. [doi:10.1093/bioinformatics/btac582](https://doi.org/10.1093/bioinformatics/btac582). Reference 21 in the scSLIDE paper.

**In one sentence.** Propeller tests whether the proportion of each cell type differs between groups of samples, judged against the spread between replicate samples.

**What they did.** Cell-type proportions vary widely between samples, and tests that ignore this often report differences that are not there. Propeller transforms each sample's proportions to even out their variability, using the logit or the arcsine square root. It then fits a linear model per cell type, sharing information across cell types, which helps most when replicates are few. In simulations the logit version detected true differences well, with few false positives, while in real data the arcsine version was less sensitive to an outlier sample. On heart development, ageing blood and COVID-19 data it found changes the original studies had not reported.

**Why scSLIDE cites it.** The scSLIDE paper cites it as an example of abundance testing between groups, and for the advice it follows when log-transforming cell-type proportions for its baseline comparison. Propeller's own transformations are the two named above.

**For village and differentiation data.** Our reading: this is the direct test when cell types are labelled and each donor line is a sample. The introduction notes that genetically different samples can be pooled and told apart by genotype, as in a village, but the paper analyses no such data.

**Read:** full text.

#### Milo: abundance testing on overlapping neighbourhoods

Dann et al. (2022). *Nature Biotechnology* 40:245-253. [doi:10.1038/s41587-021-01033-z](https://doi.org/10.1038/s41587-021-01033-z). Reference 22 in the scSLIDE paper.

**In one sentence.** Milo tests for changes in cell abundance inside many small, overlapping neighbourhoods of similar cells, so no clusters are needed.

**What they did.** Clusters fit continuous processes such as differentiation poorly, and earlier cluster-free methods either ignored variation between replicates or were built for two groups. Milo links each cell to its most similar cells, and the cells linked to a chosen cell form a neighbourhood. It counts each sample's cells in every neighbourhood and tests the counts with a model that allows for variation between replicates and for factors such as age or batch. In benchmarks it generally did better than cydar, DA-seq and cluster-based testing. In the ageing mouse thymus it found a declining precursor population that clustering had missed.

**Why scSLIDE cites it.** The scSLIDE paper cites it as an example of abundance testing between groups and as a method that uses density information. In its benchmark, Milo's table of cell counts per neighbourhood and sample was adapted to describe each sample.

**For village and differentiation data.** Our reading: Milo suits a differentiation continuum because it needs no cluster boundaries. It shows where abundance changes, not how samples relate to each other. Biological replicates are required, so each condition needs several donor lines.

**Read:** full text.

#### lochNESS: where mutant cells gather or go missing

Huang et al. (2023). *Nature* 623:772-781. [doi:10.1038/s41586-023-06548-w](https://doi.org/10.1038/s41586-023-06548-w). Reference 35 in the scSLIDE paper.

**In one sentence.** The authors profiled whole mouse embryos from 22 mutant strains, cell by cell, and devised a score, lochNESS, that shows where mutant cells are over- or under-represented.

**What they did.** Conventional checks can miss subtle defects in mutant mice. The team profiled more than 1.6 million nuclei from 101 embryos of 22 mutant and 4 wild-type strains, mostly four embryos per strain. They first tested whether cell-type proportions differed by genotype, but with so few replicates that test had limited power. LochNESS is finer: for each cell it finds the most similar cells from other embryos and compares the share of mutant cells among them with the share expected overall. It exposed changes the proportions had missed, such as two related populations shifting in opposite directions in one mutant.

**Why scSLIDE cites it.** The scSLIDE paper names it as a method that uses density information, and lochNESS is one of its benchmark comparators.

**For village and differentiation data.** Our reading: a village poses the same question with donor lines in place of mutants. The authors suggest such a score matters most when few replicates are possible. LochNESS scores cells, not samples, and the study covers one time point.

**Read:** full text.

### Describing a whole sample as a distribution

These three papers describe each sample by how its cells are spread across cell states, so that whole samples can be compared with one another.

#### GloScope: each sample as a distribution of its cells

Wang, Torous, Gong and Purdom (2024). *Genome Biology* 25:259. [doi:10.1186/s13059-024-03398-1](https://doi.org/10.1186/s13059-024-03398-1). Reference 36 in the scSLIDE paper.

**In one sentence.** GloScope treats each sample as a distribution of cells and measures how far apart the distributions of two samples are.

**What they did.** Samples are hard to compare directly: each has a different number of cells, with no pairing between them. GloScope places all cells in one low-dimensional space (a few summary axes of gene expression), estimates how each sample's cells are spread there, and computes a divergence (a number that grows as two distributions differ) for every pair of samples, without cell-type labels. The divergence it uses is the symmetrised Kullback-Leibler divergence. The Wasserstein distance, the optimal transport measure, is mentioned only as an alternative and as what PILOT uses. The authors use the method mainly to plot samples and spot batch effects and outliers, in studies of 12 to over 300 samples.

**Why scSLIDE cites it.** The scSLIDE paper names it as a density-based method and uses it as a benchmark comparator. Its main text lists GloScope under optimal transport, but its Methods say Kullback-Leibler divergence, which is what the GloScope paper uses.

**For village and differentiation data.** Our reading: a quick first look at whether donor lines or time points separate and whether any sample is an outlier. It says that samples differ, not which cells or genes are responsible.

**Read:** full text.

#### PILOT: distances between samples by optimal transport

Joodaki et al. (2024). *Molecular Systems Biology* 20:57-74. [doi:10.1038/s44320-023-00003-8](https://doi.org/10.1038/s44320-023-00003-8). Reference 81 in the scSLIDE paper.

**In one sentence.** PILOT scores how different two samples are by the least effort needed to turn one's cell-type make-up into the other's, then groups or orders the samples.

**What they did.** Standard analyses need samples labelled as case or control in advance, so they cannot reveal new groups. PILOT describes each sample by the share of its cells in each cluster (a cell type or state). Optimal transport finds the cheapest way to shift one sample's shares into another's, where moving cells between similar clusters costs less, and the total cost is the Wasserstein distance. From these distances PILOT clusters samples, or orders them with diffusion maps into a "disease progression score", then looks for cell types and genes that change along it. Across 12 data sets it recovered known sample groups better than the alternatives tested.

**Why scSLIDE cites it.** PILOT is the optimal transport comparator in the scSLIDE benchmark. Optimal transport is not a step of scSLIDE itself: this comparison is where it enters the paper.

**For village and differentiation data.** Our reading: PILOT needs clustered cells and compares composition, so in a village it would respond to shifts in cell-type abundance between donor lines. Like scSLIDE, it orders samples with diffusion maps. It was tested on disease cohorts.

**Read:** full text.

#### scLKME: describing samples through landmark cells

Yi and Stanley (2023). *bioRxiv* (preprint). [doi:10.1101/2023.11.13.566846](https://doi.org/10.1101/2023.11.13.566846). Reference 39 in the scSLIDE paper.

**In one sentence.** scLKME picks a small set of "landmark" cells and describes each sample by how close, on average, its cells lie to each landmark.

**What they did.** Earlier ways of reducing a sample to one row of numbers either clustered cells into types first or averaged all cells, which loses detail. scLKME selects landmark cells from all samples together (by geometric sketching, kernel herding or random sampling), then averages a kernel, a measure of similarity, between a sample's cells and each landmark. The resulting sample-by-landmark table holds "soft" abundance scores, with no clustering. Against four other methods it was best on two of three data sets and comparable to the best on the third. In a pregnancy cohort, the first axis of a diffusion map of the samples tracked gestational age.

**Why scSLIDE cites it.** It is the source of scSLIDE's "landmark cells" idea, and a comparator in the scSLIDE benchmark.

**For village and differentiation data.** Our reading: this is the closest relative of scSLIDE's description of a sample, and shows that samples can be compared without naming cell types. It is a preprint, tested on blood and heart samples, with 64 to 512 landmarks.

**Read:** full text.

### Ingredients scSLIDE is built from

These five papers supply building blocks that scSLIDE uses or cites: partial least squares, weighted nearest neighbours, two ways of picking representative cells, and diffusion maps.

#### Partial least squares: a regression tutorial

Geladi and Kowalski (1986). *Analytica Chimica Acta* 185:1-17. [doi:10.1016/0003-2670(86)80028-9](https://doi.org/10.1016/0003-2670%2886%2980028-9). Reference 40 in the scSLIDE paper.

**In one sentence.** This tutorial, written for chemists, explains partial least squares (PLS) regression, which predicts outcomes from measurements through summary axes chosen with those outcomes in view.

**What they did.** Finding most earlier accounts of PLS incomplete and often hard to follow, the authors build it up in steps. Multiple linear regression predicts an outcome from all measured variables at once. It has no single answer when variables outnumber the chemical samples, and can fail when variables are collinear (carry overlapping information). Principal component regression replaces the measurements with summary axes and drops the minor ones, without looking at the outcomes, so predictive information can be lost. PLS finds axes for measurements and outcomes together: in each round, each side's axes are updated from the other's scores (a sample's position along an axis). The tutorial also covers rescaling variables whose units differ and choosing how many axes to keep. No data are analysed: a simulated example is left to a companion paper.

**Why scSLIDE cites it.** The scSLIDE paper cites it for partial least squares itself. That paper describes PLS as finding a few summary axes, chosen so that positions along them vary together with the outcome variables as strongly as possible; Step 1 of this page shows this on a small example. The tutorial words it differently: its stated aim is to account for the outcomes as fully as possible while relating them usefully to the measurements. The PLS code on this page is not the tutorial's algorithm: like scSLIDE itself, it follows a later variant, canonical PLS (reference 110 in the scSLIDE paper).

**For village and differentiation data.** Our reading: the warning about principal component regression carries over. Axes drawn from gene expression alone could leave out a subtle difference between donor lines or time points. Giving PLS a label for every cell, as Step 1 of this page does, is a use the tutorial does not teach. Its outcomes are measured variables in chemistry, and it lists classification (sorting samples into named groups) among the topics it leaves out or does not cover in detail. The authors add that the best number of axes is not always clear.

**Read:** full text.

#### WNN: merging two views of the same cells

Hao et al. (2021). *Cell* 184:3573-3587.e29. [doi:10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048). Reference 10 in the scSLIDE paper.

**In one sentence.** Weighted nearest neighbour (WNN) analysis merges two kinds of measurement made on the same cells into one map of which cells are alike, deciding cell by cell how much to trust each.

**What they did.** New technologies measure two things in one cell, for example RNA and surface proteins, and each can see distinctions the other misses. For every cell, WNN finds the most similar cells in each data type and checks how well each set of neighbours predicts the cell's own profile. The better predictor gets the larger weight for that cell. One neighbour graph is then built from the weighted mix and feeds clustering, plotting or trajectory analysis. In cord blood this separated CD4 and CD8 T cells that RNA alone blurred.

**Why scSLIDE cites it.** scSLIDE uses WNN to join two descriptions of the same cells, one oriented to cell type and one label-aware. It fixes both weights at 0.5, so it does not use the cell-by-cell learning of weights that the paper introduces.

**For village and differentiation data.** Our reading: in scSLIDE the two views get equal say in every cell. The paper notes that WNN assumes its inputs do not define conflicting cell states, so check that the two views broadly agree.

**Read:** full text.

#### Geometric sketching: a subset of cells that keeps rare states

Hie et al. (2019). *Cell Systems* 8:483-493.e7. [doi:10.1016/j.cels.2019.05.003](https://doi.org/10.1016/j.cels.2019.05.003). Reference 45 in the scSLIDE paper.

**In one sentence.** Geometric sketching picks a small subset of cells that covers the whole range of cell states evenly, instead of mirroring how common each state is.

**What they did.** Very large data sets strain analysis tools, and the usual shortcut, picking cells at random, tends to lose rare cell types. The key observation is that common cell types form dense clouds in expression space, while rare ones can occupy as much room, only sparsely. The algorithm covers the data with equal-sized boxes and draws cells evenly across the boxes, which thins crowded regions and keeps sparse ones. A sketch of 2% of 665,858 mouse brain cells held on average 281 cells of the rarest type, against 31 with random sampling. One sketch revealed a rare group of inflammatory macrophages, which the authors backed with an experiment.

**Why scSLIDE cites it.** The scSLIDE paper cites it for the idea of choosing a subset of cells that covers both abundant and rare states.

**For village and differentiation data.** Our reading: a sketch is unrepresentative on purpose. That suits the choice of landmarks, because rare or passing states are covered, but abundance must be counted on all cells. The authors note that doublets (two cells captured as one) may be over-represented.

**Read:** full text.

#### Seurat v5: bridge integration and leverage-score sketching

Hao et al. (2024). *Nature Biotechnology* 42:293-304. [doi:10.1038/s41587-023-01767-y](https://doi.org/10.1038/s41587-023-01767-y). Reference 46 in the scSLIDE paper.

**In one sentence.** This paper shows how to align single-cell data sets that measure different things, and how to handle millions of cells through a well-chosen subset.

**What they did.** "Bridge integration" links two kinds of data, such as RNA and chromatin, through a data set where both were measured in the same cells, used as a dictionary. For scale, the authors compute on a "sketch" of each data set (for example 5,000 cells) and carry the result back to all cells. Cells are picked by leverage score, which is high for cells from rare populations. They merged 19 lung data sets, 1.5 million cells, in 55 minutes on one processor core, and kept a rare cell type that a random subset failed to integrate.

**Why scSLIDE cites it.** scSLIDE picks its landmark cells by this leverage-score sketching, and cites the paper for the on-disk data handling that lets it scale. The copy read here does not mention on-disk storage.

**For village and differentiation data.** Our reading: this is why a rare state in a village can still receive landmarks. A sketch is not a fair sample, so abundance should not be estimated from it.

**Read:** full text of the main paper. Its Methods are supplementary and were not retrieved.

#### Diffusion maps: ordering cells along differentiation

Haghverdi, Buettner and Theis (2015). *Bioinformatics* 31:2989-2998. [doi:10.1093/bioinformatics/btv325](https://doi.org/10.1093/bioinformatics/btv325). Reference 7 in the scSLIDE paper.

**In one sentence.** The authors adapt diffusion maps, which reduce many measurements to a few axes that need not be straight, so that single cells are laid out along the branching paths of differentiation.

**What they did.** Differentiating cells change gradually and branch into lineages, which clusters describe poorly. A diffusion map imagines a walker stepping at random between similar cells, and places cells joined by many paths close together. The authors adapt it: the walker never stays put, crowded regions are weighted down, and missing or undetected values count as uncertain. A rule of thumb suggests how close counts as similar. The authors call the result a pseudotemporal ordering of cells. They judge it by eye, and against published lineage trees, next to principal component analysis (PCA) and t-SNE, two common alternatives. On a clean simulation PCA did about as well. On a noisy, unevenly sampled one and on two mouse qPCR data sets (597 and 429 cells), the diffusion map kept the paths clearer than either. On 90 cells from human oocytes and embryos, PCA and t-SNE did almost as well.

**Why scSLIDE cites it.** scSLIDE applies diffusion maps to its sample-level table, so that samples are placed in order. Its Methods also cite the paper for the point that PCA, a straight-line method, may do poorly when variation is curved or branching. scSLIDE runs the destiny software; the paper itself describes an earlier implementation in Matlab.

**For village and differentiation data.** Our reading: a differentiating village fits the premise of gradual change, but the paper lays out cells, while scSLIDE lays out samples, which are far fewer. The authors put the 90-cell result down to the small number of cells, and we would expect the same of few samples. They note that the map does not show which way time runs, so in our view the order needs an anchor, such as the known time points. They also warn that, once "similar" is drawn widely enough, the map links all cells, so a link on the map does not prove a biological one, and that a fast transition can leave no in-between cells in a snapshot.

**Read:** full text of the main paper. Its supplementary figures, where the text places the tests of each adjustment, were not read.

### After the samples are ordered

Once the samples are in order, scSLIDE tests each gene against that order, and this paper turns the per-gene results into one number per cell type.

#### TRADE: one number for a transcriptome-wide effect

Nadig et al. (2025). *Nature Genetics* 57:1228-1237. [doi:10.1038/s41588-025-02169-3](https://doi.org/10.1038/s41588-025-02169-3). Reference 50 in the scSLIDE paper.

**In one sentence.** TRADE estimates how much a perturbation changes gene expression across the whole transcriptome, including the many small effects that never pass a significance threshold.

**What they did.** In pooled CRISPR screens read by single-cell RNA sequencing, each perturbation is seen in few cells, so many true effects go undetected. TRADE takes the usual per-gene results (an estimated change and its uncertainty), separates true effects from sampling noise, and estimates the spread of the true effects: the "transcriptome-wide impact". Two knockdowns seen in 108 and 7 cells looked alike before this correction, but afterwards only one had a clear impact. Across five large screens, significant genes explained only 13 to 36% of the impact. In the guide titration screen of Jost et al., the data set in this notebook, stronger knockdown gave larger impact.

**Why scSLIDE cites it.** scSLIDE runs TRADE on the per-gene results of its trajectory test, which gives one transcriptome-wide impact number per cell type for ranking cell types.

**For village and differentiation data.** Our reading: this ranks cell types by how strongly expression follows the sample order, without leaning on a significance cut-off. It is only as good as the per-gene test that feeds it: the authors warn that ignoring variation between donors inflates effects.

**Read:** full text.

<!-- markdownlint-disable-file MD041 -->
