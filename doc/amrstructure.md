title: Adaptive Mesh Refinement

## Introduction

This document briefly describes the AMR-related features in AGILE. The
different options can be set in the meshlist part (see [Setting parameters](par.html)) of the
**agile.par** file. For a more extensive description, you can read the article
'Parallel, grid-adaptive approaches for relativistic hydro and
magnetohydrodynamics', R. Keppens, Z. Meliani, A.J. van Marle, P. Delmont, A.
Vlasis, &amp; B. van der Holst, 2011, JCP.*
[doi:10.1016/j.jcp.2011.01.020](http://dx.doi.org/10.1016/j.jcp.2011.01.020).

AGILE uses a standard block-based, octree AMR scheme, where we have
blocks of user-controlled dimension (set by the block_nx1, block_nx2, block_nx3)
in a hierarchically nested manner. To simplify the parallelization, we gave up
flexibility to allow different sized refinement ratios between grid levels,
fixing it to 2. Also, we now use the same time step for all levels. A generic
skeleton code, generic enough to hold for any AMR code having similar
restrictions, is shown below.

```text
timeloop : do
   exit time loop by user-defined criteria
   compute timestep constraint dt over all grids
   conditionally save data
   advance all grids for one time step dt
      distinguish source- and flux-addition strategies
         select numerical scheme, geometric sources, ...
         collect fluxes at fine/coarse boundaries
      after every parallel update, on all processors
         fill all ghost cells for all grid blocks
      at the final temporal update, fix for conservation
   regrid
      quantify the error estimator on all grids
      build the new grid structure, ensure proper nesting
      load-balance the new grid structure
   update time counters
end do timeloop
```
<figure class="doc-figure">
<figcaption>
The block-AMR time-advance loop AGILE implements (<code>timeintegration()</code> in <code>src/agile.fpp</code>). All levels share one timestep <code>dt</code> and a fixed 2:1 refinement ratio; regridding and load balancing run once per step, after the update and conservative fix-up.
</figcaption>
</figure>

Some more info follows on the different aspects involved.

## Important (global) parameters

Some important global parameters are in the module
`src/mod_global_parameters.fpp`.
In particular, note that the maximum number of levels is the parameter
**nlevelshi**. The latter is default set to 20. If your want to run with more
levels, you need to change their
value and recompile. The number of levels set in the par file as
_refine_max_level_ must always be smaller or equal to **nlevelshi**.
This **nlevelshi** also returns in those parameters that are defined per
level, such as _limiter_ which needs to be set for all (default) 20 levels.

The parameter **max_blocks** controls how many blocks can be maximally in use **per MPI task**. 
It is set at runtime and can be steered by setting it in the `meshlist` of the parameter file.  
Since AGILE pre-allocates a lot of its GPU datastructures, **max_blocks** is the knob to control AGILE's memory consumption.  If you run out of memory, check how many blocks you expect and set **max_blocks** accordingly. 

## AMR criteria

This in essence describes the module `src/amr/mod_errest.fpp`, or at least its
most essential aspects.

The block-tree nature implies that a decision for refining/coarsening is to be
made on a block-by-block basis. This automated block-based regridding
procedure involves 3 steps:

```text
    (1) consider all blocks at level 1 < l < refine_max_level , with refine_max_level the maximal grid level selected;
    (2) quantify the local error $E_x$ at each gridpoint in a certain grid block;
    (3) if any point has this error exceeding a user-set threshold refine_threshold(l), refine this block (and ensure proper nesting);
    (4) if all points have their error below a user-set fraction of the threshold derefine_ratio(l) used in the previous step, coarsen the block (for l>1).
```    

The local error estimator can be one of three options, selected by
_refine_criterion_, each possibly augmented with user-defined criteria. For _refine_criterion=1_, only refinement based on `usr_refine_grid` is active. Any of
the other 3 estimators use a user-selected subset of the conserved or auxiliary
variables (or even variables that are computed dynamically at the time of
regridding), through the formula
\( E_{\mathbf{x}} = \sum_{iw}\sigma_{iw}\,E^{\mathrm{Rel}}_{\mathbf{x},iw} \).

The indices included are user-identified with the _w_refine_weight_ array, where the sum of
 weights of all variables should equal to 1. The estimated error is a weighted sum
of contributions from all variables with non-zero weight.

For _refine_criterion=3_, we select a Lohner type [R. Lohner, An adaptive finite element scheme for transient problems in CFD, Comp. Meth. App. Mech. Eng. 61, 323
(1987)] prescription as also used in the PARAMESH library or the FLASH3 code.
In our experience, it is computationally efficient as it employs only instantaneous values from
t^n. It in essence discretizes a weighted second derivative in each grid
point.

The Lohner prescription writes in formulae as

\begin{eqnarray}
  E^{\mathrm{Rel}}_{\mathbf{x},iw} &=& \sqrt{\frac{N_{iw}}{\max\!\left(D_{iw},\,\epsilon\right)}} \,,\\[2pt]
  N_{iw} &=& \sum_{i_1}\sum_{i_2}\left[\,\Delta_{i_1}\!\left(\Delta_{i_2} w_{iw}\right)\right]^2 \,,\\[2pt]
  D_{iw} &=& \sum_{i_1}\sum_{i_2}\left[\,S_{i_1}\!\left(\left|\Delta_{i_1} w_{iw}\right|\right)
             + f_{\mathrm{wave},1}\,S_{i_2}\!\left(S_{i_1}\!\left(\left|w_{iw}\right|\right)\right)\right]^2 \,,
\end{eqnarray}

where \( \Delta_i \) and \( S_i \) mean a central difference and a sum along
dimension \( i \), \( \epsilon \) a floor guarding the division, and
\( f_{\mathrm{wave},1} \) the wave-filter parameter. The wave filter is set per
level, and defaults as _amr_wavefilter(1:nlevelshi)=1.0d-2_.


## Data structures

The data structures are defined in `src/mod_physicaldata.fpp` and
`src/amr/mod_forest.fpp`, you can inspect them for learning more details.

We provide details on useful data structures. All of these are suited for any
curvilinear coordinate system, and merely reflect the tree structure of the
block-AMR. We implicitly assume a fixed refinement ratio of two. Schematic
figures for a 2D Cartesian case generalize straightforwardly to higher or
lower dimensionality.

The overall domain is considered 'rectangular', i.e. bounded by the coordinate
pairs `xprobmin1,xprobmax1`, `xprobmin2,xprobmax2` and `xprobmin3,xprobmax3`.
On the lowest grid level l=1, one controls the coarsest resolution as well as a
suitable domain decomposition, by specifying both the total number of level 1
grid cells `domain_nx1, domain_nx2, domain_nx3` along with the individual block
size per dimension `block_nx1, block_nx2, block_nx3`, which exclude the ghost
cells. The total cell number must be an integer multiple of the block size in
each dimension, so e.g. `domain_nx1 = 4*block_nx1`.

A hypothetical 2D domain is shown below, which corresponds to a domain where 4
by 3 blocks on level 1 are exploited in this domain decomposition, and where
local refinement was activated in the top-right corner: the two level-1 blocks
(3,2) and (3,3) plus (4,2) and (4,3) are resolved at level l=2, and one of
those level-2 blocks is refined once more to level l=3. Global, integer grid
indices are introduced per dimension, in a manner where knowledge of these grid
indices, combined with AMR level knowledge, instantly allows one to localize
the grid when needed. Following the figure, the grid on level l=2 indicated by
global grid indices (5,3) is indeed the fifth grid block horizontally, and the
third vertically, when the domain would be resolved fully with level l=2
blocks. The total amount of grid blocks per dimension, per level l, is stored
in `ng1(l), ng2(l), ng3(l)`, and the actual length of a grid block on level l,
per dimension, is `dg1(l), dg2(l), dg3(l)`.

<figure class="doc-figure wide">
<svg viewBox="0 0 720 400" role="img" aria-label="A 4 by 3 grid of level-1 blocks. The four blocks of the top-right corner are each subdivided into a 2 by 2 set of level-2 blocks, and one of those level-2 blocks is subdivided again into 2 by 2 level-3 blocks. Blocks carry their global integer index pair; the level-2 block (5,3) is highlighted.">
<defs>
<marker id="amrA1" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="var(--fig-accent,#B0501F)"/></marker>
</defs>
<text x="35" y="18" font-family="monospace" font-size="12" fill="var(--fig-soft,#566470)">4&#215;3 level-1 blocks, top-right corner refined to l=2 and l=3</text>
<g stroke="currentColor" stroke-width="1.3" fill="none">
<rect x="35" y="30" width="150" height="110"/><rect x="185" y="30" width="150" height="110"/>
<rect x="35" y="140" width="150" height="110"/><rect x="185" y="140" width="150" height="110"/>
<rect x="35" y="250" width="150" height="110"/><rect x="185" y="250" width="150" height="110"/><rect x="335" y="250" width="150" height="110"/><rect x="485" y="250" width="150" height="110"/>
<rect x="335" y="30" width="300" height="220"/>
</g>
<g stroke="var(--fig-diagram,#2B6C7C)" stroke-width="0.9">
<line x1="410" y1="30" x2="410" y2="250"/><line x1="485" y1="30" x2="485" y2="250"/><line x1="560" y1="30" x2="560" y2="250"/>
<line x1="335" y1="85" x2="635" y2="85"/><line x1="335" y1="140" x2="635" y2="140"/><line x1="335" y1="195" x2="635" y2="195"/>
</g>
<g stroke="var(--fig-diagram,#2B6C7C)" stroke-width="0.9">
<line x1="447.5" y1="140" x2="447.5" y2="195"/><line x1="410" y1="167.5" x2="485" y2="167.5"/>
</g>
<rect x="335" y="195" width="75" height="55" fill="var(--fig-accent-bg,#F2E4DA)" stroke="var(--fig-accent,#B0501F)" stroke-width="2.2"/>
<g font-family="monospace" font-size="12" fill="currentColor" text-anchor="middle">
<text x="110" y="90">(1,3)</text><text x="260" y="90">(2,3)</text>
<text x="110" y="200">(1,2)</text><text x="260" y="200">(2,2)</text>
<text x="110" y="310">(1,1)</text><text x="260" y="310">(2,1)</text><text x="410" y="310">(3,1)</text><text x="560" y="310">(4,1)</text>
</g>
<g font-family="monospace" font-size="9.5" fill="currentColor" text-anchor="middle">
<text x="372.5" y="171">(5,4)</text><text x="447.5" y="226">(6,3)</text>
<text x="372.5" y="226" fill="var(--fig-accent,#B0501F)">(5,3)</text>
</g>
<g font-family="monospace" font-size="7.5" fill="currentColor" text-anchor="middle">
<text x="428.75" y="156">(11,8)</text><text x="466.25" y="156">(12,8)</text><text x="428.75" y="183">(11,7)</text><text x="466.25" y="183">(12,7)</text>
</g>
<line x1="372" y1="250" x2="372" y2="376" stroke="var(--fig-accent,#B0501F)" stroke-width="1" stroke-dasharray="3 3"/>
<text x="35" y="384" font-family="monospace" font-size="10" fill="var(--fig-soft,#566470)">(5,3): 5th block across, 3rd up, counted as if the whole domain were at level 2</text>
<text x="35" y="398" font-family="monospace" font-size="9.5" fill="var(--fig-faint,#8B98A1)">ng&#94;D(l): block count per dimension at level l &#183; dg&#94;D(l): a block's length</text>
</svg>
<figcaption>
<b>Domain decomposition and global grid indices.</b> On level 1, <code>domain_nx1/block_nx1</code>, <code>domain_nx2/block_nx2</code> and <code>domain_nx3/block_nx3</code> level-1 blocks tile the rectangular domain. Local refinement replaces a level-1 block by a 2&#215;2&#215;2 set of level-2 blocks, recursively. The global indices <code>ig1(l), ig2(l), ig3(l)</code> address any block as if the whole domain were resolved at level <code>l</code>, so the highlighted level-2 block is <code>(5,3)</code> and the four level-3 blocks inside <code>(6,4)</code> are <code>(11..12, 7..8)</code>. <code>ng1(l), ng2(l), ng3(l)</code> is the level-<code>l</code> block count per dimension; <code>dg1(l), dg2(l), dg3(l)</code> a block's physical length per dimension.
</figcaption>
</figure>

The figure below reflects the tree representation of the same hypothetical
grid hierarchy, where the presence of a grid leaf at a certain grid level is
identified through a boolean variable. As indicated before, the total number
of active grid leafs _nleafs_ may change from timestep to timestep. This tree
info is stored in the structure `tree_root(ig1,ig2,ig3)`, indexed by the
level-l global grid indices, whose `%node` knows about the global grid index
through `%node%ig1`, `%node%ig2`, `%node%ig3`, the level through `%node%level`,
the processor on which it resides through the integer `%node%ipe`, and its
presence or absence through the logical `%node%leaf`.

<figure class="doc-figure wide">
<svg viewBox="0 0 720 430" role="img" aria-label="A 4 by 3 arrangement of cells matching the domain decomposition. Eight cells hold a single leaf node marked T. Four cells hold an internal node marked F with four child leaves; one of those four also has an F child that itself carries four leaves.">
<text x="25" y="18" font-family="monospace" font-size="12" fill="var(--fig-soft,#566470)">the forest: one entry per level-1 block, in domain layout</text>
<g stroke="var(--fig-faint,#8B98A1)" stroke-width="0.8" fill="none">
<rect x="25" y="30" width="165" height="125"/><rect x="190" y="30" width="165" height="125"/><rect x="355" y="30" width="165" height="125"/><rect x="520" y="30" width="165" height="125"/>
<rect x="25" y="155" width="165" height="125"/><rect x="190" y="155" width="165" height="125"/><rect x="355" y="155" width="165" height="125"/><rect x="520" y="155" width="165" height="125"/>
<rect x="25" y="280" width="165" height="125"/><rect x="190" y="280" width="165" height="125"/><rect x="355" y="280" width="165" height="125"/><rect x="520" y="280" width="165" height="125"/>
</g>
<g stroke="currentColor" stroke-width="1" fill="none">
<line x1="437.5" y1="53" x2="392.5" y2="97"/><line x1="437.5" y1="53" x2="422.5" y2="97"/><line x1="437.5" y1="53" x2="452.5" y2="97"/><line x1="437.5" y1="53" x2="482.5" y2="97"/>
<line x1="602.5" y1="53" x2="557.5" y2="97"/><line x1="602.5" y1="53" x2="587.5" y2="97"/><line x1="602.5" y1="53" x2="617.5" y2="97"/><line x1="602.5" y1="53" x2="647.5" y2="97"/>
<line x1="437.5" y1="176" x2="392.5" y2="216"/><line x1="437.5" y1="176" x2="418.5" y2="216"/><line x1="437.5" y1="176" x2="444.5" y2="216"/><line x1="437.5" y1="176" x2="470.5" y2="216"/>
<line x1="470.5" y1="228" x2="452.5" y2="262"/><line x1="470.5" y1="228" x2="466.5" y2="262"/><line x1="470.5" y1="228" x2="480.5" y2="262"/><line x1="470.5" y1="228" x2="494.5" y2="262"/>
<line x1="602.5" y1="176" x2="557.5" y2="216"/><line x1="602.5" y1="176" x2="587.5" y2="216"/><line x1="602.5" y1="176" x2="617.5" y2="216"/><line x1="602.5" y1="176" x2="647.5" y2="216"/>
</g>
<g font-family="monospace" font-size="12" text-anchor="middle">
<circle cx="107" cy="92" r="12" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="107" y="96" fill="currentColor">T</text>
<circle cx="272" cy="92" r="12" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="272" y="96" fill="currentColor">T</text>
<circle cx="107" cy="217" r="12" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="107" y="221" fill="currentColor">T</text>
<circle cx="272" cy="217" r="12" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="272" y="221" fill="currentColor">T</text>
<circle cx="107" cy="342" r="12" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="107" y="346" fill="currentColor">T</text>
<circle cx="272" cy="342" r="12" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="272" y="346" fill="currentColor">T</text>
<circle cx="437" cy="342" r="12" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="437" y="346" fill="currentColor">T</text>
<circle cx="602" cy="342" r="12" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="602" y="346" fill="currentColor">T</text>
</g>
<g font-family="monospace" font-size="11" text-anchor="middle">
<circle cx="437.5" cy="48" r="11" fill="var(--fig-accent-bg,#F2E4DA)" stroke="var(--fig-accent,#B0501F)"/><text x="437.5" y="52" fill="var(--fig-accent,#B0501F)">F</text>
<circle cx="602.5" cy="48" r="11" fill="var(--fig-accent-bg,#F2E4DA)" stroke="var(--fig-accent,#B0501F)"/><text x="602.5" y="52" fill="var(--fig-accent,#B0501F)">F</text>
<circle cx="437.5" cy="171" r="11" fill="var(--fig-accent-bg,#F2E4DA)" stroke="var(--fig-accent,#B0501F)"/><text x="437.5" y="175" fill="var(--fig-accent,#B0501F)">F</text>
<circle cx="602.5" cy="171" r="11" fill="var(--fig-accent-bg,#F2E4DA)" stroke="var(--fig-accent,#B0501F)"/><text x="602.5" y="175" fill="var(--fig-accent,#B0501F)">F</text>
</g>
<g font-family="monospace" font-size="10" text-anchor="middle" fill="currentColor">
<circle cx="392.5" cy="105" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="392.5" y="108.5">T</text>
<circle cx="422.5" cy="105" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="422.5" y="108.5">T</text>
<circle cx="452.5" cy="105" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="452.5" y="108.5">T</text>
<circle cx="482.5" cy="105" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="482.5" y="108.5">T</text>
<circle cx="557.5" cy="105" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="557.5" y="108.5">T</text>
<circle cx="587.5" cy="105" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="587.5" y="108.5">T</text>
<circle cx="617.5" cy="105" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="617.5" y="108.5">T</text>
<circle cx="647.5" cy="105" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="647.5" y="108.5">T</text>
<circle cx="392.5" cy="224" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="392.5" y="227.5">T</text>
<circle cx="418.5" cy="224" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="418.5" y="227.5">T</text>
<circle cx="444.5" cy="224" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="444.5" y="227.5">T</text>
<circle cx="557.5" cy="224" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="557.5" y="227.5">T</text>
<circle cx="587.5" cy="224" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="587.5" y="227.5">T</text>
<circle cx="617.5" cy="224" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="617.5" y="227.5">T</text>
<circle cx="647.5" cy="224" r="8" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="647.5" y="227.5">T</text>
</g>
<circle cx="470.5" cy="222" r="9" fill="var(--fig-accent-bg,#F2E4DA)" stroke="var(--fig-accent,#B0501F)"/><text x="470.5" y="225.5" font-family="monospace" font-size="10" fill="var(--fig-accent,#B0501F)" text-anchor="middle">F</text>
<g font-family="monospace" font-size="9" text-anchor="middle" fill="currentColor">
<circle cx="452.5" cy="270" r="6.5" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="452.5" y="272.5">T</text>
<circle cx="466.5" cy="270" r="6.5" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="466.5" y="272.5">T</text>
<circle cx="480.5" cy="270" r="6.5" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="480.5" y="272.5">T</text>
<circle cx="494.5" cy="270" r="6.5" fill="var(--fig-surface,#fff)" stroke="currentColor"/><text x="494.5" y="272.5">T</text>
</g>
<text x="25" y="424" font-family="monospace" font-size="10" fill="var(--fig-faint,#8B98A1)">T = leaf (%leaf = .true., holds solution data) &#183; F = internal node, refined into 2&#94;ndim children</text>
</svg>
<figcaption>
<b>The forest of octrees.</b> One <code>tree_root</code> entry per level-1 block, laid out in the domain's geometry. Each node stores its global indices <code>%ig1, %ig2, %ig3</code>, its <code>%level</code>, the owning processor <code>%ipe</code>, and <code>%leaf</code>. A leaf ("T") holds a block of solution data; an internal node ("F") has been refined and its data lives in its children. <code>nleafs</code>, the count of "T" nodes, changes at every regrid.
</figcaption>
</figure>

Various extra indices are helpful to traverse the tree structure. Local grid
indices across AMR levels are schematically given below, which are used to
identify the directional neighbours, as well as the children and parent
blocks. These are used to realize and facilitate the possible interprocessor
communication patterns.

<figure class="doc-figure wide">
<svg viewBox="0 0 720 300" role="img" aria-label="Left: a 3 by 3 coarse block neighbourhood with the centre block at offset (0,0) and neighbours at offsets -1, 0, +1 per dimension, connected to a 2 by 2 set of child blocks indexed 1..2 per dimension. Right: a 4 by 4 grid of blocks with the central 2 by 2 group highlighted as a block together with its three siblings under one parent.">
<text x="15" y="18" font-family="monospace" font-size="12" fill="var(--fig-soft,#566470)">neighbour offsets and child indices</text>
<g stroke="currentColor" stroke-width="1.1" fill="none">
<rect x="15" y="30" width="70" height="70"/><rect x="85" y="30" width="70" height="70"/><rect x="155" y="30" width="70" height="70"/>
<rect x="15" y="100" width="70" height="70"/><rect x="155" y="100" width="70" height="70"/>
<rect x="15" y="170" width="70" height="70"/><rect x="85" y="170" width="70" height="70"/><rect x="155" y="170" width="70" height="70"/>
</g>
<rect x="85" y="100" width="70" height="70" fill="var(--fig-diagram-bg,#DDEBEC)" stroke="currentColor" stroke-width="1.1"/>
<g font-family="monospace" font-size="10" fill="currentColor" text-anchor="middle">
<text x="50" y="69">(-1,1)</text><text x="120" y="69">(0,1)</text><text x="190" y="69">(1,1)</text>
<text x="50" y="139">(-1,0)</text><text x="120" y="139">(0,0)</text><text x="190" y="139">(1,0)</text>
<text x="50" y="209">(-1,-1)</text><text x="120" y="209">(0,-1)</text><text x="190" y="209">(1,-1)</text>
</g>
<g stroke="var(--fig-faint,#8B98A1)" stroke-width="0.8" stroke-dasharray="3 3">
<line x1="155" y1="100" x2="270" y2="120"/><line x1="155" y1="170" x2="270" y2="210"/>
</g>
<g stroke="var(--fig-diagram,#2B6C7C)" stroke-width="1.2" fill="none">
<rect x="270" y="120" width="90" height="90"/><line x1="315" y1="120" x2="315" y2="210"/><line x1="270" y1="165" x2="360" y2="165"/>
</g>
<g font-family="monospace" font-size="10" fill="var(--fig-diagram,#2B6C7C)" text-anchor="middle">
<text x="292" y="148">(1,2)</text><text x="338" y="148">(2,2)</text><text x="292" y="192">(1,1)</text><text x="338" y="192">(2,1)</text>
</g>
<text x="120" y="262" font-family="monospace" font-size="10" fill="var(--fig-soft,#566470)" text-anchor="middle">neighbour offset i&#94;D &#8712; {-1,0,1}</text>
<text x="315" y="234" font-family="monospace" font-size="10" fill="var(--fig-soft,#566470)" text-anchor="middle">children: 1..2 per dim</text>
<line x1="400" y1="20" x2="400" y2="280" stroke="var(--fig-faint,#8B98A1)" stroke-width="0.8"/>
<text x="430" y="18" font-family="monospace" font-size="12" fill="var(--fig-soft,#566470)">a block and its three siblings</text>
<g stroke="currentColor" stroke-width="1.1" fill="none">
<rect x="430" y="30" width="62" height="60"/><rect x="492" y="30" width="62" height="60"/><rect x="554" y="30" width="62" height="60"/><rect x="616" y="30" width="62" height="60"/>
<rect x="430" y="90" width="62" height="60"/><rect x="616" y="90" width="62" height="60"/>
<rect x="430" y="150" width="62" height="60"/><rect x="616" y="150" width="62" height="60"/>
<rect x="430" y="210" width="62" height="60"/><rect x="492" y="210" width="62" height="60"/><rect x="554" y="210" width="62" height="60"/><rect x="616" y="210" width="62" height="60"/>
</g>
<rect x="492" y="90" width="124" height="120" fill="var(--fig-accent-bg,#F2E4DA)"/>
<g stroke="var(--fig-accent,#B0501F)" stroke-width="2" fill="none"><rect x="492" y="90" width="124" height="120"/></g>
<line x1="554" y1="90" x2="554" y2="210" stroke="var(--fig-accent,#B0501F)" stroke-width="0.8" stroke-dasharray="3 3"/>
<line x1="492" y1="150" x2="616" y2="150" stroke="var(--fig-accent,#B0501F)" stroke-width="0.8" stroke-dasharray="3 3"/>
<g font-family="monospace" font-size="10" fill="currentColor" text-anchor="middle">
<text x="461" y="64">(0,3)</text><text x="523" y="64">(1,3)</text><text x="585" y="64">(2,3)</text><text x="647" y="64">(3,3)</text>
<text x="461" y="124">(0,2)</text><text x="523" y="124" fill="var(--fig-accent,#B0501F)">(1,2)</text><text x="585" y="124" fill="var(--fig-accent,#B0501F)">(2,2)</text><text x="647" y="124">(3,2)</text>
<text x="461" y="184">(0,1)</text><text x="523" y="184" fill="var(--fig-accent,#B0501F)">(1,1)</text><text x="585" y="184" fill="var(--fig-accent,#B0501F)">(2,1)</text><text x="647" y="184">(3,1)</text>
<text x="461" y="244">(0,0)</text><text x="523" y="244">(1,0)</text><text x="585" y="244">(2,0)</text><text x="647" y="244">(3,0)</text>
</g>
<text x="554" y="288" font-family="monospace" font-size="10" fill="var(--fig-soft,#566470)" text-anchor="middle">four blocks sharing one parent</text>
</svg>
<figcaption>
<b>Indices used to walk the tree.</b> A block's face neighbours are addressed by offsets <code>i1, i2, i3 &#8712; {-1,0,+1}</code> about its own <code>(0,0,0)</code>; its eight children by indices <code>1..2</code> in each of the three dimensions. The central 2&#215;2 group on the right are a block and its three siblings &#8212; the four blocks produced by refining one parent &#8212; which coarsening treats as a unit.
</figcaption>
</figure>

<figure class="doc-figure">
<svg viewBox="0 0 420 400" role="img" aria-label="Vertical layout: a 2 by 2 box of level l+1 child blocks at the top, a single level-l block with a dot in the middle, and a larger level l-1 parent block at the bottom. Arrows fan from the block up to the four children and one arrow points down from the block into the parent.">
<defs>
<marker id="amrC1" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="8" markerHeight="8" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="var(--fig-accent,#B0501F)"/></marker>
</defs>
<g stroke="currentColor" stroke-width="1.3" fill="none">
<rect x="150" y="26" width="120" height="100"/><line x1="210" y1="26" x2="210" y2="126"/><line x1="150" y1="76" x2="270" y2="76"/>
<rect x="176" y="176" width="68" height="66"/>
<rect x="110" y="292" width="200" height="90"/>
</g>
<circle cx="210" cy="209" r="6" fill="currentColor"/>
<g stroke="var(--fig-accent,#B0501F)" stroke-width="1.5" fill="none">
<path d="M204 176 L178 130" marker-end="url(#amrC1)"/><path d="M208 176 L198 130" marker-end="url(#amrC1)"/><path d="M212 176 L222 130" marker-end="url(#amrC1)"/><path d="M216 176 L242 130" marker-end="url(#amrC1)"/>
<path d="M210 242 L210 288" marker-end="url(#amrC1)"/>
</g>
<g font-family="monospace" font-size="12" fill="currentColor" text-anchor="middle">
<text x="210" y="18">children (level l+1)</text>
<text x="300" y="213" text-anchor="start">block igrid</text>
<text x="300" y="228" text-anchor="start">(level l)</text>
<text x="210" y="342">parent (level l-1)</text>
</g>
</svg>
<figcaption>
<b>Every block knows its parent and children.</b> The tree record for a block points to its parent, one level coarser, and to its <code>child(2,2,2)</code> &#8212; eight children, one level finer. Prolongation and the fine-block ghost fill read the parent; restriction (coarsening) and the conservative flux fix-up at fine/coarse boundaries read the children. Any of them may live on another processor.
</figcaption>
</figure>

The directional neighbours of a grid block, and the exchange that fills a
block's ghost layer from them, are shown below for the 3D case.

<figure class="doc-figure wide">
<svg viewBox="0 0 720 300" role="img" aria-label="Left: a 3D block drawn as an isometric cube with six arrows to six small squares labelled plus and minus x, y and z, its six face neighbours. Right: one slice through the block's ghost layer as a 3 by 3 arrangement, with arrows in from the four face neighbours and the four diagonal edge/corner neighbours.">
<defs>
<marker id="amrD1" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="currentColor"/></marker>
<marker id="amrD2" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="var(--fig-accent,#B0501F)"/></marker>
</defs>
<text x="14" y="18" font-family="monospace" font-size="12" fill="var(--fig-soft,#566470)">3D: 6 face neighbours (&#177;x, &#177;y, &#177;z)</text>
<polygon points="150,150 200,150 200,210 150,210" fill="var(--fig-accent-bg,#F2E4DA)" stroke="currentColor" stroke-width="1.3"/>
<polygon points="150,150 200,150 222,128 172,128" fill="var(--fig-diagram-bg,#DDEBEC)" stroke="currentColor" stroke-width="1.3"/>
<polygon points="200,150 222,128 222,188 200,210" fill="var(--fig-diagram-bg,#DDEBEC)" stroke="currentColor" stroke-width="1.3"/>
<g stroke="currentColor" stroke-width="1.1" fill="var(--fig-surface,#fff)">
<rect x="60" y="168" width="24" height="24"/><rect x="286" y="150" width="24" height="24"/>
<rect x="176" y="66" width="24" height="24"/><rect x="176" y="252" width="24" height="24"/>
<rect x="278" y="72" width="24" height="24"/><rect x="66" y="238" width="24" height="24"/>
</g>
<g stroke="currentColor" stroke-width="1.1" fill="none">
<path d="M148 180 L90 180" marker-end="url(#amrD1)"/><path d="M204 168 L282 162" marker-end="url(#amrD1)"/>
<path d="M175 152 L191 92" marker-end="url(#amrD1)"/><path d="M180 210 L188 248" marker-end="url(#amrD1)"/>
<path d="M210 140 L280 98" marker-end="url(#amrD1)"/><path d="M155 200 L92 236" marker-end="url(#amrD1)"/>
</g>
<g font-family="monospace" font-size="10" fill="currentColor" text-anchor="middle">
<text x="72" y="163">-x</text><text x="298" y="145">+x</text><text x="188" y="60">+z</text><text x="188" y="290">-z</text><text x="290" y="66">+y</text><text x="78" y="276">-y</text>
</g>
<line x1="380" y1="14" x2="380" y2="286" stroke="var(--fig-faint,#8B98A1)" stroke-width="0.8"/>
<text x="404" y="18" font-family="monospace" font-size="12" fill="var(--fig-soft,#566470)">one slice of the ghost-cell exchange (getbc)</text>
<g stroke="currentColor" stroke-width="1" fill="none">
<rect x="470" y="40" width="60" height="60"/><rect x="530" y="40" width="60" height="60"/><rect x="590" y="40" width="60" height="60"/>
<rect x="470" y="100" width="60" height="60"/><rect x="590" y="100" width="60" height="60"/>
<rect x="470" y="160" width="60" height="60"/><rect x="530" y="160" width="60" height="60"/><rect x="590" y="160" width="60" height="60"/>
</g>
<rect x="530" y="100" width="60" height="60" fill="var(--fig-accent-bg,#F2E4DA)" stroke="var(--fig-accent,#B0501F)" stroke-width="2.2"/>
<g stroke="var(--fig-accent,#B0501F)" stroke-width="1.1">
<path d="M560 84 L560 104" marker-end="url(#amrD2)"/><path d="M560 176 L560 156" marker-end="url(#amrD2)"/><path d="M514 130 L534 130" marker-end="url(#amrD2)"/><path d="M606 130 L586 130" marker-end="url(#amrD2)"/>
</g>
<g stroke="var(--fig-diagram,#2B6C7C)" stroke-width="1.1">
<path d="M514 84 L531 101" marker-end="url(#amrD2)"/><path d="M606 84 L589 101" marker-end="url(#amrD2)"/><path d="M514 176 L531 159" marker-end="url(#amrD2)"/><path d="M606 176 L589 159" marker-end="url(#amrD2)"/>
</g>
<g font-family="monospace" font-size="10" fill="currentColor" text-anchor="middle">
<text x="500" y="74">edge</text><text x="620" y="74">edge</text><text x="560" y="74">face</text><text x="500" y="134">face</text><text x="560" y="134" fill="var(--fig-accent,#B0501F)">block</text>
</g>
<text x="560" y="245" font-family="monospace" font-size="10" fill="var(--fig-soft,#566470)" text-anchor="middle">corners need the diagonal neighbour</text>
<text x="560" y="260" font-family="monospace" font-size="9.5" fill="var(--fig-faint,#8B98A1)" text-anchor="middle">full 3D layer: 6 faces + 12 edges + 8 corners</text>
</svg>
<figcaption>
<b>Directional neighbours and the ghost exchange (3D).</b> A block touches six neighbours across its faces. <code>getbc</code> (<code>src/mod_ghostcells_update.fpp</code>) fills the block's ghost layer from all of them and from the diagonal (edge and corner) neighbours as well &#8212; the corner ghost regions lie in no face neighbour. A neighbour at the same level is copied straight in; a coarser one is prolonged into the ghost cells; finer ones are restricted.
</figcaption>
</figure>

For parallelization, we adopted a fairly straightforward Z-order or Morton-
order space filling curve (SFC). The Morton curve threads the leaf blocks in a
recursive Z pattern &#8212; descending into each refined block in the same Z order
&#8212; which gives every leaf a single integer Morton number. The blocks are then
handed to processors in contiguous arcs of that ordering. Load-balancing is
done after every timestep, following the adaptive remeshing. When exploiting
N_p GPUs, our strategy for load balancing merely ensures that each GPU has at
least [nleafs/N_p]_int (denoting integer division) grid blocks, while the
remainder increase this number by 1 for the first as many GPUs. For the
hypothetical hierarchy above &#8212; 27 leaf blocks on 4 GPUs &#8212; the first three GPUs
hold 7 blocks each and the fourth holds 6. The grid Morton numbers of all
grids residing on processor _mype_ lie between
_Morton_start(mype),Morton_stop(mype)_. The global grid index, once you know
the grid Morton number _Morton_no_ is found from _sfc_to_igrid(Morton_no)_,
which gives the relation between the SFC and the global grid index _igrid_. The
data for the conservative variables for grid _igrid_ is then actually found
from _pw(igrid)%w_.

<figure class="doc-figure wide">
<svg viewBox="0 0 720 400" role="img" aria-label="The same 4 by 3 hierarchy as the domain-decomposition figure, with a Z-shaped polyline threading the centres of all 27 leaf blocks in Morton order. The blocks are tinted in four colours marking the four contiguous arcs assigned to GPUs 0 to 3, of sizes 7, 7, 7 and 6.">
<text x="35" y="18" font-family="monospace" font-size="12" fill="var(--fig-soft,#566470)">Morton (Z-order) curve through the 27 leaf blocks, split 7 / 7 / 7 / 6</text>
<g stroke="currentColor" stroke-width="1.3" fill="none">
<rect x="35" y="30" width="150" height="110"/><rect x="185" y="30" width="150" height="110"/>
<rect x="35" y="140" width="150" height="110"/><rect x="185" y="140" width="150" height="110"/>
<rect x="35" y="250" width="150" height="110"/><rect x="185" y="250" width="150" height="110"/><rect x="335" y="250" width="150" height="110"/><rect x="485" y="250" width="150" height="110"/>
<rect x="335" y="30" width="300" height="220"/>
</g>
<g fill="var(--fig-alt,#40704E)" fill-opacity="0.22">
<rect x="35" y="250" width="150" height="110"/><rect x="185" y="250" width="150" height="110"/><rect x="35" y="140" width="150" height="110"/><rect x="185" y="140" width="150" height="110"/><rect x="335" y="250" width="150" height="110"/><rect x="485" y="250" width="150" height="110"/><rect x="335" y="195" width="75" height="55"/>
</g>
<g fill="var(--fig-diagram,#2B6C7C)" fill-opacity="0.22">
<rect x="410" y="195" width="75" height="55"/><rect x="335" y="140" width="75" height="55"/><rect x="410" y="167.5" width="37.5" height="27.5"/><rect x="447.5" y="167.5" width="37.5" height="27.5"/><rect x="410" y="140" width="37.5" height="27.5"/><rect x="447.5" y="140" width="37.5" height="27.5"/><rect x="485" y="195" width="75" height="55"/>
</g>
<g fill="var(--fig-accent,#B0501F)" fill-opacity="0.20">
<rect x="560" y="195" width="75" height="55"/><rect x="485" y="140" width="75" height="55"/><rect x="560" y="140" width="75" height="55"/><rect x="35" y="30" width="150" height="110"/><rect x="185" y="30" width="150" height="110"/><rect x="335" y="85" width="75" height="55"/><rect x="410" y="85" width="75" height="55"/>
</g>
<g fill="var(--fig-drop,#9A5350)" fill-opacity="0.20">
<rect x="335" y="30" width="75" height="55"/><rect x="410" y="30" width="75" height="55"/><rect x="485" y="85" width="75" height="55"/><rect x="560" y="85" width="75" height="55"/><rect x="485" y="30" width="75" height="55"/><rect x="560" y="30" width="75" height="55"/>
</g>
<g stroke="var(--fig-diagram,#2B6C7C)" stroke-width="0.8">
<line x1="410" y1="30" x2="410" y2="250"/><line x1="485" y1="30" x2="485" y2="250"/><line x1="560" y1="30" x2="560" y2="250"/>
<line x1="335" y1="85" x2="635" y2="85"/><line x1="335" y1="140" x2="635" y2="140"/><line x1="335" y1="195" x2="635" y2="195"/>
<line x1="447.5" y1="140" x2="447.5" y2="195"/><line x1="410" y1="167.5" x2="485" y2="167.5"/>
</g>
<polyline points="110,305 260,305 110,195 260,195 410,305 560,305 372.5,222.5 447.5,222.5 372.5,167.5 428.75,181.25 466.25,181.25 428.75,153.75 466.25,153.75 522.5,222.5 597.5,222.5 522.5,167.5 597.5,167.5 110,85 260,85 372.5,112.5 447.5,112.5 372.5,57.5 447.5,57.5 522.5,112.5 597.5,112.5 522.5,57.5 597.5,57.5" fill="none" stroke="var(--fig-accent,#B0501F)" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>
<circle cx="110" cy="305" r="4.5" fill="var(--fig-accent,#B0501F)"/>
<text x="120" y="300" font-family="monospace" font-size="9" fill="var(--fig-soft,#566470)">Morton 0</text>
<text x="588" y="52" font-family="monospace" font-size="9" fill="var(--fig-soft,#566470)" text-anchor="end">26</text>
<g font-family="monospace" font-size="10" fill="currentColor">
<rect x="35" y="378" width="11" height="11" fill="var(--fig-alt,#40704E)" fill-opacity="0.22" stroke="var(--fig-alt,#40704E)"/><text x="52" y="387">GPU 0 (7)</text>
<rect x="150" y="378" width="11" height="11" fill="var(--fig-diagram,#2B6C7C)" fill-opacity="0.22" stroke="var(--fig-diagram,#2B6C7C)"/><text x="167" y="387">GPU 1 (7)</text>
<rect x="265" y="378" width="11" height="11" fill="var(--fig-accent,#B0501F)" fill-opacity="0.20" stroke="var(--fig-accent,#B0501F)"/><text x="282" y="387">GPU 2 (7)</text>
<rect x="380" y="378" width="11" height="11" fill="var(--fig-drop,#9A5350)" fill-opacity="0.20" stroke="var(--fig-drop,#9A5350)"/><text x="397" y="387">GPU 3 (6)</text>
</g>
</svg>
<figcaption>
<b>Morton space-filling curve and load balance.</b> The recursive Z-curve visits every leaf of the hierarchy from the previous figures, giving it a 1-D Morton number; the domain then splits into four contiguous arcs of that ordering, one per processor, sized to within one block of <code>nleafs/N_p</code>. Because the curve keeps spatial neighbours close in the ordering, each processor's blocks stay clustered, which keeps the ghost-cell communication mostly on-node. <code>Morton_start(mype)</code> / <code>Morton_stop(mype)</code> bound a processor's arc; <code>sfc_to_igrid(Morton_no)</code> maps back to the global grid index <code>igrid</code>.
</figcaption>
</figure>
