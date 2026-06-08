================================================================================
README - Master's Thesis Supplementary Files
Isabel Johanne Nilsen Wolden, NTNU 2026
Numerical Study of the Mechanical Performance of Latticed Sheet
================================================================================

This archive contains two folders:
  - scripts/             Python scripts for analytical calculations and post-processing
  - abaqus_inp_files/    Abaqus input files (.inp) for all finite element simulations

Required Python packages: numpy, matplotlib, scipy, seaborn, pathlib (standard library)
Abaqus version: 2023


================================================================================
SCRIPTS FOLDER STRUCTURE
================================================================================

analytical-framework/
    unit_cell_geometry.py
        Defines unit cell geometry, strut lengths, angles, and relative density
        as functions of the apex parameter a.

    compliance_matrix.py
        Assembles the compliance matrix for the asymmetric triangular unit cell
        and extracts effective elastic moduli (E11, E22, G12, nu12, nu21).

    effective_strength_analysis.py
        Computes effective strength by applying a force balance across the unit
        cell and finding the critical strut that reaches yield first.

    effective_fracture_toughness.py
        Estimates effective fracture toughness K_Ic using the topology factor
        approach based on strut geometry and material fracture strength.

    parametric_strength_fracture_study.py
        Runs a parametric sweep over the apex parameter a and collects
        stiffness, strength, and fracture toughness as functions of geometry.

    plot_scaling_relationships.py
        Generates log-log scaling plots of effective properties versus relative
        density. Produces the figure scaling_relationships.png.

    compare_stiffness.py
        Compares analytical stiffness predictions against FE-measured initial
        stiffness values for all nine symmetric designs.


post-processing/baseline/
    plot_baseline.py
        Plots force-displacement curves for Baseline A and B with all precut
        variants. Requires Baseline.txt data file in the same folder.

    plot_baselineA.py
        Plots force-displacement curves for Baseline A across loading rates.
        Requires BaselineA.txt.

    plot_baselineB.py
        Plots force-displacement curves for Baseline B across loading rates.
        Requires baselineB.txt.

    plot_baseline_summary.py
        Generates a multi-panel summary figure of peak force, displacement,
        energy absorption, and stiffness for all baseline configurations.


post-processing/g1/
    plot_g1_results.py
        Plots force-displacement curves for G1 variants A, B, and C at all
        three loading rates. Also generates per-variant separate figures and
        an overlay plot. Requires data files in the G1 folder.


post-processing/g2/
    plot_G2.py
        Plots force-displacement curves for G2-A and G2-B in intact and precut
        states. Requires G2 data files.


post-processing/g3/
    plot_G3.py
        Plots force-displacement curves for G3 in intact and precut states.
        Requires g3-f-d.txt.


post-processing/symmetric/
    plot_f_d_all.py
        Plots force-displacement curves for all nine symmetric designs. Produces
        separate figures for intact and precut specimens. Requires f-d-all.txt
        and f-d-precut.txt.

    plot_normalized_metrics.py
        Computes and plots normalized stiffness, peak stress, and energy density
        as functions of the apex parameter a. Reads directly from f-d-all.txt.

    plot_intact_vs_precut_overlay.py
        Generates a 3x3 grid of subplots, one per design, each showing the
        intact and precut force-displacement curves overlaid. Requires
        f-d-all.txt and f-d-precut.txt.

    plot_energy_ratio.py
        Plots the ALLKE/ALLIE kinetic-to-internal energy ratio over normalised
        time for all nine symmetric designs. Requires energy_all.txt.


post-processing/cross-design/
    create_force_energy_comparison.py
        Generates grouped bar charts comparing peak force and work to peak
        across all design families (Baseline, Symmetric, G1, G2, G3).

    visualize_complete_results.py
        Generates stiffness and notch sensitivity bar charts covering all
        design families.

    discussion_pareto_plots.py
        Generates scatter plots of peak intact force versus notch sensitivity
        and stiffness versus notch sensitivity for the discussion chapter.

    compute_all_metrics.py
        Extracts peak force, initial stiffness, work to peak, and notch
        sensitivity index from raw force-displacement text files for all designs.

    compute_stiffness.py
        Computes initial stiffness from force-displacement data using a linear
        fit over the first 10 percent of the pre-peak displacement range.

    compute_postpeak_metrics.py
        Computes post-peak ductility index and softening slope for all
        configurations.


post-processing/quasi-static/
    plot_allke_allie_all.py
        Plots ALLKE/ALLIE energy ratio histories for all design groups to
        verify quasi-static conditions. Requires data files in energy-data/.

    plot_energy_ratio_all.py
        Alternative energy ratio visualization covering all designs on a
        single figure.


appendix/
    plot_analytical_model.py
        Generates figures for the analytical cross-verification appendix,
        comparing strain-energy, compliance-matrix, and force-balance methods.

    plot_cut_force_displacement.py
        Plots force-displacement curves for the cut placement sensitivity study
        (Appendix: cut placement). Reads from cut.txt.

    plot_size_effect.py
        Plots force-displacement curves and peak force scaling for the size
        effect study (Appendix: size effect). Reads from size-effect data files.


================================================================================
ABAQUS INPUT FILES (abaqus_inp_files/)
================================================================================

All simulations use Abaqus/Explicit 2023 with an elastic-perfectly plastic
material model and ductile damage (damage initiation at 1% equivalent plastic
strain, linear evolution to deletion at u_f = 0.01 mm).

Folder structure and naming conventions are described below.


baseline/
    Simulations for Baseline designs A and B with intact and precut specimens.

    Naming convention:
      [A|B]           Design variant (A or B tessellation)
      [001|01|1|15|2] Loading rate: 001=t0.01s, 01=t0.1s, 1=t1.0s, 15=t1.5s, 2=t2.0s
      [-c]            Suffix present for precut (notched) specimens

    Examples:
      A1.inp          Baseline A, intact, t=1.0s
      A1-c.inp        Baseline A, precut, t=1.0s
      B15.inp         Baseline B, intact, t=1.5s


cut_app/
    Simulations for the cut placement sensitivity study (Appendix).
    Tests three notch locations (left, middle, right) and three cut sizes
    (through 1, 2, or 3 struts).

    Naming convention:
      loc[1|2|3]      Notch location (1=right, 2=left, 3=middle)
      [1|2|3]strut    Number of struts cut through

    Example:
      loc2-1strut.inp    Left location, cut through 1 strut


G1/
    Simulations for the G1 graded design (horizontal gradient, partial range).
    Three tessellation variants A, B, C at multiple loading rates.

    Naming convention: same as baseline with variant letter prefix.
    Example:
      G1A1.inp        G1 variant A, intact, t=1.0s
      G1C1-c.inp      G1 variant C, precut, t=1.0s


G2/
    Simulations for the G2 graded design (horizontal gradient, full range).
    Two tessellation variants A and B.

    Naming convention: same as G1.
    Example:
      G2A1.inp        G2 variant A, intact, t=1.0s
      G2B1-c.inp      G2 variant B, precut, t=1.0s


G3/
    Simulations for the G3 graded design (vertical gradient, along loading direction).

    Naming convention: same as baseline.
    Example:
      G31.inp         G3 intact, t=1.0s
      G31-c.inp       G3 precut, t=1.0s


symmetric/
    Simulations for all nine symmetric designs (apex parameter a=0.1 to 0.9).
    Both intact and precut at multiple loading rates.

    Naming convention:
      [01|02|...|09]  Apex parameter (01=a0.1, 05=a0.5, 09=a0.9)
      [001|01|1|15|2] Loading rate (same as baseline)
      [-c]            Precut specimen

    Example:
      051.inp         a=0.5, intact, t=1.0s
      051-c.inp       a=0.5, precut, t=1.0s
      0915.inp        a=0.9, intact, t=1.5s


size_effect_app/
    Simulations for the size effect study (Appendix).
    Specimens of varying matrix size n x n.

    Naming convention:
      [n]x[n]         Matrix dimension (e.g. 3x3, 6x6, 9x9, 12x12, 16x16)

    Example:
      9x9.inp         9x9 unit cell specimen


================================================================================
