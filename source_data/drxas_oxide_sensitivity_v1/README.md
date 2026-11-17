# Oxide pilot radius sensitivity

Numerical sensitivity only. All three references remain candidate; no experimental validation or promotion.

Baseline SCF=5 A, cluster=8 A. Vary one control at a time. RPATH=4 A and all other native settings retained.
Signals use native degeneracy, S02=1, sigma2=0, E0=0, dR=0. Metrics use k^2 weighting.

Review triggers (empirical, predeclared): complex NRMS >5%, amplitude NRMS >3%, amplitude-weighted circular phase RMS >0.03 rad.
These are numerical review thresholds, not accuracy guarantees. Lower-radius probes are not convergence tests.

| Reference | Setting | k window | Max per-path complex NRMS | Max amplitude NRMS | Max phase RMS (rad) | Shell complex NRMS |
|---|---|---|---:|---:|---:|---:|
| IrO2-rutile-COD2101854 | scf6_cluster8 | 3-8 | 2.312% | 0.128% | 0.023105 | 2.305% |
| IrO2-rutile-COD2101854 | scf6_cluster8 | 8-12 | 1.413% | 0.032% | 0.014124 | 1.407% |
| IrO2-rutile-COD2101854 | scf6_cluster8 | 12-16 | 1.050% | 0.016% | 0.010499 | 1.042% |
| IrO2-rutile-COD2101854 | scf6_cluster8 | 16-20 | 0.905% | 0.033% | 0.009048 | 0.866% |
| IrO2-rutile-COD2101854 | scf6_cluster8 | 3-12 | 2.138% | 0.113% | 0.021971 | 2.143% |
| IrO2-rutile-COD2101854 | scf5_cluster9 | 3-8 | 0.026% | 0.001% | 0.000262 | 0.018% |
| IrO2-rutile-COD2101854 | scf5_cluster9 | 8-12 | 0.047% | 0.003% | 0.000474 | 0.033% |
| IrO2-rutile-COD2101854 | scf5_cluster9 | 12-16 | 0.005% | 0.000% | 0.000053 | 0.005% |
| IrO2-rutile-COD2101854 | scf5_cluster9 | 16-20 | 0.053% | 0.002% | 0.000528 | 0.022% |
| IrO2-rutile-COD2101854 | scf5_cluster9 | 3-12 | 0.032% | 0.002% | 0.000327 | 0.022% |
| IrO2-rutile-COD2101854 | scf4_cluster8 | 3-8 | 0.126% | 0.070% | 0.001047 | 0.122% |
| IrO2-rutile-COD2101854 | scf4_cluster8 | 8-12 | 0.049% | 0.037% | 0.000324 | 0.039% |
| IrO2-rutile-COD2101854 | scf4_cluster8 | 12-16 | 0.057% | 0.023% | 0.000518 | 0.039% |
| IrO2-rutile-COD2101854 | scf4_cluster8 | 16-20 | 0.111% | 0.024% | 0.001095 | 0.103% |
| IrO2-rutile-COD2101854 | scf4_cluster8 | 3-12 | 0.113% | 0.064% | 0.000963 | 0.109% |
| IrO2-rutile-COD2101854 | scf5_cluster7 | 3-8 | 0.090% | 0.016% | 0.000886 | 0.086% |
| IrO2-rutile-COD2101854 | scf5_cluster7 | 8-12 | 0.049% | 0.006% | 0.000485 | 0.024% |
| IrO2-rutile-COD2101854 | scf5_cluster7 | 12-16 | 0.072% | 0.012% | 0.000716 | 0.037% |
| IrO2-rutile-COD2101854 | scf5_cluster7 | 16-20 | 0.072% | 0.008% | 0.000715 | 0.065% |
| IrO2-rutile-COD2101854 | scf5_cluster7 | 3-12 | 0.081% | 0.015% | 0.000825 | 0.077% |
| beta-PtO2-COD1530633 | scf6_cluster8 | 3-8 | 2.046% | 0.121% | 0.020412 | 2.045% |
| beta-PtO2-COD1530633 | scf6_cluster8 | 8-12 | 1.248% | 0.047% | 0.012471 | 1.241% |
| beta-PtO2-COD1530633 | scf6_cluster8 | 12-16 | 0.947% | 0.013% | 0.009464 | 0.945% |
| beta-PtO2-COD1530633 | scf6_cluster8 | 16-20 | 0.782% | 0.050% | 0.007802 | 0.779% |
| beta-PtO2-COD1530633 | scf6_cluster8 | 3-12 | 1.890% | 0.109% | 0.019366 | 1.890% |
| beta-PtO2-COD1530633 | scf5_cluster9 | 3-8 | 0.021% | 0.001% | 0.000213 | 0.016% |
| beta-PtO2-COD1530633 | scf5_cluster9 | 8-12 | 0.000% | 0.000% | 0.000004 | 0.000% |
| beta-PtO2-COD1530633 | scf5_cluster9 | 12-16 | 0.000% | 0.000% | 0.000000 | 0.000% |
| beta-PtO2-COD1530633 | scf5_cluster9 | 16-20 | 0.000% | 0.000% | 0.000000 | 0.000% |
| beta-PtO2-COD1530633 | scf5_cluster9 | 3-12 | 0.019% | 0.001% | 0.000194 | 0.014% |
| beta-PtO2-COD1530633 | scf4_cluster8 | 3-8 | 3.630% | 0.577% | 0.035743 | 3.616% |
| beta-PtO2-COD1530633 | scf4_cluster8 | 8-12 | 1.670% | 0.287% | 0.016434 | 1.665% |
| beta-PtO2-COD1530633 | scf4_cluster8 | 12-16 | 1.165% | 0.169% | 0.011516 | 1.156% |
| beta-PtO2-COD1530633 | scf4_cluster8 | 16-20 | 0.707% | 0.026% | 0.007063 | 0.707% |
| beta-PtO2-COD1530633 | scf4_cluster8 | 3-12 | 3.279% | 0.524% | 0.033405 | 3.270% |
| beta-PtO2-COD1530633 | scf5_cluster7 | 3-8 | 0.154% | 0.038% | 0.001491 | 0.152% |
| beta-PtO2-COD1530633 | scf5_cluster7 | 8-12 | 0.068% | 0.021% | 0.000654 | 0.059% |
| beta-PtO2-COD1530633 | scf5_cluster7 | 12-16 | 0.057% | 0.011% | 0.000557 | 0.042% |
| beta-PtO2-COD1530633 | scf5_cluster7 | 16-20 | 0.083% | 0.023% | 0.000794 | 0.072% |
| beta-PtO2-COD1530633 | scf5_cluster7 | 3-12 | 0.138% | 0.034% | 0.001394 | 0.136% |
| B-Ta2O5-COD2102123 | scf6_cluster8 | 3-8 | 0.116% | 0.037% | 0.001102 | 0.208% |
| B-Ta2O5-COD2102123 | scf6_cluster8 | 8-12 | 0.066% | 0.023% | 0.000621 | 0.025% |
| B-Ta2O5-COD2102123 | scf6_cluster8 | 12-16 | 0.051% | 0.021% | 0.000472 | 0.030% |
| B-Ta2O5-COD2102123 | scf6_cluster8 | 16-20 | 0.158% | 0.136% | 0.000802 | 0.080% |
| B-Ta2O5-COD2102123 | scf6_cluster8 | 3-12 | 0.102% | 0.034% | 0.001008 | 0.118% |
| B-Ta2O5-COD2102123 | scf5_cluster9 | 3-8 | 0.055% | 0.022% | 0.000510 | 0.117% |
| B-Ta2O5-COD2102123 | scf5_cluster9 | 8-12 | 0.058% | 0.005% | 0.000577 | 0.052% |
| B-Ta2O5-COD2102123 | scf5_cluster9 | 12-16 | 0.063% | 0.003% | 0.000630 | 0.017% |
| B-Ta2O5-COD2102123 | scf5_cluster9 | 16-20 | 0.082% | 0.032% | 0.000758 | 0.047% |
| B-Ta2O5-COD2102123 | scf5_cluster9 | 3-12 | 0.055% | 0.019% | 0.000514 | 0.078% |
| B-Ta2O5-COD2102123 | scf4_cluster8 | 3-8 | 1.940% | 0.127% | 0.019372 | 2.021% |
| B-Ta2O5-COD2102123 | scf4_cluster8 | 8-12 | 1.122% | 0.034% | 0.011215 | 0.968% |
| B-Ta2O5-COD2102123 | scf4_cluster8 | 12-16 | 0.860% | 0.033% | 0.008592 | 0.760% |
| B-Ta2O5-COD2102123 | scf4_cluster8 | 16-20 | 0.855% | 0.214% | 0.008489 | 0.733% |
| B-Ta2O5-COD2102123 | scf4_cluster8 | 3-12 | 1.776% | 0.112% | 0.018272 | 1.386% |
| B-Ta2O5-COD2102123 | scf5_cluster7 | 3-8 | 0.075% | 0.028% | 0.000691 | 0.139% |
| B-Ta2O5-COD2102123 | scf5_cluster7 | 8-12 | 0.055% | 0.007% | 0.000542 | 0.042% |
| B-Ta2O5-COD2102123 | scf5_cluster7 | 12-16 | 0.060% | 0.004% | 0.000601 | 0.043% |
| B-Ta2O5-COD2102123 | scf5_cluster7 | 16-20 | 0.048% | 0.017% | 0.000466 | 0.038% |
| B-Ta2O5-COD2102123 | scf5_cluster7 | 3-12 | 0.069% | 0.025% | 0.000659 | 0.085% |

## Baseline repeatability

[
  {
    "reference": "IrO2-rutile-COD2101854",
    "compared_files": 7,
    "all_input_and_raw_path_bytes_identical": true,
    "executable_fingerprints_match_source": true,
    "package_versions_match_source": true,
    "differences": []
  },
  {
    "reference": "beta-PtO2-COD1530633",
    "compared_files": 9,
    "all_input_and_raw_path_bytes_identical": true,
    "executable_fingerprints_match_source": true,
    "package_versions_match_source": true,
    "differences": []
  },
  {
    "reference": "B-Ta2O5-COD2102123",
    "compared_files": 23,
    "all_input_and_raw_path_bytes_identical": true,
    "executable_fingerprints_match_source": true,
    "package_versions_match_source": true,
    "differences": []
  }
]

## Failures

[]

The native radius-counter warning and underflow/denormal flags are preserved. Successful modules and stopping before the SCF iteration cap do not prove radius convergence.

Primary methodological references:
- https://feff.phys.washington.edu/FEFF/Docs/feff8/feff84/feff84/SCF.html
- https://feff.phys.washington.edu/feff/Docs/feff8/feff8web/node49.html
- https://xraypy.github.io/xraylarch/xafs_feffpaths.html
