# 微妙最微妙,甚深最甚深

## This repo implements Grassmann time evolving matrix product operator (GTEMPO) method for fermionic impurity problem, which includes a series of algorithms which may or may not be explicitly mentioned in our papers in recent years. Another closely related package, which uses infinite MPS technique to fully explore the time translational invariance of the IF, will also be open sourced later.

## List of algorithms (tricks) being implemented
- GTEMPO on the Keldysh contour (PRB 109, 045140 (2024))
- GTEMPO on the imaginary contour (NJP 26, 013109 (2024))
- GEMPO on the L-shaped Kadanoff contour PRB 110, 165114 (2024)()
- Fermionic QuAPI (in the QuAPI package)
- The zipup algorithm which only build the augmented density tensor on the fly (PRB 109, 045140 (2024))
- The partial-IF algorithm, which builds the MPS-IF using O(N) GMPS multiplications, each partial IF is analytically built as a GMPS with bond dimension 2 (SciPost Phys. Core 7, 063 (2024))
- The TTI-IF algorithm, which efficiently builds the MPS-IF using O(1) GMPS multiplications, by exploring the time-translational invariance. (SciPost Phys. Core 7, 063 (2024))
- Exact TTI-IF algorithm, which efficiently builds the MPS-IF using O(k) GMPS multiplications, each with bond dimension 2, where k is the number of terms in Prony algorithm (or other algorithms to approximate the hybridization as sum of exponentials), Exact TTI-IF is guaranteed to be more accurate and efficient than TTI-IF algorithm (arXiv:2510.11459, 2025)
- Measuring the current by converting the bath correlations into impurity green's functions, and than using QuAPI-similar discretization and building the current operator as a GMPS with bond dimension 2, which then allows to efficiently calculate the final observable (PRB 109, 045140 (2024))
- Deal with retarded interaction (i.e., electron-phonon interaction) (Chinese Physics Letters 42, 120701 (2025))


## Full documentation is on the way...
