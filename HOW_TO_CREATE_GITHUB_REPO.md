# Creating the GitHub repository from scratch

These steps assume you don't have the repository set up locally yet (e.g. your
old local clone is stale or you want a clean start). If you already have a
local clone of `AliJnadi/Quadratic-Stability-Control`, skip to Option B.

## Option A — Brand new repository

### 1. Create the repository on GitHub
1. Go to [github.com/new](https://github.com/new).
2. Repository name: `Quadratic-Stability-Control`.
3. Set visibility (public/private) as you prefer.
4. **Do not** initialize with a README, .gitignore, or license — we already
   have those (or will add them), and letting GitHub create its own would
   just cause a merge conflict on the first push.
5. Click **Create repository**. GitHub will show you a page with setup
   instructions and a URL like `https://github.com/AliJnadi/Quadratic-Stability-Control.git`.

### 2. Initialize the local folder and push
Unzip the delivered `Quadratic-Stability-Control_final.zip` somewhere on your
machine, then from inside that folder:

```bash
cd Quadratic-Stability-Control
git init
git add .
git commit -m "Initial commit: QSC framework, experiments, and results"
git branch -M main
git remote add origin https://github.com/AliJnadi/Quadratic-Stability-Control.git
git push -u origin main
```

If you use SSH instead of HTTPS for GitHub, use
`git@github.com:AliJnadi/Quadratic-Stability-Control.git` for the remote URL
instead.

### 3. Verify
Refresh the GitHub page — you should see all files, and the `README.md`
should render automatically on the repository's front page.

## Option B — You already have a local clone

If you have an existing local clone with its own `.git` history that you
want to keep (recommended, so you don't lose commit history), copy the
*contents* of the delivered folder into your existing clone, overwriting
files with the same name, then commit normally:

```bash
# from inside your EXISTING local clone
cp -r /path/to/delivered/Quadratic-Stability-Control/* .
cp -r /path/to/delivered/Quadratic-Stability-Control/.[^.]* . 2>/dev/null  # hidden files, if any

git status                     # review what changed
git add .
git commit -m "Fix check_stability_samples bug, add Computational Environment \
disclosure, add validation/scaling experiments, clean up README and \
references, GitHub-readiness pass"
git push
```

## Recommended repository settings (optional but useful)

- **About section** (gear icon next to "About" on the repo page): add a short
  description ("Sampling-based Quadratic Stability Control for nonlinear
  mechanical systems with parametric uncertainty") and topics like
  `robust-control`, `lmi`, `semidefinite-programming`, `matlab`, `cvx`.
- **License**: MDPI/`Actuators` articles are typically published under
  CC BY 4.0; if you want the *code* under a permissive license (MIT is
  common for research code), add a `LICENSE` file via GitHub's
  "Add file → Create new file → LICENSE" template picker. The README
  currently says "see the repository's license file" — add one so that
  sentence has something to point to, or edit the README if you'd rather
  state the license inline.
- **Releases**: once the paper is accepted, consider tagging the exact
  commit used to produce the published results (`git tag v1.0-actuators` /
  `git push --tags`) so the paper's GitHub link always resolves to the
  reviewed, reproducible code state even if you continue developing on
  `main` afterward.
- **Citation file**: GitHub recognizes a `CITATION.cff` file and adds a
  "Cite this repository" button. Worth adding once the paper has a DOI.

## After pushing: sanity-check reproducibility

Clone the repository fresh into a *different* folder and try running one
sweep end-to-end, to confirm nothing was accidentally left out (e.g. a file
only present in your MATLAB path but not actually in the repo):

```bash
git clone https://github.com/AliJnadi/Quadratic-Stability-Control.git /tmp/qsc-check
cd /tmp/qsc-check/"Matlab Codes"
# in MATLAB:
%   Bulk_Test_Samples
%   cd tests/N, make_all_plots
```
If this reproduces `tests/N/summary_N_sweep.csv` and the figures without
errors, the repository is self-contained and reproducible.
