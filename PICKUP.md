# Pickup prompt — ARS89 reproduction

Paste the block below into a fresh Claude Code session (after restart) to resume.

---

Continue the project in this directory: reproduce the numerical code from
**Adams, Ruden & Shu 1989, ApJ, 347, 959, "Eccentric Gravitational Instabilities
in Nearly Keplerian Disks"** (bibcode 1989ApJ...347..959A) in Python, using best
scientific-programming practice (vectorized array operations, no naive Python
loops over grid points, clear separation of physics / discretization / solver,
typed and documented).

The paper treats the m=1 (one-armed, eccentric) gravitational instability of a
razor-thin, nearly-Keplerian disk around a central star — the SLING mechanism,
where motion of the central star (the indirect/m=1 potential term) matters. The
core numerical task is a linear normal-mode eigenvalue problem: linearize the 2D
fluid equations (continuity + momentum) for perturbations ∝ exp(i(mφ − ωt)) with
m=1, close them with the perturbed self-gravity (Poisson integral / softened
kernel), and solve the resulting radial integro-differential equation for the
complex eigenfrequency ω and eigenfunction, subject to the paper's boundary
conditions.

Step 1 — get the actual equations before writing code:
  - Try the scanned full text: WebFetch https://articles.adsabs.harvard.edu/pdf/1989ApJ...347..959A
    (now allowlisted). It is a SCAN — if OCR comes back empty/garbled, say so.
  - Fallbacks if the scan is unreadable: (a) the companion paper Shu, Tremaine,
    Adams & Ruden 1990, ApJ 358, 495 "Sling Amplification and Eccentric
    Gravitational Instabilities in Gaseous Disks" (bibcode 1990ApJ...358..495S)
    develops the same formalism more explicitly; (b) ASK ME — I can drop a local
    PDF of ARS89 into this directory.
  - Do NOT reconstruct the equations from memory alone; confirm them against a
    source and tell me which equations you are implementing.

Step 2 — propose a short implementation plan (enter plan mode) covering: governing
equations, radial discretization, how the self-gravity kernel is handled, the
eigenvalue solver, boundary conditions, and a validation check (e.g. recover a
known limit or a figure/growth-rate from the paper). Then implement.

Environment notes:
  - .claude/settings.json already grants acceptEdits + Bash/Edit/Write/WebSearch
    and WebFetch for arxiv, IOP, OUP, A&A, and ADS hosts (incl. articles.adsabs).
  - Not a git repo yet — offer `git init` for checkpoint/rewind before long runs.
  - Once the task is well-defined, /loop can self-pace the build+validate cycle.
