# paneltests 1.0.6

* Corrected the DOIs of Galvao, Juhl, Montes-Rojas and Olmo (2017) (10.1093/jjfinec/nbx016), Little (1988) (10.1080/01621459.1988.10478722) and Swamy (1970) (10.2307/1913012).
* Demetrescu, Hosseinkouchack and Rodrigues (2023): title corrected to the published working-paper title and DOI corrected to 10.4419/96973210 (Ruhr Economic Papers 1041); the journal DOI previously given in DESCRIPTION pointed to an unrelated article and was removed.
* No changes to code.

# paneltests 1.0.5

- Resubmission addressing CRAN reviewer (Benjamin Altmann) feedback on 1.0.4:
  - Added DOI / ISBN links to all references in the Description field (Hsiao 2014, Swamy 1970, Pesaran 2004).
  - Removed the surrounding single quotes around function names in the Description text and added trailing parentheses, e.g. `'xtpretest'`-&gt; `xtpretest()`.
- Synchronised CITATION version with DESCRIPTION (1.0.5).
- No code or behavior changes.

# paneltests 1.0.4

- Resubmission to CRAN.
- Removed UTF-8 BOM (Byte Order Mark) from `inst/CITATION` and from `R/xtcsdq.R`, `R/xtmispanel.R`, `R/xtqsh.R`. The BOM in `inst/CITATION`was causing the CRAN incoming auto-check to fail with: `Invalid citation information in 'inst/CITATION': 1:1: unexpected input`.
- Synchronised the version string in `inst/CITATION` with `DESCRIPTION`.
- No code or behavior changes.

# paneltests 1.0.3

- Resubmission. Removed redundant "a suite of tools for" from the Description start. All references use `doi:10.xxx` format.

# paneltests 1.0.0

- Initial CRAN release.
- Consolidates four previously separate packages (xtpretest, xtmispanel, xtcsdq, xtqsh) into a single panel-data pre-testing suite.
- Implements:
  - Hsiao (2014) homogeneity F-tests (F1/F2/F3),
  - Swamy (1970) parameter heterogeneity test,
  - Pesaran (2004) cross-sectional dependence test (`xtpretest`),
  - Missing-data detection, mechanism testing, and imputation for unbalanced panels (`xtmispanel`),
  - Demetrescu, Hosseinkouchack and Rodrigues (2023) quantile-regression cross-sectional dependence tests (`xtcsdq`),
  - Galvao, Juhl, Montes-Rojas and Olmo (2017) panel quantile-regression slope homogeneity statistics (`xtqsh`).
