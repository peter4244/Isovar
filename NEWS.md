# isovar

## 0.0.2 (unreleased)

- `loadCopdGwas()` — reads icgcUkb COPD summary statistics.
- `annotateGwas()` — rsID-based (build-independent) join of GWAS effect estimates onto a variant table. Liftover deliberately not wired; add only if rsID coverage proves insufficient.
- `checkAlleleAlignment()` — classifies effect-allele orientation (`aligned` / `swapped` / `strand_flip` / `mismatch`) without silently flipping beta.

## 0.0.1 (unreleased)

Initial scaffold.

- `rankSplaireVariants()` — rank credible-set variants by splaire / splaireVar delta across donor / acceptor / splice-site-usage heads.
- `annotateGnomad()` — attach rsID and gnomAD v4.1 population allele frequencies via remote tabix over AWS Open Data (region-batched).
- `implicatedSpliceSites()` — derive structured splice-site objects from splaire rows, detecting the "motif shifted by indel length" pattern.
- First runs prototyped against AKR1A1 credible-set variants from HAEC-185.
