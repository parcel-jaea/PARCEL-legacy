English | [日本語](README.ja.md)

# PARCEL (legacy versions 1.0 – 1.2)

This repository keeps the past versions of PARCEL (PARallel Computing ELements), the parallel Krylov subspace
solver library written in Fortran (and CUDA from version 1.2). They are no longer developed.

**The current PARCEL is the SYCL version (PARCEL_SYCL):** <https://github.com/parcel-jaea/PARCEL-release>
(releases: <https://github.com/parcel-jaea/PARCEL-release/releases>; manual: <https://ccse2.jaea.go.jp/software/parcel/index_en.html>).

## Versions

Each version has a tag and a GitHub Release with its source tarball and SHA-256 file.
The `main` branch holds the last version (1.2); the files of each version are under its tag.

| Version | Released | Tag / download | Manual (HTML) |
|---|---|---|---|
| 1.2 | 2023-01-31 | [v1.2](https://github.com/parcel-jaea/PARCEL-legacy/releases/tag/v1.2) | [Japanese](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.2/manual/PARCEL_manual_mathjax.html) / [English](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.2/manual/PARCEL_manual_mathjax_eng.html) |
| 1.1 | 2022-01-13 | [v1.1](https://github.com/parcel-jaea/PARCEL-legacy/releases/tag/v1.1) | [Japanese](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.1/manual/PARCEL_manual_mathjax.html) / [English](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.1/manual/PARCEL_manual_mathjax_eng.html) |
| 1.01 | 2020-05-20 | [v1.01](https://github.com/parcel-jaea/PARCEL-legacy/releases/tag/v1.01) | [Japanese](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.01/manual/PARCEL_manual_mathml.html) / [English](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.01/manual/PARCEL_manual_mathml_eng.html) |
| 1.0 | 2019-03-29 | [v1.0](https://github.com/parcel-jaea/PARCEL-legacy/releases/tag/v1.0) | [Japanese](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.0/manual/PARCEL_manual_mathml.html) / [English](https://github.com/parcel-jaea/PARCEL-legacy/blob/v1.0/manual/PARCEL_manual_mathml_eng.html) |

The release dates are those announced on the page of past versions on the distribution site (the page is being closed).
How to build and use each version is in the `readme` files of that version.

The manual of each version is in `manual/` of its tarball and tag (each is a single HTML file with the figures and the math embedded, and loads nothing from the Internet).
Open it in a web browser after extracting the tarball. GitHub shows an HTML file as source code.

## Differences from the tarballs on the distribution site

The tarballs here have the same file names and top directories as those formerly on the distribution site, with these changes only:

- Build products and work files are removed (object files `*.o`, module files `*.mod`, executables, a job output file,
  and in 1.0 `SPARSE/src/tmp.f90` and `tmp.sh`).
- Lines containing personal paths of developers' accounts are removed (they were all commented out). In 1.0,
  the path of the library in `SPARSE/example5/Makefile` is replaced by `<path-to>/libjaeri_ilesf77.a`, and the working directory
  in the sample output `SPARSE/example/STDOUT_PARCEL-dev` by `<path-to>/PARCEL-dev`.
- The manual (HTML, Japanese and English) is added in `manual/`. It was distributed separately on the distribution site.
  In 1.1 and 1.2 the math, which was displayed by MathJax loaded from the Internet, is converted to MathML (the text is unchanged; the file names are kept).
- The license texts `COPYING` and `COPYING.LESSER` are added. In 1.0, `MPI/lib/.gitkeep` is added to keep the empty directory.
- File owners are not recorded, and all timestamps are set to the date of the version.

## License

GNU Lesser General Public License v3.0 (LGPL-3.0), the same as the distribution site. See `COPYING` and `COPYING.LESSER`.

## Contact

Questions: GitHub Issues of PARCEL_SYCL, <https://github.com/parcel-jaea/PARCEL-release/issues>.
These versions are not maintained; please use PARCEL_SYCL for new work.
