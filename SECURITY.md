# Security

VarSight handles genomic variants, which can be patient data. This page says what it protects against, how, and what it leaves to you.

## Reporting a problem

Please use GitHub's private vulnerability reporting on this repository (**Security** tab, then **Report a vulnerability**). Do not open a public issue for a security problem, and do not attach patient data to any report.

Only the latest release gets security fixes.

## What VarSight does

| Risk | Control | Checked by |
| - | - | - |
| Patient variants leaving the computer | Annotation has no network code. Only `pipeline/fetch.py` can open a connection, and it runs only during setup and update. A variant leaves only if a person clicks a link in the results. | `tests/test_security.py` reads every source file and fails if any other module imports a network library |
| The update check leaking anything | It runs apart from annotation, as a separate process started with fixed arguments, and takes only the data folder. It sends header requests for public files, with nothing about any variant or input file. `varsight config updates=off` turns it off | `tests/test_security.py`, `tests/test_commandline.py` |
| A licensed database reaching other people | COSMIC, OMIM, SpliceAI and PrimateAI are never downloaded for you. You supply your own file, and it stays on your computer | `tests/test_advanced.py` |
| A hostile VCF attacking whoever opens the HTML report | All text reaches the page as text, never as markup. The data can't close its script tag. Only `https` addresses become links, and they open with `noopener` | `tests/test_output.py` |
| Hidden telemetry | There is none. Setup and update requests carry only a version string. | Same test, plus the code is short enough to read |
| A tampered or corrupted download | HTTPS with certificate checks. MD5 compared where the source publishes one (UCSC, ClinVar, REVEL on Zenodo). gzip CRC checked on every archive. A size ceiling for each source. SHA-256 of every installed file is recorded, and `varsight doctor verify` compares against it. | `tests/test_fetch.py` |
| A download redirected somewhere else | A fixed list of allowed hosts. Redirects go through the same check as the first address. | `tests/test_fetch.py` |
| A bad update replacing good data | New data is built beside the current data, must not shrink by more than 10%, and must pass the known-variant self-test before it goes live. Any failure keeps the current data. The previous version is kept for `varsight rollback`. | `tests/test_manager.py` |
| SQL injection | Every query uses parameters. Databases open read-only with `trusted_schema=OFF`, and extension loading is never enabled. | `tests/test_security.py`, `tests/test_inputs.py` |
| Spreadsheet formula injection from a VCF | Excel cells are written as inline strings, so they can never be formulas. TSV cells starting with `=`, `+`, `@` or `-letter` get a leading `'`. Only `https` values become hyperlinks. | `tests/test_output.py` |
| Malformed input | Input length is capped, patterns are anchored, alleles must be plain DNA, and VCF lines have a length limit. Errors are reported per record without stopping the file. | `tests/test_inputs.py`, `tests/test_integration.py` |
| Code execution | No `eval`, `exec`, `pickle` or shell calls anywhere in `src/`. VarSight starts only two processes itself, both with fixed arguments and no shell: the update check in `updates.py`, and on Windows, when Python refuses a download's certificate, one visit by `curl.exe` from System32 to that allow-listed host, because Windows adds a trusted root only when its own client asks. Python then checks the certificate as strictly as before. Worker processes for large VCFs come from the standard `multiprocessing` module. | `tests/test_security.py` |
| Results exposed to other users | Output files are created readable only by their owner (mode 600). The data folder is mode 700. Output is never written through a symbolic link, and an existing file is not replaced without the word `replace`. | `tests/test_output.py` |
| Supply-chain attack through dependencies | There are no runtime dependencies: only the Python standard library. GitHub Actions are pinned to full commit hashes, and Dependabot proposes updates. | `tests/test_security.py`, CI |
| Patient data committed to this repository | `.gitignore` excludes VCF archives, results and data folders. | Review |

## What you still have to do

- **Treat result files as patient data.** VarSight keeps them private on disk, but where you copy them is your decision.
- **Know that links leave your computer.** Opening a Franklin, VarSome or OncoKB link sends that variant to that website. VarSight never opens links itself.
- **Encrypt the disk** if the computer holds patient data. VarSight does not encrypt files.
- **Validate before clinical use.** See [docs/validation.md](docs/validation.md). VarSight gives no clinical verdict.
- **Check licences for your use.** See [docs/data-sources.md](docs/data-sources.md).
- **Mode 600 needs a Linux filesystem.** On a Windows drive seen from WSL (`/mnt/c`, `/mnt/i`), file modes may not apply, so keep patient results on the Linux disk. On native Windows these modes do nothing at all: privacy comes from the folder you choose, so keep data and results under your own user profile.

## Not covered

- An attacker who already controls your account can read anything you can.
- A source's own server being compromised. HTTPS and checksums prove a file came intact from that server, not that its content is right. The self-test and the shrink check catch gross damage only.
- Sources that publish no checksum (RefSeq, GENCODE, HGNC, UniProt, Pfam, CIViC, ClinGen, gnomAD, HPO) rely on HTTPS, the gzip CRC where compressed, and the recorded SHA-256.
