# Contributing

## Set up

Only Python 3.10 or newer is needed.

```bash
git clone https://github.com/YOUR-USERNAME/VarSight.git
cd VarSight
./install.sh --test
```

## Before you open a pull request

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -t .
```

All tests must pass. A change to HGVS naming, normalisation or the resolver also needs:

```bash
varsight validate --sample 20000 --report validation.md
```

Say in the pull request how the rates moved, and why any group got worse.

## Rules that keep this project trustworthy

1. **No new runtime dependencies.** Standard library only. This is a security and install decision, and `tests/test_security.py` enforces it.
2. **Network code stays in `pipeline/fetch.py`.** Annotation must work with no connection.
3. **Report what sources say.** VarSight does not merge sources into one verdict, and missing data is never treated as benign.
4. **A warning is better than a guess.** If a result might be wrong, add a flag from [the list](docs/input-output.md#flags) instead of hiding it.
5. **Every new source needs a licence note** in `src/varsight/data/sources.json`, and an entry in [docs/data-sources.md](docs/data-sources.md).
6. **Never commit patient data**, even as a test case. Tests use the hand-written gene in `tests/helpers.py` and the public data in `demo/`.

## Where things are

See [docs/architecture.md](docs/architecture.md).
