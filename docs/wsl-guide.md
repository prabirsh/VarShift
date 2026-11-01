# Running VarSight on Windows with WSL

Follow these steps in order. Commands starting from step 2 go in your **WSL terminal** (Ubuntu), not PowerShell.

## 1. Open WSL

Press the Windows key, type `Ubuntu` and press Enter. If it says `Failed to start the systemd user session`, you can ignore that: VarSight doesn't need it.

Check Python:

```bash
python3 --version
```

You need 3.10 or newer. If it's missing:

```bash
sudo apt update && sudo apt install -y python3
```

## 2. Get the code

If you are using the copy already on your `I:` drive:

```bash
cd /mnt/i/VIP-Annotations/VarSight
```

Otherwise, once it's on GitHub:

```bash
git clone https://github.com/YOUR-USERNAME/VarSight.git ~/VarSight
cd ~/VarSight
```

## 3. Install and run the demo

```bash
./install.sh --test
```

This needs no sudo and no download. It takes about two minutes, and it should end with `OK` for the tests, the demo finishing, and `Installed`.

If it says `Permission denied`, run `sh install.sh --test` instead.

If it tells you to add the command to your PATH, run the line it prints. Then check:

```bash
varsight --version
```

Look at the demo results in `demo_output/`. To open them in Windows, run `explorer.exe demo_output`.

## 4. Choose where the data goes

**Use the WSL disk (`~/varsight-data`).** Lookups measured about 50 to 65 times faster there than on a Windows drive.

| | `~/varsight-data` | `/mnt/i/varsight-data` |
|---|---|---|
| Lookup speed | 27 microseconds | 1.3 to 1.8 ms |
| Physically on | Your C: drive (32 GB free when last checked) | Your I: drive (87 GB free) |
| File privacy modes | Work | May not apply |

Check you have room. You need about 25 GB free during setup:

```bash
varsight doctor --data-dir ~/varsight-data
```

## 5. Download and build the data

This is the long step: about 8 GB of transfer plus processing. My estimate is 1 to 2 hours on a typical connection, but this hasn't been timed yet.

Plug the laptop in, and keep it awake. Then run:

```bash
varsight setup --data-dir ~/varsight-data
```

**If it stops** because the connection dropped or the laptop slept, run the same command again. It resumes partial downloads, and skips anything already built.

To get only one build first:

```bash
varsight setup --data-dir ~/varsight-data --build hg19
```

To save 4.7 GB of transfer, at the cost of domain positions:

```bash
varsight setup --data-dir ~/varsight-data --skip-pfam-regions
```

At the end it runs the known-variant self-test. You want to see `11 of 11 known variants passed`.

**If a known variant fails or setup stops with an error,** please keep the message and report it. This full run hasn't been done before, so a problem here is useful to know about.

## 6. Annotate

```bash
varsight annotate "BRAF V600E" --build hg19
varsight annotate "NM_004333.6:c.1799T>A" --build hg38 --format tsv
varsight annotate "chr7:g.140453136A>T" --build hg19 --out one.xlsx
varsight annotate variants.txt --build hg19 --out list.tsv
```

For your VCF on the Windows drive:

```bash
varsight annotate /mnt/i/DCG_DRAGEN_Data/SAMPLE/SAMPLE.hard-filtered.vcf.gz --out ~/results/SAMPLE.xlsx
```

It reads the build from the VCF header. Keep results with patient data on the WSL disk, where the file privacy modes apply. To copy one to Windows:

```bash
explorer.exe ~/results
```

## 7. Measure accuracy on your own data

These two commands give you the proof for your own installation:

```bash
varsight validate --sample 20000 --report ~/varsight-validation.md
python3 tools/compare_annotated_vcf.py /mnt/i/DCG_DRAGEN_Data/SAMPLE/SAMPLE.hard-filtered.vcf.annotated.vcf.gz --differences ~/dragen-differences.txt
```

- `validate` tests 20,000 ClinVar names per build across the whole genome, and lists every mismatch group.
- `compare_annotated_vcf.py` compares against your DRAGEN annotation. It stays on your computer, but treat its output as private.

## 8. Weekly updates

Check that updating works:

```bash
varsight data update --check
```

To make it automatic, open **Windows PowerShell** (not WSL), and run this with your own WSL user name:

```powershell
cd I:\VIP-Annotations\VarSight\scripts\windows
powershell -ExecutionPolicy Bypass -File .\register-weekly-update.ps1 -RepoPath /mnt/i/VIP-Annotations/VarSight -DataDir /home/YOUR-WSL-USER/varsight-data -VarsightCommand /home/YOUR-WSL-USER/.local/bin/varsight
```

Try it straight away:

```powershell
Start-ScheduledTask -TaskName 'VarSight weekly data update'
```

Then, back in WSL, see what it did:

```bash
ls -t ~/varsight-data/reports/ | head
```

More in [updates.md](updates.md).

## 9. Publish to GitHub

Your repository has no remote yet. Before publishing:

1. Replace `YOUR-USERNAME` with your GitHub user name:

    ```bash
    grep -rl "YOUR-USERNAME" --exclude-dir=.git . | xargs sed -i 's/YOUR-USERNAME/your-real-name/g'
    ```

2. Commit and push. Git needs your name once:

    ```bash
    git config user.name "Your Name"
    git config user.email "you@example.com"
    git add -A
    git add --chmod=+x install.sh scripts/weekly-update.sh tools/build_demo_data.py tools/check_sources.py tools/compare_annotated_vcf.py
    git status
    git commit -m "Add VarSight offline variant annotation"
    git remote add origin https://github.com/your-real-name/VarSight.git
    git push -u origin main
    ```

    **Read the `git status` list before committing.** No patient VCF or results file should appear in it.

    The `--chmod` line matters: this repository was created on Windows, so without it `./install.sh` would say `Permission denied` after someone clones it.

3. On GitHub, go to **Settings**, then **Code security**, and turn on **Private vulnerability reporting**.

CI then runs on Python 3.10 to 3.14 automatically.

## Problems

| Message | What to do |
|---|---|
| `varsight: command not found` | Run the PATH line that `install.sh` printed, or use `~/.local/bin/varsight` |
| `no VarSight data in ...` | Add `--data-dir ~/varsight-data`, or run setup |
| `not enough disk space` | Free up space on C:, or choose `--data-dir /mnt/i/varsight-data`, which is slower |
| `another setup or update is running` | Wait for it. If nothing is running, run the command again: a leftover lock from a stopped run clears itself |
| `REF ... does not match` | Usually the wrong build. Check `--build` |
| `--build says hg38, but the contig lengths ... are hg19` | Your VCF is on the other build |
| Very slow lookups | Check that `--data-dir` isn't under `/mnt/`; `varsight doctor` warns about it |
| Out of disk on C: after removing data | WSL's virtual disk doesn't shrink on its own. Search for "compact WSL vhdx" |
| `expected columns are missing` | A source changed its format. Nothing was installed, and your current data is safe. Please report it |
| Git says `dubious ownership` | Run `git config --global --add safe.directory /mnt/i/VIP-Annotations/VarSight` |
