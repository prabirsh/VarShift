# Running VarSight on Ubuntu, WSL and Windows

VarSight runs the same way on all of these. Only the install step differs. Pick your system in step 1, and the rest of the guide is the same, apart from the notes marked for Windows.

## 1. Install

You need Python 3.10 or newer. Nothing else: no pip packages, no administrator rights.

### Ubuntu or WSL

In WSL, press the Windows key, type `Ubuntu` and press Enter. If it says `Failed to start the systemd user session`, you can ignore that.

```bash
python3 -V
```

If that shows an older version or nothing:

```bash
sudo apt update && sudo apt install -y python3
```

Then get the code and install. In WSL, keep the code in your Linux home folder rather than under `/mnt/`, because it starts faster there:

```bash
git clone https://github.com/prabirsh/VarSight.git ~/VarSight
cd ~/VarSight
sh install.sh
```

If it tells you to add the command to your PATH, run the line it prints.

### Windows PowerShell

Install Python from [python.org](https://www.python.org/downloads/) and tick **Add python.exe to PATH**. Then download the repository as a ZIP from GitHub and unpack it, or use `git clone`. In PowerShell, from the VarSight folder:

```powershell
powershell -ExecutionPolicy Bypass -File .\install.ps1
```

### Windows Command Prompt

The same Python, then from the VarSight folder:

```text
install.bat
```

On Windows it asks before adding the command to your PATH. **Open a new window afterwards** so the change takes effect.

### What the install does

It writes one small `varsight` launcher, runs the [demo](../demo/README.md) with no download, and offers to set up the databases. Check it:

```bash
varsight version
varsight doctor
```

Open `demo_output/variants/report.html` in a browser to see what a report looks like. It uses real data for 12 genes.

| System | To open the demo report |
| - | - |
| WSL | `explorer.exe demo_output` |
| Ubuntu desktop | `xdg-open demo_output/variants/report.html` |
| Windows | Double-click it, or `start demo_output\variants\report.html` |

## 2. Set up the databases

```bash
varsight setup
```

It asks three things, shows the size, and then waits for your yes before downloading:

```text
Which genome build do you work with?
  1  hg19 (GRCh37)
  2  hg38 (GRCh38)
  3  Both. Results then always show the position on both builds
Your choice, 1, 2 or 3: [3]

Where should the data go? [/home/you/varsight-data]

Advanced databases are optional, and can be added any time later.
   1  dbSNP                              26 GB per build, read without saving
   ...
Numbers to add now, separated by commas, or press Enter for none:

This will install: hg19 into /home/you/varsight-data
  Download   about 2.3 GB, plus 4.7 GB of Pfam read without saving
  Disk       about 5 GB kept, and keep about 15 GB free during setup (estimates)
Start? yes or no: [yes]
```

The same thing without questions, for scripts:

```bash
varsight setup build=hg38 data=/data/varsight yes
```

**Where to put the data:**

| System | Use | Why |
| - | - | - |
| WSL | `~/varsight-data`, on the Linux disk | Lookups measured 27 microseconds there, against 1.3 to 1.8 ms on a Windows drive seen from WSL (`/mnt/c`, `/mnt/i`), so about 50 to 65 times faster |
| Ubuntu | Any local disk | |
| Windows | Any local NTFS folder under your user profile | Your profile is private to you by default |

**Time:** estimated at 1 to 2 hours on a typical connection, though this hasn't been timed. Plug the laptop in and keep it awake.

**If it stops** because the connection dropped or the laptop slept, run the same command again. It resumes partial downloads and skips anything already built.

To skip the 4.7 GB Pfam read, at the cost of domain positions:

```bash
varsight setup nopfam
```

At the end it runs the known-variant self-test. You want to see `11 of 11 known variants passed`.

**If a known variant fails or setup stops with an error,** please keep the message and report it. This full run hasn't been done before, so a problem here is useful to know about.

## 3. Annotate

```bash
varsight annotate "BRAF V600E" build=hg19
varsight annotate "NM_004333.6:c.1799T>A" build=hg38
varsight annotate "chr7:g.140453136A>T" build=hg19
varsight annotate variants.txt build=hg19
varsight annotate sample.vcf.gz
```

A VCF's build is read from its header. Results go in a folder named after the input, `sample_varsight/`, holding `report.html`, `results.xlsx` and `results.tsv`.

**Always put a variant in quotes,** because of the `>`.

Paths work as usual on each system:

```bash
varsight annotate /mnt/i/runs/SAMPLE.hard-filtered.vcf.gz out=/home/you/results/SAMPLE
```

```powershell
varsight annotate "I:\runs\SAMPLE.hard-filtered.vcf.gz" "out=C:\Users\you\results\SAMPLE"
```

Keep results with patient data on a private disk: see [SECURITY.md](../SECURITY.md).

## 4. See how it decided

```bash
varsight explain "BRAF V600E" build=hg19
```

This prints every step, file and query. [How lookups work](how-lookups-work.md) explains the output.

## 5. Add advanced databases

```bash
varsight databases
varsight add hotspots
varsight add dbsnp
varsight add revel accept
varsight add cosmic_census file=Cosmic_CancerGeneCensus_v101_GRCh38.tsv.gz accept
varsight remove dbsnp
```

`varsight databases` lists each one with its size, its licence, what is installed, and how to get it. The next report gains columns only for what you added. See [data sources](data-sources.md#advanced-databases).

## 6. Stay up to date

A notice appears by itself when a newer version of a database exists, and it doesn't slow the run. When you see it:

```bash
varsight update check
varsight update
varsight rollback
```

To automate it, see [updates.md](updates.md).

## 7. Measure accuracy on your own data

These two give you the proof for your own installation:

```bash
varsight validate sample=20000 report=varsight-validation.md
varsight benchmark report=varsight-benchmark.md
```

- `validate` tests 20,000 ClinVar names per build across the whole genome, and lists every mismatch group.
- `benchmark` measures speed and memory on your computer.

To compare with an annotated VCF you already trust, such as one from DRAGEN or VEP, using `python3` on Linux and `python` on Windows:

```bash
python3 tools/compare_annotated_vcf.py SAMPLE.annotated.vcf.gz differences=differences.txt
```

It stays on your computer, but treat its output as private.

## 8. Send your changes to GitHub

The repository is at [github.com/prabirsh/VarSight](https://github.com/prabirsh/VarSight). To publish changes, from Ubuntu or WSL:

```bash
git add -A
git status
git commit -m "Describe the change"
git push
```

**Read the `git status` list before committing.** No patient VCF or results file should appear in it. `.gitignore` already leaves out VCFs, report folders, results and data folders, apart from the demo.

On GitHub, go to **Settings**, then **Code security**, and turn on **Private vulnerability reporting**. CI runs by itself on Python 3.10 to 3.14 on Linux, and on Windows.

## Problems

| Message | What to do |
| - | - |
| `varsight: command not found` | Linux: run the PATH line that `install.sh` printed, or use `~/.local/bin/varsight`. Windows: open a new window |
| `running scripts is disabled` on Windows | Use the exact command in step 1, with `-ExecutionPolicy Bypass`. It applies to that one run only |
| `Python was not found` on Windows | That's the Microsoft Store placeholder. Install Python from python.org and tick Add python.exe to PATH |
| `no VarSight data in ...` | Run `varsight setup`, or add `data=FOLDER` |
| `not enough disk space` | Free up space, or choose another folder with `data=`. In WSL a folder under `/mnt/` works but is slower |
| `another setup or update is running` | Wait for it. If nothing is running, run the command again: a leftover lock from a stopped run clears itself |
| `REF ... does not match` | Usually the wrong build. Check `build=` |
| `build= says hg38, but the contig lengths ... are hg19` | Your VCF is on the other build |
| Very slow lookups | In WSL, check that the data isn't under `/mnt/`. `varsight doctor` warns about it |
| Out of disk on C: after removing data in WSL | WSL's virtual disk doesn't shrink on its own. Search for "compact WSL vhdx" |
| `expected columns are missing` | A source changed its format. Nothing was installed and your current data is safe. Please report it |
| Git says `dubious ownership` | Run the `git config` command that the message gives |
