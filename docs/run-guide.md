# Running VarShift on Ubuntu, WSL and Windows

VarShift runs the same way on all of these. Only the install step differs. Pick your system in step 1, and the rest of the guide is the same, apart from the notes marked for Windows.

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
git clone https://github.com/prabirsh/VarShift.git ~/VarShift
cd ~/VarShift
sh install.sh
```

If it tells you to add the command to your PATH, run the line it prints.

### Windows PowerShell

Install Python from [python.org](https://www.python.org/downloads/) and tick **Add python.exe to PATH**, or in PowerShell run `winget install Python.Python.3.12` and open a new window. Then download the repository as a ZIP from GitHub and unpack it, or use `git clone`. In PowerShell, from the VarShift folder:

```powershell
powershell -ExecutionPolicy Bypass -File .\install.ps1
```

### Windows Command Prompt

The same Python, then from the VarShift folder:

```text
install.bat
```

On Windows it asks before adding the command to your PATH. **Open a new window afterwards** so the change takes effect.

### What the install does

It writes one small `varshift` launcher, runs the [demo](../demo/README.md) with no download, and offers to set up the databases. Check it:

```bash
varshift version
varshift doctor
```

Open `demo_output/variants/report.html` in a browser to see what a report looks like. It uses real data for 12 genes.

| System | To open the demo report |
| - | - |
| WSL | `explorer.exe demo_output` |
| Ubuntu desktop | `xdg-open demo_output/variants/report.html` |
| Windows | Double-click it, or `start demo_output\variants\report.html` |

### Installer options

Both installers take the same choices. Running one again is safe: it rewrites the launcher.

| To | Ubuntu or WSL (`sh install.sh ...`) | Windows (`install.ps1 ...`) |
| - | - | - |
| Put the command in another folder | `bin=DIR` (default `~/.local/bin`) | `-BinDir DIR` (default `%LOCALAPPDATA%\VarShift\bin`) |
| Skip the demo | `demo=no` | `-NoDemo` |
| Not offer setup at the end | `setup=no` | `-NoSetup` |
| Leave PATH alone | It never changes PATH; it prints the line to add | `-NoPath` |
| Run the unit tests too | `test=yes` | `-Test` |

`install.bat` passes the same `-` options on to `install.ps1`. For example, to keep everything inside one project folder and add nothing to PATH:

```bash
sh install.sh bin=/mnt/p/Projects/VarShift-bin demo=no setup=no
```

```powershell
powershell -ExecutionPolicy Bypass -File .\install.ps1 -BinDir P:\Projects\bin -NoPath -NoSetup
```

## 2. Set up the databases

```bash
varshift setup
```

It asks three things, shows the size, and then waits for your yes before downloading:

```text
Which genome build do you work with?
  1  hg19 (GRCh37)
  2  hg38 (GRCh38)
  3  Both. Results then always show the position on both builds
Your choice, 1, 2 or 3: [3]

Where should the data go? [/home/you/varshift-data]

Advanced databases are optional, and can be added any time later.
   1  dbSNP                              26 GB per build, read without saving
   ...
Numbers to add now, separated by commas, or press Enter for none:

This will install: hg19 into /home/you/varshift-data
  Download   about 1.7 GB
  Disk       about 5 GB kept, and keep about 15 GB free during setup
  Time       about 45 minutes on a 2 MB/s connection. If it stops, run it again and it carries on
Start? yes or no: [yes]
```

The same thing without questions, for scripts:

```bash
varshift setup build=hg38 data=/data/varshift yes
```

**Where to put the data:**

| System | Use | Why |
| - | - | - |
| WSL | `~/varshift-data`, on the Linux disk | Lookups measured 27 microseconds there, against 1.3 to 1.8 ms on a Windows drive seen from WSL (`/mnt/c`, `/mnt/i`) |
| Ubuntu | Any local disk | |
| Windows | Any local NTFS folder under your user profile | Your profile is private to you by default |

**If your Linux disk has no room, keep the data on a Windows drive, but build it from PowerShell.** Building from WSL on a Windows drive is about 30 times slower and, with WSL's memory tight, can stop with `Cannot allocate memory`, so `setup`, `update` and `add` refuse it unless you add the word `windowsdrive`. WSL's Linux disk is a file on C:, so `~/varshift-data` uses C: too. Set up and update in PowerShell, then point WSL at the same files:

```powershell
varshift setup build=hg19 data=I:\VIP-Annotations yes
```

```bash
varshift config data=/mnt/i/VIP-Annotations
```

Each `varshift` command then takes about 2 seconds in WSL. See [the measurements](benchmarks.md#data-on-a-windows-drive-used-from-wsl). `varshift update check` works from WSL, as it changes nothing. WSL refuses the Windows spelling `I:\VIP-Annotations`, because bash would remove the backslash.

**Which data folder is used:** `data=` on the command, else the `VARSHIFT_DATA` variable if it is set, else the folder saved by setup or `varshift config data=PATH`, else `~/varshift-data`. `varshift config` shows the folder in use, and says when `VARSHIFT_DATA` or the default chose it. If `VARSHIFT_DATA` is set, a saved folder is ignored, and `varshift config` says that too.

**Time:** about 45 minutes for one build on a 2-core laptop with a 2 MB/s connection, and an estimated 1.5 hours for both. See [the time for each step](benchmarks.md#setup-time). Plug the laptop in and keep it awake.

**While it runs,** it always says what it is doing. Every long download prints how far it has got, its speed and the time left, and every long build says once a minute that it is still going. So more than a few minutes of silence means something is wrong.

```text
[5/6] clinvar.hg19
  downloading ClinVar GRCh37 VCF (191.4 MB)
    clinvar.vcf.gz: 40% of 191 MB, 2.4 MB/s, about 48 s left
  building clinvar.hg19
    still building clinvar.hg19, 1 min so far
  finished clinvar.hg19 in 20 min
```

**If it stops,** or you stop it with Ctrl+C, run the same command again. It resumes partial downloads and skips anything already built. If the data folder isn't the one it offers, give it with `data=`.

At the end it runs the known-variant self-test. You want all of them to pass: `5 of 5 known variants passed` on hg19 alone.

**If a known variant fails or setup stops with an error,** please keep the message and report it. Few complete runs have been done so far, so a problem here is useful to know about.

## 3. The guided start

Type `varshift` on its own in a terminal, or `varshift start`:

```text
VarShift: offline variant annotation. Data: hg19, ClinVar 2026-11-01

What do you have? Type one or more, separated by ;  (or paste several lines)
  1  Protein change     BRAF V600E
  2  cDNA change        NM_004333.6:c.1799T>A
  3  Genomic position   chr7:g.140453136A>T   or   chr7-140453136-A-T
  4  rs number          rs113488022
  5  A file             variants.txt   or   sample.vcf.gz
Type a number for examples of that kind, help for what else you can type, or Enter to quit.

>
```

Type a variant the way you have it. One gets its full answer, and several get a numbered line each. A slip is asked about, and when nothing decides between candidates you choose ([how](input-output.md#when-the-text-is-a-little-off)). The data stays open, so each answer after the first is immediate. Nothing about your queries is written to disk unless you save.

| Type | To |
| - | - |
| `save` | Save the answers to Excel: `varshift-session-DATE-TIME.xlsx`, in the current folder |
| `save NAME.xlsx` | Save to that file. `.html`, `.tsv`, `.json` and `.jsonl` write those formats |
| `save FOLDER/` | Write `report.html`, `results.xlsx` and `results.tsv` in that folder |
| `list` | See the answers so far, numbered |
| A number | See that answer in full again |
| `all N` | See every candidate of answer N in full, not only the best |
| `remove N`, `clear` | Drop one answer, or all of them |
| `mode tumour`, `mode inherited`, `mode off` | Change which evidence and links come first, for the queries after it |
| `examples`, `help` | How to write a variant, and this list |
| A file name | Annotate that file into its own results folder, as `varshift annotate` would |
| Enter, or `quit` | Quit. If any answer isn't saved, it asks first |

**What is saved:** the answers you saw, with every link in the workbook's Links sheet. For each one that's the best candidate, with the others named in `also`, or the candidates you chose. A second `save` in the same minute adds `-2`, so it never overwrites. Saving to a file that exists asks before replacing it.

**No data yet?** It offers the 12-gene demo data at once, with no download, or setup. The demo is built in a temporary folder in the current directory and removed when you quit, and anything saved from it says it was made with the demo data.

Ctrl+C during a query stops only that query. At the prompt it stops the session and says how many answers weren't saved. To set the build for positions and the mode from the start:

```bash
varshift start build=hg19 mode=tumour
```

Piped or scheduled, `varshift` alone prints the command list, and `varshift start` says it needs a terminal. Neither ever waits for an answer.

## 4. Annotate

```bash
varshift annotate "BRAF V600E" build=hg19
varshift annotate "NM_004333.6:c.1799T>A" build=hg38
varshift annotate "chr7:g.140453136A>T" build=hg19
varshift annotate "BRAF V600E; KRAS G12D; TP53 R175H" build=hg19
varshift annotate variants.txt build=hg19
varshift annotate sample.vcf.gz
```

One variant gets one answer: the best candidate's card, with why it was chosen and the others named, or a choice when nothing decides between them. These change what the card shows:

```bash
varshift annotate "BRAF V600E" build=hg19 show=all
varshift annotate "KRAS G12D" build=hg19 mode=tumour links=all
varshift annotate "BRAF V600E" build=hg19 color=no
```

`show=all` prints every candidate in full, and `links=all` every link with its label and note. `mode=tumour` puts CIViC, OncoKB and cBioPortal first, and `mode=inherited` puts ClinGen, gnomAD constraint and OMIM first. `varshift config mode=tumour` saves a mode, and `varshift config color=no` saves the colour choice. Several variants separated by `;` get one numbered line each, or all go in the file that `out=` names. `varshift examples` shows every way to write a variant.

A VCF's build is read from its header. Results go in a folder named after the input, `sample_varshift/`, in the current directory, holding `report.html`, `results.xlsx` and `results.tsv`. To put them somewhere else, add `out=FOLDER`, or `out=FILE.xlsx` for one file.

**Always put a variant in quotes,** because of the `>`.

Paths work as usual on each system:

```bash
varshift annotate /mnt/i/runs/SAMPLE.hard-filtered.vcf.gz out=/home/you/results/SAMPLE
```

```powershell
varshift annotate "I:\runs\SAMPLE.hard-filtered.vcf.gz" "out=C:\Users\you\results\SAMPLE"
```

Keep results with patient data on a private disk: see [SECURITY.md](../SECURITY.md).

## 5. See how it decided

```bash
varshift explain "BRAF V600E" build=hg19
```

This prints every step, file and query. [How lookups work](how-lookups-work.md) explains the output.

## 6. Add advanced databases

```bash
varshift databases
varshift add hotspots
varshift add dbsnp
varshift add revel accept
varshift add cosmic_census file=Cosmic_CancerGeneCensus_v101_GRCh38.tsv.gz accept
varshift remove dbsnp
```

`varshift databases` lists each one with its size, its licence, what is installed, and how to get it. The next report gains columns only for what you added. See [data sources](data-sources.md#advanced-databases).

## 7. Stay up to date

A notice appears by itself when a newer version of a database exists, and it doesn't slow the run. When you see it:

```bash
varshift update check
varshift update
varshift rollback
```

To automate it, see [updates.md](updates.md).

## 8. Measure accuracy on your own data

These two give you the proof for your own installation:

```bash
varshift validate sample=20000 report=varshift-validation.md
varshift benchmark report=varshift-benchmark.md
```

- `validate` tests 20,000 ClinVar names per build across the whole genome, and lists every mismatch group.
- `benchmark` measures speed and memory on your computer.

To compare with an annotated VCF you already trust, such as one from DRAGEN or VEP, using `python3` on Linux and `python` on Windows:

```bash
python3 tools/compare_annotated_vcf.py SAMPLE.annotated.vcf.gz differences=differences.txt
```

It stays on your computer, but treat its output as private.

## 9. Send your changes to GitHub

The repository is at [github.com/prabirsh/VarShift](https://github.com/prabirsh/VarShift). To publish changes, from Ubuntu or WSL:

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
| `varshift: command not found` | Linux: run the PATH line that `install.sh` printed, or use `~/.local/bin/varshift` (or your `bin=` folder). Windows: open a new window |
| `running scripts is disabled` on Windows | Use the exact command in step 1, with `-ExecutionPolicy Bypass`. It applies to that one run only |
| `Python was not found` on Windows | That's the Microsoft Store placeholder. Install Python from python.org and tick Add python.exe to PATH, or run `winget install Python.Python.3.12`, then open a new window |
| `... is a Windows drive seen from WSL` | Run that `setup`, `update` or `add` from PowerShell, with the folder's Windows path. WSL then reads the same folder. `varshift update check` works in WSL |
| `no VarShift data in ...` | Run `varshift setup`, or add `data=FOLDER` |
| `not enough disk space` | Free up space, or choose another folder with `data=`. In WSL a folder under `/mnt/` works but is slower |
| `another setup or update is running` | Wait for it. If nothing is running, run the command again: a leftover lock from a stopped run clears itself |
| `Did you mean: ...` | VarShift couldn't read the description as written, and these are what it probably meant. Run the one you meant. At a terminal, with no `out=` or `format=`, it would have asked |
| `Not sure of the format?` | Nothing near it could be found. `varshift examples` shows every form that works |
| `varshift start needs a terminal` | It was run from a script, a pipe or a scheduled task. Use `varshift annotate` there |
| Odd characters such as `←[1;31m` around words | That console doesn't understand colour. Add `color=no`, or keep it with `varshift config color=no` |
| `the mitochondrial rCRS sequence (NC_012920.1) comes with the hg38 reference` | On hg19, mitochondrial names need hg38 set up too: `varshift setup build=hg38` |
| `REF ... does not match` | Usually the wrong build. Check `build=` |
| `build= says hg38, but the contig lengths ... are hg19` | Your VCF is on the other build |
| Very slow lookups | In WSL, check that the data isn't under `/mnt/`. `varshift doctor` warns about it |
| Out of disk on C: after removing data in WSL | WSL's virtual disk doesn't shrink on its own. Search for "compact WSL vhdx" |
| `expected columns are missing` | A source changed its format. Nothing was installed and your current data is safe. Please report it |
| Git says `dubious ownership` | Run the `git config` command that the message gives |
