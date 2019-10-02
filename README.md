# GATK-BP-bash
# Bash scripts implementing the Best Practices Workflows of the GATK team at the Broad Institute

These scripts have been developed by José R. Valverde and Lorena Magraner at CNB-CSIC

www.cnb.csic.es

These scripts are based on publicly available works:

* The gatk*.sh scripts are based on the published workflows available from the GATK team

* The scripts in NYU are based on published scripts distributed by NYU

## the deal

We'd like to properly publish this work, or at least get acknowledged for it. If you decide
to use these scripts, please do cite them in any publication, either referring to this
GitHub site, or -shuld we succeed in publishing them- the corresponding publication, which
we'll mention here. But, at the very least, we beg you to cite our work:

Valverde, J.R., Magraner, L. (2019) Bash scripts implementing GATK Best Practices 
Workflows. GitHub. https://github.com/jrvalverde/GATK-BP-bash

And, of course, we want to thank you in advance for your kindness.

## GATK-Best-Practices based workflows
 
For the GATK BP workflows, we have based the scripts on the latest available versions
published at the GATK web site. These workflows suffer from a number of "features": they 
make use of a mixture of various versions of GATK and, within them, of various sub-versions.
 
The main issue with GATK workflows is that, at face value, you need to have installed at 
least the corresponding versions (3 and 4) of GATK. Not so evident, even within the same 
version, the command line also has changes from one sub-version to another (e.g. calling the
same program may be different from one release of GATK4 to the next, some options may
disappear or change name, etc..). GATK is a software in constant flux and that shows up.
So, you should also have all the different sub-versions as well, and the same they do,
which are not documented in the workflows (you have to find out).
 
We have made our scripts to work with the latest versions of GATK available at key
points: GATK4-0.12.0 initially, and currently GATK4-1.2.0. This implies that, to run
these scripts you need to use -as of now- GATK4-1.2.0. Later versions will likely work,
but since GATK is in constant flux, it cannot be guaranteed.
 
If and when we need to run these scripts again, we will update our version of GATK to 
the latest one and re-adapt the scripts if needed. So, if you use a more recent version
of GATK, and if it fails to run, you have two options

### Fixing problems
 
* __preferably__: identify the program where the script fails, check current GATK
     documentation, pray it is up to date, and modify the script until it works (note
     that in some instances the GATK documentation may not be up to date, in our
     experience, when this happens, running the program with --help *may* give you
     updated information --though not always)
     
* if you are not in a hurry: wait to see if we get a new job to process and 
     when we update our GATK and our version of the scripts, then grab the updated one
     (this is not recommended as it mat take some time for us to get new assignments)
     
* _yes, I know, we said two, but_... you can also contact us, we cannot promise
     much: we are usually working on several projects and may not have spare time, 
     but if we can, we'll try to fix them
     
 * _a potentially better option_ is to ask us to collaborate with you on the analysis of
     your data: this would provide us with an enticement and a justification to
     continue working on this.


### So, what are these GATK scripts?
  
They are a direct translation of the latest (as of 2019) published GATK workflows to
bash, so that you can run the workflow on your computer without needing to use CROMWELL
which, in our experience, is non-trivial to install and use and adds significant overheads.

The scripts we provide implement __all the analysis steps included in the GATK BP workflows
using the same options__ (except for some oddity, actualization or version incompatibility)
as the GATK team uses in theit BP workflows. The main difference is that we do NOT split
computation at key steps: each input file is processed at once.

### Are they for me?

__BEWARE__: _they are NOT a complete translation_: the original workflows split BAM and VCF files
at several steps into many parts, process each part separately and then collect all the
results and collate (gather) them to reconstruct a single result file. They do this so
that if you want to analyze one file, you have many computers, you have CROMWELL installed
and configured to use all of them (or you pay to use their cloud solution), you can
speed up the processing of a single file.

Despite the overhead, such a speed up may be significant for a single file when you have 
many computers. It is not (and actually incurs in additional overheads) when you have 
more files to analyze than computers (or as many files as computers). So, if you have many, 
many computers (physical or virtual in the cloud), it may be worth learning to install, use 
and configure CROMWELL, specially to process life-critical diagnostic data that is being 
expected by your doctor.

If you have few files to process and you are not in a hurry to analyze them (e.g. it is a 
clinical assay), or if you only have one or a few computers to do the analysis, or if
you have many files to process, then using the CROMWELL workflows may be overkill and
it will be far easier to use our BASH translation of the CROMWELL files. You might save a
lot of work, probably a lot of time, and may find it a lot easier to configure the
analysis. Then, you are welcome to download and use these scripts.

### Other advantages

The original workflows contain provisions for __restarting calculations__. We do too and
in a simpler form: our scripts check if the output exists at each step and if it does, 
they skip the step. So, if anything goes wrong, you can just fix the problem and re-run
the script: it will not repeat any work already finished and will continue from 
wherever it left. 

This has an extra convenience: if you want to __repeat any step__ (say because a file was 
deleted accidentally or was corrupted), you can remove the files corresponding to that 
step, and only those files will be rebuilt. NOTE: only the missing ones, so later steps, 
if already done, will not be recalculated. If you want to recalculate everything after
a given step, you should remove all files after that step as well.

We also make another difference from GATK BP workflows: they use "ad-hoc", inconsistent
names for the files generated. In our scripts, each step adds a new suffix to the name
of the file generated after the previous step. This makes it easier to __follow the
progress__ of the calculation (each step will use a longer output file name) and to
control it: if you want to __repeat everything after a given step__, you can use 
__rm file.step1.step2..stepn.*__ and all files generated after stepn will be removed. Now,
simply re-run the script and those missing steps will be completed.


## NYU based best-practices workflows

As for the NYU-based scripts: we have only added additional checks to the scripts, so that
a step is only carried out if needed, i.e. if the output file is missing. Much as in the
previous description for GATK-based scripts, this allows for restarting interrupted 
jobs or for recalcularing specific steps. We do also ensure all file names follow
an add-on convention as in the former scripts: each step adds a new suffix, so that
if you want to recalculate everything after a given step you can easily remove all
subsequent files and re-run the script.

We have also updated GATK options to use the calling/command convention of the latest
GATK release (which for us has been GATK4-1.2.0 as indicated above), removing all
calling inconsistencies. If you use the same (and hopefully a later) version, you
should have no problem running them.

## Other scripts:

In `aux` you will find odds & ends of scripts that we have used at one point or
another to tweak our data. They may provide useful hints for some problems, but
they are not directly usable in general and will require detailed tweaking. They
are provided for what they're worth. We may consider at some point making them 
more general, but it is not high on our priority list.


We are working on implementing CNV measures. We are trying with a number of programs
like XHMM, CNVkit and others. We'll be adding scripts as we develop them. Same 
applies for annotation, using snpEff, ANNOVAR or other tools. Stay tuned.


## The important stuff: how do I run them?

Easy, run them. The first thing they'll do is look for a configuration file. If it
does not exist, they will copy a sample one into your folder and then you can look
at its contents and modify them to configure your environment.

The configuration file is shell script that defines several variables. All you
have to do is modify their value after the '=' sign.

Then, re-run the script and let it run...

Well, better not! If you do, it's OK, but if there is a problem, it may be
difficult to identify (GATK is too verbose). So, what we advice is that (yes,
it is OK to simply run them) instead, you pipe all output to a file so that
when the script ends, you can check if there were errors. For instance:

Runing a script (e.g.):

> bash gatk-BP2019-processing-for-variant-discovery.sh

Saving output blindly (e.g.):

> bash gatk-BP2019-processing-for-variant-discovery.sh > LOG 2>&1

Saving output and seeing it (what we prefer):

> bash gatk-BP2019-processing-for-variant-discovery.sh |& tee LOG


## Easter eggs:

In some scripts we have installed an easter egg: to make them more useful,
we detect whether we have been called with one, many or no arguments and
do the (hopefully) right thing. This allows using the same script for 
processing one file, many files or doing the default (all files in a
predefined subdirectory).

