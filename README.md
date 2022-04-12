# CBB752 Homework Repostiory
[Repo][repo] for extra credit

[repo]: https://github.com/dsm-72/cbb752

## Homework 1
Homework 1 can be found under `/hw1/`. Therein are the following:
- the provided supplementary files, renamed under the directory `/hw1/hw1_supp_files/`. 
- The original input / output files are located at `/hw1/hw1_supp_files/original`.
- The script using the original template is at `/hw1/hw1.py`

Note, when running the script with original output we get identical score matricies for the output files corresponding to `input1` and `input2`. However, on `input2` I seem to get a slightly different alignment.


```bash
# change to supp file dir where input is
$ cd hw1/hw1_supp_files/

# run script (note hw1.py is the same as hw1_template.py, just renamed following PDF instructions)
$ python3 hw1_template.py -i sample-input1.txt -s blosum62.txt
$ python3 hw1_template.py -i sample-input2.txt -s blosum62.txt

# compare output for sample 1
$ diff -E -b sample-output1.txt original/sample-output1.txt

# compare output for sample 2
$ diff -E -b sample-output2.txt original/sample-output2.txt   
119c119
<    K(LDSLSTMEDADV-F-GKA--RTFIE-APL--KNTQPIAGHP-WVK--ILD--AN--I-Y-R-C---RKFSLFESGFYGEVRAVEGHQAQSSNL--R--TV-M---EK) 
---
>    K(LDSLSTMEDADV-F-GKA--RTFIE-APL--KNTQPIAGHP-WVK--ILD--ANIY----R-C---RKFSLFESGFYGEVRAVEGHQAQSSNL--R--TV-M---EK) 
121c121
< IAET(LDKISGYNDA-LRFTGKARYRKYMERIDLSVSRAAPVEGHENRLKNGTLNGVASIIIALARPCAWLRKFSL------GE-EQ-EGH--FSSSLQTRVATVKMPDAEQ)G
\ No newline at end of file
---
> IAET(LDKISGYNDA-LRFTGKARYRKYMERIDLSVSRAAPVEGHENRLKNGTLNGVASIIIALARPCAWLRKFSL------GEEQ--EGH--FSSSLQTRVATVKMPDAEQ)G
```

Notice, that the difference in alginment seems to stem from two characters hopping places:
```bash
$ diff -E -b sample-output2.txt original/sample-output2.txt   
```
Emphasis added for snippet of sequence one:

ILD--AN--***I-Y***-R-C

ILD--AN***IY***----R-C

Emphasis added for snippet of sequence two:

------GE-***EQ***-EGH--

------GE***EQ-***-EGH--

Other than this score matrix and scores are identical. Additionally, provided file to compare to is missing a newline character.
