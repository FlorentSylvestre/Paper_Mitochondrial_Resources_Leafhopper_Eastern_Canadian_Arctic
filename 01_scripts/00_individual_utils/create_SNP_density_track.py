#!/bin/env python
import sys
from collections import defaultdict

inf, outf, wsize = sys.argv[1:4]
wsize = int(wsize)

snps_dens = defaultdict(lambda:defaultdict(int))
mpileup = open(inf, "r")
header = "Start\tScore\n"


for x, line in enumerate(mpileup):
    if line.startswith("#"):
        continue
    l = line.strip().split()
  
    if l[4] == "<*>":
        continue

    DP = list(map(int, l[9].split(":")[1].split(",")))

    if sum(DP[1:])/sum(DP) < 0.1:
        continue

    snps_dens[l[0]][int(l[1]) // wsize *wsize]  +=1
mpileup.close()

for i in snps_dens:
    with open(outf + "/" + i + ".SNPdensity", "w") as out:
            res = ""
            out.write(header)
            for x in snps_dens[i]:
                res += f"{str(x)}\t{str(snps_dens[i][x])}\n"
            out.write(res)


