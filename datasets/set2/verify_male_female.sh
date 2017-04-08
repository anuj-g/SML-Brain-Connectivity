for f in `grep " female" IACL-blsa-subjectdata-1110.csv| cut -d"," -f1` ; do ls females/"$f"_fbrCountConnMtx.csv; done | wc -l

for f in `grep " male" IACL-blsa-subjectdata-1110.csv| cut -d"," -f1` ; do ls males/"$f"_fbrCountConnMtx.csv; done | wc -l

