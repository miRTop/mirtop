
SCRIPT="$1"

function run {
    B="http://mirgenedb.org:81/static/data/$1-all.bed"
    F="http://mirgenedb.org:81/static/data/$1-$2-pri-30-30.fas"
    python $SCRIPT/prepare.py --bed ${B} --precursor30 ${F}
}

echo HSA
run hsa hg38
echo MMU
run mmu mm10
# run rno rn6
# run cpo cavPor3
# run ocu oryCun2
# run dno dasNov3
# run gga galGal4
# run dre danRer10
