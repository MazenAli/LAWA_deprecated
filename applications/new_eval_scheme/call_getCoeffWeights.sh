for f in Debugging_Run2/u_it_*.txt
do
    ./getCoeffWeigths.py < $f > ${f%.txt}_weights.txt
done

for f in Debugging_Run2/r_it_*.txt
do
    ./getCoeffWeigths.py < $f > ${f%.txt}_weights.txt
done

for f in Debugging_Run2/Ap_it_*.txt
do
    ./getCoeffWeigths.py < $f > ${f%.txt}_weights.txt
done

for f in Debugging_Run2/s_it_*.txt
do
    ./getCoeffWeigths.py < $f > ${f%.txt}_weights.txt
done