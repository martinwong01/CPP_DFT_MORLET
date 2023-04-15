
#!/bin/bash


for((i=0;i<1000;i++))
do
(time -p ./fft_bench0.exe) >> result0.txt 2>&1
(time -p ./fft_bench1.exe) >> result1.txt 2>&1
done
