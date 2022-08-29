

#for i in 10 {12..14}
#for i in 10 13
for i in {22..26}
do
    export X=$i
    printf "Launching: %d\n" $X
    sbatch --job-name="EPMC"$X run_epmcX.slurm
done
