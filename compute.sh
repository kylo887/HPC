#!/bin/bash
#To give the execute permission on "run.bash"
#<chmod +x run.bash> OR <chmod 0755 run.bash> OR <chmod 777 run.bash>
#To execute "run.bash" type <./run.bash>

# echo Strong Scaling Particle simulation is running...
# for px in {1..12}
# do
#     echo "Processor count :: $px"
#     echo "Processor count :: $px" >> sim.log
#     for dt in 0.01 0.02 0.05 0.1 0.2
#     do 
#         echo "Time step :: $dt"
#         echo "Time step :: $dt" >> sim.log
#         mpirun -np $px ./a.out 200 60 $dt >> sim.log
#         echo "------------------------" >> sim.log
#     done
# done


echo Weak Scaling Particle simulation is running...
for px in {1..12}
do
    echo "Processor count :: $px"
    echo "Processor count :: $px" >> sim.log
    for dt in 0.01 0.02 0.05 0.1 0.2
    do 
        echo "Time step :: $dt"
        echo "Time step :: $dt" >> sim.log
        NP=$px*100
        mpirun -np $px ./a.out $NP 60 $dt >> sim.log
        echo "------------------------" >> sim.log
    done
done