 #!/bin/bash
mpicc main.c
        for i in `seq 1 8`;
        do
        echo num threads: $i;
               mpiexec -n $i a.out;
                #echo '\n';
        done     
