sbatch -J 2mm -o %x_%N.o%j -e %x_%N.e%j -p gpu_edu -N 1 -n 1 -c 16 --gres=gpu:4 --time 00:20:00 --exclude=dumbo0[60,62,64-66,68-71] ./run.sh
