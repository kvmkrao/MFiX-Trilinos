gprof  -p -b /home/vkotteda/penta_tpetra/septa_gprof.exe gmon.out > profile.txt
gprof   -l   /home/vkotteda/penta_tpetra/septa_gprof.exe gmon.out > profile_line.txt
gprof   -A   /home/vkotteda/penta_tpetra/septa_gprof.exe gmon.out > profile_anotated.txt

