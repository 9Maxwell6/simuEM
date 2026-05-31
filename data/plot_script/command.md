
python3 iterations.py ../data_test/curlcurl_l2.dat ../data_test/curlcurl_l2_pc.dat -o ./curlcurl_iteration.pdf
python3 solve_time.py ../data_test/curlcurl_l2.dat ../data_test/curlcurl_l2_pc.dat -o ./curlcurl_solve_time.pdf

python3 iterations.py ../T_Omega_l2_pc.dat ../T_Omega_l2.dat -o ./T_Omega_iteration.pdf
python3 solve_time.py ../T_Omega_l2_pc.dat ../T_Omega_l2.dat -o ./T_Omega_solve_time.pdf

python3 condition_number.py ../data_test/curlcurl_l2_pc_1e-10.dat ../data_test/curlcurl_l2_1e-10.dat ../T_Omega_l2_pc.dat ../T_Omega_l2_pc_GMRES_3.dat ../T_Omega_l2_GMRES_3.dat 

