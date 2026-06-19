
python3 iterations.py ../data_test/curlcurl_l2.dat ../data_test/curlcurl_l2_pc.dat -o ./curlcurl_iteration.pdf
python3 solve_time.py ../data_test/curlcurl_l2.dat ../data_test/curlcurl_l2_pc.dat -o ./curlcurl_solve_time.pdf

python3 iterations.py ../T_Omega_l2_pc.dat ../T_Omega_l2.dat -o ./T_Omega_iteration.pdf
python3 solve_time.py ../T_Omega_l2_pc.dat ../T_Omega_l2.dat -o ./T_Omega_solve_time.pdf

python3 condition_number.py ../data_test/curlcurl_l2_pc_1e-10.dat ../data_test/curlcurl_l2_1e-10.dat ../T_Omega_l2_pc.dat ../T_Omega_l2_pc_GMRES_3.dat ../T_Omega_l2_GMRES_3.dat 

python3 iterations.py ../T_Omega_l2_pc_CG_1.dat ../T_Omega_l2_pc_CG_100.dat ../T_Omega_l2_CG.dat ../T_Omega_l2_pc_CG_1_global.dat ../T_Omega_l2_pc_CG_100_global.dat ../data_test/curlcurl_l2_1e-10.dat ../data_test/curlcurl_l2_pc_1e-10.dat 



python3 ./plot_script/iterations.py ./T_Omega_CG.dat ./T_Omega_pc_CG_1_d.dat ./T_Omega_pc_CG_1_g.dat ./T_Omega_pc_CG_1_c.dat ./T_Omega_pc_CG_1_f.dat -o ./T_Omega_iterations_CG_1.pdf


python3 ./plot_script/condition_number.py ./T_Omega_CG.dat ./T_Omega_pc_CG_1_d.dat ./T_Omega_pc_CG_1_g.dat ./T_Omega_pc_CG_1_c.dat ./T_Omega_pc_CG_1_f.dat -o ./T_Omega_condition_CG_1.pdf -n=4


python3 ./plot_script/condition_number.py  ./T_Omega_pc_CG_1_d.dat ./T_Omega_pc_CG_2_d.dat ./T_Omega_pc_CG_5_d.dat ./T_Omega_pc_CG_10_d.dat ./T_Omega_pc_CG_15_d.dat ./T_Omega_pc_CG_20_d.dat --no-fit -o ./T_Omega_condition_CG_pc_1.pdf 


python3 ./plot_script/iterations.py  ./T_Omega_pc_CG_1_d.dat ./T_Omega_pc_CG_2_d.dat ./T_Omega_pc_CG_5_d.dat ./T_Omega_pc_CG_10_d.dat ./T_Omega_pc_CG_15_d.dat ./T_Omega_pc_CG_20_d.dat  --no-fit -o ./T_Omega_iterations_CG_pc_1.pdf 


python3 ./plot_script/condition_number.py  ./T_Omega_pc_CG_1_g.dat ./T_Omega_pc_CG_2_g.dat ./T_Omega_pc_CG_5_g.dat ./T_Omega_pc_CG_10_g.dat ./T_Omega_pc_CG_15_g.dat ./T_Omega_pc_CG_20_g.dat  --no-fit -o ./T_Omega_condition_CG_pc_2.pdf 


python3 ./plot_script/iterations.py  ./T_Omega_pc_CG_1_g.dat ./T_Omega_pc_CG_2_g.dat ./T_Omega_pc_CG_5_g.dat ./T_Omega_pc_CG_10_g.dat ./T_Omega_pc_CG_15_g.dat ./T_Omega_pc_CG_20_g.dat  --no-fit -o ./T_Omega_iterations_CG_pc_2.pdf 


python3 ./plot_script/condition_number.py  ./T_Omega_pc_CG_1_c.dat ./T_Omega_pc_CG_2_c.dat ./T_Omega_pc_CG_5_c.dat ./T_Omega_pc_CG_10_c.dat ./T_Omega_pc_CG_15_c.dat ./T_Omega_pc_CG_20_c.dat  --no-fit -o ./T_Omega_condition_CG_pc_3.pdf 


python3 ./plot_script/iterations.py  ./T_Omega_pc_CG_1_c.dat ./T_Omega_pc_CG_2_c.dat ./T_Omega_pc_CG_5_c.dat ./T_Omega_pc_CG_10_c.dat ./T_Omega_pc_CG_15_c.dat ./T_Omega_pc_CG_20_c.dat  --no-fit -o ./T_Omega_iterations_CG_pc_3.pdf 


python3 ./plot_script/condition_number.py  ./T_Omega_pc_CG_1_f.dat ./T_Omega_pc_CG_2_f.dat ./T_Omega_pc_CG_5_f.dat ./T_Omega_pc_CG_10_f.dat ./T_Omega_pc_CG_15_f.dat ./T_Omega_pc_CG_20_f.dat  --no-fit -o ./T_Omega_condition_CG_pc_4.pdf 


python3 ./plot_script/iterations.py  ./T_Omega_pc_CG_1_f.dat ./T_Omega_pc_CG_2_f.dat ./T_Omega_pc_CG_5_f.dat ./T_Omega_pc_CG_10_f.dat ./T_Omega_pc_CG_15_f.dat ./T_Omega_pc_CG_20_f.dat  --no-fit -o ./T_Omega_iterations_CG_pc_4.pdf 
