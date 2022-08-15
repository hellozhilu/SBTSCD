//============================================================================
// Name        : sbts_cdp.cpp
// Author      : Zhi Lu
// Version     : 17 June, 2021
// Copyright   : Business School, University of Shanghai for Science and Technology (USST), Shanghai, China
// Description : Solution-based Tabu Search for the Capacitated Dispersion Problem
//============================================================================

#include "func_state.h"
#include "global_variables.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

using namespace std;


int main(int argc, char *argv[])
{
	//char instance_name[] = "./CDP_Instances/CDP_Instances_b02/GKD-b_11_n50_b02_m5.txt";
	char out_sol_file[] = "./rec/out_sol_test.txt";
	char out_stat_file[] = "./rec/out_stat_test.txt";
	char out_total_file[] = "./rec/out_total_test.txt";

	if (argc < 6)
	{
		cout << "cdp_ts.exe usage: instance_file Fprev A B C";
		exit(-1);
	}
	
	char *instance_name = argv[1];
	Fprev = atof(argv[2]);
	A = atoi(argv[3]);
	B = atoi(argv[4]);
	C = atoi(argv[5]);
	srand(unsigned(time(NULL)));

	Time_limit = 10.0; //10s, 300s
//	if (fabs(Fprev - 999999.0) <= PRECISION)
//		Time_limit = 300.0; //300s
	Runs = 40;
	L = pow(10,8);

	read_instance(instance_name);
//	read_case(instance_name);
	allocate_memory();
	initial_hash();

	int run_cnt = 0;
	double cost_total[Runs], time_total[Runs];
	while (run_cnt < Runs)
	{
		cout << "run_cnt=" << run_cnt << endl;
		solution_data sol_cur;
		sol_cur.ss = new int[Ins.num_v];
		Sol_best.cost = -MAXVALUE;
		Start_time = clock();

		//greedy initial solution
		greedy_construct_initialsol(sol_cur);
//		random_construct_initialsol(sol_cur);
		cout << "sol_cur.cost=" << sol_cur.cost << ", sol_cur.num_sel=" << sol_cur.num_sel << endl;

		//solution based tabu search
		solution_based_tabu_search(sol_cur);
//		attribute_based_tabu_search(sol_cur);
		cout << "run_cnt=" << run_cnt << ", Sol_best.cost=" << Sol_best.cost << ", Run_time=" << Run_time << endl << endl;

		//record information of each run
		out_sol(Sol_best, instance_name, out_sol_file);
		out_stat_results(instance_name, out_stat_file, Sol_best.cost, Run_time);
		cost_total[run_cnt] = Sol_best.cost;
		time_total[run_cnt] = Run_time;

		delete[] sol_cur.ss;
		run_cnt++;
	}

	out_total_results(instance_name, out_total_file, cost_total, time_total);
	free_memory();
	return 0;
}
