#ifndef FUNC_STATE_H
#define	FUNC_STATE_H
#include "global_variables.h"

//read instance file
void read_instance(char *instance_name);	

//read case study
void read_case(char *instance_name);

//allocate memory for some global variables
void allocate_memory();

//compute the objective of a solution 
double compute_obj(solution_data sol);

//generate a feasible initial solution in a greedy manner
void greedy_construct_initialsol(solution_data &sol);

//generate a feasible initial solution in a random manner (for analysis)
void random_construct_initialsol(solution_data &sol);

//copy solution_data structure
void copy_solution(solution_data source, solution_data &des);

//initialize the global variables: Min_dis, Sec_dis, Num_min_dis, Vec_min_dis
//time complexity: O(n * |M|), |M| is the number of vertices selected to current solution
void initialize_sup_arrays(solution_data sol);

//initial three hash vectors
void initial_hash();

//compute move gain for f(s) for each candidate move (one flip)
double compute_flip_move_gain(int node, solution_data sol);

//compute move gain for f(s) for each candidate move (two swap)
double compute_swap_move_gain(int node_x, int node_y, solution_data sol);

//void quick_sort(int *sort_array, int *node_array, int low, int high);

//verify solution
void verify_sol(solution_data sol);

//output the found best solution of each run
void out_sol(solution_data sol, char *ins_name, char *out_sol_file);

//output the found best results of each run
void out_stat_results(char *ins_name, char *out_stat_file, double best_result, double run_time);

//output the found results for all run
void out_total_results(char *ins_name, char *out_total_file, double *cost_total, double *time_total);

//free memory
void free_memory();

/******************** local search ******************/
//solution based tabu search: Add, Drop, and Swap move operators with hash function
void solution_based_tabu_search(solution_data &sol);

//traditional attribute based tabu search (for analysis)
void attribute_based_tabu_search(solution_data &sol);


#endif
