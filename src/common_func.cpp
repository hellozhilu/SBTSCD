#include "func_state.h"
#include "global_variables.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

using namespace std;


instance_data Ins;
solution_data Sol_best;

double Fprev;
double *Min_dis;
double *Sec_dis;
int *Num_min_dis;	
int *Vec_min_dis;
double Max_min_dis;
double Sec_min_dis;

int *Hash1, *Hash2, *Hash3;
int *W1, *W2, *W3;
int A, B, C;
long L;

int Runs;
double Start_time, Run_time, Time_limit;


//read instance file
void read_instance(char *instance_name)
{
	ifstream fin;
	fin.open(instance_name);
	if (fin.fail())
	{
		cout << "### error open, instance_name " << instance_name << endl;
		exit(-1);
	}

	fin >> Ins.num_v;
	fin >> Ins.capacity;
	Ins.weight = new double[Ins.num_v];

	int count_v = 0;
	while (count_v < Ins.num_v)
		fin >> Ins.weight[count_v++];
	Ins.distance = new double*[Ins.num_v];
	for (int i = 0; i < Ins.num_v; i++)
		Ins.distance[i] = new double[Ins.num_v];

	int line_count = 0;
	while (!fin.eof() && line_count < Ins.num_v)
	{
		count_v = 0;
		while (count_v < Ins.num_v)		
			fin >> Ins.distance[line_count][count_v++];
		line_count++;		
	}
	if (line_count != Ins.num_v)
	{
		cout << "### Error line_count != Ins.num_v, Ins.num_v=" << Ins.num_v
			<< ", line_count=" << line_count << endl;
		exit(-1);
	}
	cout << "finished read, running CDP" << endl;
	cout << "Instance_name=" << instance_name << ", Num_v=" << Ins.num_v
		 << ", capacity=" << Ins.capacity << endl << endl;

	fin.close();

#ifdef DEBUG_READ_INSTANCE
	cout << "num_v=" << Ins.num_v << endl;
	cout << "capacity=" << Ins.capacity << endl;
	cout << "weight=" << endl;
	for (int i = 0; i < Ins.num_v; i++)
		cout << Ins.weight[i] << " ";
	cout << endl << "distance=" << endl;;
	for (int i = 0; i < Ins.num_v; i++)
	{
		for (int j = 0; j < Ins.num_v; j++)
			cout << Ins.distance[i][j] << " ";
		cout << endl;
	}
#endif
}

//read case study
void read_case(char *instance_name)
{
	ifstream fin;
	fin.open(instance_name);
	if (fin.fail())
	{
		cout << "### error open, instance_name " << instance_name << endl;
		exit(-1);
	}

	int nb_edg;
	fin >> Ins.num_v >> nb_edg;

	Ins.distance = new double*[Ins.num_v];
	for (int i = 0; i < Ins.num_v; i++)
		Ins.distance[i] = new double[Ins.num_v];
	for (int i = 0; i < Ins.num_v; i++)
		for (int j = 0; j < Ins.num_v; j++)
			Ins.distance[i][j] = 0.0;

	int max_edg = 0;
	while (!fin.eof() && max_edg < nb_edg)
	{
		int x1, x2;
		double dist;
		fin >> x1 >> x2 >> dist;
		x1--; x2--;
		if (x1 < 0 || x2 < 0 || x1 >= Ins.num_v || x2 >= Ins.num_v)
		{
			cout << "### Error of node : x1=" << x1 << ", x2=" << x2 << endl;
			exit(-1);
		}
		Ins.distance[x1][x2] = Ins.distance[x2][x1] = dist;
		max_edg++;
	}

	fin >> Ins.capacity;
	Ins.weight = new double[Ins.num_v];
	int count_v = 0;
	while (!fin.eof() && count_v < Ins.num_v)
	{
		int x;
		fin >> x >> Ins.weight[count_v];
		count_v++;
	}

	if (max_edg != nb_edg)
	{
		cout << "### Error max_edge != nb_edg, nb_edg=" << nb_edg << ", max_edge=" << max_edg << endl;
		exit(-1);
	}
	if (count_v != Ins.num_v)
	{
		cout << "### Error count_v != Ins.num_v, Ins.num_v=" << Ins.num_v << ", count_v=" << count_v << endl;
		exit(-1);
	}
	cout << "finished read, running CDP" << endl;
	cout << "Instance_name=" << instance_name << ", Num_v=" << Ins.num_v
		 << ", capacity=" << Ins.capacity << endl << endl;

	fin.close();

#ifdef DEBUG_READ_INSTANCE
	cout << "num_v=" << Ins.num_v << endl;
	cout << "capacity=" << Ins.capacity << endl;
	cout << "weight=" << endl;
	for (int i = 0; i < Ins.num_v; i++)
		cout << Ins.weight[i] << " ";
	cout << endl << "distance=" << endl;;
	for (int i = 0; i < Ins.num_v; i++)
	{
		for (int j = 0; j < Ins.num_v; j++)
			cout << Ins.distance[i][j] << " ";
		cout << endl;
	}
#endif
}

//allocate memory for some global variables
void allocate_memory()
{
	Sol_best.ss = new int[Ins.num_v];
	Min_dis = new double[Ins.num_v];
	Sec_dis = new double[Ins.num_v];
	Num_min_dis = new int[Ins.num_v];
	Vec_min_dis = new int[Ins.num_v];
	W1 = new int[Ins.num_v];
	W2 = new int[Ins.num_v];
    W3 = new int[Ins.num_v];
	Hash1 = new int[L];
	Hash2 = new int[L];
	Hash3 = new int[L];
}


//compute the objective of a solution,
//time complexity: O(n * |M|) where |M| is the number of vertices selected to the current solution
double compute_obj(solution_data sol)
{
	double obj_min = MAXVALUE;
	for (int i = 0; i < Ins.num_v; i++)
	{
		if (sol.ss[i] == 1)
		{
			for (int j = i + 1; j < Ins.num_v; j++)
				if (sol.ss[j] == 1)
					if (Ins.distance[i][j] < obj_min - PRECISION)
						obj_min = Ins.distance[i][j];							
		}
	}
	return obj_min;
}

//generate a feasible initial solution in a greedy manner, time complexity: O(n * |M|)
void greedy_construct_initialsol(solution_data &sol)
{
	memset(sol.ss, 0, sizeof(int) * Ins.num_v);
	sol.num_sel = 0;
	sol.weight = 0;

	//choose a start point
	int select = rand() % Ins.num_v;
	sol.ss[select] = 1;
	sol.num_sel++;
	sol.weight += Ins.weight[select];
	sol.cost = 0;
	initialize_sup_arrays(sol);

	//greedy expansion
	while (sol.weight < Ins.capacity - PRECISION)
	{
		int best_array[MAXNUM];
		int best_len = 0;
		for (int i = 0; i < Ins.num_v; i++)
			if (sol.ss[i] == 0 && fabs(Min_dis[i] - Max_min_dis) <= PRECISION)
				best_array[best_len++] = i;

		if (best_len > 0)
		{
			select = best_array[rand() % best_len];
			sol.ss[select] = 1;
			sol.num_sel++;
			sol.weight += Ins.weight[select];
			initialize_sup_arrays(sol);
		}
	}
	sol.cost = compute_obj(sol);
//	verify_sol(sol); //verify solution

	//update global solution
	if (sol.cost > Sol_best.cost + PRECISION)
	{
		copy_solution(sol, Sol_best);
		Run_time = (double) (clock() - Start_time) / CLOCKS_PER_SEC;
	}
}

//generate a feasible initial solution in a random manner
void random_construct_initialsol(solution_data &sol)
{
	for (int i = 0; i < Ins.num_v; i++)
		sol.ss[i] = 0;
	sol.num_sel = 0;
	sol.weight = 0;

	//random expansion
	while (sol.weight < Ins.capacity - PRECISION)
	{
		int rand_v = rand() % Ins.num_v;
		if (sol.ss[rand_v] == 0)
		{
			sol.ss[rand_v] = 1;
			sol.num_sel++;
			sol.weight += Ins.weight[rand_v];
		}
	}
	sol.cost = compute_obj(sol);
//	verify_sol(sol); //verify solution

	//update global solution
	if (sol.cost > Sol_best.cost + PRECISION)
	{
		copy_solution(sol, Sol_best);
		Run_time = (double) (clock() - Start_time) / CLOCKS_PER_SEC;
	}
}

//copy solution_data structure
void copy_solution(solution_data source, solution_data &des)
{
	memcpy(des.ss, source.ss, sizeof(int) * Ins.num_v);
	des.num_sel = source.num_sel;
	des.cost = source.cost;
	des.weight = source.weight;
}

//initialize the global variables: Min_dis, Sec_dis, Num_min_dis, Vec_min_dis, time complexity: O(n * |M|)
void initialize_sup_arrays(solution_data sol)
{
	Max_min_dis = -MAXVALUE;
	Sec_min_dis = -MAXVALUE;
	for (int i = 0; i < Ins.num_v; i++)
	{
		Min_dis[i] = MAXVALUE;
		Sec_dis[i] = MAXVALUE;
		Num_min_dis[i] = 0;
		Vec_min_dis[i] = -1;
	}

	for (int i = 0; i < Ins.num_v; i++)
	{
		double dis_min = MAXVALUE;
		int select = -1;
		for (int j = 0; j < Ins.num_v; j++)
		{
			if (sol.ss[j] && i != j)
			{
				if (Ins.distance[i][j] < dis_min - PRECISION)
				{
					dis_min = Ins.distance[i][j];
					select = j;
					Num_min_dis[i] = 1;
				}
				else if (fabs(Ins.distance[i][j] - dis_min) <= PRECISION)
					Num_min_dis[i]++;
			}
		}
		Min_dis[i] = dis_min;
		Vec_min_dis[i] = select; //only restore the first encountered vertex with dis_min

		//calculate Max_min_dis
		if (sol.ss[i] == 0 && Min_dis[i] > Max_min_dis + PRECISION)
			Max_min_dis = Min_dis[i];
	}

	for (int i = 0; i < Ins.num_v; i++)
	{
		double dis_second_min = MAXVALUE;
//		int select = -1;
		for (int j = 0; j < Ins.num_v; j++)
		{
			if (sol.ss[j] == 1 && i != j && j != Vec_min_dis[i] && Ins.distance[i][j] < dis_second_min - PRECISION)
			{
				dis_second_min = Ins.distance[i][j];
//				select = j;
			}
		}
		Sec_dis[i] = dis_second_min;

		//calculate Sec_min_dis
		if (sol.ss[i] == 0 && Min_dis[i] > Sec_min_dis + PRECISION && fabs(Min_dis[i] - Max_min_dis) > PRECISION)
			Sec_min_dis = Min_dis[i];
	}
}

void initial_hash()
{
	/*for (int i = 0; i < Ins.num_v; i++)
	{
		W1[i] = (int) pow(i + 1, 1.3);
		W2[i] = (int) pow(i + 1, 1.5);
		W3[i] = (int) pow(i + 1, 1.8);
	}

	for (int i = 0; i < Ins.num_v; i++)
	{
		W1[i] = (int) rand() % (50 * Ins.num_v);
		W2[i] = (int) rand() % (50 * Ins.num_v);
		W3[i] = (int) rand() % (50 * Ins.num_v);
	}*/

//	int a = 300, b = 400, c = 500;
	W1[0] = A;
	W2[0] = B;
	W3[0] = C;
	for (int i = 1; i < Ins.num_v; i++)
	{
		W1[i] = (int) W1[i - 1] + A + rand() % (A / 2);
		W2[i] = (int) W1[i - 1] + B + rand() % (B / 2);
		W3[i] = (int) W1[i - 1] + C + rand() % (C / 2);
	}
}

//compute move gain for f(M) for each candidate move
double compute_flip_move_gain(int node, solution_data sol)
{
	double move_gain = 0;
	/* ADD: V\M --> M,
	 * if Min_dis[node] >= sol.cost, then move_gain = 0
	 * otherwise, move_gain is decreased
	 */
	if (sol.ss[node] == 0)
	{
		if (Min_dis[node] > sol.cost + PRECISION || fabs(Min_dis[node] - sol.cost) <= PRECISION)
			move_gain = 0;
		else
			move_gain = Min_dis[node] - sol.cost;
	}
	/* DROP: M --> V\M, only two cases:
	 * case 1: if Min_dis[node] > 0, then move_gain = 0
	 * case 2: if Min_dis[node] = 0, we evaluate for how many nodes, its Min_dis[] = 0
	 */
	else
	{
		if (Min_dis[node] > sol.cost + PRECISION)
			move_gain = 0;
		else if (fabs(Min_dis[node] - sol.cost) <= PRECISION)
		{
			double cost_min = MAXVALUE;
			/* to evaluate if there are more than one node their Min_dis[] = sol.cost
			 * if only one node it Min_dis[] = sol.cost, then cost_min = Sec_dis[i] (changed)
			 * if more than one node its Min_dis[] = sol.cost, then cost_min = Min_dis[i] (= sol.cost) (unchanged)
			 */
			for (int i = 0; i < Ins.num_v; i++)
			{
				if (sol.ss[i] == 1 && i != node)
				{
					double poten_tmp = Min_dis[i];
					if (Vec_min_dis[i] == node)
						poten_tmp = Sec_dis[i];
					if (poten_tmp < cost_min - PRECISION)
						cost_min = poten_tmp;
				}
			}
			move_gain = cost_min - sol.cost;
		}
	}
	return move_gain;
}

double compute_swap_move_gain(int node_x, int node_y, solution_data sol)
{
	double move_gain = 0;

	if (Min_dis[node_x] > sol.cost + PRECISION || fabs(Min_dis[node_x] - sol.cost) <= PRECISION)
	{
		if (Min_dis[node_y] > sol.cost + PRECISION)
			move_gain = 0;
		else if (fabs(Min_dis[node_y] - sol.cost) <= PRECISION)
		{
			double cost_min = MAXVALUE;
			for (int i = 0; i < Ins.num_v; i++)
			{
				if (sol.ss[i] == 1 && i != node_y)
				{
					double poten_tmp = Min_dis[i];
					if (Vec_min_dis[i] == node_y)
						poten_tmp = Sec_dis[i];
					if (poten_tmp < cost_min - PRECISION)
						cost_min = poten_tmp;
				}
			}

			if (Min_dis[node_x] < cost_min - PRECISION)
			{
				if (fabs(Ins.distance[node_x][node_y] - Min_dis[node_x]) > PRECISION)
					move_gain = Min_dis[node_x] - sol.cost;
				else
				{
					if (Num_min_dis[node_x] == 1)
					{
						if (Sec_dis[node_x] < cost_min - PRECISION)
							move_gain = Sec_dis[node_x] - sol.cost;
						else
							move_gain = cost_min - sol.cost;
					}
					else
						move_gain = Min_dis[node_x] - sol.cost;
				}
			}
			else
				move_gain = cost_min - sol.cost;
		}
	}
	else
	{
		if (Min_dis[node_y] > sol.cost + PRECISION)
		{
			if (fabs(Ins.distance[node_x][node_y] - Min_dis[node_x]) > PRECISION)
				move_gain = Min_dis[node_x] - sol.cost;
			else
			{
				if (Num_min_dis[node_x] == 1)
					move_gain = 0;
				else
					move_gain = Min_dis[node_x] - sol.cost;
			}
		}
		else if (fabs(Min_dis[node_y] - sol.cost) <= PRECISION)
		{
			double cost_min = MAXVALUE;
			for (int i = 0; i < Ins.num_v; i++)
			{
				if (sol.ss[i] == 1 && i != node_y)
				{
					double poten_tmp = Min_dis[i];
					if (Vec_min_dis[i] == node_y)
						poten_tmp = Sec_dis[i];
					if (poten_tmp < cost_min - PRECISION)
						cost_min = poten_tmp;
				}
			}

			if (Min_dis[node_x] < cost_min - PRECISION)
			{
				if (fabs(Ins.distance[node_x][node_y] - Min_dis[node_x]) > PRECISION)
					move_gain = Min_dis[node_x] - sol.cost;
				else
				{
					if (Num_min_dis[node_x] == 1)
					{
						if (Sec_dis[node_x] < cost_min - PRECISION)
							move_gain = Sec_dis[node_x] - sol.cost;
						else
							move_gain = cost_min - sol.cost;
					}
					else
						move_gain = Min_dis[node_x] - sol.cost;
				}
			}
			else
				move_gain = cost_min - sol.cost;
		}
	}
	return move_gain;
}

/*//from small to large
void quick_sort(int *sort_array, int *node_array, int low, int high)
{
	int l, r, sa, na;

	if (low < high)
	{
		l = low, r = high, sa = sort_array[low], na = node_array[low];

		while (l < r)
		{
			//find first > x, from right to left
			while ((l < r) && (sort_array[r] >= sa)) //large to small: only need to modify <=
				r--;
			if (l < r)
			{
				sort_array[l] = sort_array[r];
				node_array[l++] = node_array[r];
			}

			//find first <= x, from left to right
			while ((l < r) && (sort_array[l] < sa)) //large to small: only need to modify >
				l++;
			if (l < r)
			{
				sort_array[r] = sort_array[l];
				node_array[r--] = node_array[l];
			}
		}

		sort_array[l] = sa;
		node_array[l] = na;
		quick_sort(node_array, sort_array, low, l - 1);
		quick_sort(node_array, sort_array, l + 1, high);
	}
}*/

void verify_sol(solution_data sol)
{
	double weight = 0;
	double cost = MAXVALUE;
	int num_sel = 0;
	for (int i = 0; i < Ins.num_v; i++)
	{
		if (sol.ss[i])
		{
			num_sel++;
			weight += Ins.weight[i];
		}
	}
	cost = compute_obj(sol);

	//verify
	if (weight < Ins.capacity - PRECISION)
	{
		cout << "sol's weight < Capacity, an infeasible solution obtained!!!!!!" << endl;
		exit(-1);
	}
	if (fabs(sol.weight - weight) > PRECISION)
	{
		cout << "sol.weight=" << sol.weight << ", weight=" << weight << endl;
		cout << "sol.weight != weight, an infeasible solution obtained!!!!!!" << endl;
		exit(-1);
	}
	if (fabs(sol.cost - cost) > PRECISION)
	{
		cout << "sol.cost=" << sol.cost << ", cost=" << cost << endl;
		cout << "sol.cost != cost, an infeasible solution obtained!!!!!!" << endl;
		exit(-1);
	}
	if (sol.num_sel != num_sel)
	{
		cout << "sol.num_sel != num_sel, an infeasible solution obtained!!!!!!" << endl;
		exit(-1);
	}	
}

//output the found best solution of each run
void out_sol(solution_data sol, char *ins_name, char *out_sol_file)
{
	FILE *fp;
	char buff[MAXNUM];
	sprintf(buff, "%s", out_sol_file);
	fp = fopen(buff, "a+");
	fprintf(fp, "%s Num_v = %d Cap = %lf |S| = %d cost = %lf weight = %lf \n", ins_name, Ins.num_v, Ins.capacity, sol.num_sel, sol.cost, sol.weight);
	for (int i = 0; i < Ins.num_v; i++)
		if (sol.ss[i])
			fprintf(fp, "%d ", i);
	fprintf(fp, "\n");
	fclose(fp);
}

//output the found best results of each run
void out_stat_results(char *ins_name, char *out_stat_file, double best_result, double run_time)
{
	FILE *fp;
	char buff[MAXNUM];
	sprintf(buff, "%s", out_stat_file);
	fp = fopen(buff, "a+");
	fprintf(fp, "%s %d %lf %lf %lf \n", ins_name, Ins.num_v, Ins.capacity, best_result, run_time);
	fclose(fp);
}

void out_total_results(char *ins_name, char *out_total_file, double *cost_total, double *time_total)
{
	double fbest = -MAXVALUE, fworst = MAXVALUE, favg = 0.0;
	double ftime = 0.0, std = 0.0;
	int hit = 0;

	for (int i = 0; i < Runs; i++)
	{
		favg += cost_total[i];
		ftime += time_total[i];
	}
	favg /= Runs;
	ftime /= Runs;
	for (int i = 0; i < Runs; i++)
	{
		if (cost_total[i] > fbest + PRECISION || fabs(cost_total[i] - fbest) <= PRECISION)
			fbest = cost_total[i];
		if (cost_total[i] < fworst - PRECISION || fabs(cost_total[i] - fworst) <= PRECISION)
			fworst = cost_total[i];
	}
	for (int i = 0; i < Runs; i++)
	{
		if (fabs(cost_total[i] - fbest) <= PRECISION)
			hit++;
	}
	for (int i = 0; i < Runs; i++)
		std += pow((cost_total[i] - favg), 2);
	std /= Runs;
	std = sqrt(std);

	FILE *fp;
	char buff[MAXNUM];
	sprintf(buff, "%s", out_total_file);
	fp = fopen(buff, "a+");
	fprintf(fp, "%s %d %lf %lf %lf %lf %lf %d//%d %lf %lf\n", ins_name, Ins.num_v, Ins.capacity, Fprev, fbest, favg, fworst, hit, Runs, ftime, std);
	fclose(fp);
}

void free_memory()
{
	delete[] Ins.weight;
	for (int i = 0; i < Ins.num_v; i++)
		delete[] Ins.distance[i];
	delete[] Ins.distance;

	delete[] Sol_best.ss;
	delete[] Min_dis;
	delete[] Sec_dis;
	delete[] Num_min_dis;
	delete[] Vec_min_dis;

	delete[] W1;
	delete[] W2;
	delete[] W3;
	delete[] Hash1;
	delete[] Hash2;
	delete[] Hash3;
}
