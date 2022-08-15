#include "func_state.h"
#include "global_variables.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

using namespace std;


//solution based tabu search: Add, Drop, and Swap move operators with hash function
void solution_based_tabu_search(solution_data &sol)
{
	//initial hash
	int hx1 = 0, hx2 = 0, hx3 = 0;
	int hx11, hx22, hx33;

	for (int i = 0; i < L; i++)
	{
		Hash1[i] = 0;
		Hash2[i] = 0;
		Hash3[i] = 0;
	}
	for (int i = 0; i < Ins.num_v; i++)
	{
		if(sol.ss[i])
		{
			hx1 += W1[i];
			hx2 += W2[i];
			hx3 += W3[i];
		}
	}

	int iters = 0;
	int non_improve = 0;
	neighbor_data best_array[MAXNUM];
	neighbor_data tabu_best_array[MAXNUM];
	solution_data loc_best_sol;
	loc_best_sol.ss = new int[Ins.num_v];
	copy_solution(sol, loc_best_sol);
	initialize_sup_arrays(sol);

	while ((double) (clock() - Start_time) / CLOCKS_PER_SEC < Time_limit)
	{
		double best_gain = -MAXVALUE;
		double tabu_best_gain = -MAXVALUE;
		int best_len = 0;
		int tabu_best_len = 0;

		//1. Add
		for (int i = 0; i < Ins.num_v; i++)
		{
			if (sol.ss[i] == 0)
			{
				double cur_gain = compute_flip_move_gain(i, sol);

				hx11 = hx1 + W1[i];
				hx22 = hx2 + W2[i];
				hx33 = hx3 + W3[i];

				if (Hash1[hx11 % L] == 1 && Hash2[hx22 % L] == 1 && Hash3[hx33 % L] == 1)
				{
					if (cur_gain > tabu_best_gain + PRECISION)
					{
						tabu_best_array[0].gain = cur_gain;
						tabu_best_array[0].node = i;
						tabu_best_array[0].type = 1;
						tabu_best_array[0].flip = 1;
						tabu_best_gain = cur_gain;
						tabu_best_len = 1;
					}
					else if (fabs(cur_gain - tabu_best_gain) <= PRECISION && tabu_best_len < MAXNUM)
					{
						tabu_best_array[tabu_best_len].gain = cur_gain;
						tabu_best_array[tabu_best_len].node = i;
						tabu_best_array[tabu_best_len].type = 1;
						tabu_best_array[tabu_best_len].flip = 1;
						tabu_best_len++;
					}
				}
				else
				{
					if (cur_gain > best_gain + PRECISION)
					{
						best_array[0].gain = cur_gain;
						best_array[0].node = i;
						best_array[0].type = 1;
						best_array[0].flip = 1;
						best_gain = cur_gain;
						best_len = 1;
					}
					else if (fabs(cur_gain - best_gain) <= PRECISION && best_len < MAXNUM)
					{
						best_array[best_len].gain = cur_gain;
						best_array[best_len].node = i;
						best_array[best_len].type = 1;
						best_array[best_len].flip = 1;
						best_len++;
					}
				}
			}
		}

		//2. Drop
		for (int i = 0; i < Ins.num_v; i++)
		{
			if (sol.ss[i] == 1 && ((sol.weight - Ins.weight[i]) > Ins.capacity + PRECISION
				|| fabs(sol.weight - Ins.weight[i] - Ins.capacity) <= PRECISION))
			{
				double cur_gain = compute_flip_move_gain(i, sol);

				hx11 = hx1 - W1[i];
				hx22 = hx2 - W2[i];
				hx33 = hx3 - W3[i];

				if (Hash1[hx11 % L] == 1 && Hash2[hx22 % L] == 1 && Hash3[hx33 % L] == 1)
				{
					if (cur_gain > tabu_best_gain + PRECISION)
					{
						tabu_best_array[0].gain = cur_gain;
						tabu_best_array[0].node = i;
						tabu_best_array[0].type = 1;
						tabu_best_array[0].flip = 2;
						tabu_best_gain = cur_gain;
						tabu_best_len = 1;
					}
					else if (fabs(cur_gain - tabu_best_gain) <= PRECISION && tabu_best_len < MAXNUM)
					{
						tabu_best_array[tabu_best_len].gain = cur_gain;
						tabu_best_array[tabu_best_len].node = i;
						tabu_best_array[tabu_best_len].type = 1;
						tabu_best_array[tabu_best_len].flip = 2;
						tabu_best_len++;
					}
				}
				else
				{
					if (cur_gain > best_gain + PRECISION)
					{
						best_array[0].gain = cur_gain;
						best_array[0].node = i;
						best_array[0].type = 1;
						best_array[0].flip = 2;
						best_gain = cur_gain;
						best_len = 1;
					}
					else if (fabs(cur_gain - best_gain) <= PRECISION && best_len < MAXNUM)
					{
						best_array[best_len].gain = cur_gain;
						best_array[best_len].node = i;
						best_array[best_len].type = 1;
						best_array[best_len].flip = 2;
						best_len++;
					}
				}
			}
		}

		//3. Swap
		int add_array_swap[MAXNUM];
		int drop_array_swap[MAXNUM];
		int add_len_swap = 0;
		int drop_len_swap = 0;
		for (int i = 0; i < Ins.num_v; i++)
		{
			//Find add nodes with max and sec Min_dis[i]
			if (sol.ss[i] == 0 && (fabs(Min_dis[i] - Max_min_dis) <= PRECISION
				|| fabs(Min_dis[i] - Sec_min_dis) <= PRECISION))
				add_array_swap[add_len_swap++] = i;
			//Find all drop nodes with a distance equal to sol.cost
			if (sol.ss[i] == 1 && fabs(Min_dis[i] - sol.cost) <= PRECISION)
				drop_array_swap[drop_len_swap++] = i;
		}
		for (int i = 0; i < add_len_swap; i++)
		{
			int nx = add_array_swap[i];
			for (int j = 0; j < drop_len_swap; j++)
			{
				int ny = drop_array_swap[j];
				if ((sol.weight + Ins.weight[nx] - Ins.weight[ny]) < Ins.capacity - PRECISION)
					continue;
				double cur_gain = compute_swap_move_gain(nx, ny, sol);

				//SWAP
				hx11 = hx1 + (W1[nx] - W1[ny]);
				hx22 = hx2 + (W2[nx] - W2[ny]);
				hx33 = hx3 + (W3[nx] - W3[ny]);

				if (Hash1[hx11 % L] == 1 && Hash2[hx22 % L] == 1 && Hash3[hx33 % L] == 1)
				{
					if (cur_gain > tabu_best_gain + PRECISION)
					{
						tabu_best_array[0].gain = cur_gain;
						tabu_best_array[0].node_x = nx;
						tabu_best_array[0].node_y = ny;
						tabu_best_array[0].type = 2;
						tabu_best_array[0].flip = -1;
						tabu_best_gain = cur_gain;
						tabu_best_len = 1;
					}
					else if (fabs(cur_gain - tabu_best_gain) <= PRECISION && tabu_best_len < MAXNUM)
					{
						tabu_best_array[tabu_best_len].gain = cur_gain;
						tabu_best_array[tabu_best_len].node_x = nx;
						tabu_best_array[tabu_best_len].node_y = ny;
						tabu_best_array[tabu_best_len].type = 2;
						tabu_best_array[tabu_best_len].flip = -1;
						tabu_best_len++;
					}
				}
				else
				{
					if (cur_gain > best_gain + PRECISION)
					{
						best_array[0].gain = cur_gain;
						best_array[0].node_x = nx;
						best_array[0].node_y = ny;
						best_array[0].type = 2;
						best_array[0].flip = -1;
						best_gain = cur_gain;
						best_len = 1;
					}
					else if (fabs(cur_gain - best_gain) <= PRECISION && best_len < MAXNUM)
					{
						best_array[best_len].gain = cur_gain;
						best_array[best_len].node_x = nx;
						best_array[best_len].node_y = ny;
						best_array[best_len].type = 2;
						best_array[best_len].flip = -1;
						best_len++;
					}
				}
			}
		}

		//Neighboring solution transition
		int select = -1;
		int node = -1, node_x = -1, node_y = -1;
		if ((tabu_best_len > 0 && tabu_best_gain > best_gain + PRECISION
			&& (tabu_best_gain + sol.cost) > loc_best_sol.cost + PRECISION)
			|| (best_len == 0 && tabu_best_len > 0))
		{
			select = rand() % tabu_best_len;
			if (tabu_best_array[select].type == 2) //Swap
			{
				node_x = tabu_best_array[select].node_x;
				node_y = tabu_best_array[select].node_y;
				sol.cost += tabu_best_array[select].gain;
				sol.weight = sol.weight + Ins.weight[node_x] - Ins.weight[node_y];
				sol.ss[node_x] = 1;
				sol.ss[node_y] = 0;

				hx1 = hx1 + (W1[node_x] - W1[node_y]);
				hx2 = hx2 + (W2[node_x] - W2[node_y]);
				hx3 = hx3 + (W3[node_x] - W3[node_y]);
				Hash1[hx1 % L] = 1;
				Hash2[hx2 % L] = 1;
				Hash3[hx3 % L] = 1;
			}
			else if (tabu_best_array[select].type == 1) //Flip
			{
				if (tabu_best_array[select].flip == 1) //Add
				{
					node = tabu_best_array[select].node;
					sol.cost += tabu_best_array[select].gain;
					sol.num_sel++;
					sol.weight += Ins.weight[node];
					sol.ss[node] = 1;

					hx1 += W1[node];
					hx2 += W2[node];
					hx3 += W3[node];
					Hash1[hx1 % L] = 1;
					Hash2[hx2 % L] = 1;
					Hash3[hx3 % L] = 1;
				}
				else if (tabu_best_array[select].flip == 2) //Drop
				{
					node = tabu_best_array[select].node;
					sol.cost += tabu_best_array[select].gain;
					sol.num_sel--;
					sol.weight -= Ins.weight[node];
					sol.ss[node] = 0;

					hx1 -= W1[node];
					hx2 -= W2[node];
					hx3 -= W3[node];
					Hash1[hx1 % L] = 1;
					Hash2[hx2 % L] = 1;
					Hash3[hx3 % L] = 1;
				}
			}
		}
		else
		{
			select = rand() % best_len;
			if (best_array[select].type == 2) //Swap
			{
				node_x = best_array[select].node_x;
				node_y = best_array[select].node_y;
				sol.cost += best_array[select].gain;
				sol.weight = sol.weight + Ins.weight[node_x] - Ins.weight[node_y];
				sol.ss[node_x] = 1;
				sol.ss[node_y] = 0;

				hx1 = hx1 + (W1[node_x] - W1[node_y]);
				hx2 = hx2 + (W2[node_x] - W2[node_y]);
				hx3 = hx3 + (W3[node_x] - W3[node_y]);
				Hash1[hx1 % L] = 1;
				Hash2[hx2 % L] = 1;
				Hash3[hx3 % L] = 1;
			}
			else if (best_array[select].type == 1) //Flip
			{
				if (best_array[select].flip == 1) //Add
				{
					node = best_array[select].node;
					sol.cost += best_array[select].gain;
					sol.num_sel++;
					sol.weight += Ins.weight[node];
					sol.ss[node] = 1;

					hx1 += W1[node];
					hx2 += W2[node];
					hx3 += W3[node];
					Hash1[hx1 % L] = 1;
					Hash2[hx2 % L] = 1;
					Hash3[hx3 % L] = 1;
				}
				else if (best_array[select].flip == 2) //Drop
				{
					node = best_array[select].node;
					sol.cost += best_array[select].gain;
					sol.num_sel--;
					sol.weight -= Ins.weight[node];
					sol.ss[node] = 0;

					hx1 -= W1[node];
					hx2 -= W2[node];
					hx3 -= W3[node];
					Hash1[hx1 % L] = 1;
					Hash2[hx2 % L] = 1;
					Hash3[hx3 % L] = 1;
				}
			}
		}
		if (node != -1 || node_x != -1 || node_y != -1)
		{
			initialize_sup_arrays(sol);
	//		verify_sol(sol); //verify solution
		}

		//update local best solution
		if (sol.cost > loc_best_sol.cost + PRECISION)
		{
			copy_solution(sol, loc_best_sol);
			Run_time = (double) (clock() - Start_time) / CLOCKS_PER_SEC;
			non_improve = 0;
		}
		else
			non_improve++;
		if (fabs(Fprev - loc_best_sol.cost) <= PRECISION)
			break;

		if (iters % 2000 == 0)
			cout << "iters=" << iters << ", non_improve=" << non_improve
				 << ", sol.num_sel=" << sol.num_sel << ", sol=" << sol.cost
			     << ", local best=" << loc_best_sol.cost << ", Run_time=" << Run_time << endl;
		iters++;
	}
	//update global solution
	if (loc_best_sol.cost > Sol_best.cost + PRECISION)
		copy_solution(loc_best_sol, Sol_best);

	//free
	delete[] loc_best_sol.ss;
}


void attribute_based_tabu_search(solution_data &sol)
{
	int iters = 0;
	int non_improve = 0;
	neighbor_data best_array[MAXNUM];
	neighbor_data tabu_best_array[MAXNUM];
	solution_data loc_best_sol;
	loc_best_sol.ss = new int[Ins.num_v];
	copy_solution(sol, loc_best_sol);
	initialize_sup_arrays(sol);

	int tabu_tenure = 100;
	int *tabu_list = new int[Ins.num_v];
	memset(tabu_list, 0, sizeof(int) * Ins.num_v);

	while ((double) (clock() - Start_time) / CLOCKS_PER_SEC < Time_limit)
	{
		double best_gain = -MAXVALUE;
		double tabu_best_gain = -MAXVALUE;
		int best_len = 0;
		int tabu_best_len = 0;

		//1. Add
		for (int i = 0; i < Ins.num_v; i++)
		{
			if (sol.ss[i] == 0)
			{
				double cur_gain = compute_flip_move_gain(i, sol);

				if (tabu_list[i] > iters)
				{
					if (cur_gain > tabu_best_gain + PRECISION)
					{
						tabu_best_array[0].gain = cur_gain;
						tabu_best_array[0].node = i;
						tabu_best_array[0].type = 1;
						tabu_best_array[0].flip = 1;
						tabu_best_gain = cur_gain;
						tabu_best_len = 1;
					}
					else if (fabs(cur_gain - tabu_best_gain) <= PRECISION && tabu_best_len < MAXNUM)
					{
						tabu_best_array[tabu_best_len].gain = cur_gain;
						tabu_best_array[tabu_best_len].node = i;
						tabu_best_array[tabu_best_len].type = 1;
						tabu_best_array[tabu_best_len].flip = 1;
						tabu_best_len++;
					}
				}
				else
				{
					if (cur_gain > best_gain + PRECISION)
					{
						best_array[0].gain = cur_gain;
						best_array[0].node = i;
						best_array[0].type = 1;
						best_array[0].flip = 1;
						best_gain = cur_gain;
						best_len = 1;
					}
					else if (fabs(cur_gain - best_gain) <= PRECISION && best_len < MAXNUM)
					{
						best_array[best_len].gain = cur_gain;
						best_array[best_len].node = i;
						best_array[best_len].type = 1;
						best_array[best_len].flip = 1;
						best_len++;
					}
				}
			}
		}

		//2. Drop
		for (int i = 0; i < Ins.num_v; i++)
		{
			if (sol.ss[i] == 1 && ((sol.weight - Ins.weight[i]) > Ins.capacity + PRECISION
				|| fabs(sol.weight - Ins.weight[i] - Ins.capacity) <= PRECISION))
			{
				double cur_gain = compute_flip_move_gain(i, sol);

				if (tabu_list[i] > iters)
				{
					if (cur_gain > tabu_best_gain + PRECISION)
					{
						tabu_best_array[0].gain = cur_gain;
						tabu_best_array[0].node = i;
						tabu_best_array[0].type = 1;
						tabu_best_array[0].flip = 2;
						tabu_best_gain = cur_gain;
						tabu_best_len = 1;
					}
					else if (fabs(cur_gain - tabu_best_gain) <= PRECISION && tabu_best_len < MAXNUM)
					{
						tabu_best_array[tabu_best_len].gain = cur_gain;
						tabu_best_array[tabu_best_len].node = i;
						tabu_best_array[tabu_best_len].type = 1;
						tabu_best_array[tabu_best_len].flip = 2;
						tabu_best_len++;
					}
				}
				else
				{
					if (cur_gain > best_gain + PRECISION)
					{
						best_array[0].gain = cur_gain;
						best_array[0].node = i;
						best_array[0].type = 1;
						best_array[0].flip = 2;
						best_gain = cur_gain;
						best_len = 1;
					}
					else if (fabs(cur_gain - best_gain) <= PRECISION && best_len < MAXNUM)
					{
						best_array[best_len].gain = cur_gain;
						best_array[best_len].node = i;
						best_array[best_len].type = 1;
						best_array[best_len].flip = 2;
						best_len++;
					}
				}
			}
		}

		//3. Swap
		int add_array_swap[MAXNUM];
		int drop_array_swap[MAXNUM];
		int add_len_swap = 0;
		int drop_len_swap = 0;
		for (int i = 0; i < Ins.num_v; i++)
		{
			//Find add nodes with max and sec Min_dis[i]
			if (sol.ss[i] == 0 && (fabs(Min_dis[i] - Max_min_dis) <= PRECISION
				|| fabs(Min_dis[i] - Sec_min_dis) <= PRECISION))
				add_array_swap[add_len_swap++] = i;
			//Find all drop nodes with a distance equal to sol.cost
			if (sol.ss[i] == 1 && fabs(Min_dis[i] - sol.cost) <= PRECISION)
				drop_array_swap[drop_len_swap++] = i;
		}
		for (int i = 0; i < add_len_swap; i++)
		{
			int nx = add_array_swap[i];
			for (int j = 0; j < drop_len_swap; j++)
			{
				int ny = drop_array_swap[j];
				if ((sol.weight + Ins.weight[nx] - Ins.weight[ny]) < Ins.capacity - PRECISION)
					continue;
				double cur_gain = compute_swap_move_gain(nx, ny, sol);

				if (tabu_list[nx] > iters && tabu_list[ny] > iters)
				{
					if (cur_gain > tabu_best_gain + PRECISION)
					{
						tabu_best_array[0].gain = cur_gain;
						tabu_best_array[0].node_x = nx;
						tabu_best_array[0].node_y = ny;
						tabu_best_array[0].type = 2;
						tabu_best_array[0].flip = -1;
						tabu_best_gain = cur_gain;
						tabu_best_len = 1;
					}
					else if (fabs(cur_gain - tabu_best_gain) <= PRECISION && tabu_best_len < MAXNUM)
					{
						tabu_best_array[tabu_best_len].gain = cur_gain;
						tabu_best_array[tabu_best_len].node_x = nx;
						tabu_best_array[tabu_best_len].node_y = ny;
						tabu_best_array[tabu_best_len].type = 2;
						tabu_best_array[tabu_best_len].flip = -1;
						tabu_best_len++;
					}
				}
				else
				{
					if (cur_gain > best_gain + PRECISION)
					{
						best_array[0].gain = cur_gain;
						best_array[0].node_x = nx;
						best_array[0].node_y = ny;
						best_array[0].type = 2;
						best_array[0].flip = -1;
						best_gain = cur_gain;
						best_len = 1;
					}
					else if (fabs(cur_gain - best_gain) <= PRECISION && best_len < MAXNUM)
					{
						best_array[best_len].gain = cur_gain;
						best_array[best_len].node_x = nx;
						best_array[best_len].node_y = ny;
						best_array[best_len].type = 2;
						best_array[best_len].flip = -1;
						best_len++;
					}
				}
			}
		}

		//Neighboring solution transition
		int select = -1;
		int node = -1, node_x = -1, node_y = -1;
		if ((tabu_best_len > 0 && tabu_best_gain > best_gain + PRECISION
			&& (tabu_best_gain + sol.cost) > loc_best_sol.cost + PRECISION)
			|| (best_len == 0 && tabu_best_len > 0))
		{
			select = rand() % tabu_best_len;
			if (tabu_best_array[select].type == 2) //Swap
			{
				node_x = tabu_best_array[select].node_x;
				node_y = tabu_best_array[select].node_y;
				sol.cost += tabu_best_array[select].gain;
				sol.weight = sol.weight + Ins.weight[node_x] - Ins.weight[node_y];
				sol.ss[node_x] = 1;
				sol.ss[node_y] = 0;

				tabu_list[node_x] = tabu_tenure + iters;
				tabu_list[node_y] = tabu_tenure + iters;
			}
			else if (tabu_best_array[select].type == 1) //Flip
			{
				if (tabu_best_array[select].flip == 1) //Add
				{
					node = tabu_best_array[select].node;
					sol.cost += tabu_best_array[select].gain;
					sol.num_sel++;
					sol.weight += Ins.weight[node];
					sol.ss[node] = 1;

					tabu_list[node] = tabu_tenure + iters;
				}
				else if (tabu_best_array[select].flip == 2) //Drop
				{
					node = tabu_best_array[select].node;
					sol.cost += tabu_best_array[select].gain;
					sol.num_sel--;
					sol.weight -= Ins.weight[node];
					sol.ss[node] = 0;

					tabu_list[node] = tabu_tenure + iters;
				}
			}
		}
		else
		{
			select = rand() % best_len;
			if (best_array[select].type == 2) //Swap
			{
				node_x = best_array[select].node_x;
				node_y = best_array[select].node_y;
				sol.cost += best_array[select].gain;
				sol.weight = sol.weight + Ins.weight[node_x] - Ins.weight[node_y];
				sol.ss[node_x] = 1;
				sol.ss[node_y] = 0;

				tabu_list[node_x] = tabu_tenure + iters;
				tabu_list[node_y] = tabu_tenure + iters;
			}
			else if (best_array[select].type == 1) //Flip
			{
				if (best_array[select].flip == 1) //Add
				{
					node = best_array[select].node;
					sol.cost += best_array[select].gain;
					sol.num_sel++;
					sol.weight += Ins.weight[node];
					sol.ss[node] = 1;

					tabu_list[node] = tabu_tenure + iters;
				}
				else if (best_array[select].flip == 2) //Drop
				{
					node = best_array[select].node;
					sol.cost += best_array[select].gain;
					sol.num_sel--;
					sol.weight -= Ins.weight[node];
					sol.ss[node] = 0;

					tabu_list[node] = tabu_tenure + iters;
				}
			}
		}
		if (node != -1 || node_x != -1 || node_y != -1)
		{
			initialize_sup_arrays(sol);
	//		verify_sol(sol); //verify solution
		}

		//update local best solution
		if (sol.cost > loc_best_sol.cost + PRECISION)
		{
			copy_solution(sol, loc_best_sol);
			Run_time = (double) (clock() - Start_time) / CLOCKS_PER_SEC;
			non_improve = 0;
		}
		else
			non_improve++;
		if (fabs(Fprev - loc_best_sol.cost) <= PRECISION)
			break;

		if (iters % 2000 == 0)
			cout << "iters=" << iters << ", non_improve=" << non_improve
				 << ", sol.num_sel=" << sol.num_sel << ", sol=" << sol.cost
			     << ", local best=" << loc_best_sol.cost << ", Run_time=" << Run_time << endl;
		iters++;
	}
	//update global solution
	if (loc_best_sol.cost > Sol_best.cost + PRECISION)
		copy_solution(loc_best_sol, Sol_best);

	//free
	delete[] loc_best_sol.ss;
	delete[] tabu_list;
}
