#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H

//#define DEBUG_READ_INSTANCE
#define PRECISION 1.0e-6
#define MAXVALUE 99999999
#define MAXNUM 500

typedef struct instance_data
{
	int num_v;					//number of vertices
	double capacity;			//capacity constraint
	double *weight;				//weight of each vertex
	double **distance;			//distance between each pair of vertices
}instance_data;

typedef struct solution_data
{
	double cost;				//the objective function value
	double weight;				//the total weight of the solution
	int num_sel;				//the number of selected vertices of the solution
	int	*ss;					//representation of the solution, ss[i] = 1 indicates that vertex i is selected, 0 otherwise
}solution_data;

typedef struct neighbor_data
{
	double gain;				//the move gain
    int type;					//neighborhood type: Flip (Add, Drop)=1, Swap=2
    int flip;					//Flip type: Add=1, Drop=2
    int node;
    int node_x;
    int node_y;
}neighbor_data;

//global variables are all started with Capital letter
extern instance_data Ins;		//instance data
extern solution_data Sol_best;	//best solution found for each run

extern double Fprev;
extern double *Min_dis;			//Min_dis[i] denotes the minimum distance between vertex i and the vertices selected to current solution, i \in V
extern double *Sec_dis;			//Sec_dis[i] denotes the second minimum distance between vertex i and the vertices selected to current solution
extern int *Num_min_dis;		//Num_min_dis[i] denotes the number of edges (neighbors) whose distance is equal to Min_dis[i]
extern int *Vec_min_dis;		//Vec_min_dis[i] indicates any vertex v such that the distance between v and i is equal to Min_dis[i]
extern double Max_min_dis;
extern double Sec_min_dis;

extern int *Hash1, *Hash2, *Hash3;
extern int *W1, *W2, *W3;
extern int A, B, C;
extern long L;

extern int Runs;
extern double Start_time, Run_time, Time_limit;


#endif
