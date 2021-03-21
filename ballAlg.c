#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>

#define RANGE 10

typedef struct Node
{
    double **pts;
    double *center;
    double radius;
    int n_dims;
    int n_points;

    struct Node *left;
    struct Node *right;

} node;

node* newNode(double **arr)
{   
    node* aux = (node*)malloc(sizeof(node));
    aux->pts = arr;
    aux->n_points = sizeof(aux->pts);
    aux->n_dims = sizeof(aux->pts[0]);
    
    aux->center = NULL;
    aux->radius = 0;
    // Left and right child for node
    // will be initialized to null
    aux->left = NULL;
    aux->right = NULL;
    
    return aux;
}


double comp_dist(double *a, double *b , int *n_dims){
    double aux, sum = 0, two = 2;
    for(int i = 0; i < *n_dims; i++){
        aux = a[i] - b[i];
        aux = pow(aux, two);
        sum += aux;
    }
    return sqrt(sum); 
}

double furthest(double **pt_arr,int n_points){
    double *aux = pt_arr[0];
    print("%d", sizeof(aux));
    for (int i=0; i< sizeof(aux); i++){
        printf("%f ", aux[i]);
    }
    for( int i=0; i< 2; i++){
        //printf("%f ", pt_arr[0][i]);
    }

}

double **create_array_pts(int n_dims, long np)
{
    double *_p_arr;
    double **p_arr;
    _p_arr = (double *)malloc(n_dims * np * sizeof(double));
    p_arr = (double **)malloc(np * sizeof(double *));
    if ((_p_arr == NULL) || (p_arr == NULL))
    {
        printf("Error allocating array of points, exiting.\n");
        exit(4);
    }
    for (int i = 0; i < np; i++)
        p_arr[i] = &_p_arr[i * n_dims];
    return p_arr;
}

double **get_points(int argc, char *argv[], int *n_dims, long *np)
{
    double **pt_arr;
    unsigned seed;
    long i;
    int j;
    if (argc != 4)
    {
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }
    *n_dims = atoi(argv[1]);
    if (*n_dims < 2)
    {
        printf("Illegal number of dimensions (%d), must be above 1.\n", *n_dims);
        exit(2);
    }
    *np = atol(argv[2]);
    if (*np < 1)
    {
        printf("Illegal number of points (%ld), must be above 0.\n", *np);
        exit(3);
    }
    seed = atoi(argv[3]);
    srandom(seed);
    pt_arr = create_array_pts(*n_dims, *np);
    for (i = 0; i < *np; i++)
        for (j = 0; j < *n_dims; j++)
            pt_arr[i][j] = RANGE * ((double)random()) / RAND_MAX;
    return pt_arr;
}

void build_tree(node* root)
{
    
    for (int i = 0; i < root->n_points; i++)
    {
        for (int j = 0; j < root->n_dims; j++)
            printf("%f ", root->pts[i][j]);
        printf("\n");
    }
}


int main(int argc, char *argv[])
{
    int n_dims = 0;
    long n_points = 0;
    double **pts;

    double exec_time;
    exec_time = -omp_get_wtime();
    pts = get_points(argc, argv, &n_dims, &n_points);

    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_dims; j++)
            printf("%f ", pts[i][j]);
        printf("\n");
    }

    furthest(pts,n_points);

    node* root = newNode(pts);
    //build_tree(root);
    //build tree

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1lf\n", exec_time);

    //dump tree
}