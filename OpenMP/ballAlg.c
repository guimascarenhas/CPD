#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>

#define RANGE 10
#define OMP_NUM_THREADS 8

typedef struct Node
{
    double **pts;
    double *center;
    double radius;
    int n_dims;
    int n_points;
    int id;

    struct Node *left;
    struct Node *right;

} node;

node *newNode(double **arr, long n_points, int n_dims, int id)
{
    node *aux = (node *)malloc(sizeof(node));
    aux->id = id;
    aux->pts = arr;
    aux->n_points = n_points;
    aux->n_dims = n_dims;

    aux->center = NULL;
    aux->radius = 0;
    // Left and right child for node
    // will be initialized to null
    aux->left = NULL;
    aux->right = NULL;

    return aux;
}



void dump_tree(node *root){
    if(root->left != NULL)
        dump_tree(root->left);

    if(root->right != NULL)
        dump_tree(root->right);
/*
    for(int i=0; i<root->n_points; i++){
        free(root->pts[i]);
    }
    */
    free(root->pts);
    free(root->center);
    free(root);
    return;
}

double comp_dist(double *a, double *b, int *n_dims)
{
    double aux, sum = 0, two = 2;
    for (int i = 0; i < *n_dims; i++)
    {
        aux = a[i] - b[i];
        aux = pow(aux, two);
        sum += aux;
    }
    return sqrt(sum);
}

double calcRadius(double **pt_arr, int n_dims, int *indices, double *center)
{
    double dist_a = comp_dist(pt_arr[indices[0]], center, &n_dims);
    double dist_b = comp_dist(pt_arr[indices[1]], center, &n_dims);

    if (dist_a > dist_b)
        return dist_a;
    else
        return dist_b;
}

/*
Recebe o conjunto de pontos (vetor)
Retorna os indices dos dois pontos mais afastados um do outro
*/
int *furthest(double **pt_arr, long n_points, int n_dims)
{
    double *aux = pt_arr[0];
    double dist_max_b = -1;
    double dist_max_a = -1;
    int ind_a = 0, ind_b = 0;
    long i=0;
    
    /*
    Forma mais eficaz é cada thread fazer o seu cálculo de dist_max
    e depois comparar o resultado de todas elas (evita race condition)
    */
    #pragma omp parallel private(i)
    {
    double dist = 0;
    double dist_max_local = 0;
    long ind_local = 0;

    #pragma omp for nowait
    for (long i = 1; i < n_points; i++)
    {   
        dist = comp_dist(aux, pt_arr[i], &n_dims);
        if (dist > dist_max_local)
        {
            dist_max_local = dist;
            ind_local = i;
        }
    }
    #pragma omp critical
    {
        if(dist_max_local > dist_max_b){
            ind_b = ind_local;
            dist_max_b = dist_max_local;
        }
    }

    #pragma omp barrier
    /*
    printf("dist_max: %f", dist_max);
    getchar();
    dist_max = 0;
    */

    dist = 0;
    dist_max_local = 0;
    ind_local = 0;

    #pragma omp for nowait
    for (long i = 0; i < n_points; i++)
    {
        if (i != ind_b)
        {
            dist = comp_dist(pt_arr[ind_b], pt_arr[i], &n_dims);
            if (dist > dist_max_local)
            {
                dist_max_local = dist;
                ind_local = i;
            }
        }
    }

    #pragma omp critical
    {
        if(dist_max_local > dist_max_a){
            ind_a = ind_local;
            dist_max_a = dist_max_local;
        }
    }
    }

    static int return_aux[2];
    return_aux[0] = ind_a;
    return_aux[1] = ind_b;
    return return_aux;
}

double inner_product(double *x, double *y, int n_dims)
{
    double result = 0;
    int i;
    for (i = 0; i < n_dims; i++)
    {
        result += x[i] * y[i];
    }
    return result;
}

double **ort_proj(double **pt_arr, long n_points, int n_dims, int *indices)
{
    double *a = pt_arr[indices[0]];
    double *b = pt_arr[indices[1]];
    long i;
    double aux1, aux2;
    double **return_ort = (double **)malloc(n_points * sizeof(double *));

    double *x = (double *)malloc(sizeof(double) * n_dims);
    double *y = (double *)malloc(sizeof(double) * n_dims);

    #pragma omp parallel private(i)
    {
    #pragma omp for
    for (long i = 0; i < n_points; i++)
    {
        return_ort[i] = (double *)malloc(sizeof(double) * n_dims);
    }
    }

    
    for (long i = 0; i < n_points; i++)
    {
        if (i != indices[0] && i != indices[1])
        {
            for (int j = 0; j < n_dims; j++)
            {
                x[j] = pt_arr[i][j] - a[j];
                y[j] = b[j] - a[j];
                
            }

            aux1 = inner_product(x, y, n_dims);
            aux2 = inner_product(y, y, n_dims);
            aux1 = aux1 / aux2;
            for (int j = 0; j < n_dims; j++)
            {
                y[j] = aux1 * y[j];
                y[j] += a[j];
                return_ort[i][j] = y[j];
            }
        }
    }

    for (int j = 0; j < n_dims; j++)
    {
        return_ort[indices[0]][j] = a[j];
        return_ort[indices[1]][j] = b[j];
    }
    free(x);
    free(y);
    return return_ort;
}

double **create_array_pts(int n_dims, long np)
{
    double **p_arr;
    p_arr = (double **)malloc(np * sizeof(double *));
    if (p_arr == NULL)
    {
        printf("Error allocating array of points, exiting.\n");
        exit(4);
    }
    for (int i = 0; i < np; i++)
        p_arr[i] = (double *)malloc(n_dims* sizeof(double));
        
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

int cmp(const void *a, const void *b)
{
    if (*(double *)a > *(double *)b)
        return 1;
    else if (*(double *)a < *(double *)b)
        return -1;
    else
        return 0;
}

double *find(double **orto_points, double a, long n_points)
{
    for (long i = 0; i < n_points; i++)
    {
        if (orto_points[i][0] == a)
        {
            return orto_points[i];
        }
    }
}

double *avgPoint(double *a, double *b, int n_dims)
{
    double *ret_vect = (double *)malloc(sizeof(double) * n_dims);
    for (int j = 0; j < n_dims; j++)
    {
        ret_vect[j] = (a[j] + b[j]) / 2;
    }
    return ret_vect;
}

double *center(double **orto_points, long n_points, int n_dims, int *indices)
{
    double *center = (double *)malloc(sizeof(double) * n_dims);
    double *aux = (double *)malloc(sizeof(double) * n_points);
    for (int i = 0; i < n_dims; i++){
        for (long j = 0; j < n_points; j++)
        {
            aux[j] = orto_points[j][i];
        }
        if (orto_points[0][i] != orto_points[1][i])
        {
            qsort(aux, n_points, sizeof(aux[0]), cmp);
            break;
        }
    }

    /*printf("Ordered:\n");
    for (long j = 0; j < n_points; j++)
    {
        printf("%f\n", aux[j]);
    }*/
    if (n_points % 2 != 0)
    {
        center = find(orto_points, aux[(int)n_points / 2], n_points);
    }
    else if (n_points % 2 == 0)
    {
        center = avgPoint(find(orto_points, aux[(int)n_points / 2], n_points), find(orto_points, aux[(int)(n_points / 2) - 1], n_points), n_dims);
    }
    free(aux);
    return center;
}

int build_tree(node *root, int id)
{
    if (root == NULL)
    {
        return -1;
    }
    if (root->n_points > 1)
    {
        // Work on this node
        int *indices;
        indices = furthest(root->pts, root->n_points, root->n_dims);
        int half = root->n_points / 2;
        double **left = (double **)malloc(sizeof(double *) * half);
        double **right = (double **)malloc(sizeof(double *) * (half + 1));

        //Create orthogonal projection
        double **orto_points;
        orto_points = ort_proj(root->pts, root->n_points, root->n_dims, indices);

        // Calcualte center and divide (verificar numero de pontos)
        root->center = center(orto_points, root->n_points, root->n_dims, indices);
        root->radius = calcRadius(root->pts, root->n_dims, indices, root->center);

        int l_id = 0, r_id = 0;
        for (long i = 0; i < root->n_points; i++)
        {
            if (orto_points[i][0] < root->center[0])
            {   
                left[l_id] = root->pts[i];
                l_id++;
            }
            else
            {
                right[r_id] = root->pts[i];
                r_id++;
            }
        }

        // Print ID of node
        //print (radius & center)
        /*
        printf("-------------------- \n");
        printf("Radius: %f \n", root->radius);
        printf("Center: ");
        for (int i = 0; i < root->n_dims; i++)
        {
            printf("%f ", root->center[i]);
        }
        printf("\nNPoints: %d\n", root->n_points);
        */

        // call left node
        int _id = id + 1;
        node *n_left = newNode(left, l_id, root->n_dims, _id);
        root->left = n_left;
        _id = build_tree(n_left, _id);

        // call right node
        _id += 1;
        node *n_right = newNode(right, r_id, root->n_dims, _id);
        root->right = n_right;
        _id = build_tree(n_right, _id);

        return _id;
    }
    else
    {
        root->radius = 0;
        root->center = root->pts[0];
        // Print ID of node
        //print (radius & center)
        /*
        printf("--------------------end \n");
        printf("Radius: %f \n", root->radius);
        printf("Center: ");
        for (int i = 0; i < root->n_dims; i++)
        {
            printf("%f ", root->center[i]);
        }
        printf("\nNPoints: %d\n", root->n_points);
        printf("Id: %d\n", id);
        */
    }
    return id;
}

void printNode(node *root, FILE* fp)
{
    int id_left;
    int id_right;

    if (root->left == NULL)
        id_left = -1;
    else
        id_left = root->left->id;
    if (root->right == NULL)
        id_right = -1;
    else
        id_right = root->right->id;

    if (root->left != NULL)
        printNode(root->left, fp);
    if (root->right != NULL)
        printNode(root->right, fp);

    ///////////////////////////////////////////////////////////////
    fprintf(fp, "%d %d %d %f ", root->id, id_left, id_right, root->radius);
    //printf("%d %d %d %f ", root->id, id_left, id_right, root->radius);
    ////////////////////////////////////////////////////////////////
    //printf("%d %d %d %f ", root->id, id_left, id_right, root->radius);
    for (int j = 0; j < root->n_dims; j++)
    {
        fprintf(fp, "%f ", root->center[j]);
        //printf( "%f ", root->center[j]);
    }
    fprintf(fp, "\n");
    //printf("\n");
    if(root->id == 0) puts("rute");
}

/*
Recebe a raíz da árvore e o número de nós da árvore
Chama a função printNode
*/
void printTree(node *root, int n_nodes)
{
    puts("Escrita de ficheiro");
    ///////////////////////////////////////////////////////////////
    FILE *fp = NULL;
    fp = fopen("exemplogrande.tree", "w");
    if(!fp) perror("fopen");
    fprintf(fp, "%d %d\n", root->n_dims, n_nodes);
    ////////////////////////////////////////////////////////////////

    //printf("%d %d\n", root->n_dims, n_nodes);
    printNode(root, fp);
    fclose(fp);
}

int main(int argc, char *argv[])
{
    int n_dims = 0;
    long n_points = 0;
    double **pts;
    int max_id = 0;

    double exec_time;
    exec_time = -omp_get_wtime();
    pts = get_points(argc, argv, &n_dims, &n_points);

    node *root = newNode(pts, n_points, n_dims, 0);

    //build tree
    max_id=build_tree(root, 0);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1lf s\n", exec_time);

    
    if(root != NULL) printTree(root, (max_id+1));

    //dump tree
    //dump_tree(root);
}