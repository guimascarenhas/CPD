#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <math.h>

#define RANGE 10

int depht;
int last_layer_splits;
int machine_id;
double **my_pts;
long my_n_points;
int sub_tree_id;
int* machines;
int flag;


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


/*
Como o nome indica retorna um vetor que é a copia do argumento
*/
double* copy_vector(double* vector, int n_dims){
    double *aux = (double *)malloc(n_dims * sizeof(double));

    for(int i=0; i < n_dims; i++){
        aux[i] = vector[i];
    }
    return aux;
}

void printVector(double **vector, int n_points, int n_dims){
    for(int i=0; i<n_points; i++){
        for(int j=0; j<n_dims; j++){
            printf("%f ", vector[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

double log_base_2(double num){

    double aux = log(num);
    double aux2 = log(2);

    return (aux/aux2);
}

/*
Função recursiva que dá "free" nos 
vários nós da àrvore binaria
*/
void dump_tree(node *root){
    if(root->left != NULL)
        dump_tree(root->left);

    if(root->right != NULL)
        dump_tree(root->right);

    free(root->center);
    free(root);
    return;
}

/*
Recebe dois pontos e o número de dimensões destes pontos
Calcula e retorna a distancia entre eles
*/
double comp_dist(double *a, double *b, int n_dims)
{
    double aux, sum = 0, two = 2;
    for (int i = 0; i < n_dims; i++)
    {
        aux = a[i] - b[i];
        aux = pow(aux, two);
        sum += aux;
    }
    return sqrt(sum);
}

/*
Rcebe os pontos do cluster
Chama a função que computa a distância do centro a cada um desses pontos.
Retorna a maior das distâncias
*/
double calcRadius(double **pt_arr, int n_dims, double *center, int n_points)
{
    double dist_max = 0, dist = 0;
    for(int i = 0; i < n_points; i++){
        dist = comp_dist(pt_arr[i], center, n_dims);
        if(dist > dist_max){
            dist_max = dist;
        }
    }
    return dist_max;
}

/*
Rcebe coordenadas do centro e um conjunto de pontos
Chama a função que computa a distância do centro a cada um dos pontos.
Retorna a maior das distâncias (raio)
*/
int *furthest(double **pt_arr, long n_points, int n_dims)
{
    double *aux = pt_arr[0];
    double dist = 0, dist_max = 0;
    int ind_a = 0, ind_b = 0;
    for (long i = 1; i < n_points; i++)
    {
        dist = comp_dist(aux, pt_arr[i], n_dims);
        if (dist > dist_max)
        {
            dist_max = dist;
            ind_b = i;
        }
    }
    dist_max = 0;
    for (long i = 0; i < n_points; i++)
    {
        if (i != ind_b)
        {
            dist = comp_dist(pt_arr[ind_b], pt_arr[i], n_dims);
            if (dist > dist_max)
            {
                dist_max = dist;
                ind_a = i;
            }
        }
    }
    int* return_aux = (int*)malloc(2*sizeof(int));
    return_aux[0] = ind_a;
    return_aux[1] = ind_b;

    return return_aux;
}

/*
Recebe dois vetores (e o nº de dimensões dos vetores)
Calacula e retorna o produto interno dos dois
*/
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

/*
Aloca um arrar de pontos com base no numero total de pontos bem como o numero de coordenadas por pontos
*/
double **create_array_pts(int n_dims, long np)
{
    double *_p_arr;
    double **p_arr;

    _p_arr = (double *) malloc(n_dims * np * sizeof(double));
    p_arr = (double **) malloc(np * sizeof(double *));
    if((_p_arr == NULL) || (p_arr == NULL)){
        printf("Error allocating array of points, exiting.\n");
        exit(4);
    }

    for(long i = 0; i < np; i++)
        p_arr[i] = &_p_arr[i * n_dims];

    return p_arr;
}

/*
Rcebe o cluster de pontos e os indices dos pontos que irão formar a reta para o qual todos os pontos se irão projetar.
Computa a projeção na reta de todos os pontos do cluster.
Retorna as coordenadas das projeções de todos os pontos do cluster na reta
*/
double **ort_proj(double **pt_arr, long n_points, int n_dims, int *indices)
{
    double *a = pt_arr[indices[0]];
    double *b = pt_arr[indices[1]];
    double aux1, aux2;
    
    double **return_ort = create_array_pts(n_dims, n_points);

    double *x = (double *)malloc(sizeof(double) * n_dims);
    double *y = (double *)malloc(sizeof(double) * n_dims);

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
/*     if (*n_dims < 2)
    {
        printf("Illegal number of dimensions (%d), must be above 1.\n", *n_dims);
        exit(2);
    }  */
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

/*
Rcebe o cluster de pontos nas suas projecoes ortogonais.
Retorna as coordenadas do ponto central das projecoes ortogonais.
*/
double *center(double **orto_points, long n_points, int n_dims, int *indices)
{
    double *center;
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

    if (n_points % 2 != 0)
    {
        center = copy_vector(find(orto_points, aux[(int)n_points / 2], n_points), n_dims);
    }
    else if (n_points % 2 == 0)
    {
        center = avgPoint(find(orto_points, aux[(int)n_points / 2], n_points), find(orto_points, aux[(int)(n_points / 2) - 1], n_points), n_dims);
    }
    free(aux);
    return center;
}

/*
Funcao recursiva que recebendo um cluster (root) constroi a restante arvore
*/
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

        // Calcualte center and divide 
        root->center = center(orto_points, root->n_points, root->n_dims, indices);
        root->radius = calcRadius(root->pts, root->n_dims, root->center, root->n_points);

        free(indices);

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

        free(orto_points[0]);
        free(orto_points);

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

        free(left);
        free(right);
        return _id;
    }
    else
    {
        root->radius = 0;
        root->center = copy_vector(root->pts[0], root->n_dims);
    }
    return id;
}



void sendPoints(double** pts, int id, long n_points, int n_dims, int dest_id){

    MPI_Send(&n_points, 1, MPI_LONG, id, 0, MPI_COMM_WORLD);
    MPI_Send(&dest_id, 1, MPI_INT, id, 1, MPI_COMM_WORLD);

    for(long i=0; i < n_points; i++){
            MPI_Send(pts[i], n_dims, MPI_DOUBLE, id, 2, MPI_COMM_WORLD);
    }
}

double** rcvPoints(int src, int n_dims, long *n_points, int *my_id){
    MPI_Status st;
    MPI_Recv(n_points, 1, MPI_LONG, src, 0, MPI_COMM_WORLD, &st);
    MPI_Recv(my_id, 1, MPI_INT, src, 1, MPI_COMM_WORLD, &st);

    double **pts = create_array_pts(n_dims, *n_points);

    for(long i=0; i < *n_points; i++){
            MPI_Recv(pts[i], n_dims, MPI_DOUBLE, st.MPI_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return pts;
}


void printNode(node *root)
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
        printNode(root->left);
    if (root->right != NULL)
        printNode(root->right);

    printf("%d %d %d %f ", root->id, id_left, id_right, root->radius);

    for (int j = 0; j < root->n_dims; j++)
    {
        printf("%f ", root->center[j]);
    }
    printf("\n");

}

int build_tree_iter(node *root, int id, int layer)
{   
    //printf("ID: %d\n", id);
    if (root == NULL)
    {
        return -1;
    }
    if(layer > depht){
        //FAZ COISAS
        //printf("Sent to machine %d: \n", machine_id);
        machines[machine_id] = id;
        if(machine_id != 0){
            sendPoints(root->pts, machine_id, root->n_points, root->n_dims, id);
        }
        else{
            my_pts=create_array_pts(root->n_dims, root->n_points);
            for(int i=0; i<root->n_points; i++){
                for(int j=0; j<root->n_dims; j++){
                    my_pts[i][j]=root->pts[i][j];
                }
            }
            my_n_points = root->n_points;
            sub_tree_id = id;
        }
        machine_id++;
        return id;
    }
    if(layer == depht){
        if(last_layer_splits <= 0){
        //FAZ COISAS
            //printf("Sent to machine %d: \n", machine_id);
            machines[machine_id] = id;
            if(machine_id != 0){
                sendPoints(root->pts, machine_id, root->n_points, root->n_dims, id);
            }
            else{
                my_pts=create_array_pts(root->n_dims, root->n_points);
                for(int i=0; i<root->n_points; i++){
                    for(int j=0; j<root->n_dims; j++){
                        my_pts[i][j]=root->pts[i][j];
                    }
                }
                my_n_points = root->n_points;
                sub_tree_id = id;
            }
            //printVector(root->pts, root->n_points);
            machine_id++;
            return id;
        }
        else{
            last_layer_splits--;
        }
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

        // Calcualte center and divide 
        root->center = center(orto_points, root->n_points, root->n_dims, indices);
        root->radius = calcRadius(root->pts, root->n_dims, root->center, root->n_points);

        free(indices);

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

        free(orto_points[0]);
        free(orto_points);


        // call left node
        node *n_left = newNode(left, l_id, root->n_dims, (id+1));
        root->left = n_left;
        build_tree_iter(n_left, (id+1), (layer+1));

        // call right node
        node *n_right = newNode(right, r_id, root->n_dims,(id + 2*l_id));
        root->right = n_right;
        build_tree_iter(n_right, (id + 2*l_id), (layer+1));

        free(left);
        free(right);
        return id;
    }
    else
    {
        root->radius = 0;
        root->center = copy_vector(root->pts[0], root->n_dims);
    }
    return id;
}

void treeTraversal(node* root, node* new_root){

    int thrash = 0;
    if(machines[machine_id] == root->id){
        if(machine_id == 0){
            printNode(new_root);
        }
        else{
            MPI_Send(&thrash, 1, MPI_INT, machine_id, 3, MPI_COMM_WORLD);
            MPI_Recv(&thrash, 1, MPI_INT, machine_id, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }
        machine_id++;
    }
    else{

        if (root->left != NULL)
            treeTraversal(root->left, new_root);
        if (root->right != NULL)
            treeTraversal(root->right, new_root);
    
        printf("%d %d %d %f ", root->id, root->left->id, root->right->id, root->radius);
        for (int j = 0; j < root->n_dims; j++)
        {
            printf("%f ", root->center[j]);
        }
        printf("\n");

    }

}

/*
Recebe a raíz da árvore e o número de nós da árvore
Chama a função printNode
*/
void printTree(node *root, node* new_root)
{
    printf("%d %d\n", root->n_dims, 2*root->n_points - 1);

    treeTraversal(root, new_root);

}

void recvPrintOrder(node* root){

    int thrash = 0;
    MPI_Recv(&thrash, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printNode(root);
    MPI_Send(&root->id, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);
}



int main(int argc, char *argv[])
{
    int n_dims = 0;
    long n_points;
    double **pts;
    int max_id = 0;
    double exec_time;
    int id, p;
    int my_id = 0;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    machine_id = 0;
    //count = p;
    depht = floor((log_base_2(p)+1));
    last_layer_splits = p - pow(2, floor(log_base_2(p)));

    node *root = NULL;
    node *new_root = NULL;

    if(id == 0){
        //printf("depht: %d\n", depht);
        //printf("llsplits: %d\n", last_layer_splits);
        exec_time = -omp_get_wtime();
        pts = get_points(argc, argv, &n_dims, &n_points);

        root = newNode(pts, n_points, n_dims, 0);

        machines = (int*)malloc(p * sizeof(int));
        build_tree_iter(root, 0, 1);
        //printTree(root, root->n_points);

        new_root = newNode(my_pts, my_n_points, n_dims, sub_tree_id);
/*       printf("my_points: %d, n_dims: %d\n", my_n_points, n_dims);
        printVector(new_root->pts, my_n_points, n_dims); */
        build_tree(new_root, sub_tree_id);
        //printNode(root);
    }
    else{
        n_dims = atoi(argv[1]);
        pts = rcvPoints(0, n_dims, &n_points, &my_id);
        //printf("Machine %d has id %d\n", id, my_id);
        root = newNode(pts, n_points, n_dims, my_id);
        build_tree(root, my_id);
    }

    //printf("Machine %d is ready\n", id);
    MPI_Barrier(MPI_COMM_WORLD);

    if(id == 0){
        machine_id = 0;
        /* printf("machines: ");
        for(int j =0; j<p; j++){
            printf("%d ", machines[j]);
        }
        printf("\n"); */
        printTree(root, new_root);
    }
    else{
        recvPrintOrder(root);
    }


    MPI_Barrier(MPI_COMM_WORLD);

    if(id == 0){
        exec_time += omp_get_wtime();
        fprintf(stderr, "%.1lf\n", exec_time);
    }

    free(pts[0]);
    free(pts);


    MPI_Finalize();
    //printTree(root, (max_id+1));

    //dump_tree(root);
    //if(root != NULL) dump_tree(root);

    return 0;
   
}