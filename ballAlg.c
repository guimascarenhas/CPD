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

node* newNode(double **arr,long n_points,int n_dims)
{   
    node* aux = (node*)malloc(sizeof(node));
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

double comp_dist(double *a, double *b , int *n_dims){
    double aux, sum = 0, two = 2;
    for(int i = 0; i < *n_dims; i++){
        aux = a[i] - b[i];
        aux = pow(aux, two);
        sum += aux;
    }
    return sqrt(sum); 
}

double calcRadius(double **pt_arr, int n_dims, int *indices, double *center){
    double dist_a = comp_dist(pt_arr[indices[0]], center, &n_dims);
    double dist_b = comp_dist(pt_arr[indices[1]], center, &n_dims);

    if(dist_a > dist_b)
        return dist_a;
    else return dist_b; 
}

int* furthest(double **pt_arr,long n_points, int n_dims){
    double *aux = pt_arr[0];
    double dist=0, dist_max=0;
    int ind_a=0,ind_b=0;
    for (long i=1; i< n_points; i++){
        dist = comp_dist(aux,pt_arr[i],&n_dims);
        if(dist>dist_max){
            dist_max = dist;
            ind_b = i;
        }
    }
    dist_max = 0;
    for (long i=0; i< n_points; i++){
        if( i!= ind_b){
            dist = comp_dist(pt_arr[ind_b],pt_arr[i],&n_dims);
            if(dist>dist_max){
                dist_max = dist;
                ind_a = i;
            }
        }
    }
    static int return_aux[2];
    return_aux[0] = ind_a; return_aux[1] = ind_b;
    return return_aux;
}

double inner_product(double *x,double *y, int n_dims){
    double result = 0;
    int i;
    for (i = 0; i < n_dims; i++)
    {
        result += x[i]*y[i];
    }
    return result;
}

double **ort_proj(double **pt_arr,long n_points, int n_dims, int *indices){
    double *a = pt_arr[indices[0]];
    double *b = pt_arr[indices[1]];
    double aux1,aux2;
    double **return_ort = (double **)malloc( n_points * sizeof(double*));
    for (long i=0; i< n_points; i++){
        return_ort[i] = (double*)malloc(sizeof(double)*n_dims);
    }

    double* x = (double*)malloc(sizeof(double)*n_dims);
    double* y = (double*)malloc(sizeof(double)*n_dims);

    for (long i=0; i< n_points; i++){
        if( i!= indices[0] && i!= indices[1])
        {
            for (int j=0; j< n_dims; j++)
            {
                x[j] = pt_arr[i][j] - a[j];
                y[j] = b[j] - a[j];
                //printf("%d dims\n", n_dims);
            }
            //printf("%f %f \n", x[0], x[1]);
            //printf("%f %f \n", y[0], y[1]);
            aux1 = inner_product(x,y,n_dims);
            aux2 =  inner_product(y,y,n_dims);
            aux1 = aux1/aux2;
            for (int j=0; j< n_dims; j++)
            {
                y[j] = aux1 * y[j];
                y[j] += a[j];
                return_ort[i][j] = y[j];
            }
        }
    }
    for (int j=0; j< n_dims; j++)
        {
            return_ort[indices[0]][j] = a[j];
            return_ort[indices[1]][j] = b[j];
        }
    free(x);
    free(y);
    return return_ort;
}

double **create_array_pts(int n_dims, long np){
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

double **get_points(int argc, char *argv[], int *n_dims, long *np){
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

/*int cmp (const void * a, const void * b) {
    printf("%f - %f = %f\n",*(double*)a,*(double*)b, *(double*)a - *(double*)b);
   return ( *(double*)a - *(double*)b );
}*/

int cmp (const void * a, const void * b) {
    if (*(double*)a > *(double*)b)
        return 1;
    else if (*(double*)a < *(double*)b)
        return -1;
    else
        return 0;
}

double *find(double **orto_points,double a,long n_points){
    for(long i=0; i< n_points;i++){
        if(orto_points[i][0] == a){
            return orto_points[i];
        }
    }
}

double* avgPoint(double* a, double* b, int n_dims){
    double* ret_vect = (double*)malloc(sizeof(double)*n_dims);
    for (int j=0; j< n_dims; j++)
    {
        ret_vect[j] = (a[j] + b[j])/2;
    }
    return ret_vect;
}

double* center(double **orto_points, long n_points, int n_dims, int *indices){
    double *center = (double*)malloc(sizeof(double)*n_dims);
    double *aux = (double*)malloc(sizeof(double)*n_points);
    for (long j=0; j< n_points; j++)
    {
        aux[j] = orto_points[j][0];
    }
    if(orto_points[0][0] != orto_points[1][0]){
        qsort(aux,n_points,sizeof(aux[0]),cmp);
    }
    printf("Ordered:\n");
    for (long j=0; j< n_points; j++)
    {
        printf("%f\n", aux[j]);
    }
    if(n_points % 2 != 0){
        center = find(orto_points, aux[(int)n_points/2], n_points);
    }
    else if(n_points % 2 == 0){
        center = avgPoint( find(orto_points, aux[(int)n_points/2],n_points) , find(orto_points, aux[(int) (n_points/2) - 1],n_points) , n_dims);
    }
    free(aux);
    return center;
}


void build_tree(node* root){   
    if( root == NULL){
        return;
    }
    if (root->n_points > 1){
        // Work on this node
        int *indices;
        indices = furthest(root->pts,root->n_points,root->n_dims);
        int half = root->n_points/2;
        double **left = (double**)malloc(sizeof(double*) * half);
        double **right = (double**)malloc(sizeof(double*) * (half + 1));
            
        //Create orthogonal projection     
        double** orto_points;
        orto_points = ort_proj(root->pts,root->n_points,root->n_dims, indices);
        if (root->n_points == 3){
            for (int i = 0; i < root->n_points; i++)
            {
                for (int j = 0; j < root->n_dims; j++)
                    printf("%f ", orto_points[i][j]);
                printf("\n");
            }
        }

        // Calcualte center and divide (verificar numero de pontos)
        root->center = center(orto_points,root->n_points,root->n_dims, indices);
        root->radius = calcRadius(root->pts, root->n_dims, indices, root->center);
        
        int l_id = 0, r_id = 0;
        for(long i=0;i<root->n_points;i++){
            if(orto_points[i][0] < root->center[0]){
                left[l_id] = root->pts[i];
                l_id ++;
            }
            else{
                right[r_id] = root->pts[i];
                r_id ++;
            }
        }
        /*printf("Left\n");
        for(long i=0; i<l_id;i++){
            for (int j = 0; j < root->n_dims; j++)
                printf("%f ", left[i][j]);
            printf("\n");
        }

        printf("Right\n");
        for(long i=0; i<r_id;i++){
            for (int j = 0; j < root->n_dims; j++)
                printf("%f ", right[i][j]);
            printf("\n");
        }*/

        // Print ID of node
        //print (radius & center)
        printf("-------------------- \n");
        printf("Radius: %f \n", root->radius);
        printf("Center: ");
        for(int i=0; i<root->n_dims;i++){
            printf("%f ", root->center[i]);
        }
        printf("\nNPoints: %d\n", root->n_points);

        // call left node
        node* n_left = newNode(left,l_id,root->n_dims);
        root->left = n_left;
        build_tree(n_left);

        // call right node
        node* n_right = newNode(right,r_id,root->n_dims);
        root->right = n_right;
        build_tree(n_right);
    }
    else{
        root->radius = 0;
        root->center = root->pts[0];
        // Print ID of node
        //print (radius & center)
        printf("--------------------end \n");
        printf("Radius: %f \n", root->radius);
        printf("Center: ");
        for(int i=0; i<root->n_dims;i++){
            printf("%f ", root->center[i]);
        }
        printf("\nNPoints: %d\n", root->n_points);
    }
    return;
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

    //furthest(pts,n_points,n_dims);

    node* root = newNode(pts,n_points,n_dims);
    build_tree(root);
    //build tree

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1lf\n", exec_time);

    //dump tree
}