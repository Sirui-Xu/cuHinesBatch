#include <cstdio>
#include <vector>
#include <stack>
#include <assert.h>
#include <sys/time.h>
#include <omp.h>
#include <cstring>
#include <cmath>
using namespace std;
const double thres_for_parallel = 0.95;
const double thres_for_parallel_2 = 50000;
int num_procs = omp_get_num_procs();
int N;
struct Node
{
    int parent = -1;
    int end = N;
    int height = 1;
    vector<int>childrens; 
    double upper, lower, diagonal, right_side_hand;
}Nodes[6*1000*1000], Nodes_check[6*1000*1000], Nodes_copy[6*1000*1000];

void hines_backward_parallel(Node Nodes[], int start, int end);
void hines_forward_parallel(Node Nodes[], int start, int end);

bool isequal(double a, double b){
    return fabs(a - b) < 1e-5;
}

void init(FILE *fp_data){
    stack<int> parents;
    for(int i = 0; i < N; i++){
        int index;
        fscanf(fp_data, "%d %lf %lf %lf %lf %d", &index, &Nodes[i].upper, &Nodes[i].lower, &Nodes[i].right_side_hand, &Nodes[i].diagonal, &Nodes[i].parent);
        // assert(index == i);
        // assert(Nodes[i].parent >= 0 && Nodes[i].parent < i);
        Nodes[i].childrens.clear();
        if(i == 0){
            parents.push(index);
            continue;
        }else{
            Nodes[Nodes[i].parent].childrens.push_back(index);
        }
        int top = parents.top();
        while(Nodes[i].parent != top){
            Nodes[top].end = index;
            if(Nodes[top].childrens.size() == 0){
                Nodes[top].height = 1;
            }else{
                for(auto iter = Nodes[top].childrens.begin(); iter != Nodes[top].childrens.end(); ++iter){
                    Nodes[top].height = max(Nodes[top].height, Nodes[*iter].height);
                }
                Nodes[top].height += 1;
            }
            parents.pop();
            top = parents.top();
        }
        parents.push(index);
        
    }
    int top;
    while(parents.size() != 0){
        top = parents.top();
        Nodes[top].end = N;
        if(Nodes[top].childrens.size() == 0){
            Nodes[top].height = 1;
        }else{
            for(auto iter = Nodes[top].childrens.begin(); iter != Nodes[top].childrens.end(); ++iter){
                Nodes[top].height = max(Nodes[top].height, Nodes[*iter].height);
            }
            Nodes[top].height += 1;
        }
        parents.pop();
    }
}

void hines_serial(Node Nodes[], int N){
    double factor = 0;
    int parent = -1;
    for(int i = N - 1; i > 0; i--){
        factor = Nodes[i].upper / Nodes[i].diagonal;
        parent = Nodes[i].parent;
        Nodes[parent].diagonal -= factor * Nodes[i].lower;
        Nodes[parent].right_side_hand -= factor * Nodes[i].right_side_hand;
    }
    
    /*
    printf("%d\n", N);
    for(int i = 0; i < N; i++){
        printf("%d %f %f %f %f %d\n", i, Nodes[i].upper, Nodes[i].lower, Nodes[i].diagonal, Nodes[i].right_side_hand, Nodes[i].parent);
    }
    */
   
    Nodes[0].right_side_hand /= Nodes[0].diagonal;
    for(int i = 1; i < N; i++){
        parent = Nodes[i].parent;
        Nodes[i].right_side_hand -= Nodes[i].lower * Nodes[parent].right_side_hand;
        Nodes[i].right_side_hand /= Nodes[i].diagonal;
    }
}

void hines_backward(Node Nodes[], int start, int end){
    double factor = 0;
    int parent = -1;
    //#pragma omp parallel //num_threads(num_procs)
    //{
    //    #pragma omp single
    //    {
    for(int i = end - 1; i > start; i--){
        factor = Nodes[i].upper / Nodes[i].diagonal;
        parent = Nodes[i].parent;
        Nodes[parent].diagonal -= factor * Nodes[i].lower;
        Nodes[parent].right_side_hand -= factor * Nodes[i].right_side_hand;
    }
    //    }
    //}
}

void hines_backward_serial(Node Nodes[], int start, int end){
    double factor = 0;
    int parent = -1;
    int i;
    //判断返回节点的序号
    for(i = start + 1; i < end; i++){
        if(Nodes[i].childrens.size() > 1){
            hines_backward_parallel(Nodes, i, end);
            break;
        }
    }
    if(i == end) i = end - 1;
    for(int j = i; j > start; j--){
        factor = Nodes[j].upper / Nodes[j].diagonal;
        parent = Nodes[j].parent;
        Nodes[parent].diagonal -= factor * Nodes[j].lower;
        Nodes[parent].right_side_hand -= factor * Nodes[j].right_side_hand;
    }
}

void hines_forward(Node Nodes[], int start, int end){ // not for root
    int parent = -1;
    //#pragma omp parallel //num_threads(num_procs)
    //{
    //    #pragma omp single
    //    {
    for(int i = start + 1; i < end; i++){
        parent = Nodes[i].parent;
        Nodes[i].right_side_hand -= Nodes[i].lower * Nodes[parent].right_side_hand;
        Nodes[i].right_side_hand /= Nodes[i].diagonal;
    }
    //    }
    //}
}

void hines_forward_serial(Node Nodes[], int start, int end){ // not for root
    int parent = -1;
    int i;
    for(i = start + 1; i < end; i++){
        parent = Nodes[i].parent;
        Nodes[i].right_side_hand -= Nodes[i].lower * Nodes[parent].right_side_hand;
        Nodes[i].right_side_hand /= Nodes[i].diagonal;
        if(Nodes[i].childrens.size() > 1){
            hines_forward_parallel(Nodes, i, end);
            break;
        }
    }
}

void hines_backward_parallel(Node Nodes[], int start, int end){
    if(Nodes[start].height / (Nodes[start].end - start) > thres_for_parallel || (Nodes[start].end - start) < thres_for_parallel_2){
        hines_backward(Nodes, start, end);
        return;
    }
    if(Nodes[start].childrens.size() == 1){
        hines_backward_serial(Nodes, start, end);
        return;
    }
    double d = Nodes[start].diagonal, rhs = Nodes[start].right_side_hand;
    #pragma omp parallel for reduction(-:d, rhs) //num_threads(2*num_procs)//num_threads(Nodes[start].childrens.size())//schedule(dynamic) //num_threads(num_threads)
    for(int i = 0; i < Nodes[start].childrens.size(); ++i){
        double factor = 0;
        int index;
        index = Nodes[start].childrens[i];
        //#pragma omp task
        hines_backward_parallel(Nodes, index, Nodes[index].end);
        //#pragma omp taskwait
        factor = Nodes[index].upper / Nodes[index].diagonal;
        d -= factor * Nodes[index].lower;
        rhs -= factor * Nodes[index].right_side_hand;
    }
    Nodes[start].diagonal = d;
    Nodes[start].right_side_hand = rhs;
}

void hines_forward_parallel(Node Nodes[], int start, int end){ // not for root
    if(Nodes[start].height / (Nodes[start].end - start) > thres_for_parallel || (Nodes[start].end - start) < thres_for_parallel_2){
        hines_forward(Nodes, start, end);
        return;
    }
    if(Nodes[start].childrens.size() == 1){
        hines_forward_serial(Nodes, start, end);
        return;
    }
    #pragma omp parallel for //num_threads(2*num_procs)//num_threads(Nodes[start].childrens.size())//schedule(dynamic)//num_threads(num_threads)
    for(int i = 0; i < Nodes[start].childrens.size(); ++i){
        int index;
        index = Nodes[start].childrens[i];
        Nodes[index].right_side_hand -= Nodes[index].lower * Nodes[start].right_side_hand;
        Nodes[index].right_side_hand /= Nodes[index].diagonal;
        #pragma omp task
        hines_forward_parallel(Nodes, index, Nodes[index].end);
    }
}

void hines_parallel(Node Nodes[], int N){
    hines_backward_parallel(Nodes, 0, N);
    Nodes[0].right_side_hand /= Nodes[0].diagonal;
    hines_forward_parallel(Nodes, 0, N);
    #pragma omp taskwait
}

void my_memcpy(Node dst[], Node src[], int N){
    #pragma omp parallel for
    for(int i = 0; i < N; i++){
        dst[i].diagonal = src[i].diagonal;
        dst[i].end = src[i].end;
        dst[i].height = src[i].height;
        dst[i].lower = src[i].lower;
        dst[i].parent = src[i].parent;
        dst[i].right_side_hand = src[i].right_side_hand;
        dst[i].upper = src[i].upper;
        dst[i].childrens.clear();
        for(int j = 0; j < src[i].childrens.size(); j++){
            dst[i].childrens.push_back(src[i].childrens[j]);
        }
    }
}

int main(int argc, char ** argv)
{
    if(argc != 3){
        printf("The input format must be: ./serial ../data/case.txt ../presult/res.txt\n");
        return 0;
    }

    FILE *fp_data = fopen(argv[1], "r");
    if(fp_data == NULL){
        printf("The case file %s does not exists!", argv[1]);
        return 0;
    }

    //initialize
    printf("begin to initialize the tree structure...\n");
    fscanf(fp_data, "%d", &N);
    init(fp_data);
    printf("done!\n");
    fclose(fp_data);

    double start, end, time_for_parallel=0, time_for_serial=0;
    my_memcpy(Nodes_copy, Nodes, N);
    
    // serial part
    my_memcpy(Nodes_check, Nodes, N);
    //printf("The num_procs = %d\n", num_procs);

    for(int i = 0; i < 150; i++){
        start = omp_get_wtime();
        hines_serial(Nodes_check, N);
        end = omp_get_wtime();
        if(i >= 50) time_for_serial += end - start; // the first 50 iter for warmup, the latter 100 iters count.
        my_memcpy(Nodes_check, Nodes_copy, N);
    }
    printf("The time use for serial Hines algorithm is %.10f s\n", time_for_serial / 100);
    

    // parallel part
    for(int i = 0; i < 150; i++){
        start = omp_get_wtime();
        hines_parallel(Nodes, N);
        end = omp_get_wtime();
        if(i >= 50) time_for_parallel += end - start; // the first 50 iter for warmup, the latter 100 iters count.
        my_memcpy(Nodes, Nodes_copy, N);
    }

    printf("The time use for parallel Hines algorithm is %.10f s\n", time_for_parallel / 100);

    printf("The speedup is %.2f\n", time_for_serial / time_for_parallel);

    hines_parallel(Nodes, N);

    
    //check
    hines_serial(Nodes_check, N);

    for(int i = 0; i < N; i++){
        if(!(isequal(Nodes[i].upper, Nodes_check[i].upper)
        && isequal(Nodes[i].lower, Nodes_check[i].lower)
        && isequal(Nodes[i].diagonal, Nodes_check[i].diagonal)
        && isequal(Nodes[i].right_side_hand, Nodes_check[i].right_side_hand))){
            printf("parallel algorithm failed!\n");
            return 0;
        }
    }
    printf("correctness check pass!\n");
    

    //output
    FILE *fp_output = fopen(argv[2], "w");

    fprintf(fp_output, "%d\n", N);
    for(int i = 0; i < N; i++){
        fprintf(fp_output, "%d %lf %lf %lf %lf\n", i, Nodes[i].upper, Nodes[i].lower, Nodes[i].right_side_hand, Nodes[i].diagonal);
    }
    fclose(fp_output);
    printf("write back to results dir!\n");
    return 0;
}