#include <cstdio>
#include <vector>
#include <assert.h>
#include <sys/time.h>
#include <omp.h>
using namespace std;
int N;
struct Node
{
    int parent;
    vector<int>childrens; 
    double upper, lower, diagonal, right_side_hand;
}Nodes[6000000];

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

int main(int argc, char ** argv)
{
    if(argc != 3){
        printf("The input format must be: ./serial ../data/case.txt ../sresult/res.txt\n");
        return 0;
    }

    FILE *fp_data = fopen(argv[1], "r");
    if(fp_data == NULL){
        printf("The case file %s does not exists!", argv[1]);
        return 0;
    }

    fscanf(fp_data, "%d", &N);
    for(int i = 0; i < N; i++){
        int index;
        fscanf(fp_data, "%d %lf %lf %lf %lf %d", &index, &Nodes[i].upper, &Nodes[i].lower, &Nodes[i].right_side_hand, &Nodes[i].diagonal, &Nodes[i].parent);
        assert(index == i);
        if(i > 0){
            assert(Nodes[i].parent >= 0 && Nodes[i].parent < i);
            Nodes[Nodes[i].parent].childrens.push_back(index);
        }
    }
    fclose(fp_data);
    
    /*
    printf("%d\n", N);
    for(int i = 0; i < N; i++){
        printf("%d %f %f %f %f %d\n", i, Nodes[i].upper, Nodes[i].lower, Nodes[i].diagonal, Nodes[i].right_side_hand, Nodes[i].parent);
    }*/
    
    double start, end, time_for_serial=0;

    start = omp_get_wtime();
    hines_serial(Nodes, N);
    end = omp_get_wtime();
    time_for_serial = end - start;
    printf("The time use for serial Hines algorithm is %.10f s\n", time_for_serial);

    FILE *fp_output = fopen(argv[2], "w");

    fprintf(fp_output, "%d\n", N);
    for(int i = 0; i < N; i++){
        fprintf(fp_output, "%d %lf %lf %lf %lf\n", i, Nodes[i].upper, Nodes[i].lower, Nodes[i].right_side_hand, Nodes[i].diagonal);
    }
    fclose(fp_output);

}