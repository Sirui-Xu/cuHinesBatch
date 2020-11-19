### 实现思路

并行算法采用openmp实现。

- 不同子树之间并行执行深度遍历

  在初始化时保存了每个节点的子节点标号，遍历时并行循环子节点。

  1. backward阶段

     我们需要对子节点遍历、更新后才能对父节点进行更新。由于每个子节点都会让父节点进行更新，为了解决父节点更新时的冲突问题，利用到了**reduction**操作加速。

     ```c++
     double d = Nodes[start].diagonal, rhs = Nodes[start].right_side_hand;
         #pragma omp parallel for reduction(-:d, rhs) 
         for(int i = 0; i < Nodes[start].childrens.size(); ++i){
             double factor = 0;
             int index;
             index = Nodes[start].childrens[i];
             hines_backward_parallel(Nodes, index, Nodes[index].end);
             factor = Nodes[index].upper / Nodes[index].diagonal;
             d -= factor * Nodes[index].lower;
             rhs -= factor * Nodes[index].right_side_hand;
         }
         Nodes[start].diagonal = d;
         Nodes[start].right_side_hand = rhs;
     ```

  2. forward阶段

     我们对子节点的遍历是在更新完子节点后进行的，此时不需要等待子树遍历完成。因此，我们利用了**task**操作，使得子线程继续遍历，而父线程无需等待。由于是每个子节点依赖父节点的更新，子节点的更新并不存在冲突问题，因此不需要reduction操作。

     在遍历后需要加上**taskwait**回收线程。

     ```c++
     #pragma omp parallel for
         for(int i = 0; i < Nodes[start].childrens.size(); ++i){
             int index;
             index = Nodes[start].childrens[i];
             Nodes[index].right_side_hand -= Nodes[index].lower * Nodes[start].right_side_hand;
             Nodes[index].right_side_hand /= Nodes[index].diagonal;
             #pragma omp task
             hines_forward_parallel(Nodes, index, Nodes[index].end);
         }
     ```

- 设计串行循环函数、并行递归函数相互调用的递归模式

  观察到数据中存在着大量的仅有一条通路、没有其它分支的子段，如果一直递归遍历，那么函数的调用次数会非常多，代价也会很大。

  因此我们设计，遍历到某节点时，如果该节点的子节点个数较少，并行递归函数跳转到串行循环函数。而如果在串行循环函数内某节点的子节点个数较多，又跳转回并行的递归函数。

  1. backward阶段，由于是从后往前执行，串行函数必须得先判断终止节点的位置，确定串行的节点。实验中采用的子节点个数阈值为1，核心代码部分片段如下：

     ```c++
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
         for(int j = i; j > start; j--){
             factor = Nodes[j].upper / Nodes[j].diagonal;
             parent = Nodes[j].parent;
             Nodes[parent].diagonal -= factor * Nodes[j].lower;
             Nodes[parent].right_side_hand -= factor * Nodes[j].right_side_hand;
         }
     }
     
     void hines_backward_parallel(Node Nodes[], int start, int end){
         if(Nodes[start].childrens.size() == 1){
             hines_backward_serial(Nodes, start, end);
             return;
         }
       	...
     ```

  2. forward阶段，无需提前确定串行的节点，在计算过程中判断即可。实验中采用的子节点个数阈值为1，核心代码部分片段如下：

     ```c++
     void hines_forward_serial(Node Nodes[], int start, int end){ 
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
     void hines_forward_parallel(Node Nodes[], int start, int end){ 
         if(Nodes[start].childrens.size() == 1){
             hines_forward_serial(Nodes, start, end);
             return;
         }
         ...
     ```

- 设计可并行性评价指标

  设计了两个可并行性评价指标：

  1. thres1 = 子树树高 / 子树节点个数

     如果子树的树高相比子树的节点个数较少，那么这颗子树的分支可能的分叉就较多，采用并行执行的效率可能就会越高。反之，如果子树的树高与子树节点个数非常的接近，那么这颗子树调用并行的函数执行的效率必然不高。

     在实验中测得 thres1 = 0.95 表现较好

  2. thres2 = 子树节点个数

     如果子树节点个数较少，那么多开线程的代价相比计算的代价就高，此时不妨直接用串行循环来计算。

     在实验中测得令thres2 = 50000 表现较好

  子树的树高需要在初始化建树时利用动态规划算法提前计算得到，子树的节点个数需要在初始化建树时深度遍历提前计算得到。

  

### 测试性能

40核服务器上测试，串行与并行部分分别重复150次，前50次作为warmup，统计后100次的平均时间，计算speedup。

- 总体运行情况：

  | Case              | Case9  | Case10 | Case11 | Case12 |
  | ----------------- | ------ | ------ | ------ | ------ |
  | serial time / s   | 0.019  | 0.030  | 0.085  | 0.17   |
  | parallel time / s | 0.0014 | 0.0030 | 0.014  | 0.031  |
  | Speedup           | 13.6   | 10.0   | 6.1    | 5.6    |

- 对比实验（每次测试误差影响较大，因此仅选取了部分差异明显的实验组）

  1. 算法各项设计的有效性

  | Case10           | parallel time / s | speedup |
  | ---------------- | ----------------- | ------- |
  | 原始设置         | 0.0030            | 10.0    |
  | 没有串并相互调用 | 0.0033            | 8.8     |
  | 没有thres制约    | 0.0041            | 7.15    |

  2. warmup的影响

  | Case11          | parallel time / s | speedup |
  | --------------- | ----------------- | ------- |
  | warmup iter=200 | 0.014             | 6.1     |
  | warmup iter=20  | 0.016             | 5.3     |

  3. 超参数选择（仅展示部分结果）

  | Case12      | parallel time / s | speedup |
  | ----------- | ----------------- | ------- |
  | thres1=1.00 | 0.035             | 4.9     |
  | thres1=0.95 | 0.031             | 5.6     |

  4. 在小规模数据集下的有效性

  | Case5         | parallel time / s | speedup |
  | ------------- | ----------------- | ------- |
  | 原始设置      | 0.000038          | 1.01    |
  | 没有thres制约 | 0.00027           | 0.16    |

总结一下：

1. 该算法能够在最后四个大规模数据上**取得较好的结果**。

2. 从对比实验中可以看出，**串行、并行函数相互调用减轻了并行函数不断递归以及开新线程的代价；**

   通过可并行性评价指标的限制，使得一些无需并行的部分直接进入串行的函数实现，减少了一些无用的消耗。特别是对于前8个较小的case，没有限制的话会使得并行程序的速度急剧下降，**而现在我们的实现仍能在小规模数据（前8个case）上保证较好的速度。**

   **warmup对于并行程序也有一定的影响**，略微延长warmup的时间能够提升一些性能，而如果没有warmup，最开始几次的程序性能有所下降。

### 可改进的部分

1. 对不同的子节点分配不同的线程数，目前只是采用static的方式简单分配，利用dynamic schedule等等目前无法取得更好的效果，或许需要人工进行判断分配。
2. 挖掘同一路径上节点的并行性？

### 编译命令

g++ -std=c++11 -fopenmp -lm serial.cc -o serial

g++ -std=c++11 -fopenmp -lm parallel.cc -o parallel



