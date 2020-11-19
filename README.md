# cuHinesBatch

### 问题理解

算法抽象后即为一个串行的backward+forward的深度遍历。在backward过程中父节点依赖子节点的值进行更新，forward过程中子节点依赖父节点的值进行更新。我们很难挖掘父子节点之间的并行性，关键就在于如何高效地并行不同子树的遍历。




