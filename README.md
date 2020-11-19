# cuHinesBatch

### Problem Understanding
A simulation of the behavior of the Human Brain, to solve efficiently a high number of Hines systems(neurons)

It is actually a serial backward + forward depth traversal. In the backward process, the parent node updates according to the value of the child node, and the child node updates according to the value of the parent node in the forward process. It is difficult to mine the parallelism between parent and child nodes. The key is how to efficiently parallel traversal of different subtrees.


