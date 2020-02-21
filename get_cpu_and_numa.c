#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <numa.h>

int sched_getcpu(void) ;

int main(){
 int cpu = sched_getcpu();
 int node = numa_node_of_cpu(cpu);
 printf("cpu %d, node %d \n",cpu,node);
 return 0;
}
