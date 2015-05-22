/*Sample program to start ppserver on juropa*/
#include <stdio.h>
#include <mpi.h>
#include <string.h>

int main(int argc, char **argv) {
	int port, port_step, PP_Servers_Per_Node, PROCS_PER_NODE, rank, length;
	char hostname[256];

	port = atoi(argv[1]);
	PP_Servers_Per_Node = atoi(argv[2]);
	PROCS_PER_NODE = atoi(argv[3]);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(hostname, &length);

	//start ppserver only per node
	//printf("%d %s\n",rank,hostname);
	//printf("%d \n",rank%PROCS_PER_NODE);

	char cmd[256];
	char buffer[256];
	strcpy(cmd, "ppserver.py -p ");
	port_step = port + rank;
	sprintf(buffer, "%d", port_step);
	strcat(cmd, buffer);
	sprintf(buffer, "%d", 1);
	strcat(cmd, " -w ");
	strcat(cmd, buffer);
	strcat(cmd, " -t 120 -k 9999999 -s '123456' &");
	//printf("%s\n",cmd);
	if(rank%PROCS_PER_NODE < PP_Servers_Per_Node){
		//printf("%s\n",cmd);
		printf("Starting ppserver on %s %d \n", hostname, port_step);
		system(cmd);
		fflush(stdout);
	}

	MPI_Finalize();
	return 0;
}
