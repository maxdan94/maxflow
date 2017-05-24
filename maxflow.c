/*
Maximilien Danisch
May 2017
http://bit.ly/maxdan94
maximilien.danisch@gmail.com

Info:
Feel free to use these lines as you wish. This is an efficient C implementation of the Ford-Fulkerson algorithm. It should easily scale to hundreds of millions of edges if the data is not too adverserial.

To compile:
"gcc maxflow.c -O3 -o maxflow".

To execute:
"./maxflow edgelist.txt source target res.txt".
"edgelist.txt" should contain the directed edges with capcities: one edge on each line separated by spaces "n1 n2 c".
"source" and "target" are the source and target ids.
"res.txt" contains the results: "n1 n2 c f" on each line.
Will print some information in the terminal.
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

#define NLINKS 1000000000

typedef struct {
	unsigned u;//first node
	unsigned v;//second node
	int c;//capacity of the edge
	int f;//flow from u to v
} edge;

typedef struct {
	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge *edges;//list of edges

	unsigned *d;//d[i]=degree of node i
	unsigned *cd;//cumulative degree: (start with 0) length=dim+1
	unsigned *adj;//list of neighbors
	unsigned *eid;//ID of the edges
} flowgraph;

void freefg(flowgraph *g){
	free(g->edges);
	free(g->cd);
	free(g->adj);
	free(g->eid);
	free(g);
}

//compute the maximum of three unsigned
inline unsigned max3(unsigned a,unsigned b,unsigned c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the edgelist with capacities from file
flowgraph* readedgelist(char* edgelist){
	unsigned e1=NLINKS;
	flowgraph *g=malloc(sizeof(flowgraph));
	FILE *file;

	g->n=0;
	g->e=0;
	file=fopen(edgelist,"r");
	g->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%u %u %u\n", &(g->edges[g->e].u), &(g->edges[g->e].v),&(g->edges[g->e].c))==3) {//)==2){g->edges[g->e].c=1;
		g->edges[g->e].f=0;
		g->n=max3(g->n,g->edges[g->e].u,g->edges[g->e].v);
		if (g->e++==e1) {
			e1+=NLINKS;
			g->edges=realloc(g->edges,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->n++;
	g->edges=realloc(g->edges,g->e*sizeof(edge));
	return g;
}

//Building a special sparse graph structure
void mkfg(flowgraph *g){
	unsigned i,u,v;
	unsigned *d=calloc(g->n,sizeof(unsigned));
	g->cd=malloc((g->n+1)*sizeof(unsigned));
	g->adj=malloc(2*(g->e)*sizeof(unsigned));
	g->eid=malloc(2*(g->e)*sizeof(unsigned));

	for (i=0;i<g->e;i++) {
		d[g->edges[i].u]++;
		d[g->edges[i].v]++;
	}
	g->cd[0]=0;
	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+d[i-1];
		d[i-1]=0;
	}
	for (i=0;i<g->e;i++) {
		u=g->edges[i].u;
		v=g->edges[i].v;
		g->eid[g->cd[u] + d[u] ]=i;
		g->eid[g->cd[v] + d[v] ]=i;
		g->adj[g->cd[u] + d[u]++ ]=v;
		g->adj[g->cd[v] + d[v]++ ]=u;
	}
	free(d);
}


void printres(flowgraph *g,char *output){
	unsigned i;
	edge ed;
	FILE* file=fopen(output,"w");
	for (i=0;i<g->e;i++) {
		ed=g->edges[i];
		fprintf(file,"%u %u %d %d\n",ed.u,ed.v,ed.c,ed.f);
	}
	fclose(file);
}

typedef struct {
	bool d;//distance: 0 iff not infinity
	unsigned p;//predecessor
	unsigned e;//edgeid
} recover;


//finds an "augmenting path" and changes the flow accordingly
unsigned onepath(flowgraph *g,unsigned s, unsigned t) {
	unsigned n1=1,n2=0,n=0,d=0,i,j,k,u,v;
	int f=0,tmp,min;
	edge ed;

	static recover *r=NULL;
	static unsigned *l1=NULL,*l2=NULL,*l3=NULL,*l=NULL;
	if (r==NULL){
		r=malloc(g->n*sizeof(recover));
		l1=malloc(g->n*sizeof(unsigned)),l2=malloc(g->n*sizeof(unsigned)),l=malloc(g->n*sizeof(unsigned));
		for (i=0;i<g->n;i++) {
			r[i].d=1;
		}
	}

	l1[0]=s;
	r[s].d=0;

	do {
		d++;
		for (i=0;i<n1;i++) {
			u=l1[i];
			for (j=g->cd[u];j<g->cd[u+1];j++) {
				v=g->adj[j];
				k=g->eid[j];
				ed=g->edges[k];

				if (r[v].d)
					if ((ed.u==u && ed.c!=ed.f) || (ed.u==v && ed.f!=0)){
					l2[n2++]=v;
					l[n++]=v;
					r[v].d=0;
					r[v].p=u;
					r[v].e=k;
				}
			}
		}
		n1=n2,n2=0;
		l3=l2,l2=l1,l1=l3;
	} while (n1 && r[t].d);

	if (r[t].d){
		return 0;
	}

	for (i=0;i<n;i++){//resetting distances to infinity for next run
		r[l[i]].d=1;
	}

	u=t;
	min=-1;
	while (u!=s){
		ed=g->edges[r[u].e];
		v=r[u].p;
		if (ed.u==v){
			tmp=ed.c-ed.f;
			min=(min<tmp)?min:tmp;
		}
		else{
			tmp=ed.c+ed.f;
			min=(min<tmp)?min:tmp;
		}
		u=v;
	}
	u=t;
	while (u!=s){
		v=r[u].p;
		if (g->edges[r[u].e].u==v){
			g->edges[r[u].e].f-=min;
		}
		else{
			g->edges[r[u].e].f+=min;
		}
		u=v;
	}

	for (i=g->cd[s];i<g->cd[s+1];i++) {
		u=g->adj[i];
		k=g->eid[i];
		ed=g->edges[k];
		if (ed.u==s){
			f+=ed.f;
		}
		else {
			f-=ed.f;
		}
	}

	return f;

}


int main(int argc,char** argv){
	flowgraph* g;
	unsigned i, k, f1, f2;
	unsigned s=atoi(argv[2]),t=atoi(argv[3]);
	time_t t0,t1,t2;

	t1=time(NULL);
	t0=t1;

	printf("Reading edgelist from file %s\n",argv[1]);
	g=readedgelist(argv[1]);
	printf("Number of nodes: %u\n",g->n);
	printf("Number of edges: %u\n",g->e);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Building datastructure\n");

	mkfg(g);
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Computing maximum flow between nodes s=%u and t=%u\n",s,t);

	f2=0;
	do{
		printf("flow = %u\n",f2);
		f1=f2;
		f2=onepath(g,s,t);
	} while (f2!=0);

	printf("maximum flow = %u\n",f1);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Writting the results in file %s\n",argv[4]);
	printres(g,argv[4]);
	freefg(g);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	return 0;
}
