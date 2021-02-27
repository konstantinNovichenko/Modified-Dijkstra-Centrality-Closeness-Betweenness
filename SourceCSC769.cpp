//Petingi
//Delete Dijkstra/DFS
//

//Modified by Konstantin Novichenko
//for CSC 769 at College of Staten Island

#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>

using namespace std;

// FUNCTION PROVIDED BY PROFESSOR PETINGI
//node
struct edge
{ int v1;
  int v2;
};
struct node
{
    int vertex;
	int binary;
	float rel;
    node * next;
	
};

struct distanceScore {
	int index;
	int dSi;
	int diT;
};

int nods =0;
int edgs =0;
int dim = 0;
int iterat = 0;

bool * Visited;
node * glovalnode;
//void Dijkstra_A (int a, node headnodes_A[],int L[], int n);

vector<int> seq[100][1000];
//int seq[100][1000];
int numseq[100];// number of sequences-vectors

 bool cycles=false;// determine if a graph has cycles.
 //ifstream infile;
 //ofstream outfile;
 

class Graph {
private:
	int n;
    int e; // number of edges
	int A,B;	
	float rel [100][100];
	int bin [100][100];
    //int temp[100];
	int D;
	edge x;
	int counter;
	
	ifstream infile;
	
	bool * T;
	int * L;
	int * father;
	node * w;
public:  
	 node * headnodes;//linked list
	 node * headnodes_A;
	 int  weight[100][100];
	 Graph(int nods, int edgs, int dim)	
	// construtor
	 {   

		 /****************NOTE: DON'T FORGET TO CHANGE THE PATH TO THE FILE WITH THE GRID*********************************/
		 infile.open("./grid.txt");//read input file
		 // file format a  b w  where (a,b) is the edge and w is the 
		 // weight of the edge.
		 counter;
		 n=nods;
		 e = edgs;
		 D= dim;
		 counter=0;
		 headnodes= new node [n];
		 //T=new bool [n];
		 //father=new int [n];
		 //headnodes_A= new node [n];
		 // cout << "ok until here 2" << endl;
		 // headnodes is an array of nodes.
		 for (int i=0; i < n; i++)
		 {  
			 headnodes[i].vertex=i;
			 headnodes[i].next=0;
			 headnodes[i].binary=1;
			 headnodes[i].rel=0.5; 
		 }
		 // cout << "ok until here 3" << endl;
	 }

	Graph(Graph &s) //copy constructor
	
    { node *newnode, *nextn, *prev;
	  n=s.n;
      e=s.e;
      D=s.D;
	  A=s.A;
	  B=s.B;
	  x=s.x;
      headnodes= new node [n]; 
	  for (int i=0; i < n; i++)
	  {     headnodes[i].binary = s.headnodes[i].binary;
            headnodes[i].vertex = s.headnodes[i].vertex;
			headnodes[i].rel = s.headnodes[i].rel;
			headnodes[i].next = s.headnodes[i].next;
			nextn=&s.headnodes[i];
			prev=&headnodes[i];
			while (!nextn->next)
			  { newnode = new node;
                newnode->binary = nextn->next->binary;
                newnode->vertex = nextn->next->vertex;
				newnode->rel    = nextn->next->rel;
                newnode->next   = nextn->next->next;
				prev->next=newnode;
				nextn=nextn->next;
				prev=newnode;
			}
	   }
	         //cout << endl << "display"; 
			 //display();
             for (int j=0; j<100; j++)
				 for (int k=0; k<100; k++)
				 { weight[j][k]= s.weight [j][k];
	               rel [j][k]= s.rel[j][k];
				   bin [j][k] = s.bin [j][k];
			     }

	} 
	

	// FUNCTION PROVIDED BY PROFESSOR PETINGI
	void create()
	//create function
	{
		node *pre;
		node * nextn;
	    node *newnode;
		
		for (int i=0; i <e; i++)
		{  
			for (int j=0; j <e; j++)
			{
				weight[i][j] = 0;		//initializing weights to 0
				//cout<<weight[i][j]<<"\t";
			}
		}

		//cout<<"Enter the pairs of vertices E=(A,B) you wish to conect for:\n";
		//displaying edges and wheight
		//cout << "ok until here" << endl;
		

		for (int i=1; i <=e; i++)
		{  
			
			/*cout<<"Edge#"<<i<<":\nA: ";*/
			infile >>A >> B;
			
			/*cout<<"Enter wheight: ";*/
			//weight [A][B]=1;
            //cout<<weight[i][j]<<"\t";
			infile>>weight[A][B];
			weight[B][A]= weight[A][B];
			
			newnode= new node;
            newnode->vertex = B;
			newnode->binary=1;
			newnode->rel=0.5;


            if( headnodes[A].next == NULL )
			{
				newnode->next= NULL;
				
				headnodes[A].next=newnode;
			}
			else
			{
				pre= &headnodes[A];
				while( pre->next != NULL )
				{   
					 pre = pre->next;
				}
				newnode->next = NULL;
				pre->next = newnode;
			}


		//ADJACENT NODES
		
			newnode= new node;
			newnode->vertex = A;
			newnode->binary = 1;
			newnode->rel =0.5;

			if( headnodes[B].next == NULL )
			{
				newnode->next= NULL;
				headnodes[B].next=newnode;
			}
			else
			{
				pre= &headnodes[B];
				while( pre->next != NULL )
				{
					pre = pre->next;
				}
				newnode->next = NULL;
				pre->next = newnode;
			}
			
		}
     
     infile.close();

	}

	// FUNCTION PROVIDED BY PROFESSOR PETINGI
	void DFS(int father, int v)
	// DFS function
	{  
	Visited [v]=true;
	bool adjtoa=true;

    node* adjnode=headnodes[v].next;
    while (adjnode) // visit all vertices adjacent to v
	{
		if (!Visited[adjnode->vertex])
		{//if adjacent vertex to v was not visited previously
                      DFS(v,adjnode->vertex);
		}
		else if (father !=adjnode->vertex) // if the vertex adjacent to v is not the father, we have a 
		{ 
			cycles=true;
        }

		adjnode = adjnode->next;

	}
	}
	
	// FUNCTION PROVIDED BY PROFESSOR PETINGI
	void Dijkstra (int source)
		//Dijkstra function
       {
		   int infinity = 1000;
		   int s =source,Minimum,u;
		   int temp[11];			//temporary data structure which is needed to 
		   T = new bool[n];//Nodes to be visted
		   L= new int[n];//lambda
		   father = new int[n];
		   
		//INITIALIZING DATA STRUCTURES
		   for(int i =0;i<n;i++)
			{
			T[i] = false;			//T=V; all the vertices are eligible to be visited
			L[i]=infinity;			// at the beginning every vertex is at distance ?, from s
			temp[i] = infinity;		
			father[i] =-1;
			}
		   //WE ASSUME THE SOURCE IS 0 FOR EVERY GRAPH
		   L[s]=0;					// s is at distance 0 of itself
		   temp[s] =0;				// Temp has the same contents of L 
		  
		  // Let u be the vertex in T that has minimum L clearly at the beginning u=s 
		   for(int i = 0; i < n; i++)				//LOOP TROUGH ALL THE VERTICES OF THE GRAPH
		  {  
			  //cout<<endl<<"STEP "<<i<<":\n________ ";
			  Minimum = infinity;
		      for(int j = 0; j < n; j++)
			  {   if (T[j]==0){
				  if( Minimum > temp[j])
					   {
						   Minimum = temp[j];			//finding minimum L[s]
						   u = j;
					   }
			   }
			  }
			  temp[u] = infinity;				//Assigning INIFINITY to the data structure allready visited to find the next minimum L
			  //cout<<"\nU : "<<u;
			  w = &headnodes[u];				//Assigning address of headnodes[u] to find adjacent nodes for u
			  //cout<<"\tW = "<<w->vertex;
			  while( w!=NULL )					//while w adjacent to u 
			  {
				  if(T[w->vertex]==false)		// if w Exist in T, proceed
				  {
					  if (L[w->vertex]> L[u]+ weight[u][w->vertex])
					  { 
						  // if by reaching w from u has less weight
                                L[w->vertex]= L[u]+ weight[u][w->vertex]; // w is closer to s by using u;
								//cout<<"\nL[w]= L[u]+ weight(u,w)\n";
								//cout<<L[w->vertex]<<"  =   "<<L[u]<<"   +   "<<weight[u][w->vertex]<<endl;
								temp[w->vertex] = L[w->vertex];
                                father[w->vertex]=u;
								//cout<<"father[w] = "<<u<<endl;
						}
					}
			
				w = w->next;   // tranfer the address of 'temp->next' to 'temp'
			  }
			  T[headnodes[u].vertex] = true;				//Discard visited vertex u from T
		
			}


		   //DISPLAYING DATA STRUCTURES
		cout<<endl<<"--------------------------------------\nINDEX: \n";
		for(int i = 0;i<n;i++)
		{
			cout<<i<<"\t";
		}
		cout<<endl<<"--------------------------------------\nStructure T (visited): \n";
		for(int i = 0;i<n;i++)
		{
			cout<<T[i]<<"\t";
		}
		cout<<endl<<"Structure L[s](Closest Distances from S = 0): \n";
		for(int i = 0;i<n;i++)
		{
			cout<<L[i]<<"\t";
		}
	    cout<<endl<<"Structure father(routing): \n";
	    for(int i = 0;i<n;i++)
		{
			cout<<father[i]<<"\t";
		}
		cout<<"\n--------------------------------------"<<endl;
		cout<<"\t\tSHORTEST PATH: "<<n<<"->";
		for(int i = n-1;i>0;i--)
		{
			cout<<father[i]<<"->";
		}
		cout<<endl;
	}


	// Modified Dijkstra slgorithms that returns a list of all the shortes distances from the source to each vertex
	int* modifiedDijkstraCentrality(int source)		//modified Dijkstra function
	{
		int infinity = 1000;
		int s = source, Minimum, u;
		int temp[11];			//temporary data structure which is needed to 
		T = new bool[n];//Nodes to be visted
		L = new int[n];//lambda
		father = new int[n];

		//INITIALIZING DATA STRUCTURES
		for (int i = 0; i < n; i++)
		{
			T[i] = false;			//T=V; all the vertices are eligible to be visited
			L[i] = infinity;			// at the beginning every vertex is at distance ?, from s
			temp[i] = infinity;
			father[i] = -1;
		}
		//WE ASSUME THE SOURCE IS 0 FOR EVERY GRAPH
		L[s] = 0;					// s is at distance 0 of itself
		temp[s] = 0;				// Temp has the same contents of L 

	   // Let u be the vertex in T that has minimum L clearly at the beginning u=s 
		for (int i = 0; i < n; i++)				//LOOP TROUGH ALL THE VERTICES OF THE GRAPH
		{
			//cout<<endl<<"STEP "<<i<<":\n________ ";
			Minimum = infinity;
			for (int j = 0; j < n; j++)
			{
				if (T[j] == 0) {
					if (Minimum > temp[j])
					{
						Minimum = temp[j];			//finding minimum L[s]
						u = j;
					}
				}
			}
			temp[u] = infinity;				//Assigning INIFINITY to the data structure allready visited to find the next minimum L
			//cout<<"\nU : "<<u;
			w = &headnodes[u];				//Assigning address of headnodes[u] to find adjacent nodes for u
			//cout<<"\tW = "<<w->vertex;
			while (w != NULL)					//while w adjacent to u 
			{
				if (T[w->vertex] == false)		// if w Exist in T, proceed
				{
					if (L[w->vertex] > L[u] + weight[u][w->vertex])
					{
						// if by reaching w from u has less weight
						L[w->vertex] = L[u] + weight[u][w->vertex]; // w is closer to s by using u;
						//cout<<"\nL[w]= L[u]+ weight(u,w)\n";
						//cout<<L[w->vertex]<<"  =   "<<L[u]<<"   +   "<<weight[u][w->vertex]<<endl;
						temp[w->vertex] = L[w->vertex];
						father[w->vertex] = u;
						//cout<<"father[w] = "<<u<<endl;
					}
				}

				w = w->next;   // tranfer the address of 'temp->next' to 'temp'
			}
			T[headnodes[u].vertex] = true;				//Discard visited vertex u from T

		}

		// Return the list
		return L;
	}

};// end class

// A utility function used to print the solution for closeness
void printClosenessArr(std::vector<int> dist, int n)
{
	//printf("Vertex Distance from Source\n");
	for (int i = 0; i < n; ++i)
	{
		cout << i << "\t\t" << dist[i] << "\t\t" << ((float)(n - 1) /(float)dist[i]) << endl;
	}
		
}

// A utility function used to print the solution for betweeness
void printBetweennessArr(std::vector<float> dist, int n)
{
	//printf("Vertex Distance from Source\n");
	for (int i = 0; i < n; ++i)
	{
		cout << i << "\t\t" << dist[i] << endl;
	}
	
}


// Function that computes and displays Centrality Closeness points
void findCentralityCloseness(Graph g)
{
	// data for tracking the results
	std::vector<int> centralityScores;	
	int tempCentralityScore = 0;
	int bestCentralityScore = -1;
	std::vector<int> bestCentralityPointsIndices;

	for (int i = 0; i < nods; i++) // go through every node
	{		
		// Run modified Dijkstra
		int* d = g.modifiedDijkstraCentrality(i); // array of distance values from the point
		for (int j = 0; j < nods; j++)
		{
			tempCentralityScore += d[j]; // get the total distance from the vertix to other verticies
		}

		centralityScores.push_back(tempCentralityScore); // push to the list of all centrality scores

		if (bestCentralityScore > tempCentralityScore || bestCentralityScore == -1) // the first or a better result
		{		
			bestCentralityPointsIndices.clear();
			bestCentralityScore = tempCentralityScore;			
			bestCentralityPointsIndices.push_back(i);
		}
		else if (bestCentralityScore == tempCentralityScore) // as good as the other best vertices, add to the list
		{
			bestCentralityPointsIndices.push_back(i);
		}

		tempCentralityScore = 0; // reset
		
		// cleat the memory
		delete[] d;
		d = nullptr;
	}

	std::cout << "\n\n__________________________________________________________________________________\n";
	std::cout << "Closeness Centrality Scores\n\n";
	std::cout << "Vertex \t\tScore \t\tInverted Score\n";
	printClosenessArr(centralityScores, nods);

	std::cout << "\n\nThe best Centrality Closeness score is " << centralityScores[bestCentralityPointsIndices[0]] << " (inverted score: " << (float)(nods - 1) / (float)centralityScores[bestCentralityPointsIndices[0]] << ") by Verticies:\n";

	for (int i = 0; i < bestCentralityPointsIndices.size(); i++)
	{
		 cout << "#" << bestCentralityPointsIndices[i] << "\n";
	}	

}

// Function that recursively traversing the graph looking for all shortest paths.
// Sets itself as a Parent vertix
bool evaluatePath(Graph g, int indexPrev, int indexCS, int indexT, int* dT, vector<int> cP, vector<vector<int>> &p)
{	
	node tempNode = g.headnodes[indexCS]; // copy of the current node	

	while (tempNode.next != NULL) // visit all the adjacent verticies
	{
		vector<int> tempPath = cP; // current shortest path		
		tempPath.push_back(indexCS); // add this vertix to the path

		if (tempNode.vertex != indexPrev) // Check if Parent (prevent from going back in recursive calls)
		{	
			// Evaluate the distance d(t,i) + w(i, s) = d(t,s) if true, continue traversing
			if (dT[tempNode.next->vertex] + g.weight[tempNode.next->vertex][indexCS] == dT[indexCS])
			{
				if (tempNode.next->vertex == indexT) // arrived to the destination
				{
					tempPath.push_back(tempNode.next->vertex); // add the last vertix of the current shortest path to the list
					p.push_back(tempPath);	// push this shortest path to the lsit of ALL shortest paths								
				}
				else  // continue traversing
				{
					evaluatePath(g, indexCS, tempNode.next->vertex, indexT, dT, tempPath, p);
				}
			}			
		}
		
	    tempNode = *tempNode.next;	// set to the next adjacent vertix					
	}

	if (tempNode.vertex != indexT) //finished but haven't arrived to the destination
	{
		return false;
	}				
	
	return true;
}


// Function that starts traversing the graph looking for all shortest paths.
// The first vertex that sets itself as a Parent vertex
void traversingAllShortestPaths(Graph g, int indexS, int indexT, int* dT, vector<vector<int>> &p)
{	
	node tempHeadNode = g.headnodes[indexS]; // copy of the current node	
	
	while (tempHeadNode.next != NULL) // visit all the adjacent verticies
	{
		vector<int> tempP;
		tempP.push_back(indexS); // start the current shortest path	
		
	
		if (tempHeadNode.next->vertex == indexT) // check if the adjacent vertix is the destination
		{
			tempP.push_back(tempHeadNode.next->vertex);	// add the last vertix of the shortest path to the list		
			p.push_back(tempP);	// push this shortest path to the lsit of ALL shortest paths
			
		}
		// Evaluate the distance d(t,i) + w(i, s) = d(t,s) if true, continue traversing
		else if (dT[tempHeadNode.next->vertex] + g.weight[tempHeadNode.next->vertex][indexS] == dT[indexS])
		{			
			evaluatePath(g, tempHeadNode.vertex, tempHeadNode.next->vertex, indexT, dT, tempP, p);			
		}		

		tempHeadNode = *tempHeadNode.next; // set to the next adjacent vertix		
	}	

}


// Function that computes and displays the betweenness scores for each node
void findBetweenness(Graph g)
{
	std::vector<float> sigmaValues; // array total betweeness scores Sigma

	// Check all the possible paths on the graph
	for (int j = 0; j < nods - 1; j++) 
	{
		for (int k = j + 1; k < nods; k++)
		{
			
			vector<int> currentSigmaValues; // array of the current betweeness scores Sigma

			int* dS = g.modifiedDijkstraCentrality(j); // array of distance values from the point S
			int* dT = g.modifiedDijkstraCentrality(k); // array of distance values from the point T			

			// the list of all shortest paths for the current set of S and T
			vector<vector<int>> listAllSPaths;

			
			vector<vector<int>> allSPathS; // the list of shortest paths from i to S
			vector<vector<int>> allSPathT; // the list of shortest paths from i to T

			// Check all verticies
			for (int i = 0; i < nods; i++)
			{
				// Populate current sigma values with zeroes
				currentSigmaValues.push_back(0);
				
				// if the current vertix is on the shortest path, traverse and find all the shortest paths it's positioned in
				if (dT[i] + dS[i] == dS[k])
				{					
					traversingAllShortestPaths(g, i, j, dS, allSPathS); // i -> s paths					
					traversingAllShortestPaths(g, i, k, dT, allSPathT);	// i -> t paths				
				}
				
				// if the first shortest path wasn't added yet (all new paths will be unique)
				if (listAllSPaths.size() == 0)
				{
					if (allSPathS.size() != 0 && allSPathT.size() != 0) // if both parts have at least one piece of the shortest path each
					{
						// Combine all pairs of pieces of the path together and push to the list of all shortest paths
						for (int z = 0; z < allSPathS.size(); z++)
						{
							vector<int> p;
							for (int x = allSPathS[z].size() - 1; x >= 0; x--)
							{
								p.push_back(allSPathS[z][x]);
							}

							for (int kh = 0; kh < allSPathT.size(); kh++)
							{
								for (int h = 1; h < allSPathT[kh].size(); h++)
								{
									p.push_back(allSPathT[kh][h]);
								}
							}
							
							listAllSPaths.push_back(p);
						}						
					}
				}
				else // Some shortest paths exist
				{
					if (allSPathS.size() != 0 && allSPathT.size() != 0) // if both parts have at least one piece of the shortest path each
					{
						// Combine all pairs of pieces of the path together and push to the list of all shortest paths if unique
						for (int z = 0; z < allSPathS.size(); z++)
						{
							vector<int> p;
							for (int x = allSPathS[z].size() - 1; x >= 0; x--)
							{
								p.push_back(allSPathS[z][x]);
							}

							for (int kh = 0; kh < allSPathT.size(); kh++)
							{
								for (int h = 1; h < allSPathT[kh].size(); h++)
								{
									p.push_back(allSPathT[kh][h]);
								}
							}

							int unique = 0; // counter of passing the uniqueness test

							for (int ia = 0; ia < listAllSPaths.size(); ia++)
							{
								for (int ib = 0; ib < listAllSPaths[ia].size(); ib++)
								{
									if (listAllSPaths[ia][ib] != p[ib]) // unique path, update the counter
									{
										unique++;
										ib = listAllSPaths[ia].size();
									}
								}
							}

							// uniqueness test passed, add the path to the list of all shortest paths
							if (unique == listAllSPaths.size())
							{
								listAllSPaths.push_back(p);
							}
						}
					}
				}

				// reset lists with the pieces of the shortest path
				allSPathS.clear();
				allSPathT.clear();
			}			

			// Check if the S is adjacent to T and add this path to the list of sortest paths
			node tempNSource = g.headnodes[j];
			while (tempNSource.next != NULL)
			{
				if (tempNSource.next->vertex == k)
				{
					vector<int> directP;
					directP.push_back(j);
					directP.push_back(k);
					listAllSPaths.push_back(directP);
				}
				tempNSource = *tempNSource.next;
			}

			// if the list of all shortest paths is not empty
			if (listAllSPaths.size() != 0)
			{
				// get all the current value of occurence for each vertix from the list of all shortest paths
				for (int ixh = 0; ixh < listAllSPaths.size(); ixh++)
				{
					for (int yh = 1; yh < listAllSPaths[ixh].size() - 1; yh++)
					{
						currentSigmaValues[listAllSPaths[ixh][yh]]++;
					}
				}

				// if the total sigma list is empty, push into the vector array
				if (sigmaValues.size() == 0)
				{
					for (int z = 0; z < currentSigmaValues.size(); z++)
					{
						// C(B through v) = SUM(all_shortest_paths_through_v/all_the_shortest_paths)
						sigmaValues.push_back((float)currentSigmaValues[z] / (float)listAllSPaths.size());
					}
				}
				else // pass by index
				{
					for (int z = 0; z < currentSigmaValues.size(); z++)
					{
						// C(B through v) = SUM(all_shortest_paths_through_v/all_the_shortest_paths)
						sigmaValues[z] += ((float)currentSigmaValues[z] / (float)listAllSPaths.size());
					}
				}
			}
			
			// clear the dynamic variable from the memory
			delete[] dS;
			delete[] dT;
			dS = nullptr;
			dT = nullptr;
		}		
		
	}

	// Evaluate the best values
	float bestBetwSc = -1;
	vector<int> bestBetIndex;

	// Find the best scores
	for (int i = 0; i < sigmaValues.size(); i++)
	{
		if (i == 0)
		{
			bestBetwSc = sigmaValues[i];
			bestBetIndex.push_back(i);
		}
		else
		{
			if (bestBetwSc < sigmaValues[i])
			{
				bestBetwSc = sigmaValues[i];
				bestBetIndex.clear();
				bestBetIndex.push_back(i);
			}
			else if (bestBetwSc == sigmaValues[i])
			{
				bestBetIndex.push_back(i);
			}
		}
	}

	// OUTPUT
	std::cout << "\n\n__________________________________________________________________________________\n";
	std::cout << "Betweenness Centrality Scores\n\n";
	std::cout << "Vertex \t\tScore\n";
	
	printBetweennessArr(sigmaValues, nods);

	std::cout << "\n\nThe best Centrality Betweenness score is " << bestBetwSc << " by Verticies:\n";

	for (int i = 0; i < bestBetIndex.size(); i++)
	{
		cout << "#" << bestBetIndex[i] << "\n";
	}	
	cout << "\n\n";
}




int main()
{  
    //ofstream outfile;
	//outfile.open("E:/School/2020 (3) FALL/CSC 769/Project1/complete.txt"); // NOTE: specify your local path to the file

	int dim=0;
	cout<<"Please enter the number of vertices: ";
	cin >>nods;
    cout<<"Enter the number Edges: ";
    cin >>edgs;
	
    Visited=new bool [nods];
	for (int j=0; j< nods; j++){
		numseq[j]=0;
	}
    int Nmbrcmpnts=0; //we initialize the counter for the number of components

    Graph G(nods,edgs,dim);
    G.create(); 
	
	// Evaluate Closeness Centrality
	findCentralityCloseness(G);
	// Evaluate Betweenness Centrality
	findBetweenness(G);	
	

    system("pause");
	return 0; 
}



