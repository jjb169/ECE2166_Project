#include <stdio.h>
#include <chrono>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string> // std::string
#include <iostream> // std::cout
#include <algorithm>
#include <sstream> // std::stringstream
#include <list>
#include<queue>
#include <omp.h>

enum Color {WHITE, GRAY, BLACK};

//*************** BEGIN GRAPH REPRESENTATION FUNCTIONS ***************
struct edge
{
    int start;
    int destination;
    int weight;
};


void addEdge(std::vector<edge> graph[], int start, int dest, int weight)
{
    graph[start].push_back({start, dest, weight}); //add a struct directly with the use of {...}
    graph[dest].push_back({dest, start, weight}); //add a struct directly with the use of {...}
    //std::cout << "Just added for vertex " << start << " edge with dest/weight = " << dest << "/" << weight << "\n";
}

//print the adjacency list representation of graph
int printGraph(std::vector<edge> adj[], int length, int vertices[])
{
    //count var for # of edges
    int count = 0;
    for (int vertex = 0; vertex < length; vertex++) 
    {
        //make sure vertex exists/is in use
        if(vertices[vertex] == 1)
        {
            if(vertex != 0)
            {
                std::cout << "\nAdjacency list of vertex " << vertex << "\n";
                for (auto curr : adj[vertex])
                {
                    std::cout << vertex << " -- " << curr.weight << " -- " << curr.destination << "\n";
                    count++;
                }
                std::cout << "\n";
            }
        }
    }
    
    return count;
}

//actually construct the weighted directed graph
void constructGraph(std::vector<edge> graph[], char delimeter, std::fstream &file, int numVertices, int vertices[], int info[])
{
    //edge info 
    int source = 0, destination = 0, weight = 0, itemNum = 0, currRoot = -1, numEdges = 0;
    //string to hold a line of the graph (SOURCE, TARGET, WEIGHT)
    std::string row; 
    std::string curr;
    //array to hold the amount of edges for each vertex - used to determine root
    int *edgeCount = new int[numVertices]();
    
    //debugging var to find highest vertex #
    int currMax = 0;
    //loop through csv and grab all info
    if(file.is_open())
    {
        //while lines (edges) are still available, continue
        while(std::getline(file, row))
        {
            //std::cout << "Row: " << row << "\n";
            //break the row into its components
            std::stringstream items(row);
            //reset the itemNum var
            itemNum = 0;

            //take out the parts and add entry to the graph
            while(std::getline(items, curr, delimeter))
            {
                //std::cout << "Curr: " << curr << " and itemNum = " << itemNum << "\n";                
                //if itemNum == 0, it is the source vertex
                if(itemNum == 0)
                {
                    source = std::stoi(curr);
                    //debugging
                    if(source > currMax)
                        currMax = source;
                    //mark the current vertex as existing
                    vertices[source] = 1;
                    
                }
                else if(itemNum == 1) 
                {
                    //if itemNum == 1, it is the destination vertex
                    destination = std::stoi(curr);
                    //debugging
                    if(destination > currMax)
                        currMax = destination;
                    //mark the current vertex as existing
                    vertices[destination] = 1;
                }
                else if(itemNum == 2) 
                    weight = std::stoi(curr); //if itemNum == 2, it is the weight of the edge
                //don't care for anything past weight, ie time epoch
                //increment itemNum on each loop
                itemNum++;
            }
            //create a new entry for this edge
            addEdge(graph, source, destination, weight);
            numEdges++;
            
            //increment counter for the amount of edges in this node
            edgeCount[source]++;
            //if it now has more edges than the previous leader, make it the current root
            if(currRoot == -1 || edgeCount[source] > edgeCount[currRoot])
                currRoot = source;
        }
    }
    else
    {
        std::cout << "Could not open csv file, please ensure it is in the correct location with correct name and retry\n\n";
        exit(0);
    }
    //currRoot = 10; //TEMPORARY CHANGE THIS LATER ************************************************************************
    info[0] = currRoot;
    info[1] = numEdges;
    
    //debugging to find highest # of a vertex
    std::cout << "The max value of a vertex is: " << currMax << "\n";
    //free the edge count array
    free(edgeCount);

    return;
} 
//*************** END GRAPH REPRESENTATION FUNCTIONS ***************


//*************** BEGIN GRAPH CHECK METHODS (CONNECTED & CYCLIC) ***************

//******* BEGIN CONNECTED METHODS *******
void dfs(std::vector<edge> graph[], bool vis[], int curr)
{
    //mark this vertex as visited
    vis[curr] = true;
    //loop for dfs through each branch off of this vertex
    for(auto next : graph[curr])
    {
        if(!vis[next.destination])
        {
            //std::cout << "Visting from " << curr << " to " << next.destination << "\n";
            dfs(graph, vis, next.destination);
        }
    }
}

bool isConnected(std::vector<edge> graph[], int exist[], int start, int numVertices, bool notConnected[])
{
    //create two boolean arrays to mark nodes when visited
    bool *visitFirst = new bool[numVertices]();
    
    bool connected = true;
    
    //std::cout << "Starting at vertex " << start << " on original graph\n";
    
    //perform dfs on the edges starting at the root
    dfs(graph, visitFirst, start);
    
    //check if all nodes have been visited
    for (int i = 1; i <= numVertices; i++) 
    {
        //if any vertex has not been visited, graph is not connected
        if (!visitFirst[i] && exist[i])
        {
            //std::cout << "Vertex " << i << " is not connected\n";
            connected = false;
            notConnected[i] = true;
            //free(visitFirst);
            //free(visitSecond);
           // return false;
        }
    }
    
    delete[] visitFirst;
    //free(visitSecond);
    return connected;
}

//******* END CONNECTED METHODS *******

//*************** END GRAPH CHECK METHODS (CONNECTED & CYCLIC) ***************
void bfs(std::vector<edge> graph[], bool visited[], int numVertices, int root)
{   
	//create a queue to hold vertices in order seen for traversal
    std::list<int> q;
	//begin by visiting root and adding it to queue
    visited[root] = true;
    q.push_back(root);
	//while the queue is not empty, there are still vertices that need traversing
    while (!q.empty())
    {
        // dequeue front node and print it
        int v = q.front();
        q.pop_front();

		//for a given vertex, loop through its edges
        for (auto edge: graph[v])
        {
			//grab the destination for an edge to check if it has been visited
            int adjvertex = edge.destination;
            if (!visited[adjvertex])
            {
				//if not visited yet, mark it as now visited and add it to the queue to traverse
                visited[adjvertex] = true;
                q.push_back(adjvertex);
            }
        }
    }  
}
    
//*************** BEGIN MAIN FUNCTION ***************
int main(int argc, char** argv) {
    double start = 0.0;
	double end = 0.0;
	double total = 0.0;
    //variables coming from graph
    int numVertices = 0, root = -1, edges = -1;
    int info[2]; //[0] = root, [1] = # edges
    //string to hold a line of the graph (SOURCE, TARGET, WEIGHT)
    std::string row; 
    std::string curr;
    //var to hold what the delimeter between items is- 
    //ie if a graph is SOURCE/DEST/WEIGHT whether it will be shown as SOURCE DEST WEIGHT or SOURCE,DEST,WEIGHT
    char delimeter = ' ';
    //var to decide if the delimeter between values is a comma or space
    int delim = std::atoi(argv[2]);
    if(delim == 0)
        delimeter = ',';
    else if(delim == 1)
        delimeter = ' ';
	else
		delimeter = '	';
    //grab file name from cmd line
    std::string fname = argv[1];
	std::cout << "Graph file in use: " << fname << "\n";
    
    //pointer to input file containing graph data
    std::fstream file;

    //open the data file
    file.open(fname, std::ios::in);    

    //need the number of nodes in the graph, which is the first line of the csv
	//work around for the two graphs named like "201512020130.v128568730_e270234840.tsv"
    if(fname.compare("../graphs/201512020130.v128568730_e270234840.tsv") == 0)
	{
		numVertices = 128568730 + 1;
	}
	else if(fname.compare("../graphs/201512020330.v226196185_e480047894.tsv") == 0)
	{
		numVertices = 226196185 + 1;
	}
    else if(file.is_open())
    {
        //need to read until no % lead the line
        std::getline(file, row);
        //std::cout << "Row: " << row << "\n";
        int i = 0;
        while(row.at(i) == '%')
        {
            std::getline(file, row);
            //std::cout << "Row: " << row << "\n";
        }
        std::stringstream items(row);
        std::getline(items, row, delimeter);
        //grab the first line and assign value to numVertices
        //std::getline(file, row);
        //std::cout << "Row: " << row << "\n";
        numVertices = std::stoi(row) + 1; //need +1 for some graphs that start with 1 instead of 0, just defaulting to having this addition
    }
    else
    {
        std::cout << "Attempted file name: " << fname << "\n";
        std::cout << "Could not open csv file, please ensure it is in the correct location with correct name and retry\n\n";
        exit(0);
    }
    
    //helper array to see if a vertex actually exists - ie may have graph with 500 as max vertex number but only 300 total vertices in graph
    int *checkVertex = new int[numVertices]();
    //main pointer to adjacency list (will need partitioning later)
    std::vector<edge> *graph = new std::vector<edge>[numVertices];
    //debug statement
    std::cout << "# vertices in csv = " << numVertices << "\n";
    
    //construct the graph and return the root as the node with the largest amount of outgoing edges
    constructGraph(graph, delimeter, file, numVertices, checkVertex, info);
    root = info[0];
    edges = info[1];
    std::cout << "The root node for this graph is: " << root << "\n";
    std::cout << "The number of edges in this graph is: " << edges << "\n";
    
    
    //check if graph is connected
    bool *notConnected = new bool[numVertices](); //number of vertices not connected to the graph - not reached by traversal from root
    bool connected = isConnected(graph, checkVertex, root, numVertices, notConnected);
    if(!connected)
    {
        std::cout << "The initial graph is not connected!\n";
        int numOut = 0;
        for(int i = 1; i < numVertices; i++)
        {
            if(notConnected[i])
                numOut++;
        }
        std::cout << "Number of vertices not connected is: " << numOut << "\n\n";
    }
    else
        std::cout << "The initial graph is fully connected!\n\n";
    
	
    //std::cout << "Printing initial graph:\n";
    //print the initial graph
    //printGraph(graph, numVertices, checkVertex);
    
    // Start measuring time
    //auto begin = std::chrono::high_resolution_clock::now();
	bool *visited = new bool[numVertices]();
	
	//additional time vars, should prob move up top with others
	double allTime = 0.0, avgTime = 0.0;
	//going to run BFS 10 times and average the timing of the results
	for(int i = 0; i < 10; i++)
	{
		//reset this on every run
		for(int i = 0; i < numVertices; i++)
			visited[i] = false;
		printf("Begin to run and collect time for baseline....\n");
		// Begin timer and start BFS function
		start = omp_get_wtime();
		bfs(graph, visited, numVertices, root);
		//  Stop measuring time and calculate the elapsed time
		end = omp_get_wtime();
		//calculate total time
		total = end - start;
		//print execution time
		printf("Total execution time: %.8lf seconds\n", total);
		//add this runtime to the total
		allTime = total + allTime;
		
		//check for each run
		printf("Checking if all vertices have been visited...\n");
		bool allVisited = true;
		double numVisited = 0;
		//loop through all vertices and see if they are visited
		for(int i = 1; i < numVertices; i++)
		{
			//if a vertex is not connected within original graph, ignore
			//if a vertex is connected and has not been visited, change bool var to false
			if(!notConnected[i] && !visited[i])
			{
				allVisited = false;
				std::cout << "Vertex " << i << " not visited!\n";
				//break;
			}
			else if(!notConnected[i] && visited[i])
				numVisited++; //count each visited vertex
		}
		//print if all were visited or how many vertices were visited if the original graph was not fully connected
		if(allVisited)
		{
			if(connected)
			{
				printf("All %d vertices have been visited!\n\n", numVertices);
			}
			else
			{
				printf("Graph is not fully connected, %.0lf vertices reachable from the root have been visited!\n\n", numVisited);
			}
		}
		else
		{
			printf("Incomplete BFS traversal!\n\n");
			printf("Number of vertices visited: %.0lf\n\n", numVisited);
		}
	
	}
	
	//calculate the average time
	avgTime = allTime / 10;
	//print average runtime
	std::cout << "Average time across MST runs is: " << avgTime << "\n";
	
    
    delete[] checkVertex;
    delete[] graph;
    
    return 0;
}
//*************** END MAIN FUNCTION ***************