
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


//******* BEGIN CYCLIC METHODS *******
bool cyclicTraverse(std::vector<edge> graph[], std::vector<edge> mst[], int curr, int colors[])
{
    //set color of current vertex to GRAY, meaning it is being processed right now
    colors[curr] = GRAY;
    
    //traverse through all edges
    for(auto next : graph[curr])
    {
        //doing this extra for loop check due to the way I have the MST stored (ie edges in the destination slots)
        //std::cout << "Graph edge: " << curr << " --" << next.weight << "--> " << next.destination << "\n";
        bool traverse = false;
        for(auto mstCheck : mst[next.destination])
        {
            //std::cout << "MST edge: " << mstCheck.start << " --" << mstCheck.weight << "--> " << mstCheck.destination << "\n";
            if(next.destination == mstCheck.destination && next.weight == mstCheck.weight)
            {
                traverse = true;
                break;
            }
        }
        
        if(traverse)
        {
            int vertex = next.destination;
            std::cout << "Traversing from " << curr << " to " << vertex << "\n";
            //check if GRAY, if so there is a cycle
            if(colors[vertex] == GRAY)
            {
                std::cout << "Cycle found on edge " << curr << " to " << vertex << "\n";
                return true;
            }

            //if vertex not processed, but a back edge exists, return true
            if(colors[vertex] == WHITE && cyclicTraverse(graph, mst, vertex, colors))
            {
                return true;
            }
        }
    }
    
    //mark vertex as processed
    colors[curr] = BLACK;
    
    return false;
}

bool isCyclic(std::vector<edge> graph[], std::vector<edge> mst[], int numVertices)
{
    //initialize all vertex color's to white (unprocessed)
    int *colors = new int[numVertices];
    for(int i = 0; i < numVertices; i++)
        colors[i] = WHITE;
    
    //DFS all vertices
    for(int i = 0; i < numVertices; i++)
        if(colors[i] == WHITE)
            if(cyclicTraverse(graph, mst, i, colors) == true)
                return true;
    
    return false;
}
//******* END CYCLIC METHODS *******

//*************** END GRAPH CHECK METHODS (CONNECTED & CYCLIC) ***************


//*************** BEGIN BORUVKAS ALGORITHM METHODS ***************
int find(int parents[], int vert)
{
    if(parents[vert] == vert)
        return vert;
    return find(parents, parents[vert]);
}

void joinComps(int parents[], int rank[], int startOne, int startTwo)
{
    int rootOne = find(parents, startOne);
    int rootTwo = find(parents, startTwo);

    //attach smaller ranked tree under root of higher ranked one
    if(rank[rootOne] < rank[rootTwo])
        parents[rootOne] = rootTwo;
    else if(rank[rootOne] > rank[rootTwo])
        parents[rootTwo] = rootOne;
    //if ranks are the same, just make root one the main root and increment its rank
    else
    {
        parents[rootTwo] = rootOne;
        rank[rootOne]++;
    }
}
    
void minSpanningTree(std::vector<edge> graph[], std::vector<edge> result[], int root, int numVertices, int numEdges, bool notConnected[])
{
    //number of trees initially is just the number of vertices
    int numTrees = numVertices-1;
    //remove any vertices not connected
    for(int i = 1; i < numVertices; i++)
    {
        if(notConnected[i])
            numTrees--;
    }
    
    //array of vector of edges representing each tree
    std::vector<edge> *trees = new std::vector<edge>[numVertices];
    //array of each parent vertex (index holds parent value)
    int *parents = new int[numVertices]();
    //rank array to hold bigger values for bigger roots
    int *rank = new int[numVertices]();
    //vector to store smallest edge of each tree, reinitialize on every loop
    std::vector<edge> *edges = new std::vector<edge>[numVertices];
    
    //vars to hold vertex numbers
    int rootOne, rootTwo;
    
    //initialize parent array
    for(int i = 1; i < numVertices; i++)
    {
        parents[i] = i;
        rank[i] = 0;
    }
    
    //std::cout << "Number of trees to begin: " << numTrees << "\n"; 
    
    //debug var
    int edgeCount = 0;
    
    //while there are more than 1 tree (assuming entire graph is connected to begin), keep going
    while(numTrees > 1)
    {
        //clear all edges
        for(int i = 1; i < numVertices; i++)
            edges[i].clear();
        
        //loop through all edges and update the cheapest one for each component
        for(int i = 1; i < numVertices; i++)
        {
            //if the current vertex is not connected to the main graph, jump to next iteration
            if(notConnected[i])
                continue;
            //std::cout << "Looking through edges for vertex: " << i << "\n";
            for(auto curr : graph[i])
            {
                //std::cout << "Current edge " << i << " --" << curr.weight << "-- " << curr.destination << "\n";
                rootOne = i; //rootOne = edge "start"
                rootTwo = curr.destination; //rootTwo = edge "end"
                
                //find components/sets for every edge
                int setOne = find(parents, rootOne);
                int setTwo = find(parents, rootTwo);
                
                //if the two vertices found are part of the same component, ignore
                if(setOne != setTwo)
                {
                    if(edges[setOne].empty() || (edges[setOne].front().weight > curr.weight))
                    {
                        /*
                        if(!edges[setOne].empty())
                        {
                            std::cout << "S1 Edge " << rootOne << " --" << curr.weight << "-- " << rootTwo << " is less than \n";
                            std::cout << "        " << edges[setOne].front().start << " --" << edges[setOne].front().weight << "-- " << edges[setOne].front().destination << "\n";
                        }
                        */
                        edges[setOne].clear();
                        edges[setOne].push_back({rootOne, rootTwo, curr.weight});
                        //maybe push back for edges[setTwo] as well?
                        /*
                        for(auto test : edges[i])
                        {
                            std::cout << "S1 Cheapest edge for vertex " << setOne << " is " << setOne << " --" << test.weight << "-- " << test.destination << "\n";
                        }
                        */
                    }
                    if(edges[setTwo].empty() || (edges[setTwo].front().weight > curr.weight))
                    {
                        /*
                        if(!edges[setTwo].empty())
                        {
                            std::cout << "S2 Edge " << rootOne << " --" << curr.weight << "-- " << rootTwo << " is less than \n";
                            std::cout << "        " << edges[setTwo].front().start << " --" << edges[setTwo].front().weight << "-- " << edges[setTwo].front().destination << "\n";
                        }
                        */
                        edges[setTwo].clear();
                        edges[setTwo].push_back({rootOne, rootTwo, curr.weight});
                        //maybe push back for edges[setOne] as well?
                        /*
                        for(auto test : edges[i])
                        {
                            std::cout << "S2 Cheapest edge for vertex " << setTwo << " is " << setTwo << " --" << test.weight << "-- " << test.destination << "\n";
                        }
                        */
                    }
                }
            }
        }
        
        /*
        //debug - print out cheapest edges
        for(int i = 1; i < numVertices; i++)
        {
            for(auto curr : edges[i])
            {
                std::cout << "Cheapest edge for tree rooted at " << i << " is " << curr.start << " --" << curr.weight << "-- " << curr.destination << "\n";
            }
        }
        */
        
        //add the cheapest edges to the MST
        for(int i = 1; i < numVertices; i++)
        {
            //check if there is a cheapest edge for the current set
            if(!edges[i].empty())
            {
                rootOne = edges[i].front().start; //rootOne = edge "start"
                rootTwo = edges[i].front().destination; //rootTwo = edge "end"
                
                //find components/sets for every edge
                int setOne = find(parents, rootOne);
                int setTwo = find(parents, rootTwo);
                
                //make sure both setOne and setTwo are not within same component
                if(setOne != setTwo)
                {
                    //add edge to mst
                    addEdge(result, rootOne, rootTwo, edges[i].front().weight);
                    edgeCount++; //debug
                    
                    //combine the trees
                    joinComps(parents, rank, rootOne, rootTwo);
                    numTrees--;
                }
            }
        }
        //std::cout << "EdgeCount = " << edgeCount << "\n";
        //std::cout << "numTrees: " << numTrees << "\n";
        //numTrees--;
    }
    
    std::cout << "Number of edges in the MST is: " << edgeCount << "\n";
    
    delete[] parents;
    delete[] rank;
    delete[] trees;
    delete[] edges;
    
    return;
}
//*************** END BORUVKAS ALGORITHM METHODS ***************



//*************** BEGIN MAIN FUNCTION ***************
int main(int argc, char** argv) {

	//time variables
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
    else
        delimeter = ' ';
    //grab file name from cmd line
    std::string fname = argv[1]; 
    
    //pointer to input file containing graph data
    std::fstream file;

    //open the data file
    file.open(fname, std::ios::in);    

    //need the number of nodes in the graph, which is the first line of the csv
    if(file.is_open())
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
    //graph containing the resulting MST
    std::vector<edge> *result = new std::vector<edge>[numVertices];
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
    
	printf("Starting timer and running MST algorithm...\n");
    //start measuring time
	start = omp_get_wtime();	
    
	//find mst of the graph
    minSpanningTree(graph, result, root, numVertices, edges, notConnected);
    
    //get end time after ops finish
	end = omp_get_wtime();
	//calculate total time
	total = end - start;
    //print execution time
	printf("Total execution time: %.8lf seconds\n", total);
    
    //std::cout << "\nPrinting minimum spanning tree:\n";
    //print mst
    //printGraph(result, numVertices, checkVertex);
	//print total weight
	int weightTotal = 0;
	for(int i = 1; i < numVertices; i++)
	{
		for(auto curr : result[i])
        {
			weightTotal = weightTotal + curr.weight;
		}
	}
	std::cout << "Total weight for MST is: " << weightTotal << "\n";
    
    
    delete[] checkVertex;
    delete[] graph;
    delete[] result;
    
    return 0;
}
//*************** END MAIN FUNCTION ***************