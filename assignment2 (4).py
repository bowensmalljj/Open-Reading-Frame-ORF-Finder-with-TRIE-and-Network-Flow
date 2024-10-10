# Q1
class Node():
    def __init__(self):
        """
        This init function used to initialise an array of length 4
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input: -
        Return:None
        Time complexity:O(N),N is 4
        Best case analysis:O(N)
        Worst case analysis:O(N)
        Space complexity:O(N)
        Input space analysis: O(1)
        Aux space analysis:O(N)
        """
        self.childnodes = [None] * 4
        self.storedindex = []

class SuffixTrie:
    def __init__(self, string):
        """
        This init function used to initialise node for suffixtrie and store the order number of the path that contain the character
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input: String
        Return:None
        Time complexity:O(N^2),N is length of string
        Best case analysis:O(N^2)
        Worst case analysis:O(N^2)
        Space complexity:O(N)
        Input space analysis: O(N)
        Aux space analysis:O(1)
        """
        self.root = Node()
        i=0
        while i<len(string):
            current = self.root
            for j in range(i, len(string)):
                character = string[j]
                index = ord(character) - 65
                if current.childnodes[index] is None:
                    current.childnodes[index] = Node()
                    current = current.childnodes[index]
                else:
                    current = current.childnodes[index]
                current.storedindex.append(i) #store the order/sequence index number of the path that contain the character
            i+=1




class OrfFinder:
    def __init__(self, genome):
        """
        This init function used to initialise suffixtrie
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input: Genome
        Return:None
        Time complexity:O(SuffixTrie)+O(N)
        Best case analysis:O(SuffixTrie)
        Worst case analysis:O(SuffixTrie)
        Space complexity:O(N)
        Input space analysis: O(N)
        Aux space analysis:O(1)
        """
        self.string = genome
        self.suffixtrie = SuffixTrie(self.string)

    def find(self, prefix, suffix):
        """
        This find function used to find the substring that start with prefix and end with suffix
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input: prefix, suffix
        Return: substring (list of substring that start with prefix and end with suffix)
        Time complexity:O(P)*O(S)
        Best case analysis:O(1)
        Worst case analysis:O(P)*O(S)
        Space complexity:O(P)+O(S)+O(N)
        Input space analysis: prefix=O(P), suffix=O(S);P is the length of prefix and S is the length of suffix
        Aux space analysis:O(N);N is the character in substring
        """
        substrings = []
        current = self.suffixtrie.root
        #get the sequence/order number of the index number of the path that contain the prefix
        for character in prefix:
            index = ord(character) - 65
            if current.childnodes[index] is None:
                starting=[]
                break
            current = current.childnodes[index]
            starting=current.storedindex
        for index in starting:
            start=index + len(prefix)
            end=len(self.string) - len(suffix) + 1
            #from starting index to end sequence
            for idx in range(start,end):
                # check if end sequence == suffix and append substring using slicing
                if self.string[idx:idx + len(suffix)] == suffix:
                    substrings.append(self.string[index:idx + len(suffix)])

        return substrings


# ==========
# Q2

class FlowNetwork:
    """
    Reference: https://www.w3schools.com/dsa/dsa_algo_graphs_fordfulkerson.php#:~:text=The%20Ford%2DFulkerson%20algorithm%20works,the%20maximum%20flow%20is%20reached.
    https://www.geeksforgeeks.org/breadth-first-search-or-bfs-for-a-graph/
    """
    def __init__(self, vertices):
        """
        This init function used to initialise flow network and residual network (matrix or nested array)
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input: vertices (number of nodes/vertices, first node is source and last node is target)
        Return:None
        Time complexity:O(N),N is number of vertices
        Best case analysis:O(N)
        Worst case analysis:O(N)
        Space complexity:O(N^2)
        Input space analysis: O(1)
        Aux space analysis:O(N^2)
        """
        self.networkflow=[]
        self.residual_network=[]
        self.length = vertices
        for i in range(vertices):
            self.networkflow.append([None]*vertices)
            self.residual_network.append([None]*vertices)

    def insert_edge(self,start,end,flow,capacity):
        """
        This insert_edge function used to create edges from specific nodes to specific nodes with flow and capacity
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input: start (from which node), end (to which node), flow,capacity
        Return:None
        Time complexity:O(1)
        Best case analysis:O(1)
        Worst case analysis:O(1)
        Space complexity:O(1)
        Input space analysis: O(1)
        Aux space analysis:O(1)
        """
        self.networkflow[start][end]=Edge(flow,capacity)
        if self.residual_network[start][end] is None:
            #create edge if edge doesn't exist
            self.residual_network[start][end] = Edge(flow, capacity).forward()
        else:
            #update edges if edge exists
            self.residual_network[start][end].reverse_update(flow)
        if self.residual_network[end][start] is None:
            # create edge if edge doesn't exist
            self.residual_network[end][start] = Edge(flow, capacity).reverse()
        else:
            # update edges if edge exists
            self.residual_network[end][start].reverse_update(flow)


    def ford_fulkerson(self):
        """
        This ford fulkerson function used to calculate the max flow
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input:
        Return:max flow
        Time complexity:O(T)*O(BFS)*(O(N)+O(P)+O(P))+O(N);T is the number of path,N is the number of vertices,P is the list of the current shortest path
        Best case analysis:O(1)
        Worst case analysis:O(T)*O(BFS)*(O(N)+O(P)+O(P))+O(N)
        Space complexity:O(P)+O(M), M is the list of minimum flow for each edge
        Input space analysis: O(1)
        Aux space analysis:O(P)+O(M)+O(BFS)
        """
        max_flow=0
        index=0
        while index is not None:
            output = self.BFS(self.residual_network)
            #if there is list of parentnodes
            if output is not None:
                index = output[-1]
                path=[]
                min_flow=[]
                #reverse the path from destination and back to source
                for i in range(self.length):
                    if index!=0:
                        path.append(index)
                        index=output[index]
                    else:
                        path.append(0)
                        break
                index=self.length-1
                #Get minimum flow
                for i in path:
                    if self.residual_network[i][index]!=None and self.residual_network[i][index].flow!=0:
                        min_flow.append(self.residual_network[i][index].flow)
                        index=i
                index = self.length - 1
                #update max flow and update flow for each edge
                for i in path:
                    if self.residual_network[i][index] != None:
                        self.residual_network[i][index].forward_update(min(min_flow))
                        self.residual_network[index][i].reverse_update(min(min_flow))
                        #update the flow network
                        if self.networkflow[i][index] is not None:
                            self.networkflow[i][index].forward_update(-min(min_flow))
                        else:
                            self.networkflow[index][i].forward_update(min(min_flow))
                        index = i

            else:
                #count the flow from source to every neighbour and get the max flow
                for i in self.networkflow[0]:
                    if i is not None:
                        max_flow += i.flow
                return max_flow


    def BFS(self,residual_network):
        """
        This BFS function used to find the shortest path from source to target
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input:residual network
        Return:parent_vertices_nodes (list that store the parent node of each node)
        Time complexity:O(T)*O(BFS)*(O(N)+O(P)+O(P));T is the number of path,N is the number of vertices,P is the list of the current shortest path
        Best case analysis:O(1)
        Worst case analysis:O(T)*O(BFS)*(O(N)+O(P)+O(P))
        Space complexity:O(P)+O(M), M is the list of minimum flow for each edge
        Input space analysis: O(1)
        Aux space analysis:O(P)+O(M)+O(BFS)
        """
        queue = []
        visited_vertices = [False] * len(residual_network[0])
        parent_vertices = [None] * len(residual_network[0])
        #append source to queue and set it to Visited =True
        queue.append(0)
        visited_vertices[0] = True

        while len(queue) != 0:
            current_vertex = queue.pop(0)
            i = 0
            #find neighbour nodes
            while i < len(residual_network[current_vertex]):
                # If not visited before and have edges and there is still flow
                if visited_vertices[i] == False and residual_network[current_vertex][i] is not None and residual_network[current_vertex][i].flow > 0:
                    queue.append(i)
                    visited_vertices[i] = True
                    parent_vertices[i] = current_vertex
                    # check if we ge to the sink node
                    if i == len(residual_network) - 1:
                        return parent_vertices
                i += 1

class Edge:

    def __init__(self,flow,capacity):
        """
        This init function used to initialize the Edge object
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input:flow,capacity
        Return:
        Time complexity:O(1)
        Best case analysis:O(1)
        Worst case analysis:O(1)
        Space complexity:O(1)
        Input space analysis: O(1)
        Aux space analysis:O(1)
        """
        self.flow=flow
        self.capacity=capacity
    def forward(self):
        """
        This forward function used to recalculate the edge for residual network for forward edge
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input:
        Return:Edge object
        Time complexity:O(1)
        Best case analysis:O(1)
        Worst case analysis:O(1)
        Space complexity:O(1)
        Input space analysis: O(1)
        Aux space analysis:O(1)
        """
        return Edge(self.capacity - self.flow, self.capacity)

    def reverse(self):
        """
        This reverse function used to recalculate the edge for residual network for reverse edge
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input:
        Return:Edge object
        Time complexity:O(1)
        Best case analysis:O(1)
        Worst case analysis:O(1)
        Space complexity:O(1)
        Input space analysis: O(1)
        Aux space analysis:O(1)
        """
        return Edge(self.flow, self.capacity)
    def forward_update(self,flow):
        """
        This forward update function used to assign new flow to previous forward flow if there is augmenting path
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input:flow
        Return:
        Time complexity:O(1)
        Best case analysis:O(1)
        Worst case analysis:O(1)
        Space complexity:O(1)
        Input space analysis: O(1)
        Aux space analysis:O(1)
        """
        self.flow=self.flow-flow
    def reverse_update(self,flow):
        """
        This reverse update function used to assign new flow to previous reverse flow if there is augmenting path
        Written by Chong Yi Ming
        Precondition:
        Postcondition:
        Input:flow
        Return:
        Time complexity:O(1)
        Best case analysis:O(1)
        Worst case analysis:O(1)
        Space complexity:O(1)
        Input space analysis: O(1)
        Aux space analysis:O(1)
        """
        self.flow=self.flow+flow

# g = FlowNetwork(6)
# g.insert_edge(0, 1, 0,16)
# g.insert_edge(0, 2, 0,13)
# g.insert_edge(1, 2, 0,10)
# g.insert_edge(2, 1, 0,4)
# g.insert_edge(1, 3, 0,12)
# g.insert_edge(3, 2,0, 9)
# g.insert_edge(2, 4, 0,14)
# g.insert_edge(4, 3, 0, 7)
# g.insert_edge(3, 5, 0,20)
# g.insert_edge(4, 5, 0,4)
# answer=g.ford_fulkerson()
# print(answer)
